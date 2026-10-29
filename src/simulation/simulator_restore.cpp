// Checkpoint restore / overlay. Reads a rank-count-independent
// checkpoint directory (state.h5 + per-rank delta shards + shard_index.yaml +
// manifest.yaml) and overlays the recorded state onto an already-loaded
// world by global id. Writer side lives in simulator_checkpoint.cpp.

#include <H5Cpp.h>
#include <yaml-cpp/yaml.h>

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "simulation/checkpoint_masks/policy_mask_io.h"
#include "simulation/simulator.h"

namespace june {
namespace {

namespace fs = std::filesystem;

template <typename T>
std::vector<T> readVec(H5::H5File& f, const std::string& n,
                       const H5::PredType& t) {
  if (!H5Lexists(f.getId(), n.c_str(), H5P_DEFAULT)) return {};
  H5::DataSet d = f.openDataSet(n);
  hsize_t dim[1];
  d.getSpace().getSimpleExtentDims(dim);
  std::vector<T> v(dim[0]);
  if (dim[0]) d.read(v.data(), t);
  return v;
}

std::vector<std::string> readStrs(H5::H5File& f, const std::string& n) {
  if (!H5Lexists(f.getId(), n.c_str(), H5P_DEFAULT)) return {};
  H5::DataSet d = f.openDataSet(n);
  hsize_t dim[1];
  d.getSpace().getSimpleExtentDims(dim);
  std::vector<std::string> out;
  if (!dim[0]) return out;
  H5::StrType st(H5::PredType::C_S1, H5T_VARIABLE);
  std::vector<char*> raw(dim[0]);
  d.read(raw.data(), st);
  for (hsize_t i = 0; i < dim[0]; ++i) out.emplace_back(raw[i] ? raw[i] : "");
  H5::DataSpace sp = d.getSpace();
  H5Dvlen_reclaim(st.getId(), sp.getId(), H5P_DEFAULT, raw.data());
  return out;
}

// Overlay one shard's per-rank manager state (last_processed_transition_time_
// and policy frozen_states) onto the accumulator maps, keeping only entries
// for people THIS rank owns. The maps are keyed on global PersonId so a
// resume at a different rank count still routes every entry correctly.
void overlayShardManagerState(
    H5::H5File& f, const WorldState& world,
    std::unordered_map<PersonId, double>& lpt_map,
    std::unordered_map<PersonId, FrozenPersonState>& frozen_accum,
    const std::vector<std::string>& follow_names,
    std::vector<std::unordered_map<PersonId, PersonId>>& follow_accum,
    std::vector<std::unordered_set<PersonId>>& active_host_accum) {
  const auto I32 = H5::PredType::NATIVE_INT32;
  const auto U8 = H5::PredType::NATIVE_UINT8;
  const auto F64 = H5::PredType::NATIVE_DOUBLE;
  auto lp = readVec<int32_t>(f, "/epidemiology/lpt_person_id", I32);
  auto lt = readVec<double>(f, "/epidemiology/lpt_time", F64);
  for (size_t i = 0; i < lp.size(); ++i)
    if (world.getPerson(lp[i])) lpt_map[lp[i]] = lt[i];

  auto fzp = readVec<int32_t>(f, "/policy_frozen_states/person_id", I32);
  auto fzpol =
      readVec<uint8_t>(f, "/policy_frozen_states/triggering_policy_index", U8);
  auto fzh = readVec<int32_t>(
      f, "/policy_frozen_states/paused_hopped_schedule_id", I32);
  auto fzr = readVec<int32_t>(
      f, "/policy_frozen_states/paused_return_schedule_id", I32);
  auto fzv = readVec<int32_t>(f, "/policy_frozen_states/pin_venue_id", I32);
  auto fzs = readVec<int32_t>(f, "/policy_frozen_states/pin_subset_index", I32);
  for (size_t i = 0; i < fzp.size(); ++i) {
    if (!world.getPerson(fzp[i])) continue;
    FrozenPersonState st;
    st.triggering_policy_index = fzpol[i];
    st.paused_hopped_schedule_id = static_cast<int16_t>(fzh[i]);
    st.paused_return_schedule_id = static_cast<int16_t>(fzr[i]);
    st.pin_venue_id = fzv[i];
    st.pin_subset_index = fzs[i];
    frozen_accum[fzp[i]] = st;
  }

  // Follows: keep a relationship if THIS rank owns the follower, and a
  // tried-host if it owns the host. Each stochastic rule has its own named
  // subgroup; a criteria rule wrote none and rebuilds itself. Guarded so
  // checkpoints predating the follow feature still restore.
  if (H5Lexists(f.getId(), "/follow", H5P_DEFAULT) > 0) {
    for (size_t ri = 0; ri < follow_names.size(); ++ri) {
      const std::string g = "/follow/" + follow_names[ri];
      if (H5Lexists(f.getId(), g.c_str(), H5P_DEFAULT) <= 0) continue;
      auto ff = readVec<int32_t>(f, g + "/follower_id", I32);
      auto fh = readVec<int32_t>(f, g + "/host_id", I32);
      for (size_t i = 0; i < ff.size(); ++i)
        if (world.getPerson(ff[i])) follow_accum[ri][ff[i]] = fh[i];
      auto ah = readVec<int32_t>(f, g + "/active_host_id", I32);
      for (int32_t h : ah)
        if (world.getPerson(h)) active_host_accum[ri].insert(h);
    }
  }
}

// Overlay one shard's /population/ records onto world_.people owned by this
// rank. Each owned person gets every field assigned and any default-
// constructed infection / vaccine_trajectory cleared (those get reinstated by
// the per-shard infection / vaccine overlays). Returns the number of records
// this rank applied. Applicable policies are interned into policy_manager's
// set table; with no PolicyManager they must be empty.
size_t overlayShardPopulation(H5::H5File& f, WorldState& world,
                              PolicyManager* policy_manager) {
  const auto I32 = H5::PredType::NATIVE_INT32;
  const auto U8 = H5::PredType::NATIVE_UINT8;
  const auto F64 = H5::PredType::NATIVE_DOUBLE;
  auto ids = readVec<int32_t>(f, "/population/ids", I32);
  auto imm_l = readVec<double>(f, "/population/immunity_natural_level", F64);
  auto imm_a =
      readVec<double>(f, "/population/immunity_natural_acquisition_time", F64);
  auto imm_w =
      readVec<double>(f, "/population/immunity_natural_waning_rate", F64);
  auto dead = readVec<uint8_t>(f, "/population/is_dead", U8);
  auto dt = readVec<double>(f, "/population/death_time", F64);
  auto m_as =
      readMasks(f, "/population/applicable_symptom_policy_mask", ids.size());
  auto m_at =
      readMasks(f, "/population/applicable_temporal_policy_mask", ids.size());
  auto m_ps = readMasks(f, "/population/active_symptom_policy_participation",
                        ids.size());
  auto m_pt = readMasks(f, "/population/active_temporal_policy_participation",
                        ids.size());
  auto m_ds = readMasks(f, "/population/symptom_policy_decisions", ids.size());
  auto m_dt = readMasks(f, "/population/temporal_policy_decisions", ids.size());
  auto hop = readVec<int32_t>(f, "/population/hopped_schedule_id", I32);
  auto ret = readVec<int32_t>(f, "/population/return_schedule_id", I32);
  auto tsl = readVec<int32_t>(f, "/population/temp_slot_progress", I32);
  auto hop_rep = readVec<int32_t>(f, "/population/hop_repeats_remaining", I32);
  size_t n = 0;
  for (size_t i = 0; i < ids.size(); ++i) {
    Person* p = world.getPerson(ids[i]);
    if (!p) continue;  // not owned by this rank
    ++n;
    p->immunity.natural_level = imm_l[i];
    p->immunity.natural_acquisition_time = imm_a[i];
    p->immunity.natural_waning_rate = imm_w[i];
    p->is_dead = dead[i] != 0;
    p->death_time = dt[i];
    if (policy_manager) {
      policy_manager->setApplicablePolicies(*p, m_as[i], m_at[i]);
    } else if (m_as[i].any() || m_at[i].any()) {
      throw std::runtime_error(
          "checkpoint holds applicable policies but the run has no policy "
          "manager");
    }
    p->active_symptom_policy_participation = m_ps[i];
    p->active_temporal_policy_participation = m_pt[i];
    p->symptom_policy_decisions = m_ds[i];
    p->temporal_policy_decisions = m_dt[i];
    p->schedule_hop.hopped_schedule_id = static_cast<int16_t>(hop[i]);
    p->schedule_hop.return_schedule_id = static_cast<int16_t>(ret[i]);
    p->schedule_hop.temp_slot_progress = static_cast<int16_t>(tsl[i]);
    p->schedule_hop.repeats_remaining = static_cast<int16_t>(hop_rep[i]);
    // Clear any default-constructed infection/vaccine; reinstated below.
    p->infection.reset();
    p->vaccine_trajectory.reset();
  }
  return n;
}

// Overlay one shard's /infection/ records onto world_.people owned by this
// rank. Each record rebuilds the symptom trajectory from the flat
// traj_times/traj_symptom_ids arrays, then materialises the Infection via
// Infection::fromCheckpoint so the disease config is rebound to the
// currently-loaded Disease pointer. Returns the number of records consumed.
size_t overlayShardInfection(H5::H5File& f, WorldState& world,
                             Disease* disease) {
  const auto I32 = H5::PredType::NATIVE_INT32;
  const auto I64 = H5::PredType::NATIVE_INT64;
  const auto U16 = H5::PredType::NATIVE_UINT16;
  const auto F64 = H5::PredType::NATIVE_DOUBLE;
  auto ip = readVec<int32_t>(f, "/infection/person_id", I32);
  auto it_ = readVec<double>(f, "/infection/infection_time", F64);
  auto mi = readVec<double>(f, "/infection/max_infectiousness", F64);
  auto sh = readVec<double>(f, "/infection/transmission_shape", F64);
  auto ra = readVec<double>(f, "/infection/transmission_rate", F64);
  auto sft = readVec<double>(f, "/infection/transmission_shift", F64);
  auto tt = readVec<double>(f, "/infection/traj_infection_time", F64);
  auto tf = readVec<double>(f, "/infection/traj_infectiousness_factor", F64);
  auto lc = readVec<double>(f, "/infection/last_checked_time", F64);
  auto csy = readVec<uint16_t>(f, "/infection/cached_symptom_id", U16);
  auto css = readVec<double>(f, "/infection/cached_symptom_start_time", F64);
  auto toff = readVec<int64_t>(f, "/infection/traj_offsets", I64);
  auto tcnt = readVec<int32_t>(f, "/infection/traj_counts", I32);
  auto ttime = readVec<double>(f, "/infection/traj_times", F64);
  auto tsym = readVec<uint16_t>(f, "/infection/traj_symptom_ids", U16);
  size_t n = 0;
  for (size_t i = 0; i < ip.size(); ++i) {
    Person* p = world.getPerson(ip[i]);
    if (!p) continue;
    InfectionTrajectory tr;
    tr.infection_time = tt[i];
    tr.infectiousness_factor = tf[i];
    int64_t off = toff[i];
    int32_t cnt = tcnt[i];
    tr.transitions.reserve(cnt);
    for (int32_t k = 0; k < cnt; ++k)
      tr.transitions.emplace_back(ttime[off + k], tsym[off + k]);
    p->infection =
        Infection::fromCheckpoint(disease, it_[i], tr, mi[i], sh[i], ra[i],
                                  sft[i], lc[i], csy[i], css[i]);
    ++n;
  }
  return n;
}

// Overlay one shard's /vaccine/ records onto world_.people owned by this
// rank. Efficacy is NOT persisted in the shard; it is re-derived from the
// active VaccinationConfig by dose position so a resume can apply updated
// efficacy maps. Returns the number of records this rank consumed.
size_t overlayShardVaccine(H5::H5File& f, WorldState& world,
                           const VaccinationConfig& vax_cfg) {
  const auto I32 = H5::PredType::NATIVE_INT32;
  const auto I64 = H5::PredType::NATIVE_INT64;
  const auto F64 = H5::PredType::NATIVE_DOUBLE;
  auto vp = readVec<int32_t>(f, "/vaccine/person_id", I32);
  auto vname = readStrs(f, "/vaccine/vaccine_name");
  auto doff = readVec<int64_t>(f, "/vaccine/dose_offsets", I64);
  auto dcnt = readVec<int32_t>(f, "/vaccine/dose_counts", I32);
  auto dnum = readVec<int32_t>(f, "/vaccine/dose_number", I32);
  auto dad = readVec<double>(f, "/vaccine/dose_day_administered", F64);
  auto deff = readVec<double>(f, "/vaccine/dose_days_to_effective", F64);
  auto dwan = readVec<double>(f, "/vaccine/dose_days_to_waning", F64);
  auto dfin = readVec<double>(f, "/vaccine/dose_days_to_finished", F64);
  auto dwf = readVec<double>(f, "/vaccine/dose_waning_factor", F64);
  size_t n = 0;
  for (size_t i = 0; i < vp.size(); ++i) {
    Person* p = world.getPerson(vp[i]);
    if (!p) continue;
    p->vaccine_trajectory = std::make_unique<VaccineTrajectory>();
    p->vaccine_trajectory->vaccine_name = vname[i];
    auto vit = vax_cfg.vaccines.find(vname[i]);
    int64_t off = doff[i];
    int32_t cnt = dcnt[i];
    for (int32_t k = 0; k < cnt; ++k) {
      Dose d;
      d.number = dnum[off + k];
      d.day_administered = dad[off + k];
      d.days_to_effective = deff[off + k];
      d.days_to_waning = dwan[off + k];
      d.days_to_finished = dfin[off + k];
      d.waning_factor = dwf[off + k];
      if (vit != vax_cfg.vaccines.end() &&
          k < static_cast<int32_t>(vit->second.doses.size())) {
        d.infection_efficacy = vit->second.doses[k].infection_efficacy;
        d.symptom_efficacy = vit->second.doses[k].symptom_efficacy;
      }
      p->vaccine_trajectory->addDose(d);
    }
    ++n;
  }
  return n;
}

// Overlay one shard's /venue_fomite/ records onto world_.venues owned by this
// rank. Each record names (venue_id, mode, offset, count); the mode-keyed
// deque is resized + cleared + repopulated from the flat deposit arrays.
// Returns the number of records this rank consumed (for diagnostics).
size_t overlayShardFomite(H5::H5File& f, WorldState& world) {
  const auto I32 = H5::PredType::NATIVE_INT32;
  const auto I64 = H5::PredType::NATIVE_INT64;
  const auto F64 = H5::PredType::NATIVE_DOUBLE;
  auto fv = readVec<int32_t>(f, "/venue_fomite/venue_id", I32);
  auto fm = readVec<int32_t>(f, "/venue_fomite/mode_index", I32);
  auto fo = readVec<int64_t>(f, "/venue_fomite/offsets", I64);
  auto fc = readVec<int32_t>(f, "/venue_fomite/counts", I32);
  auto ftime = readVec<double>(f, "/venue_fomite/deposit_time", F64);
  auto famt = readVec<double>(f, "/venue_fomite/deposit_amount", F64);
  size_t n = 0;
  for (size_t i = 0; i < fv.size(); ++i) {
    Venue* v = world.getVenue(fv[i]);
    if (!v) continue;
    int m = fm[i];
    if (static_cast<int>(v->fomite_history.size()) <= m)
      v->fomite_history.resize(m + 1);
    auto& dq = v->fomite_history[m];
    dq.clear();
    for (int32_t k = 0; k < fc[i]; ++k)
      dq.push_back({ftime[fo[i] + k], famt[fo[i] + k]});
    ++n;
  }
  return n;
}

void overlayShardCalendarEvents(H5::H5File& f, WorldState& world,
                                CalendarEventManager::Snapshot& snap) {
  if (!f.exists("/calendar_events")) return;
  const auto I32 = H5::PredType::NATIVE_INT32;
  auto pids = readVec<int32_t>(f, "/calendar_events/person_ids", I32);
  auto eids = readVec<int32_t>(f, "/calendar_events/event_ids", I32);
  for (size_t i = 0; i < pids.size(); ++i)
    if (world.getPerson(static_cast<PersonId>(pids[i])))
      snap.active_event[static_cast<PersonId>(pids[i])] = eids[i];
}

// Open the checkpoint's manifest.yaml, enforce that the checkpoint's seed
// matches the configured seed (no silent override), and return the
// completed_day to resume from. The caller has already canonicalised cp.
// format_version comes back through the out parameter; one outside the range
// this build reads throws.
int readAndValidateManifest(const fs::path& cp, unsigned int config_seed,
                            int& format_version) {
  if (!fs::exists(cp / "manifest.yaml"))
    throw std::runtime_error(
        "restoreFromCheckpoint: no manifest.yaml (incomplete/invalid "
        "checkpoint): " +
        cp.string());
  YAML::Node man = YAML::LoadFile((cp / "manifest.yaml").string());
  format_version = man["format_version"].as<int>();
  if (format_version < 1 || format_version > 2) {
    throw std::runtime_error(
        "restoreFromCheckpoint: unsupported checkpoint format_version " +
        std::to_string(format_version) +
        " (this build reads 1 to 2): " + cp.string());
  }
  int completed_day = man["completed_day"].as<int>();
  unsigned int cp_seed = man["effective_random_seed"].as<unsigned int>();
  if (cp_seed != config_seed) {
    throw std::runtime_error(
        "restoreFromCheckpoint: effective_random_seed mismatch (checkpoint=" +
        std::to_string(cp_seed) + ", config=" + std::to_string(config_seed) +
        "). Resume with the checkpoint's seed (no silent override).");
  }
  return completed_day;
}

// Observed metric totals stored in state.h5, name order. Absent group means
// the checkpoint was written with no active metric. The rows of the unit table
// are (metric index, unit, value).
std::vector<MetricSnapshot> readStoredMetrics(const fs::path& cp) {
  H5::H5File f((cp / "state.h5").string(), H5F_ACC_RDONLY);
  if (H5Lexists(f.getId(), "/observed_metrics", H5P_DEFAULT) <= 0) return {};
  const auto names = readStrs(f, "/observed_metrics/names");
  const auto values = readVec<int64_t>(f, "/observed_metrics/values",
                                       H5::PredType::NATIVE_INT64);
  const auto scopes = readVec<int64_t>(f, "/observed_metrics/scopes",
                                       H5::PredType::NATIVE_INT64);
  const auto rows = readVec<int64_t>(f, "/observed_metrics/unit_rows",
                                     H5::PredType::NATIVE_INT64);
  if (values.size() != names.size() || scopes.size() != names.size() ||
      rows.size() % 3 != 0) {
    throw std::runtime_error(
        "restoreFromCheckpoint: observed_metrics arrays differ in length in " +
        (cp / "state.h5").string());
  }
  std::vector<MetricSnapshot> out(names.size());
  for (size_t i = 0; i < names.size(); ++i) {
    out[i].name = names[i];
    out[i].scope = static_cast<MetricScope>(scopes[i]);
    out[i].value = values[i];
  }
  for (size_t r = 0; r < rows.size(); r += 3) {
    if (rows[r] < 0 || static_cast<size_t>(rows[r]) >= names.size()) {
      throw std::runtime_error(
          "restoreFromCheckpoint: observed_metrics unit row names metric " +
          std::to_string(rows[r]) + " but the checkpoint holds " +
          std::to_string(names.size()) + " in " + (cp / "state.h5").string());
    }
    out[rows[r]].units.emplace_back(static_cast<VenueId>(rows[r + 1]),
                                    rows[r + 2]);
  }
  return out;
}

// Policy trigger state stored in state.h5. Absent group means the checkpoint
// was written with no triggered Policy. The rows of the unit table are
// (trigger index, unit, closed, anchor day, anchor value).
std::vector<PolicyManager::TriggerState> readStoredTriggers(
    const fs::path& cp) {
  H5::H5File f((cp / "state.h5").string(), H5F_ACC_RDONLY);
  if (H5Lexists(f.getId(), "/policy_triggers", H5P_DEFAULT) <= 0) return {};
  const auto keys = readStrs(f, "/policy_triggers/keys");
  const auto state =
      readVec<int64_t>(f, "/policy_triggers/state", H5::PredType::NATIVE_INT64);
  const auto scopes = readVec<int64_t>(f, "/policy_triggers/scopes",
                                       H5::PredType::NATIVE_INT64);
  const auto voluntary = readVec<int64_t>(f, "/policy_triggers/voluntary",
                                          H5::PredType::NATIVE_INT64);
  const auto rows = readVec<int64_t>(f, "/policy_triggers/unit_rows",
                                     H5::PredType::NATIVE_INT64);
  if (state.size() != 3 * keys.size() || scopes.size() != keys.size() ||
      voluntary.size() != keys.size() || rows.size() % 5 != 0) {
    throw std::runtime_error(
        "restoreFromCheckpoint: policy_triggers arrays differ in length in " +
        (cp / "state.h5").string());
  }
  std::vector<PolicyManager::TriggerState> out(keys.size());
  for (size_t i = 0; i < keys.size(); ++i) {
    out[i].key = keys[i];
    out[i].cell = {state[3 * i] != 0, state[3 * i + 1], state[3 * i + 2]};
    out[i].scope = static_cast<MetricScope>(scopes[i]);
    out[i].voluntary = voluntary[i] != 0;
  }
  for (size_t r = 0; r < rows.size(); r += 5) {
    if (rows[r] < 0 || static_cast<size_t>(rows[r]) >= keys.size()) {
      throw std::runtime_error(
          "restoreFromCheckpoint: policy_triggers unit row names trigger " +
          std::to_string(rows[r]) + " but the checkpoint holds " +
          std::to_string(keys.size()) + " in " + (cp / "state.h5").string());
    }
    out[rows[r]].units.emplace_back(
        static_cast<VenueId>(rows[r + 1]),
        TriggerCell{rows[r + 2] != 0, rows[r + 3], rows[r + 4]});
  }
  return out;
}

}  // namespace

void Simulator::validateResumeBounds(int completed_day) const {
  if (resume_from_day_ < total_days_) return;
  throw std::runtime_error(
      "checkpoint resume: requested simulation length leaves nothing to "
      "simulate. Checkpoint completed day " +
      std::to_string(completed_day) + " so the run resumes at day " +
      std::to_string(resume_from_day_) +
      ", but the configured end is total_days=" + std::to_string(total_days_) +
      " (end_date=" + config_.simulation.end_date +
      "). --days/end_date counts from the original start_date, not from "
      "the checkpoint. To run M more day(s) after this checkpoint, pass "
      "--days " +
      std::to_string(completed_day + 2) +
      " or greater (= completed_day + 1 + M).");
}

void Simulator::restoreCheckpointStateFile(const fs::path& cp) {
  H5::H5File f((cp / "state.h5").string(), H5F_ACC_RDONLY);
  auto F64 = H5::PredType::NATIVE_DOUBLE;
  auto I32 = H5::PredType::NATIVE_INT32;
  auto cst = readVec<double>(f, "/scalars/current_simulation_time", F64);
  if (!cst.empty()) current_simulation_time_ = cst[0];
  auto dtc = readVec<int32_t>(f, "/scalars/day_type_counts", I32);
  day_type_counts_.assign(dtc.begin(), dtc.end());
}

void Simulator::restoreFromCheckpoint(const std::string& checkpoint_dir) {
  const int rank = getRank();
  fs::path cp = fs::canonical(checkpoint_dir);  // resolves 'latest' symlink
  int format_version = 0;
  int completed_day = readAndValidateManifest(
      cp, config_.simulation.random_seed, format_version);

  restoreCheckpointStateFile(cp);

  // A cumulative counter cannot be rebuilt, so a checkpoint that does not
  // hold exactly the active metrics is refused rather than restarted at zero.
  // Every triggered Policy reads an active metric, so format 1, which holds
  // none, only resumes a run with no triggers either.
  if (policy_manager_) {
    auto& metrics = policy_manager_->observedMetrics();
    if (format_version == 1) {
      if (metrics.anyActive()) {
        throw std::runtime_error(
            "restoreFromCheckpoint: checkpoint is format_version 1 and holds "
            "no observed metric totals, but this run has active metrics. "
            "Resume with no triggered Policy, or restart from day 0.");
      }
    } else {
      metrics.restoreSnapshot(readStoredMetrics(cp));
      policy_manager_->restoreTriggerStates(readStoredTriggers(cp));
      // The levels are derived from the metric counts just put back.
      policy_manager_->rebuildResponseLevels();
    }
  }

  // frozen_states_ + lpt + follows are per-rank: accumulated from the shards.
  // Follows accumulate one entry per rule, matched to its shard by name.
  const auto& follow_rules = config_.coordinated_encounters.follows;
  std::vector<std::string> follow_names;
  for (const auto& r : follow_rules) follow_names.push_back(r.name);
  std::unordered_map<PersonId, double> lpt_map;
  std::unordered_map<PersonId, FrozenPersonState> frozen_accum;
  std::vector<std::unordered_map<PersonId, PersonId>> follow_accum(
      follow_rules.size());
  std::vector<std::unordered_set<PersonId>> active_host_accum(
      follow_rules.size());

  // ---- delta shards: overlay onto this rank's owned world_ by global id ----
  YAML::Node si = YAML::LoadFile((cp / "shard_index.yaml").string());
  int n_shards = si["num_ranks"].as<int>();
  size_t n_people = 0, n_inf = 0, n_vax = 0, n_fom = 0;
  CalendarEventManager::Snapshot ce_snap;
  for (int s = 0; s < n_shards; ++s) {
    fs::path sf = cp / ("delta_rank" + std::to_string(s) + ".h5");
    if (!fs::exists(sf)) continue;
    H5::H5File f(sf.string(), H5F_ACC_RDONLY);
    n_people += overlayShardPopulation(f, world_, policy_manager_.get());
    n_inf += overlayShardInfection(f, world_, disease_.get());
    n_vax += overlayShardVaccine(f, world_, config_.vaccination);
    n_fom += overlayShardFomite(f, world_);
    overlayShardManagerState(f, world_, lpt_map, frozen_accum, follow_names,
                             follow_accum, active_host_accum);
    overlayShardCalendarEvents(f, world_, ce_snap);
  }
  calendar_event_manager_.restore(std::move(ce_snap));

  if (policy_manager_) policy_manager_->setFrozenStates(frozen_accum);
  // Rebuild per-rule runtime state from the accumulated shards. A criteria rule
  // has empty accumulators and derives itself on the first slot after resume;
  // its follow_day stays -1 so that rebuild fires.
  follow_state_.assign(follow_rules.size(), FollowRuntime{});
  for (size_t ri = 0; ri < follow_rules.size(); ++ri) {
    follow_state_[ri].follower_host = std::move(follow_accum[ri]);
    follow_state_[ri].active_hosts = std::move(active_host_accum[ri]);
  }

  // ---- rebuild derived caches; set resume point ----
  if (epidemiology_) epidemiology_->restoreAfterCheckpoint(lpt_map);
  if (policy_manager_ && disease_)
    policy_manager_->initializeTransmissionModifiers(*disease_,
                                                     current_simulation_time_);
  resume_from_day_ = completed_day + 1;

  validateResumeBounds(completed_day);

  if (rank == 0) {
    std::cout << "[checkpoint] restored from " << cp.string()
              << ": day=" << completed_day
              << " sim_time=" << current_simulation_time_
              << " people=" << n_people << " infections=" << n_inf
              << " vaccinated=" << n_vax << " fomite_venues=" << n_fom
              << " — resuming at day " << resume_from_day_ << std::endl;
  }
}

}  // namespace june
