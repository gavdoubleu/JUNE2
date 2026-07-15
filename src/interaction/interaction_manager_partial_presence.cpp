// Contains the processPartialPresenceVenue family (commute-line FOI path).
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>

#include "activity/presence_window.h"
#include "activity/runtime_bin_allocator.h"
#include "epidemiology/interaction_manager.h"
#include "simulation/compartmental_model_manager.h"
#include "utils/deterministic_rng.h"
#include "utils/event_logging/event_types.h"
#include "utils/random.h"

#ifdef USE_MPI
#include <mpi.h>

#include "parallel/mpi_utils.h"
#endif

namespace june {

// Routing gate called from processVenueTransmissions (venue.cpp). Returns a
// new-infection count (forwarded from processPartialPresenceVenue) when the
// venue type is declared in SimulationConfig::partial_presence; otherwise
// returns nullopt and the caller falls through to the standard FOI path.
std::optional<int> InteractionManager::dispatchPartialPresenceIfApplicable(
    const std::vector<InteractionMember>& members, Venue* venue,
    VenueId actual_venue_id, double current_time, double delta_hours,
    std::unordered_set<PersonId>* active_infections,
    const std::unordered_set<PersonId>* visitor_ids,
    std::vector<PendingInfection>* pending_infections,
    const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
    uint8_t encounter_type_id, const CompartmentalModelManager* comp_model) {
  // Partial-presence venues (commute lines, etc.) take a different FOI path:
  // sub-interval-aware, carriage-grouped, with per-rider effective presence
  // windows. The gate is a single bit-mask test; venue types not declared in
  // SimulationConfig::partial_presence pay no cost here.
  if (!runtime_bin_allocator_ || actual_venue_id < 0 || !venue) return {};
  const uint8_t vt = venue->type_id;
  const uint64_t mask =
      simulation_config_.partial_presence.enabled_venue_type_mask;
  if (vt >= 64 || ((mask >> vt) & 1ULL) == 0) return {};
  return processPartialPresenceVenue(
      members, venue, actual_venue_id, current_time, delta_hours,
      active_infections, visitor_ids, pending_infections, visitor_data,
      encounter_type_id, comp_model);
}

void InteractionManager::accumulateOneCarriage(
    const std::vector<CarriageMember>& car, float slot_duration_min,
    double current_time, double delta_hours, int num_modes, int num_bins_needed,
    uint8_t venue_type_id, const ContactMatrix* matrix,
    const TransmissionParams& trans_params,
    std::vector<PartialPresenceSubBin>& sub_bins,
    PartialPresenceLambdaResult& result) const {
  std::vector<float> events =
      collectSubIntervalEventTimes(car, slot_duration_min);
  if (events.size() < 2) return;

  for (size_t si = 0; si + 1 < events.size(); ++si) {
    const float t0 = events[si];
    const float t1 = events[si + 1];
    const float sub_dur = t1 - t0;
    if (!(sub_dur > 0.0f)) continue;
    const double scale = static_cast<double>(sub_dur) / slot_duration_min;

    for (auto& sb : sub_bins) sb.reset(num_modes);

    // Track per-bin susceptibles in this sub-interval (rebuilt each pass).
    std::vector<std::vector<const CarriageMember*>> susc_by_bin(
        num_bins_needed);

    classifyMembersInSubInterval(car, t0, t1, scale, current_time, delta_hours,
                                 num_modes, sub_bins, susc_by_bin);

    accumulatePartialLambdaContributions(sub_bins, susc_by_bin, venue_type_id,
                                         matrix, num_bins_needed, num_modes,
                                         trans_params, result);
  }
}

void InteractionManager::accumulatePartialLambdaContributions(
    const std::vector<PartialPresenceSubBin>& sub_bins,
    const std::vector<std::vector<const CarriageMember*>>& susc_by_bin,
    uint8_t venue_type_id, const ContactMatrix* matrix, int num_bins_needed,
    int num_modes, const TransmissionParams& trans_params,
    PartialPresenceLambdaResult& result) const {
  using AccumSource = PartialPresenceAccumSource;
  auto& susc_lambda = result.susc_lambda;
  auto& susc_sources = result.susc_sources;

  for (int susc_bin = 0; susc_bin < num_bins_needed; ++susc_bin) {
    if (susc_by_bin[susc_bin].empty()) continue;

    for (int mode = 0; mode < num_modes; ++mode) {
      const ContactMatrix* mode_matrix =
          contact_matrices_.getMatrix(venue_type_id, mode);
      double mode_susc_mult =
          (mode < static_cast<int>(trans_params.modes.size()))
              ? trans_params.modes[mode].susceptibility_multiplier
              : 1.0;

      for (int inf_bin = 0; inf_bin < num_bins_needed; ++inf_bin) {
        double total_inf = sub_bins[inf_bin].total_inf_by_mode[mode];
        if (!(total_inf > 0.0)) continue;

        double contacts =
            lookupContactsForBinPair(mode_matrix, matrix, susc_bin, inf_bin);
        if (!(contacts > 0.0)) continue;

        int bin_size = sub_bins[inf_bin].total_size;
        if (susc_bin == inf_bin) bin_size = std::max(1, bin_size - 1);
        double omega = contacts / bin_size;
        double contrib = omega * total_inf * mode_susc_mult;
        if (!(contrib > 0.0)) continue;

        for (const CarriageMember* sm : susc_by_bin[susc_bin]) {
          // f_S: this susceptible's own presence cap (1.0 unless over-long).
          // Channel-symmetric with f_I (baked into total_inf via inf_sub): a
          // contact carries f_I·f_S in both directions. The source weights `w`
          // below intentionally omit f_S — it is a per-susceptible constant, so
          // it cannot change the relative infector sampling for this
          // susceptible.
          susc_lambda[sm->pid] += contrib * sm->f_presence;
          // Pick an infector for this contribution by weight-sampling
          // proportional to per-person infectiousness in this sub-interval.
          // We record one AccumSource per (susc, mode, inf_bin, sub) with
          // the FULL bin contribution; per-person sampling happens at
          // the single Bernoulli site below.
          const auto& ids = sub_bins[inf_bin].infectious_ids;
          const auto& per = sub_bins[inf_bin].inf_per_person_by_mode[mode];
          for (size_t pi = 0; pi < ids.size(); ++pi) {
            if (pi >= per.size() || !(per[pi] > 0.0)) continue;
            double w = omega * per[pi] * mode_susc_mult;
            if (!(w > 0.0)) continue;
            susc_sources[sm->pid].push_back(AccumSource{mode, ids[pi], w});
          }
        }
      }
    }
  }
}

void InteractionManager::classifyMembersInSubInterval(
    const std::vector<CarriageMember>& car, float t0, float t1, double scale,
    double current_time, double delta_hours, int num_modes,
    std::vector<PartialPresenceSubBin>& sub_bins,
    std::vector<std::vector<const CarriageMember*>>& susc_by_bin) const {
  for (const auto& m : car) {
    // Present iff [eff_board, eff_alight) covers [t0, t1).
    if (!(m.eff_board <= t0 + 1e-5f && m.eff_alight + 1e-5f >= t1)) continue;

    const int bin = m.matrix_bin;
    const bool dead = (m.person && m.person->is_dead);
    if (!dead) sub_bins[bin].total_size++;

    // Infectious?
    bool added_inf = false;
    if (m.visitor && m.visitor->is_infectious) {
      for (int mode = 0; mode < num_modes; ++mode) {
        double inf_full = (mode < VisitorInfo::MAX_MODES)
                              ? m.visitor->integrated_infectiousness[mode]
                              : 0.0;
        // f_I: this infectious rider's presence cap (1.0 unless over-long).
        double inf_sub = inf_full * scale * m.f_presence;
        if (inf_sub > 0.0) {
          if (!added_inf) {
            sub_bins[bin].infectious_ids.push_back(m.pid);
            added_inf = true;
          }
          sub_bins[bin].inf_per_person_by_mode[mode].push_back(inf_sub);
          sub_bins[bin].total_inf_by_mode[mode] += inf_sub;
        } else if (added_inf) {
          // Keep arrays aligned across modes.
          sub_bins[bin].inf_per_person_by_mode[mode].push_back(0.0);
        }
      }
    } else if (m.person && m.person->infection &&
               m.person->infection->isInfectious(current_time)) {
      const double t_end_d = current_time + delta_hours / 24.0;
      for (int mode = 0; mode < num_modes; ++mode) {
        double inf_full = m.person->infection->getIntegratedInfectiousness(
            mode, current_time, t_end_d);
        // f_I: this infectious rider's presence cap (1.0 unless over-long).
        double inf_sub = inf_full * scale * m.f_presence;
        if (inf_sub > 0.0) {
          if (!added_inf) {
            sub_bins[bin].infectious_ids.push_back(m.pid);
            added_inf = true;
          }
          sub_bins[bin].inf_per_person_by_mode[mode].push_back(inf_sub);
          sub_bins[bin].total_inf_by_mode[mode] += inf_sub;
        } else if (added_inf) {
          sub_bins[bin].inf_per_person_by_mode[mode].push_back(0.0);
        }
      }
    } else if (m.person && !m.person->infection && !dead) {
      double susc =
          m.person->getSusceptibility(current_time, disease_->getName());
      if (susc > 0.0) susc_by_bin[bin].push_back(&m);
    } else if (m.visitor && !m.visitor->is_infected &&
               m.visitor->immunity_level < 1.0) {
      susc_by_bin[bin].push_back(&m);
    }
  }
}

std::vector<float> InteractionManager::collectSubIntervalEventTimes(
    const std::vector<CarriageMember>& car, float slot_duration_min) const {
  // Slice sub-intervals over the venue's OWN [min board, max alight] range
  // (NOT [0, slot]): the windows are raw line-local offsets in this venue's
  // clock, so the boundaries are exactly the members' board/alight points.
  // Forcing 0 and slot as bounds would (a) fabricate empty edge intervals and
  // (b) wrongly truncate riders whose raw offsets exceed the slot. The
  // per-rider presence cap f_p — not a slot clamp — conserves the day budget.
  (void)slot_duration_min;
  std::vector<float> events;
  events.reserve(2 * car.size());
  for (const auto& m : car) {
    events.push_back(m.eff_board);
    events.push_back(m.eff_alight);
  }
  std::sort(events.begin(), events.end());
  events.erase(
      std::unique(events.begin(), events.end(),
                  [](float a, float b) { return std::abs(a - b) < 1e-5f; }),
      events.end());
  return events;
}

std::vector<std::vector<CarriageMember>>
InteractionManager::buildPartialPresenceCarriages(
    const std::vector<InteractionMember>& members, Venue* venue,
    VenueId actual_venue_id, const ContactMatrix* matrix, int num_bins_needed,
    uint16_t num_bins,
    const std::unordered_map<PersonId, VisitorInfo>* visitor_data) const {
  std::vector<std::vector<CarriageMember>> carriages(num_bins);

  for (const auto& m : members) {
    const uint16_t carriage =
        runtime_bin_allocator_->getBinIndex(actual_venue_id, m.id);
    // Members come from the rider table, so a missing or out-of-range carriage
    // means the two disagree. Skipping here is what used to put people on a
    // line with no force of infection in either direction, so it throws.
    if (carriage == RuntimeBinAllocator::kNoBin || carriage >= num_bins)
      throw std::runtime_error("partial-presence venue " +
                               std::to_string(actual_venue_id) + ": rider " +
                               std::to_string(m.id) + " has no carriage (got " +
                               std::to_string(static_cast<int>(carriage)) +
                               " of " + std::to_string(num_bins) + ")");

    Person* person = nullptr;
    const VisitorInfo* visitor = nullptr;
    if (!resolvePersonAndVisitor(m.id, m.array_index, visitor_data, person,
                                 visitor)) {
      throw std::runtime_error(
          "partial-presence venue " + std::to_string(actual_venue_id) +
          ": rider " + std::to_string(m.id) +
          " is neither a local person nor a visitor on this rank. The visitor "
          "exchange should have shipped every rider of a line this rank owns.");
    }

    // Window + presence cap from the allocator's global broadcast: identical
    // on every rank for the same (venue, person) pair.
    const EffectiveWindow win =
        runtime_bin_allocator_->getPresenceWindow(actual_venue_id, m.id);
    const float f_presence =
        runtime_bin_allocator_->getPresenceFactor(actual_venue_id, m.id);

    int matrix_bin = computeBinIndexForMatrix(person, venue, m.subset_index,
                                              matrix, num_bins_needed);
    if (matrix_bin < 0 || matrix_bin >= num_bins_needed) matrix_bin = 0;

    carriages[carriage].push_back(CarriageMember{
        m.id, m.array_index, m.subset_index, m.encounter_type_id, person,
        visitor, win.eff_board, win.eff_alight, f_presence, matrix_bin});
  }

  // Deterministic per-carriage order (FP-stable accumulation across ranks).
  for (auto& car : carriages) {
    std::sort(car.begin(), car.end(),
              [](const CarriageMember& a, const CarriageMember& b) {
                return a.pid < b.pid;
              });
  }
  return carriages;
}

void InteractionManager::validatePartialPresencePreconditions(
    const std::vector<InteractionMember>& members, const Venue* venue,
    VenueId actual_venue_id, uint8_t encounter_type_id) const {
  if (actual_venue_id < 0)
    throw std::runtime_error(
        "computePartialPresenceLambda: virtual encounter venues not supported");
  if (!venue)
    throw std::runtime_error("computePartialPresenceLambda: null venue");
  (void)members;
  if (encounter_type_id != kDefaultEncounterTypeId)
    throw std::runtime_error(
        "computePartialPresenceLambda: coordinated-encounter venues not "
        "supported on partial-presence types in v1");
  if (venue->parent_id >= 0)
    throw std::runtime_error(
        "computePartialPresenceLambda: parent-venue mixing not supported on "
        "partial-presence types in v1");
  if (!runtime_bin_allocator_)
    throw std::runtime_error(
        "computePartialPresenceLambda: runtime_bin_allocator_ is null (gate "
        "should have prevented this call)");
}

InteractionManager::PartialPresenceLambdaResult
InteractionManager::computePartialPresenceLambda(
    const std::vector<InteractionMember>& members, Venue* venue,
    VenueId actual_venue_id, double current_time, double delta_hours,
    const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
    uint8_t encounter_type_id) {
  PartialPresenceLambdaResult result;

  validatePartialPresencePreconditions(members, venue, actual_venue_id,
                                       encounter_type_id);

  const uint8_t venue_type_id = venue->type_id;
  const ContactMatrix* matrix = contact_matrices_.getMatrix(venue_type_id);
  if (!matrix)
    throw std::runtime_error(
        "computePartialPresenceLambda: no contact matrix for venue_type_id=" +
        std::to_string(static_cast<int>(venue_type_id)));
  const int num_bins_needed =
      std::max(1, static_cast<int>(matrix->bins.size()));

  int num_modes = disease_->numModes();
  if (num_modes == 0) num_modes = 1;
  const auto& trans_params = disease_->getTransmissionParams();

  const float slot_duration_min = static_cast<float>(delta_hours * 60.0);
  if (!(slot_duration_min > 0.0f)) return result;

  // Presence windows + per-rider caps f_p live on the allocator. It computes
  // them on each rider's home rank (raw line-local windows; f_p from the full
  // leg list) and broadcasts globally, so a cross-rank visitor's window and
  // f_p are identical to what the home rank would have computed locally.
  const uint16_t num_bins = runtime_bin_allocator_->getNumBins(actual_venue_id);
  if (num_bins == 0) return result;

  std::vector<std::vector<CarriageMember>> carriages =
      buildPartialPresenceCarriages(members, venue, actual_venue_id, matrix,
                                    num_bins_needed, num_bins, visitor_data);

  // Per-bin scratch reused across sub-intervals (cleared per sub-interval).
  std::vector<PartialPresenceSubBin> sub_bins(num_bins_needed);

  for (uint16_t c = 0; c < num_bins; ++c) {
    const auto& car = carriages[c];
    if (car.empty()) continue;
    accumulateOneCarriage(car, slot_duration_min, current_time, delta_hours,
                          num_modes, num_bins_needed, venue_type_id, matrix,
                          trans_params, sub_bins, result);
  }

  return result;
}

void InteractionManager::recordPartialPresenceCandidate(
    PersonId susc_id, PersonId infector_id, uint8_t transmission_mode_index,
    uint16_t infector_symptom_id, double current_time, uint8_t venue_type_id,
    VenueId actual_venue_id) {
  // Nothing is applied here. A rider is susceptible on every leg of their
  // journey at once, and the legs can be owned by different ranks, so infecting
  // them the moment one leg says so would let the surviving leg depend on the
  // order the venues happened to be visited in. Each leg records what it would
  // do; resolvePartialPresenceInfections then picks one winner per person.
  PendingInfection cand;
  cand.person_id = susc_id;
  cand.infector_id = infector_id;
  cand.infection_time = current_time;
  cand.venue_type_id = venue_type_id;
  cand.encounter_type_id = kDefaultEncounterTypeId;
  cand.venue_id = actual_venue_id;
  cand.infector_symptom_id = infector_symptom_id;
  cand.transmission_mode_index = transmission_mode_index;
  pp_candidates_.push_back(cand);
}

int InteractionManager::processPartialPresenceLines(
    const std::vector<VenueId>& owned_lines, double current_time,
    double delta_hours, std::unordered_set<PersonId>* active_infections,
    const std::unordered_map<PersonId, VisitorInfo>* visitor_data) {
  pp_candidates_.clear();
  if (!runtime_bin_allocator_ || owned_lines.empty()) return 0;

  std::vector<InteractionMember> members;
  for (VenueId vid : owned_lines) {
    Venue* venue = world_.getVenue(vid);
    if (!venue)
      throw std::runtime_error("partial-presence pass: rank owns line " +
                               std::to_string(vid) +
                               " but the venue is not in its world");

    const auto& riders = runtime_bin_allocator_->ridersByVenue().at(vid);
    members.clear();
    members.reserve(riders.size());
    for (const auto& r : riders) {
      // array_index is a local hint; a rider from another rank is resolved
      // through visitor_data instead, and resolvePersonAndVisitor tries both.
      size_t idx = static_cast<size_t>(-1);
      auto it = world_.person_index.find(r.pid);
      if (it != world_.person_index.end()) idx = it->second;
      members.push_back(
          InteractionMember{r.pid, idx, r.subset, kDefaultEncounterTypeId});
    }

    processPartialPresenceVenue(members, venue, vid, current_time, delta_hours,
                                active_infections, nullptr, nullptr,
                                visitor_data, kDefaultEncounterTypeId, nullptr);
  }
  return 0;  // counted once the winners are resolved
}

int InteractionManager::resolvePartialPresenceInfections(
    double current_time, std::unordered_set<PersonId>* active_infections) {
  // Legs of one journey can be owned by different ranks, so a rider's
  // candidates are scattered. Gather them all, then every rank sorts the same
  // list and reaches the same verdict; each applies only the people it owns.
  std::vector<int32_t> packed;
  packed.reserve(pp_candidates_.size() * 6);
  for (const auto& c : pp_candidates_) {
    packed.push_back(static_cast<int32_t>(c.person_id));
    packed.push_back(static_cast<int32_t>(c.venue_id));
    packed.push_back(static_cast<int32_t>(c.infector_id));
    packed.push_back(static_cast<int32_t>(c.venue_type_id));
    packed.push_back(static_cast<int32_t>(c.infector_symptom_id));
    packed.push_back(static_cast<int32_t>(c.transmission_mode_index));
  }

#ifdef USE_MPI
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size > 1) {
    int local_count = static_cast<int>(packed.size());
    std::vector<int> counts(world_size);
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    std::vector<int> displs;
    int total = 0;
    mpi_utils::computeDisplacements(counts, displs, total);
    std::vector<int32_t> global(total, 0);
    MPI_Allgatherv(packed.data(), local_count, MPI_INT, global.data(),
                   counts.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);
    packed.swap(global);
  }
#endif

  std::vector<PendingInfection> all;
  all.reserve(packed.size() / 6);
  for (size_t i = 0; i + 5 < packed.size(); i += 6) {
    PendingInfection c;
    c.person_id = static_cast<PersonId>(packed[i]);
    c.venue_id = static_cast<VenueId>(packed[i + 1]);
    c.infector_id = static_cast<PersonId>(packed[i + 2]);
    c.venue_type_id = static_cast<uint8_t>(packed[i + 3]);
    c.infector_symptom_id = static_cast<uint8_t>(packed[i + 4]);
    c.transmission_mode_index = static_cast<uint8_t>(packed[i + 5]);
    c.infection_time = current_time;
    c.encounter_type_id = kDefaultEncounterTypeId;
    all.push_back(c);
  }
  if (all.empty()) return 0;

  std::sort(all.begin(), all.end(),
            [](const PendingInfection& a, const PendingInfection& b) {
              if (a.person_id != b.person_id) return a.person_id < b.person_id;
              return a.venue_id < b.venue_id;
            });

  int applied = 0;
  PersonId last = -1;
  for (const auto& c : all) {
    if (c.person_id == last) continue;  // a lower venue id already took them
    last = c.person_id;

    Person* p = world_.getPerson(c.person_id);
    if (!p) continue;            // someone else's resident; their rank applies
    if (p->infection) continue;  // already infected, e.g. seeded this slot
    if (disease_ == nullptr) continue;

    float severity_factor = 1.0f;
    if (auto* gu = world_.getGeoUnit(p->geo_unit_id))
      severity_factor = gu->severity_factor;

    std::string venue_type_name;
    if (c.venue_type_id < world_.venue_type_names.size())
      venue_type_name = world_.venue_type_names[c.venue_type_id];

    const uint64_t seed = mix_seed(base_seed_, c.person_id,
                                   static_cast<uint64_t>(current_time * 1000),
                                   static_cast<uint64_t>(c.venue_id));
    p->infection = std::make_unique<Infection>(
        disease_, current_time, p, static_cast<unsigned int>(seed), &world_,
        venue_type_name, c.venue_id, severity_factor, c.infector_symptom_id, "",
        "", c.transmission_mode_index);

    if (event_logger_ != nullptr)
      event_logger_->logInfection(
          c.person_id, c.infector_id, c.venue_id, current_time,
          kDefaultEncounterTypeId, c.infector_symptom_id,
          c.transmission_mode_index, InfectionSource::Person);
    if (active_infections != nullptr) active_infections->insert(c.person_id);
    applied++;
  }
  return applied;
}

std::pair<int, PersonId> InteractionManager::sampleInfectorFromAccumSources(
    std::vector<PartialPresenceAccumSource>& srcs, SplitMix64& rng) const {
  if (srcs.empty()) return {0, -1};

  // Sort for determinism, then cumulative-sample.
  std::sort(srcs.begin(), srcs.end(),
            [](const PartialPresenceAccumSource& a,
               const PartialPresenceAccumSource& b) {
              if (a.mode != b.mode) return a.mode < b.mode;
              return a.infector < b.infector;
            });
  std::vector<double> cum;
  cum.reserve(srcs.size());
  double acc = 0.0;
  for (const auto& s : srcs) {
    acc += s.weighted;
    cum.push_back(acc);
  }
  int sampled = (acc > 0.0) ? sampleFromCumulative(cum, rng) : 0;
  if (sampled < 0) sampled = 0;
  if (sampled >= static_cast<int>(srcs.size())) return {0, -1};
  return {srcs[sampled].mode, srcs[sampled].infector};
}

std::vector<PersonId> InteractionManager::orderSusceptibles(
    const std::unordered_map<PersonId, double>& susc_lambda) const {
  std::vector<PersonId> ordered;
  ordered.reserve(susc_lambda.size());
  for (const auto& kv : susc_lambda) ordered.push_back(kv.first);
  std::sort(ordered.begin(), ordered.end());
  return ordered;
}

double InteractionManager::computeMemberSusceptibility(
    const Person* person, const VisitorInfo* visitor,
    double current_time) const {
  if (person) {
    return person->getSusceptibility(current_time, disease_->getName());
  }
  if (visitor) return 1.0 - visitor->immunity_level;
  return 0.0;
}

bool InteractionManager::processOnePartialSusceptible(
    PersonId susc_id, const std::unordered_map<PersonId, double>& susc_lambda,
    std::unordered_map<PersonId, std::vector<PartialPresenceAccumSource>>&
        susc_sources,
    double current_time, Venue* venue, uint8_t venue_type_id,
    VenueId actual_venue_id,
    const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
    std::unordered_set<PersonId>* active_infections,
    std::vector<PendingInfection>* pending_infections, uint64_t time_bits,
    uint64_t venue_key) {
  auto lambda_it = susc_lambda.find(susc_id);
  if (lambda_it == susc_lambda.end()) return false;
  double lambda = lambda_it->second;
  if (!(lambda > 0.0)) return false;

  Person* susc_person = world_.getPerson(susc_id);
  const VisitorInfo* visitor = nullptr;
  if (!susc_person && visitor_data) {
    auto it = visitor_data->find(susc_id);
    if (it != visitor_data->end()) visitor = &it->second;
  }

  if (!susc_person && !visitor) return false;
  double susceptibility =
      computeMemberSusceptibility(susc_person, visitor, current_time);
  if (!(susceptibility > 0.0)) return false;

  double total_risk = lambda;
  if (simulation_config_.regional_risk.enabled && venue) {
    total_risk *= venue->transmission_factor;
  }
  double prob = 1.0 - std::exp(-total_risk * susceptibility);
  if (!(prob > 1e-12)) return false;

  SplitMix64 susc_rng(mix_seed(base_seed_, susc_id, venue_key, time_bits));
  double rng_roll = uniform_dist_(susc_rng);
  if (!(rng_roll < prob)) return false;

  // Source attribution: weight-sample from accumulated AccumSource entries.
  auto src_it = susc_sources.find(susc_id);
  int sampled_mode = 0;
  PersonId infector_id = -1;
  if (src_it != susc_sources.end()) {
    std::tie(sampled_mode, infector_id) =
        sampleInfectorFromAccumSources(src_it->second, susc_rng);
  }

  uint16_t infector_symptom_id =
      resolveInfectorSymptomId(infector_id, current_time, visitor_data);

  const uint8_t transmission_mode_index = static_cast<uint8_t>(sampled_mode);
  recordPartialPresenceCandidate(susc_id, infector_id, transmission_mode_index,
                                 infector_symptom_id, current_time,
                                 venue_type_id, actual_venue_id);
  (void)active_infections;
  (void)pending_infections;
  return true;
}

int InteractionManager::processPartialPresenceVenue(
    const std::vector<InteractionMember>& members, Venue* venue,
    VenueId actual_venue_id, double current_time, double delta_hours,
    std::unordered_set<PersonId>* active_infections,
    const std::unordered_set<PersonId>* /*visitor_ids*/,
    std::vector<PendingInfection>* pending_infections,
    const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
    uint8_t encounter_type_id,
    const CompartmentalModelManager* /*comp_model*/) {
  PartialPresenceLambdaResult acc = computePartialPresenceLambda(
      members, venue, actual_venue_id, current_time, delta_hours, visitor_data,
      encounter_type_id);
  auto& susc_lambda = acc.susc_lambda;
  auto& susc_sources = acc.susc_sources;

  const uint8_t venue_type_id = venue->type_id;

  // Per-susceptible Bernoulli draws (person_id-ordered for deterministic
  // per-call work).
  std::vector<PersonId> ordered_susc = orderSusceptibles(susc_lambda);

  int new_infections = 0;
  const uint64_t time_bits = static_cast<uint64_t>(current_time * 1000);
  const uint64_t venue_key = static_cast<uint64_t>(actual_venue_id);

  for (PersonId susc_id : ordered_susc) {
    if (processOnePartialSusceptible(
            susc_id, susc_lambda, susc_sources, current_time, venue,
            venue_type_id, actual_venue_id, visitor_data, active_infections,
            pending_infections, time_bits, venue_key)) {
      new_infections++;
    }
  }

  return new_infections;
}

}  // namespace june
