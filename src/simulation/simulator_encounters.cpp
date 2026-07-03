// Simulator coordinated-encounter helpers: daily negotiation + slot injection
// + local-rank logging. Split from simulator.cpp (declared in
// simulation/simulator.h).
#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "simulation/simulator.h"

namespace june {

namespace {

// Per-def hop-trigger info: schedule to hop non-host participants into on
// injection, the schedule to return to afterwards, and the trigger
// activity's own index (stamped onto activity_index so it doesn't go
// stale once the trigger activity is coordinated_only). hop_schedule_idx
// == -1 means this def has no hop (ordinary in-place encounter).
struct EncounterHopInfo {
  int16_t hop_schedule_idx = -1;
  int16_t return_schedule_idx = -1;
  int16_t trigger_activity_idx = -1;
};

// Per-slot lookup tables derived from the coordinated-encounters config:
// encounter_type_id -> trigger activity indices, and
// encounter_type_id -> min_attendees threshold.
struct EncounterLookups {
  std::unordered_map<uint8_t, std::vector<int16_t>> trigger_activities;
  std::unordered_map<uint8_t, int> min_attendees;
  std::unordered_map<uint8_t, EncounterHopInfo> hop_info;
};

// Pass-1 result for one daily_encounter: which local participants pass the
// eligibility checks (alive + not policy-blocked at this slot), how many of
// them, and the threshold needed before this encounter gets injected.
struct EncounterEligibility {
  int encounter_idx;                     // index into daily_encounters
  std::vector<size_t> eligible_indices;  // local people passing policy
  int local_eligible;
  int min_required;
};

EncounterLookups buildEncounterLookups(
    const WorldState& world, const ScheduleConfig& schedule_config,
    const std::vector<CoordinatedEncounterDef>& encounters) {
  EncounterLookups out;
  for (const auto& def : encounters) {
    int type_id = world.getEncounterTypeIndex(def.name);
    if (type_id < 0) continue;
    std::vector<int16_t> indices;
    for (const auto& slot_name : def.trigger_slots) {
      int idx = world.getActivityIndex(slot_name);
      if (idx >= 0) indices.push_back(static_cast<int16_t>(idx));
    }
    out.trigger_activities[static_cast<uint8_t>(type_id)] = std::move(indices);
    out.min_attendees[static_cast<uint8_t>(type_id)] = def.min_attendees;

    EncounterHopInfo hop_info;
    hop_info.trigger_activity_idx = def.cached_trigger_activity_idx;
    if (def.cached_hop_schedule_idx >= 0 &&
        static_cast<size_t>(def.cached_hop_schedule_idx) <
            schedule_config.schedule_types.size()) {
      const ScheduleType& target =
          schedule_config.schedule_types[def.cached_hop_schedule_idx];
      // Only temporary schedules with flat_slots make sense as a hop
      // target here (mirrors the guard in
      // ActivityManager::maybeTriggerScheduleHop) — a permanent hop_schedule
      // is not supported by the consumeSlot0()-based onset in
      // applyEncounterInjection.
      if (target.is_temporary && !target.flat_slots.empty()) {
        hop_info.hop_schedule_idx = def.cached_hop_schedule_idx;
        hop_info.return_schedule_idx = target.return_schedule_idx;
      }
    }
    out.hop_info[static_cast<uint8_t>(type_id)] = hop_info;
  }
  return out;
}

// When both A→B and B→A proposals are accepted in the same slot, two
// encounters exist for the same pair. Collapse them by keeping the lowest
// encounter_id; then sort by encounter_id so injection order is
// deterministic across ranks (later encounters for the same person skip
// already-assigned slots).
std::vector<CoordinatedEncounter> dedupAndSortDailyEncounters(
    std::vector<CoordinatedEncounter> daily_encounters) {
  std::map<std::pair<std::set<PersonId>, int>, size_t> best;
  for (size_t i = 0; i < daily_encounters.size(); ++i) {
    auto key = std::make_pair(daily_encounters[i].participants,
                              daily_encounters[i].slot);
    auto it = best.find(key);
    if (it == best.end()) {
      best[key] = i;
    } else if (daily_encounters[i].encounter_id <
               daily_encounters[it->second].encounter_id) {
      it->second = i;
    }
  }
  std::vector<CoordinatedEncounter> deduped;
  deduped.reserve(best.size());
  for (auto& [key, idx] : best) {
    deduped.push_back(std::move(daily_encounters[idx]));
  }
  std::sort(deduped.begin(), deduped.end(),
            [](const CoordinatedEncounter& a, const CoordinatedEncounter& b) {
              return a.encounter_id < b.encounter_id;
            });
  return deduped;
}

// Pass 1: per-encounter, compute the local participants who survive the
// alive + policy-block checks. Encounters not scheduled for this slot are
// skipped. Returned vector is index-aligned with the surviving subset of
// daily_encounters (each entry stores back into the source via
// .encounter_idx).
std::vector<EncounterEligibility> computeLocalEligibility(
    const std::vector<CoordinatedEncounter>& daily_encounters,
    int time_slot_index, double current_simulation_time,
    const std::unordered_map<uint8_t, std::vector<int16_t>>&
        encounter_trigger_activities,
    const std::unordered_map<uint8_t, int>& encounter_min_attendees,
    const WorldState& world, const std::vector<PersonLocation>& locations,
    PolicyManager* policy_manager) {
  std::vector<EncounterEligibility> slot_encounters;
  for (int ei = 0; ei < static_cast<int>(daily_encounters.size()); ++ei) {
    const auto& enc = daily_encounters[ei];
    if (enc.slot != time_slot_index) continue;

    auto trig_it = encounter_trigger_activities.find(enc.encounter_type_id);

    std::vector<size_t> eligible_indices;
    for (PersonId pid : enc.participants) {
      auto it = world.person_index.find(pid);
      if (it == world.person_index.end()) continue;

      size_t array_idx = it->second;
      if (array_idx >= locations.size()) continue;

      const Person& person = world.people[array_idx];
      if (person.is_dead) continue;

      bool policy_blocked = false;
      if (policy_manager && trig_it != encounter_trigger_activities.end()) {
        for (int16_t trigger_act_idx : trig_it->second) {
          auto override = policy_manager->getOverride(
              const_cast<Person&>(world.people[array_idx]), trigger_act_idx,
              locations[array_idx].venue_id, locations[array_idx].subset_index,
              current_simulation_time, time_slot_index);
          if (override.has_value()) {
            policy_blocked = true;
            break;
          }
        }
      }
      if (!policy_blocked) eligible_indices.push_back(array_idx);
    }

    int min_required = 2;
    auto min_it = encounter_min_attendees.find(enc.encounter_type_id);
    if (min_it != encounter_min_attendees.end()) {
      min_required = min_it->second;
    }

    EncounterEligibility ee;
    ee.encounter_idx = ei;
    ee.eligible_indices = std::move(eligible_indices);
    ee.local_eligible = static_cast<int>(ee.eligible_indices.size());
    ee.min_required = min_required;
    slot_encounters.push_back(std::move(ee));
  }
  return slot_encounters;
}

// MPI eligibility exchange. Each rank only knows its local participants'
// eligibility; for encounters with any remote participants we Allgatherv
// (encounter_id, local_eligible) pairs and sum across ranks. The result
// maps encounter_id -> global eligible count, with entries only for
// encounters that had any remote participants. Serial or single-rank
// builds return an empty map.
std::unordered_map<int, int> exchangeGlobalEligibility(
    const std::vector<EncounterEligibility>& slot_encounters,
    const std::vector<CoordinatedEncounter>& daily_encounters,
    const WorldState& world, DomainManager* domain_mgr) {
  std::unordered_map<int, int> global_eligible_map;
#ifdef USE_MPI
  if (!domain_mgr) return global_eligible_map;
  // Collect (encounter_id, local_eligible) for encounters with any remote
  // participants; these are the only ones that need exchange.
  std::vector<int> local_pairs;  // flat array: [eid, count, eid, count, ...]
  for (const auto& ee : slot_encounters) {
    const auto& enc = daily_encounters[ee.encounter_idx];
    bool has_remote = false;
    for (PersonId pid : enc.participants) {
      if (world.person_index.find(pid) == world.person_index.end()) {
        has_remote = true;
        break;
      }
    }
    if (has_remote) {
      local_pairs.push_back(enc.encounter_id);
      local_pairs.push_back(ee.local_eligible);
    }
  }

  int local_count = static_cast<int>(local_pairs.size());
  int num_ranks = domain_mgr->getNumRanks();
  std::vector<int> all_counts(num_ranks);
  MPI_Allgather(&local_count, 1, MPI_INT, all_counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD);

  std::vector<int> displs(num_ranks, 0);
  int total = 0;
  for (int r = 0; r < num_ranks; ++r) {
    displs[r] = total;
    total += all_counts[r];
  }

  std::vector<int> all_pairs(total);
  MPI_Allgatherv(local_pairs.data(), local_count, MPI_INT, all_pairs.data(),
                 all_counts.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);

  for (int i = 0; i < total; i += 2) {
    global_eligible_map[all_pairs[i]] += all_pairs[i + 1];
  }
#else
  (void)slot_encounters;
  (void)daily_encounters;
  (void)world;
  (void)domain_mgr;
#endif
  return global_eligible_map;
}

// Pass 2: stamp the per-encounter venue_id / encounter_type_id onto every
// local eligible participant of any encounter that meets its min_required
// threshold. Then for virtual venues (negative IDs) register the host-rank
// ownership so the visitor exchange routes cross-rank participants
// correctly.
void applyEncounterInjection(
    const std::vector<EncounterEligibility>& slot_encounters,
    const std::vector<CoordinatedEncounter>& daily_encounters,
    const std::unordered_map<int, int>& global_eligible_map,
    const std::unordered_map<uint8_t, EncounterHopInfo>& hop_info,
    std::vector<PersonLocation>& locations, WorldState& world,
    DomainManager* domain_mgr) {
  for (size_t i = 0; i < slot_encounters.size(); ++i) {
    const auto& ee = slot_encounters[i];
    const auto& enc = daily_encounters[ee.encounter_idx];

    // Use global eligible count if available (MPI mode with remote
    // participants); otherwise use local count (serial or all-local).
    int total_eligible;
    auto ge_it = global_eligible_map.find(enc.encounter_id);
    if (ge_it != global_eligible_map.end()) {
      total_eligible = ge_it->second;
    } else {
      total_eligible = ee.local_eligible;
    }

    if (total_eligible < ee.min_required || ee.local_eligible <= 0) continue;

    auto hop_it = hop_info.find(enc.encounter_type_id);
    const EncounterHopInfo* hop =
        (hop_it != hop_info.end()) ? &hop_it->second : nullptr;

    for (size_t array_idx : ee.eligible_indices) {
      locations[array_idx].venue_id = enc.venue_id;
      locations[array_idx].encounter_type_id = enc.encounter_type_id;
      if (hop && hop->trigger_activity_idx >= 0) {
        locations[array_idx].activity_index = hop->trigger_activity_idx;
      }

      // Non-host participants of a hop-bearing encounter trip into the
      // configured schedule (e.g. day_trip/long_trip). selectVenue always
      // resolves an encounter to the host's own venue, so the host never
      // leaves home and is skipped here. The hop is started already
      // "consumed" (consumeSlot0) because this slot's venue has already
      // been stamped above, standing in for the hop target's flat_slots[0]
      // — the next timeslot resolves flat_slots[1] onward normally via
      // advanceHoppedSchedule.
      if (hop && hop->hop_schedule_idx >= 0) {
        Person& person = world.people[array_idx];
        if (person.id != enc.host_id) {
          int16_t effective_return =
              (hop->return_schedule_idx != -1)
                  ? hop->return_schedule_idx
                  : static_cast<int16_t>(person.schedule_type_id);
          person.schedule_hop =
              ScheduleHop::begin(hop->hop_schedule_idx, effective_return);
          person.schedule_hop.consumeSlot0();
        }
      }
    }

#ifdef USE_MPI
    // Register virtual venue ownership so the visitor exchange can route
    // cross-rank encounter participants to the host's rank. Physical
    // venues already have ownership via the venue ownership map; virtual
    // venues (negative IDs) need explicit registration.
    if (domain_mgr && enc.venue_id < 0) {
      int host_rank = domain_mgr->getPersonRank(enc.host_id);
      domain_mgr->setVenueRank(enc.venue_id, host_rank);
      if (host_rank == domain_mgr->getRank()) {
        domain_mgr->getDomain().registerVirtualVenue(enc.venue_id);
      }
    }
#else
    (void)domain_mgr;
#endif
  }
}

}  // namespace

void Simulator::negotiateAndLogDailyEncounters(int day, int rank) {
  if (!coordinated_encounter_manager_) return;
  coordinated_encounter_manager_->resetDaily();

  int encounter_day_type_idx = config_.schedule.getDayTypeIndex(day);

  // Phase 1: Generate proposals locally
  std::vector<EncounterProposal> local_proposals;
  coordinated_encounter_manager_->generateProposals(day, local_proposals,
                                                    encounter_day_type_idx);

  // Accumulate proposal stats for debug summary
  coordinated_encounter_manager_->accumulateProposalStats(local_proposals);

  // Phase 2: Exchange proposals across ranks, then process.
  // Keep local_proposals alive for mutual proposal detection.
  std::vector<EncounterProposal> all_proposals;
#ifdef USE_MPI
  if (domain_mgr_) {
    domain_mgr_->exchangeEncounterProposals(local_proposals, all_proposals);
  } else {
    all_proposals = local_proposals;  // copy, not move
  }
#else
  all_proposals = local_proposals;  // copy, not move
#endif

  std::vector<EncounterReply> local_replies;
  coordinated_encounter_manager_->processProposals(
      all_proposals, local_proposals, local_replies, encounter_day_type_idx);

  // Phase 3: Exchange replies back to host ranks, then finalize
  std::vector<EncounterReply> all_replies;
#ifdef USE_MPI
  if (domain_mgr_) {
    domain_mgr_->exchangeEncounterReplies(local_replies, all_replies);
  } else {
    all_replies = std::move(local_replies);
  }
#else
  all_replies = std::move(local_replies);
#endif

  // Accumulate reply stats for debug summary
  coordinated_encounter_manager_->accumulateReplyStats(all_replies);

  std::vector<CoordinatedEncounter> finalized;
  coordinated_encounter_manager_->finalizeEncounters(all_replies, finalized);

  // Accumulate finalize stats for debug summary
  coordinated_encounter_manager_->accumulateFinalizeStats(finalized);

  // Log locally-finalized encounters to HDF5 BEFORE merging remote ones, so
  // cross-rank encounters never get double-logged.
  logFinalizedEncountersLocally(finalized, rank);

  // Phase 4: Exchange finalized encounters so remote participants know
#ifdef USE_MPI
  if (domain_mgr_) {
    std::vector<CoordinatedEncounter> remote_finalized;
    domain_mgr_->exchangeFinalizedEncounters(finalized, remote_finalized);
    // Add remote encounters that involve our local people
    for (auto& enc : remote_finalized) {
      coordinated_encounter_manager_->addDailyEncounter(enc);
    }
  }
#endif
}

void Simulator::logFinalizedEncountersLocally(
    const std::vector<CoordinatedEncounter>& finalized, int rank) {
  // group_id fans one unique uint64 per real encounter across every pair-row
  // belonging to it. Rank is packed into the high 16 bits so counters stay
  // unique across MPI ranks without coordination; the low 48 bits are a
  // per-rank monotonic counter (2^48 events / rank is effectively unbounded
  // for realistic simulations).
  const uint64_t rank_prefix = static_cast<uint64_t>(rank) << 48;
  for (const auto& enc : finalized) {
    const uint64_t group_id =
        rank_prefix | (next_encounter_group_id_++ & 0x0000FFFFFFFFFFFFULL);
    for (PersonId pid : enc.participants) {
      if (pid != enc.host_id) {
        event_logger_.logCoordinatedEncounter(
            enc.host_id, pid, current_simulation_time_, enc.encounter_type_id,
            enc.slot, group_id);
      }
    }
  }
}

void Simulator::injectCoordinatedEncountersIntoSlot(int time_slot_index) {
  if (!coordinated_encounter_manager_) return;

  // Per-slot lookup tables (trigger activities + min attendees) keyed by
  // encounter_type_id. Cheap to rebuild per slot.
  EncounterLookups lookups = buildEncounterLookups(
      world_, config_.schedule, config_.coordinated_encounters.encounters);

  auto daily_encounters = dedupAndSortDailyEncounters(
      coordinated_encounter_manager_->getDailyEncounters());

  // Pass 1: compute each encounter's local eligible-participant set under
  // the alive + policy-block checks.
  std::vector<EncounterEligibility> slot_encounters = computeLocalEligibility(
      daily_encounters, time_slot_index, current_simulation_time_,
      lookups.trigger_activities, lookups.min_attendees, world_, locations_,
      policy_manager_.get());

  // Pass 1.5: in MPI mode, sum local_eligible across ranks for encounters
  // that span ranks so every rank sees the true global count.
  std::unordered_map<int, int> global_eligible_map = exchangeGlobalEligibility(
      slot_encounters, daily_encounters, world_, domain_mgr_);

  // Pass 2: stamp venue_id / encounter_type_id onto eligible participants
  // for encounters that meet their min_required threshold.
  applyEncounterInjection(slot_encounters, daily_encounters,
                          global_eligible_map, lookups.hop_info, locations_,
                          world_, domain_mgr_);
}

}  // namespace june
