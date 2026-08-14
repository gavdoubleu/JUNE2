#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "activity/coordinated_encounter_types.h"
#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/policy.h"

// Pass 1 of coordinated-encounter injection: which local participants of each
// encounter survive the alive + policy checks. Exposed here (rather than kept
// private to simulator_encounters.cpp) so the encounter tests drive the real
// production logic instead of a replica.
namespace june {
namespace encounters {

// Per-slot lookup tables derived from the coordinated-encounters config:
// encounter_type_id -> trigger activity indices, and
// encounter_type_id -> min_attendees threshold.
struct EncounterLookups {
  std::unordered_map<uint8_t, std::vector<int16_t>> trigger_activities;
  std::unordered_map<uint8_t, int> min_attendees;
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
    const WorldState& world,
    const std::vector<CoordinatedEncounterDef>& encounters);

// Per-encounter, compute the local participants who survive the alive +
// policy-block checks. Encounters not scheduled for this slot are skipped.
// Each returned entry points back at its source encounter via .encounter_idx.
//
// The policy question asked is "would a policy move this person away from the
// encounter's venue at this slot", so the venue handed to the policy is
// enc.venue_id — the venue injection is about to put them in — not their
// currently scheduled one.
std::vector<EncounterEligibility> computeLocalEligibility(
    const std::vector<CoordinatedEncounter>& daily_encounters,
    int time_slot_index, double current_simulation_time,
    const EncounterLookups& lookups, const WorldState& world,
    const std::vector<PersonLocation>& locations,
    PolicyManager* policy_manager);

}  // namespace encounters
}  // namespace june
