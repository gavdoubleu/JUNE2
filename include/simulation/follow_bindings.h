#pragma once

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/policy.h"

// Internal helpers of the follow subsystem, exposed here only so the follow
// unit tests can drive the binding logic directly on a synthetic world. The
// definitions live in simulator_encounters.cpp; nothing outside the simulator
// and its tests should depend on them.
namespace june {
namespace follow_detail {

// The host's venue of the given pool type, or -1 if it has none.
VenueId findPoolVenue(const WorldState& world, const Person& host,
                      int pool_venue_type_id);

// The host's candidate pool: co-members of its venue of the configured type, or
// its partners in the configured network. Returns global PersonIds.
std::vector<PersonId> gatherPool(const WorldState& world,
                                 const FollowConfig& fc, const Person& host);

// A host's deterministic per-(seed, host, day) probability roll.
bool hostRollsFollow(const FollowConfig& fc, uint64_t seed, PersonId host,
                     int hop_start_day);

// Rebuild the criteria bindings from scratch: every eligible person follows the
// lowest-id host in its pool. follower_excl / host_excl are the people already
// bound by an earlier rule (a follower is barred if in either; a host only if
// in host_excl).
void rebuildCriteriaBindings(
    WorldState& world, const FollowConfig& fc,
    std::unordered_map<PersonId, PersonId>& follower_host,
    std::unordered_set<PersonId>& active_hosts,
    std::vector<std::pair<PersonId, PersonId>>* remote_hosts,
    std::unordered_map<PersonId, PersonId>* new_follows,
    const std::unordered_set<PersonId>& follower_excl,
    const std::unordered_set<PersonId>& host_excl);

// Does an exception stop this follower from being mirrored onto its host for
// the slot? Any of the three firing means the follower keeps its own schedule:
// the host is at an excepted activity, the host is at an excepted venue type,
// or the follower has an activity of its own that outranks following. Pass -1
// for an activity the person does not have this slot.
bool mirrorSuppressed(const FollowConfig& fc, int16_t host_activity,
                      uint8_t host_venue_type, int16_t follower_activity);

// Does the follower's own policy outrank the mirror? The venue put to the
// policy is the host's — the venue the mirror is about to move the follower
// to — not the follower's own scheduled one, so a venue-gated policy sees
// where the follower would end up. The activity and subset stay the
// follower's, since the mirror leaves those alone. True means a policy fires
// and the follower stays where the policy put them.
bool policySuppressesMirror(PolicyManager* policy_manager, Person& follower,
                            const PersonLocation& follower_location,
                            VenueId host_venue, double current_time,
                            int time_slot_index);

// Stochastic enrolment: each host rolls once and gathers followers from its
// pool. Returns {hosts that gained followers, total local followers enrolled}.
std::pair<int, int> enrolFollowHosts(
    WorldState& world, const FollowConfig& fc, const ScheduleConfig& sched,
    uint64_t seed, bool standing, std::unordered_set<PersonId>& active_hosts,
    std::unordered_map<PersonId, PersonId>& follower_host, int current_day,
    std::vector<std::pair<PersonId, PersonId>>* remote_invites,
    std::unordered_map<PersonId, PersonId>* new_follows,
    const std::unordered_set<PersonId>& follower_excl,
    const std::unordered_set<PersonId>& host_excl);

}  // namespace follow_detail
}  // namespace june
