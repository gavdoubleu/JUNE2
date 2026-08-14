#pragma once

#include <climits>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "core/types.h"
#include "presence_window.h"

namespace june {

class Config;
class WorldState;
struct TimeSlot;

// Per-slot allocator that deals riders of partial-presence venues into
// ephemeral runtime groups. A whole line is too big to mix as one pool, so its
// riders are split into disjoint groups that can actually contact each other:
// a train's carriages, a ferry's decks. Gated on
// SimulationConfig::partial_presence; an empty config (default) makes
// allocateForSlot a one-test no-op so non-commute scenarios pay zero cost.
//
// These groups are not contact-matrix bins. A bin says what kind of person
// someone is; a group says which pocket of the venue they are physically in.
// A rider has both.
//
// Multi-presence: a rider with N legs of a single activity (e.g. a 3-leg
// train commute) is placed in N groups, one per leg, and is a member of all N
// lines. The rider table this builds IS a line's membership: a person holds a
// single PersonLocation, which cannot name four lines at once, so the FOI reads
// its passengers from here instead.
//
// Dealing is config-driven and capacity-free: per slot, the number of groups
// for a venue emerges as max(1, ceil(N_global_riders / target_group_size)),
// where target_group_size is a per-venue-type config knob (roughly vehicle
// occupancy). Assignment is a round-robin deal over a canonical hash-sort of
// the GLOBAL rider list, so groups differ in size by at most 1 and every rank
// computes bit-identical assignments for the same rider without exchanging
// indices.
//
// MPI: one Allgatherv per slot collects every (rider, leg) across ranks, along
// with its subset, venue type, window and presence factor. The venue type has
// to travel with the leg because a rank only holds types for venues in its own
// halo. Dealing is then a pure function of (seed, slot_minutes, venue_id,
// person_id, num_groups) on every rank.
class RuntimeGroupAllocator {
 public:
  static constexpr uint16_t kNoGroup = UINT16_MAX;

  // One person on one line for one slot. A commuter with four legs is four
  // Riders, one per line. Lines take their membership from these rather than
  // from the location table, because a person holds a single location but can
  // ride several lines at once.
  struct Rider {
    PersonId pid;
    SubsetIndex subset;
    uint16_t group;
  };

  RuntimeGroupAllocator(const WorldState& world, const Config& config);

  // Called once per slot, after ActivityManager::assignActivitiesFromSchedule
  // returns. Cheap no-op when partial_presence is disabled.
  //   slot:                       current TimeSlot.
  //   current_simulation_time:    days since simulation start (float).
  //   delta_hours:                slot duration in hours (for FOI clamp,
  //                                consumed downstream).
  //   locations:                  per-person slot assignments just populated.
  //
  // The allocator includes EVERY leg of every multi-leg rider on partial-
  // presence venues, even legs whose raw (t_board_min, t_alight_min) sit
  // outside [0, slot_duration_min). This is intentional: long-distance
  // commuters (e.g. Durham → London, ~3 hr) must still expose riders on
  // their destination-region legs to capture geographic disease spread.
  // Co-presence is computed on the raw line-local offsets (exact contacts);
  // each rider's total commute presence is capped at one slot-hour by the
  // per-rider presence factor f_p (see getPresenceFactor), applied downstream
  // in the FOI loop so the day budget stays bounded.
  void allocateForSlot(int time_slot_index, int day_type_idx,
                       const TimeSlot& slot, double current_simulation_time,
                       double delta_hours,
                       const std::vector<PersonLocation>& locations);

  // Runtime group index for a specific (venue, person). Returns kNoGroup if the
  // pair wasn't bucketed this slot.
  uint16_t getGroupIndex(VenueId venue_id, PersonId person_id) const {
    const uint64_t key = (static_cast<uint64_t>(venue_id) << 32) |
                         static_cast<uint64_t>(person_id);
    auto it = group_by_vid_pid_.find(key);
    return it == group_by_vid_pid_.end() ? kNoGroup : it->second;
  }

  // Number of bins allocated for a venue this slot. Returns 0 if the venue
  // isn't partial-presence or has no riders this slot. Used by the FOI loop
  // to size per-group scratch buffers.
  uint16_t getNumGroups(VenueId venue_id) const {
    auto it = num_groups_by_venue_.find(venue_id);
    return it == num_groups_by_venue_.end() ? 0 : it->second;
  }

  // Per-(venue, person) effective presence window for this slot, in minutes
  // since slot start. Computed on each rider's home rank (so the proportional
  // policy for long-distance commuters sees the rider's full leg list) and
  // broadcast globally via the same Allgatherv that distributes bin
  // assignments, so cross-rank visitors get the SAME window the home rank
  // computed  (required for MPI-deterministic FOI on cross-LGU venues).
  //
  // Returns {0, 0} if the (venue, person) wasn't bucketed this slot.
  EffectiveWindow getPresenceWindow(VenueId venue_id,
                                    PersonId person_id) const {
    const uint64_t key = (static_cast<uint64_t>(venue_id) << 32) |
                         static_cast<uint64_t>(person_id);
    auto it = windows_by_vid_pid_.find(key);
    return it == windows_by_vid_pid_.end() ? EffectiveWindow{0.0f, 0.0f}
                                           : it->second;
  }

  // Per-rider presence cap f_p = min(1, slot / T_p) for this slot, where T_p is
  // the rider's total transport minutes. f_p = 1 for journeys that fit the slot
  // (~99% of riders) and < 1 for over-long journeys, scaling that rider's FOI
  // contribution (as source AND sink) so its total modeled commute presence
  // stays ≤ one slot-hour. A pure function of the rider's own legs, computed on
  // the home rank and broadcast with the windows, so it is identical on every
  // rank (MPI bit-identity). Returns 1.0 if the (venue, person) wasn't bucketed
  // this slot.
  float getPresenceFactor(VenueId venue_id, PersonId person_id) const {
    const uint64_t key = (static_cast<uint64_t>(venue_id) << 32) |
                         static_cast<uint64_t>(person_id);
    auto it = f_presence_by_vid_pid_.find(key);
    return it == f_presence_by_vid_pid_.end() ? 1.0f : it->second;
  }

  // Who rides each line this slot. The same on every rank, and it is what the
  // FOI walks when it builds a line's groups.
  const std::unordered_map<VenueId, std::vector<Rider>>& ridersByVenue() const {
    return riders_by_venue_;
  }

  // The lines one person rides this slot, in ascending venue order. Empty for
  // someone who rides none.
  const std::vector<VenueId>& legsOf(PersonId pid) const {
    static const std::vector<VenueId> kNone;
    auto it = legs_by_person_.find(pid);
    return it == legs_by_person_.end() ? kNone : it->second;
  }

  // Seat each follower beside its host. The follower takes over the host's
  // rider entries, so it rides every leg of the host's journey, in the host's
  // group, with the host's window and presence factor. Legs of its own are
  // given up, since it travels with the host instead of on its own route.
  //
  // The pairs are exchanged across ranks before anything is applied, so every
  // rank works from the same list in the same order and the tables stay
  // bit-identical. A host that rides no line is an error: follow would have
  // placed someone on a venue the allocator never saw.
  void attachFollowers(const std::vector<std::pair<PersonId, PersonId>>& pairs);

  // True iff the SimulationConfig declares any partial-presence venue types
  // that are actually present in this world.
  bool isActive() const { return venue_type_mask_ != 0; }

  // True iff the venue's type is declared partial-presence.
  bool isPartialPresenceVenue(VenueId venue_id) const;

  // The venue's type, taken from the rider broadcast when someone rides it and
  // from the world otherwise. Ranks only hold types for venues in their halo,
  // so asking the world directly about another rank's line gives the wrong
  // answer.
  uint8_t venueTypeOf(VenueId venue_id) const;

 private:
  // Rebuilds legs_by_person_ and the per-venue rider lists from the
  // (venue, person) maps, keeping both in a canonical order.
  void reindexRiders();
  const WorldState& world_;
  const Config& config_;

  // Bit i set iff venue type id i is a partial-presence type. Copied from
  // SimulationConfig::partial_presence::enabled_venue_type_mask at ctor.
  uint64_t venue_type_mask_ = 0;

  // Global (venue, person) → bin map, identical on every rank.
  // Key layout: (uint64_t(venue_id) << 32) | uint64_t(person_id).
  std::unordered_map<uint64_t, uint16_t> group_by_vid_pid_;

  // Per-venue bin count for this slot. Used by the FOI loop to size
  // per-group scratch buffers without re-deriving from rider counts.
  std::unordered_map<VenueId, uint16_t> num_groups_by_venue_;

  // Global (venue, person) → presence window. Same broadcast as
  // group_by_vid_pid_; the FOI loop consults this for both local and
  // cross-rank visitor riders so windows are identical on every rank.
  std::unordered_map<uint64_t, EffectiveWindow> windows_by_vid_pid_;

  // Global (venue, person) → per-rider presence cap f_p. Broadcast alongside
  // the windows in the same Allgatherv so it is bit-identical on every rank.
  // The same f_p is replicated for each of a multi-leg rider's (venue, person)
  // keys (f_p is a rider-level property, not a per-venue one).
  std::unordered_map<uint64_t, float> f_presence_by_vid_pid_;

  // Global (venue, person) → the subset the rider occupies on that line, so a
  // line's members can be binned into its contact matrix without consulting
  // the location table.
  std::unordered_map<uint64_t, SubsetIndex> subset_by_vid_pid_;

  // The membership view of the maps above: who is on each line, and which
  // lines each person is on. Derived, and global like the rest.
  std::unordered_map<VenueId, std::vector<Rider>> riders_by_venue_;
  std::unordered_map<PersonId, std::vector<VenueId>> legs_by_person_;

  // Venue type of every line someone rides, as broadcast with its riders. The
  // world's own type map only covers this rank's halo.
  std::unordered_map<VenueId, uint8_t> venue_type_by_vid_;
};

}  // namespace june
