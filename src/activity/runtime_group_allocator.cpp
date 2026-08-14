#include "activity/runtime_group_allocator.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <unordered_map>

#include "core/config.h"
#include "core/world_state.h"
#include "utils/deterministic_rng.h"

#ifdef USE_MPI
#include <mpi.h>

#include "parallel/mpi_utils.h"
#endif

namespace june {

RuntimeGroupAllocator::RuntimeGroupAllocator(const WorldState& world,
                                             const Config& config)
    : world_(world), config_(config) {
  venue_type_mask_ =
      config_.simulation.partial_presence.enabled_venue_type_mask;
}

void RuntimeGroupAllocator::allocateForSlot(
    int time_slot_index, int day_type_idx, const TimeSlot& /*slot*/,
    double current_simulation_time, double delta_hours,
    const std::vector<PersonLocation>& locations) {
  // Clear group assignments from the previous slot.
  group_by_vid_pid_.clear();
  num_groups_by_venue_.clear();
  windows_by_vid_pid_.clear();
  f_presence_by_vid_pid_.clear();
  subset_by_vid_pid_.clear();
  venue_type_by_vid_.clear();
  riders_by_venue_.clear();
  legs_by_person_.clear();

  if (venue_type_mask_ == 0) return;  // feature off; zero overhead

  // Slot key for canonical RNG seeding across ranks. Sim time is in days;
  // round to integer minutes to avoid FP-summation drift between ranks.
  const uint64_t slot_minutes =
      static_cast<uint64_t>(current_simulation_time * 1440.0 + 0.5);
  const uint64_t base_seed = config_.simulation.random_seed;

  // Slot duration in minutes; input to presence_window helper.
  const float slot_duration_min = static_cast<float>(delta_hours * 60.0);

  // Membership-metadata field indices for per-leg boarding/alighting.
  // These are deterministic per world; looked up once per call (cheap).
  const int tb_field = world_.getMembershipFieldIndex("t_board_min");
  const int ta_field = world_.getMembershipFieldIndex("t_alight_min");
  const bool have_timing = (tb_field >= 0 && ta_field >= 0);

  // -------------------------------------------------------------------------
  // Phase 1: collect LOCAL (rider, leg) entries.
  //
  // For each person assigned a non-trivial activity this slot, walk every
  // activity_meta with that activity_index and add EACH leg whose venue is
  // partial-presence. A multi-leg journey (e.g. a 3-leg train commute)
  // produces 3 entries, one per leg, keyed independently by their flat
  // activity_venue index.
  //
  // We deliberately include legs whose raw (t_board_min, t_alight_min)
  // sit outside the slot window: a long-distance commuter (e.g. Durham →
  // London) must still expose riders on their destination-region legs to
  // model inter-regional disease spread. The FOI loop downstream computes
  // co-presence on these raw offsets (exact contacts) and caps each rider's
  // total commute presence at one slot-hour via the per-rider factor f_p.
  // -------------------------------------------------------------------------
  struct LocalLeg {
    uint32_t av_idx;
    VenueId venue_id;
    PersonId person_id;
    SubsetIndex subset;
    float eff_board;
    float eff_alight;
    float f_presence;
  };
  std::vector<LocalLeg> local_legs;
  local_legs.reserve(locations.size() / 64);

  // Per-rider scratch reused across iterations to avoid per-person
  // allocations. The full leg list for a single rider is small (commute
  // ≤ 4 legs), so these stay tiny.
  struct RawLeg {
    uint32_t av_idx;
    VenueId venue_id;
    SubsetIndex subset;
    float tb_min;
    float ta_min;
  };
  std::vector<RawLeg> raw_legs;
  std::vector<float> tb_buf, ta_buf;

  for (size_t i = 0; i < locations.size(); ++i) {
    const int16_t act = locations[i].activity_index;
    if (act < 0) continue;  // person inactive this slot

    const Person& person = world_.people[i];

    // Gather every partial-presence leg this rider has on this slot's
    // activity. We pass the full leg list (across all venues) to
    // computePresenceWindows so the per-rider presence cap f_p sees the
    // rider's total journey, not just one venue's piece. f_p is decided
    // here, on the rider's home rank, then broadcast with the windows.
    raw_legs.clear();
    for (int m = 0; m < person.activity_meta_count; ++m) {
      const auto& meta = world_.activity_meta[person.activity_meta_start + m];
      if (meta.activity_index != act) continue;

      auto venues = world_.getActivityVenues(meta);
      for (size_t k = 0; k < venues.size(); ++k) {
        const VenueId vid = venues[k].first;
        const uint8_t vt = world_.getVenueTypeId(vid);
        if (vt >= 64) continue;
        if (((venue_type_mask_ >> vt) & 1ULL) == 0) continue;

        const uint32_t av_idx = static_cast<uint32_t>(meta.venue_start + k);
        float tb = 0.0f, ta = slot_duration_min;
        if (have_timing) {
          tb = world_.getMembershipField(av_idx, tb_field);
          ta = world_.getMembershipField(av_idx, ta_field);
          if (tb == WorldState::kMembershipFieldAbsent ||
              ta == WorldState::kMembershipFieldAbsent) {
            // Missing per-leg metadata: degrade to full-slot presence so
            // the leg still participates in bucketing (D14: never drop
            // a leg).
            tb = 0.0f;
            ta = slot_duration_min;
          }
        }
        raw_legs.push_back({av_idx, vid, venues[k].second, tb, ta});
      }
    }
    if (raw_legs.empty()) continue;

    // No sort: windows are raw line-local offsets (each leg in its own venue's
    // clock) and f_p is order-independent, so leg order does not affect the
    // result. leg_idx remains the authoritative journey order for the future
    // multi-slot follow-up, but is not load-bearing here.
    tb_buf.resize(raw_legs.size());
    ta_buf.resize(raw_legs.size());
    for (size_t j = 0; j < raw_legs.size(); ++j) {
      tb_buf[j] = raw_legs[j].tb_min;
      ta_buf[j] = raw_legs[j].ta_min;
    }
    float f_presence = 1.0f;
    auto windows =
        computePresenceWindows(tb_buf.data(), ta_buf.data(), raw_legs.size(),
                               slot_duration_min, &f_presence);

    for (size_t j = 0; j < raw_legs.size(); ++j) {
      local_legs.push_back(LocalLeg{raw_legs[j].av_idx, raw_legs[j].venue_id,
                                    person.id, raw_legs[j].subset,
                                    windows[j].eff_board, windows[j].eff_alight,
                                    f_presence});
    }
  }

  // -------------------------------------------------------------------------
  // Phase 2: exchange (venue_id, person_id) pairs across ranks so every
  // rank has the same global rider list per partial-presence venue.
  // Multi-leg riders contribute one pair per leg (same person_id may
  // appear in multiple venues' lists; that's the point).
  // -------------------------------------------------------------------------
  // 7 int32s per leg: (vid, pid, subset, venue_type, bitcast(eff_board_f32),
  // bitcast(eff_alight_f32), bitcast(f_presence_f32)). Floats are bit-cast
  // through int32 so the existing int-typed Allgatherv ships them unchanged.
  //
  // The venue type travels with the leg because a rank only knows the types of
  // venues in its own halo. A line reachable only from another rank would
  // otherwise look untyped here, and the group deal below would quietly skip
  // it, leaving the rider tables disagreeing between ranks.
  // Every rank receives the SAME bytes the home rank packed → bit-identical
  // float values everywhere, no FP-recomputation drift. f_presence is a pure
  // function of the rider's own legs, so broadcasting it (rather than
  // recomputing per rank) keeps cross-rank visitors bit-identical to the home
  // rank.
  static_assert(
      sizeof(float) == sizeof(int32_t),
      "RuntimeGroupAllocator assumes 32-bit float for window packing");
  std::vector<int32_t> local_packed;
  local_packed.reserve(local_legs.size() * 7);
  for (const auto& lg : local_legs) {
    int32_t eb_bits, ea_bits, fp_bits;
    std::memcpy(&eb_bits, &lg.eff_board, sizeof(int32_t));
    std::memcpy(&ea_bits, &lg.eff_alight, sizeof(int32_t));
    std::memcpy(&fp_bits, &lg.f_presence, sizeof(int32_t));
    local_packed.push_back(static_cast<int32_t>(lg.venue_id));
    local_packed.push_back(static_cast<int32_t>(lg.person_id));
    local_packed.push_back(static_cast<int32_t>(lg.subset));
    local_packed.push_back(
        static_cast<int32_t>(world_.getVenueTypeId(lg.venue_id)));
    local_packed.push_back(eb_bits);
    local_packed.push_back(ea_bits);
    local_packed.push_back(fp_bits);
  }

  std::vector<int32_t> global_packed = local_packed;

#ifdef USE_MPI
  int world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size > 1) {
    int local_count = static_cast<int>(local_packed.size());
    std::vector<int> counts(world_size);
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    std::vector<int> displs;
    int total = 0;
    mpi_utils::computeDisplacements(counts, displs, total);
    global_packed.assign(total, 0);
    MPI_Allgatherv(local_packed.data(), local_count, MPI_INT,
                   global_packed.data(), counts.data(), displs.data(), MPI_INT,
                   MPI_COMM_WORLD);
  }
#endif

  // -------------------------------------------------------------------------
  // Phase 3: decode the global 5-int-per-leg stream, group by venue,
  // canonical-sort, deal into groups. Also populates windows_by_vid_pid_ and
  // f_presence_by_vid_pid_ from the broadcast float bytes (every rank sees
  // identical windows + f_p for the same rider, required for MPI determinism
  // on cross-LGU venues).
  // -------------------------------------------------------------------------
  std::unordered_map<VenueId, std::vector<PersonId>> riders_by_venue;
  windows_by_vid_pid_.reserve(global_packed.size() / 7);
  f_presence_by_vid_pid_.reserve(global_packed.size() / 7);
  subset_by_vid_pid_.reserve(global_packed.size() / 7);
  for (size_t i = 0; i + 6 < global_packed.size(); i += 7) {
    VenueId v = static_cast<VenueId>(global_packed[i]);
    PersonId p = static_cast<PersonId>(global_packed[i + 1]);
    SubsetIndex sub = static_cast<SubsetIndex>(global_packed[i + 2]);
    venue_type_by_vid_[v] = static_cast<uint8_t>(global_packed[i + 3]);
    float eb, ea, fp;
    std::memcpy(&eb, &global_packed[i + 4], sizeof(float));
    std::memcpy(&ea, &global_packed[i + 5], sizeof(float));
    std::memcpy(&fp, &global_packed[i + 6], sizeof(float));
    riders_by_venue[v].push_back(p);
    const uint64_t key =
        (static_cast<uint64_t>(v) << 32) | static_cast<uint64_t>(p);
    windows_by_vid_pid_[key] = EffectiveWindow{eb, ea};
    f_presence_by_vid_pid_[key] = fp;
    subset_by_vid_pid_[key] = sub;
  }

  // (venue_id, person_id) → group. The same person on multiple venues gets a
  // separate group per venue (multi-leg). For local riders we then look this
  // up by their leg's (vid, pid) to populate group_by_av_idx_.
  // group_by_vid_pid_ is a member so the FOI loop (which only has
  // InteractionMember = {pid, ...}) can look up groups without the flat av_idx.
  group_by_vid_pid_.reserve(global_packed.size() / 7);
  num_groups_by_venue_.reserve(riders_by_venue.size());

  for (auto& [vid, riders] : riders_by_venue) {
    const uint8_t vt = venue_type_by_vid_.at(vid);
    const int tgs = config_.simulation.partial_presence.getTargetGroupSize(vt);
    if (tgs <= 0)
      throw std::runtime_error(
          "RuntimeGroupAllocator: line " + std::to_string(vid) + " of type " +
          std::to_string(static_cast<int>(vt)) +
          " has riders but no target group size, so its groups cannot be "
          "sized");

    // Canonical shuffle: sort by hash(seed, slot_minutes, vid, pid). Same
    // hash on every rank → same order → same group assignments. PersonId as
    // a deterministic tiebreaker on the (vanishingly rare) hash collision.
    auto hash_key = [&](PersonId pid) {
      return mix_seed(
          mix_seed(base_seed, slot_minutes, static_cast<uint64_t>(vid)),
          static_cast<uint64_t>(pid));
    };
    std::sort(riders.begin(), riders.end(), [&](PersonId a, PersonId b) {
      uint64_t ha = hash_key(a), hb = hash_key(b);
      if (ha != hb) return ha < hb;
      return a < b;
    });

    const int num_groups =
        std::max(1, static_cast<int>((riders.size() + tgs - 1) / tgs));
    num_groups_by_venue_[vid] = static_cast<uint16_t>(num_groups);

    for (size_t pos = 0; pos < riders.size(); ++pos) {
      const PersonId pid = riders[pos];
      const uint16_t group = static_cast<uint16_t>(pos % num_groups);
      const uint64_t key =
          (static_cast<uint64_t>(vid) << 32) | static_cast<uint64_t>(pid);
      group_by_vid_pid_[key] = group;
    }
  }

  // -------------------------------------------------------------------------
  // Phase 4: turn the (venue, person) maps into the membership lists the FOI
  // walks.
  // -------------------------------------------------------------------------
  reindexRiders();

  // One-time visibility on first non-empty slot.
  static bool printed_first = false;
  if (!printed_first && !riders_by_venue.empty()) {
    size_t total_global_pairs = 0;
    size_t max_groups = 0;
    size_t multi_leg_local = 0;
    {
      std::unordered_map<PersonId, int> per_pid;
      for (const auto& lg : local_legs) per_pid[lg.person_id]++;
      for (auto& kv : per_pid)
        if (kv.second > 1) multi_leg_local++;
    }
    for (const auto& kv : riders_by_venue) {
      total_global_pairs += kv.second.size();
      const int tgs = config_.simulation.partial_presence.getTargetGroupSize(
          venue_type_by_vid_.at(kv.first));
      if (tgs > 0) {
        size_t nb = std::max<size_t>(1, (kv.second.size() + tgs - 1) / tgs);
        max_groups = std::max(max_groups, nb);
      }
    }
    std::cout << "[RuntimeGroupAllocator] first active slot: time_slot_index="
              << time_slot_index << " day_type=" << day_type_idx
              << " partial-presence venues=" << riders_by_venue.size()
              << " global_legs=" << total_global_pairs
              << " local_legs=" << local_legs.size()
              << " local_multi_leg_riders=" << multi_leg_local
              << " max_groups_in_a_venue=" << max_groups << std::endl;

    printed_first = true;
  }
}

uint8_t RuntimeGroupAllocator::venueTypeOf(VenueId venue_id) const {
  // A rank only knows the types of venues in its own halo, so for a line
  // somebody else is riding, trust the type that came over the wire with them.
  auto it = venue_type_by_vid_.find(venue_id);
  if (it != venue_type_by_vid_.end()) return it->second;
  return world_.getVenueTypeId(venue_id);
}

bool RuntimeGroupAllocator::isPartialPresenceVenue(VenueId venue_id) const {
  if (venue_id < 0) return false;
  const uint8_t vt = venueTypeOf(venue_id);
  return vt < 64 && ((venue_type_mask_ >> vt) & 1ULL) != 0;
}

void RuntimeGroupAllocator::reindexRiders() {
  riders_by_venue_.clear();
  legs_by_person_.clear();
  riders_by_venue_.reserve(num_groups_by_venue_.size());

  for (const auto& [key, group] : group_by_vid_pid_) {
    const VenueId vid = static_cast<VenueId>(key >> 32);
    const PersonId pid = static_cast<PersonId>(key & 0xFFFFFFFFull);
    auto sub_it = subset_by_vid_pid_.find(key);
    if (sub_it == subset_by_vid_pid_.end())
      throw std::runtime_error(
          "RuntimeGroupAllocator: rider " + std::to_string(pid) + " on venue " +
          std::to_string(vid) + " has a group but no subset");
    riders_by_venue_[vid].push_back(Rider{pid, sub_it->second, group});
    legs_by_person_[pid].push_back(vid);
  }

  // Hash-map iteration order varies between ranks, so impose one. Every rank
  // then walks each line's riders in the same order and accumulates its FOI in
  // the same sequence, which is what keeps the floating-point sums identical.
  for (auto& [vid, riders] : riders_by_venue_)
    std::sort(riders.begin(), riders.end(),
              [](const Rider& a, const Rider& b) { return a.pid < b.pid; });
  for (auto& [pid, legs] : legs_by_person_) std::sort(legs.begin(), legs.end());
}

void RuntimeGroupAllocator::attachFollowers(
    const std::vector<std::pair<PersonId, PersonId>>& pairs) {
  if (venue_type_mask_ == 0) return;

  std::vector<int32_t> packed;
  packed.reserve(pairs.size() * 2);
  for (const auto& [follower, host] : pairs) {
    packed.push_back(static_cast<int32_t>(follower));
    packed.push_back(static_cast<int32_t>(host));
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

  std::vector<std::pair<PersonId, PersonId>> all;
  all.reserve(packed.size() / 2);
  for (size_t i = 0; i + 1 < packed.size(); i += 2)
    all.emplace_back(static_cast<PersonId>(packed[i]),
                     static_cast<PersonId>(packed[i + 1]));
  if (all.empty()) return;

  // A follower binds to one host, so sorting by follower gives every rank the
  // same list with no duplicates left to worry about.
  std::sort(all.begin(), all.end());
  all.erase(std::unique(all.begin(), all.end()), all.end());

  for (const auto& [follower, host] : all) {
    const std::vector<VenueId> host_legs = legsOf(host);
    if (host_legs.empty())
      throw std::runtime_error(
          "RuntimeGroupAllocator::attachFollowers: host " +
          std::to_string(host) +
          " rides no partial-presence venue, so follower " +
          std::to_string(follower) +
          " cannot travel with it. Follow placed the follower on a line the "
          "allocator never saw.");

    // Give up the follower's own journey before taking the host's.
    for (VenueId own : legsOf(follower)) {
      const uint64_t key =
          (static_cast<uint64_t>(own) << 32) | static_cast<uint64_t>(follower);
      group_by_vid_pid_.erase(key);
      windows_by_vid_pid_.erase(key);
      f_presence_by_vid_pid_.erase(key);
      subset_by_vid_pid_.erase(key);
    }

    for (VenueId leg : host_legs) {
      const uint64_t hkey =
          (static_cast<uint64_t>(leg) << 32) | static_cast<uint64_t>(host);
      const uint64_t fkey =
          (static_cast<uint64_t>(leg) << 32) | static_cast<uint64_t>(follower);
      group_by_vid_pid_[fkey] = group_by_vid_pid_.at(hkey);
      windows_by_vid_pid_[fkey] = windows_by_vid_pid_.at(hkey);
      f_presence_by_vid_pid_[fkey] = f_presence_by_vid_pid_.at(hkey);
      subset_by_vid_pid_[fkey] = subset_by_vid_pid_.at(hkey);
    }
  }

  reindexRiders();
}

}  // namespace june
