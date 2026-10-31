#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "core/config.h"
#include "core/world_state.h"
#include "utils/deterministic_rng.h"
#include "utils/mpi_utils.h"

namespace june {

// Who is in each tier of a vaccination campaign and who would accept. Both
// depend only on attributes, ids and the run seed, so they are settled once at
// startup and never rolled again. Each person has two bits per tier: in the
// tier's group, and accepts.
//
// settleUptakeLocal works on one rank's world and needs no MPI. settleUptake
// is the collective one the simulator calls: it adds up the counts, throws one
// error with the same text on every rank, and prints a line per tier.

// Two bits per person per tier live in one 64 bit mask each.
// ponytail: a sparse list of (tier, accepted) per person has no limit, if a
// campaign ever needs more tiers.
constexpr size_t kMaxUptakeTiers = 64;

// Tags the uptake draw so it is never the same number as lotteryKey, which
// mixes the same seed, id, campaign and tier. The low 20 bits are left free for
// the tier and the person-or-unit flag.
constexpr uint64_t kUptakeDrawTag = 0xB7A4E5ULL << 20;

// The seed of the stream that settles whether one person, or one shared unit
// (keyed by its venue id), accepts. A person and a venue can share a number,
// so the flag keeps their streams apart, and the tier is in the key so no two
// tiers share a stream.
inline uint64_t uptakeDrawSeed(uint64_t seed, int64_t id,
                               uint64_t campaign_index, uint64_t tier_index,
                               bool shared_unit) {
  return mix_seed(
      seed, static_cast<uint64_t>(id), campaign_index,
      kUptakeDrawTag | (tier_index << 1) | static_cast<uint64_t>(shared_unit));
}

// The draw in [0, 1) from that stream. A pure function of its arguments, so it
// does not depend on loop order or on which rank holds the person.
inline double uptakeDraw(uint64_t seed, int64_t id, uint64_t campaign_index,
                         uint64_t tier_index, bool shared_unit) {
  SplitMix64 rng(
      uptakeDrawSeed(seed, id, campaign_index, tier_index, shared_unit));
  return uniform01(rng);
}

// One kind of failure: how many things were checked, how many failed, the
// lowest failing id and a number that belongs to it for the message.
struct UptakeCheck {
  int64_t checked = 0;
  int64_t failed = 0;
  int64_t lowest = std::numeric_limits<int64_t>::max();
  int64_t detail = 0;

  void note(bool ok, int64_t id, int64_t more = 0) {
    ++checked;
    if (ok) return;
    ++failed;
    if (id < lowest) {
      lowest = id;
      detail = more;
    }
  }
};

struct TierUptake {
  int64_t selected = 0;
  int64_t accepted = 0;
  // Shared units that made one draw for all their members
  int64_t units = 0;
  UptakeCheck unmatched;    // person with no uptake row (lowest person id)
  UptakeCheck venue_count;  // not exactly one venue (lowest person id)
  UptakeCheck split_unit;   // venue missing or elsewhere (lowest person id)
  UptakeCheck no_decider;   // unit with no decider (lowest venue id)

  int64_t declined() const { return selected - accepted; }
};

// `selected` and `accepted` hold one entry per world.people position, bit t for
// tier t. After settleUptake the counts in `tiers` are over the whole world;
// after settleUptakeLocal they are this rank's.
struct CampaignUptake {
  std::vector<uint64_t> selected;
  std::vector<uint64_t> accepted;
  std::vector<TierUptake> tiers;
};

namespace uptake_detail {

// A person who belongs to a shared unit. Sorted by venue, then person id, so a
// unit's members come out in ascending id order whatever order the world lists
// its people in.
struct Member {
  VenueId venue;
  PersonId id;
  uint32_t index;  // position in world.people
  bool in_group;
  bool decider;
};

inline void settleTier(const WorldState& world,
                       const VaccinationCampaignConfig& campaign,
                       uint64_t campaign_index, size_t t, uint64_t seed,
                       CampaignUptake& out) {
  const VaccinationTier& tier = campaign.tiers[t];
  const DecisionRules& rules = tier.decision_rules;
  const bool shared = tier.decision_unit.shared();
  const uint64_t bit = uint64_t{1} << t;
  TierUptake& result = out.tiers[t];

  const auto accept = [&](size_t i) {
    out.accepted[i] |= bit;
    ++result.accepted;
  };
  // A person decides on their own row of the table
  const auto decideAlone = [&](size_t i) {
    const Person& person = world.people[i];
    const auto probability = tier.uptake.probability(person, world);
    result.unmatched.note(probability.has_value(), person.id);
    if (probability &&
        uptakeDraw(seed, person.id, campaign_index, t, false) < *probability)
      accept(i);
  };

  std::vector<Member> members;
  for (size_t i = 0; i < world.people.size(); ++i) {
    const Person& person = world.people[i];
    if (person.is_dead) continue;
    const bool in_group =
        matchesAllCriteria(person, tier.selection_criteria, &world);
    if (in_group) {
      out.selected[i] |= bit;
      ++result.selected;
    }
    const bool in_scope =
        shared && matchesAllCriteria(person, rules.scope, &world);
    if (!in_scope) {
      if (in_group) decideAlone(i);
      continue;
    }
    // Deciders need not be in the group, and members of the group are decided
    // for. Anyone else in scope is not part of this tier's outcome.
    const bool decider = matchesAllCriteria(person, rules.deciders, &world);
    if (!in_group && !decider) continue;

    const auto venues =
        world.getActivityVenues(person, tier.decision_unit.activity_index);
    result.venue_count.note(venues.size() == 1, person.id,
                            static_cast<int64_t>(venues.size()));
    if (venues.size() != 1) continue;
    // Held whole by one rank means the venue and its people share a geo unit.
    // Geo units, not ranks, so this fails the same way at any rank count.
    const VenueId venue_id = venues[0].first;
    const Venue* venue = world.getVenue(venue_id);
    result.split_unit.note(
        venue != nullptr && venue->geo_unit_id == person.geo_unit_id, person.id,
        venue_id);
    members.push_back(
        {venue_id, person.id, static_cast<uint32_t>(i), in_group, decider});
  }
  // ponytail: this index is rebuilt for every shared tier. Share it between
  // tiers with the same activity and scope if startup time ever shows it.
  std::sort(members.begin(), members.end(),
            [](const Member& a, const Member& b) {
              return std::tie(a.venue, a.id) < std::tie(b.venue, b.id);
            });

  for (size_t lo = 0; lo < members.size();) {
    size_t hi = lo;
    while (hi < members.size() && members[hi].venue == members[lo].venue) ++hi;
    const VenueId venue_id = members[lo].venue;
    const auto first = members.begin() + static_cast<std::ptrdiff_t>(lo);
    const auto last = members.begin() + static_cast<std::ptrdiff_t>(hi);
    lo = hi;
    // A unit with nobody to decide for needs no outcome
    if (std::none_of(first, last, [](const Member& m) { return m.in_group; }))
      continue;

    // The mean is summed in ascending person id order, so its last bit does
    // not depend on how the world numbers its people.
    int64_t deciders = 0;
    int64_t answered = 0;
    double sum = 0.0;
    double lowest = 1.0;
    double highest = 0.0;
    for (auto m = first; m != last; ++m) {
      if (!m->decider) continue;
      ++deciders;
      const auto probability =
          tier.uptake.probability(world.people[m->index], world);
      result.unmatched.note(probability.has_value(), m->id);
      if (!probability) continue;
      ++answered;
      sum += *probability;
      lowest = std::min(lowest, *probability);
      highest = std::max(highest, *probability);
    }

    const bool individual = rules.if_no_decider == NoDecider::Individual;
    result.no_decider.note(deciders > 0 || individual, venue_id);
    if (deciders == 0) {
      if (individual)
        for (auto m = first; m != last; ++m)
          if (m->in_group) decideAlone(m->index);
      continue;
    }
    if (answered == 0) continue;  // every decider is already a failure

    const double probability = rules.combine == DeciderCombine::Min ? lowest
                               : rules.combine == DeciderCombine::Max
                                   ? highest
                                   : sum / static_cast<double>(answered);
    ++result.units;
    const bool unit_accepts =
        uptakeDraw(seed, venue_id, campaign_index, t, true) < probability;
    if (unit_accepts)
      for (auto m = first; m != last; ++m)
        if (m->in_group) accept(m->index);
  }
}

inline std::string tierLabel(const VaccinationCampaignConfig& campaign,
                             size_t t) {
  return "vaccination uptake: campaign '" + campaign.name + "' tier '" +
         campaign.tiers[t].name + "': ";
}

inline std::string failedOf(const UptakeCheck& check, const char* what) {
  return std::to_string(check.failed) + " of " + std::to_string(check.checked) +
         " " + what + " checked failed.";
}

// Throws for the first tier with a failure, the same text on every rank. The
// counts must already be summed over the ranks.
inline void throwOnFailure(const VaccinationCampaignConfig& campaign,
                           const CampaignUptake& global) {
  for (size_t t = 0; t < campaign.tiers.size(); ++t) {
    const TierUptake& r = global.tiers[t];
    const std::string where = tierLabel(campaign, t);
    const std::string& activity =
        campaign.tiers[t].decision_unit.shared_activity;
    if (r.venue_count.failed > 0)
      throw std::runtime_error(
          where + "person " + std::to_string(r.venue_count.lowest) + " has " +
          std::to_string(r.venue_count.detail) + " venues for activity '" +
          activity + "', but a shared decision_unit needs exactly 1. " +
          failedOf(r.venue_count, "people") +
          " Fix the world data or use 'decision_unit: person'.");
    if (r.split_unit.failed > 0)
      throw std::runtime_error(
          where + "person " + std::to_string(r.split_unit.lowest) +
          " is registered for activity '" + activity + "' to venue " +
          std::to_string(r.split_unit.detail) +
          ", which is missing or in another geo unit than the person, so the "
          "shared unit is not held whole by one rank. " +
          failedOf(r.split_unit, "people"));
    if (r.unmatched.failed > 0)
      throw std::runtime_error(
          where + "no uptake row matches person " +
          std::to_string(r.unmatched.lowest) +
          ", whose probability is needed. " + failedOf(r.unmatched, "people") +
          " Add a row with empty filters to the uptake table as a catch-all.");
    if (r.no_decider.failed > 0)
      throw std::runtime_error(
          where + "the unit at venue " + std::to_string(r.no_decider.lowest) +
          " for activity '" + activity +
          "' has people in the group but no member matching the deciders, and "
          "if_no_decider is 'error'. " +
          failedOf(r.no_decider, "units") +
          " Widen the deciders or set if_no_decider to 'individual'.");
  }
}

}  // namespace uptake_detail

// One rank's part of the pass: no MPI, counts and failures are this rank's.
// Dead people are ignored. Nothing throws on a failed check; the caller adds
// the counts of every rank and then throws, so no rank is left waiting.
inline CampaignUptake settleUptakeLocal(
    const WorldState& world, const VaccinationCampaignConfig& campaign,
    uint64_t campaign_index, uint64_t seed) {
  if (campaign.tiers.size() > kMaxUptakeTiers)
    throw std::runtime_error("vaccination campaign '" + campaign.name +
                             "' has " + std::to_string(campaign.tiers.size()) +
                             " tiers, but uptake is not built for more than " +
                             std::to_string(kMaxUptakeTiers) +
                             " tiers in one campaign");
  CampaignUptake out;
  out.selected.assign(world.people.size(), 0);
  out.accepted.assign(world.people.size(), 0);
  out.tiers.resize(campaign.tiers.size());
  for (size_t t = 0; t < campaign.tiers.size(); ++t)
    uptake_detail::settleTier(world, campaign, campaign_index, t, seed, out);
  return out;
}

// Collective: every rank calls it with its own world slice and the same
// campaign and seed. Adds up the counts over the ranks, throws if any check
// failed anywhere, and prints one line per tier on rank 0. The masks are this
// rank's; the counts in `tiers` are the whole world's.
inline CampaignUptake settleUptake(const WorldState& world,
                                   const VaccinationCampaignConfig& campaign,
                                   uint64_t campaign_index, uint64_t seed) {
  CampaignUptake out = settleUptakeLocal(world, campaign, campaign_index, seed);
  for (TierUptake& r : out.tiers) {
    r.selected = mpi_runtime::globalSum(r.selected);
    r.accepted = mpi_runtime::globalSum(r.accepted);
    r.units = mpi_runtime::globalSum(r.units);
    for (UptakeCheck* check :
         {&r.unmatched, &r.venue_count, &r.split_unit, &r.no_decider}) {
      const int64_t mine = check->lowest;
      check->checked = mpi_runtime::globalSum(check->checked);
      check->failed = mpi_runtime::globalSum(check->failed);
      check->lowest = mpi_runtime::globalMin(check->lowest);
      // Only the rank that saw the lowest failure knows its detail
      check->detail =
          mpi_runtime::globalSum(mine == check->lowest ? check->detail : 0);
    }
  }
  uptake_detail::throwOnFailure(campaign, out);

  if (mpi_runtime::state().rank == 0) {
    for (size_t t = 0; t < out.tiers.size(); ++t) {
      const TierUptake& r = out.tiers[t];
      std::cout << uptake_detail::tierLabel(campaign, t) << "selected "
                << r.selected << ", accepted " << r.accepted << ", declined "
                << r.declined() << ", shared units " << r.units << ", checked "
                << r.unmatched.checked + r.venue_count.checked +
                       r.split_unit.checked + r.no_decider.checked
                << ", failed 0" << std::endl;
    }
  }
  return out;
}

}  // namespace june
