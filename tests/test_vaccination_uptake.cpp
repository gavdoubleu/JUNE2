#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "doctest.h"
#include "epidemiology/vaccination_allocation.h"
#include "epidemiology/vaccination_uptake.h"
#include "utils/deterministic_rng.h"
#include "vaccination_uptake_fixtures.h"

using namespace june;
using namespace june::uptake_fixture;

namespace {

constexpr uint64_t kSeed = 20261028;
constexpr uint64_t kCampaignIndex = 2;

using Flags = std::map<PersonId, std::pair<uint64_t, uint64_t>>;

Flags flagsById(const WorldState& world, const CampaignUptake& result) {
  Flags flags;
  for (size_t i = 0; i < world.people.size(); ++i)
    flags[world.people[i].id] = {result.selected[i], result.accepted[i]};
  return flags;
}

bool accepted(const Flags& flags, PersonId id, size_t tier = 0) {
  return (flags.at(id).second >> tier) & 1;
}

bool selected(const Flags& flags, PersonId id, size_t tier = 0) {
  return (flags.at(id).first >> tier) & 1;
}

CampaignUptake settle(const WorldState& world,
                      const VaccinationCampaignConfig& campaign) {
  return settleUptakeLocal(world, campaign, kCampaignIndex, kSeed);
}

// The errors of the simulator, in the one text every rank must agree on
std::string errorOf(const WorldState& world,
                    const VaccinationCampaignConfig& campaign) {
  try {
    settleUptake(world, campaign, kCampaignIndex, kSeed);
  } catch (const std::runtime_error& e) {
    return e.what();
  }
  return "";
}

// What the reduce over ranks gives: counts add, the lowest failing id is the
// smallest, and its detail comes from the rank that saw it.
CampaignUptake sumRanks(const std::vector<CampaignUptake>& ranks) {
  CampaignUptake sum;
  sum.tiers.resize(ranks.at(0).tiers.size());
  for (size_t t = 0; t < sum.tiers.size(); ++t) {
    TierUptake& s = sum.tiers[t];
    const auto checks = [](TierUptake& r) {
      return std::vector<UptakeCheck*>{&r.unmatched, &r.venue_count,
                                       &r.split_unit, &r.no_decider};
    };
    for (const CampaignUptake& rank : ranks) {
      TierUptake r = rank.tiers[t];
      s.selected += r.selected;
      s.accepted += r.accepted;
      s.units += r.units;
      auto from = checks(r);
      auto into = checks(s);
      for (size_t c = 0; c < from.size(); ++c) {
        into[c]->checked += from[c]->checked;
        into[c]->failed += from[c]->failed;
        if (from[c]->lowest < into[c]->lowest) {
          into[c]->lowest = from[c]->lowest;
          into[c]->detail = from[c]->detail;
        }
      }
    }
  }
  return sum;
}

std::vector<CampaignUptake> perRank(
    const std::vector<Home>& homes, size_t ranks,
    const std::function<VaccinationCampaignConfig()>& campaign,
    std::vector<WorldState>* worlds = nullptr) {
  std::vector<CampaignUptake> out;
  for (size_t r = 0; r < ranks; ++r) {
    WorldState world =
        makeWorld(homes, [&](size_t h) { return h % ranks == r; });
    out.push_back(settle(world, campaign()));
    if (worlds) worlds->push_back(std::move(world));
  }
  return out;
}

std::string textOf(const VaccinationCampaignConfig& campaign,
                   const CampaignUptake& global) {
  try {
    uptake_detail::throwOnFailure(campaign, global);
  } catch (const std::runtime_error& e) {
    return e.what();
  }
  return "";
}

double correlation(const std::vector<double>& a, const std::vector<double>& b) {
  const double n = static_cast<double>(a.size());
  double ma = 0, mb = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    ma += a[i];
    mb += b[i];
  }
  ma /= n;
  mb /= n;
  double cov = 0, va = 0, vb = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    cov += (a[i] - ma) * (b[i] - mb);
    va += (a[i] - ma) * (a[i] - ma);
    vb += (b[i] - mb) * (b[i] - mb);
  }
  return cov / std::sqrt(va * vb);
}

VaccinationCampaignConfig childrenDecidedByParents(
    DeciderCombine combine, NoDecider none = NoDecider::Error) {
  VaccinationTier tier = children();
  shareResidence(tier, combine, none);
  tier.uptake.rows = {row({criterion("age", "<", 35.0)}, 0.9), row({}, 0.3)};
  return makeCampaign({tier});
}

}  // namespace

// ---------------------------------------------------------------------------
// A person decides alone
// ---------------------------------------------------------------------------

TEST_CASE("uptake - a person's acceptance rate matches the probability") {
  const auto homes = mixedHomes(20000);
  const WorldState world = makeWorld(homes);
  VaccinationTier tier = adults();
  tier.uptake.rows = {row({}, 0.35)};
  const auto result = settle(world, makeCampaign({tier}));
  const TierUptake& r = result.tiers[0];

  const double n = static_cast<double>(r.selected);
  const double sigma = std::sqrt(n * 0.35 * 0.65);
  CHECK(n > 30000);
  CHECK(std::abs(static_cast<double>(r.accepted) - 0.35 * n) < 4 * sigma);
  CHECK(r.declined() == r.selected - r.accepted);
  CHECK(r.units == 0);
  CHECK(r.unmatched.failed == 0);
}

TEST_CASE("uptake - probability 0 and 1 are exact") {
  const WorldState world = makeWorld(mixedHomes(2000));
  for (double p : {0.0, 1.0}) {
    CAPTURE(p);
    VaccinationTier tier = adults();
    tier.uptake.rows = {row({}, p)};
    const auto result = settle(world, makeCampaign({tier}));
    CHECK(result.tiers[0].selected > 3000);
    CHECK(result.tiers[0].accepted ==
          (p == 1.0 ? result.tiers[0].selected : 0));
  }
}

TEST_CASE(
    "uptake - first matching row wins, outsiders and the dead get no bit") {
  // ids 0 and 1 are adults aged 20 and 50, 2 is a child, 3 is a dead adult
  const std::vector<Home> homes(1, Home{{20, 50, 10, 50}});
  WorldState world = makeWorld(homes);
  world.people[3].is_dead = true;

  VaccinationTier tier = adults();
  tier.uptake.rows = {row({criterion("age", "<", 30.0)}, 0.0), row({}, 1.0)};
  auto flags = flagsById(world, settle(world, makeCampaign({tier})));
  CHECK(selected(flags, 0));
  CHECK_FALSE(accepted(flags, 0));  // the first row, not the catch-all
  CHECK(selected(flags, 1));
  CHECK(accepted(flags, 1));
  CHECK(flags.at(2) == std::make_pair(uint64_t{0}, uint64_t{0}));  // child
  CHECK(flags.at(3) == std::make_pair(uint64_t{0}, uint64_t{0}));  // dead

  // Swap the rows and the catch-all now answers for everyone
  std::swap(tier.uptake.rows[0], tier.uptake.rows[1]);
  flags = flagsById(world, settle(world, makeCampaign({tier})));
  CHECK(accepted(flags, 0));
  CHECK(accepted(flags, 1));

  const auto result = settle(world, makeCampaign({tier}));
  CHECK(result.tiers[0].selected == 2);  // the dead adult is not counted
}

TEST_CASE("uptake - a tier with no uptake accepts everyone in its group") {
  const WorldState world = makeWorld(mixedHomes(5000));
  const auto result = settle(world, makeCampaign({adults()}));
  CHECK(result.tiers[0].selected > 8000);
  CHECK(result.tiers[0].accepted == result.tiers[0].selected);
  CHECK(result.tiers[0].unmatched.failed == 0);
}

TEST_CASE("uptake - a campaign with no tiers settles nothing") {
  const WorldState world = makeWorld(mixedHomes(10));
  const auto result = settle(world, makeCampaign({}));
  CHECK(result.tiers.empty());
  for (uint64_t mask : result.selected) CHECK(mask == 0);
}

TEST_CASE("uptake - more tiers than a mask holds is refused") {
  const WorldState world = makeWorld(mixedHomes(10));
  CHECK_NOTHROW(
      settle(world, makeCampaign(std::vector<VaccinationTier>(64, adults()))));
  const auto text =
      errorOf(world, makeCampaign(std::vector<VaccinationTier>(65, adults())));
  CHECK(text.find("rollout") != std::string::npos);
  CHECK(text.find("65 tiers") != std::string::npos);
  CHECK(text.find("not built") != std::string::npos);
}

// ---------------------------------------------------------------------------
// Repeatability
// ---------------------------------------------------------------------------

TEST_CASE("uptake - same seed same result, another seed other people") {
  const WorldState world = makeWorld(mixedHomes(3000));
  VaccinationTier tier = adults();
  tier.uptake.rows = {row({}, 0.5)};
  const auto campaign = makeCampaign({tier});
  const auto first = settle(world, campaign);
  const auto again = settle(world, campaign);
  CHECK(first.accepted == again.accepted);
  CHECK(first.selected == again.selected);

  const auto other =
      settleUptakeLocal(world, campaign, kCampaignIndex, kSeed + 1);
  CHECK(other.selected == first.selected);
  CHECK(other.accepted != first.accepted);
  size_t differ = 0;
  for (size_t i = 0; i < first.accepted.size(); ++i)
    differ += first.accepted[i] != other.accepted[i];
  CHECK(differ > 1000);  // a fair coin flip each, so about half of 6400
}

TEST_CASE(
    "uptake - the order of world.people does not change anyone's result") {
  const auto homes = mixedHomes(900);
  VaccinationTier kids = children();
  shareResidence(kids, DeciderCombine::Mean);
  kids.uptake.rows = {row({criterion("age", "<", 33.0)}, 0.8), row({}, 0.4)};
  VaccinationTier grown = adults();
  grown.uptake.rows = {row({}, 0.6)};
  const auto campaign = makeCampaign({kids, grown});

  WorldState world = makeWorld(homes);
  const Flags before = flagsById(world, settle(world, campaign));
  SplitMix64 rng(7);
  shuffle_det(world.people.begin(), world.people.end(), rng);
  world.buildIndices();
  const Flags after = flagsById(world, settle(world, campaign));
  CHECK(before == after);
  CHECK(before.size() == world.people.size());
}

// ---------------------------------------------------------------------------
// A shared unit decides together
// ---------------------------------------------------------------------------

TEST_CASE("uptake - the group members of one unit share one outcome") {
  const auto homes = mixedHomes(2000);
  const WorldState world = makeWorld(homes);
  const auto result =
      settle(world, childrenDecidedByParents(DeciderCombine::Mean));
  const Flags flags = flagsById(world, result);

  size_t homes_accepting = 0, homes_declining = 0;
  PersonId id = 0;
  for (const Home& home : homes) {
    size_t child_count = 0, yes_count = 0;
    for (float age : home.ages) {
      if (age < 18 && !home.care) {
        ++child_count;
        yes_count += accepted(flags, id);
      }
      ++id;
    }
    if (child_count == 0) continue;
    CHECK((yes_count == 0 || yes_count == child_count));  // one outcome
    (yes_count > 0 ? homes_accepting : homes_declining)++;
  }
  CHECK(homes_accepting > 100);
  CHECK(homes_declining > 100);
  CHECK(result.tiers[0].units > 0);
  CHECK(result.tiers[0].unmatched.failed == 0);
}

TEST_CASE("uptake - min, mean and max pick the right unit probability") {
  // Each home: adults aged 30 (probability 0.9) and 40 (0.3), one child
  const std::vector<Home> homes(400, Home{{30, 40, 5}});
  const WorldState world = makeWorld(homes);
  std::map<DeciderCombine, Flags> outcomes;
  const std::map<DeciderCombine, double> expected = {
      {DeciderCombine::Min, 0.3},
      {DeciderCombine::Mean, (0.9 + 0.3) / 2.0},
      {DeciderCombine::Max, 0.9}};
  for (const auto& [combine, probability] : expected) {
    const auto result = settle(world, childrenDecidedByParents(combine));
    const Flags flags = flagsById(world, result);
    outcomes[combine] = flags;
    CHECK(result.tiers[0].units == 400);
    for (size_t h = 0; h < homes.size(); ++h) {
      const PersonId child = static_cast<PersonId>(3 * h + 2);
      const double u =
          uptakeDraw(kSeed, static_cast<int64_t>(h), kCampaignIndex, 0, true);
      CHECK(accepted(flags, child) == (u < probability));
    }
  }
  // The three rules really do pick different people
  CHECK(outcomes[DeciderCombine::Min] != outcomes[DeciderCombine::Mean]);
  CHECK(outcomes[DeciderCombine::Mean] != outcomes[DeciderCombine::Max]);
}

TEST_CASE("uptake - a child needs no row but a decider does") {
  // The table only has a row for adults. Children are decided for.
  const std::vector<Home> homes(50, Home{{30, 40, 5}});
  const WorldState world = makeWorld(homes);
  VaccinationTier kids = children();
  shareResidence(kids, DeciderCombine::Mean);
  kids.uptake.rows = {row({criterion("age", ">=", 18.0)}, 0.5)};
  const auto result = settle(world, makeCampaign({kids}));
  CHECK(result.tiers[0].unmatched.failed == 0);
  CHECK(result.tiers[0].unmatched.checked ==
        100);  // the two adults of 50 homes

  // The same table for children alone: a child has no row
  VaccinationTier alone = children();
  alone.uptake.rows = kids.uptake.rows;
  CHECK(settle(world, makeCampaign({alone})).tiers[0].unmatched.failed == 50);

  // A decider with no row is a failure
  VaccinationTier gap = children();
  shareResidence(gap, DeciderCombine::Mean);
  gap.uptake.rows = {row({criterion("age", ">=", 35.0)}, 0.5)};
  CHECK(settle(world, makeCampaign({gap})).tiers[0].unmatched.failed == 50);
}

TEST_CASE("uptake - people outside the scope decide on their own") {
  const auto homes = mixedHomes(700);
  const WorldState world = makeWorld(homes);
  // Everyone in the group, adults decide for the homes, care homes are outside
  VaccinationTier tier = makeTier("everyone", {});
  shareResidence(tier, DeciderCombine::Mean);
  tier.uptake.rows = {row({}, 0.5)};
  const Flags flags = flagsById(world, settle(world, makeCampaign({tier})));

  bool residents_differ = false;
  PersonId id = 0;
  for (const Home& home : homes) {
    if (home.care) {
      for (size_t k = 0; k < home.ages.size(); ++k) {
        const double u = uptakeDraw(kSeed, id + static_cast<int64_t>(k),
                                    kCampaignIndex, 0, false);
        CHECK(accepted(flags, id + static_cast<PersonId>(k)) == (u < 0.5));
      }
      residents_differ = residents_differ ||
                         accepted(flags, id) != accepted(flags, id + 1) ||
                         accepted(flags, id) != accepted(flags, id + 2);
    }
    id += static_cast<PersonId>(home.ages.size());
  }
  CHECK(residents_differ);
}

TEST_CASE(
    "uptake - a unit with no decider is an error or decides person by person") {
  // Homes of children only
  const std::vector<Home> homes(100, Home{{5, 7}});
  const WorldState world = makeWorld(homes);

  const auto strict =
      settle(world, childrenDecidedByParents(DeciderCombine::Mean));
  CHECK(strict.tiers[0].no_decider.failed == 100);
  CHECK(strict.tiers[0].no_decider.checked == 100);

  const auto campaign =
      childrenDecidedByParents(DeciderCombine::Mean, NoDecider::Individual);
  const auto individual = settle(world, campaign);
  CHECK(individual.tiers[0].no_decider.failed == 0);
  CHECK(individual.tiers[0].unmatched.failed == 0);
  CHECK(individual.tiers[0].units == 0);  // no unit drew for them
  const Flags flags = flagsById(world, individual);
  size_t siblings_differ = 0;
  for (PersonId id = 0; id < 200; ++id) {
    const double u = uptakeDraw(kSeed, id, kCampaignIndex, 0, false);
    CHECK(accepted(flags, id) == (u < 0.9));  // children are under 35: row 1
    if (id % 2 == 0)
      siblings_differ += accepted(flags, id) != accepted(flags, id + 1);
  }
  CHECK(siblings_differ > 0);
}

// ---------------------------------------------------------------------------
// The four errors
// ---------------------------------------------------------------------------

namespace {

void expectNames(const std::string& text,
                 const std::vector<std::string>& parts) {
  CAPTURE(text);
  for (const auto& part : parts) CHECK(text.find(part) != std::string::npos);
}

void expectSameTextAtEveryRankCount(
    const std::vector<Home>& homes,
    const std::function<VaccinationCampaignConfig()>& campaign,
    const std::string& whole_text) {
  for (size_t ranks : {1, 2, 3, 5}) {
    CAPTURE(ranks);
    const auto sum = sumRanks(perRank(homes, ranks, campaign));
    CHECK(textOf(campaign(), sum) == whole_text);
  }
}

}  // namespace

TEST_CASE("uptake error - a person the table cannot answer") {
  const auto campaign = [] {
    VaccinationTier tier = adults();
    tier.uptake.rows = {row({criterion("age", "<", 35.0)}, 0.5)};
    return makeCampaign({tier});
  };
  // Adults aged 30 are covered, 40 are not: the first failure is home 0's
  // second adult, person 1
  const auto homes = homesWith(Home{{30, 40}}, 3, 2);
  const WorldState world = makeWorld(homes);
  const std::string text = errorOf(world, campaign());
  expectNames(text, {"campaign 'rollout'", "tier 'adults'", "person 1",
                     "7 of 14 people checked failed"});
  expectSameTextAtEveryRankCount(homes, campaign, text);
}

TEST_CASE("uptake error - zero or several venues for the activity") {
  const auto campaign = [] {
    VaccinationTier tier = children();
    shareResidence(tier, DeciderCombine::Mean);
    // A scope on the venue type would leave out anyone with no venue at all
    tier.decision_rules.scope = {criterion("age", ">=", 0.0)};
    return makeCampaign({tier});
  };
  Home none{{30, 40, 5}};
  none.registered = false;
  auto homes = homesWith(none, 3, 2);
  WorldState world = makeWorld(homes);
  // ids: 3 homes of 3 people come first, so the first bad home starts at 9
  std::string text = errorOf(world, campaign());
  expectNames(text, {"campaign 'rollout'", "tier 'children'",
                     "person 9 has 0 venues", "6 of 21 people checked failed"});
  expectSameTextAtEveryRankCount(homes, campaign, text);

  Home several{{30, 40, 5}};
  several.extra_venues = {0};
  homes = homesWith(several, 3, 2);
  world = makeWorld(homes);
  text = errorOf(world, campaign());
  expectNames(
      text, {"campaign 'rollout'", "tier 'children'", "person 9 has 2 venues"});
  expectSameTextAtEveryRankCount(homes, campaign, text);
}

TEST_CASE("uptake error - a unit with no decider") {
  const auto campaign = [] {
    return childrenDecidedByParents(DeciderCombine::Mean);
  };
  const auto homes = homesWith(Home{{5, 7}}, 3, 2);
  const WorldState world = makeWorld(homes);
  const std::string text = errorOf(world, campaign());
  expectNames(text, {"campaign 'rollout'", "tier 'children'", "venue 3",
                     "2 of 7 units checked failed", "if_no_decider"});
  expectSameTextAtEveryRankCount(homes, campaign, text);
}

TEST_CASE("uptake error - a unit not held whole") {
  const auto campaign = [] {
    VaccinationTier tier = children();
    shareResidence(tier, DeciderCombine::Mean);
    return makeCampaign({tier});
  };
  // The residents live in another geo unit than their venue
  Home elsewhere{{30, 40, 5}};
  elsewhere.person_geo_unit = 999;
  auto homes = homesWith(elsewhere, 3, 2);
  WorldState world = makeWorld(homes);
  std::string text = errorOf(world, campaign());
  expectNames(text, {"campaign 'rollout'", "tier 'children'", "person 9",
                     "venue 3", "6 of 21 people checked failed"});
  expectSameTextAtEveryRankCount(homes, campaign, text);

  // The venue is not in the world at all
  Home missing{{30, 40, 5}};
  missing.venue_loaded = false;
  homes = homesWith(missing, 3, 2);
  world = makeWorld(homes);
  text = errorOf(world, campaign());
  expectNames(text, {"person 9", "venue 3", "missing"});
  expectSameTextAtEveryRankCount(homes, campaign, text);
}

TEST_CASE("uptake - a world that passes every check names nothing") {
  const WorldState world = makeWorld(mixedHomes(300));
  CHECK(errorOf(world, childrenDecidedByParents(DeciderCombine::Mean)) == "");
}

// ---------------------------------------------------------------------------
// Independence of the draws
// ---------------------------------------------------------------------------

TEST_CASE("uptake draw - not the lottery key, not another tier, not a venue") {
  constexpr size_t n = 200000;
  const double limit = 4.0 / std::sqrt(static_cast<double>(n));
  std::vector<double> person_draw(n), lottery(n), next_tier(n), venue_draw(n);
  for (size_t i = 0; i < n; ++i) {
    const auto id = static_cast<int64_t>(i);
    person_draw[i] = uptakeDraw(kSeed, id, kCampaignIndex, 1, false);
    lottery[i] = static_cast<double>(lotteryKey(kSeed, i, kCampaignIndex, 1)) *
                 0x1.0p-64;
    next_tier[i] = uptakeDraw(kSeed, id, kCampaignIndex, 2, false);
    venue_draw[i] = uptakeDraw(kSeed, id, kCampaignIndex, 1, true);
    CHECK(person_draw[i] >= 0.0);
    CHECK(person_draw[i] < 1.0);
  }
  CHECK(std::abs(correlation(person_draw, lottery)) < limit);
  CHECK(std::abs(correlation(person_draw, next_tier)) < limit);
  CHECK(std::abs(correlation(person_draw, venue_draw)) < limit);
  size_t same = 0;
  for (size_t i = 0; i < n; ++i)
    same += person_draw[i] == lottery[i] || person_draw[i] == next_tier[i] ||
            person_draw[i] == venue_draw[i];
  CHECK(same == 0);

  // The stream itself is not the lottery's: SplitMix64 would hide a shared seed
  // from the draws above, but not from the seeds.
  size_t shared_seed = 0;
  for (size_t i = 0; i < n; ++i)
    for (uint64_t tier : {0, 1, 2, 3})
      shared_seed +=
          uptakeDrawSeed(kSeed, static_cast<int64_t>(i), kCampaignIndex, tier,
                         false) == lotteryKey(kSeed, i, kCampaignIndex, tier);
  CHECK(shared_seed == 0);
}

// ---------------------------------------------------------------------------
// The same answer at any rank count
// ---------------------------------------------------------------------------

TEST_CASE("uptake - slices by home give the whole world's result") {
  const auto homes = mixedHomes(1200);
  const auto campaign = [] {
    VaccinationTier kids = children();
    shareResidence(kids, DeciderCombine::Mean);
    kids.uptake.rows = {row({criterion("age", "<", 33.0)}, 0.8), row({}, 0.4)};
    VaccinationTier grown = adults();
    grown.uptake.rows = {row({criterion("age", "<", 45.0)}, 0.6), row({}, 0.2)};
    VaccinationTier everyone = makeTier("everyone", {});
    shareResidence(everyone, DeciderCombine::Max);
    everyone.uptake.rows = {row({}, 0.5)};
    return makeCampaign({kids, grown, everyone});
  };
  const WorldState whole = makeWorld(homes);
  const auto expected = settle(whole, campaign());
  const Flags expected_flags = flagsById(whole, expected);
  for (const auto& r : expected.tiers) CHECK(r.selected > 0);

  for (size_t ranks : {1, 2, 3}) {
    CAPTURE(ranks);
    std::vector<WorldState> worlds;
    const auto parts = perRank(homes, ranks, campaign, &worlds);
    Flags together;
    for (size_t r = 0; r < ranks; ++r) {
      const Flags mine = flagsById(worlds[r], parts[r]);
      together.insert(mine.begin(), mine.end());
    }
    CHECK(together == expected_flags);
    const auto sum = sumRanks(parts);
    for (size_t t = 0; t < expected.tiers.size(); ++t) {
      CHECK(sum.tiers[t].selected == expected.tiers[t].selected);
      CHECK(sum.tiers[t].accepted == expected.tiers[t].accepted);
      CHECK(sum.tiers[t].units == expected.tiers[t].units);
    }
  }
}
