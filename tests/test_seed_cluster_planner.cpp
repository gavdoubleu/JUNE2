#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <vector>

#include "../include/epidemiology/seeding/seed_cluster_planner.h"

using namespace june;

namespace {

// Offers for one unit, as the ranks holding the households would make them:
// one per (candidate, budget the candidate matches), keyed by household.
// A candidate matching no budget still counts towards its household's size.
const std::vector<SeedOffer> kOffers = {
    {0x30ULL, 1, 1}, {0x30ULL, 2, 1}, {0x30ULL, 3, 1},  // three matched, size 3
    {0x10ULL, 4, 2}, {0x10ULL, 5, 2},                   // two matched, size 2
    {0x20ULL, 6, 1}, {0x20ULL, 7, 0},                   // one matched, size 2
    {0x40ULL, 8, 1}, {0x40ULL, 8, 2},                   // one matched, size 1
};

std::vector<std::vector<SeedOffer>> splitInto(const std::vector<SeedOffer>& all,
                                              size_t num_lists) {
  std::vector<std::vector<SeedOffer>> lists(num_lists);
  for (size_t i = 0; i < all.size(); ++i) {
    lists[i % num_lists].push_back(all[i]);
  }
  return lists;
}

std::vector<std::pair<PersonId, uint32_t>> named(const ClusterPlan& plan) {
  std::vector<std::pair<PersonId, uint32_t>> pairs;
  for (const auto& assignment : plan.assignments) {
    pairs.push_back({assignment.person_id, assignment.budget_index});
  }
  return pairs;
}

}  // namespace

TEST_CASE("Planner: the plan does not depend on how offers are partitioned") {
  const std::vector<int> budgets = {2, 2};

  ClusterPlan from_one_list = planClusteredSeed(splitInto(kOffers, 1), budgets);
  ClusterPlan from_three_lists =
      planClusteredSeed(splitInto(kOffers, 3), budgets);

  // Densest household first, then the next, until the total budget is met.
  CHECK(named(from_one_list) ==
        std::vector<std::pair<PersonId, uint32_t>>{{1, 0}, {2, 0}, {4, 1},
                                                   {5, 1}});
  CHECK(named(from_three_lists) == named(from_one_list));
  CHECK(from_one_list.filled_per_budget == std::vector<int>{2, 2});
}

TEST_CASE("Planner: households of equal density are ordered by key") {
  // Two households of one matched member each: only the key can separate them.
  const std::vector<SeedOffer> low_key_first = {{0x11ULL, 1, 1},
                                                {0x22ULL, 2, 1}};
  const std::vector<SeedOffer> high_key_first = {{0x22ULL, 2, 1},
                                                 {0x11ULL, 1, 1}};

  ClusterPlan offered_low_first = planClusteredSeed({low_key_first}, {1});
  ClusterPlan offered_high_first = planClusteredSeed({high_key_first}, {1});
  ClusterPlan offered_apart =
      planClusteredSeed({{high_key_first[0]}, {high_key_first[1]}}, {1});

  CHECK(named(offered_low_first) ==
        std::vector<std::pair<PersonId, uint32_t>>{{1, 0}});
  CHECK(named(offered_high_first) == named(offered_low_first));
  CHECK(named(offered_apart) == named(offered_low_first));
}

TEST_CASE("Planner: the fill respects each budget and stops at their total") {
  // One case for the first budget, none for the second.
  const std::vector<SeedOffer> offers = {
      {0x50ULL, 1, 1}, {0x50ULL, 2, 1}, {0x50ULL, 3, 2}, {0x60ULL, 4, 1}};

  ClusterPlan plan = planClusteredSeed({offers}, {1, 0});

  // The first member takes the only case; their housemates find no budget
  // still open, and the next household is never reached.
  CHECK(named(plan) == std::vector<std::pair<PersonId, uint32_t>>{{1, 0}});
  CHECK(plan.filled_per_budget == std::vector<int>{1, 0});
}

TEST_CASE("Planner: a household offered by several ranks is planned as one") {
  // 0x70 straddles three ranks. Merged it is the denser household and takes
  // the whole budget; taken a rank at a time each piece would rank below the
  // rival 0x69 and lose the fill to it.
  const std::vector<std::vector<SeedOffer>> offer_lists = {
      {{0x70ULL, 1, 1}, {0x69ULL, 4, 1}},
      {{0x70ULL, 2, 1}, {0x69ULL, 5, 1}},
      {{0x70ULL, 3, 1}}};

  ClusterPlan plan = planClusteredSeed(offer_lists, {3});

  CHECK(named(plan) == std::vector<std::pair<PersonId, uint32_t>>{
                           {1, 0}, {2, 0}, {3, 0}});
}
