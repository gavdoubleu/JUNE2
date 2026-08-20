#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <algorithm>
#include <vector>

#include "../include/epidemiology/seeding/seed_selector.h"

using namespace june;

namespace {

// Offers as they would arrive from ranks: an arbitrary key per candidate.
const std::vector<SeedOffer> kOffers = {
    {0x9A11ULL, 101}, {0x0002ULL, 102}, {0x7F30ULL, 103}, {0x0100ULL, 104},
    {0xFFFFULL, 105}, {0x0041ULL, 106}, {0x5555ULL, 107}, {0x0003ULL, 108},
};

std::vector<std::vector<SeedOffer>> splitInto(const std::vector<SeedOffer>& all,
                                              size_t num_lists) {
  std::vector<std::vector<SeedOffer>> lists(num_lists);
  for (size_t i = 0; i < all.size(); ++i) {
    lists[i % num_lists].push_back(all[i]);
  }
  return lists;
}

}  // namespace

TEST_CASE("Selector: winners do not depend on how offers are partitioned") {
  const int budget = 4;

  std::vector<PersonId> from_one_list =
      selectSeedWinners(splitInto(kOffers, 1), budget).chosen;
  std::vector<PersonId> from_two_lists =
      selectSeedWinners(splitInto(kOffers, 2), budget).chosen;
  std::vector<PersonId> from_five_lists =
      selectSeedWinners(splitInto(kOffers, 5), budget).chosen;

  CHECK(from_one_list.size() == 4);
  CHECK(from_two_lists == from_one_list);
  CHECK(from_five_lists == from_one_list);
}

TEST_CASE("Selector: a budget beyond the offers takes them all and reports the shortfall") {
  SeedSelection selection = selectSeedWinners(splitInto(kOffers, 3), 11);

  CHECK(selection.chosen.size() == kOffers.size());
  CHECK(selection.shortfall == 3);
}

TEST_CASE("Selector: a budget of zero seeds nobody and is not short") {
  SeedSelection selection = selectSeedWinners(splitInto(kOffers, 3), 0);

  CHECK(selection.chosen.empty());
  CHECK(selection.shortfall == 0);
}

TEST_CASE("Selector: a candidate offered by one rank alone can still win") {
  // The lone rank's candidate holds the lowest key of the whole pool.
  std::vector<std::vector<SeedOffer>> offer_lists = splitInto(kOffers, 2);
  offer_lists.push_back({{0x0001ULL, 999}});

  SeedSelection selection = selectSeedWinners(offer_lists, 1);

  CHECK(selection.chosen == std::vector<PersonId>{999});
}

TEST_CASE("Selector: equal keys break ties by person, not by input order") {
  // Twenty candidates sharing one key: only the tiebreak can order them.
  std::vector<SeedOffer> tied;
  for (PersonId person_id = 20; person_id >= 1; --person_id) {
    tied.push_back({0x4242ULL, person_id});
  }

  SeedSelection descending_input = selectSeedWinners({tied}, 3);
  std::reverse(tied.begin(), tied.end());
  SeedSelection ascending_input = selectSeedWinners({tied}, 3);
  SeedSelection split_input = selectSeedWinners(splitInto(tied, 4), 3);

  CHECK(descending_input.chosen == std::vector<PersonId>{1, 2, 3});
  CHECK(ascending_input.chosen == descending_input.chosen);
  CHECK(split_input.chosen == descending_input.chosen);
}

TEST_CASE("Selector: a unit held by other ranks alone is not a shortfall") {
  // The pooled offers are what separate "nobody eligible anywhere" from
  // "nobody eligible here": a rank holding none of a unit still sees the
  // offers of the ranks that do, so it reports no shortfall.
  std::vector<std::vector<SeedOffer>> offer_lists = {{}, kOffers};

  SeedSelection selection = selectSeedWinners(offer_lists, 3);

  CHECK(selection.chosen.size() == 3);
  CHECK(selection.shortfall == 0);
}
