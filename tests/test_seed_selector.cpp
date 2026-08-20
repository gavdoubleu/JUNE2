#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <algorithm>
#include <vector>

#include "../include/epidemiology/seeding/seed_selector.h"
#include "doctest.h"

using namespace june;

namespace {

// Offers as they would arrive from ranks: an arbitrary key per candidate, all
// standing against the unit's only budget.
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

std::vector<PersonId> chosenPeople(const SeedSelection& selection) {
  std::vector<PersonId> people;
  for (const auto& assignment : selection.chosen) {
    people.push_back(assignment.person_id);
  }
  return people;
}

}  // namespace

TEST_CASE("Selector: winners do not depend on how offers are partitioned") {
  const std::vector<int> budgets = {4};

  std::vector<PersonId> from_one_list =
      chosenPeople(selectSeedWinners(splitInto(kOffers, 1), budgets));
  std::vector<PersonId> from_two_lists =
      chosenPeople(selectSeedWinners(splitInto(kOffers, 2), budgets));
  std::vector<PersonId> from_five_lists =
      chosenPeople(selectSeedWinners(splitInto(kOffers, 5), budgets));

  CHECK(from_one_list.size() == 4);
  CHECK(from_two_lists == from_one_list);
  CHECK(from_five_lists == from_one_list);
}

TEST_CASE(
    "Selector: a budget beyond the offers takes them all and reports the "
    "shortfall") {
  SeedSelection selection = selectSeedWinners(splitInto(kOffers, 3), {11});

  CHECK(selection.chosen.size() == kOffers.size());
  CHECK(selection.filled_per_budget == std::vector<int>{8});
}

TEST_CASE("Selector: a budget of zero seeds nobody and is not short") {
  SeedSelection selection = selectSeedWinners(splitInto(kOffers, 3), {0});

  CHECK(selection.chosen.empty());
  CHECK(selection.filled_per_budget == std::vector<int>{0});
}

TEST_CASE("Selector: a candidate offered by one rank alone can still win") {
  // The lone rank's candidate holds the lowest key of the whole pool.
  std::vector<std::vector<SeedOffer>> offer_lists = splitInto(kOffers, 2);
  offer_lists.push_back({{0x0001ULL, 999}});

  SeedSelection selection = selectSeedWinners(offer_lists, {1});

  CHECK(chosenPeople(selection) == std::vector<PersonId>{999});
}

TEST_CASE("Selector: equal keys break ties by person, not by input order") {
  // Twenty candidates sharing one key: only the tiebreak can order them.
  std::vector<SeedOffer> tied;
  for (PersonId person_id = 20; person_id >= 1; --person_id) {
    tied.push_back({0x4242ULL, person_id});
  }

  SeedSelection descending_input = selectSeedWinners({tied}, {3});
  std::reverse(tied.begin(), tied.end());
  SeedSelection ascending_input = selectSeedWinners({tied}, {3});
  SeedSelection split_input = selectSeedWinners(splitInto(tied, 4), {3});

  CHECK(chosenPeople(descending_input) == std::vector<PersonId>{1, 2, 3});
  CHECK(chosenPeople(ascending_input) == chosenPeople(descending_input));
  CHECK(chosenPeople(split_input) == chosenPeople(descending_input));
}

TEST_CASE("Selector: a unit held by other ranks alone is not a shortfall") {
  // The pooled offers are what separate "nobody eligible anywhere" from
  // "nobody eligible here": a rank holding none of a unit still sees the
  // offers of the ranks that do, so it reports no shortfall.
  std::vector<std::vector<SeedOffer>> offer_lists = {{}, kOffers};

  SeedSelection selection = selectSeedWinners(offer_lists, {3});

  CHECK(selection.chosen.size() == 3);
  CHECK(selection.filled_per_budget == std::vector<int>{3});
}

TEST_CASE(
    "Selector: a candidate matching two budgets takes one case, not two") {
  // Overlapping target groups offer the same person against both budgets,
  // under a different key each time.
  std::vector<SeedOffer> offers = {
      {0x0010ULL, 201, 0},
      {0x0020ULL, 201, 1},
      {0x0030ULL, 202, 0},
      {0x0040ULL, 203, 1},
  };

  SeedSelection selection = selectSeedWinners({offers}, {1, 1});

  CHECK(selection.chosen.size() == 2);
  CHECK(chosenPeople(selection) == std::vector<PersonId>{201, 203});
  CHECK(selection.filled_per_budget == std::vector<int>{1, 1});
}

TEST_CASE(
    "Selector: the budget that loses a candidate refills from its next best "
    "offer") {
  // 301 is the best offer of both budgets; budget 1 must fall through to 302
  // rather than place one case fewer than it asked for.
  std::vector<SeedOffer> offers = {
      {0x0001ULL, 301, 0},
      {0x0002ULL, 301, 1},
      {0x0009ULL, 302, 1},
  };

  SeedSelection selection = selectSeedWinners({offers}, {1, 1});

  CHECK(chosenPeople(selection) == std::vector<PersonId>{301, 302});
  CHECK(selection.chosen[0].budget_index == 0);
  CHECK(selection.chosen[1].budget_index == 1);
  CHECK(selection.filled_per_budget == std::vector<int>{1, 1});
}

TEST_CASE(
    "Selector: a contested candidate goes to the budget that keys them best") {
  // 401 keys far better against budget 1 than budget 0, so budget 1 takes them
  // even though budget 0 was declared first, and budget 0 refills with 402.
  std::vector<SeedOffer> offers = {
      {0x0050ULL, 401, 0},
      {0x0005ULL, 401, 1},
      {0x0060ULL, 402, 0},
  };

  SeedSelection selection = selectSeedWinners({offers}, {1, 1});

  CHECK(chosenPeople(selection) == std::vector<PersonId>{401, 402});
  CHECK(selection.chosen[0].budget_index == 1);
  CHECK(selection.chosen[1].budget_index == 0);
}

TEST_CASE("Selector: budgets sharing too few candidates fall short honestly") {
  // Three people, both budgets asking two, and every person matches both: the
  // seed can only place three cases, and says so.
  std::vector<SeedOffer> offers;
  for (PersonId person_id = 501; person_id <= 503; ++person_id) {
    offers.push_back({0x0100ULL + person_id, person_id, 0});
    offers.push_back({0x0200ULL + person_id, person_id, 1});
  }

  SeedSelection selection = selectSeedWinners({offers}, {2, 2});

  CHECK(selection.chosen.size() == 3);
  CHECK(selection.filled_per_budget[0] + selection.filled_per_budget[1] == 3);
}

TEST_CASE(
    "Selector: overlapping budgets resolve the same however the offers are "
    "split") {
  std::vector<SeedOffer> offers;
  for (PersonId person_id = 601; person_id <= 612; ++person_id) {
    offers.push_back({0x3000ULL + (person_id * 7919) % 4096, person_id, 0});
    if (person_id % 2 == 0) {
      offers.push_back({0x3000ULL + (person_id * 104729) % 4096, person_id, 1});
    }
  }

  SeedSelection from_one_list = selectSeedWinners({offers}, {3, 2});
  SeedSelection from_four_lists =
      selectSeedWinners(splitInto(offers, 4), {3, 2});

  CHECK(from_one_list.chosen.size() == 5);
  CHECK(chosenPeople(from_four_lists) == chosenPeople(from_one_list));
  CHECK(from_four_lists.filled_per_budget == from_one_list.filled_per_budget);
}
