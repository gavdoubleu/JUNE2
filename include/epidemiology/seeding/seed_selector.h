#pragma once

#include <cstdint>
#include <vector>

#include "core/types.h"

namespace june {

// One rank's offer of one candidate against one seed budget: the candidate's
// deterministic key, derived from the run seed and the person's own identity,
// so the key does not depend on which rank holds the person.
struct SeedOffer {
  uint64_t key = 0;
  PersonId person_id = 0;
  // Which budget of the seed event this offer stands against.
  uint32_t budget_slot = 0;
};

// One person the plan seeds, and the budget the case counts against. A person
// appears at most once: a candidate matching two budgets takes one case, not
// two, whatever the target groups overlap.
struct SeedAssignment {
  PersonId person_id = 0;
  uint32_t budget_index = 0;
};

// The winners of one unit's budgets, plus how many cases each budget could
// actually place — a budget filled short of its request found nobody eligible
// anywhere, not merely nobody on this rank.
struct SeedSelection {
  std::vector<SeedAssignment> chosen;
  std::vector<int> filled_per_budget;
};

// Choose the winners of one unit's seed budgets from the offers every rank
// made. `budget_slot` on an offer is the budget it stands against, indexed
// within the unit.
// Pure: the result depends only on the multiset of offers, never on how they
// were split across the input lists, so every rank reaches the same answer
// from the same pooled offers.
SeedSelection selectSeedWinners(
    const std::vector<std::vector<SeedOffer>>& offer_lists,
    const std::vector<int>& budgets);

}  // namespace june
