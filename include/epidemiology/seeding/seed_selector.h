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

// The winners of one budget, plus the number of cases the offers could not
// cover — nobody eligible anywhere, not merely nobody on this rank.
struct SeedSelection {
  std::vector<PersonId> chosen;
  int shortfall = 0;
};

// Choose the winners of one seed budget from the offers every rank made.
// Pure: the result depends only on the multiset of offers, never on how they
// were split across the input lists, so every rank reaches the same answer
// from the same pooled offers.
SeedSelection selectSeedWinners(
    const std::vector<std::vector<SeedOffer>>& offer_lists, int budget);

}  // namespace june
