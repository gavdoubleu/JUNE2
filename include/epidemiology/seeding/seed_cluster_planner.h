#pragma once

#include <cstdint>
#include <vector>

#include "core/types.h"
#include "epidemiology/seeding/seed_selector.h"

namespace june {

// One person the plan seeds, and the budget the case counts against.
struct ClusterAssignment {
  PersonId person_id = 0;
  uint32_t budget_index = 0;
};

// The households a clustered seed fills, expressed as the people to infect.
// Every rank replays the same plan and infects only the people it holds, so a
// household straddling ranks needs no special case.
struct ClusterPlan {
  std::vector<ClusterAssignment> assignments;
  std::vector<int> filled_per_budget;
};

// Plan one unit's clustered seed from the offers every rank made.
// One offer per (candidate, budget the candidate matches): `key` is the
// candidate's household, shared by all its members, and `budget_slot` is the
// budget index plus one, or zero for a candidate matching no budget at all —
// which still counts towards its household's size, never its matched members.
// Pure: the result depends only on the multiset of offers, never on how they
// were split across the input lists.
ClusterPlan planClusteredSeed(
    const std::vector<std::vector<SeedOffer>>& offer_lists,
    const std::vector<int>& budgets);

}  // namespace june
