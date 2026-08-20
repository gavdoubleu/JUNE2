#include "epidemiology/seeding/seed_selector.h"

#include <algorithm>
#include <unordered_set>

namespace june {

SeedSelection selectSeedWinners(
    const std::vector<std::vector<SeedOffer>>& offer_lists,
    const std::vector<int>& budgets) {
  std::vector<SeedOffer> pooled;
  for (const auto& offers : offer_lists) {
    pooled.insert(pooled.end(), offers.begin(), offers.end());
  }

  std::sort(pooled.begin(), pooled.end(),
            [](const SeedOffer& a, const SeedOffer& b) {
              // The person breaks ties, so equal keys cannot be reordered by
              // which rank happened to offer them.
              if (a.key != b.key) return a.key < b.key;
              if (a.person_id != b.person_id) return a.person_id < b.person_id;
              return a.budget_slot < b.budget_slot;
            });

  SeedSelection selection;
  selection.filled_per_budget.assign(budgets.size(), 0);
  // A candidate matching two budgets is offered against both, so the first
  // budget to reach them takes them and the other refills from its next-best
  // offer. Which budget that is falls out of the keys, not out of the order
  // the budgets were declared in.
  std::unordered_set<PersonId> already_seeded;
  for (const auto& offer : pooled) {
    if (offer.budget_slot >= budgets.size()) continue;
    if (already_seeded.count(offer.person_id) != 0) continue;
    if (selection.filled_per_budget[offer.budget_slot] >=
        budgets[offer.budget_slot]) {
      continue;
    }
    selection.chosen.push_back({offer.person_id, offer.budget_slot});
    already_seeded.insert(offer.person_id);
    ++selection.filled_per_budget[offer.budget_slot];
  }
  return selection;
}

}  // namespace june
