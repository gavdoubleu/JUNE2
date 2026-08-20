#include "epidemiology/seeding/seed_selector.h"

#include <algorithm>

namespace june {

SeedSelection selectSeedWinners(
    const std::vector<std::vector<SeedOffer>>& offer_lists, int budget) {
  std::vector<SeedOffer> pooled;
  for (const auto& offers : offer_lists) {
    pooled.insert(pooled.end(), offers.begin(), offers.end());
  }

  std::sort(pooled.begin(), pooled.end(),
            [](const SeedOffer& a, const SeedOffer& b) {
              // The person breaks ties, so equal keys cannot be reordered by
              // which rank happened to offer them.
              if (a.key != b.key) return a.key < b.key;
              return a.person_id < b.person_id;
            });

  SeedSelection selection;
  for (size_t i = 0; i < pooled.size() && (int)selection.chosen.size() < budget;
       ++i) {
    selection.chosen.push_back(pooled[i].person_id);
  }
  selection.shortfall = budget - static_cast<int>(selection.chosen.size());
  return selection;
}

}  // namespace june
