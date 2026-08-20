#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace june {

// One seed budget that could not be filled: nobody eligible anywhere, not
// merely nobody on this rank. A unit can fall short against several of its
// budgets at once — overlapping target groups make its budgets compete for
// the same people — so the budget is named alongside the unit.
struct SeedShortfall {
  std::string seed_name;
  std::string geo_level;
  std::string unit_id;
  size_t budget_index = 0;
  int requested = 0;
  int available = 0;
};

// Format the shortfalls of one seeding step as one warning block, empty when
// nothing fell short. Pure, so rank 0 can emit it with no extra collective:
// every rank sees the same pooled offers and so derives the same records.
std::string formatSeedShortfallReport(
    const std::vector<SeedShortfall>& shortfalls);

}  // namespace june
