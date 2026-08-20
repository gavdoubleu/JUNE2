#pragma once

#include <string>
#include <vector>

namespace june {

// One seed budget that could not be filled: nobody eligible anywhere, not
// merely nobody on this rank.
struct SeedShortfall {
  std::string seed_name;
  std::string geo_level;
  std::string unit_id;
  int requested = 0;
  int available = 0;
};

// Format the shortfalls of one seeding step as one warning block, empty when
// nothing fell short. Pure, so rank 0 can emit it with no extra collective:
// every rank sees the same pooled offers and so derives the same records.
std::string formatSeedShortfallReport(
    const std::vector<SeedShortfall>& shortfalls);

}  // namespace june
