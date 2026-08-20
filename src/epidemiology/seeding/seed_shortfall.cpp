#include "epidemiology/seeding/seed_shortfall.h"

#include <algorithm>
#include <sstream>

namespace june {

namespace {
// A shortfall touching every unit of a national seed must not flood the log,
// and must not understate how wide it was either.
constexpr size_t kMaxReportedUnits = 10;
}  // namespace

std::string formatSeedShortfallReport(
    const std::vector<SeedShortfall>& shortfalls) {
  std::ostringstream report;
  const size_t reported = std::min(shortfalls.size(), kMaxReportedUnits);
  for (size_t unit = 0; unit < reported; ++unit) {
    const auto& shortfall = shortfalls[unit];
    report << "    [SEED SHORTFALL] seed '" << shortfall.seed_name << "' unit '"
           << shortfall.unit_id << "' (" << shortfall.geo_level
           << "): requested " << shortfall.requested << ", available "
           << shortfall.available << "\n";
  }
  if (shortfalls.size() > reported) {
    report << "    [SEED SHORTFALL] " << (shortfalls.size() - reported)
           << " further units suppressed\n";
  }
  return report.str();
}

}  // namespace june
