#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <string>
#include <vector>

#include "../include/epidemiology/seeding/seed_shortfall.h"

using namespace june;

TEST_CASE("Shortfall report: a short unit names seed, unit, level and counts") {
  std::vector<SeedShortfall> shortfalls = {
      {"february_2020", "LGU", "E06000005", 20, 12}};

  std::string report = formatSeedShortfallReport(shortfalls);

  CHECK(report.find("february_2020") != std::string::npos);
  CHECK(report.find("E06000005") != std::string::npos);
  CHECK(report.find("LGU") != std::string::npos);
  CHECK(report.find("20") != std::string::npos);
  CHECK(report.find("12") != std::string::npos);
}

namespace {

int countOccurrences(const std::string& text, const std::string& needle) {
  int count = 0;
  for (size_t at = text.find(needle); at != std::string::npos;
       at = text.find(needle, at + needle.size())) {
    ++count;
  }
  return count;
}

std::vector<SeedShortfall> manyShortfalls(int units) {
  std::vector<SeedShortfall> shortfalls;
  for (int unit = 0; unit < units; ++unit) {
    shortfalls.push_back({"february_2020", "LGU",
                          "E060000" + std::to_string(unit), 20, 12});
  }
  return shortfalls;
}

}  // namespace

TEST_CASE("Shortfall report: beyond the cap it counts what it suppressed") {
  std::string report = formatSeedShortfallReport(manyShortfalls(17));

  CHECK(countOccurrences(report, "requested") == 10);
  CHECK(report.find("7 further units suppressed") != std::string::npos);
}

TEST_CASE("Shortfall report: a unit eligible nowhere reports like any other") {
  // A mistyped unit code resolves to nobody at all. It is the shortfall most
  // worth seeing, so it gets a line of its own rather than being swallowed.
  std::vector<SeedShortfall> shortfalls = {
      {"february_2020", "LGU", "E06000005", 20, 12},
      {"february_2020", "LGU", "E06TYPO", 15, 0}};

  std::string report = formatSeedShortfallReport(shortfalls);

  CHECK(countOccurrences(report, "requested") == 2);
  CHECK(report.find("'E06TYPO' (LGU): requested 15, available 0") !=
        std::string::npos);
}

TEST_CASE("Shortfall report: nothing short emits no block") {
  CHECK(formatSeedShortfallReport({}).empty());
}
