#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "doctest.h"
#include "epidemiology/vaccination_manager.h"
#include "utils/event_logging/event_logger.h"
#include "vaccination_manager_fixtures.h"

using namespace june;
using namespace june::manager_fixture;

namespace {

// Doses given on each of the first `days` days.
std::vector<int64_t> dosesPerDay(VaccinationManager& manager,
                                 const WorldState& world, int days) {
  std::vector<int64_t> per_day;
  for (int day = 0; day < days; ++day) {
    manager.update(day);
    per_day.push_back(static_cast<int64_t>(dosesGivenOn(world, day).size()));
  }
  return per_day;
}

std::set<PersonId> idsOf(const DosesGiven& given) {
  std::set<PersonId> ids;
  for (const auto& [id, dose_index] : given) ids.insert(id);
  return ids;
}

// The `take` people of `pool` with the lowest lottery keys in tier `tier` of
// campaign 0.
std::set<PersonId> lowestKeys(const std::set<PersonId>& pool, size_t tier,
                              size_t take) {
  std::vector<std::pair<uint64_t, PersonId>> keyed;
  for (PersonId id : pool)
    keyed.emplace_back(lotteryKey(kSeed, static_cast<uint64_t>(id), 0, tier),
                       id);
  std::sort(keyed.begin(), keyed.end());
  std::set<PersonId> lowest;
  for (size_t i = 0; i < take && i < keyed.size(); ++i)
    lowest.insert(keyed[i].second);
  return lowest;
}

std::set<PersonId> withAge(const WorldState& world, bool child) {
  std::set<PersonId> ids;
  for (const Person& person : world.people)
    if ((person.age < 18) == child) ids.insert(person.id);
  return ids;
}

std::set<PersonId> minus(std::set<PersonId> pool,
                         const std::set<PersonId>& taken) {
  for (PersonId id : taken) pool.erase(id);
  return pool;
}

std::set<PersonId> plus(std::set<PersonId> a, const std::set<PersonId>& b) {
  a.insert(b.begin(), b.end());
  return a;
}

}  // namespace

// =============================================================================
// The cap on first doses
// =============================================================================

TEST_CASE("firstDoseCap is the exact floor of target times population") {
  // 0.29 is not exactly 29/100: 0.29 * 100 is 28.999999999999996 in a double,
  // and a plain floor would be one person short.
  CHECK(std::floor(0.29 * 100.0) == 28.0);
  CHECK(firstDoseCap(0.29, 100) == 29);
  CHECK(firstDoseCap(0.7, 60000000) == 42000000);
  CHECK(firstDoseCap(0.7, 630000) == 441000);
  CHECK(firstDoseCap(1.0, 12345) == 12345);
  CHECK(firstDoseCap(0.5, 5) == 2);

  // A target written with two decimals is a whole number of hundredths, so the
  // exact answer is integer arithmetic.
  for (int hundredths = 1; hundredths <= 100; ++hundredths) {
    const double target = hundredths / 100.0;
    for (int64_t population :
         {int64_t{1}, int64_t{7}, int64_t{99}, int64_t{100}, int64_t{630000},
          int64_t{60000000}, int64_t{60000001}, int64_t{1234567890}}) {
      CAPTURE(hundredths);
      CAPTURE(population);
      CHECK(firstDoseCap(target, population) == hundredths * population / 100);
    }
  }

  // A real fraction is never rounded up to the next person.
  CHECK(firstDoseCap(0.7, 3) == 2);
  CHECK(firstDoseCap(0.6999999, 10000001) == 6999999);
  CHECK(firstDoseCap(0.5, 1) == 0);
}

TEST_CASE("firstDoseCap names a target or population it cannot use") {
  CHECK_THROWS_WITH_AS(firstDoseCap(0.0, 10), doctest::Contains("coverage"),
                       std::runtime_error);
  CHECK_THROWS_AS(firstDoseCap(1.5, 10), std::runtime_error);
  CHECK_THROWS_AS(firstDoseCap(std::nan(""), 10), std::runtime_error);
  CHECK_THROWS_WITH_AS(firstDoseCap(0.5, -1), doctest::Contains("population"),
                       std::runtime_error);
}

// =============================================================================
// Campaign window and campaign day
// =============================================================================

TEST_CASE("a campaign gives doses only from its start to its end date") {
  auto campaign = vaxCampaign("window", {tierFor("all", 0)});
  campaign.start_date = "2024-01-03";  // run day 2
  campaign.end_date = "2024-01-05";    // run day 4, included
  campaign.dose_budget = DoseBudget{5, 1};
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(homesOfThree(20));
  VaccinationManager manager(world, config);

  CHECK(dosesPerDay(manager, world, 8) ==
        std::vector<int64_t>{0, 0, 5, 5, 5, 0, 0, 0});
}

TEST_CASE("campaign day 0 is the run's day 0 when there is no start date") {
  auto campaign = vaxCampaign("open", {tierFor("all", 0)});
  campaign.dose_budget = DoseBudget{5, 3};
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(homesOfThree(20));
  VaccinationManager manager(world, config);

  // The allotments of campaign days 0 to 4 of 5 doses per 3 days
  CHECK(dosesPerDay(manager, world, 5) == std::vector<int64_t>{1, 2, 2, 1, 2});
}

TEST_CASE("a start date before the run starts shifts the campaign day") {
  auto campaign = vaxCampaign("running", {tierFor("all", 0)});
  campaign.start_date = "2023-12-31";  // run day -1, so run day 0 is day 1
  campaign.dose_budget = DoseBudget{5, 3};
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(homesOfThree(20));
  VaccinationManager manager(world, config);

  std::vector<int64_t> expected;
  for (int64_t run_day = 0; run_day < 5; ++run_day)
    expected.push_back(doseAllotment(5, 3, run_day + 1));
  CHECK(expected == std::vector<int64_t>{2, 2, 1, 2, 2});
  CHECK(dosesPerDay(manager, world, 5) == expected);
}

// =============================================================================
// The dose budget
// =============================================================================

TEST_CASE("with plenty of candidates a day gives its allotment, no more") {
  auto campaign = vaxCampaign("budget", {tierFor("all", 0)});
  campaign.dose_budget = DoseBudget{7, 3};
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(homesOfThree(30));
  VaccinationManager manager(world, config);

  const auto per_day = dosesPerDay(manager, world, 6);
  for (int64_t day = 0; day < 6; ++day)
    CHECK(per_day[static_cast<size_t>(day)] == doseAllotment(7, 3, day));
  CHECK(per_day == std::vector<int64_t>{2, 2, 3, 2, 2, 3});
  // Two whole periods give exactly twice the doses
  CHECK(std::accumulate(per_day.begin(), per_day.end(), int64_t{0}) == 14);
}

TEST_CASE("doses a day cannot use are lost, not added to the next day") {
  // Six people had a first dose on day -1 and are due their second on day 1.
  // Two have had none. Day 0 can use 2 of its 5 doses.
  auto campaign =
      vaxCampaign("lost", {tierFor("second", 1), tierFor("first", 0)});
  campaign.dose_budget = DoseBudget{5, 1};
  campaign.days_to_next_dose = {2};
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(std::vector<Home>(4, Home{{30, 40}}));
  for (size_t i = 0; i < 6; ++i) giveDose(world.people[i], "TestVax", -1.0);
  VaccinationManager manager(world, config);

  const auto per_day = dosesPerDay(manager, world, 2);
  CHECK(per_day[0] == 2);
  // Banked, day 1 would give 5 + 3 = 8 of the 6 + 2 who are due
  CHECK(per_day[1] == 5);
}

// =============================================================================
// Strict priority and the lottery
// =============================================================================

TEST_CASE("a lower tier waits until the tiers above it have no one left") {
  // 10 children and 20 adults, 4 doses a day, children first
  Config config = makeConfig({vaxCampaign(
      "priority", {tierFor("children", 0, {criterion("age", "<", 18.0)}),
                   tierFor("adults", 0, {criterion("age", ">=", 18.0)})})});
  config.vaccination.campaigns[0].dose_budget = DoseBudget{4, 1};
  WorldState world = makeWorld(homesOfThree(10));
  VaccinationManager manager(world, config);

  const std::set<PersonId> children = withAge(world, true);
  const std::set<PersonId> adults = withAge(world, false);
  REQUIRE(children.size() == 10);
  REQUIRE(adults.size() == 20);

  // Day 0 and 1: only children, the four lowest keys each day
  manager.update(0);
  const auto day0 = idsOf(dosesGivenOn(world, 0));
  CHECK(day0 == lowestKeys(children, 0, 4));
  manager.update(1);
  const auto day1 = idsOf(dosesGivenOn(world, 1));
  CHECK(day1 == lowestKeys(minus(children, day0), 0, 4));

  // Day 2: the last two children, then the first two adults by their own keys
  manager.update(2);
  const auto day2 = idsOf(dosesGivenOn(world, 2));
  const auto left = minus(children, plus(day0, day1));
  REQUIRE(left.size() == 2);
  CHECK(day2 == plus(left, lowestKeys(adults, 1, 2)));

  // Day 3: four adults, never a repeat
  manager.update(3);
  const auto day3 = idsOf(dosesGivenOn(world, 3));
  CHECK(day3 == lowestKeys(minus(adults, day2), 1, 4));
}

TEST_CASE("the same seed serves the same people and another seed does not") {
  const auto serve = [](uint64_t seed) {
    Config config = makeConfig({vaxCampaign("seeded", {tierFor("all", 0)})});
    config.simulation.random_seed = seed;
    config.vaccination.campaigns[0].dose_budget = DoseBudget{10, 1};
    WorldState world = makeWorld(homesOfThree(30));
    VaccinationManager manager(world, config);
    manager.update(0);
    return idsOf(dosesGivenOn(world, 0));
  };
  CHECK(serve(1) == serve(1));
  CHECK(serve(1) != serve(2));
}

// =============================================================================
// Dose order and the coverage cap, with the daily summary rows
// =============================================================================

namespace {

struct CapRun {
  std::vector<DosesGiven> days;
  std::vector<CampaignDayEvent> rows;
  size_t vaccination_events = 0;
  int64_t total_doses = 0;
};

// 60 people, 10 doses a day, second doses ahead of first doses, a gap of two
// days and a cap of half the population (30 first doses). The campaign ends on
// run day 7.
//   day:    0   1   2   3   4   5   6   7
//   first: 10  10   0   0  10   0   0   0
//   second: 0   0  10  10   0   0  10   0
CapRun runCapScenario() {
  auto campaign =
      vaxCampaign("cap", {tierFor("second", 1), tierFor("first", 0)});
  campaign.dose_budget = DoseBudget{10, 1};
  campaign.coverage_target = 0.5;
  campaign.days_to_next_dose = {2};
  campaign.end_date = "2024-01-08";
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(homesOfThree(20));
  EventLogger logger;
  VaccinationManager manager(world, config, &logger);

  CapRun run;
  for (int day = 0; day < 10; ++day) {
    manager.update(day);
    run.days.push_back(dosesGivenOn(world, day));
  }
  run.rows = logger.getCampaignDayEvents();
  run.vaccination_events = logger.getVaccinationCount();
  run.total_doses = totalDoses(world);
  return run;
}

CampaignDayEvent rowOf(const CapRun& run, int day, uint16_t tier) {
  for (const auto& row : run.rows)
    if (row.day == day && row.tier_id == tier) return row;
  FAIL("no row for day " << day << " tier " << tier);
  return {};
}

}  // namespace

TEST_CASE("due later doses go first, and the cap stops first doses only") {
  const CapRun run = runCapScenario();
  const std::vector<int64_t> first = {10, 10, 0, 0, 10, 0, 0, 0, 0, 0};
  const std::vector<int64_t> second = {0, 0, 10, 10, 0, 0, 10, 0, 0, 0};
  for (size_t day = 0; day < run.days.size(); ++day) {
    CAPTURE(day);
    CHECK(dosesOfIndex(run.days[day], 0) == first[day]);
    CHECK(dosesOfIndex(run.days[day], 1) == second[day]);
  }
  // First doses stop exactly at floor(0.5 * 60) and second doses carry on
  // after the cap is reached on day 4 (day 6 gives ten).
  int64_t first_total = 0;
  for (const auto& day : run.days) first_total += dosesOfIndex(day, 0);
  CHECK(first_total == 30);
  CHECK(run.total_doses == 60);
}

TEST_CASE("every dose is logged once, and nothing is dosed after the end") {
  const CapRun run = runCapScenario();
  CHECK(run.vaccination_events == static_cast<size_t>(run.total_doses));
  CHECK(run.days[8].empty());
  CHECK(run.days[9].empty());
}

TEST_CASE("the daily summary has a row per tier on every active day") {
  const CapRun run = runCapScenario();
  // Run days 0 to 7, two tiers: zero rows on quiet days are rows too
  REQUIRE(run.rows.size() == 16);
  for (const auto& row : run.rows) {
    CHECK(row.campaign_id == 0);
    CHECK(row.day >= 0);
    CHECK(row.day <= 7);
    CHECK(row.allotment == 10);
  }

  // Day 0: nobody has a first dose yet, 60 could get one, 10 do
  auto row = rowOf(run, 0, 1);
  CHECK(row.candidates == 60);
  CHECK(row.doses_given == 10);
  CHECK(row.declined == 0);
  CHECK(row.started_total == 10);
  row = rowOf(run, 0, 0);
  CHECK(row.candidates == 0);
  CHECK(row.doses_given == 0);

  // Day 2: ten second doses use the whole budget, so no first doses
  row = rowOf(run, 2, 0);
  CHECK(row.candidates == 10);
  CHECK(row.doses_given == 10);
  CHECK(row.started_total == 20);
  row = rowOf(run, 2, 1);
  CHECK(row.candidates == 40);
  CHECK(row.doses_given == 0);
  CHECK(row.started_total == 20);

  // Day 4 reaches the cap, day 5 has candidates but the cap gives them nothing
  row = rowOf(run, 4, 1);
  CHECK(row.doses_given == 10);
  CHECK(row.started_total == 30);
  row = rowOf(run, 5, 1);
  CHECK(row.candidates == 30);
  CHECK(row.doses_given == 0);
  CHECK(row.started_total == 30);

  // The last active day has nothing to give but still has its rows
  row = rowOf(run, 7, 0);
  CHECK(row.candidates == 0);
  CHECK(row.doses_given == 0);
  CHECK(row.started_total == 30);
}

TEST_CASE("a campaign with no budget logs no allotment") {
  Config config = makeConfig({vaxCampaign("unlimited", {tierFor("all", 0)})});
  WorldState world = makeWorld(homesOfThree(5));
  EventLogger logger;
  VaccinationManager manager(world, config, &logger);
  manager.update(0);
  manager.update(1);

  const auto& rows = logger.getCampaignDayEvents();
  REQUIRE(rows.size() == 2);
  for (const auto& row : rows) CHECK(row.allotment == kNoAllotment);
  // Everyone is served on the first day
  CHECK(rows[0].candidates == 15);
  CHECK(rows[0].doses_given == 15);
  CHECK(rows[0].started_total == 15);
  CHECK(rows[1].candidates == 0);
  CHECK(rows[1].started_total == 15);
}

// =============================================================================
// The report at the end of a run
// =============================================================================

namespace {

// What the manager says about `campaign` after it has run for `days` days on a
// world of 15 people.
std::string coverageReport(VaccinationCampaignConfig campaign, int days) {
  Config config = makeConfig({std::move(campaign)});
  WorldState world = makeWorld(homesOfThree(5));
  VaccinationManager manager(world, config);
  for (int day = 0; day < days; ++day) manager.update(day);
  std::ostringstream report;
  manager.printCoverage(report);
  return report.str();
}

}  // namespace

TEST_CASE(
    "the end of run report says whether the coverage target was reached") {
  SUBCASE("a target that was reached") {
    // 0.4 of 15 people allows 6 first doses, and the first day gives them
    auto campaign = vaxCampaign("reached", {tierFor("all", 0)});
    campaign.coverage_target = 0.4;
    CHECK(coverageReport(campaign, 2) ==
          "vaccination campaign 'reached': 6 of 6 first doses the coverage "
          "target allows (100.0%), target reached\n");
  }

  SUBCASE("a target the budget could not reach in time") {
    // 3 doses a day for 2 days is 6 of the 15 the target allows
    auto campaign = vaxCampaign("short", {tierFor("all", 0)});
    campaign.dose_budget = DoseBudget{3, 1};
    campaign.coverage_target = 1.0;
    campaign.end_date = "2024-01-02";
    CHECK(coverageReport(campaign, 5) ==
          "vaccination campaign 'short': 6 of 15 first doses the coverage "
          "target allows (40.0%), target not reached\n");
  }

  SUBCASE("a target that nobody can reach because everyone declines") {
    VaccinationTier all = tierFor("all", 0);
    all.uptake.rows = {row({}, 0.0)};
    auto campaign = vaxCampaign("declined", {all});
    campaign.coverage_target = 0.5;
    CHECK(coverageReport(campaign, 3) ==
          "vaccination campaign 'declined': 0 of 7 first doses the coverage "
          "target allows (0.0%), target not reached\n");
  }

  SUBCASE("a campaign with no target reports who started") {
    CHECK(coverageReport(vaxCampaign("open", {tierFor("all", 0)}), 2) ==
          "vaccination campaign 'open': 15 people started a course (no "
          "coverage target)\n");
  }

  SUBCASE("a disabled vaccination config says nothing") {
    Config config = makeConfig({vaxCampaign("off", {tierFor("all", 0)})});
    config.vaccination.enabled = false;
    WorldState world = makeWorld(homesOfThree(5));
    VaccinationManager manager(world, config);
    std::ostringstream report;
    manager.printCoverage(report);
    CHECK(report.str().empty());
  }
}

// =============================================================================
// Eligibility rules a campaign keeps
// =============================================================================

TEST_CASE("a later dose waits for the gap since the last one") {
  auto campaign =
      vaxCampaign("gap", {tierFor("second", 1), tierFor("first", 0)});
  campaign.days_to_next_dose = {5};
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(homesOfThree(2));
  VaccinationManager manager(world, config);

  CHECK(dosesPerDay(manager, world, 7) ==
        std::vector<int64_t>{6, 0, 0, 0, 0, 6, 0});
}

TEST_CASE("a later dose follows the brand filter") {
  auto campaign = vaxCampaign("filter", {tierFor("second", 1)});
  campaign.courses_begun_with = {"TestVax"};
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(homesOfThree(2));
  for (size_t i = 0; i < world.people.size(); ++i)
    giveDose(world.people[i], i < 3 ? "TestVax" : "Other", -10.0);
  VaccinationManager manager(world, config);

  manager.update(0);
  CHECK(idsOf(dosesGivenOn(world, 0)) == std::set<PersonId>{0, 1, 2});
  manager.update(1);
  CHECK(dosesGivenOn(world, 1).empty());
}

TEST_CASE("with no brand filter a course goes on with the vaccine it began") {
  Config config =
      makeConfig({vaxCampaign("same brand", {tierFor("second", 1)})});
  WorldState world = makeWorld(homesOfThree(2));
  // Three people began with the campaign's vaccine and three with another
  for (size_t i = 0; i < world.people.size(); ++i)
    giveDose(world.people[i], i < 3 ? "TestVax" : "Other", -10.0);
  VaccinationManager manager(world, config);

  manager.update(0);
  CHECK(idsOf(dosesGivenOn(world, 0)) == std::set<PersonId>{0, 1, 2});
  manager.update(1);
  CHECK(dosesGivenOn(world, 1).empty());
}

TEST_CASE("a brand filter lets a booster follow another vaccine") {
  auto campaign = vaxCampaign("booster", {tierFor("third", 1)});
  campaign.courses_begun_with = {"Other"};
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(homesOfThree(2));
  for (size_t i = 0; i < world.people.size(); ++i)
    giveDose(world.people[i], i < 3 ? "TestVax" : "Other", -10.0);
  VaccinationManager manager(world, config);

  // Only the filter decides: the campaign's own vaccine is not accepted unless
  // it is listed
  manager.update(0);
  CHECK(idsOf(dosesGivenOn(world, 0)) == std::set<PersonId>{3, 4, 5});
}

TEST_CASE(
    "a decliner is never dosed and does not fall through to a later tier") {
  // Children decline for sure; the second tier would take them if they fell
  // through, since it covers everyone.
  VaccinationTier children =
      tierFor("children", 0, {criterion("age", "<", 18.0)});
  children.uptake.rows = {row({}, 0.0)};
  Config config =
      makeConfig({vaxCampaign("decline", {children, tierFor("everyone", 0)})});
  WorldState world = makeWorld(homesOfThree(10));
  EventLogger logger;
  VaccinationManager manager(world, config, &logger);

  for (int day = 0; day < 3; ++day) manager.update(day);
  const auto kids = withAge(world, true);
  const auto grown = withAge(world, false);
  std::set<PersonId> dosed;
  for (const auto& person : world.people)
    if (person.vaccine_trajectory) dosed.insert(person.id);
  CHECK(dosed == grown);
  for (PersonId id : kids) CHECK(dosed.count(id) == 0);

  // The summary counts them as waiting decliners every day, not candidates
  const auto& rows = logger.getCampaignDayEvents();
  REQUIRE(rows.size() == 6);
  for (int day = 0; day < 3; ++day) {
    CHECK(rows[2 * day].candidates == 0);
    CHECK(rows[2 * day].declined == 10);
    CHECK(rows[2 * day + 1].candidates == (day == 0 ? 20 : 0));
    CHECK(rows[2 * day + 1].declined == 0);
  }
}

TEST_CASE("a person gets one dose a day across campaigns, from the first") {
  SUBCASE("a later dose of a second campaign waits for the next day") {
    auto second = vaxCampaign("second", {tierFor("second", 1)});
    second.days_to_next_dose = {0};
    Config config =
        makeConfig({vaxCampaign("first", {tierFor("first", 0)}), second});
    WorldState world = makeWorld(homesOfThree(4));
    VaccinationManager manager(world, config);
    // 12 first doses on day 0 and, with no gap, 12 second doses on day 1
    CHECK(dosesPerDay(manager, world, 3) == std::vector<int64_t>{12, 12, 0});
    CHECK(dosesOfIndex(dosesGivenOn(world, 1), 1) == 12);
  }
  SUBCASE("the first campaign in file order claims a shared person") {
    Config config = makeConfig({vaxCampaign("a", {tierFor("t", 0)}, "TestVax"),
                                vaxCampaign("b", {tierFor("t", 0)}, "Other")});
    WorldState world = makeWorld(homesOfThree(4));
    EventLogger logger;
    VaccinationManager manager(world, config, &logger);
    manager.update(0);
    for (const auto& person : world.people) {
      REQUIRE(person.vaccine_trajectory);
      CHECK(person.vaccine_trajectory->vaccine_name == "TestVax");
      CHECK(person.vaccine_trajectory->getNumDoses() == 1);
    }
    // The second campaign saw nobody it could still dose today
    const auto& rows = logger.getCampaignDayEvents();
    REQUIRE(rows.size() == 2);
    CHECK(rows[0].campaign_id == 0);
    CHECK(rows[0].doses_given == 12);
    CHECK(rows[1].campaign_id == 1);
    CHECK(rows[1].candidates == 0);
    CHECK(rows[1].doses_given == 0);
  }
}

TEST_CASE("members of a household unit get one outcome") {
  VaccinationTier children =
      tierFor("children", 0, {criterion("age", "<", 18.0)});
  shareResidence(children, DeciderCombine::Mean);
  children.uptake.rows = {row({}, 0.5)};
  const auto campaign = vaxCampaign("household", {children});
  Config config = makeConfig({campaign});
  WorldState world = makeWorld(std::vector<Home>(60, Home{{30, 40, 5, 8}}));
  VaccinationManager manager(world, config);
  manager.update(0);

  const auto dosed = idsOf(dosesGivenOn(world, 0));
  const CampaignUptake uptake = settleUptakeLocal(world, campaign, 0, kSeed);
  size_t homes_dosed = 0;
  for (size_t home = 0; home < 60; ++home) {
    // Ids 4h + 2 and 4h + 3 are the home's two children
    const PersonId a = static_cast<PersonId>(4 * home + 2);
    const PersonId b = a + 1;
    CHECK(dosed.count(a) == dosed.count(b));
    homes_dosed += dosed.count(a);
  }
  CHECK(homes_dosed > 0);
  CHECK(homes_dosed < 60);
  // Exactly the people the uptake pass said accept
  std::set<PersonId> accepted;
  for (size_t i = 0; i < world.people.size(); ++i)
    if (uptake.accepted[i] & 1) accepted.insert(world.people[i].id);
  CHECK(dosed == accepted);
}

// =============================================================================
// Resume
// =============================================================================

TEST_CASE("a manager built before a restore carries on as the original") {
  auto campaign =
      vaxCampaign("resume", {tierFor("second", 1), tierFor("first", 0)});
  campaign.dose_budget = DoseBudget{10, 1};
  campaign.coverage_target = 0.5;
  campaign.days_to_next_dose = {2};
  Config config = makeConfig({campaign});
  const auto homes = homesOfThree(20);
  WorldState original = makeWorld(homes);
  WorldState resumed = makeWorld(homes);
  VaccinationManager uninterrupted(original, config);
  // The simulator builds its manager before a checkpoint is restored
  VaccinationManager restored(resumed, config);

  // The cap is reached on day 4, so a resumed manager that forgot who had
  // started would hand out first doses again on day 5.
  for (int day = 0; day <= 4; ++day) uninterrupted.update(day);
  for (size_t i = 0; i < original.people.size(); ++i)
    if (original.people[i].vaccine_trajectory)
      resumed.people[i].vaccine_trajectory =
          std::make_unique<VaccineTrajectory>(
              *original.people[i].vaccine_trajectory);
  REQUIRE(totalDoses(resumed) == totalDoses(original));

  int64_t first_doses_after = 0;
  for (int day = 5; day < 12; ++day) {
    CAPTURE(day);
    uninterrupted.update(day);
    restored.update(day);
    const auto expected = dosesGivenOn(original, day);
    CHECK(dosesGivenOn(resumed, day) == expected);
    first_doses_after += dosesOfIndex(expected, 0);
  }
  CHECK(first_doses_after == 0);
  CHECK(totalDoses(resumed) == totalDoses(original));
  CHECK(totalDoses(original) == 60);
}

// =============================================================================
// Disabled vaccination and bad configs
// =============================================================================

TEST_CASE("a disabled vaccination config does no uptake pass") {
  // More tiers than an uptake pass can hold: it throws if it runs at all
  std::vector<VaccinationTier> many(kMaxUptakeTiers + 1, tierFor("t", 0));
  Config config = makeConfig({vaxCampaign("big", many)});
  WorldState world = makeWorld(homesOfThree(2));

  config.vaccination.enabled = false;
  EventLogger logger;
  CHECK_NOTHROW(VaccinationManager(world, config, &logger).update(0));
  CHECK(logger.getTotalRecordCount() == 0);
  CHECK(totalDoses(world) == 0);

  config.vaccination.enabled = true;
  CHECK_THROWS_WITH_AS(VaccinationManager(world, config),
                       doctest::Contains("65 tiers"), std::runtime_error);
}

TEST_CASE("a campaign naming a missing vaccine or dose is refused by name") {
  WorldState world = makeWorld(homesOfThree(2));
  CHECK_THROWS_WITH_AS(
      VaccinationManager(
          world, makeConfig({vaxCampaign("lost", {tierFor("t", 0)}, "Nope")})),
      doctest::Contains("lost"), std::runtime_error);
  CHECK_THROWS_WITH_AS(
      VaccinationManager(
          world, makeConfig({vaxCampaign("lost", {tierFor("t", 0)}, "Nope")})),
      doctest::Contains("Nope"), std::runtime_error);
  // TestVax has three doses, numbered 0 to 2
  CHECK_THROWS_WITH_AS(
      VaccinationManager(world,
                         makeConfig({vaxCampaign("high", {tierFor("t3", 3)})})),
      doctest::Contains("tier 't3'"), std::runtime_error);
}
