#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#ifdef USE_MPI
#include <mpi.h>

#include <algorithm>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "epidemiology/vaccination_manager.h"
#include "utils/event_logging/event_logger.h"
#include "utils/mpi_utils.h"
#include "vaccination_manager_fixtures.h"

using namespace june;
using namespace june::manager_fixture;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();
  MPI_Finalize();
  return result;
}

namespace {

constexpr int kDays = 16;

int rankOf() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int sizeOf() {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

// This rank's homes only, as in the real distributed load
WorldState myWorld(const std::vector<Home>& homes) {
  const size_t rank = static_cast<size_t>(rankOf());
  const size_t size = static_cast<size_t>(sizeOf());
  return makeWorld(homes, [&](size_t h) { return h % size == rank; });
}

// Second doses first, then first doses, children ahead of adults, households
// deciding for their children. The budget is tight enough that tiers are
// partly served, and the target is reached part way through.
VaccinationCampaignConfig rollout() {
  const auto child = criterion("age", "<", 18.0);
  const auto grown = criterion("age", ">=", 18.0);
  VaccinationTier kids_second = tierFor("kids_second", 1, {child});
  shareResidence(kids_second, DeciderCombine::Mean);
  kids_second.uptake.rows = {row({}, 0.9)};
  VaccinationTier kids_first = tierFor("kids_first", 0, {child});
  shareResidence(kids_first, DeciderCombine::Mean);
  kids_first.uptake.rows = {row({criterion("age", "<", 35.0)}, 0.8),
                            row({}, 0.4)};
  VaccinationTier adults_second = tierFor("adults_second", 1, {grown});
  VaccinationTier adults_first = tierFor("adults_first", 0, {grown});
  adults_first.uptake.rows = {row({criterion("age", "<", 45.0)}, 0.6),
                              row({}, 0.3)};
  auto campaign = vaxCampaign(
      "rollout", {kids_second, kids_first, adults_second, adults_first});
  campaign.dose_budget = DoseBudget{200, 5};
  campaign.coverage_target = 0.3;
  campaign.days_to_next_dose = {3};
  return campaign;
}

struct Reference {
  std::vector<DosesGiven> days;
  std::vector<std::vector<CampaignDayEvent>> rows;  // per day
  int64_t first_doses = 0;
  int64_t cap = 0;
};

// The same rollout on the whole world, written out plainly with no MPI: sort
// every candidate of a tier by key and serve the lowest.
Reference serialReference(const std::vector<Home>& homes,
                          const VaccinationCampaignConfig& campaign) {
  WorldState world = makeWorld(homes);
  const CampaignUptake uptake = settleUptakeLocal(world, campaign, 0, kSeed);
  const size_t tiers = campaign.tiers.size();
  std::vector<bool> starts_course;
  for (const auto& tier : campaign.tiers)
    starts_course.push_back(tier.dose == 0);
  Reference ref;
  ref.cap = firstDoseCap(*campaign.coverage_target,
                         static_cast<int64_t>(world.people.size()));
  for (int day = 0; day < kDays; ++day) {
    int64_t started = 0;
    std::vector<std::vector<std::pair<uint64_t, size_t>>> candidates(tiers);
    std::vector<int64_t> declined(tiers, 0);
    for (size_t i = 0; i < world.people.size(); ++i) {
      const Person& person = world.people[i];
      const auto* trajectory = person.vaccine_trajectory.get();
      const int held = trajectory ? trajectory->getNumDoses() : 0;
      if (held > 0) ++started;
      // A gap the campaign does not list means no wait, as in the manager
      const size_t gap = static_cast<size_t>(std::max(held - 1, 0));
      const double wait = gap < campaign.days_to_next_dose.size()
                              ? campaign.days_to_next_dose[gap]
                              : 0.0;
      if (held > 0 &&
          (trajectory->getLastDoseTime() == day ||
           day < trajectory->getLastDoseTime() + wait))
        continue;
      for (size_t t = 0; t < tiers; ++t) {
        if (!((uptake.selected[i] >> t) & 1) || campaign.tiers[t].dose != held)
          continue;
        if ((uptake.accepted[i] >> t) & 1)
          candidates[t].emplace_back(
              lotteryKey(kSeed, static_cast<uint64_t>(person.id), 0, t), i);
        else
          ++declined[t];
        break;
      }
    }
    const int64_t allotment = doseAllotment(
        campaign.dose_budget->doses, campaign.dose_budget->per_days, day);
    std::vector<int64_t> counts;
    for (const auto& list : candidates)
      counts.push_back(static_cast<int64_t>(list.size()));
    const auto take = fillTiers(allotment, counts, starts_course,
                                std::max<int64_t>(0, ref.cap - started));
    int64_t started_after = started;
    ref.rows.emplace_back();
    for (size_t t = 0; t < tiers; ++t) {
      std::sort(candidates[t].begin(), candidates[t].end());
      for (int64_t k = 0; k < take[t]; ++k)
        giveDose(world.people[candidates[t][static_cast<size_t>(k)].second],
                 campaign.vaccine_type, day);
      if (starts_course[t]) started_after += take[t];
    }
    for (size_t t = 0; t < tiers; ++t)
      ref.rows.back().push_back({day, 0, static_cast<uint16_t>(t), counts[t],
                                 take[t], declined[t], allotment,
                                 started_after});
    ref.days.push_back(dosesGivenOn(world, day));
    ref.first_doses += dosesOfIndex(ref.days.back(), 0);
  }
  return ref;
}

bool same(const CampaignDayEvent& a, const CampaignDayEvent& b) {
  return a.day == b.day && a.campaign_id == b.campaign_id &&
         a.tier_id == b.tier_id && a.candidates == b.candidates &&
         a.doses_given == b.doses_given && a.declined == b.declined &&
         a.allotment == b.allotment && a.started_total == b.started_total;
}

}  // namespace

// Every rank holds only its own homes. Each day's doses, person by person, and
// the summary rows must be what one process gives the whole world.
TEST_CASE("the manager on ranks holding their own homes matches a serial run") {
  const auto homes = mixedHomes(300);
  const auto campaign = rollout();
  Config config = makeConfig({campaign});
  const Reference ref = serialReference(homes, campaign);
  // The scenario would prove little if the cap were never reached or no tier
  // were partly served
  REQUIRE(ref.first_doses == ref.cap);
  REQUIRE(ref.days[0].size() > 0);

  WorldState mine = myWorld(homes);
  std::set<PersonId> my_ids;
  for (const Person& person : mine.people) my_ids.insert(person.id);
  EventLogger logger;
  VaccinationManager manager(mine, config, &logger);

  int64_t my_doses = 0;
  for (int day = 0; day < kDays; ++day) {
    CAPTURE(day);
    manager.update(day);
    const DosesGiven local = dosesGivenOn(mine, day);
    DosesGiven expected;
    for (const auto& [id, dose_index] : ref.days[static_cast<size_t>(day)])
      if (my_ids.count(id)) expected[id] = dose_index;
    CHECK(local == expected);
    CHECK(mpi_runtime::globalSum(static_cast<int64_t>(local.size())) ==
          static_cast<int64_t>(ref.days[static_cast<size_t>(day)].size()));
    my_doses += static_cast<int64_t>(local.size());
  }

  // The end of run report counts over every rank and only rank 0 writes it
  std::ostringstream report;
  manager.printCoverage(report);
  if (rankOf() == 0) {
    CHECK(report.str() ==
          "vaccination campaign 'rollout': " + std::to_string(ref.first_doses) +
              " of " + std::to_string(ref.cap) +
              " first doses the coverage target allows (100.0%), target "
              "reached\n");
  } else {
    CHECK(report.str().empty());
  }

  // Each rank logs its own people's doses, and only rank 0 logs the summary
  CHECK(logger.getVaccinationCount() == static_cast<size_t>(my_doses));
  const auto& rows = logger.getCampaignDayEvents();
  if (rankOf() == 0) {
    std::vector<CampaignDayEvent> expected;
    for (const auto& day_rows : ref.rows)
      expected.insert(expected.end(), day_rows.begin(), day_rows.end());
    REQUIRE(rows.size() == expected.size());
    for (size_t i = 0; i < rows.size(); ++i) {
      CAPTURE(i);
      CHECK(same(rows[i], expected[i]));
    }
  } else {
    CHECK(rows.empty());
  }
}

// A disabled manager must not reduce: if rank 0 did, it would wait here for
// ranks that never call.
TEST_CASE("a disabled vaccination config makes no collective call") {
  Config config = makeConfig({rollout()});
  config.vaccination.enabled = false;
  WorldState mine = myWorld(mixedHomes(30));
  VaccinationManager manager(mine, config);
  if (rankOf() == 0)
    for (int day = 0; day < 3; ++day) manager.update(day);
  MPI_Barrier(MPI_COMM_WORLD);
  CHECK(totalDoses(mine) == 0);
}

#else
TEST_CASE("MPI vaccination manager requires MPI") {}
#endif
