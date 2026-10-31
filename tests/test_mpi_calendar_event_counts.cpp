#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#ifdef USE_MPI
#include <mpi.h>

#include <cstdint>
#include <vector>

#include "core/world_state.h"
#include "epidemiology/calendar_event.h"
#include "test_utils.h"

using namespace june;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();
  MPI_Finalize();
  return result;
}

namespace {

constexpr int kPeople = 300;

// Builds the world holding only the people this rank owns (person i goes to
// rank i % size), so the ranks together hold exactly the one-rank world.
WorldState buildShare(int rank, int size) {
  WorldState world = makeCatchmentWorld(0);
  for (int i = 0; i < kPeople; ++i) {
    if (i % size != rank) continue;
    Person& p = world.people.emplace_back();
    p.id = i;
    p.geo_unit_id = 0;
    p.age = (i % 3 == 0) ? 10.0f : 50.0f;
    if (i % 17 == 1) p.schedule_hop = ScheduleHop::begin(1, -1, 0);
  }
  world.buildIndices();
  return world;
}

CalendarEventManager::EventDayCounts countsFor(int rank, int size) {
  WorldState world = buildShare(rank, size);
  SelectionCriterion age_filter;
  age_filter.property_path = "age";
  age_filter.operator_type = ">=";
  age_filter.value = 20;
  CalendarEvent event;
  event.calendar_event_id = 1;
  event.start_day = 0;
  event.schedule_type_idx = 1;
  event.compliance_rate = 0.6f;
  event.catchment_rule_id = 0;
  event.attendee_filters = {age_filter};
  CalendarEventManager manager({{event}});
  manager.triggerEventsForDay(0, world, world.people, 77, {{0, {0}}});
  return manager.lastDayCounts().at(0);
}

}  // namespace

// Summed over ranks, the per-event counts equal what one rank holding every
// person counts, whatever the rank count.
TEST_CASE("Calendar event counts summed over ranks equal the one-rank counts") {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const auto mine = countsFor(rank, size);
  const uint64_t local[5] = {mine.in_catchment, mine.passed_filters,
                             mine.already_hopping, mine.refused, mine.hopping};
  uint64_t global[5] = {};
  MPI_Allreduce(local, global, 5, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

  const auto whole = countsFor(0, 1);
  CHECK(global[0] == whole.in_catchment);
  CHECK(global[1] == whole.passed_filters);
  CHECK(global[2] == whole.already_hopping);
  CHECK(global[3] == whole.refused);
  CHECK(global[4] == whole.hopping);
  CHECK(whole.hopping > 0);
  CHECK(whole.refused > 0);
  CHECK(whole.already_hopping > 0);
}

#else
TEST_CASE("Calendar event counts MPI test requires MPI") {}
#endif
