#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#ifdef USE_MPI
#include <mpi.h>

#include <cstdint>
#include <memory>
#include <vector>

#include "core/world_state.h"
#include "epidemiology/applicable_policies/applicable_policy_walk.h"
#include "epidemiology/policy.h"
#include "mpi_test_helpers.h"
#include "policy_test_helpers.h"
#include "venue_scope_fixtures.h"

using namespace june;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();
  MPI_Finalize();
  return result;
}

// The same deaths are split across however many ranks run, each rank counting
// only its share. The trigger must open and close on the same days whatever
// the split, and every rank must agree.
TEST_CASE("Trigger - activation days do not depend on rank count") {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  WorldState world;
  const Disease disease =
      makeStageDisease(std::make_shared<ConstantCurve>(1.0));
  PolicyManager pm(world);
  pm.observedMetrics().declare("deaths", {"death", "world", "", ""});
  TemporalPolicy policy;
  policy.name = "close";
  policy.action.replacement_activity = "residence";
  PolicyTrigger trigger;
  trigger.metric = "deaths";
  trigger.cumulative = true;
  trigger.activate_at = 10;
  trigger.calm_days = 4;
  policy.trigger = trigger;
  pm.addTemporalPolicy(policy);
  pm.observedMetrics().resolve(world, disease);
  pm.requireTriggerMetrics();

  const Person nobody{};
  // Deaths per day: 2, 3, 0, 5 reach 10 on day 3; four calm days then reopen
  // it on day 7. Deaths are numbered from 0 each day and dealt to ranks by
  // number, so the ranks' local counts differ.
  const int deaths_per_day[] = {2, 3, 0, 5, 0, 0, 0, 0, 0, 0};
  int64_t local_total = 0;
  std::vector<int64_t> in_force_next_day;
  for (int day = 0; day < 10; ++day) {
    for (int id = 0; id < deaths_per_day[day]; ++id) {
      if (id % size == rank) {
        pm.observedMetrics().onDeath(nobody);
        ++local_total;
      }
    }
    pm.observedMetrics().reduce();
    pm.advanceTriggers(day);
    in_force_next_day.push_back(
        pm.getTemporalPolicies()[0].window.contains(day + 1.0) ? 1 : 0);
  }

  // Closed from the first slot of day 4 through the last of day 7.
  const std::vector<int64_t> expected = {0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
  CHECK(in_force_next_day == expected);

  // Every rank holds the identical answer, though its own count differs.
  expectSameOnAllRanks(in_force_next_day);

  int64_t min_local = 0;
  int64_t max_local = 0;
  MPI_Allreduce(&local_total, &min_local, 1, MPI_INT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(&local_total, &max_local, 1, MPI_INT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  if (size > 1) CHECK(min_local != max_local);
}

// A school closes at two infected people and reopens at none. Each rank counts
// only its own people and holds only its own venues, yet every rank must close
// and reopen the same venues on the same days, and answer the override question
// for a venue it does not hold the same way.
TEST_CASE("Trigger - venue closures do not depend on rank count") {
  using namespace venue_scope_fixture;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  WorldState world = makeWorld(rank, size);
  const Disease disease =
      makeStageDisease(std::make_shared<ConstantCurve>(1.0));
  PolicyManager pm(world);
  pm.observedMetrics().declare(
      "cases", {"infection", "parent_venue", "primary_activity", ""});
  TemporalPolicy policy;
  policy.name = "close";
  policy.action.replacement_activity = "residence";
  policy.action.override_activities = {"primary_activity"};
  PolicyTrigger trigger;
  trigger.metric = "cases";
  trigger.activate_at = 2;
  trigger.deactivate_at = 0;
  trigger.scope = MetricScope::ParentVenue;
  policy.trigger = trigger;
  pm.addTemporalPolicy(policy);
  pm.observedMetrics().resolve(world, disease);
  pm.requireTriggerMetrics();
  pm.observedMetrics().buildVenueParents(world);
  for (Person& p : world.people)
    policy_test::importApplicable(pm, p, PolicyMask{}, 1);
  pm.resolveAll(disease);

  // Who is infected on each day. The staff member registered to the school
  // itself and the person registered to two classrooms both count at school 3;
  // the pub worker counts at the pub, a unit of its own.
  const std::vector<std::vector<PersonId>> present = {
      {kPupilA, kPupilB},
      {kPupilA},
      {},
      {kTwoClasses, kPubWorker, kSchoolStaff},
      {kOtherSchool},
  };
  const std::vector<int64_t> school_closed_next_day = {1, 1, 0, 1, 0};

  std::vector<int64_t> school_closed;
  std::vector<int64_t> pub_closed;
  std::vector<int64_t> child_overridden;
  for (int day = 0; day < static_cast<int>(present.size()); ++day) {
    for (PersonId id : present[day])
      if (id % size == rank)
        pm.observedMetrics().onPresent(world.people[id], 1);
    pm.observedMetrics().reduce();
    pm.advanceTriggers(day);
    const PolicyTrigger& t = *pm.getTemporalPolicies()[0].trigger;
    school_closed.push_back(t.unitClosed(kSchool) ? 1 : 0);
    pub_closed.push_back(t.unitClosed(kPub) ? 1 : 0);
    // Class B of the school, asked about by every rank for a person from the
    // other school, whether or not this rank holds the venue.
    child_overridden.push_back(
        pm.getOverride(world.people[kOtherSchool], kPrimary, kClassB, 0,
                       SlotVenueType::fromVenue(kClassB), day + 1.0, 0)
                .has_value()
            ? 1
            : 0);
  }

  CHECK(school_closed == school_closed_next_day);
  CHECK(child_overridden == school_closed_next_day);
  CHECK(pub_closed == std::vector<int64_t>(present.size(), 0));

  for (const auto* mine : {&school_closed, &child_overridden})
    expectSameOnAllRanks(*mine);
}

// A voluntary response counted per classroom. Each rank counts only its share
// of the infected people and holds only its own venues, yet the level at every
// classroom must be the same double on every rank and match the one worked out
// by hand from the whole count, and every rank must reach the same verdict for
// every person, including a person registered to a venue it does not hold.
TEST_CASE("Response - levels and verdicts do not depend on rank count") {
  using namespace venue_scope_fixture;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  WorldState world = makeWorld(rank, size);
  const Disease disease =
      makeStageDisease(std::make_shared<ConstantCurve>(1.0));
  PolicyManager pm(world);
  pm.setBaseSeed(77);
  pm.observedMetrics().declare("cases",
                               {"infection", "venue", "primary_activity", ""});
  TemporalPolicy policy;
  policy.name = "keep_home";
  policy.action.replacement_activity = "residence";
  policy.action.override_activities = {"primary_activity"};
  PolicyTrigger trigger;
  trigger.scope = MetricScope::Venue;
  trigger.terms = {{"cases", 0.25}};
  policy.trigger = trigger;
  pm.addTemporalPolicy(policy);
  pm.observedMetrics().resolve(world, disease);
  pm.requireTriggerMetrics();
  pm.observedMetrics().buildVenueParents(world);
  for (Person& p : world.people)
    policy_test::importApplicable(pm, p, PolicyMask{}, 1);
  pm.resolveAll(disease);

  // Class A has two cases (pupil A and the person in both classes), class B
  // two (pupil B and the person in both), class C one. Each is counted by the
  // rank that owns the person, so the ranks' local counts differ.
  const std::vector<PersonId> present = {kPupilA, kPupilB, kTwoClasses,
                                         kOtherSchool};
  for (PersonId id : present)
    if (id % size == rank) pm.observedMetrics().onPresent(world.people[id], 1);
  pm.observedMetrics().reduce();
  pm.advanceTriggers(0);

  const PolicyTrigger& t = *pm.getTemporalPolicies()[0].trigger;
  CHECK(t.unit_levels.at(kClassA) == 0.0 + 0.25 * 2.0);
  CHECK(t.unit_levels.at(kClassB) == 0.0 + 0.25 * 2.0);
  CHECK(t.unit_levels.at(kClassC) == 0.0 + 0.25 * 1.0);
  CHECK(t.unit_levels.size() == 3);

  // Every person registered to a classroom, asked at that classroom, against
  // the draw worked out here from the specification.
  struct Seat {
    PersonId person;
    VenueId classroom;
    double level;
  };
  const std::vector<Seat> seats = {{kPupilA, kClassA, 0.5},
                                   {kPupilB, kClassB, 0.5},
                                   {kTwoClasses, kClassA, 0.5},
                                   {kTwoClasses, kClassB, 0.5},
                                   {kOtherSchool, kClassC, 0.25}};
  const uint64_t kDrawKey = policyDrawKey(PolicyKind::Temporal, "keep_home");
  std::vector<int64_t> verdicts;
  for (const Seat& seat : seats) {
    const bool avoids =
        pm.getOverride(world.people[seat.person], kPrimary, seat.classroom, 0,
                       SlotVenueType::fromVenue(seat.classroom), 1.0, 0)
            .has_value();
    SplitMix64 rng(mix_seed(77, seat.person, kDrawKey, 0xA701D4ULL));
    CHECK(avoids == (seat.level > uniform01(rng)));
    verdicts.push_back(avoids ? 1 : 0);
  }

  // Identical on every rank, bit for bit.
  expectSameOnAllRanks(verdicts);
  expectSameOnAllRanks(std::vector<double>{t.unit_levels.at(kClassA),
                                           t.unit_levels.at(kClassB),
                                           t.unit_levels.at(kClassC)});
}

// The startup check sums its counts across ranks and every rank names the same
// lowest failing person, though only one of them holds that person's record.
TEST_CASE("Response - the startup check agrees on every rank") {
  using namespace venue_scope_fixture;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Every rank holds only the people it owns, as in the real partition.
  WorldState world = makeWorld(rank, size);
  std::erase_if(world.people,
                [&](const Person& p) { return p.id % size != rank; });

  const Disease disease =
      makeStageDisease(std::make_shared<ConstantCurve>(1.0));
  PolicyManager pm(world);
  pm.observedMetrics().declare("deaths", {"death", "world", "", ""});
  TemporalPolicy policy;
  policy.name = "keep_home";
  policy.action.replacement_activity = "residence";
  policy.action.override_activities = {"primary_activity"};
  PolicyTrigger trigger;
  trigger.terms = {{"deaths", 0.1}};
  trigger.decision_unit.shared_activity = "residence";
  policy.trigger = trigger;
  pm.addTemporalPolicy(policy);
  pm.observedMetrics().resolve(world, disease);
  pm.requireTriggerMetrics();
  pm.resolveAll(disease);

  CHECK_NOTHROW(pm.checkDecisionUnits());

  // Person 3 and person 4 lose their residence; the owner of 3 knows it.
  for (const PersonId id : {kPubWorker, kSchoolStaff}) {
    for (Person& p : world.people)
      if (p.id == id)
        world.activity_meta[p.activity_meta_start].venue_count = 0;
  }
  std::string msg;
  try {
    pm.checkDecisionUnits();
  } catch (const std::exception& e) {
    msg = e.what();
  }
  CHECK(msg.find("person 3 has 0 venues for activity 'residence'") !=
        std::string::npos);
  CHECK(msg.find("2 of 6 people checked failed") != std::string::npos);
}

#else
TEST_CASE("MPI policy trigger requires MPI") {}
#endif
