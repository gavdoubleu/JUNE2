#include <memory>
#include <vector>

#include "core/types.h"
#include "doctest.h"
#include "emission_fixtures.h"
#include "epidemiology/disease.h"
#include "epidemiology/emission/emission.h"

using namespace june;
using namespace emission_fixtures;

TEST_CASE("Uninfected Person emits nothing") {
  Disease disease = makeDisease(2.0);
  EmissionCalculator calculator(disease, 6.0);
  Person person{};

  Emission emission;
  emission.infectiousness_by_mode = {1.0};  // stale contents must not survive
  emission.fomite_deposits = {1.0};
  calculator.emit(person, 10.0, emission);

  CHECK(emission.infectiousness_by_mode.empty());
  CHECK(emission.fomite_deposits.empty());
}

TEST_CASE("Incubating: zero per mode, deposits as usual") {
  Disease disease = makeDisease(2.0);
  EmissionCalculator calculator(disease, 6.0);
  // Healthy (incubating) for the whole slot 10 to 10.25; mild from 11.
  Person person =
      makeInfectedPerson(disease, 9.0, {{9.0, kHealthy}, {11.0, kMild}});

  Emission emission;
  calculator.emit(person, 10.0, emission);

  std::vector<double> expected_deposits;
  calculator.fomiteSchedule().integrateDeposits(person.infection.get(), 10.0,
                                                expected_deposits);
  CHECK(emission.infectiousness_by_mode == std::vector<double>{0.0, 0.0, 0.0});
  REQUIRE(emission.fomite_deposits.size() == 3);
  CHECK(emission.fomite_deposits == expected_deposits);
  CHECK(emission.fomite_deposits[0] > 0.0);
}

TEST_CASE("Infectious: per-mode integrals and deposits") {
  Disease disease = makeDisease(2.0);
  EmissionCalculator calculator(disease, 6.0);
  // Mild from 9.5, so infectious at slot start 10; the ramp makes values
  // depend on the exact time in stage.
  Person person =
      makeInfectedPerson(disease, 9.0, {{9.0, kHealthy}, {9.5, kMild}});

  Emission emission;
  calculator.emit(person, 10.0, emission);

  const double slot_end = 10.0 + 6.0 / 24.0;
  REQUIRE(emission.infectiousness_by_mode.size() == 3);
  for (int mode = 0; mode < 3; ++mode) {
    CHECK(emission.infectiousness_by_mode[mode] ==
          person.infection->getIntegratedInfectiousness(mode, 10.0, slot_end));
  }
  CHECK(emission.infectiousness_by_mode[0] > 0.0);
  CHECK(emission.infectiousness_by_mode[1] > 0.0);

  std::vector<double> expected_deposits;
  calculator.fomiteSchedule().integrateDeposits(person.infection.get(), 10.0,
                                                expected_deposits);
  CHECK(emission.fomite_deposits == expected_deposits);
}

TEST_CASE("Onset mid-slot: the infectious part of the slot counts") {
  Disease disease = makeDisease(2.0);
  EmissionCalculator calculator(disease, 6.0);
  // Healthy at slot start 10, mild from 10.1, before the slot ends at 10.25.
  Person person =
      makeInfectedPerson(disease, 9.0, {{9.0, kHealthy}, {10.1, kMild}});

  Emission emission;
  calculator.emit(person, 10.0, emission);

  const double slot_end = 10.0 + 6.0 / 24.0;
  REQUIRE(emission.infectiousness_by_mode.size() == 3);
  for (int mode = 0; mode < 3; ++mode) {
    CHECK(emission.infectiousness_by_mode[mode] ==
          person.infection->getIntegratedInfectiousness(mode, 10.0, slot_end));
  }
  // Constant 0.3 while mild, over the 0.15 days of mild in the slot.
  CHECK(emission.infectiousness_by_mode[1] ==
        doctest::Approx(24.0 * 0.3 * 0.15));
}

TEST_CASE("Trajectory-Driven: recovered or dead at slot start emits zero") {
  // The gamma profile ignores stage, so without the t0 check a Person who is
  // recovered or dead at slot start would keep emitting.
  constexpr uint16_t kRecovered = 0;
  constexpr uint16_t kMildStage = 1;
  constexpr uint16_t kDead = 2;
  TransmissionParams transmission;
  transmission.mode = InfectiousnessMode::TRAJECTORY_DRIVEN;
  transmission.type = "gamma";
  DiseaseStageSettings stage_settings;
  stage_settings.recovered_stages = {"recovered"};
  stage_settings.fatality_stages = {"dead"};
  TrajectoryDefinition trajectory_definition;
  trajectory_definition.selection_key = "general";
  trajectory_definition.severity = 1.0;
  trajectory_definition.stages.push_back(
      {"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("TrajectoryEmission",
                  {{"recovered", -2, kRecovered},
                   {"mild", 1, kMildStage},
                   {"dead", 2, kDead}},
                  stage_settings, {trajectory_definition}, {}, transmission);
  EmissionCalculator calculator(disease, 6.0);

  auto makePerson = [&](std::vector<std::pair<double, uint16_t>> transitions) {
    InfectionTrajectory trajectory;
    trajectory.infection_time = 9.0;
    trajectory.transitions = std::move(transitions);
    Person person{};
    person.infection = Infection::fromCheckpoint(
        &disease, 9.0, trajectory, /*max_infectiousness=*/1.0,
        /*transmission_shape=*/2.0, /*transmission_rate=*/1.0,
        /*transmission_shift=*/0.0, /*last_checked_time=*/-1.0, kMildStage,
        9.0);
    return person;
  };
  auto emitAtTen = [&](const Person& person) {
    Emission emission;
    calculator.emit(person, 10.0, emission);
    REQUIRE(emission.infectiousness_by_mode.size() == 1);
    return emission.infectiousness_by_mode[0];
  };

  CHECK(emitAtTen(makePerson({{9.0, kMildStage}, {9.9, kRecovered}})) == 0.0);
  CHECK(emitAtTen(makePerson({{9.0, kMildStage}, {10.0, kDead}})) == 0.0);
  // Recovering after slot start still counts the whole slot.
  CHECK(emitAtTen(makePerson({{9.0, kMildStage}, {10.1, kRecovered}})) > 0.0);
}

TEST_CASE("Slot length fixes both the sub-bin schedule and the integrals") {
  Disease disease = makeDisease(2.0);
  EmissionCalculator six_hour(disease, 6.0);
  EmissionCalculator four_hour(disease, 4.0);
  Person person =
      makeInfectedPerson(disease, 9.0, {{9.0, kHealthy}, {9.5, kMild}});

  Emission six_hour_emission;
  Emission four_hour_emission;
  six_hour.emit(person, 10.0, six_hour_emission);
  four_hour.emit(person, 10.0, four_hour_emission);

  CHECK(six_hour.fomiteSchedule().subBinsPerMode() == std::vector<int>{3});
  CHECK(four_hour.fomiteSchedule().subBinsPerMode() == std::vector<int>{2});
  CHECK(six_hour_emission.fomite_deposits.size() == 3);
  CHECK(four_hour_emission.fomite_deposits.size() == 2);
  CHECK(four_hour_emission.infectiousness_by_mode[1] ==
        person.infection->getIntegratedInfectiousness(1, 10.0,
                                                      10.0 + 4.0 / 24.0));
  CHECK(four_hour_emission.infectiousness_by_mode[1] <
        six_hour_emission.infectiousness_by_mode[1]);
}

TEST_CASE("Symptom id before the first transition is the first stage") {
  Disease disease = makeDisease(2.0);
  // Infected at 8, exposed from 9, mild from 11.
  Person person =
      makeInfectedPerson(disease, 8.0, {{9.0, kExposed}, {11.0, kMild}});
  const InfectionTrajectory& trajectory = person.infection->getTrajectory();

  // Symptom 0 means recovered, so it is never the pre-start answer.
  CHECK(trajectory.getCurrentSymptomId(8.0) == kExposed);
  CHECK(trajectory.getCurrentSymptomId(8.99) == kExposed);
  CHECK(trajectory.getCurrentSymptomId(9.0) == kExposed);  // at a transition
  CHECK(trajectory.getCurrentSymptomId(10.0) == kExposed);
  CHECK(trajectory.getCurrentSymptomId(11.0) == kMild);
  CHECK(trajectory.getCurrentSymptomId(50.0) == kMild);  // after the last
}
