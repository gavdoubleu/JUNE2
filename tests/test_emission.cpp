#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <memory>
#include <vector>

#include "core/types.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/emission/emission.h"

using namespace june;

// Symptom ids: 0 = healthy (infected, not infectious), 1 = exposed, 2 = mild.
static constexpr uint16_t kHealthy = 0;
static constexpr uint16_t kExposed = 1;
static constexpr uint16_t kMild = 2;

// Two direct modes emitting different constants while mild, then a fomite
// mode sub-binned at `fomite_sub_bin_hours` depositing on every symptom.
static Disease makeDisease(double fomite_sub_bin_hours) {
  TransmissionParams transmission;
  transmission.mode = InfectiousnessMode::STAGE_DRIVEN;

  TransmissionMode first_direct;
  first_direct.name = "first_direct";
  first_direct.symptom_curves = {
      nullptr, nullptr, std::make_shared<LinearRampCurve>(1.0, 3.0, 1.0)};
  transmission.modes.push_back(std::move(first_direct));

  TransmissionMode second_direct;
  second_direct.name = "second_direct";
  second_direct.symptom_curves = {nullptr, nullptr,
                                  std::make_shared<ConstantCurve>(0.3)};
  transmission.modes.push_back(std::move(second_direct));

  TransmissionMode fomite;
  fomite.name = "fomite";
  fomite.type = TransmissionModeType::Fomite;
  fomite.symptom_curves = {nullptr, nullptr, nullptr};
  FomiteConfig fomite_config;
  fomite_config.sub_bin_time = fomite_sub_bin_hours;
  fomite_config.infectiousness_curve = std::make_shared<ConstantCurve>(1.0);
  fomite_config.deposition_by_symptom = {
      std::make_shared<ConstantCurve>(0.7),
      std::make_shared<ConstantCurve>(2.0),
      std::make_shared<LinearRampCurve>(1.0, 3.0, 1.0)};
  fomite.config = std::move(fomite_config);
  transmission.modes.push_back(std::move(fomite));

  std::vector<SymptomTag> symptom_tags = {
      {"healthy", -1, kHealthy}, {"exposed", 0, kExposed}, {"mild", 1, kMild}};
  TrajectoryDefinition trajectory;
  trajectory.selection_key = "general";
  trajectory.severity = 1.0;
  trajectory.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  return Disease("TestEmission", symptom_tags, {}, {trajectory}, {},
                 transmission);
}

// A Person infected at `infection_time` following exactly `transitions`
// (time, symptom id), with a cold symptom cache.
static Person makeInfectedPerson(
    const Disease& disease, double infection_time,
    std::vector<std::pair<double, uint16_t>> transitions) {
  InfectionTrajectory trajectory;
  trajectory.infection_time = infection_time;
  trajectory.transitions = std::move(transitions);
  Person person{};
  person.infection = Infection::fromCheckpoint(
      &disease, infection_time, trajectory, 1.0, 1.0, 1.0, 0.0,
      /*last_checked_time=*/-1.0, kHealthy, infection_time);
  return person;
}

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

TEST_CASE("Infected, not yet infectious: deposits only") {
  Disease disease = makeDisease(2.0);
  EmissionCalculator calculator(disease, 6.0);
  // Healthy (incubating) until 11, so not infectious at slot start 10.
  Person person =
      makeInfectedPerson(disease, 9.0, {{9.0, kHealthy}, {11.0, kMild}});
  REQUIRE_FALSE(person.infection->isInfectious(10.0));

  Emission emission;
  calculator.emit(person, 10.0, emission);

  std::vector<double> expected_deposits;
  calculator.fomiteSchedule().integrateDeposits(person.infection.get(), 10.0,
                                                expected_deposits);
  CHECK(emission.infectiousness_by_mode.empty());
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
  REQUIRE(person.infection->isInfectious(10.0));

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
