#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
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

TEST_CASE("symptomIdAt walks transitions forward from symptom 0") {
  Disease disease = makeDisease(2.0);
  // Infected at 8, exposed from 9, mild from 11.
  Person person =
      makeInfectedPerson(disease, 8.0, {{9.0, kExposed}, {11.0, kMild}});
  const Infection& infection = *person.infection;

  CHECK(infection.symptomIdAt(8.0) == 0);  // before the first transition
  CHECK(infection.symptomIdAt(8.99) == 0);
  CHECK(infection.symptomIdAt(9.0) == kExposed);  // at a transition
  CHECK(infection.symptomIdAt(10.0) == kExposed);
  CHECK(infection.symptomIdAt(11.0) == kMild);
  CHECK(infection.symptomIdAt(50.0) == kMild);  // after the last
}

TEST_CASE("symptomIdAt leaves the checkpointed symptom cache alone") {
  Disease disease = makeDisease(2.0);
  Person person =
      makeInfectedPerson(disease, 8.0, {{9.0, kExposed}, {11.0, kMild}});
  const Infection& infection = *person.infection;

  infection.symptomIdAt(10.0);

  CHECK(infection.ckptLastCheckedTime() == -1.0);
  CHECK(infection.ckptCachedSymptomId() == kHealthy);
  CHECK(infection.ckptCachedSymptomStartTime() == 8.0);
}
