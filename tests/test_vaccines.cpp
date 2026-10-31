#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "core/types.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/vaccine.h"
#include "transmission_fixtures.h"

using namespace june;

TEST_CASE("Vaccine - Efficacy Interpolation") {
  Dose dose;
  dose.number = 0;
  dose.day_administered = 10.0;
  dose.days_to_effective = 14.0;
  dose.days_to_waning = 30.0;
  dose.days_to_finished = 60.0;
  dose.waning_factor = 0.5;
  dose.infection_efficacy["TestDisease"] = {{0, 100, 0.8}};

  SUBCASE("Before administration") {
    CHECK(std::abs(dose.getEfficacy(5.0, "TestDisease", 25.0) - 0.0) < 1e-6);
  }

  SUBCASE("Ramp up - halfway to effective") {
    CHECK(std::abs(dose.getEfficacy(17.0, "TestDisease", 25.0) - 0.4) < 1e-6);
  }

  SUBCASE("Peak effective") {
    CHECK(std::abs(dose.getEfficacy(30.0, "TestDisease", 25.0) - 0.8) < 1e-6);
  }

  SUBCASE("Mid-waning") {
    // Progress = (45 - 30) / (60 - 30) = 0.5
    // Efficacy = 0.8 - 0.5 * (0.8 - 0.4) = 0.6
    CHECK(std::abs(dose.getEfficacy(55.0, "TestDisease", 25.0) - 0.6) < 1e-6);
  }

  SUBCASE("Post-waning baseline") {
    CHECK(std::abs(dose.getEfficacy(80.0, "TestDisease", 25.0) - 0.4) < 1e-6);
  }
}

namespace {

// A dose given at the start of the run. Peak infection efficacy is 0.8,
// held until day 100, then waning to half the peak by day 200.
Dose makeOnsetDose(double onset, double effective, double administered = 0.0) {
  Dose dose;
  dose.number = 0;
  dose.day_administered = administered;
  dose.days_to_onset = onset;
  dose.days_to_effective = effective;
  dose.days_to_waning = 100.0;
  dose.days_to_finished = 200.0;
  dose.waning_factor = 0.5;
  dose.infection_efficacy["TestDisease"] = {{0, 100, 0.8}};
  dose.symptom_efficacy["TestDisease"] = {{0, 100, 0.6}};
  return dose;
}

}  // namespace

TEST_CASE("Vaccine - Onset Delay") {
  SUBCASE("Onset 0 is the plain ramp from the day of the dose") {
    const Dose dose = makeOnsetDose(0.0, 14.0);
    for (double day : {0.0, 3.5, 7.0, 13.999}) {
      // The ramp as it was before onset existed, to the last bit
      CHECK(dose.getEfficacy(day, "TestDisease", 25.0) == (day / 14.0) * 0.8);
    }
    CHECK(dose.getEfficacy(14.0, "TestDisease", 25.0) == 0.8);
  }

  SUBCASE("Onset equal to days_to_effective is a step") {
    const Dose dose = makeOnsetDose(14.0, 14.0);
    CHECK(dose.getEfficacy(0.0, "TestDisease", 25.0) == 0.0);
    CHECK(dose.getEfficacy(7.0, "TestDisease", 25.0) == 0.0);
    CHECK(dose.getEfficacy(13.999, "TestDisease", 25.0) == 0.0);
    CHECK(dose.getEfficacy(14.0, "TestDisease", 25.0) == 0.8);
    // The peak is held right up to the start of waning
    CHECK(dose.getEfficacy(50.0, "TestDisease", 25.0) == 0.8);
    CHECK(dose.getEfficacy(99.999, "TestDisease", 25.0) == 0.8);
  }

  SUBCASE("A delay and then a rise") {
    const Dose dose = makeOnsetDose(7.0, 14.0);
    CHECK(dose.getEfficacy(0.0, "TestDisease", 25.0) == 0.0);
    CHECK(dose.getEfficacy(3.0, "TestDisease", 25.0) == 0.0);
    CHECK(dose.getEfficacy(7.0, "TestDisease", 25.0) == 0.0);
    CHECK(dose.getEfficacy(10.5, "TestDisease", 25.0) ==
          doctest::Approx(0.4).epsilon(1e-12));
    CHECK(dose.getEfficacy(14.0, "TestDisease", 25.0) == 0.8);
  }

  SUBCASE("Waning after the peak does not depend on onset") {
    const Dose step = makeOnsetDose(14.0, 14.0);
    const Dose ramp = makeOnsetDose(0.0, 14.0);
    const Dose delayed = makeOnsetDose(7.0, 14.0);
    for (double day : {14.0, 50.0, 100.0, 125.0, 150.0, 199.0, 200.0, 300.0}) {
      CHECK(step.getEfficacy(day, "TestDisease", 25.0) ==
            ramp.getEfficacy(day, "TestDisease", 25.0));
      CHECK(delayed.getEfficacy(day, "TestDisease", 25.0) ==
            ramp.getEfficacy(day, "TestDisease", 25.0));
    }
    // Halfway through waning: 0.8 - 0.5 * (0.8 - 0.4)
    CHECK(step.getEfficacy(150.0, "TestDisease", 25.0) ==
          doctest::Approx(0.6).epsilon(1e-12));
    CHECK(step.getEfficacy(300.0, "TestDisease", 25.0) ==
          doctest::Approx(0.4).epsilon(1e-12));
  }

  SUBCASE("Onset counts from the day of the dose") {
    const Dose dose = makeOnsetDose(14.0, 14.0, /*administered=*/10.0);
    CHECK(dose.getEfficacy(9.0, "TestDisease", 25.0) == 0.0);
    CHECK(dose.getEfficacy(23.999, "TestDisease", 25.0) == 0.0);
    CHECK(dose.getEfficacy(24.0, "TestDisease", 25.0) == 0.8);
  }

  SUBCASE("Symptom efficacy follows the same onset") {
    const Dose dose = makeOnsetDose(14.0, 14.0);
    CHECK(dose.getEfficacy(13.0, "TestDisease", 25.0, true) == 0.0);
    CHECK(dose.getEfficacy(14.0, "TestDisease", 25.0, true) == 0.6);
  }
}

TEST_CASE("Vaccine - Onset Delay Through Susceptibility") {
  Person person;
  person.age = 30;
  person.sex = Sex::MALE;
  person.immunity.natural_level = 0.0;
  person.vaccine_trajectory = std::make_unique<VaccineTrajectory>();
  person.vaccine_trajectory->addDose(makeOnsetDose(14.0, 14.0));

  // Fully susceptible the day before onset, protected from the day of onset
  CHECK(person.getSusceptibility(13.0, "TestDisease") == 1.0);
  CHECK(person.getSusceptibility(14.0, "TestDisease") ==
        doctest::Approx(1.0 - 0.8));
}

TEST_CASE("Vaccine - Hybrid Immunity Aggregation") {
  Person person;
  person.age = 30;
  person.sex = Sex::MALE;

  SUBCASE("No immunity - full susceptibility") {
    CHECK(std::abs(person.getSusceptibility(0.0, "TestDisease") - 1.0) < 1e-6);
  }

  SUBCASE("Natural immunity only - 50% protection") {
    person.immunity.natural_level = 0.5;
    person.immunity.natural_acquisition_time = 0.0;
    person.immunity.natural_waning_rate = 0.0;
    CHECK(std::abs(person.getSusceptibility(10.0, "TestDisease") - 0.5) < 1e-6);
  }

  SUBCASE("Vaccine immunity only - 80% protection") {
    person.immunity.natural_level = 0.0;
    person.vaccine_trajectory = std::make_unique<VaccineTrajectory>();

    Dose dose;
    dose.day_administered = 0.0;
    dose.days_to_effective = 0.0;
    dose.days_to_waning = 100.0;
    dose.infection_efficacy["TestDisease"] = {{0, 100, 0.8}};
    person.vaccine_trajectory->addDose(dose);

    CHECK(std::abs(person.getSusceptibility(10.0, "TestDisease") - 0.2) < 1e-6);
  }

  SUBCASE("Hybrid immunity - multiplicative") {
    person.immunity.natural_level = 0.5;
    person.immunity.natural_acquisition_time = 0.0;
    person.immunity.natural_waning_rate = 0.0;
    person.vaccine_trajectory = std::make_unique<VaccineTrajectory>();

    Dose dose;
    dose.day_administered = 0.0;
    dose.days_to_effective = 0.0;
    dose.days_to_waning = 100.0;
    dose.infection_efficacy["TestDisease"] = {{0, 100, 0.8}};
    person.vaccine_trajectory->addDose(dose);

    // natural_suscep = 0.5, vaccine_suscep = 0.2 → total = 0.1
    CHECK(std::abs(person.getSusceptibility(10.0, "TestDisease") - 0.1) < 1e-6);
  }
}

// =============================================================================
// Symptom efficacy at infection time
// =============================================================================

namespace {

// Symptom ids of the three trajectories, which each start in a symptom of
// their own so the one an infection took can be read off its first stage.
constexpr uint16_t kModerate = 1;
constexpr uint16_t kLight = 2;
constexpr uint16_t kSevere = 3;

constexpr uint64_t kInfections = 20000;

// Three trajectories that differ only in severity and in the symptom they
// start in. The least severe is listed in the middle, so moving probability to
// it cannot work by it being first or last.
//   moderate: severity 1, probability 0.3
//   light:    severity 0, probability 0.6
//   severe:   severity 2, probability 0.1
std::unique_ptr<Disease> makeShiftDisease() {
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(1.0);
  trans.symptom_id_curves = {nullptr, curve, curve, curve};

  const auto trajectory = [](const char* symptom, double severity,
                             double probability) {
    TrajectoryDefinition td;
    td.description = symptom;
    td.severity = severity;
    td.probability = probability;
    td.stages.push_back({symptom, {"constant", {{"value", 5.0}}}});
    td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});
    return td;
  };

  const std::vector<SymptomTag> tags = {
      {.name = "healthy", .value = -1, .id = 0},
      {.name = "moderate", .value = 2, .id = kModerate},
      {.name = "light", .value = 1, .id = kLight},
      {.name = "severe", .value = 3, .id = kSevere}};
  DiseaseStageSettings stage_settings;
  stage_settings.recovered_stages = {"healthy"};
  return std::make_unique<Disease>(
      "TestFlu", tags, stage_settings,
      std::vector<TrajectoryDefinition>{trajectory("moderate", 1.0, 0.3),
                                        trajectory("light", 0.0, 0.6),
                                        trajectory("severe", 2.0, 0.1)},
      OutcomeRates{}, trans);
}

Person makePerson(float age) {
  Person person;
  person.age = age;
  return person;
}

// A dose given at day 0 that protects against nothing yet
Dose makeBareDose(double onset, double effective) {
  Dose dose;
  dose.number = 0;
  dose.day_administered = 0.0;
  dose.days_to_onset = onset;
  dose.days_to_effective = effective;
  dose.days_to_waning = 1000.0;
  dose.days_to_finished = 2000.0;
  return dose;
}

Dose makeSymptomDose(const std::string& disease,
                     const std::vector<AgeEfficacy>& bands, double onset,
                     double effective) {
  Dose dose = makeBareDose(onset, effective);
  dose.symptom_efficacy[disease] = bands;
  return dose;
}

void vaccinate(Person& person, const Dose& dose) {
  person.vaccine_trajectory = std::make_unique<VaccineTrajectory>();
  person.vaccine_trajectory->addDose(dose);
}

// The trajectory each of kInfections infections took, one seed per infection.
// The same seed draws the same number whoever is infected, so two people can
// be compared seed by seed.
std::vector<uint16_t> outcomes(const Disease& disease, const Person& person,
                               double infection_time) {
  std::vector<uint16_t> taken;
  taken.reserve(kInfections);
  for (uint64_t seed = 1; seed <= kInfections; ++seed) {
    Infection infection(&disease, infection_time, &person, seed,
                        kNoTransmissionContext, nullptr, "office", 0);
    taken.push_back(
        infection.getTrajectory().getCurrentSymptomId(infection_time));
  }
  return taken;
}

double share(const std::vector<uint16_t>& taken, uint16_t symptom) {
  double count = 0.0;
  for (uint16_t t : taken) count += t == symptom;
  return count / static_cast<double>(taken.size());
}

// A share over 20000 draws has a standard deviation of at most 0.0035, so 0.015
// is more than four of them. The seeds are fixed, so a run never differs from
// the last.
constexpr double kShareTolerance = 0.015;

void checkShare(const std::vector<uint16_t>& taken, uint16_t symptom,
                double expected) {
  CAPTURE(symptom);
  CAPTURE(expected);
  CHECK(std::abs(share(taken, symptom) - expected) < kShareTolerance);
}

// Every trajectory above the least severe keeps (1 - ve) of its probability,
// the least severe takes what they lost, and nobody falls outside the three.
void checkShifted(const std::vector<uint16_t>& taken, double ve) {
  checkShare(taken, kModerate, 0.3 * (1.0 - ve));
  checkShare(taken, kSevere, 0.1 * (1.0 - ve));
  checkShare(taken, kLight, 0.6 + 0.4 * ve);
  CHECK(share(taken, kModerate) + share(taken, kLight) +
            share(taken, kSevere) ==
        doctest::Approx(1.0));
}

}  // namespace

TEST_CASE("Vaccine - Symptom Efficacy Leaves Rates Alone Without Protection") {
  const auto disease = makeShiftDisease();
  const Person unvaccinated = makePerson(30.0f);
  const std::vector<uint16_t> baseline = outcomes(*disease, unvaccinated, 0.0);

  SUBCASE("The unvaccinated take each trajectory at its own probability") {
    checkShare(baseline, kModerate, 0.3);
    checkShare(baseline, kLight, 0.6);
    checkShare(baseline, kSevere, 0.1);
  }

  SUBCASE("A dose with symptom efficacy 0 changes nothing") {
    Person person = makePerson(30.0f);
    vaccinate(person, makeSymptomDose("TestFlu", {{0, 100, 0.0}}, 0.0, 0.0));
    CHECK(outcomes(*disease, person, 0.0) == baseline);
  }

  SUBCASE("Infection efficacy alone does not move the outcome") {
    Person person = makePerson(30.0f);
    Dose dose = makeBareDose(0.0, 0.0);
    dose.infection_efficacy["TestFlu"] = {{0, 100, 1.0}};
    vaccinate(person, dose);
    CHECK(outcomes(*disease, person, 0.0) == baseline);
  }
}

TEST_CASE("Vaccine - Symptom Efficacy Moves Probability To The Least Severe") {
  const auto disease = makeShiftDisease();

  for (double ve : {0.25, 0.5, 0.9}) {
    CAPTURE(ve);
    Person person = makePerson(30.0f);
    vaccinate(person, makeSymptomDose("TestFlu", {{0, 100, ve}}, 0.0, 0.0));
    checkShifted(outcomes(*disease, person, 0.0), ve);
  }

  SUBCASE("Full efficacy leaves only the least severe trajectory") {
    Person person = makePerson(30.0f);
    vaccinate(person, makeSymptomDose("TestFlu", {{0, 100, 1.0}}, 0.0, 0.0));
    CHECK(share(outcomes(*disease, person, 0.0), kLight) == 1.0);
  }
}

TEST_CASE("Vaccine - Symptom Efficacy Is Read At Infection Time") {
  const auto disease = makeShiftDisease();
  const Person unvaccinated = makePerson(30.0f);

  SUBCASE("Nothing before onset, the peak from onset") {
    Person person = makePerson(30.0f);
    vaccinate(person, makeSymptomDose("TestFlu", {{0, 100, 1.0}}, 14.0, 14.0));
    CHECK(outcomes(*disease, person, 13.0) ==
          outcomes(*disease, unvaccinated, 13.0));
    CHECK(share(outcomes(*disease, person, 14.0), kLight) == 1.0);
    CHECK(share(outcomes(*disease, person, 100.0), kLight) == 1.0);
  }

  SUBCASE("Halfway up the ramp it is half the peak") {
    Person person = makePerson(30.0f);
    vaccinate(person, makeSymptomDose("TestFlu", {{0, 100, 0.8}}, 0.0, 20.0));
    checkShifted(outcomes(*disease, person, 10.0), 0.4);
  }
}

TEST_CASE("Vaccine - Symptom Efficacy Matches Disease And Age Band") {
  const auto disease = makeShiftDisease();
  const Person unvaccinated = makePerson(30.0f);
  const std::vector<uint16_t> baseline = outcomes(*disease, unvaccinated, 0.0);

  SUBCASE("Efficacy against another disease does not apply") {
    Person person = makePerson(30.0f);
    vaccinate(person, makeSymptomDose("OtherFlu", {{0, 100, 1.0}}, 0.0, 0.0));
    CHECK(outcomes(*disease, person, 0.0) == baseline);
  }

  const std::vector<AgeEfficacy> bands = {{0, 17, 1.0}, {18, 64, 0.5}};

  SUBCASE("The band of the person's age sets the efficacy") {
    Person child = makePerson(10.0f);
    vaccinate(child, makeSymptomDose("TestFlu", bands, 0.0, 0.0));
    CHECK(share(outcomes(*disease, child, 0.0), kLight) == 1.0);

    Person adult = makePerson(40.0f);
    vaccinate(adult, makeSymptomDose("TestFlu", bands, 0.0, 0.0));
    checkShifted(outcomes(*disease, adult, 0.0), 0.5);
  }

  SUBCASE("A person outside every band is unchanged") {
    Person elder = makePerson(80.0f);
    vaccinate(elder, makeSymptomDose("TestFlu", bands, 0.0, 0.0));
    CHECK(outcomes(*disease, elder, 0.0) == baseline);
  }
}
