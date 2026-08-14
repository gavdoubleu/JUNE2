#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <memory>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/epidemiology.h"
#include "epidemiology/interaction_manager.h"
#include "epidemiology/vaccine.h"
#include "test_utils.h"

using namespace june;

// =============================================================================
// Helpers — build minimal but valid simulation objects
// =============================================================================

struct TestFixture {
  WorldState world;
  TransmissionParams trans;
  std::vector<SymptomTag> symptom_tags;
  std::vector<TrajectoryDefinition> trajectories;
  ContactMatrixConfig cm;
  SimulationConfig sim_cfg;

  TestFixture() {
    // --- World: 1 venue "office", 2 people ---
    Venue venue;
    venue.id = 0;
    venue.type_id = 0;
    world.venue_type_names = {"office"};
    world.encounter_type_names = {"romantic_encounter", "friend_hangout"};
    world.venues.push_back(venue);

    world.people.emplace_back();
    world.people.emplace_back();
    Person& infector = world.people[0];
    infector.id = 0;
    infector.age = 30;
    infector.sex = Sex::MALE;
    infector.geo_unit_id = -1;

    Person& target = world.people[1];
    target.id = 1;
    target.age = 25;
    target.sex = Sex::FEMALE;
    target.geo_unit_id = -1;

    world.buildIndices();

    // --- Disease: stage-driven, "mild" has constant infectiousness = 1.0 ---
    trans.mode = InfectiousnessMode::STAGE_DRIVEN;
    trans.natural_immunity.level = 0.95;
    trans.natural_immunity.waning_rate = 0.001;

    auto constant_curve = std::make_shared<ConstantCurve>(1.0);
    trans.stage_curves["mild"] = constant_curve;
    trans.symptom_id_curves = {nullptr, constant_curve};

    symptom_tags = {{"healthy", -1, 0}, {"mild", 1, 1}};

    TrajectoryDefinition td;
    td.selection_key = "general";
    td.severity = 1.0;
    td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
    trajectories.push_back(td);

    // Contact matrix: tuned so omega = contacts * delta / bin_size gives
    // distinguishable infection rates for vaccinated vs unvaccinated.
    // With contacts=0.15, delta=8h, I=1.0, bin_size=1: omega ≈ 1.2
    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{0.15}};
    cm.default_matrix = default_contact_matrix;
    cm.default_characteristic_time = 1.0;
  }

  Disease makeDisease() {
    return Disease("TestFlu", symptom_tags, {}, trajectories, {}, trans);
  }

  // Create locations: both people at venue 0
  std::vector<PersonLocation> makeColocatedLocations() {
    std::vector<PersonLocation> locs;
    PersonLocation loc0;
    loc0.person_id = 0;
    loc0.venue_id = 0;
    loc0.subset_index = -1;
    loc0.activity_index = -1;
    loc0.encounter_type_id = 255;
    loc0.person_array_index = 0;
    locs.push_back(loc0);

    PersonLocation loc1;
    loc1.person_id = 1;
    loc1.venue_id = 0;
    loc1.subset_index = -1;
    loc1.activity_index = -1;
    loc1.encounter_type_id = 255;
    loc1.person_array_index = 1;
    locs.push_back(loc1);
    return locs;
  }
};

// Run N transmission trials and return the fraction of times the target got
// infected.
double measureInfectionRate(TestFixture& fix, Disease& disease, int N = 500) {
  int infections = 0;
  auto locs = fix.makeColocatedLocations();

  for (int trial = 0; trial < N; ++trial) {
    // Fresh target each trial
    fix.world.people[1].infection.reset();

    // Infector stays infected with "mild" symptom
    if (!fix.world.people[0].infection) {
      fix.world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &fix.world.people[0], trial, nullptr, "office", 0);
    }

    ParallelConfig parallel_config;
    fix.cm.allow_default_matrix = true;
    finalizeContactMatrices(fix.cm, fix.world);
    InteractionManager im(fix.world, fix.cm, fix.sim_cfg, parallel_config,
                          &disease, nullptr);
    im.processTransmissions(locs, 5.0, 8.0, nullptr);

    if (fix.world.people[1].infection) {
      infections++;
    }
  }
  return static_cast<double>(infections) / N;
}

// =============================================================================
// Integration Test 1: Vaccine → Susceptibility → Transmission
// =============================================================================

TEST_CASE("Integration: Vaccine reduces transmission probability") {
  TestFixture fix;

  // Use moderate contacts so unvaccinated rate is high but not 100%
  // omega = 0.15 * 8 / 1 = 1.2 → prob ≈ 0.70
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{0.15}};
  fix.cm.default_matrix = default_contact_matrix;

  Disease disease = fix.makeDisease();

  SUBCASE("Unvaccinated person has high infection probability") {
    // Person 1 is unvaccinated — susceptibility should be 1.0
    CHECK(fix.world.people[1].getSusceptibility(5.0, "TestFlu") ==
          doctest::Approx(1.0));

    double rate = measureInfectionRate(fix, disease);
    // With high beta * contacts * delta, unvaccinated should get infected most
    // of the time
    CHECK(rate > 0.5);
    MESSAGE("Unvaccinated infection rate: ", rate);
  }

  SUBCASE("Vaccinated person has reduced infection probability") {
    Person& target = fix.world.people[1];

    // Give target a vaccine with 80% infection efficacy for disease "TestFlu"
    target.vaccine_trajectory = std::make_unique<VaccineTrajectory>();
    target.vaccine_trajectory->vaccine_name = "TestVax";

    Dose dose;
    dose.number = 0;
    dose.day_administered = 0.0;
    dose.days_to_effective = 0.0;  // Immediate efficacy
    dose.days_to_waning = 365.0;   // No waning during test
    dose.days_to_finished = 730.0;
    dose.waning_factor = 0.5;

    // 80% infection efficacy for "TestFlu" across all ages
    AgeEfficacy all_ages;
    all_ages.min_age = 0;
    all_ages.max_age = 120;
    all_ages.efficacy = 0.8;
    dose.infection_efficacy["TestFlu"] = {all_ages};

    target.vaccine_trajectory->addDose(dose);

    // Verify susceptibility is reduced: 1.0 * (1.0 - 0.8) = 0.2
    double susc = target.getSusceptibility(5.0, "TestFlu");
    CHECK(susc == doctest::Approx(0.2).epsilon(0.01));
    MESSAGE("Vaccinated susceptibility: ", susc);

    double rate = measureInfectionRate(fix, disease);
    // Vaccinated rate should be substantially lower
    CHECK(rate < 0.5);
    MESSAGE("Vaccinated infection rate: ", rate);
  }

  SUBCASE("Vaccination reduces infection rate compared to baseline") {
    Disease disease2 = fix.makeDisease();

    // Measure unvaccinated baseline
    double unvax_rate = measureInfectionRate(fix, disease2, 1000);

    // Now vaccinate and measure again
    Person& target = fix.world.people[1];
    target.vaccine_trajectory = std::make_unique<VaccineTrajectory>();
    target.vaccine_trajectory->vaccine_name = "TestVax";

    Dose dose;
    dose.number = 0;
    dose.day_administered = 0.0;
    dose.days_to_effective = 0.0;
    dose.days_to_waning = 365.0;
    dose.days_to_finished = 730.0;
    dose.waning_factor = 0.5;
    AgeEfficacy all_ages{0, 120, 0.8};
    dose.infection_efficacy["TestFlu"] = {all_ages};
    target.vaccine_trajectory->addDose(dose);

    double vax_rate = measureInfectionRate(fix, disease2, 1000);

    MESSAGE("Unvaccinated rate: ", unvax_rate, ", Vaccinated rate: ", vax_rate);

    // The vaccinated rate should be clearly lower (at least 30% relative
    // reduction) With 80% efficacy, we expect ~80% relative reduction
    CHECK(vax_rate < unvax_rate * 0.7);
  }
}

// =============================================================================
// Integration Test 2: Natural Immunity → Re-exposure Resistance
// =============================================================================

TEST_CASE(
    "Integration: Natural immunity reduces susceptibility after recovery") {
  TestFixture fix;
  Disease disease = fix.makeDisease();

  // Simulate recovery: person was infected, recovered, immunity set
  Person& target = fix.world.people[1];
  target.immunity.natural_level = fix.trans.natural_immunity.level;  // 0.95
  target.immunity.natural_acquisition_time = 1.0;
  target.immunity.natural_waning_rate = fix.trans.natural_immunity.waning_rate;

  SUBCASE("Recovered person has very low susceptibility") {
    // Right after recovery at time=2.0
    double susc = target.getSusceptibility(2.0, "TestFlu");
    // natural level ~0.95 → susceptibility ~0.05
    CHECK(susc < 0.1);
    MESSAGE("Post-recovery susceptibility: ", susc);
  }

  SUBCASE("Recovered person has high resistance to reinfection") {
    double rate = measureInfectionRate(fix, disease);
    // With 95% natural immunity, reinfection should be much lower than
    // unvaccinated (~93%) But with very high force of infection, absolute rate
    // can still be non-trivial
    CHECK(rate < 0.40);
    MESSAGE("Reinfection rate for recovered person: ", rate);
  }

  SUBCASE("Natural immunity wanes over time") {
    double susc_early = target.getSusceptibility(2.0, "TestFlu");
    double susc_later =
        target.getSusceptibility(1000.0, "TestFlu");  // ~1000 days later

    // Waning should increase susceptibility over time
    CHECK(susc_later > susc_early);
    MESSAGE("Susceptibility at t=2: ", susc_early, ", at t=1000: ", susc_later);
  }
}

// =============================================================================
// Integration Test 3: Hybrid Immunity (Vaccine + Natural)
// =============================================================================

TEST_CASE("Integration: Hybrid immunity stacks multiplicatively") {
  TestFixture fix;
  Disease disease = fix.makeDisease();

  Person& target = fix.world.people[1];

  // Give natural immunity: 50% level
  target.immunity.natural_level = 0.5;
  target.immunity.natural_acquisition_time = 0.0;
  target.immunity.natural_waning_rate = 0.0;  // No waning for clarity

  // Give vaccine: 60% efficacy
  target.vaccine_trajectory = std::make_unique<VaccineTrajectory>();
  target.vaccine_trajectory->vaccine_name = "TestVax";
  Dose dose;
  dose.number = 0;
  dose.day_administered = 0.0;
  dose.days_to_effective = 0.0;
  dose.days_to_waning = 365.0;
  dose.days_to_finished = 730.0;
  dose.waning_factor = 0.5;
  AgeEfficacy all_ages{0, 120, 0.6};
  dose.infection_efficacy["TestFlu"] = {all_ages};
  target.vaccine_trajectory->addDose(dose);

  SUBCASE("Hybrid susceptibility is multiplicative") {
    // natural_susceptibility = 1.0 - 0.5 = 0.5
    // vaccine_susceptibility = 1.0 - 0.6 = 0.4
    // hybrid = 0.5 * 0.4 = 0.2
    double susc = target.getSusceptibility(5.0, "TestFlu");
    CHECK(susc == doctest::Approx(0.2).epsilon(0.01));
  }

  SUBCASE("Hybrid immunity provides better protection than either alone") {
    // Measure hybrid rate
    double hybrid_rate = measureInfectionRate(fix, disease, 1000);

    // Remove vaccine, measure natural-only
    target.vaccine_trajectory.reset();
    double natural_only_rate = measureInfectionRate(fix, disease, 1000);

    // Remove natural, add vaccine back, measure vaccine-only
    target.immunity.natural_level = 0.0;
    target.vaccine_trajectory = std::make_unique<VaccineTrajectory>();
    target.vaccine_trajectory->vaccine_name = "TestVax";
    Dose dose2;
    dose2.number = 0;
    dose2.day_administered = 0.0;
    dose2.days_to_effective = 0.0;
    dose2.days_to_waning = 365.0;
    dose2.days_to_finished = 730.0;
    dose2.waning_factor = 0.5;
    AgeEfficacy ae{0, 120, 0.6};
    dose2.infection_efficacy["TestFlu"] = {ae};
    target.vaccine_trajectory->addDose(dose2);
    double vax_only_rate = measureInfectionRate(fix, disease, 1000);

    MESSAGE("Hybrid: ", hybrid_rate, ", Natural-only: ", natural_only_rate,
            ", Vax-only: ", vax_only_rate);

    // Hybrid should be <= both individual protections
    CHECK(hybrid_rate <=
          natural_only_rate + 0.05);  // Small tolerance for stochasticity
    CHECK(hybrid_rate <= vax_only_rate + 0.05);
  }
}

// =============================================================================
// Integration Test 4: Seed -> Transmit -> Recover -> Immunity -> Reinfection
// =============================================================================

TEST_CASE(
    "Integration: Full infection lifecycle — seed, transmit, recover, "
    "reinfect") {
  TestFixture fix;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{0.15}};
  fix.cm.default_matrix = default_contact_matrix;

  // Use a trajectory that recovers in 5 days
  fix.trajectories.clear();
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 5.0}}}});  // 5 days mild
  td.stages.push_back(
      {"healthy", {"constant", {{"value", 100.0}}}});  // recover
  fix.trajectories.push_back(td);

  // Mark "healthy" as recovered stage
  DiseaseStageSettings stage_settings;
  stage_settings.recovered_stages = {"healthy"};

  Disease disease("TestFlu", fix.symptom_tags, stage_settings, fix.trajectories,
                  {}, fix.trans);

  // Phase 1: Infect person 0 (the infector)
  fix.world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &fix.world.people[0], 42, nullptr, "office", 0);

  // Phase 2: Transmit to person 1
  {
    auto locs = fix.makeColocatedLocations();
    ParallelConfig parallel_config;
    fix.cm.allow_default_matrix = true;
    finalizeContactMatrices(fix.cm, fix.world);
    InteractionManager im(fix.world, fix.cm, fix.sim_cfg, parallel_config,
                          &disease, nullptr);
    // Run with high enough contacts to guarantee infection
    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{100.0}};
    fix.cm.default_matrix = default_contact_matrix;
    fix.cm.allow_default_matrix = true;
    finalizeContactMatrices(fix.cm, fix.world);
    InteractionManager im_high(fix.world, fix.cm, fix.sim_cfg, parallel_config,
                               &disease, nullptr);
    im_high.processTransmissions(locs, 1.0, 8.0, nullptr);
    REQUIRE(fix.world.people[1].infection != nullptr);
  }

  // Phase 3: Simulate recovery via Epidemiology
  Epidemiology epi(fix.world, &disease);
  epi.trackInfection(1);
  std::vector<PersonLocation> empty_locs;
  epi.updateInfectionStates(
      10.0, empty_locs);  // At t=10, person 1 should be recovered

  // Verify recovery
  CHECK(fix.world.people[1].infection == nullptr);
  CHECK(fix.world.people[1].immunity.natural_level == doctest::Approx(0.95));

  // Phase 4: Attempt reinfection — person 1 should be protected
  // Reset infector
  fix.world.people[0].infection = std::make_unique<Infection>(
      &disease, 11.0, &fix.world.people[0], 99, nullptr, "office", 0);

  double susc_after_recovery =
      fix.world.people[1].getSusceptibility(11.0, "TestFlu");
  CHECK(susc_after_recovery < 0.1);  // ~5% susceptibility

  // Measure reinfection rate with moderate contacts
  ContactMatrix reinfection_contact_matrix;
  reinfection_contact_matrix.bins = {"all"};
  reinfection_contact_matrix.contacts = {{0.15}};
  fix.cm.default_matrix = reinfection_contact_matrix;
  int reinfections = 0;
  int trials = 500;
  for (int trial = 0; trial < trials; ++trial) {
    fix.world.people[1].infection.reset();
    fix.world.people[0].infection =
        std::make_unique<Infection>(&disease, 11.0, &fix.world.people[0],
                                    200 + trial, nullptr, "office", 0);

    auto locs = fix.makeColocatedLocations();
    ParallelConfig parallel_config;
    fix.cm.allow_default_matrix = true;
    finalizeContactMatrices(fix.cm, fix.world);
    InteractionManager im(fix.world, fix.cm, fix.sim_cfg, parallel_config,
                          &disease, nullptr);
    im.processTransmissions(locs, 12.0, 8.0, nullptr);
    if (fix.world.people[1].infection) reinfections++;
  }

  double reinfection_rate = static_cast<double>(reinfections) / trials;
  MESSAGE("Reinfection rate after natural immunity: ", reinfection_rate);

  // With 95% natural immunity, reinfection rate should be much lower than
  // the unvaccinated rate (which is ~0.7 with these settings)
  CHECK(reinfection_rate < 0.20);
}

// =============================================================================
// Integration Test 5: Encounter Logging & Mixing
// =============================================================================

TEST_CASE(
    "Integration: InteractionManager logs physical and virtual encounters") {
  TestFixture fix;
  Disease disease = fix.makeDisease();
  EventLogger logger;

  // Setup 10 people:
  // - Persons 0-3 are in a coordinated encounter at physical Venue 0
  // - Persons 4-6 are regular occupants at physical Venue 0 (Mixing)
  // - Persons 7-9 are in a virtual encounter (Venue -1)

  // Ensure world has enough people
  while (fix.world.people.size() < 10) {
    auto& p = fix.world.people.emplace_back();
    p.id = static_cast<PersonId>(fix.world.people.size() - 1);
  }
  fix.world.buildIndices();

  std::vector<PersonLocation> locs;

  // Group 1: Coordinated Encounter at Physical Venue 0
  for (int i = 0; i < 4; ++i) {
    PersonLocation loc;
    loc.person_id = i;
    loc.venue_id = 0;
    loc.subset_index = -1;
    loc.encounter_type_id = 0;  // "romantic_encounter" or similar index
    loc.person_array_index = i;
    locs.push_back(loc);
  }

  // Group 2: Regular Occupants at Physical Venue 0 (Should mix)
  for (int i = 4; i < 7; ++i) {
    PersonLocation loc;
    loc.person_id = i;
    loc.venue_id = 0;
    loc.subset_index = -1;
    loc.encounter_type_id = 255;  // Not an encounter
    loc.person_array_index = i;
    locs.push_back(loc);
  }

  // Group 3: Coordinated Encounter at Virtual Venue -1
  for (int i = 7; i < 10; ++i) {
    PersonLocation loc;
    loc.person_id = i;
    loc.venue_id = -1;
    loc.subset_index = -1;
    loc.encounter_type_id = 1;  // Different encounter type
    loc.person_array_index = i;
    locs.push_back(loc);
  }

  ParallelConfig parallel_config;
  fix.cm.allow_default_matrix = true;
  finalizeContactMatrices(fix.cm, fix.world);
  InteractionManager im(fix.world, fix.cm, fix.sim_cfg, parallel_config,
                        &disease, &logger);

  SUBCASE("Should log exactly 1 actual encounter per participant") {
    // Run with 1 hour delta
    im.processTransmissions(locs, 8.0, 1.0, nullptr);

    // Group 1 has 4 participants, Group 3 has 3 participants. Total = 7.
    // Since it's a weekday (Day 0 is Monday in our config), we check weekday
    // counter.
    CHECK(logger.getActualWeekdayEncounters() == 7);
    CHECK(logger.getActualWeekendEncounters() == 0);
  }

  SUBCASE("MPI Simulation: Visitors should not be logged locally") {
    std::unordered_set<PersonId> visitor_ids = {0};  // Person 0 is a visitor
    logger.clear();

    // Pass visitor_ids to simulate MPI rank processing
    im.processTransmissions(locs, 8.0, 1.0, nullptr, &visitor_ids);

    // Expected: 7 total participants, but Person 0 is a visitor.
    // So ONLY 6 should be logged on this rank.
    // CURRENTLY this will likely fail and return 7 because we don't check
    // visitor_ids in logging.
    CHECK(logger.getActualWeekdayEncounters() == 6);
  }
}
