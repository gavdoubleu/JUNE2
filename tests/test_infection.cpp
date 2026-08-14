#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "core/config.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/interaction_manager.h"
#include "test_utils.h"

using namespace june;

TEST_CASE("Disease trajectory progression") {
  // 1. Setup simple disease
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;

  auto cur = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["mild"] = cur;
  trans.symptom_id_curves = {nullptr, cur};

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.description = "Test trajectory";
  td.selection_key = "general_population";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 5.0}}}});  // 5 days mild
  td.stages.push_back(
      {"healthy", {"constant", {{"value", 100.0}}}});  // Transition to healthy
  trajectories.push_back(td);

  SymptomTag healthy;
  healthy.name = "healthy";
  healthy.id = 0;
  healthy.value = -1;
  SymptomTag mild;
  mild.name = "mild";
  mild.id = 1;
  mild.value = 1;

  DiseaseStageSettings stage_settings;
  stage_settings.recovered_stages = {"healthy"};

  Disease disease("TestFlu", {healthy, mild}, stage_settings, trajectories, {},
                  trans);

  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 123, nullptr, "office", 0);

  SUBCASE("Initial stage") {
    CHECK(p.infection->getCurrentSymptom(0.0) == "mild");
    CHECK(p.infection->getTrajectory().getCurrentSymptomId(0.0) == 1);
  }

  SUBCASE("Recovery") {
    // Trajectory: at t=0 mild(1) for 5 days. at t=5 should transition to
    // healthy(0).
    CHECK(p.infection->getTrajectory().getCurrentSymptomId(6.0) == 0);
    CHECK(p.infection->isRecovered(6.0) == true);
  }
}

TEST_CASE("Sentinel venue_id=-1 filtered by processTransmissions") {
  WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
  Venue& venue = world.venues[0];
  venue.type_id = 0;

  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto cur = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["infectious"] = cur;
  trans.symptom_id_curves = {nullptr, cur};

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general_population";
  td.stages.push_back({"infectious", {"constant", {{"value", 10.0}}}});
  td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});
  trajectories.push_back(td);

  SymptomTag healthy{.name = "healthy", .value = -1, .id = 0};
  SymptomTag sick{.name = "infectious", .value = 1, .id = 1};
  DiseaseStageSettings stage_settings;
  Disease disease("Flu", {healthy, sick}, stage_settings, trajectories, {},
                  trans);

  ContactMatrixConfig cm_config;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm_config.default_matrix = default_contact_matrix;
  SimulationConfig sim_config;
  ParallelConfig parallel_config;
  cm_config.allow_default_matrix = true;
  finalizeContactMatrices(cm_config, world, disease);
  InteractionManager im(world, cm_config, sim_config, parallel_config, &disease,
                        nullptr);

  // Infect Person 0
  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 123, nullptr, "office", 0);

  SUBCASE("Sentinel locations are skipped") {
    // Both people at venue_id=-1 with encounter_type_id=255 (unallocated)
    std::vector<PersonLocation> locs;
    locs.push_back({0, -1, -1, 0, 255, 0});
    locs.push_back({1, -1, -1, 0, 255, 1});

    int infections = im.processTransmissions(locs, 1.0, 1.0, nullptr);
    CHECK(infections == 0);
    CHECK(world.people[1].infection == nullptr);
  }

  SUBCASE("Mixed: sentinel skipped, real venue processed") {
    std::vector<PersonLocation> locs;
    // Person 0 at real venue
    locs.push_back({0, 0, -1, 0, 255, 0});
    // Person 1 at sentinel (should be filtered out)
    locs.push_back({1, -1, -1, 0, 255, 1});

    int infections = im.processTransmissions(locs, 1.0, 1.0, nullptr);
    // Person 1 is at sentinel venue, not co-located with Person 0
    CHECK(infections == 0);
    CHECK(world.people[1].infection == nullptr);
  }
}

TEST_CASE("Large venue transmission statistical correctness") {
  // Create a venue with 100 people, 50% infectious
  int num_people = 100;
  int num_infectious = 50;

  WorldState world;
  world.venue_type_names = {"office"};
  world.geo_level_names = {"city"};

  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  venue.geo_unit_id = -1;
  world.venues.push_back(venue);

  for (int i = 0; i < num_people; ++i) {
    auto& p = world.people.emplace_back();
    p.id = i;
    p.age = 30.0f;
    p.sex = (i % 2 == 0) ? Sex::MALE : Sex::FEMALE;
    p.geo_unit_id = -1;
  }
  world.buildIndices();

  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto cur = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["infectious"] = cur;
  trans.symptom_id_curves = {nullptr, cur};

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general_population";
  td.stages.push_back({"infectious", {"constant", {{"value", 100.0}}}});
  td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});
  trajectories.push_back(td);

  SymptomTag healthy{.name = "healthy", .value = -1, .id = 0};
  SymptomTag sick{.name = "infectious", .value = 1, .id = 1};
  DiseaseStageSettings stage_settings;
  Disease disease("Flu", {healthy, sick}, stage_settings, trajectories, {},
                  trans);

  // Moderate contacts so infection rate is neither 0 nor 100%
  ContactMatrixConfig cm_config;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{0.05}};
  cm_config.default_matrix = default_contact_matrix;
  SimulationConfig sim_config;
  ParallelConfig parallel_config;

  // Run multiple trials
  int total_infections = 0;
  int num_trials = 50;
  int num_susceptible = num_people - num_infectious;

  for (int trial = 0; trial < num_trials; ++trial) {
    // Reset infections
    for (auto& p : world.people) p.infection.reset();

    // Infect first num_infectious people
    for (int i = 0; i < num_infectious; ++i) {
      world.people[i].infection =
          std::make_unique<Infection>(&disease, 0.0, &world.people[i],
                                      trial * 1000 + i, nullptr, "office", 0);
    }

    cm_config.allow_default_matrix = true;
    finalizeContactMatrices(cm_config, world, disease);
    InteractionManager im(world, cm_config, sim_config, parallel_config,
                          &disease, nullptr);

    std::vector<PersonLocation> locs;
    for (int i = 0; i < num_people; ++i) {
      locs.push_back({i, 0, -1, 0, 255, static_cast<size_t>(i)});
    }

    int n = im.processTransmissions(locs, 5.0, 8.0, nullptr);
    total_infections += n;
  }

  double avg = static_cast<double>(total_infections) / num_trials;
  // With 50 infectious, 50 susceptible, contacts=0.05, delta=8h, I=1.0:
  // omega per susceptible ≈ 0.05 * 8 * 50 / 99 ≈ 0.202
  // prob ≈ 1 - exp(-0.202) ≈ 0.183
  // Expected infections per trial ≈ 50 * 0.183 ≈ 9.15
  // Allow wide tolerance for stochastic test
  CHECK(avg > 2.0);
  CHECK(avg < 30.0);
  // Verify not all susceptible got infected (would indicate a bug)
  CHECK(avg < num_susceptible);
  MESSAGE("Average infections per trial (100 people, 50 infectious): ", avg);
}

TEST_CASE("InteractionManager basic transmission") {
  WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
  Venue& venue = world.venues[0];
  venue.type_id = 0;  // office

  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto cur = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["infectious"] = cur;
  trans.symptom_id_curves = {nullptr, cur};

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.description = "Test trajectory";
  td.selection_key = "general_population";
  td.stages.push_back({"infectious", {"constant", {{"value", 10.0}}}});
  td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});
  trajectories.push_back(td);

  SymptomTag healthy;
  healthy.id = 0;
  healthy.name = "healthy";
  healthy.value = -1;
  SymptomTag sick;
  sick.id = 1;
  sick.name = "infectious";
  sick.value = 1;

  DiseaseStageSettings stage_settings;
  Disease disease("Flu", {healthy, sick}, stage_settings, trajectories, {},
                  trans);

  ContactMatrixConfig cm_config;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm_config.default_matrix =
      default_contact_matrix;  // High contacts for deterministic transmit
  SimulationConfig sim_config;

  ParallelConfig parallel_config;
  cm_config.allow_default_matrix = true;
  finalizeContactMatrices(cm_config, world, disease);
  InteractionManager im(world, cm_config, sim_config, parallel_config, &disease,
                        nullptr);

  // Infect Person 0
  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 123, nullptr, "office", 0);

  // Both at same venue
  std::vector<PersonLocation> locs;
  locs.push_back({0, 0, -1, 0, 255, 0});
  locs.push_back({1, 0, -1, 0, 255, 1});

  SUBCASE("Successful transmission") {
    im.processTransmissions(locs, 1.0, 1.0, nullptr);
    CHECK(world.people[1].infection != nullptr);
  }
}
