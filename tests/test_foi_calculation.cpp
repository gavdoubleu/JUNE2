#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <memory>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/interaction_manager.h"
#include "test_utils.h"

using namespace june;

// =============================================================================
// Force of Infection Mathematical Verification
//
// The FoI formula for a susceptible person in bin j is:
//   lambda = sum_m sum_b [ C_{j,b}^m / N_b * I_b^m * susc_mult_m ]
// Infection probability:
//   P = 1 - exp(-lambda * susceptibility)
//
// These tests verify the math by constructing scenarios with known analytical
// solutions.
// =============================================================================

TEST_CASE("FoI: Single infectious person, single susceptible, single mode") {
  // Setup: 1 infectious + 1 susceptible at same venue
  // Expected: lambda = contacts * I * delta_hours / N_bin
  //   where N_bin = max(1, 2-1) = 1 (same bin, subtract self)
  //   I = constant 1.0 (stage-driven), integrated = 1.0 * delta_hours
  //   lambda = contacts * 1.0 * delta_hours / 1
  //   prob = 1 - exp(-lambda)

  WorldState world;
  world.venue_type_names = {"office"};
  world.geo_level_names = {"city"};
  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  venue.geo_unit_id = -1;
  world.venues.push_back(venue);

  for (int i = 0; i < 2; ++i) {
    auto& p = world.people.emplace_back();
    p.id = i;
    p.age = 30.0f;
    p.sex = Sex::MALE;
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

  double contacts = 0.3;
  double delta_hours = 4.0;
  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{contacts}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim_cfg;
  ParallelConfig parallel_config;

  // Theoretical: I_integrated = curve(t_in_stage) * delta_hours
  //   ConstantCurve(1.0) evaluated over [t_in_stage, t_in_stage + delta/24]
  //   integrateStageDrivenInfectiousness returns 24 * integral(days)
  //   = 24 * 1.0 * (delta_hours/24) = delta_hours
  // lambda = contacts / N_bin * I_integrated
  //   N_bin = max(1, 2-1) = 1 (same bin, self-subtracted)
  //   lambda = 0.3 / 1 * 4.0 = 1.2
  // prob = 1 - exp(-1.2) ≈ 0.6988

  double expected_lambda = contacts * delta_hours;
  double expected_prob = 1.0 - std::exp(-expected_lambda);

  // The InteractionManager uses deterministic per-entity RNG seeded by
  // (base_seed, person_id, venue_id, time_bits). To get proper stochastic
  // sampling, we vary the base_seed across trials.
  int num_trials = 5000;
  int num_infected = 0;

  for (int trial = 0; trial < num_trials; ++trial) {
    world.people[1].infection.reset();
    if (!world.people[0].infection) {
      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);
    }

    SimulationConfig trial_cfg;
    trial_cfg.random_seed = static_cast<uint64_t>(trial);
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, disease);
    InteractionManager im(world, cm, trial_cfg, parallel_config, &disease,
                          nullptr);
    std::vector<PersonLocation> locs;
    locs.push_back({0, 0, -1, 0, 255, 0});
    locs.push_back({1, 0, -1, 0, 255, 1});

    im.processTransmissions(locs, 5.0, delta_hours, nullptr);
    if (world.people[1].infection) num_infected++;
  }

  double empirical_prob = static_cast<double>(num_infected) / num_trials;

  MESSAGE("Expected prob: ", expected_prob, ", Empirical: ", empirical_prob);
  // Allow 5% tolerance for stochastic test (with 5000 trials, std error ≈
  // 0.65%)
  CHECK(empirical_prob == doctest::Approx(expected_prob).epsilon(0.05));
}

TEST_CASE("FoI: Zero susceptibility means zero infections") {
  WorldState world;
  world.venue_type_names = {"office"};
  world.geo_level_names = {"city"};
  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  venue.geo_unit_id = -1;
  world.venues.push_back(venue);

  for (int i = 0; i < 2; ++i) {
    auto& p = world.people.emplace_back();
    p.id = i;
    p.age = 30.0f;
    p.sex = Sex::MALE;
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

  // Give person 1 full natural immunity
  world.people[1].immunity.natural_level = 1.0;
  world.people[1].immunity.natural_acquisition_time = 0.0;
  world.people[1].immunity.natural_waning_rate = 0.0;

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix = default_contact_matrix;  // Very high contacts
  SimulationConfig sim_cfg;
  ParallelConfig parallel_config;

  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 42, nullptr, "office", 0);

  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim_cfg, parallel_config, &disease, nullptr);
  std::vector<PersonLocation> locs;
  locs.push_back({0, 0, -1, 0, 255, 0});
  locs.push_back({1, 0, -1, 0, 255, 1});

  // Full immunity → susceptibility = 0 → should never get infected
  CHECK(world.people[1].getSusceptibility(5.0, "Flu") == doctest::Approx(0.0));
  im.processTransmissions(locs, 5.0, 8.0, nullptr);
  CHECK(world.people[1].infection == nullptr);
}

TEST_CASE("FoI: Higher contacts produce higher infection rate") {
  auto run_with_contacts = [](double contacts, int trials) -> double {
    WorldState world;
    world.venue_type_names = {"office"};
    world.geo_level_names = {"city"};
    Venue venue;
    venue.id = 0;
    venue.type_id = 0;
    venue.geo_unit_id = -1;
    world.venues.push_back(venue);

    for (int i = 0; i < 2; ++i) {
      auto& p = world.people.emplace_back();
      p.id = i;
      p.age = 30.0f;
      p.sex = Sex::MALE;
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

    ContactMatrixConfig cm;
    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{contacts}};
    cm.default_matrix = default_contact_matrix;
    SimulationConfig sim_cfg;
    ParallelConfig parallel_config;

    int infected_count = 0;
    for (int t = 0; t < trials; ++t) {
      world.people[1].infection.reset();
      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], t, nullptr, "office", 0);

      SimulationConfig trial_cfg;
      trial_cfg.random_seed = static_cast<uint64_t>(t);
      ParallelConfig parallel_config;
      cm.allow_default_matrix = true;
      finalizeContactMatrices(cm, world, disease);
      InteractionManager im(world, cm, trial_cfg, parallel_config, &disease,
                            nullptr);
      std::vector<PersonLocation> locs;
      locs.push_back({0, 0, -1, 0, 255, 0});
      locs.push_back({1, 0, -1, 0, 255, 1});
      im.processTransmissions(locs, 5.0, 4.0, nullptr);
      if (world.people[1].infection) infected_count++;
    }
    return static_cast<double>(infected_count) / trials;
  };

  double rate_low = run_with_contacts(0.05, 2000);
  double rate_high = run_with_contacts(0.5, 2000);

  MESSAGE("Low contacts (0.05) rate: ", rate_low,
          ", High contacts (0.5) rate: ", rate_high);
  CHECK(rate_high > rate_low);
}

TEST_CASE("FoI: No infectious people means no infections") {
  WorldState world;
  world.venue_type_names = {"office"};
  world.geo_level_names = {"city"};
  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  venue.geo_unit_id = -1;
  world.venues.push_back(venue);

  for (int i = 0; i < 5; ++i) {
    auto& p = world.people.emplace_back();
    p.id = i;
    p.age = 30.0f;
    p.sex = Sex::MALE;
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

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim_cfg;
  ParallelConfig parallel_config;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim_cfg, parallel_config, &disease, nullptr);

  // Nobody is infected
  std::vector<PersonLocation> locs;
  for (int i = 0; i < 5; ++i) {
    locs.push_back({i, 0, -1, 0, 255, static_cast<size_t>(i)});
  }

  int n = im.processTransmissions(locs, 5.0, 8.0, nullptr);
  CHECK(n == 0);
  for (auto& p : world.people) {
    CHECK(p.infection == nullptr);
  }
}

TEST_CASE("FoI: Empty venue produces no infections") {
  WorldState world;
  world.venue_type_names = {"office"};
  world.geo_level_names = {"city"};
  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  venue.geo_unit_id = -1;
  world.venues.push_back(venue);
  world.buildIndices();

  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto cur = std::make_shared<ConstantCurve>(1.0);
  trans.symptom_id_curves = {nullptr, cur};

  SymptomTag healthy{.name = "healthy", .value = -1, .id = 0};
  SymptomTag sick{.name = "infectious", .value = 1, .id = 1};
  DiseaseStageSettings stage_settings;
  std::vector<TrajectoryDefinition> trajectories;
  Disease disease("Flu", {healthy, sick}, stage_settings, trajectories, {},
                  trans);

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim_cfg;
  ParallelConfig parallel_config;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim_cfg, parallel_config, &disease, nullptr);

  std::vector<PersonLocation> locs;  // Empty
  int n = im.processTransmissions(locs, 5.0, 8.0, nullptr);
  CHECK(n == 0);
}
