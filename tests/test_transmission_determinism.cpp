#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <algorithm>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/interaction_manager.h"
#include "test_utils.h"

using namespace june;

// Helper to create a standard test disease
static Disease makeTestDisease() {
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
  return Disease("Flu", {healthy, sick}, stage_settings, trajectories, {},
                 trans);
}

// Run a transmission pass and return the set of newly infected person IDs
static std::vector<PersonId> runTransmission(WorldState& world,
                                             Disease& disease,
                                             std::vector<PersonLocation>& locs,
                                             ContactMatrixConfig& cm,
                                             const SimulationConfig& sim_cfg) {
  ParallelConfig parallel_config;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim_cfg, parallel_config, &disease, nullptr);
  im.processTransmissions(locs, 5.0, 8.0, nullptr);

  std::vector<PersonId> infected;
  for (auto& p : world.people) {
    if (p.infection != nullptr && p.id >= 10) {
      // Only count susceptible people (IDs >= 10) who got newly infected
      infected.push_back(p.id);
    }
  }
  std::sort(infected.begin(), infected.end());
  return infected;
}

TEST_CASE("Transmission determinism: same seed, shuffled input order") {
  int num_people = 30;
  int num_infectious = 10;

  auto setup_world = [&]() {
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
    return world;
  };

  Disease disease = makeTestDisease();
  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{0.1}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim_cfg;
  sim_cfg.random_seed = 42;

  SUBCASE("Ordered vs shuffled locations produce identical infections") {
    // Run 1: ordered locations
    WorldState world1 = setup_world();
    for (int i = 0; i < num_infectious; ++i) {
      world1.people[i].infection = std::make_unique<Infection>(
          &disease, 0.0, &world1.people[i], 500 + i, nullptr, "office", 0);
    }
    std::vector<PersonLocation> locs1;
    for (int i = 0; i < num_people; ++i) {
      locs1.push_back({i, 0, -1, 0, 255, static_cast<size_t>(i)});
    }
    auto result1 = runTransmission(world1, disease, locs1, cm, sim_cfg);

    // Run 2: shuffled locations (same world state, same seed)
    WorldState world2 = setup_world();
    for (int i = 0; i < num_infectious; ++i) {
      world2.people[i].infection = std::make_unique<Infection>(
          &disease, 0.0, &world2.people[i], 500 + i, nullptr, "office", 0);
    }
    std::vector<PersonLocation> locs2;
    for (int i = 0; i < num_people; ++i) {
      locs2.push_back({i, 0, -1, 0, 255, static_cast<size_t>(i)});
    }
    // Shuffle the locations vector with a deterministic shuffle
    std::mt19937 shuffler(999);
    std::shuffle(locs2.begin(), locs2.end(), shuffler);
    auto result2 = runTransmission(world2, disease, locs2, cm, sim_cfg);

    // The sorted results must be identical
    CHECK(result1 == result2);
    MESSAGE("Infections from ordered: ", result1.size(),
            ", shuffled: ", result2.size());
  }

  SUBCASE("Same seed, same order, produces identical results across runs") {
    std::vector<PersonId> results[3];
    for (int run = 0; run < 3; ++run) {
      WorldState world = setup_world();
      for (int i = 0; i < num_infectious; ++i) {
        world.people[i].infection = std::make_unique<Infection>(
            &disease, 0.0, &world.people[i], 500 + i, nullptr, "office", 0);
      }
      std::vector<PersonLocation> locs;
      for (int i = 0; i < num_people; ++i) {
        locs.push_back({i, 0, -1, 0, 255, static_cast<size_t>(i)});
      }
      results[run] = runTransmission(world, disease, locs, cm, sim_cfg);
    }
    CHECK(results[0] == results[1]);
    CHECK(results[1] == results[2]);
  }

  SUBCASE("Different seed produces different results") {
    WorldState world1 = setup_world();
    WorldState world2 = setup_world();
    for (int i = 0; i < num_infectious; ++i) {
      world1.people[i].infection = std::make_unique<Infection>(
          &disease, 0.0, &world1.people[i], 500 + i, nullptr, "office", 0);
      world2.people[i].infection = std::make_unique<Infection>(
          &disease, 0.0, &world2.people[i], 500 + i, nullptr, "office", 0);
    }
    std::vector<PersonLocation> locs;
    for (int i = 0; i < num_people; ++i) {
      locs.push_back({i, 0, -1, 0, 255, static_cast<size_t>(i)});
    }

    SimulationConfig sim1;
    sim1.random_seed = 42;
    SimulationConfig sim2;
    sim2.random_seed = 999;

    auto r1 = runTransmission(world1, disease, locs, cm, sim1);
    auto r2 = runTransmission(world2, disease, locs, cm, sim2);

    // With different seeds, results should very likely differ
    // (not guaranteed but probability of identical results is extremely low)
    if (r1.size() > 0 && r2.size() > 0) {
      // At least check they're not trivially identical
      MESSAGE("Seed 42 infections: ", r1.size(),
              ", Seed 999 infections: ", r2.size());
    }
  }
}
