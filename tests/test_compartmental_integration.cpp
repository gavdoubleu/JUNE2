#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <memory>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/infectiousness_curves.h"
#include "epidemiology/interaction_manager.h"
#include "mock_compartmental_model.h"
#include "simulation/compartmental_model_manager.h"
#include "test_utils.h"

using namespace june;

// =============================================================================
// Helpers
// =============================================================================

// Disease with one compartmental_uptake mode AND one compartmental_deposition
// mode.
static Disease makeIntegrationDisease(double susc_mult = 1.5,
                                      double deposition_rate = 0.1) {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;

  std::vector<SymptomTag> tags = {{"healthy", -1, 0}, {"infectious", 1, 1}};

  tp.symptom_id_curves = {nullptr, nullptr};

  // Mode 0: compartmental_uptake
  TransmissionMode uptake_mode;
  uptake_mode.name = "comp_uptake";
  uptake_mode.type = TransmissionModeType::CompartmentalUptake;
  uptake_mode.mode_transmissibility_multiplier = susc_mult;
  uptake_mode.symptom_curves = {nullptr, nullptr};
  CompartmentalUptakeConfig ucfg;
  ucfg.mode_index = 0;
  uptake_mode.config = ucfg;
  tp.modes.push_back(std::move(uptake_mode));

  // Mode 1: compartmental_deposition
  TransmissionMode dep_mode;
  dep_mode.name = "comp_dep";
  dep_mode.type = TransmissionModeType::CompartmentalDeposition;
  dep_mode.symptom_curves = {nullptr, nullptr};
  CompartmentalDepositionConfig dcfg;
  dcfg.mode_index = 1;
  dcfg.deposition_by_symptom.resize(2, nullptr);
  dcfg.deposition_by_symptom[1] =
      std::make_shared<ConstantCurve>(deposition_rate);
  dep_mode.config = std::move(dcfg);
  tp.modes.push_back(std::move(dep_mode));

  tp.natural_immunity.level = 0.95;
  tp.natural_immunity.waning_rate = 0.001;

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"infectious", {"constant", {{"value", 100.0}}}});
  trajectories.push_back(td);

  return Disease("IntegDisease", tags, {}, trajectories, {}, tp);
}

// Manager: node 0 owns venue_id; coupling output = coupling_output_value.
struct IntegrationManager {
  InProcessMockPlugin plugin;
  std::unique_ptr<CompartmentalModelManager> mgr;

  explicit IntegrationManager(int venue_id, float coupling_output,
                              float coupling_beta = 0.5f) {
    plugin.coupling_output_value = coupling_output;

    CompartmentalModelSteps steps;
    steps.loadSidecar =
        [coupling_beta](const std::string&) -> PluginSidecarConfig {
      PluginSidecarConfig cfg;
      cfg.plugin_so_path = "ignored.so";
      cfg.default_human_to_compartmental_model_input = coupling_beta;
      return cfg;
    };
    steps.loadPlugin = [this](const std::string&,
                              DestroyCompartmentalModelFn& out_destroy,
                              void*& out_handle) -> ICompartmentalModel* {
      out_destroy = [](ICompartmentalModel*) {};
      out_handle = nullptr;
      return &plugin;
    };
    steps.buildPartition =
        [venue_id](const std::string&, DomainManager*,
                   std::unordered_map<int, int>& venue_map) -> NodePartition {
      NodePartition p;
      p.owned_node_indices = {0};
      p.node_venue_ids = {{venue_id}};
      venue_map[venue_id] = 0;
      return p;
    };
    mgr = std::make_unique<CompartmentalModelManager>("any.yaml", nullptr,
                                                      std::move(steps));
  }
};

// Run one slot of the bidirectional coupling sequence manually.
// Returns number of new infections from processTransmissions.
static int runOneSlot(InteractionManager& im, CompartmentalModelManager& mgr,
                      WorldState& world, const Disease& disease,
                      const std::vector<PersonLocation>& locs, double t,
                      double delta_hours) {
  const double dt_days = delta_hours / 24.0;

  mgr.advance(static_cast<float>(dt_days), static_cast<float>(t));

  const CompartmentalModelManager* comp_model = mgr.isActive() ? &mgr : nullptr;
  int new_infections = im.processTransmissions(
      locs, t, delta_hours, nullptr, nullptr, nullptr, nullptr, comp_model);

  mgr.computeDepositionWriteback(locs, world, disease, t, t + dt_days);

  mgr.maybeSnapshot(static_cast<float>(t));

  return new_infections;
}

// =============================================================================
// Tests
// =============================================================================

TEST_CASE("Integration: per-slot call sequence is correct over 3 slots") {
  const int n_slots = 3;
  const double delta_hours = 1.0;

  WorldState world;
  Venue v;
  v.id = 10;
  v.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(v);
  for (int i = 0; i < 2; ++i) {
    auto& p = world.people.emplace_back();
    p.id = i;
    p.age = 30;
    p.sex = Sex::MALE;
    p.geo_unit_id = -1;
  }
  world.buildIndices();

  Disease disease = makeIntegrationDisease();
  IntegrationManager im_mgr(10, /*coupling_output=*/0.0f);
  REQUIRE(im_mgr.mgr->isActive());

  ContactMatrixConfig cm;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease);

  std::vector<PersonLocation> locs = {
      {0, 10, -1, -1, 255, 0},
      {1, 10, -1, -1, 255, 1},
  };

  for (int s = 0; s < n_slots; ++s)
    runOneSlot(im, *im_mgr.mgr, world, disease, locs, s * delta_hours / 24.0,
               delta_hours);

  CHECK(im_mgr.plugin.advance_count == n_slots);
  CHECK(im_mgr.plugin.snapshot_count == n_slots);
  CHECK(im_mgr.plugin.coupling_inputs_call_count == n_slots);
}

TEST_CASE(
    "Integration: non-zero coupling → susceptibles infected (uptake path)") {
  // lambda = coupling_output * coupling_beta * susc_mult
  //        = 5.0 * 0.5 * 1.5 = 3.75
  // prob = 1 - exp(-3.75) ≈ 0.976 — extremely likely over a few trials
  const float coupling_output = 5.0f;
  const float coupling_beta = 0.5f;
  const double susc_mult = 1.5;

  WorldState world;
  Venue v;
  v.id = 10;
  v.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(v);
  auto& p0 = world.people.emplace_back();
  p0.id = 0;
  p0.age = 30;
  p0.sex = Sex::MALE;
  p0.geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeIntegrationDisease(susc_mult);
  IntegrationManager im_mgr(10, coupling_output, coupling_beta);
  REQUIRE(im_mgr.mgr->isActive());

  // Run up to 20 seeds — with prob ≈ 0.976 per slot, almost certain to infect.
  bool got_infected = false;
  for (int seed = 1; seed <= 20 && !got_infected; ++seed) {
    world.people[0].infection.reset();

    ContactMatrixConfig cm;
    SimulationConfig sim;
    sim.random_seed = static_cast<uint64_t>(seed);
    ParallelConfig par;
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, disease);
    InteractionManager im(world, cm, sim, par, &disease);

    std::vector<PersonLocation> locs = {{0, 10, -1, -1, 255, 0}};
    runOneSlot(im, *im_mgr.mgr, world, disease, locs, 0.0, 1.0);

    if (world.people[0].infection) got_infected = true;
  }
  CHECK(got_infected);
}

TEST_CASE("Integration: infected person → non-zero deposition write-back") {
  const double deposition_rate = 0.1;
  const double delta_days = 1.0 / 24.0;
  const float coupling_beta = 0.5f;
  const double expected_dep =
      24.0 * deposition_rate * delta_days * coupling_beta;

  WorldState world;
  Venue v;
  v.id = 10;
  v.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(v);
  auto& p = world.people.emplace_back();
  p.id = 0;
  p.age = 30;
  p.sex = Sex::MALE;
  p.geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeIntegrationDisease(1.5, deposition_rate);
  IntegrationManager im_mgr(10, 0.0f, coupling_beta);
  REQUIRE(im_mgr.mgr->isActive());

  // Seed the person as infected at t=0
  world.people[0].infection =
      std::make_unique<Infection>(&disease, 0.0, &world.people[0], 42);

  ContactMatrixConfig cm;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease);

  std::vector<PersonLocation> locs = {{0, 10, -1, -1, 255, 0}};
  runOneSlot(im, *im_mgr.mgr, world, disease, locs, 0.0, 1.0);

  REQUIRE(im_mgr.plugin.coupling_inputs_call_count == 1);
  REQUIRE(im_mgr.plugin.last_coupling_inputs.size() == 1);
  CHECK(im_mgr.plugin.last_coupling_inputs[0] ==
        doctest::Approx(static_cast<float>(expected_dep)).epsilon(0.01));
}

TEST_CASE(
    "Integration: zero coupling → no compartmental infections, deposition "
    "still runs") {
  WorldState world;
  Venue v;
  v.id = 10;
  v.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(v);
  auto& p = world.people.emplace_back();
  p.id = 0;
  p.age = 30;
  p.sex = Sex::MALE;
  p.geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeIntegrationDisease();
  IntegrationManager im_mgr(10, /*coupling_output=*/0.0f);

  ContactMatrixConfig cm;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease);

  std::vector<PersonLocation> locs = {{0, 10, -1, -1, 255, 0}};

  for (int seed = 1; seed <= 50; ++seed) {
    world.people[0].infection.reset();
    sim.random_seed = static_cast<uint64_t>(seed);
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, disease);
    InteractionManager im2(world, cm, sim, par, &disease);
    runOneSlot(im2, *im_mgr.mgr, world, disease, locs, 0.0, 1.0);
    CHECK(world.people[0].infection == nullptr);
  }

  // writeCouplingInputs should still have been called each slot
  CHECK(im_mgr.plugin.coupling_inputs_call_count == 50);
  // No infected person → zero deposition
  if (!im_mgr.plugin.last_coupling_inputs.empty())
    CHECK(im_mgr.plugin.last_coupling_inputs[0] == doctest::Approx(0.0f));
}

TEST_CASE(
    "Integration: inactive manager → no coupling calls, processTransmissions "
    "safe") {
  WorldState world;
  Venue v;
  v.id = 10;
  v.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(v);
  auto& p = world.people.emplace_back();
  p.id = 0;
  p.age = 30;
  p.sex = Sex::MALE;
  p.geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeIntegrationDisease();
  CompartmentalModelManager inactive_mgr("", nullptr);
  REQUIRE_FALSE(inactive_mgr.isActive());

  ContactMatrixConfig cm;
  SimulationConfig sim;
  ParallelConfig par;
  std::vector<PersonLocation> locs = {{0, 10, -1, -1, 255, 0}};

  for (int s = 0; s < 3; ++s) {
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, disease);
    InteractionManager im(world, cm, sim, par, &disease);
    REQUIRE_NOTHROW(runOneSlot(im, inactive_mgr, world, disease, locs,
                               s * 1.0 / 24.0, 1.0));
  }
  // Inactive manager → never infected by compartmental source
  CHECK(world.people[0].infection == nullptr);
}

// =============================================================================
// Slice 6: different venue types produce different FOI from identical node
// output
// =============================================================================

TEST_CASE(
    "Integration: per-venue-type foi_scale produces different infection "
    "probability") {
  // household foi_scale=0.5, church foi_scale=0.1 — same node output
  // Expected: household infection prob > church infection prob
  const float node_output = 4.0f;
  const double susc_mult = 1.5;

  auto count_infections = [&](const std::string& venue_type_name,
                              float foi_scale, int num_trials) -> int {
    InProcessMockPlugin plugin;
    plugin.coupling_output_value = node_output;

    CompartmentalModelSteps steps;
    steps.loadSidecar = [foi_scale, &venue_type_name](
                            const std::string&) -> PluginSidecarConfig {
      PluginSidecarConfig cfg;
      cfg.plugin_so_path = "ignored.so";
      CouplingMatrix cm;
      cm.values = {foi_scale};
      cfg.output_foi_matrix.matrices[venue_type_name] = cm;
      cfg.output_foi_matrix.default_value = 1.0f;
      return cfg;
    };
    steps.loadPlugin = [&plugin](const std::string&,
                                 DestroyCompartmentalModelFn& d,
                                 void*& h) -> ICompartmentalModel* {
      d = [](ICompartmentalModel*) {};
      h = nullptr;
      return &plugin;
    };
    steps.buildPartition =
        [](const std::string&, DomainManager*,
           std::unordered_map<int, int>& vm) -> NodePartition {
      NodePartition p;
      p.owned_node_indices = {0};
      p.node_venue_ids = {{10}};
      vm[10] = 0;
      return p;
    };

    auto mgr = std::make_unique<CompartmentalModelManager>("any.yaml", nullptr,
                                                           std::move(steps));
    mgr->resolveVenueTypes({venue_type_name});

    WorldState world;
    Venue v;
    v.id = 10;
    v.type_id = 0;
    world.venue_type_names = {venue_type_name};
    world.venues.push_back(v);
    world.people.emplace_back();
    world.people[0].id = 0;
    world.people[0].age = 30;
    world.people[0].sex = Sex::MALE;
    world.people[0].geo_unit_id = -1;
    world.buildIndices();

    Disease disease = makeIntegrationDisease(susc_mult);
    int infected = 0;
    for (int trial = 0; trial < num_trials; ++trial) {
      world.people[0].infection.reset();
      std::vector<PersonLocation> locs = {{0, 10, -1, -1, 255, 0}};
      ContactMatrixConfig cm;
      SimulationConfig sim;
      sim.random_seed = static_cast<uint64_t>(trial + 1);
      ParallelConfig par;
      cm.allow_default_matrix = true;
      finalizeContactMatrices(cm, world, disease);
      InteractionManager im(world, cm, sim, par, &disease);
      runOneSlot(im, *mgr, world, disease, locs, 0.0, 1.0);
      if (world.people[0].infection) ++infected;
    }
    return infected;
  };

  const int n = 4000;
  int hh_infected = count_infections("household", 0.5f, n);
  int church_infected = count_infections("church", 0.1f, n);

  // household: lambda = 4.0 * 0.5 * 1.5 = 3.0  → prob ≈ 0.95
  // church:    lambda = 4.0 * 0.1 * 1.5 = 0.6  → prob ≈ 0.45
  // Household should have significantly more infections (at least 30% more)
  CHECK(hh_infected > church_infected);
  double hh_prob = static_cast<double>(hh_infected) / n;
  double church_prob = static_cast<double>(church_infected) / n;
  CHECK(hh_prob > church_prob + 0.3);
}
