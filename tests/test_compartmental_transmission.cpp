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

// ---------------------------------------------------------------------------
// Build a disease with one standard mode and one compartmental_uptake mode.
// ---------------------------------------------------------------------------
static Disease makeDiseaseWithCompartmentalUptake(double susc_mult = 1.5) {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;

  std::vector<SymptomTag> tags = {{"healthy", -1, 0}, {"infectious", 1, 1}};

  auto direct_curve = std::make_shared<ConstantCurve>(1.0);
  tp.symptom_id_curves = {nullptr, direct_curve};
  tp.stage_curves["infectious"] = direct_curve;

  TransmissionMode direct;
  direct.name = "direct";
  direct.symptom_curves = {nullptr, direct_curve};
  tp.modes.push_back(std::move(direct));

  TransmissionMode uptake;
  uptake.name = "comp_uptake";
  uptake.type = TransmissionModeType::CompartmentalUptake;
  uptake.mode_transmissibility_multiplier = susc_mult;
  uptake.symptom_curves = {nullptr, nullptr};
  CompartmentalUptakeConfig ucfg;
  ucfg.mode_index = 1;
  uptake.config = ucfg;
  tp.modes.push_back(std::move(uptake));

  tp.natural_immunity.level = 0.95;
  tp.natural_immunity.waning_rate = 0.001;

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"infectious", {"constant", {{"value", 100.0}}}});
  trajectories.push_back(td);

  return Disease("TestComp", tags, {}, trajectories, {}, tp);
}

// ---------------------------------------------------------------------------
// Build a minimal manager via step injection.
//   plugin: InProcessMockPlugin owned by the caller (must outlive the manager).
//   venue_id: the single venue to map to node 0.
//   coupling_output: value the plugin returns from getCouplingOutputs.
// ---------------------------------------------------------------------------
static CompartmentalModelManager makeManager(InProcessMockPlugin& plugin,
                                             VenueId venue_id,
                                             float coupling_output = 0.0f) {
  plugin.coupling_output_value = coupling_output;

  CompartmentalModelSteps steps;
  steps.loadSidecar = [](const std::string&) -> PluginSidecarConfig {
    PluginSidecarConfig cfg;
    cfg.plugin_so_path = "ignored.so";
    return cfg;
  };
  steps.loadPlugin = [&plugin](const std::string&,
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
    p.node_venue_ids = {{static_cast<int>(venue_id)}};
    venue_map[venue_id] = 0;
    return p;
  };

  return CompartmentalModelManager("any.yaml", nullptr, std::move(steps));
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("InfectionSource::Compartmental is distinct from Person and Fomite") {
  CHECK(InfectionSource::Compartmental != InfectionSource::Person);
  CHECK(InfectionSource::Compartmental != InfectionSource::Fomite);
}

TEST_CASE("Compartmental uptake: zero coupling → no compartmental infections") {
  WorldState world;
  Venue venue;
  venue.id = 10;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);
  world.people.emplace_back();
  world.people[0].id = 0;
  world.people[0].age = 30;
  world.people[0].sex = Sex::MALE;
  world.people[0].geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeDiseaseWithCompartmentalUptake();
  ContactMatrixConfig cm;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease);

  std::vector<PersonLocation> locs;
  PersonLocation loc;
  loc.person_id = 0;
  loc.venue_id = 10;
  loc.subset_index = -1;
  loc.activity_index = -1;
  loc.encounter_type_id = 255;
  loc.person_array_index = 0;
  locs.push_back(loc);

  InProcessMockPlugin mock_plugin;
  auto mgr = makeManager(mock_plugin, 10, /*coupling_output=*/0.0f);
  REQUIRE(mgr.isActive());
  CHECK(mgr.venueToLocalNodeIndex(10) == 0);

  int new_infections = im.processTransmissions(locs, 5.0, 1.0, nullptr, nullptr,
                                               nullptr, nullptr, &mgr);
  CHECK(new_infections == 0);
  CHECK(world.people[0].infection == nullptr);
}

TEST_CASE("Compartmental uptake: lambda = coupling_output * susc_mult") {
  // lambda = coupling_output * susc_mult
  // The coupling matrix (beta) is the plugin's responsibility to pre-scale.
  // prob = 1 - exp(-lambda * susceptibility)
  // lambda = node_output * foi_scale(default=1.0) * susc_mult
  const float coupling_output = 2.0f;
  const double susc_mult = 1.5;
  const double expected_lambda = coupling_output * 1.0 * susc_mult;
  const double expected_prob = 1.0 - std::exp(-expected_lambda);

  WorldState world;
  Venue venue;
  venue.id = 10;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);
  world.people.emplace_back();
  world.people[0].id = 0;
  world.people[0].age = 30;
  world.people[0].sex = Sex::MALE;
  world.people[0].geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeDiseaseWithCompartmentalUptake(susc_mult);
  InProcessMockPlugin mock_plugin;
  auto mgr = makeManager(mock_plugin, 10, coupling_output);
  REQUIRE(mgr.isActive());

  ContactMatrixConfig cm;
  SimulationConfig sim;
  ParallelConfig par;

  // Monte Carlo over many trials to verify probability
  int num_trials = 8000;
  int num_infected = 0;
  for (int trial = 0; trial < num_trials; ++trial) {
    world.people[0].infection.reset();

    std::vector<PersonLocation> locs;
    PersonLocation loc;
    loc.person_id = 0;
    loc.venue_id = 10;
    loc.subset_index = -1;
    loc.activity_index = -1;
    loc.encounter_type_id = 255;
    loc.person_array_index = 0;
    locs.push_back(loc);

    sim.random_seed = static_cast<uint64_t>(trial + 1);
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, disease);
    InteractionManager im(world, cm, sim, par, &disease);
    mgr.advance(0.0f, 0.0f);  // invalidate buffer cache each trial
    im.processTransmissions(locs, 5.0, 1.0, nullptr, nullptr, nullptr, nullptr,
                            &mgr);

    if (world.people[0].infection) ++num_infected;
  }

  double observed_prob = static_cast<double>(num_infected) / num_trials;
  // Allow ±3% absolute tolerance (Monte Carlo noise)
  CHECK(std::abs(observed_prob - expected_prob) < 0.03);
}

TEST_CASE(
    "Compartmental uptake: infection source is "
    "InfectionSource::Compartmental") {
  WorldState world;
  Venue venue;
  venue.id = 10;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);
  world.people.emplace_back();
  world.people[0].id = 0;
  world.people[0].age = 30;
  world.people[0].sex = Sex::MALE;
  world.people[0].geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeDiseaseWithCompartmentalUptake();
  InProcessMockPlugin mock_plugin;
  auto mgr = makeManager(mock_plugin, 10, /*coupling_output=*/1000.0f);

  ContactMatrixConfig cm;
  SimulationConfig sim;
  ParallelConfig par;

  // Try up to 1000 seeds until we get an infection (coupling=1000 → near
  // certain)
  bool found_compartmental = false;
  for (int seed = 1; seed <= 1000 && !found_compartmental; ++seed) {
    world.people[0].infection.reset();
    std::vector<PersonLocation> locs;
    PersonLocation loc;
    loc.person_id = 0;
    loc.venue_id = 10;
    loc.subset_index = -1;
    loc.activity_index = -1;
    loc.encounter_type_id = 255;
    loc.person_array_index = 0;
    locs.push_back(loc);

    sim.random_seed = static_cast<uint64_t>(seed);
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, disease);
    InteractionManager im(world, cm, sim, par, &disease);
    mgr.advance(0.0f, 0.0f);
    im.processTransmissions(locs, 5.0, 1.0, nullptr, nullptr, nullptr, nullptr,
                            &mgr);

    if (world.people[0].infection) {
      found_compartmental = true;
    }
  }
  CHECK(found_compartmental);
}

TEST_CASE("Compartmental uptake: venue not in node → zero lambda") {
  WorldState world;
  Venue venue;
  venue.id = 10;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);
  world.people.emplace_back();
  world.people[0].id = 0;
  world.people[0].age = 30;
  world.people[0].sex = Sex::MALE;
  world.people[0].geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeDiseaseWithCompartmentalUptake();
  // Manager maps venue 99 (not venue 10) → node 0; coupling=1000 but shouldn't
  // matter
  InProcessMockPlugin mock_plugin;
  auto mgr = makeManager(mock_plugin, 99, /*coupling_output=*/1000.0f);
  REQUIRE(mgr.venueToLocalNodeIndex(10) == -1);

  ContactMatrixConfig cm;
  SimulationConfig sim;
  ParallelConfig par;

  int infections = 0;
  for (int seed = 1; seed <= 200; ++seed) {
    world.people[0].infection.reset();
    std::vector<PersonLocation> locs;
    PersonLocation loc;
    loc.person_id = 0;
    loc.venue_id = 10;
    loc.subset_index = -1;
    loc.activity_index = -1;
    loc.encounter_type_id = 255;
    loc.person_array_index = 0;
    locs.push_back(loc);

    sim.random_seed = static_cast<uint64_t>(seed);
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, disease);
    InteractionManager im(world, cm, sim, par, &disease);
    mgr.advance(0.0f, 0.0f);
    im.processTransmissions(locs, 5.0, 1.0, nullptr, nullptr, nullptr, nullptr,
                            &mgr);
    if (world.people[0].infection) ++infections;
  }
  // No node mapped → no compartmental infections; no direct infectors either
  CHECK(infections == 0);
}

// =============================================================================
// Slice 4: output_foi_scale applied — halved scale gives proportionally lower
// lambda
// =============================================================================

TEST_CASE(
    "Compartmental uptake: output_foi_scale=0.5 halves effective lambda") {
  // With foi_scale=1.0: lambda = node_output * 1.0 * susc_mult
  // With foi_scale=0.5: lambda = node_output * 0.5 * susc_mult
  const float node_output = 3.0f;
  const double susc_mult = 1.0;

  auto measure_prob = [&](float foi_scale) -> double {
    InProcessMockPlugin mock_plugin;
    CompartmentalModelSteps steps;
    steps.loadSidecar = [foi_scale](const std::string&) -> PluginSidecarConfig {
      PluginSidecarConfig cfg;
      cfg.plugin_so_path = "ignored.so";
      cfg.output_foi_matrix.default_value = foi_scale;
      return cfg;
    };
    steps.loadPlugin = [&](const std::string&,
                           DestroyCompartmentalModelFn& out_destroy,
                           void*& out_handle) -> ICompartmentalModel* {
      out_destroy = [](ICompartmentalModel*) {};
      out_handle = nullptr;
      return &mock_plugin;
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
    mock_plugin.coupling_output_value = node_output;

    CompartmentalModelManager mgr("any.yaml", nullptr, std::move(steps));
    REQUIRE(mgr.isActive());

    WorldState world;
    Venue v;
    v.id = 10;
    v.type_id = 0;
    world.venue_type_names = {"office"};
    world.venues.push_back(v);
    world.people.emplace_back();
    world.people[0].id = 0;
    world.people[0].age = 30;
    world.people[0].sex = Sex::MALE;
    world.people[0].geo_unit_id = -1;
    world.buildIndices();

    Disease disease = makeDiseaseWithCompartmentalUptake(susc_mult);

    int num_infected = 0;
    const int num_trials = 6000;
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
      mgr.advance(0.0f, 0.0f);
      im.processTransmissions(locs, 5.0, 1.0, nullptr, nullptr, nullptr,
                              nullptr, &mgr);
      if (world.people[0].infection) ++num_infected;
    }
    return static_cast<double>(num_infected) / num_trials;
  };

  double prob_full = measure_prob(1.0f);
  double prob_half = measure_prob(0.5f);

  double expected_full = 1.0 - std::exp(-node_output * 1.0 * susc_mult);
  double expected_half = 1.0 - std::exp(-node_output * 0.5 * susc_mult);

  CHECK(std::abs(prob_full - expected_full) < 0.03);
  CHECK(std::abs(prob_half - expected_half) < 0.03);
}

// =============================================================================
// Slice 5: comp_model=nullptr skips compartmental uptake entirely
// =============================================================================

TEST_CASE(
    "Compartmental uptake: comp_model=nullptr → no compartmental infections") {
  WorldState world;
  Venue venue;
  venue.id = 10;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);
  world.people.emplace_back();
  world.people[0].id = 0;
  world.people[0].age = 30;
  world.people[0].sex = Sex::MALE;
  world.people[0].geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeDiseaseWithCompartmentalUptake();
  ContactMatrixConfig cm;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease);

  std::vector<PersonLocation> locs = {{0, 10, -1, -1, 255, 0}};

  // No manager — compartmental uptake path must be skipped
  int new_infections = im.processTransmissions(locs, 5.0, 1.0, nullptr, nullptr,
                                               nullptr, nullptr, nullptr);
  CHECK(new_infections == 0);
  CHECK(world.people[0].infection == nullptr);
}
