#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <memory>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/infectiousness_curves.h"
#include "mock_compartmental_model.h"
#include "simulation/compartmental_model_manager.h"
#include "transmission_fixtures.h"

using namespace june;

// computeDepositionWriteback is declared in compartmental_model_manager.h

// ---------------------------------------------------------------------------
// Disease with one compartmental_deposition mode.
// ---------------------------------------------------------------------------
static Disease makeDiseaseWithDeposition(double deposition_rate = 0.03) {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  std::vector<SymptomTag> tags = {{"healthy", -1, 0}, {"infectious", 1, 1}};
  tp.symptom_id_curves = {nullptr, nullptr};

  TransmissionMode dep_mode;
  dep_mode.name = "comp_dep";
  dep_mode.type = TransmissionModeType::CompartmentalDeposition;
  dep_mode.symptom_curves = {nullptr, nullptr};
  CompartmentalDepositionConfig dcfg;
  dcfg.mode_index = 0;
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

  return Disease("TestDep", tags, {}, trajectories, {}, tp);
}

// Manager that:
//   - uses a caller-supplied InProcessMockPlugin (for ownership/inspection)
//   - maps every venue in node_venues to local node index 0
static std::unique_ptr<CompartmentalModelManager> makeManager(
    InProcessMockPlugin& plugin, const std::vector<int>& node_venues,
    float coupling_beta = 0.001f) {
  CompartmentalModelSteps steps;
  steps.loadSidecar =
      [coupling_beta](const std::string&) -> PluginSidecarConfig {
    PluginSidecarConfig cfg;
    cfg.plugin_so_path = "ignored.so";
    cfg.default_human_to_compartmental_model_input = coupling_beta;
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
      [node_venues](const std::string&, DomainManager*,
                    std::unordered_map<int, int>& venue_map) -> NodePartition {
    NodePartition p;
    p.owned_node_indices = {0};
    p.node_venue_ids = {node_venues};
    for (int vid : node_venues) venue_map[vid] = 0;
    return p;
  };
  return std::make_unique<CompartmentalModelManager>("any.yaml", nullptr,
                                                     std::move(steps));
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("Write-back: inactive manager → no plugin call") {
  WorldState world;
  world.people.emplace_back();
  world.people[0].id = 0;
  world.buildIndices();
  Disease disease = makeDiseaseWithDeposition();
  CompartmentalModelManager inactive_mgr("", nullptr);
  REQUIRE_FALSE(inactive_mgr.isActive());
  std::vector<PersonLocation> locs;
  REQUIRE_NOTHROW(inactive_mgr.computeDepositionWriteback(
      locs, world, disease, 5.0, 5.0 + 1.0 / 24.0));
}

TEST_CASE(
    "Write-back: disease has no deposition modes → writeCouplingInputs not "
    "called") {
  WorldState world;
  Venue v;
  v.id = 10;
  v.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(v);
  world.people.emplace_back();
  world.people[0].id = 0;
  world.buildIndices();

  // Direct-only disease (no compartmental_deposition)
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(1.0);
  tp.symptom_id_curves = {nullptr, curve};
  TransmissionMode direct;
  direct.name = "direct";
  direct.symptom_curves = {nullptr, curve};
  tp.modes.push_back(std::move(direct));
  std::vector<SymptomTag> tags = {{"healthy", -1, 0}, {"infectious", 1, 1}};
  Disease no_dep("NoDep", tags, {}, {}, {}, tp);

  InProcessMockPlugin plugin;
  auto mgr = makeManager(plugin, {10});
  REQUIRE(mgr->isActive());

  std::vector<PersonLocation> locs = {{0, 10, -1, -1, 255, 0}};
  mgr->computeDepositionWriteback(locs, world, no_dep, 5.0, 5.0 + 1.0 / 24.0);
  // No deposition modes → writeCouplingInputs receives all zeros or is not
  // called Either way, plugin should show 0 non-zero entries
  for (float v : plugin.last_coupling_inputs) CHECK(v == doctest::Approx(0.0f));
}

TEST_CASE("Write-back: uninfected person → zero deposition per node") {
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

  Disease disease = makeDiseaseWithDeposition(0.03);
  InProcessMockPlugin plugin;
  auto mgr = makeManager(plugin, {10});
  REQUIRE(mgr->isActive());

  // person NOT infected
  std::vector<PersonLocation> locs = {{0, 10, -1, -1, 255, 0}};
  mgr->computeDepositionWriteback(locs, world, disease, 5.0, 5.0 + 1.0 / 24.0);

  REQUIRE(plugin.coupling_inputs_call_count == 1);
  CHECK(plugin.last_coupling_inputs.size() == 1);
  CHECK(plugin.last_coupling_inputs[0] == doctest::Approx(0.0f));
}

TEST_CASE("Write-back: infected person at mapped venue → positive deposition") {
  const double deposition_rate = 0.03;
  const double delta_days = 1.0 / 24.0;  // 1 hour
  // ConstantCurve.integrate(a, b) = value * (b - a)
  // integrateFomiteDeposition-style: 24 * integral_in_days = 24 * dep *
  // delta_days
  const double expected_dep = 24.0 * deposition_rate * delta_days;

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

  Disease disease = makeDiseaseWithDeposition(deposition_rate);
  InProcessMockPlugin plugin;
  auto mgr = makeManager(plugin, {10});
  REQUIRE(mgr->isActive());

  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 42, kNoTransmissionContext);

  std::vector<PersonLocation> locs = {{0, 10, -1, -1, 255, 0}};
  double t0 = 5.0, t1 = t0 + delta_days;
  mgr->computeDepositionWriteback(locs, world, disease, t0, t1);

  REQUIRE(plugin.coupling_inputs_call_count == 1);
  REQUIRE(plugin.last_coupling_inputs.size() == 1);
  // Deposition scaled by coupling_beta (default matrix value = 0.001)
  double expected_scaled = expected_dep * mgr->getCouplingMatrix().getDefault();
  CHECK(plugin.last_coupling_inputs[0] ==
        doctest::Approx(static_cast<float>(expected_scaled)).epsilon(0.01));
}

TEST_CASE(
    "Write-back: person at unmapped venue → excluded from node deposition") {
  WorldState world;
  Venue v1;
  v1.id = 10;
  v1.type_id = 0;
  Venue v2;
  v2.id = 20;
  v2.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(v1);
  world.venues.push_back(v2);
  world.people.emplace_back();
  world.people[0].id = 0;
  world.people[0].age = 30;
  world.people[0].sex = Sex::MALE;
  world.people[0].geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeDiseaseWithDeposition(0.03);
  InProcessMockPlugin plugin;
  // node 0 owns venue 10; venue 20 is not mapped
  auto mgr = makeManager(plugin, {10});
  REQUIRE(mgr->venueToLocalNodeIndex(10) == 0);
  REQUIRE(mgr->venueToLocalNodeIndex(20) == -1);

  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 42, kNoTransmissionContext);

  // person at venue 20 (unmapped)
  std::vector<PersonLocation> locs = {{0, 20, -1, -1, 255, 0}};
  mgr->computeDepositionWriteback(locs, world, disease, 5.0, 5.0 + 1.0 / 24.0);

  REQUIRE(plugin.coupling_inputs_call_count == 1);
  CHECK(plugin.last_coupling_inputs[0] == doctest::Approx(0.0f));
}
