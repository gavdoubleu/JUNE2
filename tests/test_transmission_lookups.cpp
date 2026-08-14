#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <memory>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/interaction_manager.h"
#include "loaders/config_loader.h"
#include "test_utils.h"

using namespace june;

TEST_CASE("InteractionManager - Venue Matrix Lookups via Resolved Config") {
  // 1. Load actual config that parses string-based venue names mapping to betas
  ContactMatrixConfig cm =
      ConfigLoader::loadContactMatrices("tests/configs/contact_matrices.yaml");

  // 2. Setup a minimal WorldState with those venue types
  WorldState world;
  Venue venue;
  venue.id = 0;

  // office is venue type 1
  world.venue_type_names = {"household", "office"};
  venue.type_id = 1;
  world.venues.push_back(venue);

  world.people.emplace_back();
  world.people.emplace_back();
  world.people[0].id = 0;
  world.people[0].age = 30;
  world.people[0].sex = Sex::MALE;
  world.people[0].geo_unit_id = -1;
  world.people[1].id = 1;
  world.people[1].age = 25;
  world.people[1].sex = Sex::FEMALE;
  world.people[1].geo_unit_id = -1;
  world.buildIndices();

  // 3. Set up low default contacts — if the matrix lookup is skipped, this
  // default matrix's contacts get used instead of the office matrix's (5.0).
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{0.000001}};
  cm.default_matrix = default_contact_matrix;

  // 4. Resolve every (venue type, mode) into the table transmission reads.
  // household has no matrix of its own here, so accept the default for it.
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world);

  // 5. Set up disease
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  // No fomite modes configured (fomite_configs is empty by default)
  auto cur = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["mild"] = cur;
  trans.symptom_id_curves = {nullptr, cur};

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 10.0}}}});
  trajectories.push_back(td);

  SymptomTag healthy{"healthy", -1, 0};
  SymptomTag mild{"mild", 1, 1};
  Disease disease("TestFlu", {healthy, mild}, {}, trajectories, {}, trans);

  SimulationConfig sim_config;
  ParallelConfig parallel_config;
  InteractionManager im(world, cm, sim_config, parallel_config, &disease,
                        nullptr);

  // 6. Test transmission in "office" venue (venue.id = 0, type = 1)
  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 0, nullptr, "office", 0);

  std::vector<PersonLocation> locs;
  PersonLocation loc0, loc1;
  loc0.person_id = 0;
  loc0.venue_id = 0;
  loc0.subset_index = -1;
  loc0.encounter_type_id = 255;
  loc0.person_array_index = 0;
  loc1.person_id = 1;
  loc1.venue_id = 0;
  loc1.subset_index = -1;
  loc1.encounter_type_id = 255;
  loc1.person_array_index = 1;
  locs.push_back(loc0);
  locs.push_back(loc1);

  // Verify matrix lookup was successful during transmission calculations
  int infections = 0;
  for (int t = 0; t < 1000; ++t) {
    world.people[1].infection.reset();
    im.processTransmissions(locs, 5.0, 24.0, nullptr);
    if (world.people[1].infection) infections++;
  }

  CHECK(infections > 50);
}

TEST_CASE(
    "InteractionManager - Virtual Coordinated Encounter Matrix Resolution "
    "Regression") {
  // 1. Load config with romantic_encounter matrix
  ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
      "tests/configs/romantic_regression.yaml");

  // 2. Set up a world with BOTH a physical venue at venue_type_id 0 AND a
  // virtual encounter at encounter_type_id 0. The integer id collision is
  // the exact precondition that triggered the aliasing bug: passing
  // encounter_type_id through the venue-indexed matrix arrays would have
  // pulled the office matrix (C=9.0) instead of romantic_encounter
  // (C=1.0). The difference is what this test detects.
  WorldState world;
  world.venue_type_names = {"office"};                  // venue_type_id 0
  world.encounter_type_names = {"romantic_encounter"};  // enc_type_id 0

  // Create a "virtual venue" object (negative ID) to mirror how the
  // HDF5Loader stages coordinated encounter venues (type_id = 255).
  Venue v_coordinated;
  v_coordinated.id = -1001;
  v_coordinated.type_id = 255;
  world.venues.push_back(v_coordinated);

  world.people.emplace_back();
  world.people.emplace_back();
  world.people[0].id = 0;
  world.people[0].age = 30;
  world.people[0].sex = Sex::MALE;
  world.people[0].geo_unit_id = -1;
  world.people[1].id = 1;
  world.people[1].age = 25;
  world.people[1].sex = Sex::FEMALE;
  world.people[1].geo_unit_id = -1;
  world.buildIndices();

  // 3. Wire the virtual encounter → matrix mapping the same way
  // CoordinatedEncounterConfig::resolve does in production. The hot path
  // looks up virtual matrices by encounter_type_id, not by name, so this
  // map is the canonical entry point.
  cm.virtual_matrix_names[0] = "romantic_encounter";
  cm.virtual_encounter_type_ids = {0};

  // 4. Resolve config
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world);

  // Direct pointer assertions: the lookup keyed by encounter_type_id must
  // resolve to the romantic_encounter matrix (C=1.0), NOT the office
  // matrix (C=9.0) that sits at the colliding venue_type_id. Before the
  // fix, venue-keyed and encounter-keyed lookups were conflated and this
  // would return the office matrix.
  const ContactMatrix& by_enc = cm.getVirtualBinStructure(0);
  REQUIRE(!by_enc.contacts.empty());
  REQUIRE(!by_enc.contacts[0].empty());
  CHECK(by_enc.contacts[0][0] == doctest::Approx(1.0));

  // And the venue getter should still return the office matrix at the
  // same integer id — these are two independent lookups now, and the
  // test proves it.
  const ContactMatrix& by_venue = cm.getBinStructure(static_cast<uint8_t>(0));
  REQUIRE(!by_venue.contacts.empty());
  REQUIRE(!by_venue.contacts[0].empty());
  CHECK(by_venue.contacts[0][0] == doctest::Approx(9.0));

  // 4. Setup disease with high infectiousness
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto cur = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["mild"] = cur;
  trans.symptom_id_curves = {nullptr, cur};

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 10.0}}}});
  trajectories.push_back(td);

  SymptomTag healthy{"healthy", -1, 0};
  SymptomTag mild{"mild", 1, 1};
  Disease disease("TestFlu", {healthy, mild}, {}, trajectories, {}, trans);

  SimulationConfig sim_config;
  ParallelConfig parallel_config;
  InteractionManager im(world, cm, sim_config, parallel_config, &disease,
                        nullptr);

  // 5. Test transmission with romantic_encounter (encounter_type_id = 0)
  // Use the virtual venue ID
  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 0, nullptr, "romantic_encounter", 0);

  std::vector<PersonLocation> locs;
  PersonLocation loc0, loc1;
  loc0.person_id = 0;
  loc0.venue_id = -1001;
  loc0.subset_index = -1;
  loc0.encounter_type_id = 0;
  loc0.person_array_index = 0;
  loc1.person_id = 1;
  loc1.venue_id = -1001;
  loc1.subset_index = -1;
  loc1.encounter_type_id = 0;
  loc1.person_array_index = 1;
  locs.push_back(loc0);
  locs.push_back(loc1);

  // Verify that the romantic_encounter contact matrix (contacts=1.0) is used
  // for virtual encounters, not the default matrix.
  // FoI formula: lambda = delta_hours * C / bin_size * infectiousness *
  // susc_mult With C=1.0 (from matrix), delta=24h, bin_size=1: lambda = 24 →
  // prob ≈ 1.0
  int infections = 0;
  for (int t = 0; t < 100; ++t) {
    world.people[1].infection.reset();
    im.processTransmissions(locs, 5.0, 24.0, nullptr);
    if (world.people[1].infection) infections++;
  }

  CHECK(infections > 30);
}

TEST_CASE(
    "ContactMatrixConfig - encounter type with no matrix of its own is "
    "refused at load unless the default is opted into") {
  // An encounter type declared in the world but never given a
  // virtual_contact_matrix entry (missing config, typo'd name, ...) has no
  // contact structure of its own. Borrowing the scenario-wide default for it
  // has to be spelled out in the scenario, because the resulting physics is
  // not what the scenario says about that encounter type.
  auto build = [](WorldState& world, ContactMatrixConfig& cm) {
    world.venue_type_names = {"office"};
    world.encounter_type_names = {"unconfigured_encounter"};
    world.buildIndices();
    // Held at a virtual venue, so it needs a matrix of its own; it just has
    // not been given one. CoordinatedEncounterConfig::resolve records this
    // in production.
    cm.virtual_encounter_type_ids = {0};
  };

  SUBCASE("strict by default: load refuses and names the encounter type") {
    ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
        "tests/configs/romantic_regression.yaml");
    WorldState world;
    build(world, cm);
    REQUIRE(cm.allow_default_matrix == false);
    CHECK_THROWS_WITH_AS(finalizeContactMatrices(cm, world),
                         doctest::Contains("unconfigured_encounter"),
                         std::runtime_error);
  }

  SUBCASE("opted in: resolves to the default matrix") {
    ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
        "tests/configs/romantic_regression.yaml");
    WorldState world;
    build(world, cm);
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world);

    const ContactMatrix& resolved = cm.getVirtualMatrix(0, 0);
    REQUIRE(!resolved.contacts.empty());
    REQUIRE(!resolved.contacts[0].empty());
    CHECK(resolved.contacts[0][0] == doctest::Approx(1.0));
  }
}
