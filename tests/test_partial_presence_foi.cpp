#define DOCTEST_CONFIG_IMPLEMENT
#include <cmath>
#include <memory>
#include <vector>

#include "activity/runtime_group_allocator.h"
#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/interaction_manager.h"
#include "test_utils.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace june;

// The RuntimeGroupAllocator calls MPI_Comm_size during allocateForSlot when
// the build has USE_MPI defined, even with a single process. Initialise MPI
// here so a serial doctest run doesn't abort.
int main(int argc, char** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  int res = context.run();
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return res;
}

// =============================================================================
// Partial-presence FOI integration tests.
//
// Each test asserts a physics property of the route-dynamics pipeline:
//   1. Partial-overlap exposure during shared window only:
//      a susceptible's λ comes ONLY from sub-intervals where an infectious
//      rider is co-present in the same runtime bin.
//   2. Bin isolation:
//      a susceptible bucketed alone (target_group_size = 1, num_bins >= 2)
//      accrues zero λ even when other riders on the same venue are infectious.
//   3. Gate non-regression:
//      declaring partial_presence on a venue type that no venue in the
//      scenario uses leaves the fast-path FOI behavior identical to the
//      empty-partial-presence baseline.
//   4. Default-sourced contact structure:
//      a partial-presence venue type with no matrix of its own resolves
//      against the scenario default and accrues that default's contact rate,
//      instead of failing the moment such a venue transmits.
//
// Synthetic worlds with anonymous IDs — see
// feedback_tests_name_physics_not_scenarios memory: no Durham, no
// Chester-le-Street naming in test code.
// =============================================================================

namespace {

// Build a minimal disease with constant unit infectiousness so integrated
// infectiousness over delta_hours = delta_hours. Caller owns the unique_ptr;
// returning by value would force a copy that Disease doesn't expose.
// (Mirrors the inline pattern in test_foi_calculation.cpp.)
std::unique_ptr<Disease> makeUnitConstantDisease() {
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
  return std::make_unique<Disease>(
      "Flu", std::vector<SymptomTag>{healthy, sick}, stage_settings,
      trajectories, OutcomeRates{}, trans);
}

// Build a single-bin, single-mode contact matrix for a partial-presence
// venue type and register it in the config, then run the same finalize chain
// Simulator runs so the integer-keyed lookups transmission reads are filled.
void registerSimpleContactMatrix(ContactMatrixConfig& cm,
                                 const std::string& venue_type_name,
                                 const WorldState& world,
                                 double contacts_per_slot) {
  ContactMatrix mat;
  mat.bins = {"rider"};
  mat.contacts = {{contacts_per_slot}};
  mat.proportion_physical = {{0.0}};
  cm.matrices[venue_type_name] = mat;
  cm.mode_matrices[venue_type_name]["respiratory"] = mat;
  if (cm.mode_names.empty()) cm.mode_names = {"respiratory"};
  finalizeContactMatrices(cm, world);
}

// Add a person to the world. Sets only the fields the partial-presence
// pipeline reads; tests overlay infection / immunity as needed.
void addAdult(WorldState& world, PersonId pid) {
  Person& p = world.people.emplace_back();
  p.id = pid;
  p.age = 30.0f + static_cast<float>(pid);
  p.sex = Sex::MALE;
  p.geo_unit_id = 0;
}

// Common geo unit + activity name registration. Called once per test world.
void registerScaffolding(WorldState& world) {
  world.geo_level_names = {"city"};
  GeographicalUnit gu;
  gu.id = 0;
  gu.name = "test_unit";
  gu.level_id = 0;
  gu.parent_id = -1;
  world.geo_units.push_back(gu);
  world.activity_names = {"commute"};
}

// Compute the integrated infectiousness from a single infectious rider over
// a fraction `scale` of the slot. With ConstantCurve(1.0) the full-slot
// integrated infectiousness equals delta_hours (see test_foi_calculation
// comments). Scaling by sub_interval / slot_duration yields the per-sub
// contribution.
double expectedInfContribution(double delta_hours, double sub_dur_min,
                               double slot_dur_min) {
  return delta_hours * (sub_dur_min / slot_dur_min);
}

}  // namespace

// -----------------------------------------------------------------------------
// Test 1: Partial-overlap exposure during shared window only.
//
// Synthetic world: one partial-presence venue ("route_line"), two riders.
//   Rider A (id 0) — infectious, present [0, 20)
//   Rider B (id 1) — susceptible, present [10, 20)
// Slot duration = 60 min. target_group_size large enough that both riders
// land in a single bin. Single matrix bin.
//
// Expected λ_B:
//   - sub-interval [0, 10): B not present → no contribution.
//   - sub-interval [10, 20): both present, sub_dur = 10.
//       scale = 10/60 = 1/6
//       inf_full_A = delta_hours (= 1.0 for slot_min=60, delta_h=1)
//       inf_sub_A  = 1/6
//       N_present (own bin) = max(1, 2 - 1) = 1
//       omega = contacts / N = contacts / 1
//       contrib = omega * inf_sub_A = contacts * (1/6)
//   - sub-interval [20, 60): neither present → no contribution.
//
// With contacts = 6.0 we expect λ_B = 1.0 exactly.
// -----------------------------------------------------------------------------
TEST_CASE(
    "partial-presence FOI: susceptible accrues lambda only during shared "
    "window") {
  WorldState world;
  registerScaffolding(world);

  Config config;
  SimulationConfig& sim_cfg = config.simulation;
  sim_cfg.random_seed = 12345;

  const uint8_t line_type_id = addPartialPresenceVenueType(
      world, sim_cfg, "route_line", /*target_group_size=*/100);
  const VenueId line = addPartialPresenceVenue(world, line_type_id);

  addAdult(world, /*pid=*/0);  // infectious
  addAdult(world, /*pid=*/1);  // susceptible

  std::vector<PartialPresenceLegSpec> legs;
  legs.push_back(PartialPresenceLegSpec{0, line, 0.0f, 20.0f});
  legs.push_back(PartialPresenceLegSpec{1, line, 10.0f, 20.0f});
  addPartialPresenceLegs(world, /*activity_index=*/0, legs);
  world.buildIndices();

  // Contact matrix: single rider bin, contacts = 6.0 → tidy expected λ = 1.0.
  ContactMatrixConfig cm;
  registerSimpleContactMatrix(cm, "route_line", world, /*contacts=*/6.0);

  auto disease = makeUnitConstantDisease();
  ParallelConfig parallel_config;

  // Make rider A infectious. Constant trajectory means they stay infectious
  // throughout the slot.
  world.people[0].infection =
      std::make_unique<Infection>(disease.get(), 0.0, &world.people[0],
                                  /*seed=*/7, &world, "route_line", line);

  // Allocator + interaction manager wiring.
  RuntimeGroupAllocator allocator(world, config);
  InteractionManager im(world, cm, sim_cfg, parallel_config, disease.get(),
                        nullptr);
  im.setRuntimeGroupAllocator(&allocator);

  // Drive allocator for this slot.
  std::vector<PersonLocation> locs =
      makePartialPresenceLocations(world, /*activity_index=*/0);
  TimeSlot slot;  // contents unused by allocator
  const double delta_hours = 1.0;
  allocator.allocateForSlot(/*time_slot_index=*/0, /*day_type_idx=*/0, slot,
                            /*current_simulation_time=*/0.0, delta_hours, locs);

  // Sanity: both riders bucketed, single bin.
  REQUIRE(allocator.getNumGroups(line) == 1);
  REQUIRE(allocator.getGroupIndex(line, 0) == 0);
  REQUIRE(allocator.getGroupIndex(line, 1) == 0);

  // Build the InteractionMember list (one per location, both on `line`).
  std::vector<InteractionMember> members;
  for (const auto& l : locs) {
    InteractionMember m;
    m.id = l.person_id;
    m.array_index = l.person_array_index;
    m.subset_index = l.subset_index;
    m.encounter_type_id = l.encounter_type_id;
    members.push_back(m);
  }

  Venue* venue_ptr = world.getVenue(line);
  REQUIRE(venue_ptr != nullptr);

  // Drive the λ accumulator directly (no Bernoulli, no infection writes).
  auto acc = im.computePartialPresenceLambda(
      members, venue_ptr, line, /*current_time=*/5.0, delta_hours,
      /*visitor_data=*/nullptr, /*encounter_type_id=*/255);

  // Expected: λ_B = contacts * (overlap / slot) * (full_inf / N_bin)
  //               = 6.0 * (10 / 60) * (1.0 / 1) = 1.0
  // λ_A = 0 (A is infectious, not susceptible).
  const double expected_lambda_b = 6.0 * (10.0 / 60.0) * 1.0;

  REQUIRE(acc.susc_lambda.count(/*B=*/1) == 1);
  CHECK(acc.susc_lambda.at(1) ==
        doctest::Approx(expected_lambda_b).epsilon(1e-6));
  // A must not appear as a susceptible — they are the infectious rider.
  CHECK(acc.susc_lambda.count(/*A=*/0) == 0);

  // Source attribution: B's λ must be attributable entirely to A.
  REQUIRE(acc.susc_sources.count(1) == 1);
  const auto& sources_b = acc.susc_sources.at(1);
  REQUIRE_FALSE(sources_b.empty());
  double total_weight = 0.0;
  for (const auto& s : sources_b) {
    CHECK(s.infector == /*A=*/0);
    total_weight += s.weighted;
  }
  CHECK(total_weight == doctest::Approx(expected_lambda_b).epsilon(1e-6));
}

// -----------------------------------------------------------------------------
// Test 2: Bin isolation.
//
// Same two-rider setup, but target_group_size = 1 so the allocator buckets
// each rider into their OWN runtime bin (num_bins = 2). The susceptible
// rider must accrue zero λ because their bin contains only themselves —
// the infectious rider is in a different bin and cannot contribute.
//
// Control: with target_group_size = 100 (the test above) both riders share
// a bin and λ_B > 0. This test inverts the bucketing and asserts λ_B == 0.
// -----------------------------------------------------------------------------
TEST_CASE("partial-presence FOI: bin isolation yields zero cross-bin lambda") {
  WorldState world;
  registerScaffolding(world);

  Config config;
  SimulationConfig& sim_cfg = config.simulation;
  sim_cfg.random_seed = 12345;

  const uint8_t line_type_id = addPartialPresenceVenueType(
      world, sim_cfg, "route_line", /*target_group_size=*/1);
  const VenueId line = addPartialPresenceVenue(world, line_type_id);

  addAdult(world, /*pid=*/0);  // infectious
  addAdult(world, /*pid=*/1);  // susceptible

  // Full-slot presence — eliminates partial-overlap as a confound. Bin
  // isolation is the only thing under test.
  std::vector<PartialPresenceLegSpec> legs;
  legs.push_back(PartialPresenceLegSpec{0, line, 0.0f, 60.0f});
  legs.push_back(PartialPresenceLegSpec{1, line, 0.0f, 60.0f});
  addPartialPresenceLegs(world, /*activity_index=*/0, legs);
  world.buildIndices();

  ContactMatrixConfig cm;
  registerSimpleContactMatrix(cm, "route_line", world, /*contacts=*/6.0);

  auto disease = makeUnitConstantDisease();
  ParallelConfig parallel_config;
  world.people[0].infection =
      std::make_unique<Infection>(disease.get(), 0.0, &world.people[0],
                                  /*seed=*/7, &world, "route_line", line);

  RuntimeGroupAllocator allocator(world, config);
  InteractionManager im(world, cm, sim_cfg, parallel_config, disease.get(),
                        nullptr);
  im.setRuntimeGroupAllocator(&allocator);

  std::vector<PersonLocation> locs =
      makePartialPresenceLocations(world, /*activity_index=*/0);
  TimeSlot slot;
  const double delta_hours = 1.0;
  allocator.allocateForSlot(0, 0, slot, 0.0, delta_hours, locs);

  // tgs=1, 2 riders → 2 bins, one rider in each. Different bins is the
  // load-bearing invariant for this test; assert it explicitly.
  REQUIRE(allocator.getNumGroups(line) == 2);
  const uint16_t bin_a = allocator.getGroupIndex(line, 0);
  const uint16_t bin_b = allocator.getGroupIndex(line, 1);
  REQUIRE(bin_a != bin_b);

  std::vector<InteractionMember> members;
  for (const auto& l : locs) {
    InteractionMember m;
    m.id = l.person_id;
    m.array_index = l.person_array_index;
    m.subset_index = l.subset_index;
    m.encounter_type_id = l.encounter_type_id;
    members.push_back(m);
  }

  Venue* venue_ptr = world.getVenue(line);
  auto acc = im.computePartialPresenceLambda(
      members, venue_ptr, line, /*current_time=*/5.0, delta_hours,
      /*visitor_data=*/nullptr, /*encounter_type_id=*/255);

  // B's bin contains only B — no infectious co-rider → no λ accumulated.
  // The FOI loop never inserts B into susc_lambda when there is no
  // contribution, so we expect the map either to lack B or to map B to 0.
  if (acc.susc_lambda.count(1) > 0) {
    CHECK(acc.susc_lambda.at(1) == doctest::Approx(0.0).epsilon(1e-9));
  }
  // No source attribution either.
  if (acc.susc_sources.count(1) > 0) {
    CHECK(acc.susc_sources.at(1).empty());
  }
}

// -----------------------------------------------------------------------------
// Test 3: Gate non-regression.
//
// The gate in processVenueTransmissions consults
// SimulationConfig::partial_presence::enabled_venue_type_mask. Declaring
// partial-presence on a venue type that no venue in the scenario uses must
// NOT change FOI behavior on the regular fast path.
//
// Setup: a tiny "office" scenario with 1 infectious + 1 susceptible.
//   Run A — partial_presence config is empty.
//   Run B — partial_presence declares "route_line" with target_group_size 100
//           (but the world has no route_line venues).
// With identical seeds, the resulting infection state must be identical
// in both runs.
// -----------------------------------------------------------------------------
TEST_CASE(
    "partial-presence FOI: declaring partial-presence on an unrelated type "
    "leaves fast-path behavior unchanged") {
  auto run_one = [&](bool declare_unrelated_partial_presence, int seed) {
    WorldState world;
    world.geo_level_names = {"city"};
    GeographicalUnit gu;
    gu.id = 0;
    gu.name = "test_unit";
    gu.level_id = 0;
    gu.parent_id = -1;
    world.geo_units.push_back(gu);

    world.venue_type_names = {"office"};
    Venue v;
    v.id = 0;
    v.type_id = 0;  // office — NOT partial-presence
    v.geo_unit_id = 0;
    v.parent_id = -1;
    world.venues.push_back(v);

    for (int i = 0; i < 2; ++i) {
      Person& p = world.people.emplace_back();
      p.id = i;
      p.age = 30.0f;
      p.sex = Sex::MALE;
      p.geo_unit_id = 0;
    }
    world.buildIndices();

    ContactMatrixConfig cm;
    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{0.3}};
    cm.default_matrix = default_contact_matrix;

    Config config;
    SimulationConfig& sim_cfg = config.simulation;
    sim_cfg.random_seed = static_cast<uint64_t>(seed);

    if (declare_unrelated_partial_presence) {
      // Declares a partial-presence venue type that doesn't appear in this
      // world. The gate must skip the office venue (type_id 0) because the
      // mask bit for type_id 0 stays clear.
      sim_cfg.partial_presence.target_group_size_by_name["route_line"] = 100;
      // Reserve a phantom type_id 7 — past office's id 0, no venue uses it.
      sim_cfg.partial_presence.target_group_size_by_type_id.resize(8, 0);
      sim_cfg.partial_presence.target_group_size_by_type_id[7] = 100;
      sim_cfg.partial_presence.enabled_venue_type_mask = (1ULL << 7);
    }

    ParallelConfig parallel_config;
    auto disease = makeUnitConstantDisease();

    world.people[0].infection =
        std::make_unique<Infection>(disease.get(), 0.0, &world.people[0],
                                    /*seed=*/seed, &world, "office", 0);

    // Wire an allocator only when partial_presence is declared (mirrors how
    // Simulator wires it conditionally).
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, *disease);
    RuntimeGroupAllocator allocator(world, config);
    InteractionManager im(world, cm, sim_cfg, parallel_config, disease.get(),
                          nullptr);
    if (declare_unrelated_partial_presence) {
      im.setRuntimeGroupAllocator(&allocator);
      TimeSlot slot;
      allocator.allocateForSlot(0, 0, slot, /*current_time=*/0.0,
                                /*delta_hours=*/4.0, /*locations=*/{});
    }

    std::vector<PersonLocation> locs;
    locs.push_back({0, 0, -1, 0, 255, 0});
    locs.push_back({1, 0, -1, 0, 255, 1});
    im.processTransmissions(locs, /*current_time=*/5.0, /*delta_hours=*/4.0,
                            nullptr);
    return world.people[1].infection != nullptr;
  };

  // Run identical scenarios under both partial_presence configurations across
  // a small seed sweep. Outcomes must match exactly per seed — the gate is
  // pure on (venue_type_id, mask) and the mask never selects the office venue.
  for (int seed = 1; seed <= 25; ++seed) {
    const bool baseline = run_one(/*declare=*/false, seed);
    const bool gated = run_one(/*declare=*/true, seed);
    CHECK_MESSAGE(baseline == gated, "seed=", seed);
  }
}

// -----------------------------------------------------------------------------
// Test 4: A partial-presence venue type whose contact structure comes from the
// scenario default, not from an entry of its own.
//
// The full-presence path has always tolerated this: a venue type with no
// matrix of its own picks up the scenario default. The partial-presence path
// used a lookup that could not see the default at all, so the same scenario
// killed the run the first time such a venue transmitted.
//
// The physics asserted here is the same partial-overlap property as Test 1 —
// λ accrues only over the shared window — with the contact rate sourced from
// the default matrix. Two riders, windows [0,20) and [10,20), slot 60 min:
//   λ_B = contacts * (overlap / slot) * (full_inf / N_bin)
//       = 6.0 * (10 / 60) * (1.0 / 1) = 1.0
// -----------------------------------------------------------------------------
TEST_CASE(
    "partial-presence FOI: venue type with only a default contact matrix "
    "resolves and accrues the default's contact rate") {
  auto buildWorld = [](WorldState& world, Config& config) -> VenueId {
    registerScaffolding(world);
    const uint8_t line_type_id = addPartialPresenceVenueType(
        world, config.simulation, "route_line", /*target_group_size=*/100);
    const VenueId line = addPartialPresenceVenue(world, line_type_id);
    addAdult(world, /*pid=*/0);  // infectious
    addAdult(world, /*pid=*/1);  // susceptible
    std::vector<PartialPresenceLegSpec> legs;
    legs.push_back(PartialPresenceLegSpec{0, line, 0.0f, 20.0f});
    legs.push_back(PartialPresenceLegSpec{1, line, 10.0f, 20.0f});
    addPartialPresenceLegs(world, /*activity_index=*/0, legs);
    world.buildIndices();
    return line;
  };

  // The only contact structure in the scenario: no per-venue-type entry.
  auto defaultOnlyMatrices = []() {
    ContactMatrixConfig cm;
    ContactMatrix mat;
    mat.bins = {"rider"};
    mat.contacts = {{6.0}};
    mat.proportion_physical = {{0.0}};
    cm.default_matrix = mat;
    return cm;
  };

  SUBCASE("leaning on the default is refused at load unless stated") {
    WorldState world;
    Config config;
    buildWorld(world, config);
    ContactMatrixConfig cm = defaultOnlyMatrices();
    REQUIRE(cm.allow_default_matrix == false);
    CHECK_THROWS_WITH_AS(finalizeContactMatrices(cm, world),
                         doctest::Contains("route_line"), std::runtime_error);
  }

  SUBCASE("stated: lambda accrues over the shared window at the default rate") {
    WorldState world;
    Config config;
    SimulationConfig& sim_cfg = config.simulation;
    sim_cfg.random_seed = 12345;
    const VenueId line = buildWorld(world, config);

    ContactMatrixConfig cm = defaultOnlyMatrices();
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world);

    auto disease = makeUnitConstantDisease();
    ParallelConfig parallel_config;
    world.people[0].infection =
        std::make_unique<Infection>(disease.get(), 0.0, &world.people[0],
                                    /*seed=*/7, &world, "route_line", line);

    RuntimeGroupAllocator allocator(world, config);
    InteractionManager im(world, cm, sim_cfg, parallel_config, disease.get(),
                          nullptr);
    im.setRuntimeGroupAllocator(&allocator);

    std::vector<PersonLocation> locs =
        makePartialPresenceLocations(world, /*activity_index=*/0);
    TimeSlot slot;
    const double delta_hours = 1.0;
    allocator.allocateForSlot(/*time_slot_index=*/0, /*day_type_idx=*/0, slot,
                              /*current_simulation_time=*/0.0, delta_hours,
                              locs);
    REQUIRE(allocator.getNumGroups(line) == 1);

    std::vector<InteractionMember> members;
    for (const auto& l : locs) {
      InteractionMember m;
      m.id = l.person_id;
      m.array_index = l.person_array_index;
      m.subset_index = l.subset_index;
      m.encounter_type_id = l.encounter_type_id;
      members.push_back(m);
    }

    Venue* venue_ptr = world.getVenue(line);
    REQUIRE(venue_ptr != nullptr);

    auto acc = im.computePartialPresenceLambda(
        members, venue_ptr, line, /*current_time=*/5.0, delta_hours,
        /*visitor_data=*/nullptr, /*encounter_type_id=*/255);

    const double expected_lambda_b = 6.0 * (10.0 / 60.0) * 1.0;
    REQUIRE(acc.susc_lambda.count(/*B=*/1) == 1);
    CHECK(acc.susc_lambda.at(1) ==
          doctest::Approx(expected_lambda_b).epsilon(1e-6));
    CHECK(acc.susc_lambda.count(/*A=*/0) == 0);
  }
}
