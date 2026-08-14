#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <memory>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/epidemiology.h"
#include "epidemiology/interaction_manager.h"
#include "test_utils.h"

using namespace june;

// =============================================================================
// Helper: Build a Disease with 2 modes — direct (mode 0) and fomite (mode 1)
// =============================================================================
static Disease makeDiseaseWithFomite(double deposition_rate = 1.0,
                                     double fomite_infectiousness = 1.0,
                                     double fomite_susc_mult = 1.0) {
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;

  // Symptom tags: 0 = healthy, 1 = mild
  std::vector<SymptomTag> symptom_tags = {
      {"healthy", -1, 0},
      {"mild", 1, 1},
  };

  // Direct mode curves: mild has constant infectiousness = 1.0
  auto constant_curve = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["mild"] = constant_curve;
  trans.symptom_id_curves = {nullptr, constant_curve};

  // Mode 0: "direct" — Standard mode
  TransmissionMode direct_mode;
  direct_mode.name = "direct";
  direct_mode.symptom_curves = {nullptr, constant_curve};

  // Mode 1: "fomite_env" — Fomite mode
  TransmissionMode fomite_mode;
  fomite_mode.name = "fomite_env";
  fomite_mode.type = TransmissionModeType::Fomite;
  fomite_mode.mode_transmissibility_multiplier = fomite_susc_mult;
  fomite_mode.symptom_curves = {nullptr, nullptr};

  FomiteConfig fcfg;
  fcfg.mode_index = 1;
  fcfg.max_age = 2.0;
  fcfg.infectiousness_curve =
      std::make_shared<ConstantCurve>(fomite_infectiousness);
  fcfg.deposition_by_symptom = {
      nullptr,
      std::make_shared<ConstantCurve>(deposition_rate),
  };
  fomite_mode.config = std::move(fcfg);

  trans.modes.push_back(std::move(direct_mode));
  trans.modes.push_back(std::move(fomite_mode));

  // Natural immunity (default)
  trans.natural_immunity.level = 0.95;
  trans.natural_immunity.waning_rate = 0.001;

  // Trajectories: everyone gets "mild" for 100 hours
  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  trajectories.push_back(td);

  return Disease("TestFomite", symptom_tags, {}, trajectories, {}, trans);
}

// Helper: Initialize venue fomite_history for all venues
static void initVenueFomiteHistory(WorldState& world, int num_fomite_modes) {
  for (auto& venue : world.venues) {
    venue.fomite_history.assign(num_fomite_modes, {});
  }
}

// =============================================================================
// Tests
// =============================================================================

TEST_CASE("Fomite Deposition") {
  WorldState world;
  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);

  world.people.emplace_back();
  Person& p1 = world.people[0];
  p1.id = 0;
  p1.age = 30;
  p1.sex = Sex::MALE;
  p1.geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeDiseaseWithFomite(2.0 /* deposition_rate */);
  initVenueFomiteHistory(world, 1);  // 1 fomite mode

  // Infect person 0
  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 42, nullptr, "office", 0);

  // Place person at venue
  std::vector<PersonLocation> locs;
  PersonLocation loc;
  loc.person_id = 0;
  loc.venue_id = 0;
  loc.subset_index = -1;
  loc.activity_index = -1;
  loc.encounter_type_id = 255;
  loc.person_array_index = 0;
  locs.push_back(loc);

  ContactMatrixConfig cm_config;
  SimulationConfig sim_config;
  ParallelConfig parallel_config;
  cm_config.allow_default_matrix = true;
  finalizeContactMatrices(cm_config, world, disease);
  InteractionManager im(world, cm_config, sim_config, parallel_config, &disease,
                        nullptr);

  double delta_hours = 1.0;
  im.processTransmissions(locs, 5.0, delta_hours, nullptr);

  // Fomite history for mode 0 (local index) should have entries
  CHECK(!world.venues[0].fomite_history[0].empty());
  if (!world.venues[0].fomite_history[0].empty()) {
    double deposited = world.venues[0].fomite_history[0].back().amount;
    // deposition_rate=2.0, delta_hours=1.0 → amount ≈ 2.0
    CHECK(deposited > 0.0);
    CHECK(deposited == doctest::Approx(2.0).epsilon(0.5));
  }
}

TEST_CASE("No Deposition When Healthy") {
  WorldState world;
  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);

  world.people.emplace_back();
  Person& p1 = world.people[0];
  p1.id = 0;
  p1.age = 30;
  p1.sex = Sex::MALE;
  p1.geo_unit_id = -1;
  // p1 is NOT infected
  world.buildIndices();

  Disease disease = makeDiseaseWithFomite();
  initVenueFomiteHistory(world, 1);

  std::vector<PersonLocation> locs;
  PersonLocation loc;
  loc.person_id = 0;
  loc.venue_id = 0;
  loc.subset_index = -1;
  loc.activity_index = -1;
  loc.encounter_type_id = 255;
  loc.person_array_index = 0;
  locs.push_back(loc);

  ContactMatrixConfig cm_config;
  SimulationConfig sim_config;
  ParallelConfig parallel_config;
  cm_config.allow_default_matrix = true;
  finalizeContactMatrices(cm_config, world, disease);
  InteractionManager im(world, cm_config, sim_config, parallel_config, &disease,
                        nullptr);

  im.processTransmissions(locs, 5.0, 1.0, nullptr);

  // No infected person → no deposition
  CHECK(world.venues[0].fomite_history[0].empty());
}

TEST_CASE("Fomite Max Age Pruning") {
  WorldState world;
  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);
  world.buildIndices();

  // Disease with max_age = 2 days
  Disease disease = makeDiseaseWithFomite();

  // Manually set up fomite history
  world.venues[0].fomite_history.assign(1, {});
  // Old event: deposited at time 0.0 days
  world.venues[0].fomite_history[0].push_back({0.0, 5.0});
  // Recent event: deposited at time 4.5 days
  world.venues[0].fomite_history[0].push_back({4.5, 3.0});

  CHECK(world.venues[0].fomite_history[0].size() == 2);

  // Create Epidemiology to call updateVenueFomites
  Epidemiology epi(world, &disease);

  // current_time = 5.0 days → age of old event = 5.0 days > max_age(2 days)
  // age of recent event = 0.5 days < max_age(2 days)
  epi.updateVenueFomites(5.0, 1.0);

  // Old event should be pruned, recent should remain
  CHECK(world.venues[0].fomite_history[0].size() == 1);
  CHECK(world.venues[0].fomite_history[0].front().amount ==
        doctest::Approx(3.0));
}

TEST_CASE("Fomite Transmission To Susceptible") {
  WorldState world;
  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);

  // Only a susceptible person (no infector present)
  world.people.emplace_back();
  Person& p1 = world.people[0];
  p1.id = 0;
  p1.age = 30;
  p1.sex = Sex::MALE;
  p1.geo_unit_id = -1;
  world.buildIndices();

  // High fomite infectiousness + high susceptibility multiplier to make
  // infection near-certain
  Disease disease = makeDiseaseWithFomite(
      1.0,   // deposition_rate (irrelevant here)
      10.0,  // fomite infectiousness (high)
      50.0   // fomite susceptibility multiplier (very high)
  );
  initVenueFomiteHistory(world, 1);

  // Pre-load a large fomite deposit
  world.venues[0].fomite_history[0].push_back({4.0, 100.0});

  std::vector<PersonLocation> locs;
  PersonLocation loc;
  loc.person_id = 0;
  loc.venue_id = 0;
  loc.subset_index = -1;
  loc.activity_index = -1;
  loc.encounter_type_id = 255;
  loc.person_array_index = 0;
  locs.push_back(loc);

  ContactMatrixConfig cm_config;
  SimulationConfig sim_config;
  ParallelConfig parallel_config;

  // Run multiple trials to ensure infection happens
  bool infected = false;
  for (int trial = 0; trial < 20 && !infected; ++trial) {
    world.people[0].infection.reset();
    cm_config.allow_default_matrix = true;
    finalizeContactMatrices(cm_config, world, disease);
    InteractionManager im(world, cm_config, sim_config, parallel_config,
                          &disease, nullptr);
    im.processTransmissions(locs, 5.0, 8.0, nullptr);
    if (world.people[0].infection) infected = true;
  }
  CHECK(infected);
}

TEST_CASE("Fomite-Only Infection (No Direct Infector Present)") {
  // Verify that fomite transmission works even when no infectious person is
  // present at the venue — only deposited material drives infection.
  WorldState world;
  Venue venue;
  venue.id = 0;
  venue.type_id = 0;
  world.venue_type_names = {"office"};
  world.venues.push_back(venue);

  world.people.emplace_back();
  Person& p1 = world.people[0];
  p1.id = 0;
  p1.age = 30;
  p1.sex = Sex::MALE;
  p1.geo_unit_id = -1;
  world.buildIndices();

  Disease disease = makeDiseaseWithFomite(1.0, 10.0, 50.0);
  initVenueFomiteHistory(world, 1);

  // Pre-load fomite deposit (no infectious person at venue)
  world.venues[0].fomite_history[0].push_back({4.0, 100.0});

  std::vector<PersonLocation> locs;
  PersonLocation loc;
  loc.person_id = 0;
  loc.venue_id = 0;
  loc.subset_index = -1;
  loc.activity_index = -1;
  loc.encounter_type_id = 255;
  loc.person_array_index = 0;
  locs.push_back(loc);

  ContactMatrixConfig cm_config;
  SimulationConfig sim_config;
  ParallelConfig parallel_config;

  bool infected = false;
  for (int trial = 0; trial < 50 && !infected; ++trial) {
    world.people[0].infection.reset();
    cm_config.allow_default_matrix = true;
    finalizeContactMatrices(cm_config, world, disease);
    InteractionManager im(world, cm_config, sim_config, parallel_config,
                          &disease, nullptr);
    im.processTransmissions(locs, 5.0, 8.0, nullptr);
    if (world.people[0].infection) {
      infected = true;
    }
  }
  // Person should be infected purely via fomite, with no direct infector
  CHECK(infected);
}

// =============================================================================
// Regression test: stale bin data must not leak across venue calls
//
// Reproduces the scenario where InteractionManager processes venues with
// different bin counts through a shared bins_buffer_. If bins are not fully
// cleared between calls, higher-indexed bins retain stale infectious data
// from earlier venues and corrupt later venues' transmission calculations.
//
// Sequence: venue 0 (3 bins, infectious person in bin 2)
//         → venue 1 (1 bin, healthy person only)
//         → venue 2 (3 bins, susceptible person in bin 2, NO infection source)
//
// If bin 2 retains stale infectiousness from venue 0, venue 2's susceptible
// person will see phantom infectious pressure and may become falsely infected.
// =============================================================================
TEST_CASE("No Stale Bin Data Across Venues With Different Bin Counts") {
  // --- World setup: 3 venues, 2 venue types ---
  WorldState world;
  // venue type 0 = "school" (will have 3 age-based bins)
  // venue type 1 = "home"   (no contact matrix → 1 bin)
  world.venue_type_names = {"school", "home"};

  Venue v0;
  v0.id = 0;
  v0.type_id = 0;  // school (3 bins)
  world.venues.push_back(v0);

  Venue v1;
  v1.id = 1;
  v1.type_id = 1;  // home (1 bin)
  world.venues.push_back(v1);

  Venue v2;
  v2.id = 2;
  v2.type_id = 0;  // school (3 bins)
  world.venues.push_back(v2);

  // Person 0: age 70 → bin 2, infectious, at venue 0
  // Person 1: age 40 → bin 1, healthy, at venue 1 (buffer venue)
  // Person 2: age 70 → bin 2, susceptible, at venue 2
  // Person 3: age 25 → bin 0, at venue 2 (needed for bin_size denominator)
  for (int i = 0; i < 4; ++i) {
    world.people.emplace_back();
    world.people[i].id = i;
    world.people[i].sex = Sex::MALE;
    world.people[i].geo_unit_id = -1;
  }
  world.people[0].age = 70;
  world.people[1].age = 40;
  world.people[2].age = 70;
  world.people[3].age = 25;
  world.buildIndices();

  // --- Disease: direct-contact mode with high infectiousness ---
  Disease disease = makeDiseaseWithFomite(
      1.0,  // deposition_rate (irrelevant — no fomite history at venue 2)
      0.0,  // fomite_infectiousness = 0 (isolate direct-contact path)
      0.0   // fomite_susc_mult = 0
  );
  initVenueFomiteHistory(world, 1);

  // Infect person 0 (at venue 0)
  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 42, nullptr, "school", 0);

  // --- Contact matrix for "school": 3 age-based bins with high contacts ---
  ContactMatrixConfig cm_config;
  cm_config.default_beta = 1.0;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{10.0}};
  cm_config.default_matrix = default_contact_matrix;

  ContactMatrix school_matrix;
  school_matrix.bins = {"0-30", "31-60", "61+"};
  // 3x3 contact matrix with high within-bin and cross-bin contacts
  school_matrix.contacts = {
      {10.0, 5.0, 5.0},  // from bin 0
      {5.0, 10.0, 5.0},  // from bin 1
      {5.0, 5.0, 10.0},  // from bin 2
  };
  school_matrix.characteristic_time = 1.0;
  cm_config.matrices["school"] = std::move(school_matrix);
  cm_config.betas["school"] = 1.0;
  // No matrix for "home" → defaults to 1 bin

  cm_config.resolve(world);

  // --- Locations: all 4 people at their respective venues ---
  std::vector<PersonLocation> locs;
  auto makeLoc = [](PersonId pid, VenueId vid, size_t arr_idx) {
    PersonLocation loc;
    loc.person_id = pid;
    loc.venue_id = vid;
    loc.subset_index = -1;
    loc.activity_index = -1;
    loc.encounter_type_id = 255;
    loc.person_array_index = arr_idx;
    return loc;
  };
  locs.push_back(makeLoc(0, 0, 0));  // person 0 → venue 0 (school, 3 bins)
  locs.push_back(makeLoc(1, 1, 1));  // person 1 → venue 1 (home, 1 bin)
  locs.push_back(makeLoc(2, 2, 2));  // person 2 → venue 2 (school, 3 bins)
  locs.push_back(makeLoc(3, 2, 3));  // person 3 → venue 2 (school, 3 bins)

  // --- Run multiple trials: if stale data leaks, high beta should cause
  //     near-certain false infection of person 2 ---
  SimulationConfig sim_config;
  ParallelConfig parallel_config;

  int false_infections = 0;
  const int num_trials = 50;
  for (int trial = 0; trial < num_trials; ++trial) {
    // Reset person 2 (susceptible target)
    world.people[2].infection.reset();
    // Reset person 3 (bystander)
    world.people[3].infection.reset();
    // Ensure person 0 stays infectious
    if (!world.people[0].infection) {
      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], 42, nullptr, "school", 0);
    }
    // Clear fomite history to isolate direct-contact path
    for (auto& venue : world.venues) {
      for (auto& hist : venue.fomite_history) hist.clear();
    }

    cm_config.allow_default_matrix = true;
    finalizeContactMatrices(cm_config, world, disease);
    InteractionManager im(world, cm_config, sim_config, parallel_config,
                          &disease, nullptr);
    im.processTransmissions(locs, 5.0, 8.0, nullptr);

    if (world.people[2].infection) {
      false_infections++;
    }
  }

  // Person 2 is at venue 2 which has NO infectious person and NO fomite
  // history. Any infection of person 2 is a false positive caused by stale
  // bin data leaking from venue 0 through the shared bins_buffer_.
  CHECK(false_infections == 0);
}
