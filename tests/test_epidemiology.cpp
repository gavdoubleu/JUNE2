#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/epidemiology.h"
#include "transmission_fixtures.h"

using namespace june;

TEST_CASE("Epidemiology - Fomite Pruning") {
  WorldState world;

  world.venues.emplace_back();
  world.venues.back().id = 0;
  world.venues.back().geo_unit_id = 0;
  world.venues.back().type_id = 0;

  // Build a minimal disease with a fomite config (max_age = 2.0 days)
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto constant_curve = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["mild"] = constant_curve;
  trans.symptom_id_curves = {nullptr, constant_curve};

  TransmissionMode direct_mode;
  direct_mode.name = "direct";
  direct_mode.type = TransmissionModeType::Fomite;
  direct_mode.symptom_curves = {nullptr, constant_curve};
  FomiteConfig fcfg;
  fcfg.mode_index = 0;
  fcfg.max_age = 2.0;
  fcfg.infectiousness_curve = std::make_shared<ConstantCurve>(1.0);
  direct_mode.config = std::move(fcfg);
  trans.modes.push_back(std::move(direct_mode));

  std::vector<SymptomTag> symptom_tags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 10.0}}}});
  trajectories.push_back(td);

  Disease disease("TestDisease", symptom_tags, {}, trajectories, {}, trans);
  Epidemiology epi(world, &disease);

  // Initialize venue fomite history for 1 mode
  world.venues[0].fomite_history.assign(1, {});

  // Pre-populate with old event (time=0.0)
  world.venues[0].fomite_history[0].push_back({0.0, 5.0});
  CHECK(world.venues[0].fomite_history[0].size() == 1);

  // current_time=3.0 days → age = 3.0 days > max_age(2.0 days) → should be
  // pruned
  epi.updateVenueFomites(3.0, 1.0);
  CHECK(world.venues[0].fomite_history[0].empty());

  // Add a recent event (time=2.5 days) and call again at time=3.0
  // age = 0.5 days < max_age(2.0 days) → should survive
  world.venues[0].fomite_history[0].push_back({2.5, 3.0});
  epi.updateVenueFomites(3.0, 1.0);
  CHECK(world.venues[0].fomite_history[0].size() == 1);
  CHECK(world.venues[0].fomite_history[0].front().amount ==
        doctest::Approx(3.0));
}

TEST_CASE("Epidemiology - updateInfectionStates progression") {
  WorldState world;

  // Build disease with 3-stage trajectory: mild -> severe -> healthy
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  trans.natural_immunity.level = 0.95;
  trans.natural_immunity.waning_rate = 0.001;

  auto constant_curve = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["mild"] = constant_curve;
  trans.stage_curves["severe"] = constant_curve;
  trans.symptom_id_curves = {nullptr, constant_curve, constant_curve};

  std::vector<SymptomTag> symptom_tags = {
      {"healthy", -1, 0}, {"mild", 1, 1}, {"severe", 2, 2}};

  DiseaseStageSettings stage_settings;
  stage_settings.recovered_stages = {"healthy"};

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 3.0}}}});  // 3 days mild
  td.stages.push_back(
      {"severe", {"constant", {{"value", 2.0}}}});  // 2 days severe
  td.stages.push_back(
      {"healthy", {"constant", {{"value", 100.0}}}});  // recover
  trajectories.push_back(td);

  Disease disease("TestDisease", symptom_tags, stage_settings, trajectories, {},
                  trans);
  Epidemiology epi(world, &disease);

  // Create 3 people, infect at staggered times
  for (int i = 0; i < 3; ++i) {
    auto& p = world.people.emplace_back();
    p.id = i;
    p.age = 30.0f;
    p.geo_unit_id = -1;
  }
  world.buildIndices();

  // Infect person 0 at t=0, person 1 at t=1, person 2 at t=2
  for (int i = 0; i < 3; ++i) {
    world.people[i].infection = std::make_unique<Infection>(
        &disease, static_cast<double>(i), &world.people[i], 100 + i,
        kNoTransmissionContext, nullptr, "office", 0);
    epi.trackInfection(i);
  }

  std::vector<PersonLocation> empty_locs;

  SUBCASE("All 3 active after infection") {
    CHECK(epi.getActiveInfections().size() == 3);
  }

  SUBCASE("Person 0 transitions through stages correctly") {
    // At t=1.0: person 0 should be in "mild" (infected at t=0, mild lasts 3
    // days)
    CHECK(world.people[0].infection->getCurrentSymptom(1.0) == "mild");

    // At t=4.0: person 0 should be in "severe" (mild ends at t=3)
    CHECK(world.people[0].infection->getCurrentSymptom(4.0) == "severe");

    // At t=6.0: person 0 should be "healthy" (severe ends at t=5)
    CHECK(world.people[0].infection->isRecovered(6.0));
  }

  SUBCASE("updateInfectionStates removes recovered people") {
    // Run until person 0 should have recovered (t=6)
    epi.updateInfectionStates(6.0, empty_locs);
    // Person 0 infected at t=0, recovers at t=5 — should be removed
    CHECK(epi.getActiveInfections().count(0) == 0);
    CHECK(world.people[0].infection == nullptr);
    // Person 1 infected at t=1, recovers at t=6 — should also be removed
    CHECK(epi.getActiveInfections().count(1) == 0);
    // Person 2 infected at t=2, recovers at t=7 — still active
    CHECK(epi.getActiveInfections().count(2) == 1);
  }

  SUBCASE("Recovery sets natural immunity") {
    epi.updateInfectionStates(6.0, empty_locs);
    // Person 0 should now have natural immunity
    CHECK(world.people[0].immunity.natural_level == doctest::Approx(0.95));
    CHECK(world.people[0].immunity.natural_acquisition_time ==
          doctest::Approx(6.0));
    CHECK(world.people[0].immunity.natural_waning_rate ==
          doctest::Approx(0.001));
  }

  SUBCASE("All recovered after enough time") {
    epi.updateInfectionStates(10.0, empty_locs);
    CHECK(epi.getActiveInfections().empty());
    for (int i = 0; i < 3; ++i) {
      CHECK(world.people[i].infection == nullptr);
      CHECK(world.people[i].immunity.natural_level > 0.0);
    }
  }
}

TEST_CASE("Epidemiology - updateInfectionStates with fatal trajectory") {
  WorldState world;

  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto constant_curve = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["mild"] = constant_curve;
  trans.stage_curves["dead"] = constant_curve;
  trans.symptom_id_curves = {nullptr, constant_curve, constant_curve};

  std::vector<SymptomTag> symptom_tags = {
      {"healthy", -1, 0}, {"mild", 1, 1}, {"dead", 3, 2}};

  DiseaseStageSettings stage_settings;
  stage_settings.fatality_stages = {"dead"};

  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 2.0}}}});
  td.stages.push_back({"dead", {"constant", {{"value", 100.0}}}});
  trajectories.push_back(td);

  Disease disease("TestDisease", symptom_tags, stage_settings, trajectories, {},
                  trans);
  Epidemiology epi(world, &disease);

  auto& p = world.people.emplace_back();
  p.id = 0;
  p.age = 80.0f;
  p.geo_unit_id = -1;
  world.buildIndices();

  p.infection = std::make_unique<Infection>(
      &disease, 0.0, &p, 42, kNoTransmissionContext, nullptr, "office", 0);
  epi.trackInfection(0);

  std::vector<PersonLocation> empty_locs;
  epi.updateInfectionStates(5.0, empty_locs);

  CHECK(epi.getActiveInfections().empty());
  CHECK(p.is_dead);
  CHECK(p.death_time >= 0.0);
  CHECK(p.infection == nullptr);
}

TEST_CASE("Epidemiology - Track Infection") {
  WorldState world;

  // Build a minimal disease
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto constant_curve = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["mild"] = constant_curve;
  trans.symptom_id_curves = {nullptr, constant_curve};

  std::vector<SymptomTag> symptom_tags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 10.0}}}});
  trajectories.push_back(td);

  Disease disease("TestDisease", symptom_tags, {}, trajectories, {}, trans);
  Epidemiology epi(world, &disease);

  world.people.emplace_back();
  world.people.back().id = 1;
  world.person_index[1] = 0;

  epi.trackInfection(1);
  CHECK(epi.getActiveInfections().size() == 1);
  CHECK(epi.getActiveInfections().count(1) == 1);

  epi.untrackInfection(1);
  CHECK(epi.getActiveInfections().size() == 0);
}
