#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <stdexcept>
#include <string>

#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/infection_seed.h"
#include "epidemiology/transmission/infection_context.h"
#include "loaders/disease_loader.h"

using namespace june;

namespace {

static WorldState buildOnePersonWorld() {
  WorldState world;
  Person& person = world.people.emplace_back();
  person.id = 0;
  person.age = 30;
  world.buildIndices();
  return world;
}

static SelectionCriterion contextCriterion(const std::string& fact,
                                           const std::string& value) {
  SelectionCriterion criterion;
  criterion.property_path = fact;
  criterion.operator_type = "==";
  criterion.value = value;
  return criterion;
}

// A plague-shaped disease whose outcome table is the given rows.
static Disease buildDisease(
    std::vector<OutcomeRow> rows,
    std::vector<TrajectoryDefinition> trajectories = {}) {
  SymptomTag recovered{.name = "recovered", .value = -3, .id = 0};
  SymptomTag pneumonic{.name = "primary_pneumonic", .value = 2, .id = 1};

  TransmissionParams transmission;
  TransmissionMode animal_bite;
  animal_bite.name = "animal_bite";
  TransmissionMode respiratory;
  respiratory.name = "respiratory";
  transmission.modes = {animal_bite, respiratory};

  OutcomeRates rates;
  rates.rows = std::move(rows);
  return Disease("Plague", {recovered, pneumonic}, DiseaseStageSettings{},
                 trajectories, rates, transmission);
}

static OutcomeRow rowWith(SelectionCriterion criterion) {
  OutcomeRow row;
  row.criteria = {std::move(criterion)};
  row.probabilities = {{"mild", 1.0}};
  return row;
}

}  // namespace

TEST_CASE("an outcome row filtering on a known transmission mode resolves") {
  WorldState world = buildOnePersonWorld();
  Disease disease =
      buildDisease({rowWith(contextCriterion("transmission_mode", "respiratory"))});
  CHECK_NOTHROW(disease.resolve(world));
}

TEST_CASE("an outcome row filtering on a known infector symptom resolves") {
  WorldState world = buildOnePersonWorld();
  Disease disease = buildDisease(
      {rowWith(contextCriterion("infector_symptom", "primary_pneumonic"))});
  CHECK_NOTHROW(disease.resolve(world));
}

TEST_CASE("an outcome row naming a transmission mode the disease lacks is refused, naming the row and value") {
  WorldState world = buildOnePersonWorld();
  OutcomeRow default_row;
  default_row.probabilities = {{"mild", 1.0}};
  Disease disease = buildDisease(
      {default_row, rowWith(contextCriterion("transmission_mode", "rat_flee_bite"))});
  std::string message;
  try {
    disease.resolve(world);
  } catch (const std::runtime_error& error) {
    message = error.what();
  }
  CHECK(message.find("row 1") != std::string::npos);
  CHECK(message.find("rat_flee_bite") != std::string::npos);
}

TEST_CASE("an outcome row naming a symptom the disease lacks is refused, naming the row and value") {
  WorldState world = buildOnePersonWorld();
  Disease disease = buildDisease(
      {rowWith(contextCriterion("infector_symptom", "primary_pneumonik"))});
  std::string message;
  try {
    disease.resolve(world);
  } catch (const std::runtime_error& error) {
    message = error.what();
  }
  CHECK(message.find("row 0") != std::string::npos);
  CHECK(message.find("primary_pneumonik") != std::string::npos);
}

TEST_CASE("a context filter outside an outcome table is still refused") {
  WorldState world = buildOnePersonWorld();
  InfectionSeedEvent seed;
  seed.name = "seed_filtering_on_mode";
  seed.attribute_filters = {contextCriterion("transmission_mode", "respiratory")};
  InfectionSeedConfig config;
  config.seeds = {seed};
  CHECK_THROWS_AS(config.resolve(world), std::runtime_error);
}

TEST_CASE("the plague outcome table, filtered on mode and infector symptom, resolves") {
  WorldState world = buildOnePersonWorld();
  Disease disease =
      DiseaseLoader::loadFromYAML("configs/config_plague/disease_plague.yaml");
  CHECK_NOTHROW(disease.resolve(world));
}

TEST_CASE(
    "an absent infector symptom fails both == and != rows; a row with no "
    "criteria matches") {
  WorldState world = buildOnePersonWorld();
  Disease disease = buildDisease({});
  InfectionContext context = buildInfectionContext(
      TransmissionRecord{InfectionSource::Person, kNoSymptomId, 1}, disease);

  SelectionCriterion symptom_is =
      contextCriterion("infector_symptom", "recovered");
  SelectionCriterion symptom_is_not = symptom_is;
  symptom_is_not.operator_type = "!=";
  const Person& person = world.people[0];
  CHECK_FALSE(
      filtering::matchesCriteria(person, &world, {symptom_is}, context));
  CHECK_FALSE(
      filtering::matchesCriteria(person, &world, {symptom_is_not}, context));
  CHECK(filtering::matchesCriteria(person, &world, {}, context));
}

TEST_CASE("an absent transmission mode fails both == and != rows") {
  WorldState world = buildOnePersonWorld();
  Disease disease = buildDisease({});
  InfectionContext context = buildInfectionContext(
      TransmissionRecord{InfectionSource::Person, 1, kNoModeIndex}, disease);

  SelectionCriterion mode_is =
      contextCriterion("transmission_mode", "animal_bite");
  SelectionCriterion mode_is_not = mode_is;
  mode_is_not.operator_type = "!=";
  const Person& person = world.people[0];
  CHECK_FALSE(filtering::matchesCriteria(person, &world, {mode_is}, context));
  CHECK_FALSE(
      filtering::matchesCriteria(person, &world, {mode_is_not}, context));
}

// A one-stage trajectory, so the chosen trajectory is visible as its symptom.
static TrajectoryDefinition trajectoryInto(const std::string& selection_key,
                                           const std::string& symptom) {
  TrajectoryDefinition trajectory;
  trajectory.selection_key = selection_key;
  TrajectoryStage stage;
  stage.symptom_tag = symptom;
  stage.completion_time.type = "constant";
  stage.completion_time.params = {{"value", 1.0}};
  trajectory.stages = {stage};
  return trajectory;
}

TEST_CASE(
    "an undeclared seed takes the default row, not the (recovered, first mode) "
    "row") {
  WorldState world = buildOnePersonWorld();
  OutcomeRow recovered_bite_row;
  recovered_bite_row.criteria = {
      contextCriterion("infector_symptom", "recovered"),
      contextCriterion("transmission_mode", "animal_bite")};
  recovered_bite_row.probabilities = {{"recovered_bite", 1.0}};
  OutcomeRow default_row;
  default_row.probabilities = {{"default", 1.0}};
  Disease disease =
      buildDisease({recovered_bite_row, default_row},
                   {trajectoryInto("recovered_bite", "recovered"),
                    trajectoryInto("default", "primary_pneumonic")});

  Infection infection(
      &disease, 0.0, &world.people[0], 42,
      TransmissionRecord{InfectionSource::Seed, kNoSymptomId, kNoModeIndex},
      &world);
  CHECK(disease.getSymptomName(
            infection.getTrajectory().transitions.at(0).second) ==
        "primary_pneumonic");
}

TEST_CASE("an infection_source == seed row matches a seed but not a fomite") {
  WorldState world = buildOnePersonWorld();
  Disease disease = buildDisease({});
  InfectionContext seed_context = buildInfectionContext(
      TransmissionRecord{InfectionSource::Seed, kNoSymptomId, kNoModeIndex},
      disease);
  InfectionContext fomite_context = buildInfectionContext(
      TransmissionRecord{InfectionSource::Fomite, kNoSymptomId, 1}, disease);

  SelectionCriterion source_is_seed =
      contextCriterion("infection_source", "seed");
  const Person& person = world.people[0];
  CHECK(filtering::matchesCriteria(person, &world, {source_is_seed},
                                   seed_context));
  CHECK_FALSE(filtering::matchesCriteria(person, &world, {source_is_seed},
                                         fomite_context));
}

TEST_CASE(
    "an infection_source != seed row excludes seeds and matches Person and "
    "Fomite infections") {
  WorldState world = buildOnePersonWorld();
  Disease disease = buildDisease({});
  SelectionCriterion source_is_not_seed =
      contextCriterion("infection_source", "seed");
  source_is_not_seed.operator_type = "!=";
  const Person& person = world.people[0];
  auto matchesFrom = [&](InfectionSource source) {
    InfectionContext context = buildInfectionContext(
        TransmissionRecord{source, kNoSymptomId, kNoModeIndex}, disease);
    return filtering::matchesCriteria(person, &world, {source_is_not_seed},
                                      context);
  };
  CHECK_FALSE(matchesFrom(InfectionSource::Seed));
  CHECK(matchesFrom(InfectionSource::Person));
  CHECK(matchesFrom(InfectionSource::Fomite));
}

TEST_CASE(
    "a row combining source, mode and age matches only when all three hold") {
  WorldState world;
  Person& younger = world.people.emplace_back();
  younger.id = 0;
  younger.age = 30;
  Person& older = world.people.emplace_back();
  older.id = 1;
  older.age = 70;
  world.buildIndices();

  SelectionCriterion aged_sixty_or_over;
  aged_sixty_or_over.property_path = "age";
  aged_sixty_or_over.operator_type = ">=";
  aged_sixty_or_over.value = 60;
  OutcomeRow combined_row;
  combined_row.criteria = {contextCriterion("infection_source", "seed"),
                           contextCriterion("transmission_mode", "animal_bite"),
                           aged_sixty_or_over};
  combined_row.probabilities = {{"severe", 1.0}};
  OutcomeRow default_row;
  default_row.probabilities = {{"mild", 1.0}};
  Disease disease = buildDisease({combined_row, default_row});
  REQUIRE(disease.resolve(world).empty());

  const uint8_t animal_bite = 0;
  const uint8_t respiratory = 1;
  auto takesCombinedRow = [&](const Person& person, InfectionSource source,
                              uint8_t mode) {
    InfectionContext context = buildInfectionContext(
        TransmissionRecord{source, kNoSymptomId, mode}, disease);
    return disease.getOutcomeRates().getRate(person, &world, "severe",
                                             context) == 1.0;
  };
  CHECK(takesCombinedRow(world.people[1], InfectionSource::Seed, animal_bite));
  CHECK_FALSE(
      takesCombinedRow(world.people[0], InfectionSource::Seed, animal_bite));
  CHECK_FALSE(
      takesCombinedRow(world.people[1], InfectionSource::Person, animal_bite));
  CHECK_FALSE(
      takesCombinedRow(world.people[1], InfectionSource::Seed, respiratory));
}

TEST_CASE(
    "a seed-only row above the default takes seeds; a Person infection falls "
    "through to the default") {
  WorldState world = buildOnePersonWorld();
  OutcomeRow seed_row = rowWith(contextCriterion("infection_source", "seed"));
  seed_row.probabilities = {{"severe", 1.0}};
  OutcomeRow default_row;
  default_row.probabilities = {{"mild", 1.0}};
  Disease disease = buildDisease({seed_row, default_row});
  REQUIRE(disease.resolve(world).empty());

  auto severeRateFrom = [&](InfectionSource source) {
    InfectionContext context = buildInfectionContext(
        TransmissionRecord{source, kNoSymptomId, kNoModeIndex}, disease);
    return disease.getOutcomeRates().getRate(world.people[0], &world, "severe",
                                             context);
  };
  CHECK(severeRateFrom(InfectionSource::Seed) == 1.0);
  CHECK(severeRateFrom(InfectionSource::Person) == 0.0);
}

TEST_CASE(
    "an outcome row naming an unknown infection source is refused, naming the "
    "row and value") {
  WorldState world = buildOnePersonWorld();
  OutcomeRow default_row;
  default_row.probabilities = {{"mild", 1.0}};
  Disease disease = buildDisease(
      {default_row, rowWith(contextCriterion("infection_source", "sead"))});
  std::string message;
  try {
    disease.resolve(world);
  } catch (const std::runtime_error& error) {
    message = error.what();
  }
  CHECK(message.find("row 1") != std::string::npos);
  CHECK(message.find("sead") != std::string::npos);
}

TEST_CASE(
    "each infection_source name matches only infections from that source") {
  WorldState world = buildOnePersonWorld();
  Disease disease = buildDisease({});
  const std::vector<std::pair<std::string, InfectionSource>> sources = {
      {"person", InfectionSource::Person},
      {"fomite", InfectionSource::Fomite},
      {"compartmental", InfectionSource::Compartmental},
      {"seed", InfectionSource::Seed}};
  for (const auto& [row_name, row_source] : sources) {
    for (const auto& [context_name, context_source] : sources) {
      CAPTURE(row_name);
      CAPTURE(context_name);
      InfectionContext context = buildInfectionContext(
          TransmissionRecord{context_source, kNoSymptomId, kNoModeIndex},
          disease);
      CHECK(filtering::matchesCriteria(
                world.people[0], &world,
                {contextCriterion("infection_source", row_name)},
                context) == (row_source == context_source));
    }
  }
}
