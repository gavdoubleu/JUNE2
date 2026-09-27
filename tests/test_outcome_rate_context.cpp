#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <stdexcept>
#include <string>

#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/infection_seed.h"
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
static Disease buildDisease(std::vector<OutcomeRow> rows) {
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
  return Disease("Plague", {recovered, pneumonic}, DiseaseStageSettings{}, {},
                 rates, transmission);
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
