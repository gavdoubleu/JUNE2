#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <string>
#include <vector>

#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "utils/filtering.h"

using namespace june;

namespace {

// One nation over one fine unit, one person in it: enough for every path that
// resolves against a world without needing venues or activities.
static WorldState buildSinglePersonWorld() {
  WorldState world;
  world.geo_level_names = {"XLGU", "SGU"};

  GeographicalUnit nation;
  nation.id = 0;
  nation.name = "Wales";
  nation.level_id = 0;
  nation.parent_id = -1;
  world.geo_units.push_back(nation);

  GeographicalUnit output_area;
  output_area.id = 10;
  output_area.name = "W00000001";
  output_area.level_id = 1;
  output_area.parent_id = 0;
  world.geo_units.push_back(output_area);

  Person& person = world.people.emplace_back();
  person.id = 0;
  person.age = 59.5f;
  person.sex = Sex::FEMALE;
  person.geo_unit_id = 10;

  world.buildIndices();
  return world;
}

static SelectionCriterion makeCriterion(const std::string& property_path,
                                        const std::string& operator_type,
                                        PropertyValue value) {
  SelectionCriterion criterion;
  criterion.property_path = property_path;
  criterion.operator_type = operator_type;
  criterion.value = std::move(value);
  return criterion;
}

static bool evaluateResolved(SelectionCriterion criterion,
                             const WorldState& world) {
  criterion.resolveOrThrow(world, "test");
  return criterion.evaluate(world.people.front(), &world);
}

}  // namespace

TEST_CASE("numeric operators compare an int threshold against fractional age") {
  WorldState world = buildSinglePersonWorld();
  CHECK(evaluateResolved(makeCriterion("age", ">", 59), world));
  CHECK_FALSE(evaluateResolved(makeCriterion("age", "<", 59), world));
  CHECK(evaluateResolved(makeCriterion("age", ">=", 59), world));
  CHECK_FALSE(evaluateResolved(makeCriterion("age", "<=", 59), world));
  CHECK(evaluateResolved(makeCriterion("age", "<=", 60), world));
  CHECK_FALSE(evaluateResolved(makeCriterion("age", "==", 59), world));
  CHECK(evaluateResolved(makeCriterion("age", "!=", 59), world));
}

TEST_CASE("equality on sex goes through the interned code") {
  WorldState world = buildSinglePersonWorld();
  CHECK(evaluateResolved(makeCriterion("sex", "==", std::string("female")),
                         world));
  CHECK_FALSE(
      evaluateResolved(makeCriterion("sex", "==", std::string("M")), world));
  CHECK(
      evaluateResolved(makeCriterion("sex", "!=", std::string("male")), world));
}

TEST_CASE("in matches a geographical unit id list") {
  WorldState world = buildSinglePersonWorld();
  CHECK(evaluateResolved(
      makeCriterion("geo_unit_id", "in", std::vector<int32_t>{3, 10}), world));
  CHECK_FALSE(evaluateResolved(
      makeCriterion("geo_unit_id", "in", std::vector<int32_t>{3, 11}), world));
}

TEST_CASE("contains matches a substring of a string value") {
  WorldState world = buildSinglePersonWorld();
  CHECK(evaluateResolved(makeCriterion("sex", "contains", std::string("fem")),
                         world));
  CHECK_FALSE(evaluateResolved(
      makeCriterion("sex", "contains", std::string("xyz")), world));
}

TEST_CASE("boolean predicates support == and != only") {
  WorldState world = buildSinglePersonWorld();
  CHECK(evaluateResolved(makeCriterion("is_alive", "==", true), world));
  CHECK(evaluateResolved(makeCriterion("is_alive", "!=", false), world));

  SelectionCriterion greater = makeCriterion("is_alive", ">", true);
  greater.resolve(world);
  CHECK_FALSE(greater.evaluate(world.people.front(), &world));
}

TEST_CASE("an unsupported operator matches nobody and is refused at load") {
  WorldState world = buildSinglePersonWorld();
  SelectionCriterion unsupported = makeCriterion("age", "=~", 59);
  unsupported.resolve(world);
  CHECK_FALSE(unsupported.evaluate(world.people.front(), &world));
  CHECK_THROWS_WITH(unsupported.resolveOrThrow(world, "test"),
                    doctest::Contains("operator '=~' is not supported"));
}

TEST_CASE("ancestor geography honours == != and in") {
  WorldState world = buildSinglePersonWorld();
  CHECK(evaluateResolved(
      makeCriterion("geo_unit.XLGU", "==", std::string("Wales")), world));
  CHECK_FALSE(evaluateResolved(
      makeCriterion("geo_unit.XLGU", "!=", std::string("Wales")), world));
  CHECK(evaluateResolved(
      makeCriterion("geo_unit.XLGU", "in", std::vector<std::string>{"Wales"}),
      world));
}

TEST_CASE("a != expression parsed from text excludes the named value") {
  WorldState world = buildSinglePersonWorld();
  std::vector<SelectionCriterion> criteria =
      filtering::parseConjunctiveExpression("age>=18 AND sex!=female");
  for (SelectionCriterion& criterion : criteria)
    criterion.resolveOrThrow(world, "test");
  CHECK_FALSE(
      filtering::matchesCriteria(world.people.front(), &world, criteria));
}

namespace {

// Two outcome rows keyed on infection context, then a catch-all on age.
static OutcomeRates buildContextOutcomeRates() {
  const std::vector<std::string> headers = {
      "filter.age", "filter.infector_symptom", "filter.transmission_mode",
      "death"};
  const std::vector<std::pair<int, std::string>> filter_columns =
      filtering::findFilterColumns(headers);

  auto addRow = [&](const std::vector<std::string>& fields, double death) {
    OutcomeRow row;
    row.criteria = filtering::parseCriteriaFromRow(fields, filter_columns);
    row.probabilities["death"] = death;
    return row;
  };

  OutcomeRates rates;
  rates.rows.push_back(addRow({"18-99", "pneumonic", "respiratory"}, 0.9));
  rates.rows.push_back(addRow({"18-99", "bubonic", "animal_bite"}, 0.4));
  rates.rows.push_back(addRow({"0-17", "", ""}, 0.1));
  return rates;
}

}  // namespace

TEST_CASE("outcome rates keyed on infection context resolve and pick rows") {
  WorldState world = buildSinglePersonWorld();
  OutcomeRates rates = buildContextOutcomeRates();
  REQUIRE_NOTHROW(rates.resolve(world));

  const Person& adult = world.people.front();
  CHECK(rates.getRate(adult, &world, "death",
                      InfectionContext{"pneumonic", "respiratory"}) ==
        doctest::Approx(0.9));
  CHECK(rates.getRate(adult, &world, "death",
                      InfectionContext{"bubonic", "animal_bite"}) ==
        doctest::Approx(0.4));
  CHECK(rates.getRate(adult, &world, "death",
                      InfectionContext{"bubonic", "respiratory"}) == 0.0);
}

TEST_CASE("an empty infection context fails both == and !=") {
  WorldState world = buildSinglePersonWorld();
  const Person& person = world.people.front();
  const InfectionContext seeded{};
  const InfectionContext bitten{"bubonic", "animal_bite"};

  std::vector<SelectionCriterion> not_respiratory = {
      makeCriterion("transmission_mode", "!=", std::string("respiratory"))};
  std::vector<SelectionCriterion> is_bubonic = {
      makeCriterion("infector_symptom", "==", std::string("bubonic"))};
  for (SelectionCriterion& criterion : not_respiratory)
    criterion.resolveOrThrow(world, "test");
  for (SelectionCriterion& criterion : is_bubonic)
    criterion.resolveOrThrow(world, "test");

  CHECK_FALSE(
      filtering::matchesCriteria(person, &world, not_respiratory, seeded));
  CHECK_FALSE(filtering::matchesCriteria(person, &world, is_bubonic, seeded));
  CHECK(filtering::matchesCriteria(person, &world, not_respiratory, bitten));
  CHECK(filtering::matchesCriteria(person, &world, is_bubonic, bitten));
}

TEST_CASE("unresolved criteria evaluate without a world") {
  WorldState world = buildSinglePersonWorld();
  const Person& person = world.people.front();
  const std::vector<SelectionCriterion> criteria = {
      makeCriterion("age", ">=", 18),
      makeCriterion("transmission_mode", "==", std::string("respiratory"))};

  CHECK(filtering::matchesCriteria(person, nullptr, criteria,
                                   InfectionContext{"", "respiratory"}));
  CHECK_FALSE(filtering::matchesCriteria(person, nullptr, criteria,
                                         InfectionContext{"", "animal_bite"}));
}

TEST_CASE("a world-dependent criterion without a world matches nobody") {
  WorldState world = buildSinglePersonWorld();
  const std::vector<SelectionCriterion> criteria = {
      makeCriterion("properties.occupation", "==", std::string("farmer"))};
  CHECK_FALSE(
      filtering::matchesCriteria(world.people.front(), nullptr, criteria));
}
