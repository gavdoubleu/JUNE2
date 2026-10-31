#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "loaders/config_loader_detail.h"
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

static WorldState buildListPropertyWorld() {
  WorldState world;
  world.person_property_names = {"comorbidities"};
  world.person_property_list_value_registries["comorbidities"] = {
      {"crd", "cancer"}, {"flu"}, {}};

  for (int i = 0; i < 4; ++i) {
    Person& person = world.people.emplace_back();
    person.id = i;
    person.properties_start = static_cast<uint32_t>(i);
    person.properties_count = 1;
    world.person_properties.push_back(i < 3 ? i : -1);
  }
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

TEST_CASE("list-valued properties use exact membership") {
  WorldState world = buildListPropertyWorld();
  const auto& with_cancer = world.people[0];
  const auto& empty = world.people[2];
  const auto& missing = world.people[3];

  auto property = world.getPersonProperty(with_cancer, "comorbidities");
  REQUIRE(property.has_value());
  REQUIRE(std::holds_alternative<std::vector<std::string>>(*property));
  CHECK(std::get<std::vector<std::string>>(*property) ==
        std::vector<std::string>{"crd", "cancer"});

  SelectionCriterion equals_cancer =
      makeCriterion("properties.comorbidities", "==", "cancer");
  SelectionCriterion contains_cancer =
      makeCriterion("properties.comorbidities", "contains", "cancer");
  SelectionCriterion contains_fragment =
      makeCriterion("properties.comorbidities", "contains", "can");
  SelectionCriterion not_cancer =
      makeCriterion("properties.comorbidities", "!=", "cancer");
  for (SelectionCriterion* criterion : {&equals_cancer, &contains_cancer,
                                         &contains_fragment, &not_cancer})
    criterion->resolveOrThrow(world, "test");

  CHECK(equals_cancer.evaluate(with_cancer, &world));
  CHECK(contains_cancer.evaluate(with_cancer, &world));
  CHECK_FALSE(contains_fragment.evaluate(with_cancer, &world));
  CHECK_FALSE(not_cancer.evaluate(with_cancer, &world));
  CHECK(not_cancer.evaluate(empty, &world));
  CHECK_FALSE(not_cancer.evaluate(missing, &world));
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

TEST_CASE("infection context criteria refuse other operators and non-strings") {
  WorldState world = buildSinglePersonWorld();
  SelectionCriterion ordered =
      makeCriterion("infector_symptom", ">", std::string("bubonic"));
  CHECK_THROWS_WITH(ordered.resolveOrThrow(world, "test"),
                    doctest::Contains("supports only == !="));

  SelectionCriterion numeric = makeCriterion("transmission_mode", "==", 3);
  CHECK_THROWS_WITH(numeric.resolveOrThrow(world, "test"),
                    doctest::Contains("compares against a name"));

  SelectionCriterion valid =
      makeCriterion("transmission_mode", "!=", std::string("respiratory"));
  CHECK_NOTHROW(valid.resolveOrThrow(world, "test"));
}

// Person 0 lives in a household, person 1 in a care home, and person 2 is
// registered to no residence at all.
static WorldState buildResidenceTypeWorld() {
  WorldState world;
  world.activity_names = {"residence"};
  world.venue_type_names = {"household", "care_home"};
  for (VenueId id = 0; id < 2; ++id) {
    Venue venue{};
    venue.id = id;
    venue.type_id = static_cast<uint8_t>(id);
    venue.geo_unit_id = 0;
    world.venues.push_back(venue);
  }
  for (PersonId id = 0; id < 3; ++id) {
    Person& person = world.people.emplace_back();
    person.id = id;
    person.age = 30.0f;
    person.geo_unit_id = -1;
    person.activity_meta_start =
        static_cast<uint32_t>(world.activity_meta.size());
    if (id == 2) continue;
    world.activity_meta.push_back(
        {0, static_cast<uint32_t>(world.activity_venues.size()), 1});
    world.activity_venues.push_back({static_cast<VenueId>(id), 0});
    ++person.activity_meta_count;
  }
  world.buildIndices();
  return world;
}

TEST_CASE("not-equal on an activity's venue type excludes the type itself") {
  const WorldState world = buildResidenceTypeWorld();
  const auto matches = [&](const char* op, const char* type, PersonId id) {
    SelectionCriterion criterion =
        makeCriterion("activities.residence.venue_type", op, std::string(type));
    return criterion.evaluate(world.people[id], &world);
  };

  SUBCASE("a type the world has") {
    CHECK_FALSE(matches("==", "care_home", 0));
    CHECK(matches("==", "care_home", 1));
    CHECK_FALSE(matches("==", "care_home", 2));
    CHECK(matches("!=", "care_home", 0));
    // A care home resident does not live somewhere that is not a care home.
    CHECK_FALSE(matches("!=", "care_home", 1));
    // With no residence there is nothing to compare, so neither operator holds.
    CHECK_FALSE(matches("!=", "care_home", 2));
  }

  SUBCASE("a type the world does not have") {
    CHECK_FALSE(matches("==", "hospice", 0));
    CHECK_FALSE(matches("==", "hospice", 1));
    CHECK(matches("!=", "hospice", 0));
    CHECK(matches("!=", "hospice", 1));
    CHECK_FALSE(matches("!=", "hospice", 2));
  }
}

TEST_CASE("is_alive written in YAML selects the living or the dead") {
  WorldState world = buildSinglePersonWorld();
  Person& dead = world.people.emplace_back();
  dead.id = 1;
  dead.age = 70.0f;
  dead.sex = Sex::MALE;
  dead.geo_unit_id = 10;
  dead.is_dead = true;
  world.buildIndices();
  const Person& living = world.people[0];

  auto fromYaml = [&](const std::string& op, const std::string& value) {
    std::vector<SelectionCriterion> criteria;
    config_detail::parseSelectionCriteria(
        YAML::Load("[{property: is_alive, operator: \"" + op +
                   "\", value: " + value + "}]"),
        criteria);
    REQUIRE(criteria.size() == 1);
    return criteria.front();
  };
  auto selects = [&](const std::string& op, const std::string& value,
                     const Person& person) {
    SelectionCriterion criterion = fromYaml(op, value);
    criterion.resolveOrThrow(world, "test");
    return criterion.evaluate(person, &world);
  };

  CHECK(selects("==", "true", living));
  CHECK_FALSE(selects("==", "true", dead));
  CHECK(selects("!=", "false", living));
  CHECK_FALSE(selects("!=", "false", dead));
  CHECK(selects("==", "false", dead));
  CHECK_FALSE(selects("==", "false", living));
  CHECK(selects("!=", "true", dead));
  CHECK_FALSE(selects("!=", "true", living));

  SelectionCriterion greater = fromYaml(">", "true");
  CHECK_THROWS_WITH(
      greater.resolveOrThrow(world, "test"),
      doctest::Contains("'is_alive' compares with == or != against true or false"));
  SelectionCriterion maybe = fromYaml("==", "maybe");
  CHECK_THROWS_WITH(
      maybe.resolveOrThrow(world, "test"),
      doctest::Contains("'is_alive' compares with == or != against true or false"));
}
