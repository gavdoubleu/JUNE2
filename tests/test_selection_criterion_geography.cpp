#include <yaml-cpp/yaml.h>

#include <iostream>
#include <sstream>
#include <string>

#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/infection_seed.h"
#include "loaders/config_loader.h"
#include "loaders/config_loader_detail.h"
#include "loaders/policy_loader.h"
#include "test_utils.h"
#include "utils/filtering.h"

using namespace june;

namespace {

// Two-level geography mirroring the GB world's shape in miniature: named
// nations at the coarse level, output-area-like units at the fine level, and
// every person sitting at the fine level.
static WorldState buildNationWorld() {
  WorldState world;
  world.geo_level_names = {"XLGU", "SGU"};

  auto addUnit = [&world](GeoUnitId id, const std::string& name,
                          uint8_t level_id, GeoUnitId parent_id) {
    GeographicalUnit unit;
    unit.id = id;
    unit.name = name;
    unit.level_id = level_id;
    unit.parent_id = parent_id;
    world.geo_units.push_back(unit);
  };

  addUnit(0, "Scotland", 0, -1);
  addUnit(1, "London", 0, -1);
  addUnit(2, "Wales", 0, -1);
  addUnit(10, "S00000001", 1, 0);
  addUnit(11, "E00000001", 1, 1);
  addUnit(12, "W00000001", 1, 2);

  auto addPerson = [&world](PersonId id, GeoUnitId geo_unit_id) {
    Person& person = world.people.emplace_back();
    person.id = id;
    person.geo_unit_id = geo_unit_id;
  };

  addPerson(0, 10);  // Scotland
  addPerson(1, 11);  // London
  addPerson(2, 12);  // Wales

  world.buildIndices();
  return world;
}

// Resolve `criterion` against `world`, returning whatever it wrote to stderr.
static std::string resolveCapturingStderr(SelectionCriterion& criterion,
                                          const WorldState& world) {
  std::ostringstream captured;
  std::streambuf* previous = std::cerr.rdbuf(captured.rdbuf());
  criterion.resolveOrThrow(world, "test");
  std::cerr.rdbuf(previous);
  return captured.str();
}

static SelectionCriterion xlguCriterion(const std::string& op,
                                        PropertyValue value) {
  SelectionCriterion criterion;
  criterion.property_path = "geo_unit.XLGU";
  criterion.operator_type = op;
  criterion.value = std::move(value);
  return criterion;
}

TEST_CASE("activity venue criteria resolve global venue IDs") {
  WorldState world;
  world.activity_names = {"residence"};
  world.venue_type_names = {"household", "hospital"};

  Person& person = world.people.emplace_back();
  person.id = 7;
  person.age = 40;
  person.sex = Sex::FEMALE;
  person.activity_meta_start = 0;
  person.activity_meta_count = 1;
  world.activity_meta.push_back({0, 0, 1});

  // Venue 200 is a global venue that is not present in this rank's local
  // venue vector. This is the MPI ownership shape that exposed the bug.
  world.activity_venues.push_back({200, 0});
  world.buildIndices();
  world.setGlobalVenueType(200, 1);

  SelectionCriterion criterion;
  criterion.property_path = "activities.residence.venue_type";
  criterion.operator_type = "==";
  criterion.value = std::string("hospital");
  criterion.resolveOrThrow(world, "test");

  CHECK(criterion.evaluate(world.people.front(), &world));
}

TEST_CASE("activity venue criteria select the matching outcome row globally") {
  WorldState world;
  world.activity_names = {"residence"};
  world.venue_type_names = {"household", "hospital"};

  Person& person = world.people.emplace_back();
  person.id = 7;
  person.activity_meta_start = 0;
  person.activity_meta_count = 1;
  world.activity_meta.push_back({0, 0, 1});
  world.activity_venues.push_back({200, 0});
  world.buildIndices();
  world.setGlobalVenueType(200, 1);

  SelectionCriterion hospital;
  hospital.property_path = "activities.residence.venue_type";
  hospital.operator_type = "==";
  hospital.value = std::string("hospital");

  OutcomeRow hospital_row;
  hospital_row.criteria = {hospital};
  hospital_row.probabilities = {{"symptomatic", 1.0}};

  OutcomeRow fallback_row;
  fallback_row.probabilities = {{"asymptomatic", 1.0}};

  OutcomeRates rates;
  rates.rows = {hospital_row, fallback_row};
  rates.resolve(world);

  CHECK(rates.getRate(person, &world, "symptomatic") == doctest::Approx(1.0));
}

TEST_CASE("activity venue string operators see global venue types") {
  WorldState world;
  world.activity_names = {"residence"};
  world.venue_type_names = {"household", "hospital"};

  Person& person = world.people.emplace_back();
  person.id = 7;
  person.activity_meta_start = 0;
  person.activity_meta_count = 1;
  world.activity_meta.push_back({0, 0, 1});
  world.activity_venues.push_back({200, 0});
  world.buildIndices();
  world.setGlobalVenueType(200, 1);

  SelectionCriterion criterion;
  criterion.property_path = "activities.residence.venue_type";
  criterion.operator_type = "contains";
  criterion.value = std::string("hosp");
  criterion.resolveOrThrow(world, "test");

  CHECK(criterion.evaluate(person, &world));
}

}  // namespace

TEST_CASE("activity venue criteria use global types for cross-rank venues") {
  WorldState world;
  world.activity_names = {"residence"};
  world.venue_type_names = {"household", "hospital"};

  Person& person = world.people.emplace_back();
  person.id = 7;
  person.age = 40;
  person.sex = Sex::FEMALE;
  person.activity_meta_start = 0;
  person.activity_meta_count = 1;
  world.activity_meta.push_back({0, 0, 1});
  // Venue 200 is deliberately absent from this rank's local venues.
  world.activity_venues.push_back({200, 0});
  world.buildIndices();
  world.setGlobalVenueType(200, 1);

  SelectionCriterion criterion;
  criterion.property_path = "activities.residence.venue_type";
  criterion.operator_type = "==";
  criterion.value = std::string("hospital");
  criterion.resolveOrThrow(world, "test");

  CHECK(criterion.evaluate(world.people.front(), &world));
}

TEST_CASE(
    "a unit name found at no level is an error unless absent units are "
    "allowed") {
  WorldState world = buildNationWorld();

  SUBCASE("not allowed: resolving throws") {
    SelectionCriterion criterion = xlguCriterion("==", std::string("Atlantis"));
    CHECK_THROWS_AS(criterion.resolveOrThrow(world, "test"),
                    std::runtime_error);
  }

  SUBCASE("allowed: recorded, and matches nobody") {
    SelectionCriterion criterion = xlguCriterion("==", std::string("Atlantis"));
    criterion.allow_absent_geo_units = true;
    CHECK_NOTHROW(criterion.resolveOrThrow(world, "test"));
    CHECK(criterion.absentGeoUnitNames() ==
          std::vector<std::string>{"Atlantis"});
    for (const Person& person : world.people) {
      CHECK_FALSE(criterion.evaluate(person, &world));
    }
  }

  SUBCASE("allowed: the names that do exist still select") {
    SelectionCriterion criterion =
        xlguCriterion("in", std::vector<std::string>{"Atlantis", "Wales"});
    criterion.allow_absent_geo_units = true;
    criterion.resolveOrThrow(world, "test");
    CHECK(criterion.absentGeoUnitNames() ==
          std::vector<std::string>{"Atlantis"});
    CHECK_FALSE(criterion.evaluate(world.people[0], &world));
    CHECK(criterion.evaluate(world.people[2], &world));
  }
}

TEST_CASE(
    "a unit name at a different level is an error even when absent units are "
    "allowed") {
  WorldState world = buildNationWorld();
  // S00000001 exists, but at SGU: that is a misspelt level, not a missing
  // place.
  SelectionCriterion criterion = xlguCriterion("==", std::string("S00000001"));
  criterion.allow_absent_geo_units = true;
  CHECK_THROWS_AS(criterion.resolveOrThrow(world, "test"), std::runtime_error);
}

TEST_CASE(
    "an outcome table reports rows naming units the world lacks, and those "
    "rows match nobody") {
  WorldState world = buildNationWorld();

  OutcomeRates rates;
  OutcomeRow elsewhere;
  elsewhere.criteria = {xlguCriterion("==", std::string("Atlantis"))};
  elsewhere.probabilities = {{"mild", 0.9}};
  OutcomeRow everyone;
  everyone.probabilities = {{"mild", 0.4}};
  rates.rows = {elsewhere, everyone};

  const std::vector<std::string> absent = rates.resolve(world, {}, {});
  REQUIRE(absent.size() == 1);
  CHECK(absent[0].find("row 0") != std::string::npos);
  CHECK(absent[0].find("Atlantis") != std::string::npos);

  // First match wins, so a row that matched anyone would give 0.9.
  CHECK(rates.getRate(world.people[0], &world, "mild") == doctest::Approx(0.4));
}

TEST_CASE("an outcome table row the world cannot answer is an error") {
  WorldState world = buildNationWorld();
  SelectionCriterion unknown;
  unknown.property_path = "properties.no_such_property";
  unknown.operator_type = "==";
  unknown.value = std::string("x");

  OutcomeRates rates;
  OutcomeRow row;
  row.criteria = {unknown};
  rates.rows = {row};
  CHECK_THROWS_AS(rates.resolve(world, {}, {}), std::runtime_error);
}

TEST_CASE(
    "an infection seed attribute filter the world cannot answer is an error") {
  WorldState world = buildNationWorld();
  SelectionCriterion unknown;
  unknown.property_path = "properties.no_such_property";
  unknown.operator_type = "==";
  unknown.value = std::string("x");

  InfectionSeedEvent seed;
  seed.name = "seed_with_bad_filter";
  seed.attribute_filters = {unknown};
  InfectionSeedConfig config;
  config.seeds = {seed};
  CHECK_THROWS_AS(config.resolve(world), std::runtime_error);
}

TEST_CASE("geo_unit.<LEVEL> == name selects people under that ancestor") {
  WorldState world = buildNationWorld();

  SelectionCriterion criterion;
  criterion.property_path = "geo_unit.XLGU";
  criterion.operator_type = "==";
  criterion.value = std::string("Scotland");
  criterion.resolveOrThrow(world, "test");

  CHECK(criterion.evaluate(world.people[0], &world));
  CHECK_FALSE(criterion.evaluate(world.people[1], &world));

  // Dense ids, and no world passed: same answers.
  CHECK(criterion.evaluate(world.people[0]));
  CHECK_FALSE(criterion.evaluate(world.people[1]));
}

TEST_CASE(
    "geo_unit.<LEVEL> in a list of names selects people under any of them") {
  WorldState world = buildNationWorld();

  SelectionCriterion criterion;
  criterion.property_path = "geo_unit.XLGU";
  criterion.operator_type = "in";
  criterion.value = std::vector<std::string>{"London", "Wales"};
  criterion.resolveOrThrow(world, "test");

  CHECK_FALSE(criterion.evaluate(world.people[0], &world));
  CHECK(criterion.evaluate(world.people[1], &world));
  CHECK(criterion.evaluate(world.people[2], &world));
}

TEST_CASE("geo_unit.<LEVEL> != excludes, and over a list means not-in") {
  WorldState world = buildNationWorld();

  SUBCASE("single name") {
    SelectionCriterion criterion;
    criterion.property_path = "geo_unit.XLGU";
    criterion.operator_type = "!=";
    criterion.value = std::string("Scotland");
    criterion.resolveOrThrow(world, "test");

    CHECK_FALSE(criterion.evaluate(world.people[0], &world));
    CHECK(criterion.evaluate(world.people[1], &world));
  }

  SUBCASE("list of names") {
    SelectionCriterion criterion;
    criterion.property_path = "geo_unit.XLGU";
    criterion.operator_type = "!=";
    criterion.value = std::vector<std::string>{"Scotland", "Wales"};
    criterion.resolveOrThrow(world, "test");

    CHECK_FALSE(criterion.evaluate(world.people[0], &world));
    CHECK(criterion.evaluate(world.people[1], &world));
    CHECK_FALSE(criterion.evaluate(world.people[2], &world));
  }
}

TEST_CASE(
    "geo_unit.<LEVEL> at the person's own level matches the unit itself") {
  WorldState world = buildNationWorld();

  SelectionCriterion criterion;
  criterion.property_path = "geo_unit.SGU";
  criterion.operator_type = "==";
  criterion.value = std::string("S00000001");
  criterion.resolveOrThrow(world, "test");

  CHECK(criterion.evaluate(world.people[0], &world));
  CHECK_FALSE(criterion.evaluate(world.people[1], &world));
}

TEST_CASE("a person with no ancestor at the level matches neither == nor !=") {
  WorldState world = buildNationWorld();
  // A unit sitting at the fine level with no parent: absent at XLGU.
  GeographicalUnit stateless;
  stateless.id = 13;
  stateless.name = "Z00000001";
  stateless.level_id = 1;
  stateless.parent_id = -1;
  world.geo_units.push_back(stateless);
  Person& person = world.people.emplace_back();
  person.id = 3;
  person.geo_unit_id = 13;
  world.buildIndices();

  SelectionCriterion equals;
  equals.property_path = "geo_unit.XLGU";
  equals.operator_type = "==";
  equals.value = std::string("Scotland");
  equals.resolveOrThrow(world, "test");

  SelectionCriterion differs;
  differs.property_path = "geo_unit.XLGU";
  differs.operator_type = "!=";
  differs.value = std::string("Scotland");
  differs.resolveOrThrow(world, "test");

  CHECK_FALSE(equals.evaluate(world.people[3], &world));
  CHECK_FALSE(differs.evaluate(world.people[3], &world));
}

TEST_CASE("the no-ancestor warning fires only for units people live in") {
  SUBCASE("a correctly specified world stays silent") {
    WorldState world = buildNationWorld();
    // The nations have no ancestor at SGU, but nobody lives in one directly.
    SelectionCriterion criterion;
    criterion.property_path = "geo_unit.SGU";
    criterion.operator_type = "==";
    criterion.value = std::string("S00000001");

    CHECK(resolveCapturingStderr(criterion, world).empty());
  }

  SUBCASE("a unit people live in with no ancestor at the level warns") {
    WorldState world = buildNationWorld();
    GeographicalUnit stateless;
    stateless.id = 13;
    stateless.name = "Z00000001";
    stateless.level_id = 1;
    stateless.parent_id = -1;
    world.geo_units.push_back(stateless);
    Person& person = world.people.emplace_back();
    person.id = 3;
    person.geo_unit_id = 13;
    world.buildIndices();

    SelectionCriterion criterion;
    criterion.property_path = "geo_unit.XLGU";
    criterion.operator_type = "==";
    criterion.value = std::string("Scotland");

    const std::string warning = resolveCapturingStderr(criterion, world);
    CHECK(warning.find("1 inhabited geographical units") != std::string::npos);
  }
}

TEST_CASE(
    "a geo_unit criterion this world cannot answer fails loudly at load") {
  WorldState world = buildNationWorld();

  auto criterionFor = [](const std::string& path, const std::string& op,
                         PropertyValue value) {
    SelectionCriterion criterion;
    criterion.property_path = path;
    criterion.operator_type = op;
    criterion.value = std::move(value);
    return criterion;
  };

  SUBCASE("unknown level name") {
    auto criterion =
        criterionFor("geo_unit.REGION", "==", std::string("Scotland"));
    CHECK_THROWS_AS(criterion.resolveOrThrow(world, "test"),
                    std::runtime_error);
  }

  SUBCASE("unknown unit name") {
    auto criterion =
        criterionFor("geo_unit.XLGU", "==", std::string("Scotlnad"));
    CHECK_THROWS_AS(criterion.resolveOrThrow(world, "test"),
                    std::runtime_error);
  }

  SUBCASE("name that exists only at another level") {
    auto criterion =
        criterionFor("geo_unit.XLGU", "==", std::string("S00000001"));
    CHECK_THROWS_AS(criterion.resolveOrThrow(world, "test"),
                    std::runtime_error);
  }

  SUBCASE("one bad name in an otherwise good list") {
    auto criterion = criterionFor(
        "geo_unit.XLGU", "in", std::vector<std::string>{"Scotland", "Mercia"});
    CHECK_THROWS_AS(criterion.resolveOrThrow(world, "test"),
                    std::runtime_error);
  }

  SUBCASE("ambiguous name at the level") {
    GeographicalUnit duplicate;
    duplicate.id = 3;
    duplicate.name = "Scotland";
    duplicate.level_id = 0;
    duplicate.parent_id = -1;
    world.geo_units.push_back(duplicate);
    world.buildIndices();

    auto criterion =
        criterionFor("geo_unit.XLGU", "==", std::string("Scotland"));
    CHECK_THROWS_AS(criterion.resolveOrThrow(world, "test"),
                    std::runtime_error);
  }

  SUBCASE("a raw unit id instead of a name") {
    auto criterion = criterionFor("geo_unit.XLGU", "==", int32_t(0));
    CHECK_THROWS_AS(criterion.resolveOrThrow(world, "test"),
                    std::runtime_error);
  }

  SUBCASE("an operator geography cannot answer") {
    auto criterion =
        criterionFor("geo_unit.XLGU", ">", std::string("Scotland"));
    CHECK_THROWS_AS(criterion.resolveOrThrow(world, "test"),
                    std::runtime_error);
  }
}

TEST_CASE("a sparse geo id space gives the same answers") {
  WorldState world;
  world.geo_level_names = {"XLGU", "SGU"};

  auto addUnit = [&world](GeoUnitId id, const std::string& name,
                          uint8_t level_id, GeoUnitId parent_id) {
    GeographicalUnit unit;
    unit.id = id;
    unit.name = name;
    unit.level_id = level_id;
    unit.parent_id = parent_id;
    world.geo_units.push_back(unit);
  };

  addUnit(1000000, "Scotland", 0, -1);
  addUnit(2000000, "London", 0, -1);
  addUnit(3000000, "S00000001", 1, 1000000);
  addUnit(4000000, "E00000001", 1, 2000000);

  Person& scot = world.people.emplace_back();
  scot.id = 0;
  scot.geo_unit_id = 3000000;
  Person& londoner = world.people.emplace_back();
  londoner.id = 1;
  londoner.geo_unit_id = 4000000;
  world.buildIndices();

  SelectionCriterion criterion;
  criterion.property_path = "geo_unit.XLGU";
  criterion.operator_type = "==";
  criterion.value = std::string("Scotland");
  criterion.resolveOrThrow(world, "test");

  CHECK(criterion.evaluate(world.people[0], &world));
  CHECK_FALSE(criterion.evaluate(world.people[1], &world));

  // Sparse ids force the fallback keying; the answer must not depend on it, nor
  // on the caller having a WorldState in hand (ActivityManager has none).
  CHECK(criterion.evaluate(world.people[0]));
  CHECK_FALSE(criterion.evaluate(world.people[1]));

  // A unit this world does not carry is no ancestor of anything.
  Person stranger;
  stranger.id = 2;
  stranger.geo_unit_id = 3500000;
  CHECK_FALSE(criterion.evaluate(stranger));
}

TEST_CASE("a policies.yaml applies_to accepts a list of unit names") {
  WorldState world = buildNationWorld();
  world.activity_names = {"residence", "primary_activity"};
  world.buildIndices();

  ScopedTestFiles files{"june_test_geo_policies"};
  const auto path =
      files.write("policies.yaml",
                  "policies:\n"
                  "  temporal_policies:\n"
                  "    - name: \"english_and_welsh_lockdown\"\n"
                  "      start_date: \"2020-03-23\"\n"
                  "      end_date: \"2020-05-15\"\n"
                  "      override_activities: [\"primary_activity\"]\n"
                  "      replacement: \"residence\"\n"
                  "      compliance_rate: 1.0\n"
                  "      applies_to:\n"
                  "        - property: \"geo_unit.XLGU\"\n"
                  "          operator: \"in\"\n"
                  "          value: [\"London\", \"Wales\"]\n");

  PolicyManager policy_manager(world);
  PolicyLoader::loadPolicies(policy_manager, path.string(), "2020-01-01");
  REQUIRE(policy_manager.getTemporalPolicyCount() == 1);

  const TemporalPolicy& policy = policy_manager.getTemporalPolicies()[0];
  CHECK_FALSE(policy.appliesTo(world.people[0], &world));
  CHECK(policy.appliesTo(world.people[1], &world));
  CHECK(policy.appliesTo(world.people[2], &world));
}

TEST_CASE("a schedules.yaml selection accepts a list of unit names") {
  YAML::Node node = YAML::Load(
      "- property: \"geo_unit.XLGU\"\n"
      "  operator: \"in\"\n"
      "  value: [\"London\", \"Wales\"]\n");

  std::vector<SelectionCriterion> criteria;
  config_detail::parseSelectionCriteria(node, criteria);

  REQUIRE(criteria.size() == 1);
  REQUIRE(std::holds_alternative<std::vector<std::string>>(criteria[0].value));
  CHECK(std::get<std::vector<std::string>>(criteria[0].value).size() == 2);
}

TEST_CASE("selection shape validation preserves absent and empty semantics") {
  std::vector<SelectionCriterion> criteria;

  const YAML::Node absent = YAML::Load("name: fallback")["selection"];
  CHECK_FALSE(absent.IsDefined());
  CHECK_NOTHROW(config_detail::parseSelectionCriteria(absent, criteria));

  const YAML::Node empty = YAML::Load("selection: []")["selection"];
  CHECK_NOTHROW(config_detail::parseSelectionCriteria(empty, criteria));
  CHECK(criteria.empty());

  for (const std::string& value : {"null", "everyone", "{property: age}"}) {
    const YAML::Node malformed = YAML::Load("selection: " + value)["selection"];
    CHECK_THROWS_WITH(
        config_detail::parseSelectionCriteria(malformed, criteria),
        doctest::Contains("selection must be a sequence"));
  }
}

TEST_CASE("every selection loader rejects a present malformed selection") {
  ScopedTestFiles files{"june_malformed_selection"};

  const auto schedules = files.write("schedules.yaml",
                                     "day_type_cycle: [weekday]\n"
                                     "schedule_types:\n"
                                     "  fallback:\n"
                                     "    selection: {}\n");
  CHECK_THROWS_WITH(ConfigLoader::loadSchedule(schedules.string()),
                    doctest::Contains("selection must be a sequence"));
  const auto valid_schedules = files.write("valid_schedules.yaml",
                                           "day_type_cycle: [weekday]\n"
                                           "schedule_types:\n"
                                           "  fallback:\n"
                                           "    selection: []\n");
  CHECK(ConfigLoader::loadSchedule(valid_schedules.string())
            .schedule_types.size() == 1);

  const auto preferences = files.write("preferences.yaml",
                                       "profiles:\n"
                                       "  - name: bad\n"
                                       "    selection: null\n");
  CHECK_THROWS_WITH(ConfigLoader::loadActivityPreferences(preferences.string()),
                    doctest::Contains("selection must be a sequence"));
  const auto valid_preferences = files.write("valid_preferences.yaml",
                                             "profiles:\n"
                                             "  - name: fallback\n"
                                             "    selection: []\n");
  CHECK(ConfigLoader::loadActivityPreferences(valid_preferences.string())
            .profiles.size() == 1);

  const auto vaccines = files.write("vaccines.yaml",
                                    "enabled: true\n"
                                    "campaigns:\n"
                                    "  bad:\n"
                                    "    vaccine_type: V\n"
                                    "    tiers:\n"
                                    "      - {name: t, dose: 0, selection: "
                                    "{property: age}}\n");
  CHECK_THROWS_WITH(ConfigLoader::loadVaccination(vaccines.string()),
                    doctest::Contains("selection must be a sequence"));
  const auto valid_vaccines = files.write("valid_vaccines.yaml",
                                          "enabled: true\n"
                                          "vaccines:\n"
                                          "  V:\n"
                                          "    doses:\n"
                                          "      - {days_to_effective: 14, "
                                          "days_to_waning: 100, "
                                          "days_to_finished: 200}\n"
                                          "campaigns:\n"
                                          "  fallback:\n"
                                          "    vaccine_type: V\n"
                                          "    tiers:\n"
                                          "      - {name: t, dose: 0, "
                                          "selection: []}\n");
  CHECK(
      ConfigLoader::loadVaccination(valid_vaccines.string()).campaigns.size() ==
      1);

  const auto policy_applies =
      files.write("policy_applies.yaml",
                  "policies:\n"
                  "  temporal_policies:\n"
                  "    - name: bad\n"
                  "      applies_to: {property: age}\n");
  WorldState policy_world = buildNationWorld();
  PolicyManager policy_manager(policy_world);
  CHECK_THROWS_WITH(PolicyLoader::loadPolicies(
                        policy_manager, policy_applies.string(), "2020-01-01"),
                    doctest::Contains("selection must be a sequence"));
  const auto valid_policy_applies = files.write("valid_policy_applies.yaml",
                                                "policies:\n"
                                                "  temporal_policies:\n"
                                                "    - name: everyone\n"
                                                "      applies_to: []\n");
  WorldState valid_policy_world = buildNationWorld();
  PolicyManager valid_policy_manager(valid_policy_world);
  CHECK_NOTHROW(PolicyLoader::loadPolicies(
      valid_policy_manager, valid_policy_applies.string(), "2020-01-01"));

  const auto policy_exempt = files.write("policy_exempt.yaml",
                                         "policies:\n"
                                         "  temporal_policies:\n"
                                         "    - name: bad\n"
                                         "      exempt:\n"
                                         "        - activity: leisure\n"
                                         "          selection: everyone\n");
  WorldState exempt_world = buildNationWorld();
  PolicyManager exempt_manager(exempt_world);
  CHECK_THROWS_WITH(PolicyLoader::loadPolicies(
                        exempt_manager, policy_exempt.string(), "2020-01-01"),
                    doctest::Contains("selection must be a sequence"));
  const auto valid_policy_exempt = files.write("valid_policy_exempt.yaml",
                                               "policies:\n"
                                               "  temporal_policies:\n"
                                               "    - name: everyone\n"
                                               "      exempt:\n"
                                               "        - activity: leisure\n"
                                               "          selection: []\n");
  WorldState valid_exempt_world = buildNationWorld();
  PolicyManager valid_exempt_manager(valid_exempt_world);
  CHECK_NOTHROW(PolicyLoader::loadPolicies(
      valid_exempt_manager, valid_policy_exempt.string(), "2020-01-01"));

  const auto follow = files.write("follow.yaml",
                                  "coordinated_encounters:\n"
                                  "  encounters: []\n"
                                  "  follows:\n"
                                  "    - name: bad\n"
                                  "      eligibility: {property: age}\n");
  CHECK_THROWS_WITH(ConfigLoader::loadCoordinatedEncounters(follow.string()),
                    doctest::Contains("follow.eligibility must be a list"));
  const auto valid_follow = files.write("valid_follow.yaml",
                                        "coordinated_encounters:\n"
                                        "  encounters: []\n"
                                        "  follows:\n"
                                        "    - name: everyone\n"
                                        "      eligibility: []\n");
  CHECK(ConfigLoader::loadCoordinatedEncounters(valid_follow.string())
            .follows.size() == 1);
}

TEST_CASE(
    "a temporal policy refuses to resolve a filter this world cannot answer") {
  WorldState world = buildNationWorld();
  world.activity_names = {"residence", "primary_activity"};
  world.buildIndices();

  TemporalPolicy policy;
  policy.name = "misspelt_filter";
  policy.action.replacement_activity = "residence";
  SelectionCriterion misspelt;
  misspelt.property_path = "geo_unit.XLGU";
  misspelt.operator_type = "==";
  misspelt.value = std::string("Scotlnad");
  policy.applies_to.push_back(misspelt);

  CHECK_THROWS_AS(policy.resolve(world), std::runtime_error);
}

TEST_CASE("a list value on a non-geographical property must be whole numbers") {
  // Regression: a `catch (...)` fallback used to reinterpret any list that
  // failed int conversion as a list of unit names. Only geo_unit.<LEVEL>
  // reads those, so every other property got a criterion that loaded clean
  // and then matched nobody — a config mistake reading as "nobody qualifies".

  SUBCASE("a list of strings is rejected") {
    YAML::Node node = YAML::Load(
        "- property: \"properties.region\"\n"
        "  operator: \"in\"\n"
        "  value: [\"North East\", \"Yorkshire\"]\n");

    std::vector<SelectionCriterion> criteria;
    CHECK_THROWS_AS(config_detail::parseSelectionCriteria(node, criteria),
                    std::runtime_error);
  }

  SUBCASE("a list of floats is rejected") {
    YAML::Node node = YAML::Load(
        "- property: \"age\"\n"
        "  operator: \"in\"\n"
        "  value: [1.5, 2.5]\n");

    std::vector<SelectionCriterion> criteria;
    CHECK_THROWS_AS(config_detail::parseSelectionCriteria(node, criteria),
                    std::runtime_error);
  }

  SUBCASE("a list of ints still parses") {
    YAML::Node node = YAML::Load(
        "- property: \"geo_unit_id\"\n"
        "  operator: \"in\"\n"
        "  value: [10, 11]\n");

    std::vector<SelectionCriterion> criteria;
    config_detail::parseSelectionCriteria(node, criteria);

    REQUIRE(criteria.size() == 1);
    CHECK(std::holds_alternative<std::vector<int32_t>>(criteria[0].value));
  }
}

TEST_CASE(
    "the policy loader rejects the same list values as the schedule loader") {
  WorldState world = buildNationWorld();
  world.activity_names = {"residence", "primary_activity"};
  world.buildIndices();

  ScopedTestFiles files{"june_test_bad_list_policies"};
  const auto path =
      files.write("policies.yaml",
                  "policies:\n"
                  "  temporal_policies:\n"
                  "    - name: \"regional_lockdown\"\n"
                  "      start_date: \"2020-03-23\"\n"
                  "      override_activities: [\"primary_activity\"]\n"
                  "      replacement: \"residence\"\n"
                  "      applies_to:\n"
                  "        - property: \"properties.region\"\n"
                  "          operator: \"in\"\n"
                  "          value: [\"North East\", \"Yorkshire\"]\n");

  PolicyManager policy_manager(world);
  CHECK_THROWS_AS(
      PolicyLoader::loadPolicies(policy_manager, path.string(), "2020-01-01"),
      std::runtime_error);
}
