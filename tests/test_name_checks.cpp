// Names a config writes must exist in the world. Each test sets one name
// wrong and expects a load time error that lists the valid names, then spells
// it correctly and expects a clean load. Worlds are synthetic.

#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/policy.h"

using namespace june;

namespace {

WorldState namedWorld() {
  WorldState world;
  world.activity_names = {"residence", "work", "leisure"};
  world.venue_type_names = {"home", "office", "pub"};
  world.network_names = {"friendships"};
  world.subset_type_names = {"worker", "patient", "resident"};
  world.person_property_names = {"ethnicity", "work_mode", "comorbidities"};
  world.buildIndices();
  return world;
}

ContactMatrix matrixWithBins(std::vector<std::string> bins) {
  ContactMatrix matrix;
  matrix.bins = std::move(bins);
  matrix.contacts.assign(matrix.bins.size(),
                         std::vector<double>(matrix.bins.size(), 1.0));
  return matrix;
}

CoordinatedEncounterConfig oneEncounter(const std::string& network,
                                        const std::string& trigger,
                                        const std::string& venue) {
  CoordinatedEncounterConfig config;
  config.enabled = true;
  CoordinatedEncounterDef def;
  def.name = "social_encounters";
  def.network = network;
  def.trigger_slots = {trigger};
  if (!venue.empty()) def.allowed_venues = {venue};
  config.encounters = {def};
  return config;
}

}  // namespace

TEST_CASE("a multi bin contact matrix names its bins in this world") {
  WorldState world = namedWorld();

  SUBCASE("a misspelt bin is refused with the world's subsets listed") {
    ContactMatrixConfig matrices;
    matrices.mode_matrices["hospital"]["respiratory"] =
        matrixWithBins({"worker", "patients"});
    CHECK_THROWS_WITH_AS(
        matrices.resolve(world),
        doctest::Contains("contact_matrices.yaml: venue type 'hospital', mode "
                          "'respiratory': bin 'patients' matches no subset "
                          "name in this world, no sex and no age band. "
                          "Subsets in this world: worker, patient, resident"),
        std::runtime_error);
  }

  SUBCASE("subset, sex and age band bins load") {
    ContactMatrixConfig matrices;
    matrices.mode_matrices["hospital"]["respiratory"] =
        matrixWithBins({"worker", "patient"});
    matrices.mode_matrices["pub"]["respiratory"] =
        matrixWithBins({"male", "female"});
    matrices.mode_matrices["office"]["respiratory"] =
        matrixWithBins({"0-17", "18-64", "65+"});
    matrices.matrices["home"] = matrixWithBins({"[0,17]", "[18,99]"});
    CHECK_NOTHROW(matrices.resolve(world));
    const ContactMatrix& office = matrices.mode_matrices["office"]["respiratory"];
    CHECK(office.age_to_bin[10] == 0);
    CHECK(office.age_to_bin[30] == 1);
    CHECK(office.age_to_bin[80] == 2);
    CHECK(matrices.mode_matrices["pub"]["respiratory"].male_bin == 0);
  }

  SUBCASE("a single bin matrix may carry any name") {
    ContactMatrixConfig matrices;
    matrices.mode_matrices["office"]["respiratory"] =
        matrixWithBins({"made_up_name"});
    CHECK_NOTHROW(matrices.resolve(world));
  }

  SUBCASE("the default matrix follows the same rule") {
    ContactMatrixConfig matrices;
    matrices.default_matrix = matrixWithBins({"worker", "typo"});
    CHECK_THROWS_WITH_AS(matrices.resolve(world),
                         doctest::Contains("the default matrix: bin 'typo'"),
                         std::runtime_error);
  }
}

TEST_CASE("a coordinated encounter names a network, activities and venues") {
  WorldState world = namedWorld();
  ContactMatrixConfig matrices;

  SUBCASE("a misspelt network") {
    auto config = oneEncounter("frendships", "leisure", "pub");
    CHECK_THROWS_WITH_AS(
        config.resolve(world, matrices),
        doctest::Contains("encounter 'social_encounters' network "
                          "'frendships' is not a network in this world. "
                          "Networks: friendships"),
        std::runtime_error);
  }
  SUBCASE("a misspelt trigger slot") {
    auto config = oneEncounter("friendships", "leisur", "pub");
    CHECK_THROWS_WITH_AS(
        config.resolve(world, matrices),
        doctest::Contains("trigger_slots names 'leisur', which is not an "
                          "activity in this world. Activities: residence, "
                          "work, leisure"),
        std::runtime_error);
  }
  SUBCASE("a misspelt allowed venue") {
    auto config = oneEncounter("friendships", "leisure", "pubb");
    CHECK_THROWS_WITH_AS(
        config.resolve(world, matrices),
        doctest::Contains("allowed_venues names 'pubb', which is not a venue "
                          "type in this world. Venue types: home, office, pub"),
        std::runtime_error);
  }
  SUBCASE("correct spelling loads, and no allowed venues is fine") {
    auto config = oneEncounter("friendships", "leisure", "pub");
    CHECK_NOTHROW(config.resolve(world, matrices));
    auto virtual_style = oneEncounter("friendships", "leisure", "");
    CHECK_NOTHROW(virtual_style.resolve(world, matrices));
  }
}

TEST_CASE("a policy names activities that exist") {
  WorldState world = namedWorld();

  SUBCASE("a misspelt override activity") {
    PolicyAction action;
    action.override_activities = {"leisur"};
    CHECK_THROWS_WITH_AS(
        action.resolve(world, "stay_in"),
        doctest::Contains("policy 'stay_in': override_activities names "
                          "'leisur', which is not an activity in this world. "
                          "Activities: residence, work, leisure"),
        std::runtime_error);
  }
  SUBCASE("a misspelt replacement") {
    PolicyAction action;
    action.override_activities = {"leisure"};
    action.replacement_activity = "residense";
    CHECK_THROWS_WITH_AS(
        action.resolve(world, "stay_in"),
        doctest::Contains("policy 'stay_in': replacement activity "
                          "'residense' is not an activity in this world"),
        std::runtime_error);
  }
  SUBCASE("correct names load, and star or an empty list mean all") {
    PolicyAction named;
    named.override_activities = {"leisure", "work"};
    named.replacement_activity = "residence";
    CHECK_NOTHROW(named.resolve(world, "p"));
    CHECK(named.shouldOverride(2));
    CHECK_FALSE(named.shouldOverride(0));

    PolicyAction star;
    star.override_activities = {"*"};
    star.replacement_activity = "residence";
    CHECK_NOTHROW(star.resolve(world, "p"));
    CHECK(star.override_all);

    PolicyAction empty;
    empty.replacement_activity = "residence";
    CHECK_NOTHROW(empty.resolve(world, "p"));
    CHECK(empty.override_all);
  }
}

TEST_CASE("output.summary_properties names properties the world carries") {
  WorldState world = namedWorld();
  SimulationConfig sim;
  sim.save_population_summary = true;
  sim.summary_properties_explicit = true;

  SUBCASE("a misspelt property") {
    sim.summary_properties = {"ethnicity", "work_mod"};
    CHECK_THROWS_WITH_AS(
        sim.resolve(world),
        doctest::Contains("output.summary_properties: 'work_mod' is not a "
                          "person property in this world. Properties: "
                          "ethnicity, work_mode, comorbidities"),
        std::runtime_error);
  }
  SUBCASE("more than four names are refused, not cut") {
    world.person_property_names = {"a", "b", "c", "d", "e"};
    world.buildIndices();
    sim.summary_properties = {"a", "b", "c", "d", "e"};
    CHECK_THROWS_WITH_AS(
        sim.resolve(world),
        doctest::Contains("lists 5 properties, at most 4 are stored"),
        std::runtime_error);
  }
  SUBCASE("correct names load and are all kept") {
    sim.summary_properties = {"ethnicity", "work_mode", "comorbidities"};
    CHECK_NOTHROW(sim.resolve(world));
    CHECK(sim.summaryPropertiesFor(world) == sim.summary_properties);
  }
  SUBCASE("the defaults are used only for properties the world carries") {
    SimulationConfig defaults;
    defaults.save_population_summary = true;
    CHECK_NOTHROW(defaults.resolve(world));
    CHECK(defaults.summaryPropertiesFor(world) ==
          std::vector<std::string>{"ethnicity", "work_mode"});
  }
}

TEST_CASE("partial presence prints one line for each type the world lacks") {
  WorldState world = namedWorld();
  SimulationConfig sim;
  sim.partial_presence.target_group_size_by_name = {{"pub", 10},
                                                    {"tube_line", 20}};

  std::ostringstream captured;
  std::streambuf* original = std::cout.rdbuf(captured.rdbuf());
  sim.resolve(world);
  std::cout.rdbuf(original);

  const std::string text = captured.str();
  CHECK(text.find("partial_presence: 'tube_line' is not in this world, "
                  "ignored") != std::string::npos);
  CHECK(text.find("'pub' is not in this world") == std::string::npos);
  CHECK(sim.partial_presence.target_group_size_by_type_id[2] == 10);
}

namespace {

ScheduleConfig oneSlotSchedule(const std::function<void(TimeSlot&)>& edit,
                               const std::function<void(ScheduleType&)>& tweak =
                                   [](ScheduleType&) {}) {
  ScheduleConfig config;
  ScheduleType sched;
  sched.name = "worker";
  TimeSlot slot;
  slot.name = "day";
  slot.allowed_activities = {"work", "leisure"};
  edit(slot);
  sched.flat_slots = {slot};
  tweak(sched);
  ScheduleType other;
  other.name = "away";
  config.schedule_types = {sched, other};
  return config;
}

}  // namespace

TEST_CASE("a schedule names activities, venue types and schedules that exist") {
  WorldState world = namedWorld();

  SUBCASE("hop_on_activity activity") {
    auto config = oneSlotSchedule(
        [](TimeSlot& s) { s.hop_on_activity = {{"wrk", "away"}}; });
    CHECK_THROWS_WITH_AS(config.resolveSlots(world),
                         doctest::Contains("hop_on_activity names 'wrk', which "
                                           "is not an activity in this world"),
                         std::runtime_error);
  }
  SUBCASE("hop_on_activity target schedule") {
    auto config = oneSlotSchedule(
        [](TimeSlot& s) { s.hop_on_activity = {{"work", "awya"}}; });
    CHECK_THROWS_WITH_AS(
        config.resolveSlots(world),
        doctest::Contains("targets schedule 'awya', which is not a schedule "
                          "type. Schedule types: worker, away"),
        std::runtime_error);
  }
  SUBCASE("property_hop_dispatch activity") {
    auto config = oneSlotSchedule([](TimeSlot& s) {
      s.property_hop_dispatch["wrk"] = {"work_mode", "mode_{value}"};
    });
    CHECK_THROWS_WITH_AS(
        config.resolveSlots(world),
        doctest::Contains("property_hop_dispatch names 'wrk', which is not an "
                          "activity in this world"),
        std::runtime_error);
  }
  SUBCASE("force_hybrid_activities") {
    auto config = oneSlotSchedule(
        [](TimeSlot&) {},
        [](ScheduleType& t) { t.force_hybrid_activities = {"leisur"}; });
    CHECK_THROWS_WITH_AS(
        config.resolveSlots(world),
        doctest::Contains("schedule type 'worker' force_hybrid_activities "
                          "names 'leisur', which is not an activity in this "
                          "world. Activities: residence, work, leisure"),
        std::runtime_error);
  }
  SUBCASE("linked_activities") {
    auto config = oneSlotSchedule(
        [](TimeSlot&) {},
        [](ScheduleType& t) { t.linked_activities = {"leisur"}; });
    CHECK_THROWS_WITH_AS(
        config.resolveSlots(world),
        doctest::Contains("schedule type 'worker' linked_activities names "
                          "'leisur'"),
        std::runtime_error);
  }
  SUBCASE("specified_activity type and venue type") {
    auto bad_type = oneSlotSchedule([](TimeSlot& s) {
      SpecifiedActivity spec;
      spec.type = "wrk";
      s.specified_activity = spec;
    });
    CHECK_THROWS_WITH_AS(
        bad_type.resolveSlots(world),
        doctest::Contains("specified_activity.type names 'wrk'"),
        std::runtime_error);

    auto bad_venue = oneSlotSchedule([](TimeSlot& s) {
      SpecifiedActivity spec;
      spec.type = "work";
      spec.venue_type = "ofice";
      s.specified_activity = spec;
    });
    CHECK_THROWS_WITH_AS(
        bad_venue.resolveSlots(world),
        doctest::Contains("specified_activity.venue_type names 'ofice', which "
                          "is not a venue type in this world. Venue types: "
                          "home, office, pub"),
        std::runtime_error);
  }
  SUBCASE("correct spelling loads") {
    auto config = oneSlotSchedule(
        [](TimeSlot& s) {
          s.hop_on_activity = {{"work", "away"}};
          s.property_hop_dispatch["work"] = {"work_mode", "mode_{value}"};
          SpecifiedActivity spec;
          spec.type = "work";
          spec.venue_type = "office";
          s.specified_activity = spec;
        },
        [](ScheduleType& t) {
          t.force_hybrid_activities = {"leisure"};
          t.linked_activities = {"work", "leisure"};
        });
    CHECK_NOTHROW(config.resolveSlots(world));
  }
}
