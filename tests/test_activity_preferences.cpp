#include "core/config.h"
#include "doctest.h"
#include "test_utils.h"

using namespace june;

namespace {

SelectionCriterion ageAtLeast(int years) {
  SelectionCriterion c;
  c.property_path = "age";
  c.operator_type = ">=";
  c.value = years;
  return c;
}

// Two activities and three venue types, so a profile can be told apart by
// which activity and which venue type it speaks about.
WorldState makeWorld() {
  WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
  world.activity_names = {"leisure", "work"};
  world.venue_type_names = {"pub", "gym", "office"};
  world.people[0].age = 70.0f;
  world.people[1].age = 20.0f;
  return world;
}

PreferenceProfile makeProfile(const std::string& activity, int priority,
                              std::vector<SelectionCriterion> selection,
                              std::unordered_map<std::string, double> weights) {
  PreferenceProfile p;
  p.name = "profile";
  p.activity = activity;
  p.priority = priority;
  p.selection_criteria = std::move(selection);
  p.preference_weights = std::move(weights);
  return p;
}

}  // namespace

TEST_CASE("ActivityPreferences - a profile for an activity weights that activity") {
  WorldState world = makeWorld();
  ActivityPreferenceConfig prefs;
  prefs.profiles.push_back(
      makeProfile("leisure", 10, {ageAtLeast(65)}, {{"pub", 5.0}, {"gym", 0.25}}));
  prefs.resolve(world);

  const int leisure = world.getActivityIndex("leisure");
  const int pub = world.getVenueTypeIndex("pub");
  const int gym = world.getVenueTypeIndex("gym");
  const int office = world.getVenueTypeIndex("office");

  CHECK(prefs.getWeight(world.people[0], leisure, pub, &world) == 5.0);
  CHECK(prefs.getWeight(world.people[0], leisure, gym, &world) == 0.25);
  // A venue type the profile does not list keeps weight 1.
  CHECK(prefs.getWeight(world.people[0], leisure, office, &world) == 1.0);
  // The selection still gates membership.
  CHECK(prefs.getWeight(world.people[1], leisure, pub, &world) == 1.0);
}

TEST_CASE("ActivityPreferences - a profile does not weight other activities") {
  WorldState world = makeWorld();
  ActivityPreferenceConfig prefs;
  prefs.profiles.push_back(
      makeProfile("leisure", 10, {ageAtLeast(65)}, {{"pub", 5.0}}));
  prefs.resolve(world);

  const int work = world.getActivityIndex("work");
  const int pub = world.getVenueTypeIndex("pub");
  CHECK(prefs.getWeight(world.people[0], work, pub, &world) == 1.0);
}

TEST_CASE("ActivityPreferences - higher priority wins among matching profiles") {
  WorldState world = makeWorld();
  ActivityPreferenceConfig prefs;
  prefs.profiles.push_back(
      makeProfile("leisure", 20, {ageAtLeast(65)}, {{"pub", 7.0}}));
  prefs.profiles.push_back(
      makeProfile("leisure", 10, {ageAtLeast(18)}, {{"pub", 2.0}}));
  prefs.resolve(world);

  const int leisure = world.getActivityIndex("leisure");
  const int pub = world.getVenueTypeIndex("pub");
  CHECK(prefs.getWeight(world.people[0], leisure, pub, &world) == 7.0);
  CHECK(prefs.getWeight(world.people[1], leisure, pub, &world) == 2.0);
}

TEST_CASE("ActivityPreferences - an empty activity applies to every activity") {
  WorldState world = makeWorld();
  ActivityPreferenceConfig prefs;
  prefs.profiles.push_back(makeProfile("", 0, {}, {{"pub", 3.0}}));
  prefs.resolve(world);

  const int pub = world.getVenueTypeIndex("pub");
  for (const char* activity : {"leisure", "work"}) {
    CHECK(prefs.getWeight(world.people[0], world.getActivityIndex(activity), pub,
                          &world) == 3.0);
  }
}

TEST_CASE("ActivityPreferences - a profile naming an unknown activity is rejected") {
  WorldState world = makeWorld();
  ActivityPreferenceConfig prefs;
  prefs.profiles.push_back(makeProfile("liesure", 10, {}, {{"pub", 5.0}}));
  CHECK_THROWS_WITH_AS(prefs.resolve(world),
                       doctest::Contains("liesure"), std::runtime_error);
}

TEST_CASE("ActivityPreferences - a weight for an unknown venue type is rejected") {
  WorldState world = makeWorld();
  ActivityPreferenceConfig prefs;
  prefs.profiles.push_back(makeProfile("leisure", 10, {}, {{"cinemaa", 5.0}}));
  CHECK_THROWS_WITH_AS(prefs.resolve(world), doctest::Contains("cinemaa"),
                       std::runtime_error);
}
