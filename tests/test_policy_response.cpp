#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/epidemiology.h"
#include "epidemiology/policy.h"
#include "loaders/policy_loader.h"
#include "policy_test_helpers.h"
#include "test_utils.h"
#include "transmission_fixtures.h"
#include "utils/deterministic_rng.h"
#include "venue_scope_fixtures.h"

using namespace june;
using namespace june::policy_test;
namespace vs = june::venue_scope_fixture;

namespace {

constexpr int16_t kResidence = 0;
constexpr int16_t kPrimary = 1;
constexpr uint8_t kHouseholdType = 0;
constexpr uint8_t kSchoolType = 1;
constexpr uint8_t kClassroomType = 2;
constexpr uint64_t kSeed = 2024;

// Mild lasts the whole run, so a case follows one stage unless it moves people.
Disease makeDisease() {
  return makeMildDeadDisease("ResponseDisease", 100.0, 1000.0);
}

// Households of one adult and `kids` pupils each. Pupil j of household h is
// registered to one classroom of a single school, dealt round robin, so
// siblings can sit in different classrooms of the same school.
//
//   venues 0 .. H-1     households (no parent)
//   venue  H            the school (no parent)
//   venues H+1 .. H+C   classrooms (parent: the school)
struct Shape {
  int households = 0;
  int kids = 0;
  int classrooms = 0;

  VenueId school() const { return households; }
  PersonId adult(int h) const { return h * (kids + 1); }
  PersonId kid(int h, int j) const { return h * (kids + 1) + 1 + j; }
  VenueId classroomOf(int h, int j) const {
    return school() + 1 + (h * kids + j) % classrooms;
  }
};

WorldState makeSchoolWorld(const Shape& shape) {
  WorldState world;
  world.activity_names = {"residence", "primary_activity"};
  world.venue_type_names = {"household", "school", "classroom"};

  const auto add_venue = [&](VenueId id, uint8_t type, VenueId parent) {
    Venue v{};
    v.id = id;
    v.type_id = type;
    v.geo_unit_id = 0;
    v.parent_id = parent;
    world.venues.push_back(v);
  };
  for (int h = 0; h < shape.households; ++h) add_venue(h, kHouseholdType, -1);
  add_venue(shape.school(), kSchoolType, -1);
  for (int c = 0; c < shape.classrooms; ++c)
    add_venue(shape.school() + 1 + c, kClassroomType, shape.school());

  for (int h = 0; h < shape.households; ++h) {
    for (int slot = 0; slot <= shape.kids; ++slot) {
      Person& p = world.people.emplace_back();
      p.id = shape.adult(h) + slot;
      p.age = slot == 0 ? 40.0f : 10.0f;
      p.geo_unit_id = -1;
      p.activity_meta_start = static_cast<uint32_t>(world.activity_meta.size());
      const auto add = [&](int16_t activity, VenueId venue) {
        world.activity_meta.push_back(
            {activity, static_cast<uint32_t>(world.activity_venues.size()), 1});
        world.activity_venues.push_back({venue, 0});
        ++p.activity_meta_count;
      };
      add(kResidence, h);
      if (slot > 0) add(kPrimary, shape.classroomOf(h, slot - 1));
    }
  }
  world.buildIndices();
  return world;
}

const std::string kMetrics =
    "observed_metrics:\n"
    "  deaths: {source: death, scope: world}\n"
    "  infected: {source: infection, scope: world}\n"
    "  class_cases: {source: infection, scope: venue, activity: "
    "primary_activity}\n"
    "  school_cases: {source: infection, scope: parent_venue, activity: "
    "primary_activity}\n"
    "  class_deaths: {source: death, scope: venue, activity: "
    "primary_activity}\n";

std::string temporalFile(const std::string& policies) {
  return kMetrics + "policies:\n  temporal_policies:\n" + policies;
}

// One temporal policy that sends primary_activity home, with the given trigger
// and extra policy lines (already indented).
std::string keepHome(const std::string& name, const std::string& trigger,
                     const std::string& extra = "") {
  return "    - name: " + name +
         "\n      override_activities: [primary_activity]\n"
         "      replacement: residence\n" +
         extra + "      trigger: " + trigger + "\n";
}

std::string response(const std::string& terms,
                     const std::string& unit = "person",
                     const std::string& base = "0.0") {
  return "{response: {base: " + base + ", terms: [" + terms +
         "]}, decision_unit: " + unit + "}";
}

const std::string kHousehold = "{shared_venue: residence}";
const std::string kDeaths = "{metric: deaths, weight: 0.1}";

// The draw a unit is held to, written out from the specification: a pure
// function of the base seed, the unit key and the policy index. The temporal
// policy at position i draws with index i + kTemporalOffset.
constexpr uint32_t kTemporalOffset =
    PolicyManager::kTemporalComplianceIndexOffset;

bool expectedAvoid(uint64_t key, uint32_t draw_index, double level) {
  SplitMix64 rng(mix_seed(kSeed, key, draw_index, 0xA701D4ULL));
  return level > uniform01(rng);
}

// A world with the policies loaded, driven the way the simulator drives it:
// count, reduce, advance.
struct Responding {
  WorldState world;
  Disease disease = makeDisease();
  PolicyManager pm;

  Responding(WorldState w, const std::string& file_text)
      : world(std::move(w)), pm(world) {
    pm.setBaseSeed(kSeed);
    REQUIRE(loadPolicyText(pm, file_text) == "");
    pm.observedMetrics().resolve(world, disease);
    pm.requireTriggerMetrics();
    pm.observedMetrics().buildVenueParents(world);
    pm.resolveAll(disease);
    applyEveryPolicy(pm, world);
  }

  void endDay(int day, int world_deaths = 0, std::vector<PersonId> present = {},
              std::vector<PersonId> died = {}) {
    for (int i = 0; i < world_deaths; ++i)
      pm.observedMetrics().onDeath(world.people[0]);
    for (PersonId id : died) pm.observedMetrics().onDeath(world.people[id]);
    for (PersonId id : present)
      pm.observedMetrics().onPresent(world.people[id], 1);
    pm.observedMetrics().reduce();
    pm.advanceTriggers(day);
  }

  const PolicyTrigger& trigger(size_t i = 0) const {
    return *pm.getTemporalPolicies()[i].trigger;
  }

  bool avoids(PersonId who, VenueId at, double time) {
    return pm
        .getOverride(world.people[who], kPrimary, at, 0,
                     SlotVenueType::fromVenue(at), time, 0)
        .has_value();
  }

  bool suppressed(PersonId who, VenueId at, uint8_t type, double time) const {
    return pm.suppressesParticipation(world.people[who], kPrimary,
                                      SlotVenueType::known(type, at), time);
  }

  // Which households keep their first pupil home at `time`.
  std::vector<bool> households(const Shape& shape, double time) {
    std::vector<bool> out;
    for (int h = 0; h < shape.households; ++h)
      out.push_back(avoids(shape.kid(h, 0), shape.classroomOf(h, 0), time));
    return out;
  }
};

int countTrue(const std::vector<bool>& v) {
  return static_cast<int>(std::count(v.begin(), v.end(), true));
}

// True when every unit avoiding in `a` also avoids in `b`.
bool subset(const std::vector<bool>& a, const std::vector<bool>& b) {
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i] && !b[i]) return false;
  return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Loader
// ---------------------------------------------------------------------------

TEST_CASE("PolicyLoader - a response loads and the policy is in force") {
  WorldState world;
  PolicyManager pm(world);
  REQUIRE(loadPolicyText(
              pm, temporalFile(keepHome(
                      "keep_home",
                      response(kDeaths + ", {metric: infected, weight: 0.02}",
                               kHousehold, "0.05")))) == "");
  const auto& policy = pm.getTemporalPolicies()[0];
  REQUIRE(policy.trigger.has_value());
  const PolicyTrigger& t = *policy.trigger;
  CHECK(t.voluntary());
  CHECK(t.base == 0.05);
  REQUIRE(t.terms.size() == 2);
  CHECK(t.terms[0].metric == "deaths");
  CHECK(t.terms[1].weight == 0.02);
  CHECK(t.decision_unit.shared_activity == "residence");
  CHECK(t.scope == MetricScope::World);
  // Nothing holds it off: the date window alone decides.
  CHECK(policy.window.contains(0.0));
  // Before the first day boundary the level is the base.
  CHECK(t.levelAtWorld() == 0.05);
}

TEST_CASE("PolicyLoader - the base is optional and means zero") {
  WorldState world;
  PolicyManager pm(world);
  REQUIRE(loadPolicyText(
              pm, temporalFile(keepHome(
                      "keep_home", "{response: {terms: [" + kDeaths +
                                       "]}, decision_unit: person}"))) == "");
  CHECK(pm.getTemporalPolicies()[0].trigger->base == 0.0);
}

TEST_CASE("PolicyLoader - every bad response is refused with a reason") {
  struct Bad {
    std::string trigger;
    const char* expected;
  };
  const std::string terms = "terms: [" + kDeaths + "]";
  const std::vector<Bad> cases = {
      {"{activate_at: 5, calm_days: 3, response: {" + terms +
           "}, decision_unit: person}",
       "never both; found 'activate_at' next to 'response'"},
      {"{metric: deaths, response: {" + terms + "}, decision_unit: person}",
       "never both; found 'metric' next to 'response'"},
      {"{response: {" + terms + "}}", "needs 'decision_unit'"},
      {"{response: {" + terms + "}, decision_unit: household}",
       "'decision_unit' must be 'person'"},
      {"{response: {" + terms + "}, decision_unit: {other: residence}}",
       "'decision_unit' must be 'person'"},
      {"{response: {" + terms +
           "}, decision_unit: {shared_venue: residence, extra: 1}}",
       "'decision_unit' must be 'person'"},
      {"{response: {" + terms + "}, decision_unit: {shared_venue: ''}}",
       "'shared_venue' needs an activity name"},
      {"{response: {base: 1.5, " + terms + "}, decision_unit: person}",
       "'base' must be from 0 to 1"},
      {"{response: {base: -0.1, " + terms + "}, decision_unit: person}",
       "'base' must be from 0 to 1"},
      {"{response: {base: lots, " + terms + "}, decision_unit: person}",
       "'base' must be a finite number"},
      {"{response: {terms: [{metric: deaths, weight: 0}]}, decision_unit: "
       "person}",
       "must be above 0"},
      {"{response: {terms: [{metric: deaths, weight: -0.2}]}, decision_unit: "
       "person}",
       "must be above 0"},
      {"{response: {terms: [{metric: deaths, weight: heavy}]}, decision_unit: "
       "person}",
       "'weight' must be a finite number"},
      {"{response: {terms: [{metric: nope, weight: 0.1}]}, decision_unit: "
       "person}",
       "unknown metric 'nope'"},
      {"{response: {terms: []}, decision_unit: person}",
       "non empty 'terms' list"},
      {"{response: {base: 0.1}, decision_unit: person}",
       "non empty 'terms' list"},
      {"{response: {terms: [" + kDeaths + ", " + kDeaths +
           "]}, decision_unit: person}",
       "appears in two terms"},
      {"{response: {terms: [{metric: deaths}]}, decision_unit: person}",
       "needs both 'metric' and 'weight'"},
      {"{response: {terms: [{metric: deaths, weight: 0.1, colour: red}]}, "
       "decision_unit: person}",
       "unknown field 'colour' in a term"},
      {"{response: {" + terms + ", colour: red}, decision_unit: person}",
       "unknown field 'colour' in 'response'"},
      {"{response: [" + kDeaths + "], decision_unit: person}",
       "'response' must be a mapping"},
      // A world term next to a venue term is a known limitation.
      {"{response: {terms: [" + kDeaths +
           ", {metric: class_cases, weight: 0.1}]}, decision_unit: person}",
       "differ in scope or activity"},
      {"{response: {terms: [{metric: class_cases, weight: 0.1}, "
       "{metric: school_cases, weight: 0.1}]}, decision_unit: person}",
       "differ in scope or activity"},
      {"{metric: deaths, activate_at: 5, calm_days: 3, decision_unit: person}",
       "'decision_unit' belongs to the 'response' form"},
  };
  for (const auto& c : cases) {
    WorldState world;
    PolicyManager pm(world);
    const std::string msg =
        loadPolicyText(pm, temporalFile(keepHome("keep_home", c.trigger)));
    INFO(c.trigger);
    CHECK(has(msg, c.expected));
    CHECK(has(msg, "temporal_policy 'keep_home' trigger: "));
  }
}

TEST_CASE("PolicyLoader - a policy with a response refuses what it cannot do") {
  const std::string trigger = response(kDeaths);
  const std::string symptom_head =
      kMetrics +
      "policies:\n  symptom_policies:\n    - name: stay_home\n"
      "      symptoms: [mild]\n      override_activities: [primary_activity]\n";
  struct Bad {
    std::string text;
    const char* expected;
  };
  const std::vector<Bad> cases = {
      // The key must be absent: even a value of 1 is refused.
      {temporalFile(
           keepHome("keep_home", trigger, "      compliance_rate: 1.0\n")),
       "'compliance_rate' must be absent"},
      {temporalFile(
           keepHome("keep_home", trigger, "      compliance_rate: 0.5\n")),
       "even a value of 1"},
      {temporalFile(keepHome("keep_home", trigger,
                             "      transmission_effects_file: masks.csv\n")),
       "transmission effects on a policy with a 'response' are not built yet"},
      {symptom_head + "      trigger: " +
           response("{metric: class_cases, weight: 0.1}") + "\n",
       "is counted per venue, which is not built yet for a symptom policy"},
      {symptom_head + "      trigger: " + response(kDeaths, kHousehold) + "\n",
       "'shared_venue' is not built yet for a symptom policy"},
      {symptom_head + "      replacement_schedule: travel\n      trigger: " +
           trigger + "\n",
       "'replacement_schedule' is not built yet on a policy with a 'response'"},
      {symptom_head +
           "      follow_up_policy: other\n      trigger: " + trigger + "\n",
       "'follow_up_policy' is not built yet on a policy with a 'response'"},
      {symptom_head + "      compliance_rate: 1\n      trigger: " + trigger +
           "\n",
       "'compliance_rate' must be absent"},
  };
  for (const auto& c : cases) {
    WorldState world;
    PolicyManager pm(world);
    INFO(c.text);
    CHECK(has(loadPolicyText(pm, c.text), c.expected));
  }
}

TEST_CASE(
    "PolicyLoader - a venue response must be allowed to act on its activity") {
  const std::string venue_trigger =
      response("{metric: class_cases, weight: 0.1}");
  const auto text = [&](const std::string& activity_line) {
    return kMetrics + "policies:\n  temporal_policies:\n" +
           "    - name: keep_home\n" + activity_line +
           "      replacement: residence\n      trigger: " + venue_trigger +
           "\n";
  };
  for (const char* lines : {"      override_activities: [leisure]\n",
                            "      override_activities: [residence]\n"}) {
    WorldState world;
    PolicyManager pm(world);
    INFO(lines);
    CHECK(has(loadPolicyText(pm, text(lines)),
              "'override_activities' must include 'primary_activity'"));
  }
  for (const char* lines :
       {"      override_activities: [primary_activity, leisure]\n",
        "      override_activities: ['*']\n", ""}) {
    WorldState world;
    PolicyManager pm(world);
    INFO(lines);
    CHECK(loadPolicyText(pm, text(lines)) == "");
  }
}

TEST_CASE("PolicyManager - resolving refuses what only the world can tell") {
  Disease disease = makeDisease();
  WorldState world = makeSchoolWorld({4, 1, 2});
  PolicyManager pm(world);

  SUBCASE("a decision unit naming an unknown activity") {
    REQUIRE(loadPolicyText(
                pm, temporalFile(keepHome(
                        "keep_home",
                        response(kDeaths, "{shared_venue: nowhere}")))) == "");
    pm.observedMetrics().resolve(world, disease);
    pm.requireTriggerMetrics();
    CHECK(has(thrown([&] { pm.resolveAll(disease); }),
              "decision_unit names the unknown activity 'nowhere'"));
  }
  SUBCASE("a voluntary policy named as a follow up") {
    const std::string text = kMetrics +
                             "policies:\n  symptom_policies:\n"
                             "    - name: first\n      symptoms: [mild]\n"
                             "      override_activities: [primary_activity]\n"
                             "      follow_up_policy: second\n"
                             "    - name: second\n      symptoms: [mild]\n"
                             "      override_activities: [primary_activity]\n"
                             "      trigger: " +
                             response(kDeaths) + "\n";
    REQUIRE(loadPolicyText(pm, text) == "");
    pm.observedMetrics().resolve(world, disease);
    pm.requireTriggerMetrics();
    CHECK(has(thrown([&] { pm.resolveAll(disease); }),
              "names the voluntary policy 'second' as its follow_up_policy"));
  }
}

// ---------------------------------------------------------------------------
// The level
// ---------------------------------------------------------------------------

TEST_CASE("Response - the level is the base plus the weighted counts") {
  const Shape shape{4, 1, 2};
  Responding r(
      makeSchoolWorld(shape),
      temporalFile(keepHome(
          "keep_home", response(kDeaths + ", {metric: infected, weight: 0.25}",
                                "person", "0.05"))));
  CHECK(r.trigger().levelAtWorld() == 0.05);  // before any day: the base

  // 2 deaths and 1 infected: 0.05 + 0.1 * 2 + 0.25 * 1, in term order.
  r.endDay(0, 2, {shape.kid(0, 0)});
  CHECK(r.trigger().levelAtWorld() == 0.05 + 0.1 * 2.0 + 0.25 * 1.0);

  // Deaths keep adding; the infected count is current and falls to 0.
  r.endDay(1, 1);
  CHECK(r.trigger().levelAtWorld() == 0.05 + 0.1 * 3.0 + 0.25 * 0.0);
}

TEST_CASE("Response - the level is clamped at 1 and never below 0") {
  const Shape shape{4, 1, 2};
  Responding r(
      makeSchoolWorld(shape),
      temporalFile(keepHome(
          "keep_home", response("{metric: deaths, weight: 0.4}", "person"))));
  CHECK(r.trigger().levelAtWorld() == 0.0);
  r.endDay(0, 2);
  CHECK(r.trigger().levelAtWorld() == 0.0 + 0.4 * 2.0);
  r.endDay(1, 1);  // 3 deaths: 1.2 before the clamp
  CHECK(r.trigger().levelAtWorld() == 1.0);
  r.endDay(2, 100);
  CHECK(r.trigger().levelAtWorld() == 1.0);

  // Level 1 sends everyone: every draw is below 1.
  for (int h = 0; h < shape.households; ++h)
    CHECK(r.avoids(shape.kid(h, 0), shape.classroomOf(h, 0), 3.0));
}

TEST_CASE("Response - day zero already uses the base") {
  const Shape shape{400, 1, 4};
  Responding r(
      makeSchoolWorld(shape),
      temporalFile(keepHome("keep_home", response(kDeaths, "person", "1.0"))));
  // Nothing has been reduced yet. With a base of 1 everyone is out from the
  // first slot, where a threshold trigger would still be off.
  CHECK(countTrue(r.households(shape, 0.0)) == shape.households);
}

// ---------------------------------------------------------------------------
// The draw
// ---------------------------------------------------------------------------

TEST_CASE(
    "Response - a person avoids exactly when the level beats their draw") {
  const Shape shape{300, 1, 4};
  Responding r(makeSchoolWorld(shape),
               temporalFile(keepHome("keep_home", response(kDeaths))));
  for (int day = 0; day < 11; ++day) {
    const double level = std::min(1.0, 0.0 + 0.1 * static_cast<double>(day));
    const double time = day + 1.0;
    int avoiding = 0;
    for (int h = 0; h < shape.households; ++h) {
      const PersonId kid = shape.kid(h, 0);
      const bool expected = expectedAvoid(kid, kTemporalOffset, level);
      CHECK(r.avoids(kid, shape.classroomOf(h, 0), time) == expected);
      avoiding += expected;
    }
    INFO("day " << day << " level " << level);
    if (day == 0) CHECK(avoiding == 0);
    if (day == 10) CHECK(avoiding == shape.households);
    r.endDay(day, 1);  // one death a day, so the next day's level is higher
  }
}

TEST_CASE("Response - a shared decision unit draws once per venue") {
  const Shape shape{300, 2, 4};
  Responding r(
      makeSchoolWorld(shape),
      temporalFile(keepHome("keep_home", response(kDeaths, kHousehold))));
  r.endDay(0, 4);
  const double level = 0.0 + 0.1 * 4.0;
  for (int h = 0; h < shape.households; ++h) {
    const bool expected =
        expectedAvoid(static_cast<uint64_t>(h), kTemporalOffset, level);
    for (int j = 0; j < shape.kids; ++j) {
      INFO("household " << h << " child " << j);
      CHECK(r.avoids(shape.kid(h, j), shape.classroomOf(h, j), 1.0) ==
            expected);
    }
  }
}

TEST_CASE(
    "Response - siblings agree when they share a unit and differ when not") {
  const Shape shape{300, 3, 5};
  const auto disagreeing = [&](const std::string& unit) {
    Responding r(makeSchoolWorld(shape),
                 temporalFile(keepHome("keep_home", response(kDeaths, unit))));
    r.endDay(0, 5);  // level 0.5
    int households_split = 0;
    for (int h = 0; h < shape.households; ++h) {
      bool any = false, all = true;
      for (int j = 0; j < shape.kids; ++j) {
        const bool a = r.avoids(shape.kid(h, j), shape.classroomOf(h, j), 1.0);
        any = any || a;
        all = all && a;
      }
      if (any && !all) ++households_split;
    }
    return households_split;
  };
  CHECK(disagreeing(kHousehold) == 0);
  // Three pupils each, level one half: a household splits with probability
  // 3/4, so about 225 of 300 do and the bound below is loose.
  CHECK(disagreeing("person") > 100);
}

TEST_CASE("Response - the share of units avoiding approaches the level") {
  const Shape shape{3000, 1, 8};
  Responding r(
      makeSchoolWorld(shape),
      temporalFile(keepHome("keep_home", response(kDeaths, kHousehold))));
  r.endDay(0, 3);
  const double level = 0.0 + 0.1 * 3.0;
  const int avoiding = countTrue(r.households(shape, 1.0));
  // Binomial(3000, 0.3): mean 900, standard deviation about 25. The bound is 5
  // standard deviations, so a correct draw fails it about once in 1.7 million.
  const double mean = shape.households * level;
  const double sd = std::sqrt(shape.households * level * (1.0 - level));
  INFO("avoiding " << avoiding << " expected " << mean << " sd " << sd);
  CHECK(std::abs(avoiding - mean) < 5.0 * sd);
}

TEST_CASE("Response - two policies draw independently") {
  const Shape shape{1000, 1, 8};
  const std::string trigger = response(kDeaths, kHousehold);
  Responding first(makeSchoolWorld(shape),
                   temporalFile(keepHome("keep_home", trigger)));
  // The same policy listed second, behind one that nobody complies with.
  Responding second(
      makeSchoolWorld(shape),
      temporalFile("    - name: nobody\n"
                   "      override_activities: [primary_activity]\n"
                   "      compliance_rate: 0\n" +
                   keepHome("keep_home", trigger)));
  first.endDay(0, 5);
  second.endDay(0, 5);
  const auto a = first.households(shape, 1.0);
  const auto b = second.households(shape, 1.0);
  for (int h = 0; h < shape.households; ++h)
    CHECK(b[static_cast<size_t>(h)] ==
          expectedAvoid(static_cast<uint64_t>(h), kTemporalOffset + 1, 0.5));

  int both = 0;
  for (size_t i = 0; i < a.size(); ++i) both += a[i] && b[i];
  // Independent halves overlap in a quarter: 250 of 1000, standard deviation
  // about 14. Identical draws would overlap in 500.
  CHECK(std::abs(both - 250) < 5 * 14);
  CHECK(a != b);
}

TEST_CASE("Response - an exempt person and a non exempt one share a unit") {
  const Shape shape{200, 2, 4};
  Responding r(makeSchoolWorld(shape),
               temporalFile(keepHome("keep_home", response(kDeaths, kHousehold),
                                     "      exempt:\n"
                                     "        - activity: primary_activity\n"
                                     "          selection:\n"
                                     "            - property: age\n"
                                     "              operator: \">=\"\n"
                                     "              value: 12\n")));
  for (int h = 0; h < shape.households; ++h)
    r.world.people[static_cast<size_t>(shape.kid(h, 1))].age = 12.0f;
  r.endDay(0, 5);

  int avoiding_units = 0;
  for (int h = 0; h < shape.households; ++h) {
    const bool unit_avoids =
        expectedAvoid(static_cast<uint64_t>(h), kTemporalOffset, 0.5);
    avoiding_units += unit_avoids;
    CHECK(r.avoids(shape.kid(h, 0), shape.classroomOf(h, 0), 1.0) ==
          unit_avoids);
    // The sibling aged 12 is exempt, whatever the unit decided.
    CHECK_FALSE(r.avoids(shape.kid(h, 1), shape.classroomOf(h, 1), 1.0));
    // And the exemption is asked of encounters and follows too.
    CHECK_FALSE(r.suppressed(shape.kid(h, 1), shape.classroomOf(h, 1),
                             kClassroomType, 1.0));
  }
  CHECK(avoiding_units > 50);
}

// ---------------------------------------------------------------------------
// Monotone, current, permanent
// ---------------------------------------------------------------------------

TEST_CASE("Response - raising a cumulative metric never frees a unit") {
  const Shape shape{500, 1, 8};
  Responding r(makeSchoolWorld(shape),
               temporalFile(keepHome(
                   "keep_home",
                   response("{metric: deaths, weight: 0.05}", kHousehold))));
  std::vector<bool> previous = r.households(shape, 0.5);
  CHECK(countTrue(previous) == 0);
  int last = 0;
  for (int day = 0; day < 25; ++day) {
    r.endDay(day, day % 3);  // 0, 1, 2 deaths in turn, never fewer
    const auto now = r.households(shape, day + 1.5);
    INFO("day " << day);
    CHECK(subset(previous, now));
    CHECK(countTrue(now) >= last);
    last = countTrue(now);
    previous = now;
  }
  CHECK(last > 0);
  // Deaths stop: the total stays, so avoidance is permanent.
  for (int day = 25; day < 40; ++day) {
    r.endDay(day);
    CHECK(r.households(shape, day + 1.5) == previous);
  }
}

TEST_CASE("Response - a current metric lifts avoidance when the count falls") {
  const Shape shape{500, 1, 8};
  Responding r(makeSchoolWorld(shape),
               temporalFile(keepHome(
                   "keep_home",
                   response("{metric: infected, weight: 0.1}", kHousehold))));
  const auto present = [&](int n) {
    std::vector<PersonId> out;
    for (int i = 0; i < n; ++i) out.push_back(shape.adult(i));
    return out;
  };
  r.endDay(0, 0, present(5));
  const auto high = r.households(shape, 1.5);
  r.endDay(1, 0, present(2));
  const auto lower = r.households(shape, 2.5);
  r.endDay(2, 0, present(0));
  const auto none = r.households(shape, 3.5);
  r.endDay(3, 0, present(5));
  const auto again = r.households(shape, 4.5);

  CHECK(countTrue(high) > 0);
  CHECK(subset(lower, high));
  CHECK(countTrue(lower) < countTrue(high));
  CHECK(countTrue(none) == 0);
  // Nothing is re-rolled: the same level brings back the same units.
  CHECK(again == high);
}

TEST_CASE("Response - a death count at a venue keeps that venue's level") {
  const Shape shape{200, 1, 4};
  Responding r(makeSchoolWorld(shape),
               temporalFile(keepHome(
                   "keep_home", response("{metric: class_deaths, weight: 0.5}",
                                         kHousehold))));
  // The first pupil of classroom 0 dies. The death counts once, at that
  // classroom, and the level there reaches only the pupils registered to it.
  std::vector<PersonId> classmates;
  for (int h = 0; h < shape.households; ++h)
    if (shape.classroomOf(h, 0) == shape.school() + 1)
      classmates.push_back(shape.kid(h, 0));
  REQUIRE(classmates.size() > 20);
  r.endDay(0, 0, {}, {classmates[0]});
  const VenueId first = shape.school() + 1;
  for (int day = 1; day < 20; ++day) r.endDay(day);

  for (int h = 0; h < shape.households; ++h) {
    const VenueId room = shape.classroomOf(h, 0);
    const bool expected =
        room == first &&
        expectedAvoid(static_cast<uint64_t>(h), kTemporalOffset, 0.5);
    CHECK(r.avoids(shape.kid(h, 0), room, 21.0) == expected);
  }
}

// ---------------------------------------------------------------------------
// Who is affected
// ---------------------------------------------------------------------------

namespace {

// The shared fixture world, with a response whose weight is 1: a count of one
// at a unit puts its level at 1, and no count leaves it at 0. Verdicts then
// do not depend on any draw.
std::string fixtureFile(const std::string& metric) {
  return temporalFile(
      keepHome("keep_home", response("{metric: " + metric + ", weight: 1.0}")));
}

}  // namespace

TEST_CASE("Response - only people the metric counts at a venue avoid it") {
  Responding r(vs::makeWorld(), fixtureFile("class_cases"));
  r.endDay(0, 0, {vs::kPupilA});  // one case in class A

  // A member of class A avoids it; class B has no case, so its members do not.
  CHECK(r.avoids(vs::kPupilA, vs::kClassA, 1.0));
  CHECK_FALSE(r.avoids(vs::kPupilB, vs::kClassB, 1.0));
  // Visitors and staff scheduled at class A are not members, so unlike a
  // closure the level does not reach them.
  CHECK_FALSE(r.avoids(vs::kPubWorker, vs::kClassA, 1.0));
  CHECK_FALSE(r.avoids(vs::kSchoolStaff, vs::kClassA, 1.0));
  CHECK_FALSE(r.avoids(vs::kOtherSchool, vs::kClassA, 1.0));
  // Pupil B is not registered to class A either.
  CHECK_FALSE(r.avoids(vs::kPupilB, vs::kClassA, 1.0));
  // A venue with no case is untouched, even one the pupil visits for leisure.
  CHECK_FALSE(r.avoids(vs::kPupilA, vs::kPub, 1.0));
}

TEST_CASE(
    "Response - a person with several venues responds to the scheduled one") {
  Responding r(vs::makeWorld(), fixtureFile("class_cases"));
  r.endDay(0, 0, {vs::kPupilA});  // a case in class A only

  // Registered to both classrooms: avoids A today, attends B.
  CHECK(r.avoids(vs::kTwoClasses, vs::kClassA, 1.0));
  CHECK_FALSE(r.avoids(vs::kTwoClasses, vs::kClassB, 1.0));
  // No venue at all (a traveller in transit) is never avoided.
  CHECK_FALSE(r.pm.getOverride(r.world.people[vs::kTwoClasses], kPrimary, -1, 0,
                               SlotVenueType::absent(), 1.0, 0)
                  .has_value());
}

TEST_CASE(
    "Response - a parent venue level reaches every child venue's members") {
  Responding r(vs::makeWorld(), fixtureFile("school_cases"));
  // One case, in class A. At parent venue scope it counts at the school.
  r.endDay(0, 0, {vs::kPupilA});

  CHECK(r.avoids(vs::kPupilA, vs::kClassA, 1.0));
  // A pupil of class B, which has no case, shares the school's level.
  CHECK(r.avoids(vs::kPupilB, vs::kClassB, 1.0));
  // A person registered to two classrooms counts once at the school.
  CHECK(r.avoids(vs::kTwoClasses, vs::kClassB, 1.0));
  // Staff registered to the school building itself are members of its unit.
  CHECK(r.avoids(vs::kSchoolStaff, vs::kSchool, 1.0));
  // The other school has no case, and a pub worker is not counted at a school.
  CHECK_FALSE(r.avoids(vs::kOtherSchool, vs::kClassC, 1.0));
  CHECK_FALSE(r.avoids(vs::kPubWorker, vs::kClassA, 1.0));
}

TEST_CASE("Response - the date window still applies") {
  const Shape shape{50, 1, 4};
  Responding r(makeSchoolWorld(shape),
               temporalFile(keepHome("keep_home",
                                     response("{metric: deaths, weight: 1.0}"),
                                     "      start_date: 2020-02-05\n"
                                     "      end_date: 2020-02-08\n")));
  r.endDay(0, 1);  // level 1 for everyone
  const VenueId room = shape.classroomOf(0, 0);
  CHECK_FALSE(r.avoids(shape.kid(0, 0), room, 1.5));
  CHECK(r.avoids(shape.kid(0, 0), room, 2.0));
  CHECK(r.avoids(shape.kid(0, 0), room, 4.9));
  CHECK_FALSE(r.avoids(shape.kid(0, 0), room, 5.0));
}

// ---------------------------------------------------------------------------
// No latch
// ---------------------------------------------------------------------------

TEST_CASE(
    "Response - a voluntary policy latches nothing and a plain one does") {
  const Shape shape{20, 1, 4};
  Responding voluntary(
      makeSchoolWorld(shape),
      temporalFile(keepHome("keep_home", response(kDeaths, kHousehold))));
  voluntary.endDay(0, 5);
  voluntary.households(shape, 1.0);
  for (const Person& p : voluntary.world.people) {
    CHECK(voluntary.pm.decidedPolicies(p, PolicyKind::Temporal).none());
    CHECK(voluntary.pm.participatingPolicies(p, PolicyKind::Temporal).none());
  }

  Responding plain(
      makeSchoolWorld(shape),
      temporalFile("    - name: stay_home\n"
                   "      override_activities: [primary_activity]\n"
                   "      replacement: residence\n"
                   "      compliance_rate: 1.0\n"));
  plain.households(shape, 1.0);
  for (int h = 0; h < shape.households; ++h) {
    const Person& p = plain.world.people[static_cast<size_t>(shape.kid(h, 0))];
    CHECK(plain.pm.decidedPolicies(p, PolicyKind::Temporal) == PolicyMask{0b1});
    CHECK(plain.pm.participatingPolicies(p, PolicyKind::Temporal) ==
          PolicyMask{0b1});
  }
}

TEST_CASE("Response - encounters and follows are told the same verdict") {
  const Shape shape{300, 1, 4};
  Responding r(
      makeSchoolWorld(shape),
      temporalFile(keepHome("keep_home", response(kDeaths, kHousehold))));
  r.endDay(0, 5);
  int avoiding = 0;
  for (int h = 0; h < shape.households; ++h) {
    const PersonId kid = shape.kid(h, 0);
    const VenueId room = shape.classroomOf(h, 0);
    // Asked first: the question must not decide anything for the schedule.
    const bool asked = r.suppressed(kid, room, kClassroomType, 1.0);
    CHECK(r.pm.decidedPolicies(r.world.people[static_cast<size_t>(kid)],
                               PolicyKind::Temporal)
              .none());
    CHECK(asked == r.avoids(kid, room, 1.0));
    avoiding += asked;
  }
  CHECK(avoiding > 100);
  CHECK(avoiding < 200);
}

// ---------------------------------------------------------------------------
// Symptom policies
// ---------------------------------------------------------------------------

TEST_CASE("Response - symptomatic people stay home more as deaths rise") {
  const Shape shape{400, 1, 4};
  const std::string text =
      kMetrics +
      "policies:\n  symptom_policies:\n"
      "    - name: sick_stay_home\n      symptoms: [mild]\n"
      "      override_activities: [primary_activity]\n"
      "      replacement: residence\n      trigger: " +
      response(kDeaths) + "\n";
  Responding r(makeSchoolWorld(shape), text);

  // Every other pupil is sick; the rest are well.
  for (int h = 0; h < shape.households; h += 2) {
    Person& p = r.world.people[static_cast<size_t>(shape.kid(h, 0))];
    p.infection = std::make_unique<Infection>(&r.disease, 0.0, &p, 42,
                                              kNoTransmissionContext, nullptr,
                                              "household", 0);
  }
  const auto sick_avoiding = [&](double time) {
    int sick = 0;
    for (int h = 0; h < shape.households; ++h) {
      const PersonId kid = shape.kid(h, 0);
      const bool out = r.avoids(kid, shape.classroomOf(h, 0), time);
      if (h % 2 == 0) {
        sick += out;
      } else {
        CHECK_FALSE(out);  // the well never avoid under a symptom policy
      }
    }
    return sick;
  };
  const auto draw = [&](double level) {
    int expected = 0;
    for (int h = 0; h < shape.households; h += 2)
      expected +=
          expectedAvoid(static_cast<uint64_t>(shape.kid(h, 0)), 0, level);
    return expected;
  };

  CHECK(sick_avoiding(1.0) == 0);  // no deaths yet: level 0
  int previous = 0;
  for (int day = 0; day < 10; ++day) {
    r.endDay(day, 1);
    const double level =
        std::min(1.0, 0.0 + 0.1 * static_cast<double>(day + 1));
    const int now = sick_avoiding(day + 1.5);
    INFO("day " << day << " level " << level);
    CHECK(now == draw(level));
    CHECK(now >= previous);
    previous = now;
  }
  CHECK(previous == shape.households / 2);  // level 1: every sick pupil
  // No decision was latched on the way.
  for (const Person& p : r.world.people)
    CHECK(r.pm.decidedPolicies(p, PolicyKind::Symptom).none());
}

// ---------------------------------------------------------------------------
// Startup check of the decision unit
// ---------------------------------------------------------------------------

TEST_CASE("Response - a good world passes the startup check") {
  const Shape shape{50, 2, 4};
  SUBCASE("a shared decision unit with a venue metric") {
    Responding r(makeSchoolWorld(shape),
                 temporalFile(keepHome(
                     "keep_home", response("{metric: class_cases, weight: 0.1}",
                                           kHousehold))));
    CHECK_NOTHROW(r.pm.checkDecisionUnits());
  }
  SUBCASE("a shared decision unit with a world metric") {
    Responding r(
        makeSchoolWorld(shape),
        temporalFile(keepHome("keep_home", response(kDeaths, kHousehold))));
    CHECK_NOTHROW(r.pm.checkDecisionUnits());
  }
}

TEST_CASE("Response - a decision unit of one person needs no startup check") {
  const Shape shape{50, 2, 4};
  const auto lose_residence = [&](Responding& r, PersonId who) {
    const Person& p = r.world.people[static_cast<size_t>(who)];
    r.world.activity_meta[p.activity_meta_start].venue_count = 0;
  };
  SUBCASE("a response whose members decide alone") {
    Responding r(makeSchoolWorld(shape),
                 temporalFile(keepHome("keep_home", response(kDeaths))));
    lose_residence(r, shape.adult(3));
    CHECK_NOTHROW(r.pm.checkDecisionUnits());
  }
  SUBCASE("a policy that is not voluntary") {
    Responding r(
        makeSchoolWorld(shape),
        temporalFile("    - name: plain\n"
                     "      override_activities: [primary_activity]\n"));
    lose_residence(r, shape.adult(3));
    CHECK_NOTHROW(r.pm.checkDecisionUnits());
  }
}

TEST_CASE(
    "Response - a person without exactly one venue fails the startup check") {
  const Shape shape{50, 2, 4};
  const auto residence_meta = [&](Responding& r, PersonId who) -> auto& {
    const Person& p = r.world.people[static_cast<size_t>(who)];
    return r.world.activity_meta[p.activity_meta_start];  // residence is first
  };

  SUBCASE("no residence, for a pupil the venue metric counts") {
    Responding r(makeSchoolWorld(shape),
                 temporalFile(keepHome(
                     "keep_home", response("{metric: class_cases, weight: 0.1}",
                                           kHousehold))));
    residence_meta(r, shape.kid(7, 1)).venue_count = 0;
    residence_meta(r, shape.kid(30, 0)).venue_count = 0;
    const std::string msg = thrown([&] { r.pm.checkDecisionUnits(); });
    // The lowest failing id is named, with the activity and the venue count.
    CHECK(has(msg, "person " + std::to_string(shape.kid(7, 1)) +
                       " has 0 venues for activity 'residence'"));
    CHECK(has(msg, "exactly 1 is required"));
    CHECK(has(msg, "2 of 100 people checked failed"));
  }
  SUBCASE("two residences") {
    Responding r(
        makeSchoolWorld(shape),
        temporalFile(keepHome("keep_home", response(kDeaths, kHousehold))));
    residence_meta(r, shape.adult(3)).venue_count = 2;
    const std::string msg = thrown([&] { r.pm.checkDecisionUnits(); });
    CHECK(has(msg, "person " + std::to_string(shape.adult(3)) +
                       " has 2 venues for activity 'residence'"));
    CHECK(has(msg, "1 of 150 people checked failed"));
  }
  SUBCASE("an adult the venue metric cannot count is not checked") {
    Responding r(makeSchoolWorld(shape),
                 temporalFile(keepHome(
                     "keep_home", response("{metric: class_cases, weight: 0.1}",
                                           kHousehold))));
    residence_meta(r, shape.adult(3)).venue_count =
        0;  // registered to no class
    CHECK(thrown([&] { r.pm.checkDecisionUnits(); }) == "");
  }
  SUBCASE("a person outside applies_to is not checked") {
    Responding r(
        makeSchoolWorld(shape),
        temporalFile(keepHome("keep_home", response(kDeaths, kHousehold),
                              "      applies_to:\n"
                              "        - property: age\n"
                              "          operator: \"<\"\n"
                              "          value: 18\n")));
    residence_meta(r, shape.adult(3)).venue_count = 0;  // an adult
    CHECK(thrown([&] { r.pm.checkDecisionUnits(); }) == "");
    residence_meta(r, shape.kid(3, 0)).venue_count = 0;  // a child
    CHECK(has(thrown([&] { r.pm.checkDecisionUnits(); }),
              "1 of 100 people checked failed"));
  }
}

// ---------------------------------------------------------------------------
// Checkpoint
// ---------------------------------------------------------------------------

TEST_CASE(
    "Response - a restored run reads the same levels and gives the same "
    "verdicts") {
  const Shape shape{300, 2, 6};
  const std::string text = temporalFile(
      keepHome("keep_home", response("{metric: class_cases, weight: 0.2}, "
                                     "{metric: class_deaths, weight: 0.05}",
                                     kHousehold, "0.1")));
  // A voluntary response needs terms of one scope, so both are per venue.
  Responding before(makeSchoolWorld(shape), text);
  const auto pupils = [&](int from, int to) {
    std::vector<PersonId> out;
    for (int h = from; h < to; ++h)
      for (int j = 0; j < shape.kids; ++j) out.push_back(shape.kid(h, j));
    return out;
  };
  before.endDay(0, 0, pupils(0, 60), {shape.kid(1, 0), shape.kid(2, 1)});
  before.endDay(1, 0, pupils(30, 90));
  REQUIRE(before.trigger().unit_levels.size() > 0);

  // What a checkpoint keeps: the metric values per venue, current ones too.
  const auto metrics = before.pm.observedMetrics().snapshot();
  bool current_has_units = false;
  for (const auto& m : metrics)
    if (m.name == "class_cases") current_has_units = !m.units.empty();
  CHECK(current_has_units);
  const auto states = before.pm.triggerStates();
  REQUIRE(states.size() == 1);
  CHECK(states[0].voluntary);

  Responding after(makeSchoolWorld(shape), text);
  after.pm.observedMetrics().restoreSnapshot(metrics);
  after.pm.restoreTriggerStates(states);
  after.pm.rebuildResponseLevels();

  // The first resumed day: the same level at every venue, bit for bit.
  CHECK(after.trigger().unit_levels == before.trigger().unit_levels);
  for (int h = 0; h < shape.households; ++h)
    for (int j = 0; j < shape.kids; ++j)
      CHECK(after.avoids(shape.kid(h, j), shape.classroomOf(h, j), 2.0) ==
            before.avoids(shape.kid(h, j), shape.classroomOf(h, j), 2.0));

  // And the two runs stay together.
  for (int day = 2; day < 6; ++day) {
    before.endDay(day, 0, pupils(day * 10, day * 10 + 40));
    after.endDay(day, 0, pupils(day * 10, day * 10 + 40));
    CHECK(after.trigger().unit_levels == before.trigger().unit_levels);
    CHECK(after.households(shape, day + 1.5) ==
          before.households(shape, day + 1.5));
  }
}

TEST_CASE("Response - restoring a trigger read in the other form is refused") {
  const Shape shape{10, 1, 2};
  const std::string threshold = temporalFile(
      keepHome("keep_home", "{metric: deaths, activate_at: 3, calm_days: 5}"));
  const std::string voluntary =
      temporalFile(keepHome("keep_home", response(kDeaths)));
  Responding t(makeSchoolWorld(shape), threshold);
  Responding v(makeSchoolWorld(shape), voluntary);
  CHECK(has(thrown([&] { v.pm.restoreTriggerStates(t.pm.triggerStates()); }),
            "is a threshold trigger in the checkpoint and a voluntary "
            "response trigger in this run"));
  CHECK(has(thrown([&] { t.pm.restoreTriggerStates(v.pm.triggerStates()); }),
            "is a voluntary response trigger in the checkpoint and a "
            "threshold trigger in this run"));
}
