#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <string>

#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/policy.h"
#include "loaders/policy_loader.h"
#include "policy_test_helpers.h"
#include "test_utils.h"
#include "venue_scope_fixtures.h"

using namespace june;
using namespace june::policy_test;

namespace {

constexpr int16_t kResidence = 0;
constexpr int16_t kLeisure = 1;
constexpr VenueId kHome = 0;
constexpr VenueId kPub = 1;

// Mild ends after two days, so a case sees a person die.
Disease makeDisease() {
  return makeMildDeadDisease("TriggerDisease", 2.0, 100.0);
}

// One person whose residence is the home venue; a leisure closure sends them
// there.
WorldState buildWorld() {
  WorldState world = TestWorldFactory::createMinimalWorld(1, 2);
  world.activity_names = {"residence", "leisure"};
  world.venue_type_names = {"household", "pub"};
  world.venues[0].type_id = 0;
  world.venues[1].type_id = 1;
  world.buildIndices();
  for (auto& person : world.people) {
    person.activity_meta_start =
        static_cast<uint32_t>(world.activity_meta.size());
    person.activity_meta_count = 1;
    world.activity_meta.push_back(
        {kResidence, static_cast<uint32_t>(world.activity_venues.size()), 1});
    world.activity_venues.push_back({kHome, 0});
  }
  return world;
}

const std::string kMetrics =
    "observed_metrics:\n"
    "  deaths: {source: death, scope: world}\n"
    "  infected: {source: infection, scope: world}\n";

// A venue scope metric for tests that load a policy file against an empty
// world, where the activity does not resolve.
const std::string kVenueMetric =
    "  pub_cases: {source: infection, scope: venue, activity: leisure}\n";

// A policy file with one temporal policy carrying `trigger`. `window` is extra
// date lines for the policy, already indented.
std::string policyText(const std::string& trigger,
                       const std::string& window = "",
                       const std::string& extra_metrics = "") {
  return kMetrics + extra_metrics +
         "policies:\n"
         "  temporal_policies:\n"
         "    - name: close_leisure\n" +
         window +
         "      override_activities: [leisure]\n"
         "      replacement: residence\n"
         "      compliance_rate: 1.0\n"
         "      trigger: " +
         trigger + "\n";
}

PolicyTrigger cumulativeTrigger(int64_t activate_at, int64_t calm_days,
                                int64_t tolerance = 0,
                                int64_t tolerance_pct = 0) {
  PolicyTrigger t;
  t.cumulative = true;
  t.activate_at = activate_at;
  t.calm_days = calm_days;
  t.tolerance = tolerance;
  t.tolerance_pct = tolerance_pct;
  return t;
}

// The person a world scope death is counted for; a world total ignores them.
const Person nobody{};

// A loaded world with one triggered leisure closure, driven day by day the way
// the simulator drives it: count, reduce, advance.
struct TriggeredWorld {
  WorldState world = buildWorld();
  Disease disease = makeDisease();
  PolicyManager pm{world};

  explicit TriggeredWorld(const std::string& text) {
    REQUIRE(loadPolicyText(pm, text) == "");
    pm.observedMetrics().resolve(world, disease);
    pm.requireTriggerMetrics();
    pm.resolveAll(disease);
    applyTemporalPolicyZero(pm, world);
  }

  void endDay(int day, int deaths) {
    for (int i = 0; i < deaths; ++i) pm.observedMetrics().onDeath(nobody);
    pm.observedMetrics().reduce();
    pm.advanceTriggers(day);
  }

  bool overridden(double time) {
    return pm
        .getOverride(world.people[0], kLeisure, kPub, 0,
                     SlotVenueType::fromVenue(kPub), time, 0)
        .has_value();
  }
};

const std::string kDeathsTrigger =
    "{metric: deaths, activate_at: 3, calm_days: 5}";

}  // namespace

// ---------------------------------------------------------------------------
// State machine
// ---------------------------------------------------------------------------

TEST_CASE("PolicyTrigger - cumulative metric follows the 100 deaths example") {
  PolicyTrigger t = cumulativeTrigger(100, 200, 5);
  t.step(49, 99);
  CHECK_FALSE(t.closed);
  t.step(50, 100);
  CHECK(t.closed);
  t.step(60, 105);  // growth of 5 is within the tolerance
  CHECK(t.closed);
  t.step(249, 105);  // 199 days since closing
  CHECK(t.closed);
  t.step(250, 105);  // 200 days
  CHECK_FALSE(t.closed);
  CHECK(t.anchor_value == 105);  // the baseline is the value at reopening
  t.step(251, 204);
  CHECK_FALSE(t.closed);
  t.step(252, 205);  // 100 more than the baseline
  CHECK(t.closed);
}

TEST_CASE("PolicyTrigger - growth beyond the tolerance restarts the clock") {
  PolicyTrigger t = cumulativeTrigger(100, 200, 5);
  t.step(50, 100);
  t.step(100, 106);  // 6 more than 5 allowed: anchor moves to (100, 106)
  CHECK(t.anchor_day == 100);
  CHECK(t.anchor_value == 106);
  t.step(299, 106);
  CHECK(t.closed);
  t.step(300, 106);
  CHECK_FALSE(t.closed);

  // Growth counts from the anchor, not from the day before, so a slow creep
  // still moves it.
  PolicyTrigger creep = cumulativeTrigger(100, 200, 5);
  creep.step(50, 100);
  creep.step(60, 104);
  CHECK(creep.anchor_day == 50);
  creep.step(70, 108);
  CHECK(creep.anchor_day == 70);
}

TEST_CASE("PolicyTrigger - the percentage tolerance is of the anchor value") {
  // The first closure already has a non zero percentage tolerance: 5 percent
  // of the 100 it closed at, not 5 percent of the baseline of 0.
  PolicyTrigger t = cumulativeTrigger(100, 10, 0, 5);
  t.step(1, 100);
  CHECK(t.closed);
  t.step(5, 105);  // exactly the 5 allowed
  CHECK(t.anchor_day == 1);
  t.step(11, 105);
  CHECK_FALSE(t.closed);

  // A rise beyond it moves the anchor, and the next tolerance is 5 percent of
  // the new anchor value, rounded down.
  PolicyTrigger moved = cumulativeTrigger(100, 10, 0, 5);
  moved.step(1, 100);
  moved.step(5, 106);
  CHECK(moved.anchor_day == 5);
  CHECK(moved.anchor_value == 106);
  moved.step(14, 111);  // 5 more than 106, and 5 percent of 106 is 5
  CHECK(moved.anchor_day == 5);
  moved.step(15, 111);
  CHECK_FALSE(moved.closed);

  // 2 percent of an anchor of 400 deaths is a tolerance of 8.
  PolicyTrigger deaths = cumulativeTrigger(400, 30, 0, 2);
  deaths.step(10, 400);
  deaths.step(20, 408);
  CHECK(deaths.anchor_day == 10);
  deaths.step(21, 409);
  CHECK(deaths.anchor_day == 21);

  // The larger of the absolute and the percentage tolerance applies.
  PolicyTrigger both = cumulativeTrigger(400, 30, 10, 2);
  both.step(10, 400);
  both.step(20, 410);
  CHECK(both.anchor_day == 10);
}

TEST_CASE("PolicyTrigger - a current metric has hysteresis") {
  PolicyTrigger t;
  t.activate_at = 10;
  t.deactivate_at = 4;
  const int64_t levels[] = {9, 10, 7, 5, 4, 9, 10};
  const bool expected[] = {false, true, true, true, false, false, true};
  for (int day = 0; day < 7; ++day) {
    t.step(day, levels[day]);
    INFO("day " << day << " level " << levels[day]);
    CHECK(t.closed == expected[day]);
  }

  // A level of 0 is reachable, so a deactivate_at of 0 is a real threshold.
  PolicyTrigger zero;
  zero.activate_at = 10;
  zero.deactivate_at = 0;
  zero.step(0, 10);
  zero.step(1, 1);
  CHECK(zero.closed);
  zero.step(2, 0);
  CHECK_FALSE(zero.closed);
}

// ---------------------------------------------------------------------------
// Loader
// ---------------------------------------------------------------------------

TEST_CASE("PolicyLoader - a good trigger loads and the policy starts off") {
  WorldState world;
  SUBCASE("cumulative metric") {
    PolicyManager pm(world);
    REQUIRE(loadPolicyText(pm, policyText("{metric: deaths, activate_at: 100, "
                                          "calm_days: 200, tolerance: 5, "
                                          "tolerance_pct: 2}")) == "");
    const auto& trigger = pm.getTemporalPolicies()[0].trigger;
    REQUIRE(trigger.has_value());
    CHECK(trigger->cumulative);
    CHECK(trigger->activate_at == 100);
    CHECK(trigger->calm_days == 200);
    CHECK(trigger->tolerance == 5);
    CHECK(trigger->tolerance_pct == 2);
    CHECK_FALSE(pm.getTemporalPolicies()[0].window.contains(5.0));
  }
  SUBCASE("current metric, deactivate_at 0 allowed") {
    PolicyManager pm(world);
    REQUIRE(loadPolicyText(pm, policyText("{metric: infected, activate_at: 10, "
                                          "deactivate_at: 0}")) == "");
    const auto& trigger = pm.getTemporalPolicies()[0].trigger;
    REQUIRE(trigger.has_value());
    CHECK_FALSE(trigger->cumulative);
    CHECK(*trigger->deactivate_at == 0);
  }
  SUBCASE("symptom policy") {
    PolicyManager pm(world);
    const std::string text =
        kMetrics +
        "policies:\n"
        "  symptom_policies:\n"
        "    - name: stay_home\n"
        "      symptoms: [mild]\n"
        "      override_activities: [leisure]\n"
        "      trigger: {metric: deaths, activate_at: 3, calm_days: 5}\n";
    REQUIRE(loadPolicyText(pm, text) == "");
    REQUIRE(pm.getSymptomPolicies()[0].trigger.has_value());
    CHECK_FALSE(pm.getSymptomPolicies()[0].window.contains(5.0));
  }
  SUBCASE("a policy with no trigger is in force as before") {
    PolicyManager pm(world);
    REQUIRE(loadPolicyText(
                pm, kMetrics + "policies:\n  temporal_policies:\n"
                               "    - name: close_leisure\n"
                               "      override_activities: [leisure]\n") == "");
    CHECK_FALSE(pm.getTemporalPolicies()[0].trigger.has_value());
    CHECK(pm.getTemporalPolicies()[0].window.contains(5.0));
    pm.requireTriggerMetrics();
    CHECK_FALSE(pm.observedMetrics().anyActive());
  }
}

TEST_CASE("PolicyLoader - every bad trigger is refused with a reason") {
  struct Bad {
    const char* trigger;
    const char* expected;
  };
  const Bad cases[] = {
      {"{metric: nope, activate_at: 5, calm_days: 3}", "unknown metric 'nope'"},
      {"{activate_at: 5, calm_days: 3}", "'metric' is required"},
      {"{metric: deaths, activate_at: 5, calm_days: 3, deactivate_at: 1}",
       "'deactivate_at' is not valid for the cumulative metric 'deaths'"},
      {"{metric: infected, activate_at: 5, deactivate_at: 5}",
       "must be below 'activate_at'"},
      {"{metric: infected, activate_at: 5, deactivate_at: 8}",
       "must be below 'activate_at'"},
      {"{metric: infected, activate_at: 5}",
       "'deactivate_at' is required for the current metric 'infected'"},
      {"{metric: deaths, activate_at: 5}",
       "'calm_days' is required for the cumulative metric 'deaths'"},
      {"{metric: infected, activate_at: 5, deactivate_at: 1, calm_days: 3}",
       "only valid for a cumulative metric"},
      {"{metric: deaths, calm_days: 3}", "'activate_at' is required"},
      {"{metric: deaths, activate_at: 5, calm_days: 3, response: []}",
       "never both"},
      {"{metric: deaths, response: {}}",
       "never both; found 'metric' next to 'response'"},
      {"{metric: deaths}", "needs the threshold fields"},
      {"{metric: deaths, activate_at: 5, calm_days: 3, decision_unit: person}",
       "'decision_unit' belongs to the 'response' form"},
      {"{metric: deaths, activate_at: 0, calm_days: 3}",
       "'activate_at' must be at least 1"},
      {"{metric: deaths, activate_at: -3, calm_days: 3}",
       "'activate_at' must be at least 1"},
      {"{metric: deaths, activate_at: 5, calm_days: 0}",
       "'calm_days' must be at least 1"},
      {"{metric: deaths, activate_at: 5, calm_days: -1}",
       "'calm_days' must be at least 1"},
      {"{metric: infected, activate_at: 5, deactivate_at: -1}",
       "'deactivate_at' must not be negative"},
      {"{metric: deaths, activate_at: 5, calm_days: 3, tolerance: -1}",
       "'tolerance' must not be negative"},
      {"{metric: deaths, activate_at: 5, calm_days: 3, tolerance_pct: -1}",
       "whole percent from 0 to 100"},
      {"{metric: deaths, activate_at: 5, calm_days: 3, tolerance_pct: 101}",
       "whole percent from 0 to 100"},
      {"{metric: deaths, activate_at: 2.5, calm_days: 3}",
       "'activate_at' must be a whole number"},
      {"{metric: deaths, activate_at: ten, calm_days: 3}",
       "'activate_at' must be a whole number"},
      {"{metric: deaths, activate_at: 5, calm_days: 3, colour: red}",
       "unknown field 'colour'"},
  };
  WorldState world;
  for (const auto& c : cases) {
    PolicyManager pm(world);
    const std::string msg =
        loadPolicyText(pm, policyText(c.trigger, "", kVenueMetric));
    INFO(c.trigger);
    CHECK(has(msg, c.expected));
    CHECK(has(msg, "temporal_policy 'close_leisure' trigger: "));
  }
}

// ---------------------------------------------------------------------------
// In force
// ---------------------------------------------------------------------------

TEST_CASE("ActiveWindow - a trigger holds the policy off whatever the dates") {
  ActiveWindow window{0.0, 10.0};
  CHECK(window.contains(5.0));
  window.trigger_state = false;
  CHECK_FALSE(window.contains(5.0));
}

TEST_CASE("Trigger - the effect starts at the first slot of the next day") {
  TriggeredWorld w(policyText(kDeathsTrigger));
  CHECK_FALSE(w.overridden(0.5));  // starts off

  w.endDay(0, 2);
  CHECK_FALSE(w.overridden(1.0));

  // Day 1 is still running: nothing changes until its boundary.
  CHECK_FALSE(w.overridden(1.9));
  w.endDay(1, 1);  // the third death closes the trigger
  CHECK(w.overridden(2.0));
  CHECK(w.overridden(2.9));

  // Closed on day 1, so 5 calm days reopen it at the end of day 6.
  for (int day = 2; day <= 5; ++day) w.endDay(day, 0);
  CHECK(w.overridden(6.0));
  w.endDay(6, 0);
  CHECK_FALSE(w.overridden(7.0));
}

TEST_CASE("Trigger - the state evolves while the window is not open") {
  // The window opens on day 20 (2020-02-23 against a start of 2020-02-03).
  const std::string window = "      start_date: \"2020-02-23\"\n";

  SUBCASE("a threshold crossed before the window takes effect when it opens") {
    TriggeredWorld w(
        policyText("{metric: deaths, activate_at: 3, calm_days: 100}", window));
    w.endDay(0, 3);
    CHECK_FALSE(w.overridden(5.0));  // closed, but the window is not open
    for (int day = 1; day < 20; ++day) w.endDay(day, 0);
    CHECK_FALSE(w.overridden(19.9));
    CHECK(w.overridden(20.0));
  }
  SUBCASE("a trigger that reopened before the window does not hold it") {
    TriggeredWorld w(policyText(kDeathsTrigger, window));
    w.endDay(0, 3);
    for (int day = 1; day < 19; ++day) w.endDay(day, 0);  // reopens at day 5
    w.endDay(19, 0);
    CHECK_FALSE(w.overridden(20.0));
  }
}

TEST_CASE("Trigger - a policy without one ignores the daily step") {
  WorldState world = buildWorld();
  Disease disease = makeDisease();
  PolicyManager pm(world);
  TemporalPolicy plain;
  plain.name = "plain";
  plain.action.override_activities.insert("leisure");
  plain.action.replacement_activity = "residence";
  plain.action.compliance_rate = 1.0;
  pm.addTemporalPolicy(plain);
  pm.observedMetrics().resolve(world, disease);
  pm.requireTriggerMetrics();
  pm.resolveAll(disease);
  CHECK_FALSE(pm.observedMetrics().anyActive());
  pm.advanceTriggers(3);
  CHECK(pm.getTemporalPolicies()[0].window.contains(4.0));
  CHECK(pm.triggerStates().empty());
}

// ---------------------------------------------------------------------------
// Checkpoint
// ---------------------------------------------------------------------------

TEST_CASE("Trigger - closed state survives a round trip and resumes the same") {
  const std::string text = policyText(
      "{metric: deaths, activate_at: 3, calm_days: 5, tolerance: 1}");
  TriggeredWorld before(text);
  before.endDay(0, 2);
  before.endDay(1, 1);  // closes at day 1
  before.endDay(2, 1);  // growth of 1 is within the tolerance
  REQUIRE(before.overridden(3.0));

  const auto states = before.pm.triggerStates();
  REQUIRE(states.size() == 1);
  CHECK(states[0].key == "temporal/close_leisure");
  CHECK(states[0].cell.closed);
  CHECK(states[0].cell.anchor_day == 1);
  CHECK(states[0].cell.anchor_value == 3);

  // A fresh run restored the way the simulator does: metric totals, then
  // trigger state.
  TriggeredWorld after(text);
  after.pm.observedMetrics().restoreSnapshot(
      before.pm.observedMetrics().snapshot());
  after.pm.restoreTriggerStates(states);
  CHECK(after.overridden(3.0));

  // Both reopen on the same day, the end of day 6.
  for (int day = 3; day <= 5; ++day) {
    before.endDay(day, 0);
    after.endDay(day, 0);
    CHECK(before.overridden(day + 1.0));
    CHECK(after.overridden(day + 1.0));
  }
  before.endDay(6, 0);
  after.endDay(6, 0);
  CHECK_FALSE(before.overridden(7.0));
  CHECK_FALSE(after.overridden(7.0));
}

TEST_CASE("Trigger - restoring a different set of triggers is refused") {
  const std::string text = policyText(kDeathsTrigger);
  TriggeredWorld with(text);
  with.endDay(0, 3);
  const auto states = with.pm.triggerStates();

  SUBCASE("stored triggers, none in this run") {
    WorldState world;
    PolicyManager none(world);
    CHECK(has(thrown([&] { none.restoreTriggerStates(states); }),
              "cannot be rebuilt"));
  }
  SUBCASE("a trigger in this run, none stored") {
    CHECK(has(thrown([&] { with.pm.restoreTriggerStates({}); }),
              "cannot be rebuilt"));
  }
  SUBCASE("same count, different policy") {
    auto renamed = states;
    renamed[0].key = "temporal/other";
    CHECK(has(thrown([&] { with.pm.restoreTriggerStates(renamed); }),
              "temporal/other"));
  }
}

TEST_CASE("Trigger - reordered triggered policies restore by name") {
  const auto policy = [](const std::string& name, int activate_at) {
    return "    - name: " + name +
           "\n"
           "      override_activities: [leisure]\n"
           "      replacement: residence\n"
           "      compliance_rate: 1.0\n"
           "      trigger: {metric: deaths, activate_at: " +
           std::to_string(activate_at) + ", calm_days: 5}\n";
  };
  const std::string header = kMetrics + "policies:\n  temporal_policies:\n";
  // TriggeredWorld makes only the first temporal policy applicable.
  TriggeredWorld before(header + policy("closes", 3) + policy("stays_open", 100));
  applyEveryPolicy(before.pm, before.world);
  before.endDay(0, 3);
  const auto states = before.pm.triggerStates();
  REQUIRE(states.size() == 2);
  REQUIRE(states[0].key == "temporal/closes");
  REQUIRE(states[0].cell.closed);
  REQUIRE_FALSE(states[1].cell.closed);
  CHECK(before.overridden(1.0));

  TriggeredWorld after(header + policy("stays_open", 100) + policy("closes", 3));
  applyEveryPolicy(after.pm, after.world);
  after.pm.observedMetrics().restoreSnapshot(
      before.pm.observedMetrics().snapshot());
  after.pm.restoreTriggerStates(states);
  const auto restored = after.pm.triggerStates();
  REQUIRE(restored.size() == 2);
  CHECK(restored[0].key == "temporal/stays_open");
  CHECK_FALSE(restored[0].cell.closed);
  CHECK(restored[1].key == "temporal/closes");
  CHECK(restored[1].cell.closed);
  CHECK(restored[1].cell.anchor_day == states[0].cell.anchor_day);
  CHECK(after.overridden(1.0));
}

// ---------------------------------------------------------------------------
// Triggers counted per venue and per parent venue
// ---------------------------------------------------------------------------

namespace vs = june::venue_scope_fixture;

namespace {

// A policy file with one temporal policy that closes venues, triggered by a
// metric counted over primary_activity. `policy_lines` are extra policy fields,
// already indented.
std::string venuePolicyText(
    const std::string& scope,
    const std::string& trigger =
        "{metric: cases, activate_at: 2, "
        "deactivate_at: 0}",
    const std::string& policy_lines =
        "      override_activities: [primary_activity]\n",
    const std::string& source = "infection") {
  return "observed_metrics:\n"
         "  cases: {source: " +
         source + ", scope: " + scope +
         ", activity: primary_activity}\n"
         "policies:\n"
         "  temporal_policies:\n"
         "    - name: close_venue\n" +
         policy_lines +
         "      replacement: residence\n"
         "      trigger: " +
         trigger + "\n";
}

// The fixture world with one venue triggered Policy, driven the way the
// simulator drives it: count, reduce, advance.
struct VenueTriggered {
  WorldState world = vs::makeWorld();
  Disease disease = makeDisease();
  PolicyManager pm{world};

  explicit VenueTriggered(const std::string& text) {
    REQUIRE(loadPolicyText(pm, text) == "");
    pm.observedMetrics().resolve(world, disease);
    pm.requireTriggerMetrics();
    pm.observedMetrics().buildVenueParents(world);
    applyTemporalPolicyZero(pm, world);
    pm.resolveAll(disease);
  }

  void endDay(int day, std::vector<PersonId> present = {},
              std::vector<PersonId> died = {}) {
    for (PersonId id : present)
      pm.observedMetrics().onPresent(world.people[id], 1);
    for (PersonId id : died) pm.observedMetrics().onDeath(world.people[id]);
    pm.observedMetrics().reduce();
    pm.advanceTriggers(day);
  }

  bool closed(VenueId unit) const {
    return pm.getTemporalPolicies()[0].trigger->unitClosed(unit);
  }

  bool overridden(PersonId who, VenueId at, double time,
                  int16_t activity = vs::kPrimary) {
    return pm
        .getOverride(world.people[who], activity, at, 0,
                     SlotVenueType::fromVenue(at), time, 0)
        .has_value();
  }

  bool suppressed(PersonId who, VenueId at, uint8_t type, double time,
                  int16_t activity = vs::kPrimary) const {
    return pm.suppressesParticipation(world.people[who], activity,
                                      SlotVenueType::known(type, at), time);
  }
};

}  // namespace

TEST_CASE("VenueTrigger - two venues are in different states at once") {
  VenueTriggered t(venuePolicyText("venue"));

  // Class A has two people infected and class B one: only A reaches 2.
  t.endDay(0, {vs::kPupilA, vs::kTwoClasses});
  CHECK(t.closed(vs::kClassA));
  CHECK_FALSE(t.closed(vs::kClassB));
  CHECK(t.overridden(vs::kPupilA, vs::kClassA, 1.0));
  CHECK_FALSE(t.overridden(vs::kPupilB, vs::kClassB, 1.0));

  // Class B reaches 2 and closes. A is down to 1, which is above the reopening
  // level of 0, so it stays closed: the hysteresis is per venue.
  t.endDay(1, {vs::kPupilB, vs::kTwoClasses});
  CHECK(t.closed(vs::kClassA));
  CHECK(t.closed(vs::kClassB));

  // A is clear and reopens while B, at 1, stays closed. Opposite states.
  t.endDay(2, {vs::kPupilB});
  CHECK_FALSE(t.closed(vs::kClassA));
  CHECK(t.closed(vs::kClassB));
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 3.0));
  CHECK(t.overridden(vs::kPupilB, vs::kClassB, 3.0));
  // The other school's classroom was never counted.
  CHECK_FALSE(t.overridden(vs::kOtherSchool, vs::kClassC, 3.0));
}

TEST_CASE("VenueTrigger - a venue closes, reopens and closes again") {
  VenueTriggered t(venuePolicyText("venue"));
  const std::vector<PersonId> both = {vs::kPupilA, vs::kTwoClasses};

  t.endDay(0, both);  // 2: closes
  CHECK(t.overridden(vs::kPupilA, vs::kClassA, 1.0));
  t.endDay(1, {vs::kPupilA});  // 1: stays closed, hysteresis
  CHECK(t.overridden(vs::kPupilA, vs::kClassA, 2.0));
  t.endDay(2);  // 0: reopens
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 3.0));
  t.endDay(3, {vs::kPupilA});  // 1: below 2, still open
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 4.0));
  t.endDay(4, both);  // 2: closes a second time
  CHECK(t.overridden(vs::kPupilA, vs::kClassA, 5.0));
  t.endDay(5);
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 6.0));
}

TEST_CASE("VenueTrigger - a closure takes effect from the next day") {
  VenueTriggered t(venuePolicyText("venue"));
  for (PersonId id : {vs::kPupilA, vs::kTwoClasses})
    t.pm.observedMetrics().onPresent(t.world.people[id], 1);
  // Counted but not yet reduced and advanced: still open.
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 0.9));
  t.pm.observedMetrics().reduce();
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 0.9));
  t.pm.advanceTriggers(0);
  CHECK(t.overridden(vs::kPupilA, vs::kClassA, 1.0));
}

TEST_CASE("VenueTrigger - everyone scheduled at a closed venue is overridden") {
  VenueTriggered t(venuePolicyText("venue"));
  t.endDay(0, {vs::kPupilA, vs::kTwoClasses});

  // A visitor whose own registration is elsewhere, scheduled at class A today.
  CHECK(t.overridden(vs::kOtherSchool, vs::kClassA, 1.0));
  CHECK(t.overridden(vs::kPubWorker, vs::kClassA, 1.0));
  // Only the activity the policy names: a leisure visit is untouched.
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 1.0, vs::kLeisure));
  // The override sends people to their residence.
  const auto loc = t.pm.getOverride(
      t.world.people[vs::kOtherSchool], vs::kPrimary, vs::kClassA, 0,
      SlotVenueType::fromVenue(vs::kClassA), 1.0, 0);
  REQUIRE(loc.has_value());
  CHECK(loc->venue_id == vs::kHome);
  CHECK(loc->activity_index == vs::kResidence);
}

TEST_CASE(
    "VenueTrigger - an empty activity list closes the venue for every "
    "activity") {
  VenueTriggered t(venuePolicyText("venue",
                                   "{metric: cases, activate_at: 2, "
                                   "deactivate_at: 0}",
                                   ""));
  t.endDay(0, {vs::kPupilA, vs::kTwoClasses});
  CHECK(t.overridden(vs::kPupilA, vs::kClassA, 1.0, vs::kPrimary));
  CHECK(t.overridden(vs::kPupilA, vs::kClassA, 1.0, vs::kLeisure));
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassB, 1.0, vs::kLeisure));
}

TEST_CASE(
    "VenueTrigger - an exempt person stays, a non-exempt one goes, same venue "
    "same day") {
  VenueTriggered t(venuePolicyText(
      "venue", "{metric: cases, activate_at: 2, deactivate_at: 0}",
      "      override_activities: [primary_activity]\n"
      "      exempt:\n"
      "        - activity: primary_activity\n"
      "          selection:\n"
      "            - property: age\n"
      "              operator: \">=\"\n"
      "              value: 60\n"));
  t.world.people[vs::kTwoClasses].age = 61.0f;  // the exempt one
  t.endDay(0, {vs::kPupilA, vs::kTwoClasses});
  REQUIRE(t.closed(vs::kClassA));

  CHECK(t.overridden(vs::kPupilA, vs::kClassA, 1.0));
  CHECK_FALSE(t.overridden(vs::kTwoClasses, vs::kClassA, 1.0));
  // The same gates apply to the encounter and follow checks.
  CHECK(t.suppressed(vs::kPupilA, vs::kClassA, vs::kClassroomType, 1.0));
  CHECK_FALSE(
      t.suppressed(vs::kTwoClasses, vs::kClassA, vs::kClassroomType, 1.0));
}

TEST_CASE("VenueTrigger - the venue type gate still narrows what closes") {
  VenueTriggered t(venuePolicyText(
      "parent_venue", "{metric: cases, activate_at: 2, deactivate_at: 0}",
      "      override_activities: [primary_activity]\n"
      "      override_venue_types: [classroom]\n"));
  t.endDay(0, {vs::kPupilA, vs::kPupilB});
  REQUIRE(t.closed(vs::kSchool));

  // Classrooms of the closed school go, the school building itself is not a
  // classroom, so the gate lets its staff through.
  CHECK(t.overridden(vs::kPupilA, vs::kClassA, 1.0));
  CHECK(t.overridden(vs::kPupilB, vs::kClassB, 1.0));
  CHECK_FALSE(t.overridden(vs::kSchoolStaff, vs::kSchool, 1.0));
}

TEST_CASE(
    "VenueTrigger - closing a parent venue closes every child and the parent "
    "itself") {
  VenueTriggered t(venuePolicyText("parent_venue"));
  t.endDay(0, {vs::kPupilA, vs::kPupilB});  // two pupils of school 3

  CHECK(t.closed(vs::kSchool));
  // Children the infected pupils do not attend are closed too.
  CHECK(t.overridden(vs::kTwoClasses, vs::kClassA, 1.0));
  CHECK(t.overridden(vs::kTwoClasses, vs::kClassB, 1.0));
  // A person scheduled at the school itself.
  CHECK(t.overridden(vs::kSchoolStaff, vs::kSchool, 1.0));
  // The other school and an unrelated venue stay open.
  CHECK_FALSE(t.overridden(vs::kOtherSchool, vs::kClassC, 1.0));
  CHECK_FALSE(t.overridden(vs::kPubWorker, vs::kPub, 1.0));
  // The child is not itself the unit that closed.
  CHECK_FALSE(t.closed(vs::kClassA));
}

TEST_CASE(
    "VenueTrigger - one person in two children counts once at the parent") {
  VenueTriggered t(venuePolicyText("parent_venue"));
  // The person registered to both classrooms is the only one infected. Per
  // venue that is one at each classroom; per parent it is one, not two.
  t.endDay(0, {vs::kTwoClasses});
  CHECK_FALSE(t.closed(vs::kSchool));
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 1.0));
  // A second infected pupil takes the school to two.
  t.endDay(1, {vs::kTwoClasses, vs::kPupilB});
  CHECK(t.closed(vs::kSchool));
}

TEST_CASE("VenueTrigger - a venue with no parent is its own unit") {
  VenueTriggered t(venuePolicyText(
      "parent_venue", "{metric: cases, activate_at: 1, deactivate_at: 0}"));
  t.endDay(0, {vs::kPubWorker});  // registered to the pub for primary_activity
  CHECK(t.closed(vs::kPub));
  CHECK(t.overridden(vs::kPubWorker, vs::kPub, 1.0));
  CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 1.0));
  CHECK_FALSE(t.closed(vs::kSchool));
}

TEST_CASE(
    "VenueTrigger - a cumulative venue metric keeps its own baseline per "
    "venue") {
  VenueTriggered t(
      venuePolicyText("venue", "{metric: cases, activate_at: 1, calm_days: 2}",
                      "      "
                      "override_activities: [primary_activity]\n",
                      "death"));
  t.endDay(0, {}, {vs::kPupilA});  // class A: 1 death, closes
  CHECK(t.closed(vs::kClassA));
  CHECK_FALSE(t.closed(vs::kClassB));
  t.endDay(1);
  CHECK(t.closed(vs::kClassA));
  t.endDay(2, {}, {vs::kPupilB});  // A calm for 2 days: reopens; B closes
  CHECK_FALSE(t.closed(vs::kClassA));
  CHECK(t.closed(vs::kClassB));
  // A's baseline is now 1, so one more death is growth of 1 and closes it.
  t.endDay(3, {}, {vs::kTwoClasses});
  CHECK(t.closed(vs::kClassA));
  CHECK(t.closed(vs::kClassB));
}

TEST_CASE("VenueTrigger - the scheduled venue decides, not the pin venue") {
  VenueTriggered t(venuePolicyText("venue"));
  t.endDay(0, {vs::kPupilA, vs::kTwoClasses});
  Person& traveller = t.world.people[vs::kPupilB];

  // In transit: pinned at closed class A but occupying no venue this slot.
  auto in_transit = t.pm.getOverride(traveller, vs::kPrimary, vs::kClassA, 0,
                                     SlotVenueType::absent(), 1.0, 0);
  CHECK_FALSE(in_transit.has_value());
  // Scheduled at closed class A while pinned somewhere open.
  auto scheduled =
      t.pm.getOverride(traveller, vs::kPrimary, vs::kClassB, 0,
                       SlotVenueType::fromVenue(vs::kClassA), 1.0, 0);
  CHECK(scheduled.has_value());
}

TEST_CASE("VenueTrigger - encounters and follows are asked about the venue") {
  VenueTriggered t(venuePolicyText("venue"));
  t.endDay(0, {vs::kPupilA, vs::kTwoClasses});

  CHECK(t.suppressed(vs::kPupilB, vs::kClassA, vs::kClassroomType, 1.0));
  CHECK_FALSE(t.suppressed(vs::kPupilB, vs::kClassB, vs::kClassroomType, 1.0));
  // No venue named, or a virtual encounter: nothing to be closed.
  CHECK_FALSE(t.pm.suppressesParticipation(
      t.world.people[vs::kPupilB], vs::kPrimary,
      SlotVenueType::known(vs::kClassroomType), 1.0));
  CHECK_FALSE(t.pm.suppressesParticipation(
      t.world.people[vs::kPupilB], vs::kPrimary, SlotVenueType::absent(), 1.0));
}

TEST_CASE(
    "PolicyLoader - a trigger counted per venue loads and starts with every "
    "venue open") {
  for (const char* scope : {"venue", "parent_venue"}) {
    INFO(scope);
    VenueTriggered t(venuePolicyText(scope));
    const auto& trigger = t.pm.getTemporalPolicies()[0].trigger;
    REQUIRE(trigger.has_value());
    CHECK(trigger->perUnit());
    CHECK(trigger->closedUnits().empty());
    // In force from the start; the set of closed venues decides where it acts.
    CHECK(t.pm.getTemporalPolicies()[0].window.contains(0.5));
    CHECK_FALSE(t.overridden(vs::kPupilA, vs::kClassA, 0.5));
  }
}

TEST_CASE("PolicyLoader - every bad venue trigger is refused with a reason") {
  struct Bad {
    std::string text;
    const char* expected;
  };
  const std::string metric =
      "observed_metrics:\n"
      "  cases: {source: infection, scope: venue, activity: "
      "primary_activity}\n";
  const std::string trig = "{metric: cases, activate_at: 2, deactivate_at: 0}";
  const std::vector<Bad> cases = {
      {metric +
           "policies:\n  symptom_policies:\n    - name: stay_home\n"
           "      symptoms: [mild]\n      override_activities: "
           "[primary_activity]\n"
           "      trigger: " +
           trig + "\n",
       "is counted per venue, which is not built yet for a symptom policy"},
      {venuePolicyText("venue", trig,
                       "      override_activities: [primary_activity]\n"
                       "      compliance_rate: 0.7\n"),
       "'compliance_rate' must be 1 on a policy triggered per venue"},
      {venuePolicyText("venue", trig,
                       "      override_activities: [primary_activity]\n"
                       "      transmission_effects_file: effects.csv\n"),
       "transmission effects on a policy triggered per venue are not built "
       "yet"},
      {venuePolicyText("venue", trig, "      override_activities: [leisure]\n"),
       "'override_activities' must include 'primary_activity'"},
      {venuePolicyText("parent_venue", trig,
                       "      override_activities: [leisure, residence]\n"),
       "'override_activities' must include 'primary_activity'"},
  };
  for (const auto& c : cases) {
    WorldState world;
    PolicyManager pm(world);
    INFO(c.text);
    CHECK(has(loadPolicyText(pm, c.text), c.expected));
  }

  SUBCASE("a list that includes the activity, '*' and an empty list all load") {
    for (const char* lines :
         {"      override_activities: [primary_activity, leisure]\n",
          "      override_activities: ['*']\n", ""}) {
      WorldState world;
      PolicyManager pm(world);
      INFO(lines);
      CHECK(loadPolicyText(pm, venuePolicyText("venue", trig, lines)) == "");
    }
  }
}

TEST_CASE(
    "VenueTrigger - closed venues and a reopened baseline survive a round "
    "trip") {
  const std::string text = venuePolicyText(
      "venue", "{metric: cases, activate_at: 1, calm_days: 2}",
      "      override_activities: [primary_activity]\n", "death");
  VenueTriggered before(text);
  before.endDay(0, {}, {vs::kPupilA});  // A closes
  before.endDay(1);
  before.endDay(2, {}, {vs::kPupilB});  // A reopens with baseline 1, B closes
  REQUIRE_FALSE(before.closed(vs::kClassA));
  REQUIRE(before.closed(vs::kClassB));

  const auto states = before.pm.triggerStates();
  REQUIRE(states.size() == 1);
  CHECK(states[0].key == "temporal/close_venue");
  CHECK(states[0].scope == MetricScope::Venue);
  REQUIRE(states[0].units.size() == 2);
  CHECK(states[0].units[0].first == vs::kClassA);
  CHECK_FALSE(states[0].units[0].second.closed);
  CHECK(states[0].units[0].second.anchor_value == 1);  // the baseline
  CHECK(states[0].units[1].first == vs::kClassB);
  CHECK(states[0].units[1].second.closed);

  VenueTriggered after(text);
  after.pm.observedMetrics().restoreSnapshot(
      before.pm.observedMetrics().snapshot());
  after.pm.restoreTriggerStates(states);
  CHECK_FALSE(after.closed(vs::kClassA));
  CHECK(after.closed(vs::kClassB));
  CHECK(after.overridden(vs::kPupilB, vs::kClassB, 3.0));

  // From here the two runs agree every day: A closes again on a death (its
  // baseline survived), B's anchor moves with the death, both reopen together.
  const std::vector<std::vector<PersonId>> deaths = {
      {vs::kTwoClasses}, {}, {}, {}};
  for (size_t i = 0; i < deaths.size(); ++i) {
    const int day = 3 + static_cast<int>(i);
    before.endDay(day, {}, deaths[i]);
    after.endDay(day, {}, deaths[i]);
    for (VenueId v : {vs::kClassA, vs::kClassB}) {
      INFO("day " << day << " venue " << v);
      CHECK(before.closed(v) == after.closed(v));
    }
    CHECK(before.overridden(vs::kPupilA, vs::kClassA, day + 1.0) ==
          after.overridden(vs::kPupilA, vs::kClassA, day + 1.0));
  }
  // Day 3: both closed. Reopened together two calm days later, at day 5.
  CHECK_FALSE(before.closed(vs::kClassA));
  CHECK_FALSE(before.closed(vs::kClassB));
}

TEST_CASE("VenueTrigger - restoring under another scope is refused") {
  const std::string venue_text = venuePolicyText(
      "venue", "{metric: cases, activate_at: 1, calm_days: 2}",
      "      override_activities: [primary_activity]\n", "death");
  const std::string world_text =
      "observed_metrics:\n"
      "  cases: {source: death, scope: world}\n"
      "policies:\n  temporal_policies:\n    - name: close_venue\n"
      "      override_activities: [primary_activity]\n"
      "      replacement: residence\n"
      "      trigger: {metric: cases, activate_at: 1, calm_days: 2}\n";

  VenueTriggered per_venue(venue_text);
  per_venue.endDay(0, {}, {vs::kPupilA});
  VenueTriggered world_total(world_text);
  world_total.endDay(0, {}, {vs::kPupilA});

  CHECK(has(thrown([&] {
              world_total.pm.restoreTriggerStates(per_venue.pm.triggerStates());
            }),
            "another scope than the checkpoint's"));
  CHECK(has(thrown([&] {
              per_venue.pm.restoreTriggerStates(world_total.pm.triggerStates());
            }),
            "another scope than the checkpoint's"));
}
