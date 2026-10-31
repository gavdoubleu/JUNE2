
#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/calendar_event.h"
#include "test_utils.h"

using namespace june;

namespace {

CalendarEvent makeCatchmentEvent(int32_t id, int16_t sched_idx,
                                 float compliance = 1.0f,
                                 int16_t duration = 1) {
  CalendarEvent e;
  e.calendar_event_id = id;
  e.start_day = 0;
  e.schedule_type_idx = sched_idx;
  e.compliance_rate = compliance;
  e.duration_days = duration;
  e.catchment_rule_id = 0;
  return e;
}

}  // namespace

// =============================================================================
// Trigger: hop fields and active-event state
// =============================================================================

TEST_CASE("trigger sets hop fields for catchment-rule event") {
  WorldState world = makeCatchmentWorld();
  CalendarEventManager manager({{makeCatchmentEvent(1, 1)}});
  manager.triggerEventsForDay(0, world, world.people, 123, {{0, {0}}});

  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 1);
  CHECK(world.people[0].schedule_hop.return_schedule_id == -1);
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 0);
  CHECK(world.people[0].schedule_hop.repeats_remaining == 0);
  CHECK(manager.hasActiveEvent(world.people[0].id));
  CHECK(manager.stats().triggered == 1);
}

TEST_CASE("trigger sets hop_repeats_remaining from duration_days") {
  WorldState world = makeCatchmentWorld();
  CalendarEventManager manager({{makeCatchmentEvent(1, 1, 1.0f, 3)}});
  manager.triggerEventsForDay(0, world, world.people, 123, {{0, {0}}});
  CHECK(world.people[0].schedule_hop.repeats_remaining == 2);
}

TEST_CASE("trigger skips a person already on a hopped schedule") {
  WorldState world = makeCatchmentWorld();
  world.people[0].schedule_hop.hopped_schedule_id = 3;  // already mid-hop
  CalendarEventManager manager({{makeCatchmentEvent(1, 5)}});
  manager.triggerEventsForDay(0, world, world.people, 123, {{0, {0}}});

  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 3);  // unchanged
  CHECK_FALSE(manager.hasActiveEvent(world.people[0].id));
  CHECK(manager.stats().triggered == 0);
  CHECK(manager.stats().skipped_collision == 1);
}

TEST_CASE("compliance rate 1.0 always triggers, 0.0 never triggers") {
  SUBCASE("compliance 1.0 triggers") {
    WorldState world = makeCatchmentWorld();
    CalendarEventManager manager({{makeCatchmentEvent(1, 1, 1.0f)}});
    manager.triggerEventsForDay(0, world, world.people, 123, {{0, {0}}});
    CHECK(manager.stats().triggered == 1);
    CHECK(manager.stats().skipped_compliance == 0);
  }
  SUBCASE("compliance 0.0 never triggers") {
    WorldState world = makeCatchmentWorld();
    CalendarEventManager manager({{makeCatchmentEvent(1, 1, 0.0f)}});
    manager.triggerEventsForDay(0, world, world.people, 123, {{0, {0}}});
    CHECK(manager.stats().triggered == 0);
    CHECK(manager.stats().skipped_compliance == 1);
  }
}

// =============================================================================
// Catchment-rule attendee resolution
// =============================================================================

TEST_CASE("catchment-rule event triggers people in specified geo_units") {
  WorldState world;
  world.geo_level_names = {"sgu"};
  for (GeoUnitId gid : {10, 11, 12}) {
    GeographicalUnit gu;
    gu.id = gid; gu.parent_id = -1; gu.level_id = 0;
    world.geo_units.push_back(gu);
  }
  world.venue_type_names = {"fair"};
  Venue v; v.id = 0; v.type_id = 0; v.geo_unit_id = 10;
  world.venues.push_back(v);
  world.activity_names = {"Fair_accommodation"};
  world.schedule_type_names = {"regular", "Fair_day_trip"};

  for (int i = 0; i < 3; ++i) {
    Person& p = world.people.emplace_back();
    p.id = i; p.geo_unit_id = 10 + i; p.age = 30.0f;
  }
  world.buildIndices();

  CalendarEvent event;
  event.calendar_event_id = 7;
  event.start_day = 0;
  event.schedule_type_idx = 1;
  event.compliance_rate = 1.0f;
  event.catchment_rule_id = 7;

  CalendarEventManager manager({{event}});
  manager.triggerEventsForDay(0, world, world.people, 42, {{7, {10, 11}}});

  CHECK(manager.stats().triggered == 2);
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 1);
  CHECK(world.people[1].schedule_hop.hopped_schedule_id == 1);
  CHECK(world.people[2].schedule_hop.hopped_schedule_id == -1);
}

TEST_CASE("attendee_filters on catchment event exclude non-matching people") {
  WorldState world;
  world.geo_level_names = {"sgu"};
  GeographicalUnit gu; gu.id = 0; gu.parent_id = -1; gu.level_id = 0;
  world.geo_units.push_back(gu);
  world.venue_type_names = {"fair"};
  Venue v; v.id = 0; v.type_id = 0; v.geo_unit_id = 0;
  world.venues.push_back(v);
  world.activity_names = {"Fair_accommodation"};
  world.schedule_type_names = {"regular", "Fair_day_trip"};
  world.person_property_names = {"age", "sex"};

  for (int i = 0; i < 3; ++i) {
    Person& p = world.people.emplace_back();
    p.id = i; p.geo_unit_id = 0;
    p.age = 10.0f + 20.0f * i;  // 10, 30, 50
  }
  world.buildIndices();

  SelectionCriterion age_filter;
  age_filter.property_path = "age";
  age_filter.operator_type = ">=";
  age_filter.value = 20;

  CalendarEvent event;
  event.calendar_event_id = 1;
  event.start_day = 0;
  event.schedule_type_idx = 1;
  event.compliance_rate = 1.0f;
  event.catchment_rule_id = 0;
  event.attendee_filters = {age_filter};

  CalendarEventManager manager({{event}});
  manager.triggerEventsForDay(0, world, world.people, 42, {{0, {0}}});

  CHECK(manager.stats().triggered == 2);
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == -1);
  CHECK(world.people[1].schedule_hop.hopped_schedule_id == 1);
  CHECK(world.people[2].schedule_hop.hopped_schedule_id == 1);
}

TEST_CASE("sweepCompletedHops removes entry once hop becomes inactive") {
  WorldState world = makeCatchmentWorld();
  CalendarEventManager manager({{makeCatchmentEvent(1, 1, 1.0f, 1)}});

  manager.triggerEventsForDay(0, world, world.people, 123, {{0, {0}}});
  REQUIRE(manager.hasActiveEvent(world.people[0].id));

  world.people[0].schedule_hop.clear();  // isActive() == false

  manager.triggerEventsForDay(99, world, world.people, 123, {});
  CHECK_FALSE(manager.hasActiveEvent(world.people[0].id));
}

// =============================================================================
// Per-event counts and the start-up rule check
// =============================================================================

TEST_CASE("event whose attendee filters select nobody records zero passed and zero triggered") {
  WorldState world = makeCatchmentWorld(5);
  for (auto& p : world.people) p.age = 10.0f;

  SelectionCriterion age_filter;
  age_filter.property_path = "age";
  age_filter.operator_type = ">=";
  age_filter.value = 90;

  CalendarEvent event = makeCatchmentEvent(7, 1);
  event.attendee_filters = {age_filter};
  CalendarEventManager manager({{event}});
  manager.triggerEventsForDay(0, world, world.people, 42, {{0, {0}}});

  REQUIRE(manager.lastDayCounts().size() == 1);
  const auto& c = manager.lastDayCounts()[0];
  CHECK(c.calendar_event_id == 7);
  CHECK(c.in_catchment == 5);
  CHECK(c.passed_filters == 0);
  CHECK(c.already_hopping == 0);
  CHECK(c.refused == 0);
  CHECK(c.hopping == 0);
  CHECK(manager.stats().triggered == 0);
}

TEST_CASE("per-event counts add up: catchment, filters, collisions, refusals, hops") {
  WorldState world = makeCatchmentWorld(200);
  for (size_t i = 0; i < world.people.size(); ++i)
    world.people[i].age = (i % 2 == 0) ? 10.0f : 50.0f;
  // Ten of the adults are already mid-hop.
  int mid_hop = 0;
  for (size_t i = 1; i < world.people.size() && mid_hop < 10; i += 2, ++mid_hop)
    world.people[i].schedule_hop = ScheduleHop::begin(1, -1, 0);

  SelectionCriterion age_filter;
  age_filter.property_path = "age";
  age_filter.operator_type = ">=";
  age_filter.value = 20;

  CalendarEvent event = makeCatchmentEvent(3, 1, 0.5f);
  event.attendee_filters = {age_filter};
  CalendarEventManager manager({{event}});
  manager.triggerEventsForDay(0, world, world.people, 99, {{0, {0}}});

  REQUIRE(manager.lastDayCounts().size() == 1);
  const auto& c = manager.lastDayCounts()[0];
  CHECK(c.in_catchment == 200);
  CHECK(c.passed_filters == 100);
  CHECK(c.already_hopping == 10);
  CHECK(c.in_catchment >= c.passed_filters);
  CHECK(c.passed_filters == c.already_hopping + c.refused + c.hopping);
  CHECK(c.refused > 0);
  CHECK(c.hopping > 0);
  CHECK(c.hopping == manager.stats().triggered);
  CHECK(c.refused == manager.stats().skipped_compliance);
  CHECK(c.already_hopping == manager.stats().skipped_collision);

  // A day with no events leaves nothing to report.
  manager.triggerEventsForDay(5, world, world.people, 99, {{0, {0}}});
  CHECK(manager.lastDayCounts().empty());
}

TEST_CASE("an event naming a catchment rule that is not loaded throws at start-up") {
  CalendarEvent event = makeCatchmentEvent(4, 1);
  event.catchment_rule_id = 999;
  CalendarEventManager manager({{event}});

  const std::unordered_map<int32_t, std::vector<GeoUnitId>> rules = {{0, {0}}};
  try {
    manager.checkCatchmentRules(rules, "rules.csv");
    FAIL("expected a throw");
  } catch (const std::runtime_error& e) {
    CHECK(std::string(e.what()) ==
          "calendar event 4 names catchment_rule_id 999, which is not in "
          "rules.csv");
  }

  // No rules file configured at all is also refused, with its own wording.
  CHECK_THROWS_WITH_AS(manager.checkCatchmentRules({}, ""),
                       doctest::Contains("none is configured"),
                       std::runtime_error);

  // A known rule id passes.
  CalendarEventManager ok({{makeCatchmentEvent(5, 1)}});
  CHECK_NOTHROW(ok.checkCatchmentRules(rules, "rules.csv"));
}

TEST_CASE("event report line format") {
  CalendarEventManager::EventDayCounts c;
  c.calendar_event_id = 1;
  c.catchment_rule_id = 2;
  c.in_catchment = 10;
  c.passed_filters = 8;
  c.already_hopping = 1;
  c.refused = 2;
  c.hopping = 5;
  CHECK(formatCalendarEventLine(c, "2020-03-14") ==
        "[CalendarEvent] id=1 2020-03-14 catchment_rule=2: 10 in catchment, 8 "
        "pass filters, 1 already hopping, 2 refused, 5 hopping\n");
}
