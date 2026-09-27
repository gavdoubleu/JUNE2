#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "activity/activity_manager.h"
#include "core/config.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/policy.h"
#include "test_utils.h"
#include "transmission_fixtures.h"

using namespace june;

// Helper: pre-resolve a TimeSlot's allowed_activity_indices from world
static void resolveSlotIndices(TimeSlot& slot, const WorldState& world) {
  slot.allowed_activity_indices.clear();
  for (const auto& act : slot.allowed_activities) {
    int idx = world.getActivityIndex(act);
    if (idx >= 0)
      slot.allowed_activity_indices.push_back(static_cast<int16_t>(idx));
  }
}

// =============================================================================
// Cycle 1: hop_on_activity config resolution
// =============================================================================

TEST_CASE("hop_on_activity resolves to schedule index") {
  WorldState world = TestWorldFactory::createMinimalWorld(1, 1);
  world.activity_names = {"residence", "work", "none", "dead"};

  ScheduleConfig sched_config;
  sched_config.day_type_cycle = {"workday"};
  sched_config.day_type_names = {"workday"};

  ScheduleType regular;
  regular.name = "regular";
  ScheduleType work_sched;
  work_sched.name = "work_schedule";

  TimeSlot slot;
  slot.name = "morning";
  slot.allowed_activities = {"work", "residence"};
  slot.hop_on_activity["work"] = "work_schedule";
  regular.slots_by_day_type["workday"].push_back(slot);

  sched_config.schedule_types.push_back(regular);
  sched_config.schedule_types.push_back(work_sched);

  sched_config.resolveSlots(world);

  const TimeSlot& resolved = sched_config.schedule_types[0].slots_by_day_type["workday"][0];

  CHECK(resolved.hop_schedule_by_activity_idx[1] == 1);  // "work" -> "work_schedule" (idx 1)
  CHECK(resolved.hop_schedule_by_activity_idx[0] == -1); // "residence" -> no hop
}

// =============================================================================
// Cycle 2: activity triggers immediate hop — person assigned from flat_slots[0]
// =============================================================================

TEST_CASE("activity triggers immediate hop to new schedule slot 0") {
  // activity layout: residence=0, work=1, special_work=2, none=3, dead=4
  WorldState world = TestWorldFactory::createMinimalWorld(1, 2);
  world.activity_names = {"residence", "work", "special_work", "none", "dead"};
  world.venue_type_names = {"office", "special_office"};
  world.venues[0].type_id = 0;  // office (work venue)
  world.venues[1].type_id = 1;  // special_office (special_work venue)
  world.buildIndices();

  // Person has work (venue 0) and special_work (venue 1)
  world.people[0].activity_meta_start =
      static_cast<uint32_t>(world.activity_meta.size());
  world.people[0].activity_meta_count = 2;
  world.activity_meta.push_back(
      {1, static_cast<uint32_t>(world.activity_venues.size()), 1});  // work
  world.activity_venues.push_back({0, 0});
  world.activity_meta.push_back(
      {2, static_cast<uint32_t>(world.activity_venues.size()), 1});  // special_work
  world.activity_venues.push_back({1, 0});

  // Schedule 0: "regular" — has a slot that hops to temp_sched on "work"
  ScheduleType regular;
  regular.name = "regular";
  TimeSlot trigger_slot;
  trigger_slot.name = "morning";
  trigger_slot.allowed_activities = {"work"};
  trigger_slot.hop_on_activity["work"] = "temp_sched";
  resolveSlotIndices(trigger_slot, world);
  regular.slots_by_day_type["workday"].push_back(trigger_slot);

  // Schedule 1: "temp_sched" (temporary) — flat_slots[0] has special_work only
  ScheduleType temp_sched;
  temp_sched.name = "temp_sched";
  temp_sched.is_temporary = true;
  TimeSlot temp_slot0;
  temp_slot0.name = "temp_morning";
  temp_slot0.allowed_activities = {"special_work"};
  resolveSlotIndices(temp_slot0, world);
  temp_sched.flat_slots.push_back(temp_slot0);

  Config config;
  config.schedule.day_type_cycle = {"workday"};
  config.schedule.day_type_names = {"workday"};
  config.schedule.schedule_types.push_back(regular);
  config.schedule.schedule_types.push_back(temp_sched);
  config.schedule.default_schedule_type = "regular";
  config.performance.precompute_schedules = false;

  // Resolve so hop_schedule_by_activity_idx is populated
  config.resolve(world);

  // Precomputed: person 0 deterministically assigned "work" (idx 1) at venue 0
  world.num_day_types = 1;
  world.precomputed_schedules.resize(1);
  world.precomputed_schedules[0].push_back(ScheduleEntry(1, 0, 0, true));  // work
  world.schedule_starts.assign(1, 0);
  world.schedule_counts.assign(1, 1);
  world.people[0].schedule_computed = true;
  world.people[0].schedule_type_id = 0;  // "regular"

  world.schedule_type_names = {"regular", "temp_sched"};
  ActivityManager manager(world, config);
  manager.assignScheduleTypes();

  std::vector<PersonLocation> locations(1);
  manager.assignActivitiesFromSchedule(0, 0, locations);

  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 1);   // hopped to temp_sched
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 1);   // advanced past slot 0
  CHECK(locations[0].activity_index == 2);           // special_work (from flat_slots[0])
}

// =============================================================================
// Cycle 3: non-triggering activity leaves person unhopped
// =============================================================================

TEST_CASE("non-triggering activity does not hop") {
  // Same setup as Cycle 2 but precomputed entry forces "residence" (idx 0)
  WorldState world = TestWorldFactory::createMinimalWorld(1, 2);
  world.activity_names = {"residence", "work", "special_work", "none", "dead"};
  world.venue_type_names = {"home", "office", "special_office"};
  world.venues[0].type_id = 0;
  world.venues[1].type_id = 1;
  world.buildIndices();

  // Person has residence (venue 0) and work (venue 1)
  world.people[0].activity_meta_start =
      static_cast<uint32_t>(world.activity_meta.size());
  world.people[0].activity_meta_count = 2;
  world.activity_meta.push_back(
      {0, static_cast<uint32_t>(world.activity_venues.size()), 1});  // residence
  world.activity_venues.push_back({0, 0});
  world.activity_meta.push_back(
      {1, static_cast<uint32_t>(world.activity_venues.size()), 1});  // work
  world.activity_venues.push_back({1, 0});

  ScheduleType regular;
  regular.name = "regular";
  TimeSlot trigger_slot;
  trigger_slot.name = "morning";
  trigger_slot.allowed_activities = {"work", "residence"};
  trigger_slot.hop_on_activity["work"] = "temp_sched";
  resolveSlotIndices(trigger_slot, world);
  regular.slots_by_day_type["workday"].push_back(trigger_slot);

  ScheduleType temp_sched;
  temp_sched.name = "temp_sched";
  temp_sched.is_temporary = true;
  TimeSlot temp_slot0;
  temp_slot0.name = "temp_morning";
  temp_slot0.allowed_activities = {"special_work"};
  resolveSlotIndices(temp_slot0, world);
  temp_sched.flat_slots.push_back(temp_slot0);

  Config config;
  config.schedule.day_type_cycle = {"workday"};
  config.schedule.day_type_names = {"workday"};
  config.schedule.schedule_types.push_back(regular);
  config.schedule.schedule_types.push_back(temp_sched);
  config.schedule.default_schedule_type = "regular";
  config.performance.precompute_schedules = false;
  config.resolve(world);

  // Precomputed: person assigned "residence" (idx 0) — no hop should trigger
  world.num_day_types = 1;
  world.precomputed_schedules.resize(1);
  world.precomputed_schedules[0].push_back(ScheduleEntry(0, 0, 0, true));  // residence
  world.schedule_starts.assign(1, 0);
  world.schedule_counts.assign(1, 1);
  world.people[0].schedule_computed = true;
  world.people[0].schedule_type_id = 0;
  world.schedule_type_names = {"regular", "temp_sched"};

  ActivityManager manager(world, config);
  manager.assignScheduleTypes();

  std::vector<PersonLocation> locations(1);
  manager.assignActivitiesFromSchedule(0, 0, locations);

  CHECK(world.people[0].schedule_hop.hopped_schedule_id == -1);  // no hop
  CHECK(locations[0].activity_index == 0);           // stayed at residence
}

// =============================================================================
// Cycle 4: person on temp schedule advances temp_slot_progress each call
// =============================================================================

TEST_CASE("person on temp schedule advances through flat_slots") {
  // 3-slot temp schedule; 2 calls advance progress without exhausting all slots
  WorldState world = TestWorldFactory::createMinimalWorld(1, 2);
  world.activity_names = {"residence", "task_a", "task_b", "none", "dead"};
  world.venue_type_names = {"home", "venue_a", "venue_b"};
  world.venues[0].type_id = 1;  // venue_a
  world.venues[1].type_id = 2;  // venue_b
  world.buildIndices();

  world.people[0].activity_meta_start =
      static_cast<uint32_t>(world.activity_meta.size());
  world.people[0].activity_meta_count = 2;
  world.activity_meta.push_back(
      {1, static_cast<uint32_t>(world.activity_venues.size()), 1});  // task_a
  world.activity_venues.push_back({0, 0});
  world.activity_meta.push_back(
      {2, static_cast<uint32_t>(world.activity_venues.size()), 1});  // task_b
  world.activity_venues.push_back({1, 0});

  ScheduleType temp_sched;
  temp_sched.name = "temp_sched";
  temp_sched.is_temporary = true;
  TimeSlot slot_a, slot_b, slot_a2;
  slot_a.name = "step_a";
  slot_a.allowed_activities = {"task_a"};
  resolveSlotIndices(slot_a, world);
  slot_b.name = "step_b";
  slot_b.allowed_activities = {"task_b"};
  resolveSlotIndices(slot_b, world);
  slot_a2.name = "step_a2";
  slot_a2.allowed_activities = {"task_a"};
  resolveSlotIndices(slot_a2, world);
  temp_sched.flat_slots.push_back(slot_a);
  temp_sched.flat_slots.push_back(slot_b);
  temp_sched.flat_slots.push_back(slot_a2);  // 3rd slot keeps person hopped after Call 2

  Config config;
  config.schedule.day_type_cycle = {"workday"};
  config.schedule.day_type_names = {"workday"};
  config.schedule.schedule_types.push_back(temp_sched);
  config.schedule.default_schedule_type = "temp_sched";
  config.resolve(world);

  // Manually put the person in a hopped state at progress=0
  world.people[0].schedule_hop.hopped_schedule_id = 0;   // temp_sched is index 0
  world.people[0].schedule_hop.return_schedule_id = -1;
  world.people[0].schedule_hop.temp_slot_progress = 0;
  world.people[0].schedule_computed = true;
  world.people[0].schedule_type_id = 0;
  world.schedule_type_names = {"temp_sched"};
  // Precomputed schedule must exist (but hopped path bypasses it)
  world.num_day_types = 1;
  world.precomputed_schedules.resize(1);
  world.precomputed_schedules[0].push_back(ScheduleEntry(0, -1, -1, false));
  world.schedule_starts.assign(1, 0);
  world.schedule_counts.assign(1, 1);

  ActivityManager manager(world, config);
  manager.assignScheduleTypes();

  std::vector<PersonLocation> locations(1);

  // Call 1: uses flat_slots[0] (task_a)
  manager.assignActivitiesFromSchedule(0, 0, locations);
  CHECK(locations[0].activity_index == 1);                      // task_a
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 1);
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 0);               // still hopped

  // Call 2: uses flat_slots[1] (task_b), still before end
  manager.assignActivitiesFromSchedule(0, 0, locations);
  CHECK(locations[0].activity_index == 2);                      // task_b
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 2);
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 0);               // still hopped (return happens at boundary)
}

// =============================================================================
// Cycle 5: auto-return to original schedule after last flat_slot
// =============================================================================

TEST_CASE("temp schedule auto-returns to original schedule after last slot") {
  WorldState world = TestWorldFactory::createMinimalWorld(1, 2);
  world.activity_names = {"residence", "task_a", "none", "dead"};
  world.venue_type_names = {"home", "venue_a"};
  world.venues[0].type_id = 0;  // home
  world.venues[1].type_id = 1;  // venue_a
  world.buildIndices();

  // Person has residence (venue 0) and task_a (venue 1)
  world.people[0].activity_meta_start =
      static_cast<uint32_t>(world.activity_meta.size());
  world.people[0].activity_meta_count = 2;
  world.activity_meta.push_back(
      {0, static_cast<uint32_t>(world.activity_venues.size()), 1});  // residence
  world.activity_venues.push_back({0, 0});
  world.activity_meta.push_back(
      {1, static_cast<uint32_t>(world.activity_venues.size()), 1});  // task_a
  world.activity_venues.push_back({1, 0});

  // Schedule 0: "regular" (the original)
  ScheduleType regular;
  regular.name = "regular";
  TimeSlot reg_slot;
  reg_slot.name = "day";
  reg_slot.allowed_activities = {"residence"};
  resolveSlotIndices(reg_slot, world);
  regular.slots_by_day_type["workday"].push_back(reg_slot);

  // Schedule 1: "temp_sched" — 1 flat_slot only
  ScheduleType temp_sched;
  temp_sched.name = "temp_sched";
  temp_sched.is_temporary = true;
  TimeSlot temp_slot;
  temp_slot.name = "step_a";
  temp_slot.allowed_activities = {"task_a"};
  resolveSlotIndices(temp_slot, world);
  temp_sched.flat_slots.push_back(temp_slot);

  Config config;
  config.schedule.day_type_cycle = {"workday"};
  config.schedule.day_type_names = {"workday"};
  config.schedule.schedule_types.push_back(regular);
  config.schedule.schedule_types.push_back(temp_sched);
  config.schedule.default_schedule_type = "regular";
  config.resolve(world);

  world.schedule_type_names = {"regular", "temp_sched"};
  world.num_day_types = 1;
  world.precomputed_schedules.resize(1);
  world.precomputed_schedules[0].push_back(ScheduleEntry(0, 0, 0, true));
  world.schedule_starts.assign(1, 0);
  world.schedule_counts.assign(1, 1);
  world.people[0].schedule_computed = true;
  world.people[0].schedule_type_id = 0;  // "regular" is the original

  // Put person in hopped state: on temp_sched (idx 1), 1 slot remaining
  world.people[0].schedule_hop.hopped_schedule_id = 1;
  world.people[0].schedule_hop.return_schedule_id = -1;  // default: return to original
  world.people[0].schedule_hop.temp_slot_progress = 0;

  ActivityManager manager(world, config);
  manager.assignScheduleTypes();

  std::vector<PersonLocation> locations(1);

  // Executing the 1 temp slot triggers return
  manager.assignActivitiesFromSchedule(0, 0, locations);

  CHECK(world.people[0].schedule_hop.hopped_schedule_id == -1);                     // returned
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 0);                      // reset
  CHECK(world.people[0].cached_schedule_type_->name == "regular");     // back to original
}

// =============================================================================
// Cycle 6: return_schedule sends person to a specified schedule, not the original
// =============================================================================

TEST_CASE("temp schedule returns to specified return_schedule") {
  WorldState world = TestWorldFactory::createMinimalWorld(1, 2);
  world.activity_names = {"residence", "task_a", "none", "dead"};
  world.venue_type_names = {"home", "venue_a"};
  world.venues[0].type_id = 0;
  world.venues[1].type_id = 1;
  world.buildIndices();

  world.people[0].activity_meta_start =
      static_cast<uint32_t>(world.activity_meta.size());
  world.people[0].activity_meta_count = 2;
  world.activity_meta.push_back(
      {0, static_cast<uint32_t>(world.activity_venues.size()), 1});  // residence
  world.activity_venues.push_back({0, 0});
  world.activity_meta.push_back(
      {1, static_cast<uint32_t>(world.activity_venues.size()), 1});  // task_a
  world.activity_venues.push_back({1, 0});

  // Schedule 0: "original"
  ScheduleType original;
  original.name = "original";
  TimeSlot orig_slot;
  orig_slot.name = "day";
  orig_slot.allowed_activities = {"residence"};
  resolveSlotIndices(orig_slot, world);
  original.slots_by_day_type["workday"].push_back(orig_slot);

  // Schedule 1: "special_return" — the intended return target
  ScheduleType special_return;
  special_return.name = "special_return";
  TimeSlot sr_slot;
  sr_slot.name = "day";
  sr_slot.allowed_activities = {"residence"};
  resolveSlotIndices(sr_slot, world);
  special_return.slots_by_day_type["workday"].push_back(sr_slot);

  // Schedule 2: "temp_sched" — 1 flat_slot, returns to "special_return"
  ScheduleType temp_sched;
  temp_sched.name = "temp_sched";
  temp_sched.is_temporary = true;
  temp_sched.return_schedule = "special_return";
  TimeSlot temp_slot;
  temp_slot.name = "step";
  temp_slot.allowed_activities = {"task_a"};
  resolveSlotIndices(temp_slot, world);
  temp_sched.flat_slots.push_back(temp_slot);

  Config config;
  config.schedule.day_type_cycle = {"workday"};
  config.schedule.day_type_names = {"workday"};
  config.schedule.schedule_types.push_back(original);
  config.schedule.schedule_types.push_back(special_return);
  config.schedule.schedule_types.push_back(temp_sched);
  config.schedule.default_schedule_type = "original";
  config.resolve(world);

  world.schedule_type_names = {"original", "special_return", "temp_sched"};
  world.num_day_types = 1;
  world.precomputed_schedules.resize(1);
  world.precomputed_schedules[0].push_back(ScheduleEntry(0, 0, 0, true));
  world.schedule_starts.assign(1, 0);
  world.schedule_counts.assign(1, 1);
  world.people[0].schedule_computed = true;
  world.people[0].schedule_type_id = 0;  // "original"

  // Hop to temp_sched; return_schedule_id must be set from return_schedule_idx
  world.people[0].schedule_hop.hopped_schedule_id = 2;  // temp_sched
  int16_t return_idx = config.schedule.schedule_types[2].return_schedule_idx;
  world.people[0].schedule_hop.return_schedule_id = return_idx;
  world.people[0].schedule_hop.temp_slot_progress = 0;

  ActivityManager manager(world, config);
  manager.assignScheduleTypes();

  std::vector<PersonLocation> locations(1);
  manager.assignActivitiesFromSchedule(0, 0, locations);

  REQUIRE(return_idx == 1);  // "special_return" is index 1
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == -1);
  CHECK(world.people[0].cached_schedule_type_->name == "special_return");
}

// =============================================================================
// Cycle 7: hop_repeats_remaining — multi-day event stays hopped across days
// =============================================================================
// NOTE: with the monotonic-progress fix, temp_slot_progress is no longer reset
// to 0 on day-boundary wrap — it keeps incrementing so findLastNonNullVenueOnHop
// can scan across boundaries via k % n.  Mid-hop assertions reflect this.
//

TEST_CASE("temp schedule repeats N times before returning when hop_repeats_remaining > 0") {
  WorldState world = TestWorldFactory::createMinimalWorld(1, 2);
  world.activity_names = {"residence", "task_a", "none", "dead"};
  world.venue_type_names = {"home", "venue_a"};
  world.venues[0].type_id = 0;
  world.venues[1].type_id = 1;
  world.buildIndices();

  world.people[0].activity_meta_start =
      static_cast<uint32_t>(world.activity_meta.size());
  world.people[0].activity_meta_count = 2;
  world.activity_meta.push_back(
      {0, static_cast<uint32_t>(world.activity_venues.size()), 1});
  world.activity_venues.push_back({0, 0});
  world.activity_meta.push_back(
      {1, static_cast<uint32_t>(world.activity_venues.size()), 1});
  world.activity_venues.push_back({1, 0});

  ScheduleType regular;
  regular.name = "regular";
  TimeSlot reg_slot;
  reg_slot.name = "day";
  reg_slot.allowed_activities = {"residence"};
  resolveSlotIndices(reg_slot, world);
  regular.slots_by_day_type["workday"].push_back(reg_slot);

  ScheduleType temp_sched;
  temp_sched.name = "temp_sched";
  temp_sched.is_temporary = true;
  TimeSlot temp_slot;
  temp_slot.name = "step_a";
  temp_slot.allowed_activities = {"task_a"};
  resolveSlotIndices(temp_slot, world);
  temp_sched.flat_slots.push_back(temp_slot);  // 1 slot per "day"

  Config config;
  config.schedule.day_type_cycle = {"workday"};
  config.schedule.day_type_names = {"workday"};
  config.schedule.schedule_types.push_back(regular);
  config.schedule.schedule_types.push_back(temp_sched);
  config.schedule.default_schedule_type = "regular";
  config.resolve(world);

  world.schedule_type_names = {"regular", "temp_sched"};
  world.num_day_types = 1;

  world.people[0].schedule_hop.hopped_schedule_id = 1;   // temp_sched
  world.people[0].schedule_hop.return_schedule_id = -1;
  world.people[0].schedule_hop.temp_slot_progress = 0;
  world.people[0].schedule_hop.repeats_remaining = 2; // 2 remaining loops after first

  ActivityManager manager(world, config);
  manager.assignScheduleTypes();

  std::vector<PersonLocation> locations(1);

  // Day 1: exhausts the 1-slot schedule, hop_repeats_remaining: 2 -> 1, stays hopped
  manager.assignActivitiesFromSchedule(0, 0, locations);
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 1);   // still hopped
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 1);   // monotonic: not reset
  CHECK(world.people[0].schedule_hop.repeats_remaining == 1);

  // Day 2: decrement to 0, still loops (>0 check was pre-decrement)
  manager.assignActivitiesFromSchedule(0, 0, locations);
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 1);
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 2);   // monotonic: not reset
  CHECK(world.people[0].schedule_hop.repeats_remaining == 0);

  // Day 3: hop_repeats_remaining is 0, so this exhaustion triggers return
  manager.assignActivitiesFromSchedule(0, 0, locations);
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == -1);   // returned
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 0);
}

// =============================================================================
// Cycle 8: day-boundary wrap with transit slot — temp_slot_progress must stay
// monotonically increasing so findLastNonNullVenueOnHop can scan back across
// the wrap via (k % n).  A reset to 0 leaves k = -2 < 0 → empty scan →
// wrong home fallback when a policy freeze fires on a transit slot.
// =============================================================================

TEST_CASE("multi-day hop keeps monotonic temp_slot_progress across day-boundary wrap") {
  // activities: residence=0, lodging=1, no_venue=2, none=3, dead=4
  WorldState world = TestWorldFactory::createMinimalWorld(1, 2);
  world.activity_names = {"residence", "lodging", "no_venue", "none", "dead"};
  world.venue_type_names = {"home", "guest_house"};
  world.venues[0].type_id = 0;
  world.venues[1].type_id = 1;
  world.buildIndices();

  // Person: residence at venue 0, lodging at venue 1
  world.people[0].activity_meta_start =
      static_cast<uint32_t>(world.activity_meta.size());
  world.people[0].activity_meta_count = 2;
  world.activity_meta.push_back(
      {0, static_cast<uint32_t>(world.activity_venues.size()), 1});
  world.activity_venues.push_back({0, 0});
  world.activity_meta.push_back(
      {1, static_cast<uint32_t>(world.activity_venues.size()), 1});
  world.activity_venues.push_back({1, 0});

  // 2-slot daily cycle: [transit (no_venue), overnight (lodging)]
  ScheduleType temp_sched;
  temp_sched.name = "temp_sched";
  temp_sched.is_temporary = true;
  TimeSlot transit_slot, lodge_slot;
  transit_slot.name = "transit";
  transit_slot.allowed_activities = {"no_venue"};
  resolveSlotIndices(transit_slot, world);
  lodge_slot.name = "at_lodge";
  lodge_slot.allowed_activities = {"lodging"};
  resolveSlotIndices(lodge_slot, world);
  temp_sched.flat_slots.push_back(transit_slot);
  temp_sched.flat_slots.push_back(lodge_slot);

  Config config;
  config.schedule.day_type_cycle = {"workday"};
  config.schedule.day_type_names = {"workday"};
  config.schedule.schedule_types.push_back(temp_sched);
  config.schedule.default_schedule_type = "temp_sched";
  config.resolve(world);

  world.schedule_type_names = {"temp_sched"};
  world.num_day_types = 1;

  world.people[0].schedule_hop.hopped_schedule_id = 0;
  world.people[0].schedule_hop.return_schedule_id = -1;
  world.people[0].schedule_hop.temp_slot_progress = 0;
  world.people[0].schedule_hop.repeats_remaining = 1;  // 2 full day-cycles
  world.people[0].schedule_computed = true;
  world.people[0].schedule_type_id = 0;

  ActivityManager manager(world, config);
  manager.assignScheduleTypes();
  std::vector<PersonLocation> locations(1);

  // Day 1 slot 0: transit → venue = -1
  manager.assignActivitiesFromSchedule(0, 0, locations);
  CHECK(locations[0].venue_id == -1);
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 1);

  // Day 1 slot 1: lodging → venue = 1; day-boundary wrap occurs here
  manager.assignActivitiesFromSchedule(0, 0, locations);
  CHECK(locations[0].venue_id == 1);
  CHECK(world.people[0].schedule_hop.repeats_remaining == 0);
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 0);
  // Monotonic: progress must be 2 so findLastNonNullVenueOnHop starts at
  // k = 2-2 = 0 → s = 0 % 2 = 0 (transit, skipped) → k = -1 stop; but the
  // key correctness is the Day 2 transit slot below.
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 2);  // RED before fix (was 0)

  // Day 2 slot 0: transit again after wrap; progress must be 3 so that a
  // findLastNonNullVenueOnHop scan starts at k=1, s=1%2=1 (lodging) → returns
  // lodge venue (1) rather than home (bug: scan started at k=-1, empty → home).
  manager.assignActivitiesFromSchedule(0, 0, locations);
  CHECK(locations[0].venue_id == -1);
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 3);  // RED before fix (was 1)
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == 0);

  // Day 2 slot 1: lodging, then hop ends (repeats exhausted)
  manager.assignActivitiesFromSchedule(0, 0, locations);
  CHECK(locations[0].venue_id == 1);
  CHECK(world.people[0].schedule_hop.hopped_schedule_id == -1);
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 0);  // reset on hop end only
}

// =============================================================================
// Cycle 9: per-day-type participation — back-scan must replay the venue the
// forward path assigned, not one resolved under the wrong day type.
//
// Regression for the latent day_type_idx divergence: findLastNonNullVenueOnHop
// used to pass day_type_idx = -1, falling back to empty participation, so for a
// temp schedule whose overnight slot picks between two activities by day-type
// participation it back-scanned the WRONG activity → wrong venue → policy pin at
// a place never visited.  With both paths sharing resolveHopSlot and a real
// per-slot day type, the freeze pin must land on the day-0 overnight venue.
// =============================================================================

TEST_CASE("back-scan pins the venue forward path assigned under per-day-type "
          "participation") {
  // activities: residence=0, lodge_a=1, lodge_b=2, no_venue=3, none=4, dead=5
  WorldState world = TestWorldFactory::createMinimalWorld(1, 3);
  world.activity_names = {"residence", "lodge_a", "lodge_b",
                          "no_venue",  "none",    "dead"};
  world.venue_type_names = {"home", "guest_a", "guest_b"};
  world.venues[0].type_id = 0;  // home   -> venue 0
  world.venues[1].type_id = 1;  // lodge_a -> venue 1
  world.venues[2].type_id = 2;  // lodge_b -> venue 2
  world.buildIndices();

  // Person venue mappings: residence@0, lodge_a@1, lodge_b@2
  world.people[0].activity_meta_start =
      static_cast<uint32_t>(world.activity_meta.size());
  world.people[0].activity_meta_count = 3;
  world.activity_meta.push_back(
      {0, static_cast<uint32_t>(world.activity_venues.size()), 1});
  world.activity_venues.push_back({0, 0});
  world.activity_meta.push_back(
      {1, static_cast<uint32_t>(world.activity_venues.size()), 1});
  world.activity_venues.push_back({1, 0});
  world.activity_meta.push_back(
      {2, static_cast<uint32_t>(world.activity_venues.size()), 1});
  world.activity_venues.push_back({2, 0});

  // 2-slot daily cycle: [transit (no_venue), overnight (lodge_a OR lodge_b)].
  // Day type dtA picks lodge_a, dtB picks lodge_b via participation = 1.0.
  ScheduleType temp_sched;
  temp_sched.name = "temp_sched";
  temp_sched.is_temporary = true;
  temp_sched.participation_by_day_type["dtA"] = {{"lodge_a", 1.0}};
  temp_sched.participation_by_day_type["dtB"] = {{"lodge_b", 1.0}};
  TimeSlot transit_slot, overnight_slot;
  transit_slot.name = "transit";
  transit_slot.allowed_activities = {"no_venue"};
  resolveSlotIndices(transit_slot, world);
  overnight_slot.name = "overnight";
  overnight_slot.allowed_activities = {"lodge_a", "lodge_b"};
  resolveSlotIndices(overnight_slot, world);
  temp_sched.flat_slots.push_back(transit_slot);
  temp_sched.flat_slots.push_back(overnight_slot);

  // A second schedule the freeze policy hops the person into when it fires.
  ScheduleType freeze_sched;
  freeze_sched.name = "freeze_sched";

  Config config;
  config.schedule.day_type_cycle = {"dtA", "dtB"};
  config.schedule.day_type_names = {"dtA", "dtB"};
  config.schedule.schedule_types.push_back(temp_sched);
  config.schedule.schedule_types.push_back(freeze_sched);
  config.schedule.default_schedule_type = "temp_sched";
  config.resolve(world);

  world.schedule_type_names = {"temp_sched", "freeze_sched"};
  world.num_day_types = 2;

  world.people[0].schedule_hop.hopped_schedule_id = 0;
  world.people[0].schedule_hop.return_schedule_id = -1;
  world.people[0].schedule_hop.temp_slot_progress = 0;
  world.people[0].schedule_hop.repeats_remaining = 1;  // stay hopped across day boundary
  world.people[0].schedule_computed = true;
  world.people[0].schedule_type_id = 0;

  ActivityManager manager(world, config);
  manager.assignScheduleTypes();
  std::vector<PersonLocation> locations(1);

  // --- Day 0 (dtA): run both slots with NO policy attached, building hop
  //     history [transit -> none, overnight -> lodge_a (venue 1)].
  manager.setCurrentDay(0);
  manager.assignActivitiesFromSchedule(0, 0, locations);  // transit
  CHECK(locations[0].venue_id == -1);
  manager.assignActivitiesFromSchedule(1, 0, locations);  // overnight -> dtA
  CHECK(locations[0].venue_id == 1);  // lodge_a venue, NOT lodge_b
  CHECK(world.people[0].schedule_hop.temp_slot_progress == 2);

  // --- Build a symptom freeze policy that fires on the no_venue transit slot
  //     and pins the person at their last real overnight venue.
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(1.0);
  trans.stage_curves["sick"] = curve;
  trans.symptom_id_curves = {nullptr, curve};
  std::vector<SymptomTag> symptom_tags = {{"healthy", -1, 0}, {"sick", 1, 1}};
  DiseaseStageSettings stage_settings;
  stage_settings.recovered_stages = {"healthy"};
  std::vector<TrajectoryDefinition> trajectories;
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"sick", {"constant", {{"value", 100.0}}}});
  trajectories.push_back(td);
  Disease disease("TestDisease", symptom_tags, stage_settings, trajectories, {},
                  trans);

  PolicyManager pm(world);
  SymptomPolicy freeze;
  freeze.name = "sick_traveller_freeze";
  freeze.trigger_symptoms = {"sick"};
  freeze.action.override_activities = {"no_venue"};
  freeze.action.replacement_activity = "residence";
  freeze.action.replacement_schedule = "freeze_sched";
  freeze.action.compliance_rate = 1.0;
  pm.addSymptomPolicy(freeze);
  pm.resolveAll(disease);

  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 42, kNoTransmissionContext, nullptr,
      "guest_a", 0);
  world.people[0].applicable_symptom_policy_mask = 1;

  // --- Day 1 (dtB): policy now active. Transit slot resolves no_venue, the
  //     freeze fires and pins the person at the back-scanned overnight venue.
  manager.setPolicyManager(&pm);
  manager.setCurrentTime(1.0);  // within the sick window -> policy triggers
  manager.setCurrentDay(1);
  manager.assignActivitiesFromSchedule(0, 1, locations);  // transit, dtB

  // The last real venue before this transit was day-0 overnight = lodge_a
  // (venue 1).  The buggy -1 day_type back-scan instead resolved lodge_b
  // (venue 2) via empty participation falling through to the last activity.
  CHECK(locations[0].venue_id == 1);
}
