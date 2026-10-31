#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/config.h"
#include "core/types.h"

namespace june {

class WorldState;

// A calendar-triggered event. On `start_day`, attendees drawn from the
// catchment rule hop onto the temporary schedule `schedule_type_idx`. Venue
// assignment is handled by OnTheFlyVenueAllocator at activity-selection time.
struct CalendarEvent {
  int32_t calendar_event_id = -1;
  int start_day = -1;              // sim day (0-based) the event triggers on
  int16_t schedule_type_idx = -1;  // index into ScheduleConfig::schedule_types
  float compliance_rate = 1.0f;    // P(an eligible attendee actually hops)
  int16_t duration_days = 1;       // number of schedule-hop days (default 1)
  int32_t catchment_rule_id = -1;
  GeoUnitId hosting_geo_unit_id = -1;
  std::string venue_type_name;
  std::vector<SelectionCriterion> attendee_filters;
};

// Owns calendar events and the thin per-person active-event-id state.
// Venue resolution is delegated to OnTheFlyVenueAllocator via
// getActiveHostingGeoUnit, which ActivityManager calls at assignment time.
class CalendarEventManager {
 public:
  CalendarEventManager() = default;
  explicit CalendarEventManager(
      std::vector<std::vector<CalendarEvent>> events_by_day);

  // events_by_id_ points into events_by_day_; move only.
  CalendarEventManager(const CalendarEventManager&) = delete;
  CalendarEventManager& operator=(const CalendarEventManager&) = delete;
  CalendarEventManager(CalendarEventManager&&) = default;
  CalendarEventManager& operator=(CalendarEventManager&&) = default;

  // Fire all events scheduled for `day`: for each, resolve attendees via the
  // catchment rule, then for each compliant, non-colliding attendee set the
  // hop fields and record the active calendar_event_id. `base_seed` drives
  // compliance rolls and venue hashing. `catchment_rules` maps
  // catchment_rule_id -> geo_unit list.
  void triggerEventsForDay(
      int day, const WorldState& world, std::vector<Person>& people,
      uint64_t base_seed,
      const std::unordered_map<int32_t, std::vector<GeoUnitId>>&
          catchment_rules = {});

  // Erase active-event entries for persons whose hop has completed (i.e.
  // !schedule_hop.isActive()). Called at the top of triggerEventsForDay so
  // stale entries are cleaned up once per day without coupling to
  // ActivityManager.
  void sweepCompletedHops(
      const std::unordered_map<PersonId, size_t>& person_index,
      const std::vector<Person>& people);

  // Checkpoint seam: serialisable state.
  struct Snapshot {
    std::unordered_map<PersonId, int32_t> active_event;
  };
  Snapshot snapshot_for_checkpoint() const;
  void restore(Snapshot snapshot);

  // Trigger diagnostics.
  struct Stats {
    size_t triggered = 0;
    size_t skipped_collision = 0;
    size_t skipped_compliance = 0;
  };
  const Stats& stats() const { return stats_; }

  // What one event did on the last triggerEventsForDay call, on this rank only.
  // All counts are integers so a sum across ranks does not depend on the rank
  // count. `in_catchment` is everyone the rule's geo units hold, `passed_filters`
  // those who also pass the attendee filters; each of those lands in exactly
  // one of already_hopping, refused or hopping.
  struct EventDayCounts {
    int32_t calendar_event_id = -1;
    int32_t catchment_rule_id = -1;
    uint64_t in_catchment = 0;
    uint64_t passed_filters = 0;
    uint64_t already_hopping = 0;
    uint64_t refused = 0;
    uint64_t hopping = 0;
  };
  // One entry per event scheduled for the day, in file order; empty on a day
  // with no events. Replaced at the start of every triggerEventsForDay call.
  const std::vector<EventDayCounts>& lastDayCounts() const {
    return last_day_counts_;
  }

  bool hasActiveEvent(PersonId person_id) const {
    return active_event_.find(person_id) != active_event_.end();
  }

  // Returns the hosting_geo_unit_id of the person's active event, or nullopt
  // if the person has no active event or the event has no hosting geo-unit.
  std::optional<GeoUnitId> getActiveHostingGeoUnit(PersonId person_id) const;

  // Every event's catchment_rule_id must name a rule. Throws a descriptive
  // error otherwise, so an event that could reach nobody is caught at start-up.
  // `rules_source` is the rules file path, or empty if none was configured.
  void checkCatchmentRules(
      const std::unordered_map<int32_t, std::vector<GeoUnitId>>&
          catchment_rules,
      const std::string& rules_source) const;

  // Distinct hosting geo units across all events (dups harmless). Used by the
  // halo prototype to precompute hosting-strategy OTF pools.
  std::vector<GeoUnitId> hostingGeoUnits() const;

 private:
  // events_by_day_[day] -> events starting that day.
  std::vector<std::vector<CalendarEvent>> events_by_day_;
  // O(1) lookup from calendar_event_id to its definition (pointer into
  // events_by_day_).
  std::unordered_map<int32_t, const CalendarEvent*> events_by_id_;
  // person -> opaque active calendar_event_id (never a venue).
  std::unordered_map<PersonId, int32_t> active_event_;
  Stats stats_;
  std::vector<EventDayCounts> last_day_counts_;

  // Attendees for a catchment-rule event: gathered from people_by_geo_unit for
  // each geo_unit in the catchment rule, filtered by event.attendee_filters.
  // `in_catchment` receives the number of people examined before filtering.
  std::vector<PersonId> attendeesForCatchmentEvent(
      const CalendarEvent& event, const WorldState& world,
      const std::unordered_map<int32_t, std::vector<GeoUnitId>>&
          catchment_rules,
      uint64_t& in_catchment) const;
};

// One report line for an event's day, e.g.
// "[CalendarEvent] id=1 2020-03-14 catchment_rule=1: 10 in catchment, ..."
// `counts` hold the totals summed over every rank.
std::string formatCalendarEventLine(
    const CalendarEventManager::EventDayCounts& counts,
    const std::string& date);

static_assert(!std::is_copy_constructible_v<CalendarEventManager>);

}  // namespace june
