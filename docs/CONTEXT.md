# JUNE2

C++ epidemiological simulation engine. Consumes a serialized world (people,
venues, activities) produced by the MAY world-builder and runs the
disease-transmission/scheduling loop over it.

## Language

**Activity**:
A named category of thing a person can do (e.g. work, school, residence,
`Fair_accommodation`). Identified by a constant `activity_name`/
`activity_index`, shared across all people and all occurrences of that
category.
_Avoid_: activity type (ambiguous with venue type), event (see **Calendar Event**).

**Venue**:
A physical place a person can be assigned to for an activity (e.g. a
specific school, a specific guest-house). Has a `venue_id` and a
`venue.type` drawn from a single global registry capped at 256 distinct
types.

**Subset**:
A fixed group/role within a Venue (e.g. a class within a school, a room's
occupants within a guest-house). Membership is pre-baked at world-build
time and not modified during the simulation run. `subset_index` is local
to its Venue, so the pair `(venue_id, subset_index)` identifies a Subset
uniquely and globally — regardless of which Activity reaches it. This
holds even when a Venue hosts multiple Subsets serving different Calendar
Events (see **Calendar Event**).
_Avoid_: group (too generic), cohort.

**Calendar Event**:
A scheduled occurrence (not an Activity) that some people attend on a
specific date, triggering a temporary schedule hop. Examples include fairs
requiring accommodation away from residence. Identified by a
`calendar_event_id`. Attendees are drawn at trigger time from
`people_by_geo_unit` using a **Catchment rule**, requiring no pre-baked
membership.
_Avoid_: feast, fair (too narrow — Calendar Event is the general concept).

**Catchment rule**:
A geo-unit eligibility list (`catchment_rule_id → [geo_unit_id, …]`)
resolved against `WorldState::people_by_geo_unit` at
calendar-event-trigger time to select Calendar Event attendees. Not a
Subset: membership is never pre-baked; it is recomputed each time the
Calendar Event fires.
_Avoid_: Subset (pre-baked, world-build-time — the opposite of this).

**Coordinated Encounter**:
A live-negotiated Activity assignment (`CoordinatedEncounterDef`, defined in
`coordinated_encounters.yaml`): a host proposes to network partners at a
`trigger_slots` Activity, invitees are drawn from a named network and
liveness-checked, and accepting participants are assigned a Venue from
`allowed_venues`. Replaces the ordinary pre-baked-venue-map resolution for
Activities marked `coordinated_only_activities` in `schedules.yaml` (that map
has no liveness check, so it kept sending visitors to dead contacts).
Optionally sets `hop_schedule`, which moves non-host participants onto a
temporary Schedule Hop on acceptance (e.g. a multi-slot trip) — the host
always stays at their own Venue. No `hop_schedule` means an in-place
encounter: participants are assigned the Venue directly, without hopping.
_Avoid_: "negotiated activity", "live encounter".

**Schedule Hop**:
A temporary or permanent departure from a Person's assigned schedule type,
owned by the `ScheduleHop` struct on Person; active when
`hopped_schedule_id != -1`. A _temporary_ hop advances monotonically through
the hopped schedule's `flat_slots` (the counter is never reset on day-boundary
wrap, so backward scans reach earlier repeats via `% n`); `repeats_remaining`
counts full-cycle repeats still to run (0 = final/only). Auto-return fires when
a full cycle completes with none remaining. Whether a hop is temporary is a
property of its ScheduleType (`is_temporary`), not of the hop fields. All
temporary hops begin at progress = 0, via two call-site onset patterns:
_immediate_ (`maybeTriggerScheduleHop` runs slot 0 then advances, leaving
progress = 1) and _deferred_ (`triggerEventsForDay` stops after `begin()`; slot
0 runs on the first advance). _Permanent_ hops (a property-dispatched permanent
schedule, or a policy freeze-in-place swap) set `hopped_schedule_id` without
auto-return.
_Avoid_: "temp schedule", "hopped state".

**Logical day** (of a hop slot):
The calendar day on which a specific slot within a Schedule Hop was originally
scheduled to execute: `hop_start_day + k / n`, where `k` is
`temp_slot_progress` at time of execution and `n` is `flat_slots.size()`.
Equals `current_sim_day_` on the forward pass but differs when
`findLastNonNullVenueOnHop` re-resolves an earlier slot on a later calendar
day. OTF daily venue seeds must use the logical day so forward and backward
resolution always agree.
_Avoid_: confusing with `current_sim_day_` (the wall-clock day at call time).

**`current_sim_day_`** (simulation wall-clock day):
The calendar day the simulation is currently advancing through at the moment a
function executes. Identical to the logical day for forward-pass slot
resolution, but not for `findLastNonNullVenueOnHop`, which re-resolves
previously-executed slots from a later wall-clock position.
_Avoid_: treating as synonymous with logical day inside hop backward scans.

## Relationships

- A **Person** has, per **Activity**, a list of candidate **Venue**
  (+ **Subset**) references — populated at world-load, never grown during
  the simulation.
- A **Calendar Event** is not itself an **Activity** — it is calendar data
  that triggers a schedule hop into a designated **Activity** (e.g.
  `Fair_accommodation`). The specific **Venue** is resolved at trigger
  time: attendees drawn from `people_by_geo_unit` using a **Catchment
  rule**; no pre-baked membership required.

**Venue assignment strategy**:
What determines which Venue an attendee occupies during a calendar-triggered
schedule hop. The eligible Venue pool is the venues of `venue_type_name`
located in the geo-unit determined by the OTF rule (hosting or resident,
depending on the rule's strategy). One Venue is then chosen from that pool for
a specific Person by deterministic hash-select.
_Avoid_: venue resolution strategy (resolution is the act of calling the
strategy, not the strategy itself).

**Visitor**:
A person attending a venue owned by a different MPI rank for a given time
slot. The person's disease state is sent to the venue's owning rank via MPI
exchange; transmission is computed there. If the visitor is infected, a
pending infection is routed back to the person's home rank and applied. The
person's `active_event_` entry (if any) and all simulation state are held
exclusively on their home rank.
_Avoid_: ghost, halo (no replication of person state across ranks).

## Flagged ambiguities

- An earlier session used the term "Feast" for this concept before
  settling on **Fair** (2026-06-19) — if "Feast" appears in older notes
  or branches, read it as **Calendar Event**.
- "Fair" was used as the canonical term from 2026-06-19 until 2026-06-29,
  when it was replaced by **Calendar Event** to reflect the broader concept
  (Fair is one example). If "Fair" appears in older notes, ADRs, or
  branches, read it as **Calendar Event**.
