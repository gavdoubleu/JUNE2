# june_events

Python library for reading `simulation_events.h5`, the event-log output of the
JUNE2 simulation engine (see [JUNE2 simulation engine](../../docs/CONTEXT.md)).
Provides functions to load, inspect, decode, and join the file's tables;
leaves plotting and analysis logic to the caller.

## Language

**Event type**:
One of the named tables under `events/` in `simulation_events.h5` (e.g.
`infections`, `deaths`, `symptom_changes`, `coordinated_encounters`) — one row
per occurrence, always carrying a `time` field. Not every file has every
event type.
_Avoid_: event group (HDF5 term, not a domain concept here).

**Raw event table**:
The unmodified contents of one `events/*` dataset, loaded as a DataFrame with
exactly the columns the file stores (byte-string columns undecoded, id
columns undecoded). Produced by `io/`. Distinct from an **Enrichment**, which
augments a raw event table with data from elsewhere.
_Avoid_: events dataframe (ambiguous with enriched output).

**Lookup table**:
A `lookups/*` table (`people`, `venues`, `people_properties/*`,
`population_summary`) that maps an id (`person_id`, `venue_id`) to static
metadata about that entity — one row per id. Loaded by `io/`, joined onto
event tables by `enrich/`. A duplicate id is treated as corruption, not a
legitimate multi-row entity: `enrich/` raises on it by default. Callers who
need to force the merge through anyway (e.g. to inspect corrupt data rather
than fail outright) can opt out via an explicit parameter — the fan-out is
never silent.
`people_properties/*` is a special case: each property (`ethnicity`,
`comorbidities`, etc.) is its own same-length dataset with no `person_id`
column of its own — row *i* corresponds to row *i* of `lookups/people` by
construction (both written from the same ordered iteration in the engine's
event writer), not by an id join. Unlike a duplicate-id in a joined
**Lookup table**, a length mismatch here has no legitimate reading and
nothing meaningful to force through — `io/` raises unconditionally, no
escape hatch.

**Registry**:
An ordered list of strings under `metadata/registries/*` (e.g.
`encounter_types`, `symptoms`) that a uint-coded column indexes into. The
column stores the index, not the string, to keep row size small.

**Sentinel value**:
A reserved id/index value meaning "not applicable" rather than a real
reference. Two conventions exist in the current schema: `venue_id == -999`
(`INFECTION_SEED_VENUE_ID` in the engine) marks a seed infection with no
venue; `255` marks an unset uint8 registry index. Two distinct engine
constants share that `255` value — `kDefaultEncounterTypeId` (`encounter_type_id`,
meaning "ordinary non-coordinated encounter", not an error) and `kNoSymptomId`
(`infector_symptom_id`, meaning "no infector to read a symptom from" —
seed/fomite/compartmental infections). `decode/sentinels.py` keeps them as
separate named constants (`DEFAULT_ENCOUNTER_TYPE_ID`, `NO_SYMPTOM_ID`)
mirroring `include/core/types.h`, not one shared value, since a future engine
change could split their values apart. `old_symptom_id`/`new_symptom_id`
carry no sentinel — always a real logged transition. Every `255` registry
sentinel must be excluded before a registry decode (see **Registry decode**),
not treated as data. `-999` is different: `lookups/venues` carries a real row
for `venue_id == -999` (`type == "infection_seed"`), so joining on it is not
an error case — it resolves to a meaningful label, not `NaN`.

**Registry decode**:
Translating a registry-indexed column (e.g. `encounter_type_id`) into its
string values via the matching `metadata/registries/*` list, mapping the
sentinel to a fixed placeholder (e.g. `"unknown"`). Handled generically by
`decode/`, parameterised by column, registry, and sentinel — not one function
per field. Distinct from **No-match value**: a `255` in the file is a
recorded fact ("this event has no venue-side encounter type"); a `NaN` from
an upstream **Enrichment** means the join found nothing to report — the file
never asserted anything. `decode_registry_column` labels these separately
(`unset_label` vs `no_match_label`) rather than collapsing them. A third
case — an index present in the column but absent from the loaded registry —
is neither: it is not a **Sentinel value** (the engine did assert something)
and not a **No-match value** (no **Enrichment** was involved). It signals a
genuinely incomplete registry, e.g. a multi-domain merge that only copied
`metadata/` from one input file. `decode_registry_column` raises rather than
silently labelling it, naming the id_column/registry/offending index.

**No-match value**:
A `NaN` in a registry-indexed column caused by an upstream **Enrichment**
(e.g. `enrich_with_state_at_time` when an id had no prior state logged
before `time`) rather than a value written by the engine. Not a **Sentinel
value** — the file itself never stores `NaN` in these columns; it only
appears after a join.

**Enrichment**:
Joining a raw event table against one or more lookup tables (by `person_id`,
`venue_id`, etc.) to attach metadata columns, or against another event table
via a time-aware join (e.g. an infector's symptom state at the moment of
transmission). Produced by `enrich/`; never mutates its inputs.

**Introspection**:
Reading a file's structure (which groups/datasets exist, row counts, byte
sizes, registry contents) without loading any event or lookup data. Exists
because different `simulation_events.h5` files vary — different registries,
different venue types, some missing `coordinated_encounters` or
`people_properties` entirely. Lives in `inspect/`, separate from `io/`
because it answers "what's in this file" rather than "give me this file's
data".

**Seam**:
The one file per subfolder (`io/__init__.py`, `decode/__init__.py`, etc.)
that all calls from outside the subfolder must go through. Internal files
within a subfolder are free to be restructured as long as the seam's
signatures hold.

**Enriched event table**:
The output of `load_enriched_events()`: one **Raw event table** with its
known **Registry** columns decoded and `people`/`venues` **Lookup tables**
joined on, in one call. Lives at top level (`june_events/load_enriched.py`),
not inside `io/`, `decode/`, or `enrich/`, because it orchestrates all three
rather than belonging to one. Distinct from a manually chained
**Enrichment** in that its registry-column mapping and join flags are fixed
defaults the caller can opt out of, not steps the caller assembles.

## Relationships

- A **Raw event table** has zero or more columns that are **Registry**
  indices, decodable via `decode/` against the matching **Registry**.
- An **Enrichment** takes a **Raw event table** plus one or more
  **Lookup tables** (or another enriched event table, for the time-aware
  case) and returns a new, wider DataFrame — the inputs are untouched.
- **Introspection** never touches **Raw event table** or **Lookup table**
  contents — dataset shape/dtype only, read directly from HDF5 metadata.

## Example dialogue

> **Dev:** "Why isn't `load_infections()` returning venue names directly?"
> **Maintainer:** "`io/` only ever returns a **Raw event table** — undecoded,
> unjoined. You call `enrich_with_venues(raw_infections, venues_lookup)`
> yourself. Keeps memory bounded and lets you skip the join entirely for
> files where you don't need it."
>
> **Dev:** "What about `coordinated_encounters` — where's the loader?"
> **Maintainer:** "There isn't one yet. It's the largest **Event type** by
> far (up to ~22GB on the biggest runs) and needs a query interface, not a
> DataFrame loader — that's a later phase, not part of this library yet."

## Flagged ambiguities

- `-999` and `255` are both "sentinel values" but are not interchangeable:
  `-999` is a real `venue_id`'s reserved value (int32), `255` is a reserved
  *registry index* (uint8). A generic decode/sentinel helper must take the
  sentinel as a parameter rather than hardcoding one.
- The `tests/fixtures/simulation_events_fixture.h5` real-data fixture predates
  the engine's uint8_t narrowing of `infector_symptom_id` (#16): it still
  stores `infector_symptom_id` as `uint16` and `0` (not `255`) for no-infector
  rows. Rerun `tests/fixtures/build_fixture.py` against a current
  `simulation_events.h5` before trusting fixture-based assertions about
  `infector_symptom_id` decode.
- `coordinated_encounters` is deliberately out of scope for this library's
  current phase (see ADR [0004](../../docs/adr/0004-size-checked-chunking-for-raw-table-loads.md)
  for why the raw-table loader is chunked at all) — do not add a bare
  `load_coordinated_encounters()` without first designing the query
  interface it actually needs.
- The engine also writes an `activities` **Registry** and an `activity_index`
  column, but that column lives on `/lookups/person_activities`
  (`PersonActivityRecord`, `event_types.h`), not on any **Raw event table** —
  this library has no loader for that lookup table yet, so `activities` has
  no decode path. Same "not built yet" status as `coordinated_encounters`,
  not a gap in `load_enriched_events`.
