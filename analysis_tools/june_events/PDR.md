# june_events — PDR

Library of composable Python functions for reading `simulation_events.h5`.
Plotting/cosmetics are left to the caller; example notebooks are written
separately by the user. See [CONTEXT.md](./CONTEXT.md) for terminology.

## Goal

Give notebooks a small set of reliable building blocks — load a table, load
a lookup, decode a registry column, join event to lookup, inspect a file's
shape — that work across `simulation_events.h5` files that vary in size and
schema (registries, venue types, which optional tables exist).

## Scope (this phase)

- `io/`, `decode/`, `enrich/`, `inspect/` subpackages, each with a seam
  `__init__.py` that all outside callers go through.
- Event types: `infections`, `deaths`, `symptom_changes`.
- Lookup tables: `people`, `venues` (dedicated wrappers with auto byte-decode);
  any other `lookups/*` dataset (`people_properties/*`, `population_summary`,
  `population_networks/*`) reachable via the generic raw-table loader.
- Registry decode: generic, parameterised by column/registry/sentinel — not
  bespoke per field.
- Enrichment: id-based joins (event → people, event → venues) and one
  generalised time-aware join (state-at-event-time via `merge_asof`).
- Introspection: `inspect_file()` — structure, sizes, registries — with no
  data materialised.

## Out of scope (this phase)

- `coordinated_encounters` — any access at all, including metadata peek.
  Up to ~22GB / 691M rows on observed runs; needs a query interface, not a
  DataFrame loader. Revisit as its own phase.
- `world_state.h5` / geo enrichment (venue and person coordinates, names,
  levels). Separate concern from `simulation_events.h5`.
- Plotting.
- pip packaging (`pyproject.toml`) — plain importable package for now.
- Automated tests — deferred to implementation time.

## Package layout

`__init__.py` re-exports each subpackage's `__all__`, so callers can
`from june_events import load_raw_table, decode_registry_column, ...`
without knowing which seam a function lives under.

```
analysis_tools/june_events/
  __init__.py
  CONTEXT.md
  PDR.md
  load_enriched.py         # seam: load_enriched_events
  io/
    __init__.py        # seam: load_raw_table, load_people_lookup, load_venues_lookup
    raw_tables.py
    lookups.py
  decode/
    __init__.py         # seam: decode_registry_column, load_registry
    registries.py
    sentinels.py
  enrich/
    __init__.py          # seam: enrich_with_people, enrich_with_venues, enrich_with_state_at_time
    joins.py
    temporal.py
  inspect/
    __init__.py           # seam: inspect_file
    types.py               # FileSummary, DatasetSummary dataclasses
    scan.py
```

## Function sketches

`io/raw_tables.py`
```python
def load_raw_table(
    path: str,
    dataset_path: str,             # e.g. "events/infections"
    chunk_threshold_bytes: int = 500_000_000,
    chunk_rows: int = 5_000_000,
) -> Optional[pd.DataFrame]:
    """None + logged warning if dataset_path absent. Single read below
    threshold, row-chunked read+concat above it."""
```

`io/lookups.py`
```python
def load_people_lookup(path: str, include_properties: bool = True) -> Optional[pd.DataFrame]
def load_venues_lookup(path: str) -> Optional[pd.DataFrame]
    # both auto-decode dtype.kind == 'S' columns to str
```

`decode/sentinels.py`
```python
SEED_VENUE_ID = -999              # INFECTION_SEED_VENUE_ID in the engine
NO_INFECTOR_ID = -1               # kInvalidPersonId, as used for infector_id
DEFAULT_ENCOUNTER_TYPE_ID = 255   # kDefaultEncounterTypeId (encounter_type_id)
NO_SYMPTOM_ID = 255               # kNoSymptomId (infector_symptom_id)
```

`decode/registries.py`
```python
def load_registry(path: str, registry_name: str) -> Optional[list[str]]
def decode_registry_column(
    df: pd.DataFrame,
    id_column: str,
    registry: list[str],
    unset_value: int | None = None,
    unset_label: str = "unknown",
    no_match_label: str = "not_recorded",
) -> pd.Series
    # NaN (from an upstream Enrichment join, e.g. enrich_with_state_at_time
    # with no prior state) is distinct from unset_value: NaN -> no_match_label,
    # unset_value -> unset_label
```

`enrich/joins.py`
```python
def enrich_with_people(event_df, people_lookup, prefix="person_") -> pd.DataFrame
def enrich_with_venues(event_df, venues_lookup, prefix="venue_") -> pd.DataFrame
    # left join on event_df["person_id"] / event_df["venue_id"] only (no
    # id_column param this phase, e.g. no infector-side enrichment yet);
    # venues_lookup carries a real venue_id == SEED_VENUE_ID row
    # (type == "infection_seed"), so seed infections join naturally to a
    # meaningful label rather than needing NaN-by-design special-casing;
    # prefix applies to lookup's non-key columns only, not the id column itself
```

`enrich/temporal.py`
```python
def enrich_with_state_at_time(
    event_df: pd.DataFrame,
    state_change_df: pd.DataFrame,
    id_column: str,             # e.g. "infector_id"
    time_column: str = "time",
    state_column: str = "new_symptom_id",
) -> pd.DataFrame
    # backward merge_asof; generalises the infector-symptom-at-infection
    # pattern to any (id, time, state) triple. state_change_df's own person
    # identifier is assumed to always be "person_id" (per CONTEXT.md's
    # Lookup table convention) and is matched against event_df[id_column]
    # internally — no separate parameter for it. Rows where the id had no
    # state_change_df entry at or before time_column correctly get NaN
    # (e.g. the id's current state was never logged as an explicit
    # transition before this instant) rather than a fallback default.
```

`inspect/types.py`
```python
@dataclass
class DatasetSummary:
    path: str
    n_rows: int
    dtype: str
    nbytes: int

@dataclass
class FileSummary:
    path: str
    datasets: list[DatasetSummary]
    registries: dict[str, list[str]]
```

`inspect/scan.py`
```python
def inspect_file(path: str) -> FileSummary
```

`load_enriched.py`
```python
DEFAULT_REGISTRY_COLUMNS = {
    "encounter_type_id": "encounter_types",
    "infector_symptom_id": "symptoms",
    "old_symptom_id": "symptoms",
    "new_symptom_id": "symptoms",
}

def load_enriched_events(
    path: str,
    dataset_path: str,                 # e.g. "events/infections"
    *,
    with_people: bool = True,
    with_venues: bool = True,
    registry_columns: dict[str, str] | None = None,  # id_column -> registry_name
    people_prefix: str = "person_",
    venues_prefix: str = "venue_",
) -> Optional[pd.DataFrame]:
    """One-call load + decode + enrich for notebooks that don't need to
    care about file size (see io.load_raw_table's chunk_threshold_bytes
    for the underlying cutoff). Orchestrates io/decode/enrich; not a seam
    of any one of them.

    Order: load_raw_table -> decode each registry_columns entry present
    in the table (skipped with a warning if the file lacks that registry;
    decoded value written to a new column, id_column's "_id" suffix
    stripped, original id column untouched) -> enrich_with_people if
    with_people and "person_id" present -> enrich_with_venues if
    with_venues and "venue_id" present.

    registry_columns defaults to DEFAULT_REGISTRY_COLUMNS. Does not
    include enrich_with_state_at_time (needs a second event table the
    caller must choose and load themselves) or transmission_mode_index
    (no matching metadata/registries/* entry exists in any inspected
    file).

    None + logged warning if dataset_path absent, matching
    load_raw_table.
    """
```

## Data model notes (from inspecting 5 real runs, 2026-07-10)

- `events/{infections,deaths,symptom_changes}` stay well under the chunking
  threshold even on multi-GB runs (e.g. 14.7M-row `symptom_changes` ≈ 0.3GB).
- `events/coordinated_encounters` ranges 227M–691M rows (~7–22GB) across the
  five runs inspected — confirms it needs its own phase, not a DataFrame load.
- `lookups/people` and `lookups/venues` are fixed-schema; `people_properties/*`
  and `population_networks/*` are dynamic per config (rat-related properties
  only appear in plague configs) — the generic raw-table loader handles these
  without needing per-field code.
- Root-level HDF5 attrs are empty on every file inspected — there is no
  schema-version marker in the file itself. `inspect_file()` is the only way
  to discover what a given file actually contains.
- Sentinels confirmed against `include/core/types.h` /
  `include/utils/event_logging/event_types.h` (as of the "Consistent
  sentinels" rollout, #16): `INFECTION_SEED_VENUE_ID = -999` (int32
  `venue_id`); `255` (uint8 registry-index unset) covers two distinct engine
  constants — `kDefaultEncounterTypeId` (`encounter_type_id`) and
  `kNoSymptomId` (`infector_symptom_id`, narrowed uint16_t -> uint8_t in the
  same rollout, written directly by the engine for seed/fomite/compartmental
  infections that have no infector to read a symptom from).
  `old_symptom_id`/`new_symptom_id` are always a real logged transition, no
  sentinel, so `decode_registry_column` is never called with an unset value
  on those two columns.
- `chunk_threshold_bytes = 500_000_000` / `chunk_rows = 5_000_000` confirmed as
  final defaults (Phase 2) — kept as a safety net per ADR-0004 even though no
  in-scope table on any inspected run reaches the threshold.
- `include_properties=True` on `load_people_lookup` attaches each
  `lookups/people_properties/*` dataset as a column by row position, not by an
  id join — confirmed against `src/utils/event_logging/event_writer_lookups.cpp`
  (`writePersonLookupTable`), which writes `lookups/people` and every
  `people_properties/*` dataset from the same ordered `world.people` iteration.

## Open questions

None outstanding. `enrich_with_state_at_time`'s generalisation was confirmed
against two call sites on real data: infector-symptom-at-infection
(infections × symptom_changes, id_column="infector_id", >99% match against
the engine's own infector_symptom_id where a prior state was logged; the
remainder are correctly NaN — the infector's current state was never logged
as an explicit transition before that instant) and deaths × symptom_changes
(id_column="person_id", 100% coverage — every death is preceded by a logged
symptom transition).
