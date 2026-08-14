# JUNE 2 — User Guide

This guide is for users who want to **run** JUNE 2, without having to edit the code.
Everything you need to control lives in two directories:

```
configs/<scenario>/      — YAML + CSV that describe your simulation
data/                    — Population-scale CSVs (rates, contact data, geography)
worlds/                  — Pre-built population HDF5 files (the "world")
```

---

## 1. Running a simulation

### 1.1 Build once

```bash
./rebuild.sh                # local build
```

The executable lands at `build/disease_sim`.

### 1.2 Pick a world and a config

A run needs **one world file** and **one `simulation.yaml`** (which pulls in
all the other YAML files via `config_paths:`).

```bash
# Serial
./build/disease_sim worlds/world_2021.h5 \
    --config configs/config_2021/simulation.yaml

# Parallel (MPI domain decomposition)
mpirun -n 4 ./build/disease_sim worlds/world_2021.h5 \
    --config configs/config_2021/simulation.yaml
```

> Output: `simulation_events.h5` (HDF5 event log) + the contents of
> `output/` (partition map, stats CSVs, etc.).

### 1.3 The "configs only" rule

YAML and CSV are loaded at runtime. **You do not have to rebuild after editing them.**
If you change anything in `src/` or `include/`, you must `./rebuild.sh`.

---

## 2. The shape of a config directory

A typical scenario folder (e.g. `configs/config_2021/`):

```
simulation.yaml               ← entry point; references everything else
disease.yaml                  ← what the pathogen does
infection_seeds.yaml          ← who gets infected on day 0 (and later seedings)
contact_matrices.yaml         ← per-venue contact rates and β (transmission rate)
schedules.yaml                ← weekly time-use templates
activity_preferences.yaml     ← which leisure venue each person prefers
coordinated_encounters.yaml   ← coordinated multi-person encounter recipes
policies.yaml                 ← lockdowns, isolation, behavioural overrides
vaccines.yaml                 ← vaccine products and rollout campaigns
parallel.yaml                 ← MPI domain decomposition
performance.yaml              ← schedule pre-compute toggles
bulk_seeds.csv                ← optional: many seeding events in one file
```

Only `simulation.yaml` is required by name; everything else is referenced
through `simulation.yaml`'s `config_paths:` block, so you can rename or
relocate sub-files freely as long as the path matches.

---

## 3. `simulation.yaml` — the entry point

Minimal example:

```yaml
time:
    start_date: "2024-01-01"
    end_date:   "2024-01-31"

config_paths:
    disease_file:               "configs/config_2021/disease.yaml"
    contact_matrices_file:      "configs/config_2021/contact_matrices.yaml"
    schedules_file:             "configs/config_2021/schedules.yaml"
    vaccines_file:              "configs/config_2021/vaccines.yaml"
    activity_preferences_file:  "configs/config_2021/activity_preferences.yaml"
    coordinated_encounters_file:"configs/config_2021/coordinated_encounters.yaml"
    performance_file:           "configs/config_2021/performance.yaml"
    parallel_file:              "configs/config_2021/parallel.yaml"
    policies_file:              "configs/config_2021/policies.yaml"
    infection_seeds_file:       "configs/config_2021/infection_seeds.yaml"

random_seed: 12345

output:
    stats_interval_days:   1
    flush_interval_days:   7
    max_event_buffer_size: 100000
    save_coordinated_encounters: false

regional_risk:
    enabled: false
    regional_risk_file: "data/regional_risk.csv"
```

| Field | What it controls |
|---|---|
| `time.start_date` / `end_date` | First and last simulated day. Must be ISO `YYYY-MM-DD`. |
| `random_seed` | Master RNG seed. Same seed + same world ⇒ bit-identical output (across MPI rank counts too). |
| `output.stats_interval_days` | How often summary stats are written. |
| `output.flush_interval_days` | How often the event buffer is flushed to HDF5. |
| `output.save_coordinated_encounters` | Write coordinated encounters to `/events/coordinated_encounters`. Default `false` — it's the largest dataset in `simulation_events.h5` and rarely needed for analysis; encounters still happen and still drive transmission either way, only the H5 record is skipped. |
| `regional_risk` | Optional per-geo multipliers on transmission and severity. |

> **Paths** are resolved relative to the working directory you launched
> from (typically the repo root). All sample configs assume that.

---

## 4. The disease — `disease.yaml`

This file defines **the pathogen's biology and clinical course**. It has
four blocks: `settings`, `symptom_tags`, `outcome_rates_csv`, `trajectories`,
`transmission`, and `immunity`.

### 4.1 Symptom tags

Every state a person can be in. Each tag has a numeric `value` used to
order them (higher = sicker). Stages with the same name across trajectories
share a tag.

```yaml
symptom_tags:
  - {name: healthy,        value: -1}
  - {name: exposed,        value:  0}
  - {name: asymptomatic,   value:  1}
  - {name: mild,           value:  2}
  - {name: severe,         value:  3}
  - {name: hospitalised,   value:  4}
  - {name: intensive_care, value:  5}
  - {name: dead_home,      value:  6}
  - {name: dead_hospital,  value:  7}
  - {name: recovered,      value: -2}
```

`settings` then groups these into roles the simulator needs to recognise:
which tags are terminal (`fatality_stage`), which mean recovery, which
require a hospital bed (`hospitalised_stage` / `intensive_care_stage`).
Add a new disease state by adding a tag here AND wiring it into the
relevant settings list AND using it in a trajectory.

### 4.2 Outcome rates CSV

A per-(age, sex, population) probability table that decides which
trajectory a newly infected person follows. Defined in
`data/infection_outcome_rates_<disease>.csv`:

```csv
filter.age,filter.sex,asymptomatic,mild,severe,hospital,icu,home_ifr,hospital_ifr,icu_ifr
0-4,male,0.456,0.533,0.0100,0.0009,9.1e-05,1.4e-05,0.0,0.0
...
```

- `filter.*` columns identify the bin (age band, sex, population type).
- The remaining columns are probabilities of each outcome for that bin
  (must sum to ≤ 1.0; the residual goes to `severe` if you declared it
  in `calculated_outcomes`).
- The `outcome_rates_csv:` block in `disease.yaml` tells the loader how
  to interpret these columns:
  - `populations:` lets you define population-typed sub-tables (e.g.
    `gp` general, `ch` care-home) using a `selection:` rule on person
    attributes. The CSV columns must then be prefixed `gp_*`, `ch_*`.
  - `outcomes:` maps each CSV column to the symptom tag it triggers.
  - `calculated_outcomes:` defines outcomes whose probability is
    `1 - sum(others)` rather than a CSV column.

### 4.3 Trajectories

A trajectory is the **clinical path** a person walks once their outcome
is sampled. Each has a `selection_key` (matching an outcome in the CSV),
a `severity` (0–1, used by vaccines to scale benefit), and an ordered
list of `stages`. Each stage has a `symptom_tag` and a `completion_time`
distribution (how long until the next stage).

```yaml
- description: "exposed => mild => recovered"
  selection_key: mild
  severity: 0.2
  stages:
    - symptom_tag: exposed
      completion_time: {type: beta, a: 2.29, b: 19.05, loc: 0.39, scale: 39.8}
    - symptom_tag: mild
      completion_time: {type: constant, value: 20.0}
    - symptom_tag: recovered
      completion_time: {type: constant, value: 0.0}
```

Distribution `type` can be `constant`, `normal`, `lognormal`, `beta`,
`gamma`, `exponweib`. Parameters match SciPy conventions.

To **add a new clinical pathway**, add a row in the outcome-rates CSV
(or a calculated outcome), then add a corresponding trajectory whose
`selection_key` matches.

### 4.4 Transmission

Two modes are supported, chosen by the top-level `mode:` key:

- **Trajectory-Driven** (default). A single global infectiousness
  curve (gamma PDF, etc.) parameterised by `max_infectiousness`,
  `shape`, `rate`, `shift`. Each symptom stage's relative
  infectiousness comes from `infectiousness_factor` on the trajectory.
- **Stage-Driven**. Each symptom tag has its own curve in
  `stage_curves:` and the clock resets when the person enters that
  stage. Useful for diseases where infectiousness changes character
  (e.g. plague's bubonic→pneumonic shift).

Within each mode you list `modes:` (the **transmission modes** —
`respiratory`, `physical_contact`, `fomite`, etc.). Each gets its own
`mode_transmissibility_multiplier`. **The per-mode `beta` is not here** — it's
in `contact_matrices.yaml`, per venue.

### 4.5 Immunity

```yaml
immunity:
  level: 0.95
  waning_rate: 0.001
```

Probability of becoming immune on recovery and the rate at which immunity
decays per day.

---

## 5. Contact matrices — `contact_matrices.yaml`

This file says **how often people meet at each venue** and **how
transmissive that venue is**, per transmission mode. It is the single
biggest knob for outbreak intensity.

```yaml
contact_matrices:
  household:
    modes:
      respiratory:
        beta: 0.0673
        bins: &b [Kids, Young Adults, Adults, Old Adults]
        contacts:
          - [1.63, 0.78, 0.78, 0.78]
          - [0.55, 0.45, 0.45, 0.45]
          - [0.39, 0.44, 0.44, 0.44]
          - [0.35, 0.45, 0.45, 0.71]
      physical_contact:
        beta: 0.134
        bins: *b
        contacts:
          - [0.84, 0.52, 0.52, 0.52]
          - [...]
```

Each venue declares one block per transmission mode:

- `beta` — venue-and-mode transmission rate (probability per
  contact-second; treat it as a tunable knob).
- `bins` — labelled subgroups within the venue. Both modes must use
  the same bins.
- `contacts[i][j]` — daily contacts an `i`-bin person has with a
  `j`-bin person.

The simulator multiplies these: `effective[i][j] = beta * contacts[i][j]`.

> **Don't double-count.** The `mode_transmissibility_multiplier` in
> `disease.yaml` is folded into per-venue `beta` here. Leave it at 1.0
> in `disease.yaml` unless you know what you're doing.

A `default_contacts_matrix:` block is used as the fallback for any
venue not explicitly listed.

---

## 6. Infection seeding — `infection_seeds.yaml` + `bulk_seeds.csv`

Seeds are the **initial spark**. JUNE 2 supports four types:

```yaml
global_parameters:
  default_seed_strength: 1.0
  base_cases_per_capita: 0.0002

infection_seeds:
  # 1. Uniform — random across the whole population
  - name: "basic"
    type: "uniform"
    date: "2024-01-01 08:00"
    parameters:
      cases_per_capita_multiplier: 15.0
      attribute_filters:                # optional restriction
        properties.ethnicity: "W"

  # 2. Exact — N cases in named geographic units, by age group
  - name: "targeted"
    type: "exact"
    date: "2024-01-02 12:00"
    geo_level: "MGU"
    parameters:
      age_groups: ["0-17", "18-64", "65+"]
      units:
        "E02004292": [5, 10, 2]   # 5 kids, 10 adults, 2 seniors
        "E02001234": 20           # 20 total

  # 3. Clustered — like exact, but cases cluster within households
  - name: "household_burst"
    type: "clustered"
    date: "2024-01-03 18:00"
    geo_level: "SGU"
    parameters:
      age_groups: ["18-30"]
      units: {"E01000001": 5}
```

Or for many events at once, point to a CSV:

```yaml
bulk_csv: "configs/config_2021/bulk_seeds.csv"
```

`bulk_seeds.csv` columns:

```
name,date,type,geo_level,geo_unit,trajectory_key,start_symptom,cases,cases_per_capita,
filter.properties.ethnicity,filter.sex,filter.age_groups
```

`trajectory_key` lets you force a specific clinical path on seeded
cases (e.g. seed plague cases directly into the `pneumonic_perished`
trajectory). `start_symptom` lets you skip the latent period.

> For **outbreak modelling**, prefer `exact` or `clustered` seeding in
> high-density geographic units over `uniform-with-filter` — a
> uniform seed wastes most cases on people the filter rejects.

---

## 7. Schedules — `schedules.yaml` (and optional `schedules.csv`)

A schedule is a **weekly time-use template**. Every person is assigned
exactly one schedule type at world load.

```yaml
day_type_cycle: [weekday, weekday, weekday, weekday, weekday, weekend, weekend]
default_schedule_type: "has_primary_activity"

schedule_types:
  has_primary_activity:
    priority: 0
    selection: []                       # empty = matches everyone (fallback)
    weekday:
      - {name: morning,           start: "08:00", end: "09:00", allowed_activities: [residence]}
      - {name: primary_activity,  start: "09:00", end: "17:00", allowed_activities: [primary_activity, residence]}
      - {name: evening,           start: "18:00", end: "21:00", allowed_activities: [leisure, residence],
         coordinated_only_activities: [sex]}
      - {name: night,             start: "21:00", end: "08:00", allowed_activities: [residence]}
    weekend: [...]
    participation:
      primary_activity: {weekday: 1.0, weekend: 0.0}
      leisure:          {weekday: 0.3, weekend: 0.5}
      residence:        {weekday: 1.0, weekend: 1.0}
```

- `priority` + `selection` — first matching schedule wins. Use
  `selection: []` for the catch-all fallback.
- Each time slot lists `allowed_activities` (what the person *may*
  do) and optionally `coordinated_only_activities` (slots that exist
  only for coordinated encounters — e.g. `sex` doesn't roll on its
  own; it's only realised when a partner pair commits).
- `participation` is the per-day probability the person engages in
  that activity at all. `0.3` weekday leisure ⇒ 30% chance to leave
  the house for leisure on a weekday.

For complex worlds (e.g. medieval) you can drive schedule assignment
from a `schedules.csv` with `filter.*` columns and `schedule.<name>`
probability columns. See `configs/config_plague/schedules.csv`.

### Activity preferences — `activity_preferences.yaml`

Once a person is going to do `leisure`, *which* venue category? This
file holds priority-ranked profiles with `selection:` criteria and
relative `weights:` per venue type:

```yaml
profiles:
  - name: "Young Males (Leisure)"
    activity: "leisure"
    priority: 10
    selection:
      - {property: age, operator: "<", value: 30}
      - {property: sex, operator: "==", value: "M"}
    weights: {gym: 5.0, pub: 3.0, cinema: 2.0, grocery: 0.5}
```

The first profile whose `selection` matches wins; weights are
normalised to a categorical distribution.

---

## 8. Networks — loaded from the world

Persistent social networks (e.g. `friendships`) live in the **world
HDF5 file** built by [MAY](https://github.com/mtcorread/MAY). They
are static across a run: every person's network ties are read at
startup and accessed via `WorldState::getNetworkPartners(person,
network_name)`. There are no formation or dissolution mechanics in
the public engine.

To use a network in a scenario, reference it by name from a
coordinated encounter (see §9):

```yaml
- name: "social_encounters"
  network: "friendships"        # must match a network type in the world HDF5
  ...
```

If the named network is absent from the world file, the encounter's
partner list will simply be empty.

---

## 9. Coordinated encounters — `coordinated_encounters.yaml`

These are activities where **two or more specific people must commit
to meet** (sex, group encounters, etc.) — distinct from the implicit
mass-action contacts that contact matrices model.

```yaml
coordinated_encounters:
  enabled: true
  log_commitments: true

  encounters:
    - name: "social_encounters"
      priority: 10
      enabled: true
      network: "friendships"             # which network supplies partners
      trigger_slots: ["leisure"]         # which schedule slot triggers it
      allowed_venues: ["pub", "cinema"]  # real venues for the meet-up
      daily_max_distribution:
        type: poisson
        mean: 1.5
      proposal_probability: 0.05
      invite_distribution: {type: binomial, p: 0.3}
      acceptance_probability: 0.8
      is_virtual: false
```

Key fields:

| Field | Purpose |
|---|---|
| `network` | The network in the world HDF5 whose ties supply partner candidates. |
| `trigger_slots` | Schedule slot names that activate this encounter. |
| `proposal_probability` / `acceptance_probability` | Multipliers on the offer/accept stages. |
| `allowed_venues` | Real venue names. Empty + `is_virtual: true` ⇒ encounter happens at a synthetic venue. |
| `virtual_contact_matrix` | When `is_virtual: true`, which key in `contact_matrices.yaml` supplies the contact rates and β. |

---

## 10. Policies — `policies.yaml`

Two kinds of policies, evaluated in order:

```yaml
policies:
  symptom_policies:                    # behavioural responses to illness
    - name: "isolate_when_severe"
      symptoms: ["severe"]
      override_activities: ["primary_activity", "leisure"]
      replacement: "residence"
      compliance_rate: 0.8

    - name: "hospitalize"
      symptoms: ["hospitalised", "intensive_care"]
      override_activities: "*"
      replacement: "medical_facility"
      compliance_rate: 1.0

  temporal_policies:                   # date-bounded restrictions
    - name: "first_lockdown"
      start_date: "2024-03-23"
      end_date:   "2024-05-15"
      override_activities: ["primary_activity", "leisure"]
      replacement: "residence"
      compliance_rate: 0.85
      applies_to:                      # optional selection
        - {property: age, operator: "<", value: 65}
```

- Symptom policies fire when a person has any of the listed symptom
  tags. They override the activities the schedule asked for.
- Temporal policies are active in date windows (inclusive). If both
  apply, the **first matching** policy wins — order them by
  specificity.
- `override_activities: "*"` overrides everything.
- `compliance_rate` is rolled per slot, per person.
- `applies_to` uses the same selection grammar as everywhere else
  (`property`, `operator`, `value`).
- Use `follow_up_policy` to chain policies (e.g. mild ⇒ severe).

---

## 11. Vaccines — `vaccines.yaml`

```yaml
enabled: true

vaccines:
  Pfizer:
    doses:
      - days_to_effective: 5
        days_to_waning:    6
        days_to_finished:  7
        waning_factor:     1.0
        infection_efficacy:
          covid19: {"0-100": 0.52}
        symptom_efficacy:
          covid19: {"0-100": 0.52}
      - {... dose 2 ...}

campaigns:
  "Care Home Pfizer":
    start_date:        "2024-01-01"
    end_date:          "2024-01-20"
    days_to_next_dose: [0, 14]
    daily_coverage:    0.3
    vaccine_type:      "Pfizer"
    dose_sequence:     [0, 1]
    selection:
      - {property: "activities.residence.venue_type", operator: "==", value: "care_home"}
```

- A **vaccine** is a sequence of doses. Each dose has its own ramp-up
  (`days_to_effective`), waning window, and per-disease, per-age-band
  efficacy against infection and symptoms. Efficacy multiplies the
  trajectory's `severity` to compute symptom protection.
- A **campaign** rolls out a vaccine to a selected sub-population at
  `daily_coverage` until end-date. `last_dose_type_filter` lets a
  booster campaign target only people whose previous dose was a
  specific brand.

---

## 12. Parallel — `parallel.yaml`

If you run with MPI, this file controls domain decomposition.

```yaml
parallel:
  enabled: true
  partitioning:
    level: "MGU"                          # geography unit to partition on
    centroids_file: "data/domain_decomposition/2021/mgu_centroids.csv"
    adjacency_file: "data/domain_decomposition/2021/mgu_adjacency_graph.json"
    metis: {imbalance_tolerance: 0.05}
  chunked_loading:
    person_metadata_chunk_size: 100000
    geo_unit_chunk_size:        1000
  communication:
    buffer_size_mb: 256
  output:
    save_partition: true
    partition_file: "output/partition.json"
```

- `level:` must match the granularity of your centroid/adjacency files
  (`SGU` ⊂ `MGU` ⊂ `LGU` ⊂ `XLGU`). For full-England runs use `MGU`;
  for the 1911 Durham world there is currently only one `MGU`, so
  partitioning is done at `SGU`.
- Output is **bit-identical across rank counts** (with the same
  random seed). If `mpirun -n 1` and `mpirun -n 4` disagree, that's a
  bug — file an issue.

---

## 13. Performance — `performance.yaml`

```yaml
performance:
  precompute_schedules: true
  deterministic_activities: ["residence", "primary_activity"]
  track_active_infections_only: true
```

- `deterministic_activities` are pre-computed once at world load and
  cached. `leisure` etc. stay stochastic at runtime.
- `track_active_infections_only` skips state updates for healthy
  people. Default on; turn off only if you need full event traces.

---

## 14. The `data/` directory

Reference data shared across scenarios. Edit (or replace) these CSVs
to retune the model without touching code.

| File | Purpose | Schema |
|---|---|---|
| `infection_outcome_rates_<disease>.csv` | Per (age, sex, pop) outcome probabilities | `filter.age,filter.sex,asymptomatic,mild,severe,hospital,icu,home_ifr,hospital_ifr,icu_ifr` |
| `bulk_seeds_<disease>.csv` | Bulk seeding events | see §6 |
| `regional_risk.csv` | Per-geo transmission/severity multipliers | `geo_unit,transmission_factor,severity_factor` |
| `domain_decomposition/<world>/` | METIS centroids + adjacency | CSV + JSON |

Convention: every column starting `filter.` is a row-selection key
(equality match); the rest are values.

---

## 15. The `worlds/` directory

Pre-built HDF5 population files. **You don't edit these directly** —
they're produced by [MAY](https://github.com/mtcorread/MAY), the
companion world-builder repository.

| World | Use |
|---|---|
| `worlds/world_2021.h5` | Full England 2021 census (~60M, partitioned at MGU) |
| `worlds/medieval.h5`   | Medieval plague world |
| `worlds/1911_durham_world.h5` | 1911 Durham (1 MGU + 268 SGU parishes) |

Pick a world that matches your scenario's geography. The world
defines the population, households, workplaces, schools, geography
hierarchy, and pre-existing networks (e.g. friendships). The
configs tell that population *what to do*.

---

## 16. Writing a new scenario from scratch

The fastest path:

1. **Copy** an existing scenario folder closest to your use case
   (`cp -r configs/config_2021 configs/config_myrun`).
2. **Edit `simulation.yaml`** — dates, output cadence, and update
   the `config_paths:` to point inside `configs/config_myrun/`.
3. **Pick a disease** — either edit `disease.yaml` and the matching
   `data/infection_outcome_rates_<disease>.csv`, or swap in plague.
   Make sure `disease.name` matches the CSV file you're using.
4. **Tune contacts** — `contact_matrices.yaml` is your single
   biggest knob on outbreak intensity. Start there before touching
   anything else.
5. **Choose seeds** — start with one `uniform` seed; switch to
   `exact`/`clustered` once you want geographic targeting.
6. **(Optional) policies, vaccines, networks** — add as needed.

---

## 17. Checkpoint / restart

A run can be snapshotted at configurable points and later resumed. The
continuation is **bit-identical to the uninterrupted run, regardless of
the MPI rank count** used for either the checkpoint or the resume
(set-identical events; the only cross-rank difference is the documented
rank-stamped `coordinated_encounters.group_id`).

### 17.1 Enabling checkpoints

Add a `checkpoint:` block to `simulation.yaml`:

```yaml
checkpoint:
  enabled: true
  output_dir: checkpoints/   # resolved under the run directory
  every_n_days: 30           # interval cadence
  on_dates: null             # OR an explicit list of "YYYY-MM-DD"
  keep_last: 3               # 0 = keep all; otherwise rotate oldest
```

**Cadence is mutually exclusive.** If `on_dates` is a non-null,
non-empty list it takes precedence and `every_n_days` is **ignored**
(a one-line notice is printed at startup). Use `null` to mark a field
absent. With neither set, no checkpoints are written.

- `every_n_days: N` → checkpoint at the end of every `N`-th completed
  day (`(day+1) % N == 0`).
- `on_dates: ["2026-03-01", "2026-06-15"]` → checkpoint at the end of
  exactly those calendar days. Dates outside the simulation window are
  warned about at startup and never fire.

Checkpoint frequency has **zero effect on results** — it is purely an
operational knob (the snapshot is read-only and the RNG is stateless).

### 17.2 What a checkpoint looks like

```
runs/<run-id>/checkpoints/
  checkpoint_20260301_day045/
    delta_rank0.h5 …          # per-rank mutable-state shard
    shard_index.yaml          # geo_unit -> (shard, offset, count)
    state.h5                  # global scalars, applied seeds, …
    manifest.yaml             # written LAST = commit marker
  latest -> checkpoint_20260301_day045
```

Writes are atomic: a checkpoint is built in `*.tmp/` and renamed only
once complete, and `manifest.yaml` is written last — a crash mid-write
never leaves a selectable, corrupt checkpoint. The 25 GB world is **not**
copied; the checkpoint stores only the small mutable delta and is
restored on top of the normal (fast) world load.

### 17.3 Resuming

```bash
# resume the newest checkpoint of a previous run (any rank count)
mpirun -np 8 ./build/disease_sim \
  --config configs/config_2021/simulation.yaml \
  --world worlds/world_2021.h5 \
  --restart-from runs/20260301-120000/checkpoints/latest
```

`--restart-from` accepts a checkpoint directory or the `latest`
symlink. The resume is a **new run directory** (the original is left
untouched, so the two can be diffed). The checkpoint's recorded
effective seed is authoritative; passing a conflicting `--seed` is
rejected (no silent override). Start-of-simulation infection seeding is
not re-applied on resume.

**Run length on resume.** You may combine `--restart-from` with
`--days` (or change `end_date`) to extend/shorten the continued run,
**but `--days N`/`end_date` is anchored to the original `start_date`,
not to the checkpoint.** It sets the absolute end (day index from
start); the resumed run fills `[completed_day + 1, total_days)`. To run
**M more days** after a checkpoint that completed at day **D**, pass
`--days (D + 1 + M)`. Requesting an end at or before the resume day is
**rejected with a descriptive error** (not a silent no-op).

### 17.4 Limitations

- **Compartmental-model scenarios are not supported.** If a
  `compartmental_model_sidecar` is configured together with
  checkpointing, the run aborts with a clear error — the external ODE
  plugin's internal state is opaque to the engine. (Plugin
  serialize/deserialize hooks are a planned follow-up.)
- Checkpoints are valid only against the exact world they were taken
  from.

### 17.5 Verifying determinism

`tests/checkpoint_determinism_check.py` drives the binary through a
`(write_ranks -> resume_ranks)` matrix and asserts set-identical events
versus an uninterrupted baseline. It is wired as the CTest
`test_checkpoint_determinism` (skips cleanly without the world file):

```bash
ctest -R test_checkpoint_determinism --output-on-failure
```

---

## 18. Quick troubleshooting

| Symptom | Most likely cause |
|---|---|
| Output differs across rank counts | Bug or non-deterministic seeding — never accept partial determinism. |
| HDF5 read aborts with "file is locked" | Two `disease_sim` runs against the same world in parallel — don't. Run sequentially or copy the world. |
| Many "encounter_type=255" events | Normal: sentinel for non-coordinated venue infections (typically ~87%). |
| Low/zero infections despite seeds | Check `contact_matrices.yaml` betas, not seed counts. β is the dominant knob. |
| `--restart-from` aborts on seed mismatch | The checkpoint's effective seed is authoritative; drop the conflicting `--seed`. |
| Checkpointing aborts immediately | A `compartmental_model_sidecar` is set — checkpointing is unsupported for compartmental scenarios. |

---

*Anything you can't reach through these files is a missing config
hook, not a reason to patch C++. File an issue or extend the
loader; don't bake a value into the source.*