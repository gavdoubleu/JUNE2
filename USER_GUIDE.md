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
./build/disease_sim --world worlds/world_2021.h5 \
    --config configs/config_2021/simulation.yaml

# Parallel (MPI domain decomposition)
mpirun -n 4 ./build/disease_sim --world worlds/world_2021.h5 \
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
policies.yaml                 ← lockdowns, isolation, behavioural overrides, observed metrics
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

```

| Field | What it controls |
|---|---|
| `time.start_date` / `end_date` | First and last simulated day. Must be ISO `YYYY-MM-DD`. |
| `random_seed` | Master RNG seed. Same seed + same world ⇒ bit-identical output (across MPI rank counts too). |
| `output.stats_interval_days` | How often summary stats are written. |
| `output.flush_interval_days` | How often the event buffer is flushed to HDF5. |
| `output.save_coordinated_encounters` | Write coordinated encounters to `/events/coordinated_encounters`. Default `false` — it's the largest dataset in `simulation_events.h5` and rarely needed for analysis; encounters still happen and still drive transmission either way, only the H5 record is skipped. |

Every row of `/events/coordinated_encounters` carries a `group_id` that is
shared by the rows of one encounter. It is built from the encounter itself, so
it is the same at any rank count: the simulation day in bits 62 to 51 (up to
4095), the host's running number for that day in bits 50 to 36 (up to 32767),
and the host's person id in bits 35 to 0. The value is always positive as a
signed 64 bit integer. A day, a count or a person id that does not fit stops
the run with an error that names the value and the limit.

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
- The remaining columns are probabilities of each outcome for that bin.
  Each row must sum to 1.0 (within 0.01), or the loader stops with an
  error naming the row. Nothing is filled in for you, so there is no
  residual outcome.
- `outcome_rates_csv:` in `disease.yaml` is either a bare file path or
  a mapping with a `file:` key:

  ```yaml
  outcome_rates_csv: data/infection_outcome_rates_<disease>.csv
  # or
  outcome_rates_csv:
    file: data/infection_outcome_rates_<disease>.csv
  ```

  The path is relative to the directory of `disease.yaml`. A table with
  `filter.*` columns needs nothing else.
- `populations:`, `sexes:` and `outcomes:` inside the mapping form are
  only used to check a wide-layout table, one whose columns are named
  `<population>_<outcome>_<sex>`. When at least one such column is
  present, each population and sex group must sum to 1.0 on its own.
  Tables with `filter.*` columns do not need them.

### 4.3 Trajectories

A trajectory is the **clinical path** a person walks once their outcome
is sampled. Each has a `selection_key` (matching an outcome in the CSV),
a `severity` (a number that only ranks the trajectories; a vaccine uses it to
find the least severe one, see section 11), and an ordered list of `stages`.
Each stage has a `symptom_tag` and a `completion_time` distribution (how long
until the next stage).

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

**How a person is put in a bin.** The first rule that applies wins:

1. If the person's subset in the venue has the same name as a bin, that
   bin (for example `resident` or `worker`).
2. If the bins include `male` and `female`, the bin for the person's sex.
3. If the bins are age bands, the band holding the person's age in whole
   years. Bands are written `min-max`, `min+` or `[min,max]`, and ages
   above 99 count as 99.

Anyone who matches nothing goes to the first bin. For example, a care
home visitor in a `[resident, worker]` matrix is binned as a resident. A
matrix with a single bin puts everyone in it.

A `default_contacts_matrix:` block is the fallback for any venue type
(or mode) with no matrix of its own. Using it is refused at load, with
the affected pairs listed, unless you set `allow_default_matrix: true`
at the top level of `contact_matrices.yaml`. The setting only says you
accept the default for those venue types.

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
- A schedule type with no `selection` matches everyone at its priority.
  A type that should be reachable only from `schedules.csv` needs a
  selection that can never match, for example
  `{property: age, operator: "<", value: 0}`.

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

Enabling or disabling an encounter changes the random draws of every
encounter listed after it, because the draw key counts the enabled
encounters in order. Expect different (but still reproducible) results
for those encounters after such an edit.

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
- Temporal policies are active in half-open date windows (`start_date`
  included, `end_date` excluded).
- Policies are tried in order, and each makes its own compliance roll for
  each person. A person who refuses a policy is treated as if it did not
  apply to them, and the next policy that covers them rolls separately.
  So overlapping policies at 0.85 and 0.8 hold a person with probability
  1 - 0.15 x 0.2 = 0.97, not 0.85. To give a subgroup its own rate, make
  `applies_to` disjoint. The same holds for symptom policies.
- `override_activities: "*"` overrides everything.
- `compliance_rate` is a sticky per-person decision for each policy. The
  draw is keyed by the policy's kind and name, not its position, so
  reordering, adding or removing other policies leaves it unchanged;
  renaming the policy changes it.
- Names must be unique within each kind.
- Limits: at most 128 policies of each kind per run, and at most 32 of
  each kind may apply to any one person. A person matching more (for
  example, 33 temporal policies with no `applies_to`) stops the run at
  startup with an error naming the person, kind and count. Narrow
  `applies_to`, or merge policies, to stay under it.
- `applies_to` uses the same selection grammar as everywhere else
  (`property`, `operator`, `value`).
- Use `follow_up_policy` to chain policies (e.g. mild ⇒ severe).
- Either kind can carry a `trigger:` that ties it to an observed count, such as
  deaths or cases in a classroom, instead of (or as well as) a date. See
  [Observed metrics](#observed-metrics) below.

### Policies with a `replacement_schedule`

A policy normally swaps the activity a person was going to do for its
`replacement`. Setting `replacement_schedule` instead moves a person who is
already on a temporary schedule hop (a trip) onto the named schedule. It works
the same in symptom and temporal policies:

```yaml
temporal_policies:
  - name: "curtail_trips"
    start_date: "2020-03-23"
    end_date: "2020-05-11"
    override_activities: ["leisure"]
    replacement_schedule: "stay_put"     # a temporary schedule
    compliance_rate: 0.7
```

This is not intuitive, so note the rules:

- The policy acts only on people who are on a hop. Everyone else is left
  alone, even though `override_activities` names an activity they are doing and
  `replacement` would normally send them home. In the example, a person who
  stays in town for leisure is unaffected, and 70 percent of the people who set
  off on a trip are cut short.
- The person is moved at the first slot where the policy covers them. For that
  slot they are placed at the `replacement` activity (residence unless you say
  otherwise). From the next slot they run the replacement schedule from the
  position they had reached on the trip.
- For a temporal policy the swap is final. Nothing is restored when the window
  ends, and the trip finishes on the replacement schedule's timetable. A
  symptom policy instead pins the person where they are and puts the trip back
  when the symptoms end.
- The named schedule must exist in `schedules.yaml`. An unknown name stops the
  run at load time.

### Transmission effects

A policy can reference a CSV of transmission effects with
`transmission_effects_file`. A policy with that field and no location action
changes transmission only; it does not send anyone home. The path may be
relative to `policies.yaml`.

```yaml
policies:
  temporal_policies:
    - name: masks
      start_date: "2024-03-23"
      end_date: "2024-05-15"
      compliance_rate: 0.8
      transmission_effects_file: masks.csv
    - name: distancing
      start_date: "2024-03-23"
      end_date: "2024-05-15"
      compliance_rate: 1.0
      transmission_effects_file: distancing.csv
```

`masks.csv` can select people using the existing `filter.*` columns, including
`filter.properties.<name>`:

```csv
filter.properties.risk_group,scope,transmission_mode,effect_channel,multiplier
high,person,respiratory,source_infectiousness,0.60
high,person,respiratory,target_susceptibility,0.70
```

`distancing.csv` selects venues with the same column convention:

```csv
filter.venue_type,scope,transmission_mode,effect_channel,multiplier
school,venue,respiratory,contact_intensity,0.80
```

The other channels are `environmental_risk` (for example, cleaning on a
fomite mode) and the two person channels shown above. Each row is one effect;
all distinct matching rows and active policies multiply. A missing effect uses
`1.0`, and duplicate definitions are rejected. Venue effects apply once per
venue and require a temporal policy with `compliance_rate: 1`; symptom policies
can carry person effects. Source/target person effects share their policy's
sticky compliance decision with any location action.

### Observed metrics

A policy can follow the epidemic instead of the calendar. That takes two
steps: declare a count in an `observed_metrics:` block at the top level of
`policies.yaml` (next to `policies:`), then give a policy a `trigger:` that
names it.

```yaml
observed_metrics:
  deaths:               {source: death,     scope: world}
  icu_now:              {source: stage,     scope: world, stage: intensive_care}
  primary_venue_cases:  {source: infection, scope: venue,        activity: primary_activity}
  primary_parent_cases: {source: infection, scope: parent_venue, activity: primary_activity}
```

The name on the left is only a label. A policy refers to the metric by it, and
the engine reads nothing into the words, so choose a name that says what is
counted. These four count:

- `deaths`: every death in the whole population so far.
- `icu_now`: the number of people in the `intensive_care` stage right now, across
  the whole population.
- `primary_venue_cases`: for each venue separately, the infected people
  registered to that venue for `primary_activity`. A classroom's count comes from
  its pupils and any staff registered to it, an office's from its workers, a care
  home's from the people who have it as their primary activity, and so on.
  Counts are never mixed between venues.
- `primary_parent_cases`: the same people, pooled by parent venue. A school's
  count adds up its classrooms, a company's adds up its offices, and a venue with
  no parent counts on its own.

| Field | What it controls |
|---|---|
| `source` | What is counted. One of the four sources below. |
| `scope` | Where it is counted: `world`, `venue` or `parent_venue`. |
| `activity` | An activity name known to the world, such as `primary_activity`. Required for `venue` and `parent_venue`, and not allowed for `world`. |
| `stage` | A symptom tag from `disease.yaml`. Required by `stage` and `stage_admission`, not allowed with the other sources. |

| Source | Kind | Counts |
|---|---|---|
| `infection` | current | People infected right now, from the moment of exposure, so exposed people are included. A person leaves the count when they recover or die. |
| `stage` | current | Infected people whose current stage is the named one. |
| `death` | cumulative | Deaths so far. |
| `stage_admission` | cumulative | Entries into the named stage, including the stage a person is placed in when infected. A person moving from the stage to itself is not an entry, and `stage_admission` on `exposed` counts every new infection. |

A current metric goes up and down, and a cumulative one only rises. The source
decides which, and the trigger fields a metric needs depend on it.

`world` is a single total. `venue` keeps a separate count for each venue and
counts a person toward every venue the world registers to them for the
metric's `activity`. `parent_venue` counts them toward the parent of each of
those venues instead (the school for a classroom, the company for an office),
or toward the venue itself when it has no parent. A person counts at most once
per venue or parent venue, whatever the source, so someone registered to two
classrooms of one school counts once at that school. People who die keep their
registrations, so a death counts at the venues the person was registered to.

A metric cannot be limited to one venue type. It counts at every venue that
people are registered to for its activity, and only the policy's
`override_venue_types` narrows where the policy acts. A metric on
`primary_activity` therefore keeps a count for every office and care home too.
A policy that names only classrooms never uses those counts. Venue scope suits
activities with one home venue, such as school or work. For an activity like
`leisure`, where a person is registered to several venues, a venue count is a
weaker signal.

Counts are read once a day, after the last slot, and a change applies from the
first slot of the next day. A death is counted at the end of the slot it happens
in, so each day's total holds every death up to midnight. Counts are sums of whole
numbers across ranks, so results do not depend on the number of MPI processes.

A metric is only counted when a trigger reads it. Declared metrics that nothing
reads are still checked at startup (an unknown activity or stage is an error)
but cost nothing, and a run with no triggers does no counting at all.

### Triggering a policy from a metric

Add `trigger:` to a symptom or temporal policy. With a world scope metric, the
policy is in force when its date window holds and its trigger is closed.
Closing a trigger puts the policy in force and reopening it lifts the policy.
`start_date` and `end_date` are optional on a triggered policy, and without
them the trigger alone decides. (Venue scope metrics and responses work per
venue, as the next two sections describe.)

```yaml
observed_metrics:
  icu_now: {source: stage, scope: world, stage: intensive_care}
  deaths:  {source: death, scope: world}

policies:
  temporal_policies:
    - name: "stay_home_when_icu_is_busy"
      override_activities: ["primary_activity", "leisure"]
      replacement: "residence"
      compliance_rate: 0.8
      trigger:
        metric: icu_now          # a current metric
        activate_at: 50          # close at 50 or more
        deactivate_at: 20        # reopen at 20 or fewer

    - name: "stay_home_after_deaths"
      override_activities: ["primary_activity", "leisure"]
      replacement: "residence"
      trigger:
        metric: deaths           # a cumulative metric
        activate_at: 100         # close after 100 deaths
        calm_days: 14            # reopen after 14 days without a rise
        tolerance: 2             # a rise of 2 or fewer still counts as calm
```

`icu_now` is the number of people in intensive care right now, so it goes up and
down. `deaths` is the running total of deaths, so it only rises.

The fields you write depend on whether the metric is current or cumulative.
The loader rejects any other mix.

| Field | Current metric | Cumulative metric |
|---|---|---|
| `activate_at` | Required, a whole number of at least 1. Closes at this value or above. | Required, at least 1. Closes once the metric has grown by this much since the baseline. |
| `deactivate_at` | Required, from 0 up to below `activate_at`. Reopens at this value or below. | Not allowed. |
| `calm_days` | Not allowed. | Required, at least 1. |
| `tolerance` | Not allowed. | Optional whole number, default 0. |
| `tolerance_pct` | Not allowed. | Optional whole percent from 0 to 100, default 0. |

For a current metric the gap between `activate_at` and `deactivate_at` stops
the policy from flickering on and off around one value.

For a cumulative metric `activate_at` is growth and not a level. The baseline is
0 at the start of the run and the current value after each reopening. Once
closed, the trigger keeps a reference point: the value at closing, and later the
value each time the metric rose by more than the tolerance. When `calm_days`
days pass since the reference point without such a rise, the trigger reopens.
The tolerance is `tolerance` or `tolerance_pct` percent of the reference value,
whichever is larger, so a reference of 400 deaths with `tolerance_pct: 2`
tolerates 8 more. Rises are measured against the reference point and not
against yesterday, so a slow trickle that adds up past the tolerance still
moves the reference forward and restarts the clock.

For example, `activate_at: 100`, `calm_days: 200` and `tolerance: 5` closes at
100 deaths, reopens after 200 days with at most 5 further deaths, and closes
again after another 100.

A trigger starts open, so the policy starts out of force. The trigger is first
evaluated at the end of day 0, so a scenario that is already above
`activate_at` at the start closes from the second day. The trigger keeps
updating while the date window is closed, so a threshold crossed before
`start_date` takes effect when the window opens.

On a world scope metric the trigger switches the whole policy on and off. A
symptom policy still only affects people with its listed symptoms while it is in
force.

### Closing venues one at a time

When the metric has `scope: venue` or `scope: parent_venue`, the trigger keeps
its own state for every venue (or parent venue) it has seen, so one classroom
can be closed while the next stays open.

```yaml
observed_metrics:
  primary_venue_cases:  {source: infection, scope: venue,        activity: primary_activity}
  primary_parent_cases: {source: infection, scope: parent_venue, activity: primary_activity}

policies:
  temporal_policies:
    - name: "close_classroom"
      override_activities: ["primary_activity"]
      override_venue_types: ["classroom"]
      replacement: "residence"
      trigger:
        metric: primary_venue_cases    # infected people per venue
        activate_at: 2
        deactivate_at: 0

    - name: "close_school"
      override_activities: ["primary_activity"]
      override_venue_types: ["school", "classroom"]
      replacement: "residence"
      trigger:
        metric: primary_parent_cases   # infected people per parent venue
        activate_at: 5
        deactivate_at: 0
```

`close_classroom` closes a single classroom once 2 of the people registered to it
are infected, and reopens it when none are. `close_school` uses the count pooled
over a school's classrooms, so it closes every classroom of that school at 5 and
reopens them at 0. Both metrics also count offices and other venues that use
`primary_activity`. These policies leave those alone because
`override_venue_types` names only classrooms and schools.

`override_venue_types` limits a policy to venues of those types, using the venue
type names in your world. Without it (and without `exempt_venue_types`) the
policy acts on every venue the metric closes, and with no `override_activities`
it acts on every activity. Check that `override_activities` and
`override_venue_types` match where the metric counts, because a policy only acts
on the venues its metric closes.

A closed venue sends everyone scheduled there for a listed activity to the
replacement: pupils, staff and visitors, whatever household they come from. The
policy's `applies_to` and exemptions still apply. A closed venue also
hosts no coordinated encounter and takes in no follower mirroring a host. A
closed parent venue closes the parent and all of its children. A person who is
scheduled at no venue in a slot is not affected.

The loader refuses a venue triggered policy that:

- is a symptom policy (write a temporal policy instead);
- has a `compliance_rate` other than 1, since a closed venue closes for
  everyone;
- uses `transmission_effects_file`;
- lists `override_activities` without the metric's activity (an empty list or
  `"*"` means every activity).

Policies are checked in order and the first match wins, as everywhere else in
this file.

### Voluntary avoidance

A trigger can also model people choosing to stay away from a venue that is still
officially open. Instead of the threshold fields it takes a `response`, which
replaces `compliance_rate`.

```yaml
observed_metrics:
  primary_venue_cases: {source: infection, scope: venue, activity: primary_activity}

policies:
  temporal_policies:
    - name: "keep_children_home"
      override_activities: ["primary_activity"]
      override_venue_types: ["classroom"]
      replacement: "residence"
      applies_to:
        - {property: age, operator: "<", value: 18}
      trigger:
        response:
          base: 0.02
          terms:
            - {metric: primary_venue_cases, weight: 0.05}
        decision_unit: {shared_venue: residence}
```

Here `primary_venue_cases` counts the people who are infected right now, exposed
ones included, separately for each venue. A person counts toward the venue they
are registered to for `primary_activity`: for a classroom that is its pupils and
any staff registered to it, and for an office its workers. The policy only
acts on classrooms because of `override_venue_types`, so the offices' counts go
unused.

Each day the engine works out an avoidance level:

```
level = base + weight_1 * metric_1 + weight_2 * metric_2 + ...     (capped at 1)
```

Every decision unit has one fixed random draw between 0 and 1, set by the run
seed, the unit and the policy. A unit avoids exactly when the level is above its
draw. At a level of 0.3 about 30% of units avoid. If the level rises to 0.5 the
same 30% still avoid, plus more, and if it falls they return in the reverse
order. Nothing is re-rolled during the run, and the draws do not depend on the
number of ranks.

| Field | What it controls |
|---|---|
| `base` | Optional level before any count, from 0 to 1. Default 0. |
| `terms` | A non-empty list of `{metric, weight}`. Each weight is above 0, and a metric appears once. All terms must share one scope and one activity. |
| `decision_unit` | Required. `person` makes each person decide alone. `{shared_venue: <activity>}` makes everyone registered to the same venue for that activity decide together, for example `residence` so that people from one household avoid together. |

For a venue or parent venue metric the level is worked out per venue. Venues the
metric has not counted sit at `base`. Only people the metric counts at a venue
can avoid it, meaning those registered to it for the metric's activity. A
visitor, or someone scheduled at a venue they are not registered to, is not
reached, which differs from a closure. The policy's `override_activities` must
include the metric's activity.

#### Who a response reaches

A response acts on the person who is scheduled at the venue, in their own slot.
In the example above that is the child, who is sent to the replacement. A
parent's `primary_activity` slot is at a workplace, not a classroom, so nothing
changes for them, and nobody decides on a child's behalf. Whether a child
avoids depends only on the level at their venue and the draw of their decision
unit.

With `shared_venue: residence` every member of a household has the same draw,
which acts as that household's tolerance for risk. Each child compares the level
at their own venue with it:

- Siblings in one classroom, or at venues with equal levels, always match. One
  is never home while the other is at school.
- Siblings in different classrooms or at different schools each use their own
  venue's level against the same draw. The child at the busier venue goes home
  first, and the sibling follows once their own level also passes the draw.
- With a `parent_venue` metric the level is the school's, so siblings at one
  school match and siblings at different schools can differ.
- With `decision_unit: person` every child has their own draw, so siblings can
  differ even at equal levels.

Without `applies_to`, the policy also reaches adults who are registered to and
scheduled at a classroom, such as teachers. They follow the same rule with their
own household's draw. Add `applies_to` as above to limit the policy to children.

`applies_to` selects on properties of the person, such as age or a property the
world carries (`properties.<name>`). It cannot select on a person's role within a
venue, such as the `student` and `worker` subsets of a classroom, so age stands
in for "pupil" here. If your world has a person property for occupation, use that
instead.

With a cumulative metric such as deaths the level only rises, so avoidance is
permanent. With a current metric it falls when the venue calms down and people
return.

On day 0 the level is `base`, because nothing has been counted yet. With a
`base` above 0 avoidance therefore starts in the first slot, where a threshold
trigger would still be open.

A shared decision unit needs every person it covers to have exactly one venue
for the shared activity. That is checked once at startup for the people the
policy applies to and the metric can count. If anyone fails, the run stops with
an error naming the lowest failing person id, the activity and the number of
venues found. Fix the world data or use `decision_unit: person`.

Symptom policies can use a response too, for example symptomatic people staying
home more as the total number of deaths rises:

```yaml
observed_metrics:
  deaths: {source: death, scope: world}

policies:
  symptom_policies:
    - name: "stay_home_as_deaths_rise"
      symptoms: ["mild", "severe"]
      override_activities: ["primary_activity", "leisure"]
      replacement: "residence"
      trigger:
        response:
          terms:
            - {metric: deaths, weight: 0.002}
        decision_unit: person
```

A symptom policy with a response needs a world scope metric and
`decision_unit: person`, and cannot set `replacement_schedule` or
`follow_up_policy`. A symptom policy cannot name a policy with a response as its
`follow_up_policy` either.

The loader refuses any policy with a response that also sets `compliance_rate`
(even 1) or `transmission_effects_file`, or that mixes the threshold fields
(`metric`, `activate_at`, `deactivate_at`, `calm_days`, `tolerance`,
`tolerance_pct`) into the same trigger. `decision_unit` belongs to the response
form only.

---

## 11. Vaccines — `vaccines.yaml`

`vaccines.yaml` has three top level keys: `enabled`, `vaccines` and `campaigns`.
A vaccine says how well each dose protects and for how long. A campaign says who
is offered which vaccine, how many doses there are each day and in what order
people are served. If the file is missing, vaccination is off and a warning is
printed. `enabled: false` turns it off too.

### Vaccines and doses

```yaml
enabled: true

vaccines:
  Shot:
    doses:
      - days_to_onset:     14
        days_to_effective: 28
        days_to_waning:    180
        days_to_finished:  360
        waning_factor:     0.5
        infection_efficacy:
          covid19: {"0-17": 0.6, "18-100": 0.7}
        symptom_efficacy:
          covid19: {"0-100": 0.7}
      - days_to_effective: 14
        days_to_waning:    180
        days_to_finished:  360
        infection_efficacy:
          covid19: 0.9
        symptom_efficacy:
          covid19: 0.9
      - days_to_effective: 7
        days_to_waning:    180
        days_to_finished:  360
        infection_efficacy:
          covid19: 0.95
        symptom_efficacy:
          covid19: 0.95
  OtherShot:
    doses:
      - {days_to_effective: 14, days_to_waning: 180, days_to_finished: 360,
         infection_efficacy: {covid19: 0.5}, symptom_efficacy: {covid19: 0.5}}
      - {days_to_effective: 14, days_to_waning: 180, days_to_finished: 360,
         infection_efficacy: {covid19: 0.8}, symptom_efficacy: {covid19: 0.8}}
```

- A vaccine is a list of doses. The first in the list is dose 0, the next is
  dose 1, and so on. Tiers ask for a dose by this number.
- `days_to_effective`, `days_to_waning` and `days_to_finished` are required. If
  one is missing the load stops and the message names the vaccine and the dose.
- `days_to_onset` is optional and defaults to 0. It must be between 0 and
  `days_to_effective`.
- `waning_factor` is optional and defaults to 1.
- Every time is in days since that dose was given.

Protection follows the dose like this:

| Days since the dose | Protection |
|---|---|
| Before `days_to_onset` | None |
| Onset to `days_to_effective` | A straight rise to the peak |
| `days_to_effective` to `days_to_waning` | The peak |
| `days_to_waning` to `days_to_finished` | A straight fall to peak times `waning_factor` |
| After `days_to_finished` | Peak times `waning_factor` |

If onset equals `days_to_effective` there is no rise, only a step. With onset
14 and effective 14 a 70% dose gives 0 on day 13 and 70% on day 14. With onset 0
and effective 14 it gives 35% on day 7. With `waning_factor` left at 1 protection
never falls.

`infection_efficacy` and `symptom_efficacy` are set per disease and per age band.
The key under them is the disease's name. A band is `"min-max"` in whole years,
both ends included. A plain number covers ages 0 to 100. A person outside every
band, or a disease with no entry, gets no protection from that dose. A person with
several doses gets the highest value of any dose, not the sum. This applies to
infection and symptom efficacy separately.

Any other key in a dose is an error, so a misspelt `days_to_onset` cannot leave
the onset at 0. So is a band that is not written `min-max`, or an efficacy
outside 0 to 1.

### What symptom efficacy does

Infection efficacy lowers susceptibility. The infection pressure on a person is
multiplied by 1 minus the infection efficacy on that day.

Symptom efficacy acts at the moment of infection, on the outcome. Each trajectory
in `disease.yaml` has a probability for that person and a `severity`. Every
trajectory more severe than the least severe one loses the fraction
`symptom_efficacy` of its probability. All of it moves to the least severe
trajectory, so the total stays the same. A 70% symptom efficacy takes 70% of the
probability of every more severe trajectory, whatever its severity. If two
trajectories share the lowest severity, the first one listed receives the
probability and the other is left alone. The efficacy used is the one on the day
of infection.

The two stack. With 70% infection efficacy and 70% symptom efficacy, the
infection pressure is 0.3 of its unvaccinated value and an infection carries 0.3
of the chance of a more severe outcome. At low pressure that is 0.3 x 0.3 = 0.09
of the unvaccinated chance, so the reduction in more severe outcomes per exposure
is 91%, not 70%.

### A campaign

```yaml
campaigns:
  "Autumn rollout":
    vaccine_type:      "Shot"
    start_date:        "2024-09-02"
    end_date:          "2024-12-01"
    days_to_next_dose: [28]
    dose_budget:       {doses: 70000, per_days: 7}
    coverage_target:   0.7
    tiers:
      - name: everyone_dose_1
        dose: 1
      - name: everyone_dose_0
        dose: 0
```

| Key | Required | What it does |
|---|---|---|
| `vaccine_type` | yes | The name of a vaccine in this file |
| `tiers` | yes | The ordered list of groups to serve |
| `start_date` | no | First day of the campaign, `YYYY-MM-DD`. Absent: the first day of the run |
| `end_date` | no | Last day of the campaign, included. Absent: the end of the run |
| `days_to_next_dose` | no | Days to wait after each dose |
| `courses_begun_with` | no | The vaccines a course may have begun with, for a booster campaign |
| `dose_budget` | no | Doses per period. Absent: no limit |
| `coverage_target` | no | Stop starting new courses at this share of the population |

Any other key is an error, and the message lists the valid ones. A campaign with
no `tiers` is an error.

### Dose budget

`dose_budget: {doses: D, per_days: P}` gives D doses every P days. Both are whole
numbers above 0. On campaign day k the campaign may give `floor((k+1) x D / P)`
minus `floor(k x D / P)` doses. Over any whole number of periods that adds up to
exactly D per period, with no drift. For `{doses: 10, per_days: 3}`
the days give 3, 3, 4, then 3, 3, 4 again.

Unspent doses are lost. If 40 people are due and the budget for the day is 100,
40 are given and the other 60 are not carried to tomorrow. The budget is shared by
every tier and every dose of the course. With no `dose_budget` there is no limit:
everyone who is due and accepts is dosed that day.

Campaign day k is the run day minus the run day of `start_date`. With no
`start_date` it is the run day, so campaign day 0 is run day 0.

### Coverage target

`coverage_target` caps first doses at `floor(target x population at day 0)`. A
target of 0.7 in a world of 600000 people allows 420000 first doses. Second and
later doses of courses already started carry on after the cap is reached. Tiers
for a later dose are not cut by the cap, so they still use the budget.

Coverage counts people who have had a first dose under this campaign's vaccine,
including people given it by another campaign of the same vaccine and people who
have since died. If uptake, eligibility or the budget makes the target
unreachable, the run goes on and the end of run line says so (see Output).

### Tiers

```yaml
    tiers:
      - name: over_65_dose_1
        dose: 1
        selection:
          - {property: age, operator: ">=", value: 65}
        uptake: 0.95
      - name: over_65_dose_0
        dose: 0
        selection:
          - {property: age, operator: ">=", value: 65}
        uptake: 0.8
      - name: everyone_else_dose_1
        dose: 1
      - name: everyone_else_dose_0
        dose: 0
```

Each tier has:

- `name`: required, and different from every other tier in the campaign.
- `dose`: required, a whole number from 0. The vaccine must have that dose. 0 is
  the first.
- `selection`: the usual selection syntax. Empty or absent means everyone.
- `uptake`: see Uptake. Absent means everyone accepts.
- `decision_unit`: see Decision unit. Absent means `person`.

Tiers are served in the order listed, strictly. A tier gets doses only after every
tier above it has taken all its candidates for the day. The coverage cap is the one
exception (see Coverage target). A person's tier on a given day is the first tier
whose group they are in and whose `dose` is the number of doses they already hold.
A person with no doses is served by a `dose: 0` tier, and a person with one dose by
a `dose: 1` tier.

To put second doses before new first doses, list the `dose: 1` tier above the
`dose: 0` tier, as above.

Inside a tier that the budget serves only in part, who is served comes from a
seeded lottery. The lottery order depends on the run seed and on the positions of
the campaign and the tier in the file. It does not change from day to day, so
someone passed over today is next in line tomorrow.

A campaign can have at most 64 tiers.

### Later doses

`days_to_next_dose[n-1]` is the number of days to wait after dose n before the
next one, counted from the day dose n was given. In `[28]` a second dose is due 28
days after the first. A missing entry means no wait. A person is never dosed twice
on the same day.

A later dose continues the vaccine the course began with. A campaign whose
`vaccine_type` differs from the vaccine of a person's first dose does not serve
that person. To serve people who began with another vaccine, a booster campaign
lists them in `courses_begun_with`:

```yaml
campaigns:
  "Booster":
    vaccine_type: "Shot"
    courses_begun_with: ["Shot", "OtherShot"]
    tiers:
      - {name: booster, dose: 2}
```

The list is checked against the vaccine of the person's first dose, and it is
not checked for people who hold no doses. Dose 2 is the third entry in the `doses`
list of `Shot`, so `Shot` needs at least three doses for the loader to accept this.
A booster tier names the dose people are due next, which after two doses is dose 2.

### Uptake

A tier's `uptake` is the probability that a person or a unit accepts the dose. It
takes one of three forms.

A number, for everyone in the tier:

```yaml
        uptake: 0.9
```

A list of rows. The first row that matches wins, and a row with no `selection` is
the catch all:

```yaml
        uptake:
          - selection:
              - {property: properties.ethnicity, operator: "==", value: "W"}
              - {property: age, operator: ">=", value: 65}
            probability: 0.95
          - selection:
              - {property: age, operator: "<", value: 18}
            probability: 0.6
          - probability: 0.8
```

A CSV file:

```yaml
        uptake: {file: uptake_first_dose.csv}
```

```csv
filter.age,filter.properties.ethnicity,probability
65-120,,0.95
18-64,W,0.8
18-64,,0.7
,,0.5
```

The file has `filter.*` columns and one `probability` column, and no others. A
filter cell can be a value, a range such as `18-64`, or a comparison such as
`<18` or `>=65`. An empty cell matches anyone. Rows are tried in order, and the
first match wins. The path is used as written if that file exists, otherwise it is
read relative to the folder of the YAML file. The file is read at startup.

Any attribute the world records can be a filter, for example
`properties.ethnicity`.

Uptake is settled once, at startup, for each decision unit in each tier. It
depends only on the person's attributes, their id and the run seed. A decline is
final for that tier. A person who declined is not asked again, and the decline
does not pass them on to a later tier. A tier for a later dose has its own uptake.
Use it for the share of people who finish a course they started.

If no row matches a person whose probability is needed, startup stops and the
error names the tier and the person. End the table with a row that has no filters
if everyone should be covered.

### Decision unit

`decision_unit` is `person` (the default) or a shared unit:

```yaml
        decision_unit:
          shared_venue: residence
          scope:
            - {property: "activities.residence.venue_type", operator: "==", value: "household"}
          deciders:
            - {property: age, operator: ">=", value: 18}
          combine: mean
          if_no_decider: error
```

Everyone in `scope` who is registered to the same venue for the named activity
decides together. The unit makes one draw, so all its members in the tier's group
get the same answer.

- `shared_venue`: the activity name. Required.
- `scope`: who is inside a unit. People outside it decide alone. Required. Write
  it with `==` on the venue type you want in, as above. This keeps care homes and
  halls out.
- `deciders`: the members whose uptake rows set the unit's probability. Required.
  A decider does not have to be in the tier's group. A child can be vaccinated on
  a parent's decision, and the parent is not offered a dose by that tier.
- `combine`: how the deciders' probabilities become one: `min`, `mean` or `max`.
  The default is `mean`.
- `if_no_decider`: what a unit with no decider does. `error` stops the run at
  startup and is the default. `individual` has each member decide alone.

In a shared tier the uptake table is read for the deciders. People outside `scope`
use their own row.

Each person in a unit must be registered to exactly one venue for the activity. A
unit must also be held whole by one rank, which means its venue and its members are
in the same geo unit. If either is not true, startup fails and the message names the
person and the venue.

### Several campaigns

Campaigns run in file order. Each has its own budget and its own cap. A person gets
at most one dose a day, from the first campaign that reaches them. Nothing is
shared between campaigns except that coverage counts a vaccine's first doses
whichever campaign gave them.

### Restart

A checkpoint stores each person's doses and nothing about a campaign. Each day the
manager counts again who has started a course. A resumed run gives the same doses
as the run it continues.

### Output

- `/events/vaccinations`: one row per dose.
- `/events/campaign_days`: one row per campaign, tier and active day, with these
  columns: `day` (the run day), `campaign_id`, `tier_id`, `candidates` (people in
  this tier today who accepted), `doses_given`, `declined` (people in this tier
  today who declined), `allotment` (the day's budget, or `-1` when there is no
  budget) and `started_total` (people who have started a course by the end of the
  day). A person is in a tier on a day when it is their tier that day and their
  dose is due. A run that logs no rows has no `/events/campaign_days` in its file at
  all.
- `/metadata/registries/campaigns`: the campaign names, in file order. This is
  what `campaign_id` indexes.
- `/metadata/registries/campaign_tiers/<campaign_id>`: the tier names of that
  campaign, in order, for the campaign's number. This is what `tier_id` indexes.

At startup, rank 0 prints one line for each tier:

```
vaccination uptake: campaign '<campaign>' tier '<tier>': selected <n>, accepted <n>, declined <n>, shared units <n>, checked <n>, failed 0
```

At the end of the run it prints one line for each campaign:

```
vaccination campaign '<campaign>': <n> of <cap> first doses the coverage target allows (<p>%), target reached
```

The last words are `target not reached` when the cap was not reached. A campaign
with no coverage target prints `<n> people started a course (no coverage target)`.

### A worked example

A two dose vaccine with a 14 day onset step and 70% efficacy. The budget is 10000
doses every 30 days, and first doses stop at 70% of the population. Children are
served first, then adults. Within each group second doses come before first doses.
Children's first doses are decided by the adults in their home.

```yaml
enabled: true

vaccines:
  Shot:
    doses:
      - days_to_onset:     14
        days_to_effective: 14
        days_to_waning:    30
        days_to_finished:  60
        infection_efficacy:
          covid19: 0.7
        symptom_efficacy:
          covid19: 0.7
      - days_to_onset:     14
        days_to_effective: 14
        days_to_waning:    30
        days_to_finished:  60
        infection_efficacy:
          covid19: 0.7
        symptom_efficacy:
          covid19: 0.7

campaigns:
  "National rollout":
    vaccine_type:      "Shot"
    start_date:        "2024-07-01"
    days_to_next_dose: [28]
    dose_budget:       {doses: 10000, per_days: 30}
    coverage_target:   0.7
    tiers:
      - name: children_dose_1
        dose: 1
        selection:
          - {property: age, operator: "<", value: 18}
        uptake: 0.95
      - name: children_dose_0
        dose: 0
        selection:
          - {property: age, operator: "<", value: 18}
        uptake:
          - selection:
              - {property: age, operator: "<", value: 40}
            probability: 0.6
          - probability: 0.85
        decision_unit:
          shared_venue: residence
          scope:
            - {property: "activities.residence.venue_type", operator: "==", value: "household"}
          deciders:
            - {property: age, operator: ">=", value: 18}
          combine: mean
          if_no_decider: individual
      - name: adults_dose_1
        dose: 1
        selection:
          - {property: age, operator: ">=", value: 18}
        uptake: 0.95
      - name: adults_dose_0
        dose: 0
        selection:
          - {property: age, operator: ">=", value: 18}
        uptake:
          - selection:
              - {property: age, operator: ">=", value: 65}
            probability: 0.9
          - probability: 0.7
```

Each dose gives nothing for 14 days, then 70%. `waning_factor` is left out, so
protection never falls.

### Errors you will meet

A load error starts with `Error loading vaccination file '<path>': `, which the
table leaves out. A startup error stops the run before day 0. In the table
`<...>` stands for your own names and numbers.

| Mistake | Message | Fix |
|---|---|---|
| A typo in a campaign key, such as `dose_budjet` | `vaccination campaign '<campaign>': unknown field 'dose_budjet'. Fields are start_date, end_date, vaccine_type, days_to_next_dose, courses_begun_with, dose_budget, coverage_target, tiers.` | Use one of the listed names |
| No `vaccine_type` | `vaccination campaign '<campaign>': needs 'vaccine_type'` | Name a vaccine from `vaccines:` |
| A typo in a dose key, such as `days_to_onsett` | `vaccine '<vaccine>' dose <n>: unknown field 'days_to_onsett'. Fields are days_to_onset, days_to_effective, days_to_waning, days_to_finished, waning_factor, infection_efficacy, symptom_efficacy.` | Use one of the listed names |
| An age band not written `min-max`, such as `0_17` | `vaccine '<vaccine>' dose <n> infection_efficacy <disease> band '0_17' must be written min-max in whole years, such as 0-17` | Write the band as `0-17` |
| An efficacy outside 0 to 1 | `vaccine '<vaccine>' dose <n> infection_efficacy <disease> efficacy must be a number from 0 to 1, got '1.5'` | Use a number from 0 to 1 |
| A dose with no timing field | `vaccine '<vaccine>' dose <n>: missing required field 'days_to_waning'` | Add the field to that dose |
| Onset after the peak | `vaccine '<vaccine>' dose <n>: days_to_onset must be between 0 and days_to_effective, got 20 with days_to_effective 14` | Lower `days_to_onset` or raise `days_to_effective` |
| A tier asks for a dose the vaccine does not have | `vaccination campaign '<campaign>' tier '<tier>' dose 2 but vaccine 'Shot' has 2 doses (0 to 1)` | Use a dose in that range, or add the dose to the vaccine |
| A probability outside 0 to 1 | `vaccination campaign '<campaign>' tier '<tier>' uptake must be a number from 0 to 1, got '1.5'` | Use a number from 0 to 1 |
| A person that no uptake row matches | `vaccination uptake: campaign '<campaign>' tier '<tier>': no uptake row matches person <id>, whose probability is needed. <n> of <m> people checked failed. Add a row with empty filters to the uptake table as a catch-all.` | End the table with a row that has no `selection` |
| A shared unit with no decider | `vaccination uptake: campaign '<campaign>' tier '<tier>': the unit at venue <id> for activity '<activity>' has people in the group but no member matching the deciders, and if_no_decider is 'error'. <n> of <m> units checked failed. Widen the deciders or set if_no_decider to 'individual'.` | Widen `deciders`, or set `if_no_decider: individual` |
| A shared unit not held whole by one rank | `vaccination uptake: campaign '<campaign>' tier '<tier>': person <id> is registered for activity '<activity>' to venue <id>, which is missing or in another geo unit than the person, so the shared unit is not held whole by one rank. <n> of <m> people checked failed.` | Fix the world data, or use `decision_unit: person` |

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
```

The loader reads only `enabled`, `partitioning` and `chunked_loading`.

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
```

- `deterministic_activities` are pre-computed once at world load and
  cached. `leisure` etc. stay stochastic at runtime.

---

## 14. The `data/` directory

Reference data shared across scenarios. Edit (or replace) these CSVs
to retune the model without touching code.

| File | Purpose | Schema |
|---|---|---|
| `infection_outcome_rates_<disease>.csv` | Per (age, sex, pop) outcome probabilities | `filter.age,filter.sex,asymptomatic,mild,severe,hospital,icu,home_ifr,hospital_ifr,icu_ifr` |
| `bulk_seeds_<disease>.csv` | Bulk seeding events | see §6 |
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
(set-identical events, every field included).

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
- A checkpoint stores the running totals of every metric a trigger reads
  and the state of every trigger. A cumulative total cannot be rebuilt, so
  a resume must use the same triggered policies (same kind, name and form,
  counted over the same scope, in any order) and read the same set of
  metrics. Anything else stops the run with an error that lists both sets.
  Declared metrics that no trigger reads are not stored and may change
  freely.
- Untriggered policies may change between checkpoint and resume: fewer,
  more, reordered, or with edited selection criteria. Policies are
  matched by name within their kind (symptom or temporal), which is why
  names must be unique per kind (a duplicate is refused at load, on any
  run). Who each policy applies to is recomputed from the current config.
  A person's compliance decision carries over when the policy still
  exists and still applies to them; otherwise it is dropped. New policies
  start undecided. A person held in place by a removed policy is released
  at resume. Rank 0 prints what was dropped, added or moved; nothing is
  printed when the list is unchanged.
- Draws are keyed by policy name, so a policy that only moved position
  draws as it would have in an uninterrupted run. A new or renamed policy
  draws afresh.
- A checkpoint whose `format_version` (in its `manifest.yaml`) differs from
  the one this build writes is always refused, and the error names both
  versions. Start the run again from day 0.

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
| Load error `observed_metrics '<name>': unknown activity` or `unknown stage` | The activity is not in the world, or the stage is not a symptom tag in `disease.yaml`. The stage error lists the known tags. |
| Load error `... trigger: 'deactivate_at' is required` (or `'calm_days' is required`) | A current metric needs `deactivate_at`, and a cumulative metric needs `calm_days` and rejects `deactivate_at`. |
| Load error `'compliance_rate' must be absent on a policy with a 'response'` | Remove the key, even a value of 1. The avoidance level replaces it. |
| Startup error `... exactly 1 is required` naming a person | That person has no venue, or several, for the `shared_venue` activity. Fix the world data or use `decision_unit: person`. |
| A triggered policy never seems to act | Check the metric name and scope, and that `override_activities` and `override_venue_types` match where the metric counts. A trigger is first evaluated at the end of day 0. |
| Resume fails with `the checkpoint holds metrics [...] but this run has [...]` | The run reads a different set of metrics than the one that wrote the checkpoint. Resume with the same triggered policies. |

---

*Anything you can't reach through these files is a missing config
hook, not a reason to patch C++. File an issue or extend the
loader; don't bake a value into the source.*
