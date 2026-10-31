// Config files are checked strictly: a key no loader reads is an error, and a
// file that was named but cannot be used is an error. These tests put one
// stray key at every structural level of a minimal valid file and expect the
// loader to name it, and they load every shipped config directory so a dead key
// cannot creep back in.

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "activity/on_the_fly_venue_allocator.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/infection_seed.h"
#include "epidemiology/policy.h"
#include "loaders/config_loader.h"
#include "loaders/disease_loader.h"
#include "loaders/policy_loader.h"
#include "test_utils.h"

using namespace june;
namespace fs = std::filesystem;

namespace {

constexpr const char* kStray = "zz_stray";

using Path = std::vector<std::string>;

// `yaml` with a stray key added to the mapping at `path`. A path step that is
// all digits indexes a list; anything else is a mapping key. An empty path is
// the root.
std::string withStrayKey(const std::string& yaml, const Path& path) {
  YAML::Node root = YAML::Load(yaml);
  YAML::Node node = root;
  for (const auto& step : path) {
    const bool is_index =
        !step.empty() && std::all_of(step.begin(), step.end(), [](char c) {
          return std::isdigit(static_cast<unsigned char>(c));
        });
    REQUIRE_MESSAGE(node.IsDefined(),
                    "no such level in the test file: " << step);
    // reset() rebinds the handle; assigning would overwrite the parent.
    if (is_index && node.IsSequence()) {
      node.reset(node[static_cast<size_t>(std::stoul(step))]);
    } else {
      REQUIRE_MESSAGE(node[step].IsDefined(),
                      "no such level in the test file: " << step);
      node.reset(node[step]);
    }
  }
  REQUIRE_MESSAGE(node.IsMap(), "level is not a mapping");
  node[kStray] = 1;
  return YAML::Dump(root);
}

std::string describe(const Path& path) {
  std::string joined = "root";
  for (const auto& step : path) joined += "/" + step;
  return joined;
}

struct LoaderCase {
  std::string name;
  std::string yaml;
  std::vector<Path> levels;
  // Loads the file at the given path (may write sibling files first).
  std::function<void(const ScopedTestFiles&, const std::string&)> load;
};

void expectEveryLevelRejected(const LoaderCase& c) {
  {
    ScopedTestFiles files{"june_strict_valid"};
    const auto file = files.write("config.yaml", c.yaml);
    INFO("the unmodified " << c.name << " file must load");
    CHECK_NOTHROW(c.load(files, file.string()));
  }
  for (const auto& level : c.levels) {
    ScopedTestFiles files{"june_strict_stray"};
    const auto file = files.write("config.yaml", withStrayKey(c.yaml, level));
    INFO(c.name << ": stray key at " << describe(level));
    CHECK_THROWS_WITH_AS(c.load(files, file.string()),
                         doctest::Contains("unknown field"),
                         std::runtime_error);
    CHECK_THROWS_WITH_AS(c.load(files, file.string()),
                         doctest::Contains(kStray), std::runtime_error);
  }
}

// ---------------------------------------------------------------------------
// Minimal valid files
// ---------------------------------------------------------------------------

const char* kSimulationYaml = R"(
time:
  start_date: "2020-03-02"
  end_date: "2020-03-09"
config_paths:
  disease_file: disease.yaml
random_seed: 1
verbose: false
output:
  stats_interval_days: 1
partial_presence:
  enabled_venue_types:
    train_line: {target_group_size: 100}
checkpoint:
  enabled: false
)";

const char* kPerformanceYaml = R"(
performance:
  precompute_schedules: true
  deterministic_activities: [residence]
  hybrid_activities: [primary_activity]
  stochastic_activities: [leisure]
)";

const char* kParallelYaml = R"(
parallel:
  enabled: true
  partitioning:
    level: MGU
    centroids_file: centroids.csv
    adjacency_file: adjacency.json
    metis:
      imbalance_tolerance: 0.05
  chunked_loading:
    person_metadata_chunk_size: 1000
    geo_unit_chunk_size: 10
)";

const char* kMultiModeMatricesYaml = R"(
contact_matrices:
  tube:
    modes:
      airborne:
        bins: [all]
        beta: 1.0
        contacts: [[1.0]]
allow_default_matrix: true
default_contacts_matrix:
  modes:
    airborne:
      bins: [all]
      beta: 1.0
      contacts: [[1.0]]
)";

const char* kSingleModeMatricesYaml = R"(
contact_matrices:
  pub:
    bins: [all]
    beta: 1.0
    contacts: [[1.0]]
default_contacts_matrix:
  bins: [all]
  contacts: [[1.0]]
)";

const char* kScheduleYaml = R"(
day_type_cycle: [weekday, weekend]
default_schedule_type: adult
schedule_types:
  adult:
    priority: 1
    is_temporary: false
    selection:
      - {property: age, operator: ">", value: 17}
    weekday:
      - name: day
        start: "08:00"
        end: "18:00"
        allowed_activities: [residence, leisure]
        specified_activity: {type: primary_activity, index: 0, venue_type: office}
        hop_on_activity:
          commute: other_schedule
          leisure:
            property: work_mode
            schedule_template: "{value}_schedule"
    weekend:
      - name: day
        start: "08:00"
        end: "18:00"
        allowed_activities: [residence]
    force_hybrid_activities: [leisure]
    linked_activities: [leisure]
    participation:
      leisure:
        weekday: 0.5
        weekend: 0.7
)";

const char* kPreferencesYaml = R"(
profiles:
  - name: adults
    activity: leisure
    priority: 1
    selection:
      - {property: age, operator: ">", value: 17}
    weights:
      pub: 2.0
)";

const char* kSeedsYaml = R"(
global_parameters:
  base_cases_per_capita: 0.001
  default_seed_strength: 1.0
infection_seeds:
  - name: first
    type: uniform
    date: "2020-03-02 08:00"
    trajectory_key: general
    parameters:
      seed_strength: 1.0
      cases_per_capita_multiplier: 1.0
      attribute_filters:
        sex: male
  - name: second
    type: exact
    date: "2020-03-03 08:00"
    geo_level: MGU
    parameters:
      age_groups: ["0-17"]
      units:
        E02001234: 3
)";

const char* kEncountersYaml = R"(
coordinated_encounters:
  enabled: true
  encounters:
    - name: social
      enabled: true
      network: friendships
      trigger_slots: [evening]
      allowed_venues: [pub]
      proposal_probability: 0.1
      acceptance_probability: 0.5
      invite_distribution: {type: poisson, mean: 2}
      daily_max_distribution: {type: fixed, count: 3}
      min_attendees: 2
      priority: 1
  follows:
    - name: toddlers
      pool_venue_type: household
      network: parents
      encounter_type: follow
      probability: 1.0
      log: false
      establishment: criteria
      span: hop
      eligibility:
        - {property: age, operator: "<", value: 5}
      host_eligibility:
        - {property: age, operator: ">", value: 17}
      activity_exceptions: [leisure]
      venue_exceptions: [school]
      follower_activity_exceptions: [residence]
)";

const char* kOnTheFlyYaml = R"(
rules:
  visit:
    strategy: hosting_geo_unit
    venue_type: pub
    venue_stability: daily
    geo_unit_level: MGU
activity_rules:
  leisure: visit
)";

const char* kPoliciesYaml = R"(
policies:
  temporal_policies:
    - name: close
      start_date: "2020-03-09"
      end_date: "2020-04-09"
      override_activities: [leisure]
      override_venue_types: [pub]
      exempt_venue_types: [school]
      replacement: residence
      compliance_rate: 0.5
      applies_to:
        - {property: age, operator: ">", value: 17}
      exempt:
        - activity: leisure
          selection:
            - {property: age, operator: "<", value: 5}
  symptom_policies:
    - name: isolate
      symptoms: [mild]
      start_time: 0
      override_activities: [leisure]
      compliance_rate: 0.9
      inherit_compliance: true
      inherit_refusal: false
)";

const char* kVaccinesYaml = R"(
enabled: true
vaccines:
  Jab:
    doses:
      - days_to_onset: 1
        days_to_effective: 5
        days_to_waning: 6
        days_to_finished: 7
        waning_factor: 1.0
        infection_efficacy:
          covid19:
            "0-100": 0.5
campaigns:
  first:
    vaccine_type: Jab
    start_date: "2020-03-09"
    tiers:
      - name: everyone
        dose: 0
)";

const char* kTrajectoryDiseaseYaml = R"(
disease:
  name: strict_test
  settings:
    fatality_stage:
      - name: dead
    recovered_stage:
      - name: recovered
  symptom_tags:
    - {name: recovered, value: -2}
    - {name: healthy, value: -1}
    - {name: exposed, value: 0}
    - {name: mild, value: 1}
    - {name: dead, value: 2}
  trajectories:
    - description: mild course
      selection_key: mild
      probability: 1.0
      severity: 0.0
      infectiousness_factor: 1.0
      stages:
        - symptom_tag: exposed
          completion_time: {type: constant, value: 2.0}
        - symptom_tag: mild
          completion_time: {type: lognormal, s: 0.5, loc: 0.0, scale: 2.0}
        - symptom_tag: recovered
          completion_time: {type: constant, value: 0.0}
  transmission:
    mode: Trajectory-Driven
    type: gamma
    max_infectiousness: {type: lognormal, s: 0.5, loc: 0.0, scale: 1.0}
    shape: {type: constant, value: 1.56}
    rate: {type: constant, value: 0.53}
    shift: {type: constant, value: -2.12}
    modes:
      - {name: airborne, mode_transmissibility_multiplier: 1.0}
  immunity:
    level: 0.9
    waning_rate: 0.01
)";

const char* kStageDrivenDiseaseYaml = R"(
disease:
  name: strict_stage_test
  settings:
    recovered_stage:
      - name: recovered
  symptom_tags:
    - {name: recovered, value: -2}
    - {name: healthy, value: -1}
    - {name: exposed, value: 0}
    - {name: mild, value: 1}
  outcome_rates_csv:
    file: rates.csv
    populations:
      gp: {name: general_population}
    sexes: [male, female]
    outcomes:
      mild: {symptom_tag: mild}
  trajectories:
    - description: mild course
      selection_key: mild
      stages:
        - symptom_tag: exposed
          completion_time: {type: constant, value: 2.0}
        - symptom_tag: mild
          completion_time: {type: constant, value: 3.0}
  transmission:
    mode: Stage-Driven
    modes:
      - name: airborne
        mode_transmissibility_multiplier: 1.0
      - name: surface
        type: fomite
        max_age: 10
        sub_bin_time: 0
        fomite_curve: {type: constant, value: 1.0}
        deposition_stages:
          mild: {type: constant, value: 1.0}
      - name: air_pool
        type: compartmental_deposition
        deposition_stages:
          mild: {type: linear_ramp, start_value: 0.0, end_value: 1.0, ramp_duration: 2.0}
      - name: air_uptake
        type: compartmental_uptake
    stage_curves:
      airborne:
        mild: {type: gamma, max_infectiousness: 1.0, shape: 1.5, rate: 0.5, shift: 0.0}
        exposed: {type: exponential_decay, initial_value: 1.0, decay_rate: 0.5, delay: 0.0}
)";

void writeRatesCsv(const ScopedTestFiles& files) {
  files.write("rates.csv", "mild\n");
}

void loadPoliciesFile(const std::string& file) {
  WorldState world;
  PolicyManager manager(world);
  PolicyLoader::loadPolicies(manager, file, "2020-03-02");
}

std::vector<LoaderCase> loaderCases() {
  std::vector<LoaderCase> cases;

  cases.push_back({"simulation.yaml",
                   kSimulationYaml,
                   {{},
                    {"time"},
                    {"config_paths"},
                    {"output"},
                    {"partial_presence"},
                    {"partial_presence", "enabled_venue_types", "train_line"},
                    {"checkpoint"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     ConfigLoader::loadSimulation(f);
                   }});

  cases.push_back({"performance.yaml",
                   kPerformanceYaml,
                   {{}, {"performance"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     ConfigLoader::loadPerformance(f);
                   }});

  cases.push_back({"parallel.yaml",
                   kParallelYaml,
                   {{},
                    {"parallel"},
                    {"parallel", "partitioning"},
                    {"parallel", "partitioning", "metis"},
                    {"parallel", "chunked_loading"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     ConfigLoader::loadParallel(f);
                   }});

  cases.push_back({"contact_matrices.yaml (modes)",
                   kMultiModeMatricesYaml,
                   {{},
                    {"contact_matrices", "tube"},
                    {"contact_matrices", "tube", "modes", "airborne"},
                    {"default_contacts_matrix"},
                    {"default_contacts_matrix", "modes", "airborne"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     ConfigLoader::loadContactMatrices(f);
                   }});

  cases.push_back(
      {"contact_matrices.yaml (single matrix)",
       kSingleModeMatricesYaml,
       {{}, {"contact_matrices", "pub"}, {"default_contacts_matrix"}},
       [](const ScopedTestFiles&, const std::string& f) {
         ConfigLoader::loadContactMatrices(f);
       }});

  cases.push_back(
      {"schedules.yaml",
       kScheduleYaml,
       {{},
        {"schedule_types", "adult"},
        {"schedule_types", "adult", "selection", "0"},
        {"schedule_types", "adult", "weekday", "0"},
        {"schedule_types", "adult", "weekday", "0", "specified_activity"},
        {"schedule_types", "adult", "weekday", "0", "hop_on_activity",
         "leisure"},
        {"schedule_types", "adult", "participation", "leisure"}},
       [](const ScopedTestFiles&, const std::string& f) {
         ConfigLoader::loadSchedule(f);
       }});

  cases.push_back({"activity_preferences.yaml",
                   kPreferencesYaml,
                   {{}, {"profiles", "0"}, {"profiles", "0", "selection", "0"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     ConfigLoader::loadActivityPreferences(f);
                   }});

  cases.push_back({"infection_seeds.yaml",
                   kSeedsYaml,
                   {{},
                    {"global_parameters"},
                    {"infection_seeds", "0"},
                    {"infection_seeds", "0", "parameters"},
                    {"infection_seeds", "1"},
                    {"infection_seeds", "1", "parameters"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     InfectionSeedConfigLoader::loadFromFile(f);
                   }});

  cases.push_back(
      {"coordinated_encounters.yaml",
       kEncountersYaml,
       {{},
        {"coordinated_encounters"},
        {"coordinated_encounters", "encounters", "0"},
        {"coordinated_encounters", "encounters", "0", "invite_distribution"},
        {"coordinated_encounters", "encounters", "0", "daily_max_distribution"},
        {"coordinated_encounters", "follows", "0"},
        {"coordinated_encounters", "follows", "0", "eligibility", "0"}},
       [](const ScopedTestFiles&, const std::string& f) {
         ConfigLoader::loadCoordinatedEncounters(f);
       }});

  cases.push_back({"on_the_fly_venues.yaml",
                   kOnTheFlyYaml,
                   {{}, {"rules", "visit"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     OnTheFlyVenueAllocator allocator(f);
                   }});

  cases.push_back({"policies.yaml",
                   kPoliciesYaml,
                   {{},
                    {"policies"},
                    {"policies", "temporal_policies", "0"},
                    {"policies", "temporal_policies", "0", "applies_to", "0"},
                    {"policies", "temporal_policies", "0", "exempt", "0"},
                    {"policies", "symptom_policies", "0"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     loadPoliciesFile(f);
                   }});

  cases.push_back({"vaccines.yaml",
                   kVaccinesYaml,
                   {{},
                    {"vaccines", "Jab"},
                    {"vaccines", "Jab", "doses", "0"},
                    {"campaigns", "first"},
                    {"campaigns", "first", "tiers", "0"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     ConfigLoader::loadVaccination(f);
                   }});

  cases.push_back({"disease.yaml (trajectory driven)",
                   kTrajectoryDiseaseYaml,
                   {{},
                    {"disease"},
                    {"disease", "settings"},
                    {"disease", "symptom_tags", "0"},
                    {"disease", "trajectories", "0"},
                    {"disease", "trajectories", "0", "stages", "0"},
                    {"disease", "transmission"},
                    {"disease", "transmission", "modes", "0"},
                    {"disease", "immunity"}},
                   [](const ScopedTestFiles&, const std::string& f) {
                     DiseaseLoader::loadFromYAML(f);
                   }});

  cases.push_back(
      {"disease.yaml (stage driven)",
       kStageDrivenDiseaseYaml,
       {{"disease", "outcome_rates_csv"},
        {"disease", "transmission"},
        {"disease", "transmission", "modes", "0"},
        {"disease", "transmission", "modes", "1"},
        {"disease", "transmission", "modes", "1", "fomite_curve"},
        {"disease", "transmission", "modes", "2"},
        {"disease", "transmission", "modes", "3"},
        {"disease", "transmission", "stage_curves", "airborne", "mild"}},
       [](const ScopedTestFiles& files, const std::string& f) {
         writeRatesCsv(files);
         DiseaseLoader::loadFromYAML(f);
       }});

  return cases;
}

}  // namespace

TEST_CASE("every loader names a stray key at every structural level") {
  for (const auto& c : loaderCases()) {
    SUBCASE(c.name.c_str()) { expectEveryLevelRejected(c); }
  }
}

TEST_CASE("a stray key in a selection criterion is named") {
  ScopedTestFiles files{"june_strict_selection"};
  const auto file = files.write("prefs.yaml", R"(
profiles:
  - name: adults
    activity: leisure
    selection:
      - {property: age, operator: ">", value: 17, comparison: strict}
)");
  CHECK_THROWS_WITH_AS(ConfigLoader::loadActivityPreferences(file.string()),
                       doctest::Contains("unknown field 'comparison'"),
                       std::runtime_error);
}

// ---------------------------------------------------------------------------
// The two typos that used to change a run without a word
// ---------------------------------------------------------------------------

TEST_CASE("a misspelt policies_file in config_paths is rejected") {
  ScopedTestFiles files{"june_strict_paths"};
  const auto file = files.write("simulation.yaml", R"(
config_paths:
  disease_file: disease.yaml
  policies_fil: policies.yaml
)");
  CHECK_THROWS_WITH_AS(ConfigLoader::loadSimulation(file.string()),
                       doctest::Contains("unknown field 'policies_fil'"),
                       std::runtime_error);
}

TEST_CASE("a misspelt compliance_rate in a policy is rejected") {
  ScopedTestFiles files{"june_strict_policy"};
  const auto file = files.write("policies.yaml", R"(
policies:
  temporal_policies:
    - name: close
      start_date: "2020-03-09"
      override_activities: [leisure]
      complience_rate: 0.85
)");
  CHECK_THROWS_WITH_AS(loadPoliciesFile(file.string()),
                       doctest::Contains("unknown field 'complience_rate'"),
                       std::runtime_error);
  // The message says which policy it was found in.
  CHECK_THROWS_WITH_AS(loadPoliciesFile(file.string()),
                       doctest::Contains("temporal_policy 'close'"),
                       std::runtime_error);
}

// ---------------------------------------------------------------------------
// Distributions
// ---------------------------------------------------------------------------

namespace {

// A disease whose only trajectory stage takes `dist` (a YAML flow mapping) as
// its completion time.
void loadDiseaseWithDistribution(const std::string& dist) {
  ScopedTestFiles files{"june_strict_dist"};
  const auto file = files.write("disease.yaml", R"(
disease:
  name: dist_test
  symptom_tags:
    - {name: recovered, value: -2}
    - {name: healthy, value: -1}
    - {name: exposed, value: 0}
  trajectories:
    - description: only course
      selection_key: general
      probability: 1.0
      stages:
        - symptom_tag: exposed
          completion_time: )" + dist + "\n");
  DiseaseLoader::loadFromYAML(file.string());
}

}  // namespace

TEST_CASE("each distribution type loads with its documented parameters") {
  CHECK_NOTHROW(loadDiseaseWithDistribution("{type: constant, value: 2.0}"));
  CHECK_NOTHROW(loadDiseaseWithDistribution("{type: normal}"));
  CHECK_NOTHROW(
      loadDiseaseWithDistribution("{type: normal, loc: 3, scale: 1}"));
  CHECK_NOTHROW(loadDiseaseWithDistribution(
      "{type: lognormal, s: 0.5, loc: 0, scale: 2}"));
  CHECK_NOTHROW(loadDiseaseWithDistribution(
      "{type: beta, a: 2, b: 3, loc: 0, scale: 5}"));
  CHECK_NOTHROW(loadDiseaseWithDistribution(
      "{type: exponweib, a: 1.5, c: 2, loc: 0, scale: 4}"));
  CHECK_NOTHROW(
      loadDiseaseWithDistribution("{type: gamma, shape: 2, scale: 3}"));
  CHECK_NOTHROW(
      loadDiseaseWithDistribution("{type: gamma, shape: 2, scale: 3, loc: 1}"));
}

TEST_CASE("a distribution missing a required parameter names it and its type") {
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: lognormal, loc: 0, scale: 2}"),
      doctest::Contains("lognormal distribution is missing 's'"),
      std::runtime_error);
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: constant}"),
      doctest::Contains("constant distribution is missing 'value'"),
      std::runtime_error);
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: beta, a: 2, b: 3, loc: 0}"),
      doctest::Contains("beta distribution is missing 'scale'"),
      std::runtime_error);
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: exponweib, a: 1, c: 2, scale: 4}"),
      doctest::Contains("exponweib distribution is missing 'loc'"),
      std::runtime_error);
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: gamma, scale: 3}"),
      doctest::Contains("gamma distribution is missing 'shape'"),
      std::runtime_error);
}

TEST_CASE("a distribution with an unknown parameter names it and its type") {
  // A normal distribution used to read scal as no scale at all and carry on
  // with a scale of 1.
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: normal, loc: 3, scal: 1}"),
      doctest::Contains("unknown parameter 'scal' for a normal distribution"),
      std::runtime_error);
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: constant, value: 1, scale: 2}"),
      doctest::Contains(
          "unknown parameter 'scale' for a constant distribution"),
      std::runtime_error);
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: gamma, shape: 2, scale: 3, rate: 1}"),
      doctest::Contains("unknown parameter 'rate' for a gamma distribution"),
      std::runtime_error);
}

TEST_CASE("a distribution with an unknown or missing type is rejected") {
  CHECK_THROWS_WITH_AS(loadDiseaseWithDistribution("{type: weibull, scale: 2}"),
                       doctest::Contains("unknown distribution type 'weibull'"),
                       std::runtime_error);
  CHECK_THROWS_WITH_AS(loadDiseaseWithDistribution("{value: 2}"),
                       doctest::Contains("must have a 'type'"),
                       std::runtime_error);
}

TEST_CASE("a non numeric distribution parameter is rejected, not skipped") {
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: constant, value: abc}"),
      doctest::Contains(
          "parameter 'value' of a constant distribution must be a "
          "number"),
      std::runtime_error);
  CHECK_THROWS_WITH_AS(
      loadDiseaseWithDistribution("{type: normal, loc: 3, scale: [1, 2]}"),
      doctest::Contains("parameter 'scale' of a normal distribution must be a "
                        "number"),
      std::runtime_error);
}

TEST_CASE("an unset transmission distribution still means zero") {
  // shape, rate and shift are optional; leaving them out keeps the empty type
  // sampleFromDistribution reads as 0.
  ScopedTestFiles files{"june_strict_unset"};
  const auto file = files.write("disease.yaml", R"(
disease:
  name: unset_test
  symptom_tags:
    - {name: recovered, value: -2}
    - {name: healthy, value: -1}
    - {name: exposed, value: 0}
  transmission:
    mode: Trajectory-Driven
    type: gamma
    max_infectiousness: {type: constant, value: 1.0}
)");
  Disease disease = DiseaseLoader::loadFromYAML(file.string());
  CHECK(disease.getTransmissionParams().shift.type.empty());
}

// ---------------------------------------------------------------------------
// Errors are not swallowed
// ---------------------------------------------------------------------------

TEST_CASE("a malformed performance file throws and names the file") {
  ScopedTestFiles files{"june_strict_perf"};
  const auto file = files.write(
      "performance.yaml", "performance:\n  deterministic_activities: x: y\n");
  CHECK_THROWS_WITH_AS(ConfigLoader::loadPerformance(file.string()),
                       doctest::Contains("performance.yaml"),
                       std::runtime_error);

  const auto bad_value = files.write(
      "bad_value.yaml", "performance:\n  precompute_schedules: sometimes\n");
  CHECK_THROWS_WITH_AS(ConfigLoader::loadPerformance(bad_value.string()),
                       doctest::Contains("bad_value.yaml"), std::runtime_error);

  const auto bad_list =
      files.write("bad_list.yaml", "performance:\n  hybrid_activities: 7\n");
  CHECK_THROWS_WITH_AS(ConfigLoader::loadPerformance(bad_list.string()),
                       doctest::Contains("bad_list.yaml"), std::runtime_error);
}

TEST_CASE("a performance file that is configured but missing throws") {
  ScopedTestFiles files{"june_strict_perf_missing"};
  const auto missing = files.path("not_there.yaml").string();
  CHECK_THROWS_WITH_AS(ConfigLoader::loadPerformance(missing),
                       doctest::Contains("not_there.yaml"), std::runtime_error);
}

TEST_CASE("no performance file configured means the defaults") {
  PerformanceConfig defaults;
  PerformanceConfig loaded = ConfigLoader::loadPerformance("");
  CHECK(loaded.precompute_schedules == defaults.precompute_schedules);
  CHECK(loaded.deterministic_activities == defaults.deterministic_activities);
}

TEST_CASE("a malformed activity preferences file throws and names the file") {
  ScopedTestFiles files{"june_strict_prefs"};
  const auto bad_weight = files.write("bad_weight.yaml", R"(
profiles:
  - name: adults
    activity: leisure
    weights: {grocery: abc}
)");
  CHECK_THROWS_WITH_AS(
      ConfigLoader::loadActivityPreferences(bad_weight.string()),
      doctest::Contains("bad_weight.yaml"), std::runtime_error);

  const auto bad_syntax =
      files.write("bad_syntax.yaml", "profiles:\n  - name: a: b\n");
  CHECK_THROWS_WITH_AS(
      ConfigLoader::loadActivityPreferences(bad_syntax.string()),
      doctest::Contains("bad_syntax.yaml"), std::runtime_error);
}

TEST_CASE(
    "an activity preferences file that is configured but missing throws") {
  ScopedTestFiles files{"june_strict_prefs_missing"};
  const auto missing = files.path("not_there.yaml").string();
  CHECK_THROWS_WITH_AS(ConfigLoader::loadActivityPreferences(missing),
                       doctest::Contains("not_there.yaml"), std::runtime_error);
}

TEST_CASE("no activity preferences file configured means no profiles") {
  CHECK(ConfigLoader::loadActivityPreferences("").profiles.empty());
}

TEST_CASE("a vaccines file that is configured but missing throws") {
  ScopedTestFiles files{"june_strict_vaccines_missing"};
  const auto missing = files.path("not_there.yaml").string();
  CHECK_THROWS_WITH_AS(ConfigLoader::loadVaccination(missing),
                       doctest::Contains("not_there.yaml"), std::runtime_error);
}

TEST_CASE("no vaccines file configured means vaccination is off") {
  CHECK_FALSE(ConfigLoader::loadVaccination("").enabled);
}

TEST_CASE("default_seed_strength is read and checked") {
  ScopedTestFiles files{"june_strict_seed_strength"};
  const auto set = files.write("set.yaml", R"(
global_parameters:
  default_seed_strength: 2.5
infection_seeds:
  - name: a
    type: uniform
    date: "2020-03-02 08:00"
    parameters:
      cases_per_capita_multiplier: 1.0
  - name: b
    type: uniform
    date: "2020-03-03 08:00"
    parameters:
      seed_strength: 0.5
      cases_per_capita_multiplier: 1.0
)");
  const auto config = InfectionSeedConfigLoader::loadFromFile(set.string());
  CHECK(config.global_params.default_seed_strength == doctest::Approx(2.5));
  REQUIRE(config.seeds.size() == 2);
  CHECK(config.seeds[0].seed_strength == doctest::Approx(2.5));
  CHECK(config.seeds[1].seed_strength == doctest::Approx(0.5));

  const auto negative = files.write("negative.yaml", R"(
global_parameters:
  default_seed_strength: -1.0
)");
  CHECK_THROWS_WITH_AS(
      InfectionSeedConfigLoader::loadFromFile(negative.string()),
      doctest::Contains("default_seed_strength"), std::runtime_error);
}

// ---------------------------------------------------------------------------
// Every shipped config loads
// ---------------------------------------------------------------------------

namespace {

// Load the files of one scenario that need no world. Files whose data lives
// outside the repository (a CSV on a data drive) are pointed at a stub of the
// same shape, so the keys are still checked.
void loadScenario(const std::string& simulation_file) {
  INFO("scenario " << simulation_file);
  const SimulationConfig sim = ConfigLoader::loadSimulation(simulation_file);

  CHECK_NOTHROW(ConfigLoader::loadSchedule(sim.schedules_file));
  CHECK_NOTHROW(ConfigLoader::loadContactMatrices(sim.contact_matrices_file));
  CHECK_NOTHROW(ConfigLoader::loadPerformance(sim.performance_file));
  CHECK_NOTHROW(ConfigLoader::loadParallel(sim.parallel_file));
  CHECK_NOTHROW(
      ConfigLoader::loadActivityPreferences(sim.activity_preferences_file));
  CHECK_NOTHROW(ConfigLoader::loadVaccination(sim.vaccines_file));
  CHECK_NOTHROW(
      ConfigLoader::loadCoordinatedEncounters(sim.coordinated_encounters_file));

  if (!sim.on_the_fly_venues_file.empty()) {
    CHECK_NOTHROW(OnTheFlyVenueAllocator(sim.on_the_fly_venues_file));
  }
  if (!sim.policies_file.empty()) {
    CHECK_NOTHROW(loadPoliciesFile(sim.policies_file));
  }
  if (!sim.infection_seeds_file.empty()) {
    // A bulk seed table that lives on a data drive is left out of a copy, so
    // the keys of the file itself are still checked.
    ScopedTestFiles seed_files{"june_strict_shipped_seeds"};
    std::string seeds_file = sim.infection_seeds_file;
    YAML::Node seeds = YAML::LoadFile(seeds_file);
    if (seeds["bulk_csv"] && !fs::exists(seeds["bulk_csv"].as<std::string>())) {
      seeds.remove("bulk_csv");
      seeds_file = seed_files.write("seeds.yaml", YAML::Dump(seeds)).string();
    }
    CHECK_NOTHROW(InfectionSeedConfigLoader::loadFromFile(seeds_file));
  }

  // The disease file names an outcome table relative to itself. When the table
  // is not in the repository, load a copy that points at an empty one.
  ScopedTestFiles files{"june_strict_shipped_disease"};
  YAML::Node root = YAML::LoadFile(sim.disease_file);
  YAML::Node csv = root["disease"]["outcome_rates_csv"];
  std::string disease_file = sim.disease_file;
  if (csv) {
    const std::string rel =
        csv.IsScalar() ? csv.as<std::string>() : csv["file"].as<std::string>();
    const fs::path table = fs::path(sim.disease_file).parent_path() / rel;
    if (!fs::exists(table)) {
      files.write("rates.csv", "mild\n");
      if (csv.IsScalar()) {
        root["disease"]["outcome_rates_csv"] = "rates.csv";
      } else {
        csv["file"] = "rates.csv";
      }
      disease_file = files.write("disease.yaml", YAML::Dump(root)).string();
    }
  }
  CHECK_NOTHROW(DiseaseLoader::loadFromYAML(disease_file));
}

}  // namespace

TEST_CASE("every shipped scenario loads without an unknown key") {
  const std::vector<std::string> scenarios = {
      "configs/config_1911/simulation.yaml",
      "configs/config_2021/simulation.yaml",
      "configs/config_plague/simulation.yaml",
      "tests/configs/fomite_simulation.yaml",
      "tests/configs/sex_age_bins/simulation.yaml",
  };
  for (const auto& scenario : scenarios) {
    SUBCASE(scenario.c_str()) { loadScenario(scenario); }
  }
}

TEST_CASE("every standalone test config loads without an unknown key") {
  CHECK_NOTHROW(
      ConfigLoader::loadContactMatrices("tests/configs/contact_matrices.yaml"));
  CHECK_NOTHROW(ConfigLoader::loadContactMatrices(
      "tests/configs/contact_matrices_modes_default.yaml"));
  CHECK_NOTHROW(ConfigLoader::loadContactMatrices(
      "tests/configs/fomite_contact_matrices.yaml"));
  CHECK_NOTHROW(ConfigLoader::loadContactMatrices(
      "tests/configs/romantic_regression.yaml"));
  CHECK_NOTHROW(ConfigLoader::loadContactMatrices(
      "tests/configs/contact_matrices_reversed_modes.yaml"));
  CHECK_NOTHROW(ConfigLoader::loadVaccination("tests/configs/vaccines.yaml"));
  CHECK_NOTHROW(
      DiseaseLoader::loadFromYAML("tests/configs/simulator_disease.yaml"));
  CHECK_NOTHROW(InfectionSeedConfigLoader::loadFromFile(
      "tests/configs/simulator_infection_seeds.yaml"));
}

namespace {

std::string replaceOnce(std::string text, const std::string& from,
                        const std::string& to) {
  const size_t at = text.find(from);
  REQUIRE_MESSAGE(at != std::string::npos, "test text lacks: " << from);
  text.replace(at, from.size(), to);
  return text;
}

void loadDiseaseText(const std::string& yaml) {
  ScopedTestFiles files;
  writeRatesCsv(files);
  const auto file = files.write("disease.yaml", yaml);
  DiseaseLoader::loadFromYAML(file.string());
}

}  // namespace

TEST_CASE("disease names must be symptom tags and transmission modes") {
  const std::string stage = kStageDrivenDiseaseYaml;
  CHECK_NOTHROW(loadDiseaseText(stage));

  SUBCASE("a stage_curves stage that is not a symptom tag") {
    CHECK_THROWS_WITH_AS(
        loadDiseaseText(replaceOnce(stage, "        mild: {type: gamma",
                                    "        mld: {type: gamma")),
        doctest::Contains("stage_curves.airborne names stage 'mld', which is "
                          "not a symptom tag. Tags: recovered, healthy, "
                          "exposed, mild"),
        std::runtime_error);
  }
  SUBCASE("a stage_curves mode that is not a transmission mode") {
    CHECK_THROWS_WITH_AS(
        loadDiseaseText(replaceOnce(stage, "      airborne:\n        mild",
                                    "      airbone:\n        mild")),
        doctest::Contains("stage_curves names mode 'airbone', which is not a "
                          "transmission mode. Modes: airborne, surface"),
        std::runtime_error);
  }
  SUBCASE("a deposition stage that is not a symptom tag") {
    CHECK_THROWS_WITH_AS(
        loadDiseaseText(replaceOnce(stage, "          mild: {type: constant",
                                    "          mld: {type: constant")),
        doctest::Contains("fomite mode 'surface' deposition_stages names "
                          "stage 'mld', which is not a symptom tag"),
        std::runtime_error);
  }
  SUBCASE("a transmission mode that is neither of the two") {
    CHECK_THROWS_WITH_AS(
        loadDiseaseText(replaceOnce(stage, "mode: Stage-Driven",
                                    "mode: Stage-driven")),
        doctest::Contains("transmission.mode 'Stage-driven' is not a mode. "
                          "Modes: Stage-Driven, Trajectory-Driven"),
        std::runtime_error);
  }
  SUBCASE("a trajectory stage without a symptom_tag") {
    CHECK_THROWS_WITH_AS(
        loadDiseaseText(replaceOnce(
            stage, "        - symptom_tag: mild\n          completion_time",
            "        - completion_time")),
        doctest::Contains("is missing 'symptom_tag'"), std::runtime_error);
  }
  SUBCASE("a trajectory without stages") {
    CHECK_THROWS_WITH_AS(
        loadDiseaseText(replaceOnce(
            stage,
            "      stages:\n"
            "        - symptom_tag: exposed\n"
            "          completion_time: {type: constant, value: 2.0}\n"
            "        - symptom_tag: mild\n"
            "          completion_time: {type: constant, value: 3.0}\n",
            "")),
        doctest::Contains("is missing 'stages'"), std::runtime_error);
  }
}
