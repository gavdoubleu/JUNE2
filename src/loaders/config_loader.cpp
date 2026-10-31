#include "loaders/config_loader.h"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "loaders/config_loader_detail.h"
#include "utils/filtered_csv.h"
#include "utils/filtering.h"

namespace june {

namespace {

using ::june::config_detail::parseSelectionCriteria;
using ::june::config_detail::rejectUnknownKeys;
using ::june::config_detail::rejectUnknownKeysAmong;

// Load a YAML file whose path was written down in the config, so a missing or
// unreadable file is an error that names it.
YAML::Node loadConfiguredFile(const std::string& filename,
                              const std::string& what) {
  if (!std::filesystem::exists(filename)) {
    throw std::runtime_error(what + " file not found: " + filename + ".");
  }
  try {
    return YAML::LoadFile(filename);
  } catch (const YAML::Exception& e) {
    throw std::runtime_error(what + " file " + filename + ": " + e.what());
  }
}

// Parse a YAML sequence of `{name, start, end, allowed_activities, ...}`
// entries into a vector of TimeSlot. Shared by the per-schedule-type loop
// in loadSchedule (both flat_slots and per-day-type slot lists).
std::vector<TimeSlot> parseSlotList(const YAML::Node& slot_list_node,
                                    const std::string& where,
                                    const std::string& filename) {
  std::vector<TimeSlot> slots;
  for (const auto& slot : slot_list_node) {
    rejectUnknownKeys(slot,
                      {"name", "start", "end", "allowed_activities",
                       "coordinated_only_activities", "specified_activity",
                       "hop_on_activity"},
                      filename + ": ", " in a slot of '" + where + "'");
    TimeSlot ts;
    ts.name = slot["name"].as<std::string>();
    ts.start = slot["start"].as<std::string>();
    ts.end = slot["end"].as<std::string>();
    ts.allowed_activities =
        slot["allowed_activities"].as<std::vector<std::string>>();

    if (slot["coordinated_only_activities"])
      ts.coordinated_only_activities =
          slot["coordinated_only_activities"].as<std::vector<std::string>>();

    if (slot["specified_activity"]) {
      rejectUnknownKeys(
          slot["specified_activity"], {"type", "index", "venue_type"},
          filename + ": ",
          " in 'specified_activity' of a slot of '" + where + "'");
      SpecifiedActivity spec_act;
      spec_act.type = slot["specified_activity"]["type"].as<std::string>();
      spec_act.index = slot["specified_activity"]["index"].as<int>();
      if (slot["specified_activity"]["venue_type"]) {
        spec_act.venue_type =
            slot["specified_activity"]["venue_type"].as<std::string>();
      }
      ts.specified_activity = spec_act;
    }

    if (slot["hop_on_activity"]) {
      for (const auto& hop_kv : slot["hop_on_activity"]) {
        std::string act_name = hop_kv.first.as<std::string>();
        if (hop_kv.second.IsScalar()) {
          // Direct hop: activity_name → schedule_name
          ts.hop_on_activity[act_name] = hop_kv.second.as<std::string>();
        } else if (hop_kv.second.IsMap()) {
          // Property-dispatched hop: schedule name resolved at runtime
          rejectUnknownKeys(
              hop_kv.second, {"property", "schedule_template"}, filename + ": ",
              " in hop_on_activity '" + act_name + "' of '" + where + "'");
          TimeSlot::PropertyDispatchHop dispatch;
          dispatch.property_name = hop_kv.second["property"].as<std::string>();
          dispatch.schedule_name_template =
              hop_kv.second["schedule_template"].as<std::string>();
          ts.property_hop_dispatch[act_name] = std::move(dispatch);
        }
      }
    }

    slots.push_back(ts);
  }
  return slots;
}

// Read the required `config_paths:` block onto the SimulationConfig file
// fields.
void parseSimulationConfigPaths(const YAML::Node& paths,
                                SimulationConfig& config,
                                const std::string& filename) {
  rejectUnknownKeys(
      paths,
      {"disease_file", "contact_matrices_file", "schedules_file",
       "vaccines_file", "activity_preferences_file",
       "coordinated_encounters_file", "performance_file", "parallel_file",
       "policies_file", "infection_seeds_file", "compartmental_model_sidecar",
       "calendar_events_file", "calendar_event_catchment_rules_file",
       "on_the_fly_venues_file"},
      filename + ": ", " in 'config_paths'");
  if (paths["disease_file"])
    config.disease_file = paths["disease_file"].as<std::string>();
  if (paths["contact_matrices_file"])
    config.contact_matrices_file =
        paths["contact_matrices_file"].as<std::string>();
  if (paths["schedules_file"])
    config.schedules_file = paths["schedules_file"].as<std::string>();
  if (paths["vaccines_file"])
    config.vaccines_file = paths["vaccines_file"].as<std::string>();
  if (paths["activity_preferences_file"])
    config.activity_preferences_file =
        paths["activity_preferences_file"].as<std::string>();
  if (paths["coordinated_encounters_file"])
    config.coordinated_encounters_file =
        paths["coordinated_encounters_file"].as<std::string>();
  if (paths["performance_file"])
    config.performance_file = paths["performance_file"].as<std::string>();
  if (paths["parallel_file"])
    config.parallel_file = paths["parallel_file"].as<std::string>();
  if (paths["policies_file"])
    config.policies_file = paths["policies_file"].as<std::string>();
  if (paths["infection_seeds_file"])
    config.infection_seeds_file =
        paths["infection_seeds_file"].as<std::string>();
  if (paths["compartmental_model_sidecar"])
    config.compartmental_model_sidecar =
        paths["compartmental_model_sidecar"].as<std::string>();
  if (paths["calendar_events_file"])
    config.calendar_events_file =
        paths["calendar_events_file"].as<std::string>();
  if (paths["calendar_event_catchment_rules_file"])
    config.calendar_event_catchment_rules_file =
        paths["calendar_event_catchment_rules_file"].as<std::string>();
  if (paths["on_the_fly_venues_file"])
    config.on_the_fly_venues_file =
        paths["on_the_fly_venues_file"].as<std::string>();
}

// Set `mode` from output[key] when the key is present. The writers only know
// "none", "all" and "infected_only", and any other spelling would save no
// table and say nothing, so it throws here instead.
void readLookupMode(const YAML::Node& output, const std::string& key,
                    std::string& mode) {
  if (!output[key]) return;
  const std::string value = output[key].as<std::string>();
  if (value != "none" && value != "all" && value != "infected_only") {
    throw std::runtime_error("output." + key + " is '" + value +
                             "'; valid values are 'none', 'all' and "
                             "'infected_only'.");
  }
  mode = value;
}

// Read the `output:` block onto the SimulationConfig output fields.
// `stats_interval_days` is required when the block is present; everything
// else is optional.
void parseSimulationOutput(const YAML::Node& output, SimulationConfig& config,
                           const std::string& filename) {
  rejectUnknownKeys(
      output,
      {"stats_interval_days", "save_full_person_details",
       "save_person_activities", "compression_level", "flush_interval_days",
       "max_event_buffer_size", "save_population_summary",
       "save_coordinated_encounters", "summary_properties"},
      filename + ": ", " in 'output'");
  config.stats_interval_days = output["stats_interval_days"].as<int>();

  readLookupMode(output, "save_full_person_details",
                 config.save_full_person_details);
  readLookupMode(output, "save_person_activities",
                 config.save_person_activities);
  if (output["compression_level"]) {
    config.compression_level = output["compression_level"].as<int>();
  }
  if (output["flush_interval_days"]) {
    config.flush_interval_days = output["flush_interval_days"].as<int>();
  }
  if (output["max_event_buffer_size"]) {
    config.max_event_buffer_size = output["max_event_buffer_size"].as<int>();
  }
  if (output["save_population_summary"]) {
    config.save_population_summary =
        output["save_population_summary"].as<bool>();
  }
  if (output["save_coordinated_encounters"]) {
    config.save_coordinated_encounters =
        output["save_coordinated_encounters"].as<bool>();
  }
  if (output["summary_properties"]) {
    config.summary_properties =
        output["summary_properties"].as<std::vector<std::string>>();
    config.summary_properties_explicit = true;
  }
}

// Read the `partial_presence.enabled_venue_types:` map onto the per-venue
// target_group_size table.
void parseSimulationPartialPresence(const YAML::Node& pp,
                                    SimulationConfig& config,
                                    const std::string& filename) {
  rejectUnknownKeys(pp, {"enabled_venue_types"}, filename + ": ",
                    " in 'partial_presence'");
  if (pp["enabled_venue_types"] && pp["enabled_venue_types"].IsMap()) {
    for (auto it = pp["enabled_venue_types"].begin();
         it != pp["enabled_venue_types"].end(); ++it) {
      std::string name = it->first.as<std::string>();
      rejectUnknownKeys(it->second, {"target_group_size"}, filename + ": ",
                        " in partial_presence venue type '" + name + "'");
      int tgs = it->second["target_group_size"]
                    ? it->second["target_group_size"].as<int>()
                    : 0;
      config.partial_presence.target_group_size_by_name[name] = tgs;
    }
  }
}

// Parse a single contact-matrix YAML node into a ContactMatrix. Multiplies
// the `contacts` rows by the optional `beta` scaling factor at parse time
// (matching the historical behaviour of the inline lambda this replaced).
ContactMatrix parseContactMatrix(const YAML::Node& matrix_node,
                                 const std::string& filename,
                                 const std::string& within) {
  rejectUnknownKeys(matrix_node, {"bins", "beta", "contacts"}, filename + ": ",
                    within);
  ContactMatrix cm;
  if (matrix_node["bins"]) {
    cm.bins = matrix_node["bins"].as<std::vector<std::string>>();
  }
  double beta = 1.0;
  if (matrix_node["beta"]) {
    beta = matrix_node["beta"].as<double>();
  }
  if (matrix_node["contacts"]) {
    for (const auto& row : matrix_node["contacts"]) {
      std::vector<double> row_vec;
      for (const auto& val : row) row_vec.push_back(val.as<double>() * beta);
      cm.contacts.push_back(row_vec);
    }
  }
  return cm;
}

// Read the `contact_matrices:` block: for each venue type, dispatch on
// whether the entry declares `modes:` (multi-mode) or is a single contact
// matrix. Records mode names in YAML insertion order on config.mode_names.
void parseContactMatricesList(const YAML::Node& cm_node,
                              ContactMatrixConfig& config,
                              const std::string& filename) {
  std::vector<std::string> seen_modes_ordered;
  auto addMode = [&seen_modes_ordered](const std::string& m) {
    for (const auto& s : seen_modes_ordered)
      if (s == m) return;
    seen_modes_ordered.push_back(m);
  };

  for (const auto& matrix_kv : cm_node) {
    std::string venue_type = matrix_kv.first.as<std::string>();
    const auto& matrix_node = matrix_kv.second;

    if (matrix_node["modes"]) {
      rejectUnknownKeys(matrix_node, {"modes"}, filename + ": ",
                        " in venue type '" + venue_type + "'");
      // Multi-mode: parse each mode entry. Nothing is copied into the flat
      // `matrices` map — a type that declares modes is described by those
      // modes, and mirroring an arbitrary one of them into a mode-free slot
      // is how a mode-specific matrix used to end up answering mode-free
      // questions.
      for (const auto& mode_kv : matrix_node["modes"]) {
        std::string mode_name = mode_kv.first.as<std::string>();
        addMode(mode_name);
        config.mode_matrices[venue_type][mode_name] = parseContactMatrix(
            mode_kv.second, filename,
            " in mode '" + mode_name + "' of venue type '" + venue_type + "'");
      }
    } else {
      // Single-mode fallback: store as "default" mode and in flat map
      addMode("default");
      ContactMatrix cm = parseContactMatrix(
          matrix_node, filename, " in venue type '" + venue_type + "'");
      config.mode_matrices[venue_type]["default"] = cm;
      config.matrices[venue_type] = cm;
    }
  }

  if (!seen_modes_ordered.empty()) {
    config.mode_names = seen_modes_ordered;
  }
}

// Read the required `parallel.partitioning:` block. Throws if any of the
// three mandatory keys (level / centroids_file / adjacency_file) is missing,
// embedding the source `filename` in the error so the user knows which
// parallel.yaml triggered it. The optional metis.imbalance_tolerance is
// applied if present.
void parseParallelPartitioning(const YAML::Node& part,
                               const std::string& filename,
                               ParallelConfig& config) {
  rejectUnknownKeys(part,
                    {"level", "centroids_file", "adjacency_file", "metis"},
                    filename + ": ", " in 'parallel.partitioning'");
  if (!part["level"]) {
    throw std::runtime_error(
        "Parallel config '" + filename +
        "' is missing required key 'parallel.partitioning.level' "
        "(e.g. 'SGU', 'MGU', 'LGU').");
  }
  config.partition_level = part["level"].as<std::string>();

  if (!part["centroids_file"]) {
    throw std::runtime_error(
        "Parallel config '" + filename +
        "' is missing required key 'parallel.partitioning.centroids_file'.");
  }
  config.centroids_file = part["centroids_file"].as<std::string>();

  if (!part["adjacency_file"]) {
    throw std::runtime_error(
        "Parallel config '" + filename +
        "' is missing required key 'parallel.partitioning.adjacency_file'.");
  }
  config.adjacency_file = part["adjacency_file"].as<std::string>();

  if (part["metis"]) {
    rejectUnknownKeys(part["metis"], {"imbalance_tolerance"}, filename + ": ",
                      " in 'parallel.partitioning.metis'");
  }
  if (part["metis"] && part["metis"]["imbalance_tolerance"]) {
    config.metis_imbalance_tolerance =
        part["metis"]["imbalance_tolerance"].as<double>();
  }
}

// Read the `checkpoint:` block. Cadence is mutually exclusive: on_dates
// (non-null, non-empty) takes precedence over every_n_days. A null YAML
// value leaves the corresponding optional empty.
void parseSimulationCheckpoint(const YAML::Node& cp, SimulationConfig& config,
                               const std::string& filename) {
  rejectUnknownKeys(
      cp, {"enabled", "output_dir", "every_n_days", "on_dates", "keep_last"},
      filename + ": ", " in 'checkpoint'");
  if (cp["enabled"]) {
    config.checkpoint.enabled = cp["enabled"].as<bool>();
  }
  if (cp["output_dir"]) {
    config.checkpoint.output_dir = cp["output_dir"].as<std::string>();
  }
  if (cp["every_n_days"] && !cp["every_n_days"].IsNull()) {
    config.checkpoint.every_n_days = cp["every_n_days"].as<int>();
  }
  if (cp["on_dates"] && !cp["on_dates"].IsNull()) {
    config.checkpoint.on_dates = cp["on_dates"].as<std::vector<std::string>>();
  }
  if (cp["keep_last"]) {
    config.checkpoint.keep_last = cp["keep_last"].as<int>();
  }
}

// Parse one entry from `schedule_types:` (a named schedule). Reads the
// flag/scalar fields, the optional `selection` block (via
// parseSelectionCriteria), the flat_slots / per-day-type slot lists (via
// parseSlotList), the per-schedule activity overrides, and the
// participation-rate table. `day_type_names` lists the day-type keys to
// look up on the YAML node.
ScheduleType parseScheduleType(const std::string& name,
                               const YAML::Node& type_node,
                               const std::vector<std::string>& day_type_names,
                               const std::string& filename) {
  std::vector<std::string> type_fields = {
      "priority",          "is_temporary", "return_schedule",
      "flat_slots",        "selection",    "force_hybrid_activities",
      "linked_activities", "participation"};
  type_fields.insert(type_fields.end(), day_type_names.begin(),
                     day_type_names.end());
  rejectUnknownKeysAmong(type_node, type_fields, filename + ": ",
                         " in schedule type '" + name + "'");
  ScheduleType sched_type;
  sched_type.name = name;

  if (type_node["priority"]) {
    sched_type.priority = type_node["priority"].as<int>();
  }

  if (type_node["is_temporary"]) {
    sched_type.is_temporary = type_node["is_temporary"].as<bool>();
  }
  if (type_node["return_schedule"]) {
    sched_type.return_schedule = type_node["return_schedule"].as<std::string>();
  }
  if (type_node["flat_slots"]) {
    sched_type.flat_slots =
        parseSlotList(type_node["flat_slots"], name, filename);
  }

  if (type_node["selection"].IsDefined()) {
    parseSelectionCriteria(type_node["selection"],
                           sched_type.selection_criteria);
  }

  for (const auto& dt_name : day_type_names) {
    if (type_node[dt_name]) {
      sched_type.slots_by_day_type[dt_name] =
          parseSlotList(type_node[dt_name], name, filename);
    }
  }

  if (type_node["force_hybrid_activities"]) {
    sched_type.force_hybrid_activities =
        type_node["force_hybrid_activities"].as<std::vector<std::string>>();
  }

  if (type_node["linked_activities"]) {
    sched_type.linked_activities =
        type_node["linked_activities"].as<std::vector<std::string>>();
  }

  if (type_node["participation"]) {
    for (const auto& kv : type_node["participation"]) {
      std::string activity = kv.first.as<std::string>();
      rejectUnknownKeysAmong(kv.second, day_type_names, filename + ": ",
                             " in participation of '" + activity +
                                 "' in schedule type '" + name + "'");
      for (const auto& dt_name : day_type_names) {
        if (kv.second[dt_name]) {
          sched_type.participation_by_day_type[dt_name][activity] =
              kv.second[dt_name].as<double>();
        }
      }
    }
  }

  return sched_type;
}

}  // namespace

namespace config_detail {

PropertyValue parseCriterionSequenceValue(const YAML::Node& value_node,
                                          const std::string& property_path) {
  if (SelectionCriterion::comparesAgainstUnitNames(property_path)) {
    try {
      return value_node.as<std::vector<std::string>>();
    } catch (const YAML::Exception&) {
      throw std::runtime_error("selection criterion '" + property_path +
                               "': list value must be geographical unit names");
    }
  }

  try {
    return value_node.as<std::vector<int32_t>>();
  } catch (const YAML::Exception&) {
    throw std::runtime_error(
        "selection criterion '" + property_path +
        "': list value must be whole numbers (only 'geo_unit.<LEVEL>' "
        "compares against a list of names)");
  }
}

// Parse a YAML sequence of `{property, operator, value}` entries into a vector
// of SelectionCriterion. Scalar values use bool -> int -> double -> string;
// sequence values are validated by the local parseCriterionSequenceValue
// helper.
// Declared in loaders/config_loader_detail.h so the vaccination-loader TU can
// reuse it.
void parseSelectionCriteria(const YAML::Node& selection_node,
                            std::vector<SelectionCriterion>& out) {
  if (!selection_node.IsDefined()) return;
  if (!selection_node.IsSequence()) {
    throw MalformedSelectionError(
        "selection must be a sequence of {property, operator, value} "
        "entries; use [] for a catch-all selection");
  }

  for (const auto& criterion_node : selection_node) {
    rejectUnknownKeys(criterion_node, {"property", "operator", "value"},
                      "selection criterion: ");
    SelectionCriterion criterion;
    criterion.property_path = criterion_node["property"].as<std::string>();
    criterion.operator_type = criterion_node["operator"].as<std::string>();

    const auto& value_node = criterion_node["value"];
    if (value_node.IsSequence()) {
      criterion.value =
          parseCriterionSequenceValue(value_node, criterion.property_path);
    } else if (value_node.IsScalar()) {
      const std::string scalar = value_node.as<std::string>();
      if (scalar == "true" || scalar == "false") {
        criterion.value = value_node.as<bool>() ? 1 : 0;
      } else {
        try {
          criterion.value = value_node.as<int>();
        } catch (...) {
          try {
            criterion.value = value_node.as<double>();
          } catch (...) {
            criterion.value = scalar;
          }
        }
      }
    }
    out.push_back(std::move(criterion));
  }
}

void rejectUnknownKeys(const YAML::Node& node,
                       std::initializer_list<const char*> fields,
                       const std::string& where, const std::string& within) {
  rejectUnknownKeysAmong(node,
                         std::vector<std::string>(fields.begin(), fields.end()),
                         where, within);
}

void rejectUnknownKeysAmong(const YAML::Node& node,
                            const std::vector<std::string>& fields,
                            const std::string& where,
                            const std::string& within) {
  if (node.IsDefined() && !node.IsNull() && !node.IsMap()) {
    throw std::runtime_error(where + "expected a mapping of fields" + within +
                             ".");
  }
  for (const auto& kv : node) {
    const std::string key = kv.first.as<std::string>();
    if (std::find(fields.begin(), fields.end(), key) != fields.end()) {
      continue;
    }
    std::string known;
    for (const auto& field : fields)
      known += (known.empty() ? "" : ", ") + field;
    throw std::runtime_error(where + "unknown field '" + key + "'" + within +
                             ". Fields are " + known + ".");
  }
}

}  // namespace config_detail

Config ConfigLoader::loadAll(const std::string& simulation_file) {
  Config config;

  // 1. Load the master simulation config first
  config.simulation = loadSimulation(simulation_file);

  // 2. Load all other configs using paths specified in the simulation config
  config.schedule = loadSchedule(config.simulation.schedules_file);
  config.contact_matrices =
      loadContactMatrices(config.simulation.contact_matrices_file);
  config.performance = loadPerformance(config.simulation.performance_file);
  config.parallel = loadParallel(config.simulation.parallel_file);
  config.activity_preferences =
      loadActivityPreferences(config.simulation.activity_preferences_file);
  config.vaccination = loadVaccination(config.simulation.vaccines_file);
  config.coordinated_encounters =
      loadCoordinatedEncounters(config.simulation.coordinated_encounters_file);

  return config;
}

SimulationConfig ConfigLoader::loadSimulation(const std::string& filename) {
  YAML::Node root = YAML::LoadFile(filename);
  SimulationConfig config;
  rejectUnknownKeys(root,
                    {"time", "config_paths", "random_seed", "verbose", "output",
                     "partial_presence", "checkpoint"},
                    filename + ": ");

  if (root["time"]) {
    auto time = root["time"];
    rejectUnknownKeys(time, {"start_date", "end_date"}, filename + ": ",
                      " in 'time'");
    config.start_date = time["start_date"].as<std::string>();
    config.end_date = time["end_date"].as<std::string>();
  }

  if (!root["config_paths"] || !root["config_paths"].IsMap()) {
    throw std::runtime_error(
        "simulation.yaml requires a top-level 'config_paths:' mapping");
  }
  parseSimulationConfigPaths(root["config_paths"], config, filename);

  if (root["random_seed"]) {
    config.random_seed = root["random_seed"].as<unsigned int>();
  }
  if (root["verbose"]) {
    config.verbose = root["verbose"].as<bool>();
  }

  if (root["output"]) {
    parseSimulationOutput(root["output"], config, filename);
  }

  if (root["partial_presence"]) {
    parseSimulationPartialPresence(root["partial_presence"], config, filename);
  }

  if (root["checkpoint"]) {
    parseSimulationCheckpoint(root["checkpoint"], config, filename);
  }

  return config;
}

ScheduleConfig ConfigLoader::loadSchedule(const std::string& filename) {
  YAML::Node root = YAML::LoadFile(filename);
  ScheduleConfig config;
  rejectUnknownKeys(root,
                    {"day_type_cycle", "schedule_types", "schedule_csv",
                     "default_schedule_type"},
                    filename + ": ");

  // Day type cycle
  if (root["day_type_cycle"]) {
    config.day_type_cycle =
        root["day_type_cycle"].as<std::vector<std::string>>();
  }

  // Derive ordered unique day type names from cycle
  for (const auto& dt : config.day_type_cycle) {
    if (std::find(config.day_type_names.begin(), config.day_type_names.end(),
                  dt) == config.day_type_names.end()) {
      config.day_type_names.push_back(dt);
    }
  }

  // Schedule types
  if (root["schedule_types"]) {
    // Optional CSV for probabilistic schedule assignment
    if (root["schedule_csv"]) {
      config.csv_path = root["schedule_csv"].as<std::string>();
    }

    // Load default schedule type
    if (root["default_schedule_type"]) {
      config.default_schedule_type =
          root["default_schedule_type"].as<std::string>();
    }

    // Parse each schedule type
    for (const auto& type_kv : root["schedule_types"]) {
      config.schedule_types.push_back(
          parseScheduleType(type_kv.first.as<std::string>(), type_kv.second,
                            config.day_type_names, filename));
    }

    // Sort schedule types by priority (highest first)
    std::sort(config.schedule_types.begin(), config.schedule_types.end(),
              [](const ScheduleType& a, const ScheduleType& b) {
                return a.priority > b.priority;
              });
  }

  return config;
}

ContactMatrixConfig ConfigLoader::loadContactMatrices(
    const std::string& filename) {
  YAML::Node root = YAML::LoadFile(filename);
  ContactMatrixConfig config;
  rejectUnknownKeys(
      root,
      {"contact_matrices", "allow_default_matrix", "default_contacts_matrix"},
      filename + ": ");

  if (root["contact_matrices"]) {
    parseContactMatricesList(root["contact_matrices"], config, filename);
  }

  if (root["allow_default_matrix"]) {
    config.allow_default_matrix = root["allow_default_matrix"].as<bool>();
  }

  if (!root["default_contacts_matrix"]) {
    throw std::runtime_error(
        "Contact matrix config '" + filename +
        "' is missing required key 'default_contacts_matrix'.");
  }
  const auto& default_node = root["default_contacts_matrix"];

  if (default_node["modes"]) {
    rejectUnknownKeys(default_node, {"modes"}, filename + ": ",
                      " in 'default_contacts_matrix'");
    std::unordered_map<std::string, ContactMatrix> mode_matrices;
    for (const auto& mode_kv : default_node["modes"]) {
      const std::string mode_name = mode_kv.first.as<std::string>();
      mode_matrices[mode_name] = parseContactMatrix(
          mode_kv.second, filename,
          " in mode '" + mode_name + "' of 'default_contacts_matrix'");
    }

    std::vector<std::string> missing_modes;
    for (const auto& mode_name : config.mode_names) {
      if (!mode_matrices.count(mode_name)) missing_modes.push_back(mode_name);
    }
    if (!missing_modes.empty()) {
      std::string missing_list;
      for (size_t i = 0; i < missing_modes.size(); ++i) {
        if (i) missing_list += ", ";
        missing_list += missing_modes[i];
      }
      throw std::runtime_error(
          "Contact matrix config '" + filename +
          "': 'default_contacts_matrix' is missing mode(s) [" + missing_list +
          "] present in 'contact_matrices'.");
    }

    config.default_mode_matrices = std::move(mode_matrices);
  } else {
    config.default_matrix = parseContactMatrix(default_node, filename,
                                               " in 'default_contacts_matrix'");
  }

  return config;
}

PerformanceConfig ConfigLoader::loadPerformance(const std::string& filename) {
  PerformanceConfig config;

  // No performance_file in config_paths means the defaults. A path that was
  // written down must load.
  if (filename.empty()) return config;

  YAML::Node root = loadConfiguredFile(filename, "Performance");
  try {
    rejectUnknownKeys(root, {"performance"}, filename + ": ");

    if (root["performance"]) {
      const auto& perf = root["performance"];
      rejectUnknownKeys(perf,
                        {"precompute_schedules", "deterministic_activities",
                         "hybrid_activities", "stochastic_activities"},
                        filename + ": ", " in 'performance'");

      if (perf["precompute_schedules"]) {
        config.precompute_schedules = perf["precompute_schedules"].as<bool>();
      }

      if (perf["deterministic_activities"]) {
        config.deterministic_activities =
            perf["deterministic_activities"].as<std::vector<std::string>>();
      }

      if (perf["hybrid_activities"]) {
        config.hybrid_activities =
            perf["hybrid_activities"].as<std::vector<std::string>>();
      }

      if (perf["stochastic_activities"]) {
        config.stochastic_activities =
            perf["stochastic_activities"].as<std::vector<std::string>>();
      }
    }
  } catch (const YAML::Exception& e) {
    throw std::runtime_error("Performance file " + filename + ": " + e.what());
  }

  return config;
}

ParallelConfig ConfigLoader::loadParallel(const std::string& filename) {
  ParallelConfig config;
  YAML::Node root = YAML::LoadFile(filename);

  rejectUnknownKeys(root, {"parallel"}, filename + ": ");
  if (!root["parallel"]) {
    throw std::runtime_error("Parallel config file '" + filename +
                             "' is missing the top-level 'parallel:' block.");
  }
  YAML::Node parallel = root["parallel"];
  rejectUnknownKeys(parallel, {"enabled", "partitioning", "chunked_loading"},
                    filename + ": ", " in 'parallel'");

  if (parallel["enabled"]) {
    config.enabled = parallel["enabled"].as<bool>();
  }

  if (!parallel["partitioning"]) {
    throw std::runtime_error(
        "Parallel config '" + filename +
        "' is missing the required 'parallel.partitioning' block "
        "(must specify level, centroids_file, adjacency_file).");
  }
  parseParallelPartitioning(parallel["partitioning"], filename, config);

  if (parallel["chunked_loading"]) {
    YAML::Node chunked = parallel["chunked_loading"];
    rejectUnknownKeys(chunked,
                      {"person_metadata_chunk_size", "geo_unit_chunk_size"},
                      filename + ": ", " in 'parallel.chunked_loading'");
    if (chunked["person_metadata_chunk_size"]) {
      config.person_metadata_chunk_size =
          chunked["person_metadata_chunk_size"].as<size_t>();
    }
    if (chunked["geo_unit_chunk_size"]) {
      config.geo_unit_chunk_size = chunked["geo_unit_chunk_size"].as<size_t>();
    }
  }

  return config;
}

ActivityPreferenceConfig ConfigLoader::loadActivityPreferences(
    const std::string& filename) {
  ActivityPreferenceConfig config;

  // No activity_preferences_file in config_paths means every weight is 1.0. A
  // path that was written down must load.
  if (filename.empty()) return config;

  YAML::Node root = loadConfiguredFile(filename, "Activity preferences");
  try {
    rejectUnknownKeys(root, {"profiles"}, filename + ": ");
    if (root["profiles"]) {
      for (const auto& profile_node : root["profiles"]) {
        PreferenceProfile profile;
        profile.name = profile_node["name"]
                           ? profile_node["name"].as<std::string>()
                           : "unnamed";
        rejectUnknownKeys(
            profile_node,
            {"name", "activity", "priority", "selection", "weights"},
            filename + ": ", " in profile '" + profile.name + "'");
        profile.activity = profile_node["activity"]
                               ? profile_node["activity"].as<std::string>()
                               : "";
        profile.priority =
            profile_node["priority"] ? profile_node["priority"].as<int>() : 0;

        // Load selection criteria
        if (profile_node["selection"].IsDefined()) {
          parseSelectionCriteria(profile_node["selection"],
                                 profile.selection_criteria);
        }

        // Load preference weights
        if (profile_node["weights"]) {
          for (const auto& weight_kv : profile_node["weights"]) {
            profile.preference_weights[weight_kv.first.as<std::string>()] =
                weight_kv.second.as<double>();
          }
        }
        config.profiles.push_back(profile);
      }
    }

    // Sort by priority
    std::sort(config.profiles.begin(), config.profiles.end(),
              [](const PreferenceProfile& a, const PreferenceProfile& b) {
                return a.priority > b.priority;
              });

  } catch (const YAML::Exception& e) {
    throw std::runtime_error("Activity preferences file " + filename + ": " +
                             e.what());
  }

  return config;
}

// loadVaccination is defined in config_loader_vaccination.cpp.
// loadCoordinatedEncounters is defined in config_loader_encounters.cpp.

}  // namespace june
