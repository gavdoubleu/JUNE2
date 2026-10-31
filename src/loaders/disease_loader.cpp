#include "loaders/disease_loader.h"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <unordered_set>

#include "loaders/config_loader_detail.h"
#include "utils/filtered_csv.h"

namespace june {

using ::june::config_detail::rejectUnknownKeys;

namespace {

void logCurveRescale(const std::string& context_label, const char* type_name,
                     double max_inf, double factor) {
  std::cout << "[DEBUG] Curve (" << context_label << "): type=" << type_name
            << "  old_peak=" << max_inf * factor
            << "  max_infectiousness=" << max_inf
            << "  rescale_factor=" << factor
            << "  (set max_infectiousness: " << max_inf * factor
            << " in YAML to match previous behaviour)\n";
}

}  // namespace

namespace {

// =============================================================================
// Parse Distribution Parameters from YAML
// =============================================================================

// The parameters each distribution type takes. `required` ones are the ones
// sampleFromDistribution cannot do without; the rest fall back to its defaults.
struct DistributionSpec {
  const char* type;
  std::vector<std::string> required;
  std::vector<std::string> optional;
};

const std::vector<DistributionSpec>& distributionSpecs() {
  static const std::vector<DistributionSpec> specs = {
      {"constant", {"value"}, {}},
      {"normal", {}, {"loc", "scale"}},
      {"lognormal", {"s", "loc", "scale"}, {}},
      {"beta", {"a", "b", "loc", "scale"}, {}},
      {"exponweib", {"a", "c", "loc", "scale"}, {}},
      {"gamma", {"shape", "scale"}, {"loc"}},
  };
  return specs;
}

std::string joinNames(const std::vector<std::string>& names) {
  std::string joined;
  for (const auto& name : names) joined += (joined.empty() ? "" : ", ") + name;
  return joined;
}

DistributionParams parseDistribution(const YAML::Node& dist_node,
                                     const std::string& context) {
  DistributionParams dist;

  if (!dist_node["type"]) {
    throw std::runtime_error(context +
                             ": distribution must have a 'type' field");
  }

  dist.type = dist_node["type"].as<std::string>();

  const DistributionSpec* spec = nullptr;
  std::vector<std::string> types;
  for (const auto& candidate : distributionSpecs()) {
    types.push_back(candidate.type);
    if (dist.type == candidate.type) spec = &candidate;
  }
  if (!spec) {
    throw std::runtime_error(context + ": unknown distribution type '" +
                             dist.type + "'. Types are " + joinNames(types) +
                             ".");
  }

  std::vector<std::string> params = spec->required;
  params.insert(params.end(), spec->optional.begin(), spec->optional.end());

  for (auto it = dist_node.begin(); it != dist_node.end(); ++it) {
    const std::string key = it->first.as<std::string>();
    if (key == "type") continue;
    if (std::find(params.begin(), params.end(), key) == params.end()) {
      throw std::runtime_error(
          context + ": unknown parameter '" + key + "' for a " + dist.type +
          " distribution. Parameters are " + joinNames(params) + ".");
    }
    try {
      dist.params[key] = it->second.as<double>();
    } catch (const YAML::Exception&) {
      throw std::runtime_error(context + ": parameter '" + key + "' of a " +
                               dist.type +
                               " distribution must be a number, got '" +
                               YAML::Dump(it->second) + "'.");
    }
  }

  for (const auto& key : spec->required) {
    if (!dist.params.count(key)) {
      throw std::runtime_error(context + ": " + dist.type +
                               " distribution is missing '" + key +
                               "'. Parameters are " + joinNames(params) + ".");
    }
  }

  return dist;
}

// Throws unless `name` is one of the disease's symptom tags. `where` says
// which block named it.
void requireSymptomTag(const std::string& name, const std::string& where,
                       const std::vector<SymptomTag>& symptom_tags) {
  for (const auto& tag : symptom_tags) {
    if (tag.name == name) return;
  }
  std::string tags;
  for (const auto& tag : symptom_tags) {
    if (!tags.empty()) tags += ", ";
    tags += tag.name;
  }
  throw std::runtime_error("disease: " + where + " names stage '" + name +
                           "', which is not a symptom tag. Tags: " + tags);
}

// =============================================================================
// Parse Trajectory Definitions from YAML
// =============================================================================

TrajectoryStage parseTrajectoryStage(const YAML::Node& stage_node,
                                     const std::string& traj_label) {
  rejectUnknownKeys(stage_node, {"symptom_tag", "completion_time"}, "",
                    "a stage of " + traj_label);
  if (!stage_node["symptom_tag"]) {
    throw std::runtime_error("disease: a stage of " + traj_label +
                             " is missing 'symptom_tag'");
  }
  TrajectoryStage stage;
  stage.symptom_tag = stage_node["symptom_tag"].as<std::string>();
  if (stage_node["completion_time"]) {
    stage.completion_time = parseDistribution(
        stage_node["completion_time"],
        traj_label + " stage '" + stage.symptom_tag + "' completion_time");
  }
  return stage;
}

TrajectoryDefinition parseOneTrajectory(const YAML::Node& traj_node,
                                        size_t index) {
  rejectUnknownKeys(
      traj_node,
      {"description", "selection_key", "probability", "severity",
       "infectiousness_factor", "start_stage", "stages"},
      "disease: ", " in trajectory[" + std::to_string(index) + "]");
  TrajectoryDefinition traj;

  if (traj_node["description"]) {
    traj.description = traj_node["description"].as<std::string>();
  }

  if (traj_node["selection_key"]) {
    traj.selection_key = traj_node["selection_key"].as<std::string>();
  }

  if (traj_node["probability"]) {
    traj.probability = traj_node["probability"].as<double>();
  }

  traj.severity =
      traj_node["severity"] ? traj_node["severity"].as<double>() : 0.0;

  traj.infectiousness_factor =
      traj_node["infectiousness_factor"]
          ? traj_node["infectiousness_factor"].as<double>()
          : 1.0;

  if (traj_node["start_stage"]) {
    traj.start_stage = traj_node["start_stage"].as<std::string>();
  }

  const std::string traj_label =
      traj.description.empty() ? ("trajectory[" + std::to_string(index) + "]")
                               : ("trajectory \"" + traj.description + "\"");
  if (!traj_node["stages"]) {
    throw std::runtime_error("disease: " + traj_label + " is missing 'stages'");
  }
  for (const auto& stage_node : traj_node["stages"]) {
    traj.stages.push_back(parseTrajectoryStage(stage_node, traj_label));
  }

  return traj;
}

std::vector<TrajectoryDefinition> parseTrajectories(
    const YAML::Node& trajectories_node) {
  std::vector<TrajectoryDefinition> trajectories;

  if (!trajectories_node || !trajectories_node.IsSequence()) {
    std::cerr << "Warning: No trajectories defined or invalid format"
              << std::endl;
    return trajectories;
  }

  size_t index = 0;
  for (const auto& traj_node : trajectories_node) {
    trajectories.push_back(parseOneTrajectory(traj_node, index++));
  }

  return trajectories;
}

}  // namespace

// =============================================================================
// Load Outcome Rates from filter-column CSV
// =============================================================================

OutcomeRates DiseaseLoader::loadOutcomeRatesFromCSV(
    const std::string& csv_path) {
  OutcomeRates rates;
  csv::FilteredTable table = csv::loadFilteredCSV(csv_path);

  int row_num = 0;
  for (const auto& r : table.rows) {
    ++row_num;
    OutcomeRow row;
    row.criteria = r.criteria;
    for (const auto& name : table.value_columns) {
      auto it = r.values.find(name);
      if (it == r.values.end() || it->second.empty()) continue;
      try {
        row.probabilities[name] = std::stod(it->second);
      } catch (const std::exception&) {
        throw std::runtime_error("Outcome rates CSV '" + csv_path + "' row " +
                                 std::to_string(row_num) + " column '" + name +
                                 "' has non-numeric value '" + it->second +
                                 "'");
      }
    }
    rates.rows.push_back(std::move(row));
  }

  return rates;
}

// =============================================================================
// Section helpers for loadFromYAML
// =============================================================================

namespace {

std::shared_ptr<InfectiousnessCurve> parseCurve(
    const YAML::Node& curve_node, const std::string& context_label,
    bool verbose);
void parseDepositionStages(
    const YAML::Node& mode_node,
    std::vector<std::shared_ptr<InfectiousnessCurve>>& deposition_by_symptom,
    const std::string& mode_type_prefix, const std::string& mode_name,
    const std::vector<SymptomTag>& symptom_tags, bool verbose);
void loadTransmissionTrajectoryDriven(const YAML::Node& trans_node,
                                      TransmissionParams& transmission);
void loadTransmissionStageDriven(const YAML::Node& trans_node,
                                 TransmissionParams& transmission,
                                 const std::vector<SymptomTag>& symptom_tags,
                                 bool verbose);
void loadTransmissionStageDrivenFlat(
    const YAML::Node& trans_node, TransmissionParams& transmission,
    const std::vector<SymptomTag>& symptom_tags, bool verbose);
void parseStageDrivenModes(const YAML::Node& modes_node,
                           TransmissionParams& transmission,
                           const std::vector<SymptomTag>& symptom_tags,
                           bool verbose);
void attachStageCurvesToModes(const YAML::Node& stage_curves_node,
                              TransmissionParams& transmission,
                              const std::vector<SymptomTag>& symptom_tags,
                              bool verbose);
void finalizeStageDrivenMultiMode(TransmissionParams& transmission,
                                  const std::vector<SymptomTag>& symptom_tags);
void loadModeCompartmentalUptake(const YAML::Node& mode_node,
                                 TransmissionMode& tmode, int mode_idx);
void loadModeCompartmentalDeposition(
    const YAML::Node& mode_node, TransmissionMode& tmode, int mode_idx,
    const std::vector<SymptomTag>& symptom_tags, bool verbose);

void loadModeFomite(const YAML::Node& mode_node, TransmissionMode& tmode,
                    int mode_idx, const std::vector<SymptomTag>& symptom_tags,
                    bool verbose) {
  tmode.type = TransmissionModeType::Fomite;
  rejectUnknownKeys(
      mode_node,
      {"name", "mode_transmissibility_multiplier", "type", "max_age",
       "sub_bin_time", "fomite_curve", "deposition_stages"},
      "disease: ", " in fomite mode '" + tmode.name + "'");
  FomiteConfig fcfg;
  fcfg.mode_index = mode_idx;
  fcfg.max_age =
      mode_node["max_age"] ? mode_node["max_age"].as<double>() : 14.0;
  fcfg.sub_bin_time =
      mode_node["sub_bin_time"] ? mode_node["sub_bin_time"].as<double>() : 0.0;
  if (mode_node["fomite_curve"]) {
    fcfg.infectiousness_curve =
        parseCurve(mode_node["fomite_curve"],
                   "fomite / " + tmode.name + " / fomite_curve", verbose);
  } else {
    fcfg.infectiousness_curve = std::make_shared<ConstantCurve>(0.0);
  }
  parseDepositionStages(mode_node, fcfg.deposition_by_symptom, "fomite",
                        tmode.name, symptom_tags, verbose);
  tmode.config = std::move(fcfg);
  std::cout << "[DiseaseLoader] Fomite mode '" << tmode.name
            << "' registered at index " << mode_idx << std::endl;
}

void loadTransmission(const YAML::Node& config,
                      TransmissionParams& transmission,
                      const std::vector<SymptomTag>& symptom_tags,
                      bool verbose) {
  if (!config["transmission"]) return;
  auto trans_node = config["transmission"];

  // Detect mode (default to Trajectory-Driven)
  std::string mode_str = trans_node["mode"]
                             ? trans_node["mode"].as<std::string>()
                             : "Trajectory-Driven";
  if (mode_str != "Stage-Driven" && mode_str != "Trajectory-Driven") {
    throw std::runtime_error("disease: transmission.mode '" + mode_str +
                             "' is not a mode. Modes: Stage-Driven, "
                             "Trajectory-Driven");
  }
  transmission.mode = (mode_str == "Stage-Driven")
                          ? InfectiousnessMode::STAGE_DRIVEN
                          : InfectiousnessMode::TRAJECTORY_DRIVEN;
  if (transmission.mode == InfectiousnessMode::TRAJECTORY_DRIVEN) {
    rejectUnknownKeys(trans_node,
                      {"mode", "type", "max_infectiousness", "shape", "rate",
                       "shift", "modes"},
                      "disease: ", " in 'transmission' (Trajectory-Driven)");
  } else {
    rejectUnknownKeys(trans_node, {"mode", "modes", "stage_curves"},
                      "disease: ", " in 'transmission' (Stage-Driven)");
  }

  if (transmission.mode == InfectiousnessMode::TRAJECTORY_DRIVEN) {
    loadTransmissionTrajectoryDriven(trans_node, transmission);
  } else {
    loadTransmissionStageDriven(trans_node, transmission, symptom_tags,
                                verbose);
  }
  // Mode indices are narrowed to uint8_t like symptom ids; kNoModeIndex
  // reserves 255 as the "absent" sentinel, so a real index can't reach it.
  if (transmission.modes.size() >= kNoModeIndex) {
    throw std::runtime_error(
        "DiseaseLoader::loadTransmission: too many transmission modes (" +
        std::to_string(transmission.modes.size()) +
        "); indices are stored as uint8_t and must stay below " +
        std::to_string(static_cast<int>(kNoModeIndex)));
  }
}

void loadTransmissionStageDriven(const YAML::Node& trans_node,
                                 TransmissionParams& transmission,
                                 const std::vector<SymptomTag>& symptom_tags,
                                 bool verbose) {
  if (trans_node["modes"]) {
    // Multi-mode: parse modes list with per-mode stage_curves
    parseStageDrivenModes(trans_node["modes"], transmission, symptom_tags,
                          verbose);
    if (trans_node["stage_curves"]) {
      attachStageCurvesToModes(trans_node["stage_curves"], transmission,
                               symptom_tags, verbose);
    }
    finalizeStageDrivenMultiMode(transmission, symptom_tags);
  } else {
    loadTransmissionStageDrivenFlat(trans_node, transmission, symptom_tags,
                                    verbose);
  }
}

void loadTransmissionStageDrivenFlat(
    const YAML::Node& trans_node, TransmissionParams& transmission,
    const std::vector<SymptomTag>& symptom_tags, bool verbose) {
  // Flat stage_curves at top level: single Standard mode.
  std::map<std::string, std::shared_ptr<InfectiousnessCurve>> stage_curves;
  if (trans_node["stage_curves"]) {
    for (auto it = trans_node["stage_curves"].begin();
         it != trans_node["stage_curves"].end(); ++it) {
      std::string stage_name = it->first.as<std::string>();
      requireSymptomTag(stage_name, "stage_curves", symptom_tags);
      stage_curves[stage_name] =
          parseCurve(it->second, "stage_curves / " + stage_name, verbose);
    }
  }

  // Populate the ID-indexed vector used for runtime lookups.
  transmission.symptom_id_curves.resize(symptom_tags.size(), nullptr);
  for (const auto& tag : symptom_tags) {
    auto it = stage_curves.find(tag.name);
    if (it != stage_curves.end()) {
      transmission.symptom_id_curves[tag.id] = it->second;
    }
  }

  // Single Standard mode.
  TransmissionMode tmode;
  tmode.name = "default";
  tmode.symptom_curves = transmission.symptom_id_curves;
  transmission.modes.push_back(std::move(tmode));
}

void finalizeStageDrivenMultiMode(TransmissionParams& transmission,
                                  const std::vector<SymptomTag>& symptom_tags) {
  // Fill zero symptom_curves for modes without explicit stage_curves
  // (fomite, compartmental modes have no person-to-person infectiousness).
  std::vector<std::shared_ptr<InfectiousnessCurve>> zero_curves(
      symptom_tags.size(), std::make_shared<ConstantCurve>(0.0));
  for (auto& tmode : transmission.modes) {
    if (tmode.symptom_curves.empty()) {
      tmode.symptom_curves = zero_curves;
    }
  }

  // Populate flat symptom_id_curves from first Standard mode.
  transmission.symptom_id_curves.resize(symptom_tags.size(), nullptr);
  for (const auto& tmode : transmission.modes) {
    if (tmode.type == TransmissionModeType::Standard) {
      transmission.symptom_id_curves = tmode.symptom_curves;
      break;
    }
  }
}

void attachStageCurvesToModes(const YAML::Node& stage_curves_node,
                              TransmissionParams& transmission,
                              const std::vector<SymptomTag>& symptom_tags,
                              bool verbose) {
  for (const auto& mode_kv : stage_curves_node) {
    std::string mname = mode_kv.first.as<std::string>();
    const bool is_mode =
        std::any_of(transmission.modes.begin(), transmission.modes.end(),
                    [&](const TransmissionMode& m) { return m.name == mname; });
    if (!is_mode) {
      std::string modes;
      for (const auto& m : transmission.modes) {
        if (!modes.empty()) modes += ", ";
        modes += m.name;
      }
      throw std::runtime_error(
          "disease: transmission.stage_curves names mode '" + mname +
          "', which is not a transmission mode. Modes: " + modes);
    }
    std::vector<std::shared_ptr<InfectiousnessCurve>> mode_curves(
        symptom_tags.size(), nullptr);
    for (auto it = mode_kv.second.begin(); it != mode_kv.second.end(); ++it) {
      std::string stage_name = it->first.as<std::string>();
      requireSymptomTag(stage_name, "stage_curves." + mname, symptom_tags);
      auto curve = parseCurve(it->second, mname + " / " + stage_name, verbose);
      for (const auto& tag : symptom_tags) {
        if (tag.name == stage_name) {
          mode_curves[tag.id] = curve;
          break;
        }
      }
    }
    // Find mode index by name and store symptom_curves on the mode.
    for (auto& tmode : transmission.modes) {
      if (tmode.name == mname) {
        for (const auto& tag : symptom_tags) {
          if (tag.id < mode_curves.size() && !mode_curves[tag.id]) {
            if (tag.value > 0) {
              std::cerr << "[DiseaseLoader] Warning: mode '" << mname
                        << "' has no stage_curve for symptom '" << tag.name
                        << "'; using zero infectiousness." << std::endl;
            }
            mode_curves[tag.id] = std::make_shared<ConstantCurve>(0.0);
          }
        }
        tmode.symptom_curves = mode_curves;
        break;
      }
    }
  }
}

void parseStageDrivenModes(const YAML::Node& modes_node,
                           TransmissionParams& transmission,
                           const std::vector<SymptomTag>& symptom_tags,
                           bool verbose) {
  for (const auto& mode_node : modes_node) {
    TransmissionMode tmode;
    tmode.name =
        mode_node["name"] ? mode_node["name"].as<std::string>() : "default";
    tmode.mode_transmissibility_multiplier =
        mode_node["mode_transmissibility_multiplier"]
            ? mode_node["mode_transmissibility_multiplier"].as<double>()
            : 1.0;

    std::string mode_type =
        mode_node["type"] ? mode_node["type"].as<std::string>() : "";
    int mode_idx = static_cast<int>(transmission.modes.size());

    if (mode_type.empty()) {
      rejectUnknownKeys(mode_node, {"name", "mode_transmissibility_multiplier"},
                        "disease: ", " in mode '" + tmode.name + "'");
    }
    if (mode_type == "fomite") {
      loadModeFomite(mode_node, tmode, mode_idx, symptom_tags, verbose);
    } else if (mode_type == "compartmental_uptake") {
      loadModeCompartmentalUptake(mode_node, tmode, mode_idx);
    } else if (mode_type == "compartmental_deposition") {
      loadModeCompartmentalDeposition(mode_node, tmode, mode_idx, symptom_tags,
                                      verbose);
    } else if (!mode_type.empty()) {
      throw std::runtime_error(
          "disease: unknown type '" + mode_type + "' in mode '" + tmode.name +
          "'. Types are fomite, compartmental_uptake, compartmental_deposition "
          "(leave 'type' out for an ordinary mode).");
    }

    transmission.modes.push_back(std::move(tmode));
  }
}

void loadModeCompartmentalUptake(const YAML::Node& mode_node,
                                 TransmissionMode& tmode, int mode_idx) {
  rejectUnknownKeys(
      mode_node, {"name", "mode_transmissibility_multiplier", "type"},
      "disease: ", " in compartmental_uptake mode '" + tmode.name + "'");
  tmode.type = TransmissionModeType::CompartmentalUptake;
  CompartmentalUptakeConfig ucfg;
  ucfg.mode_index = mode_idx;
  tmode.config = ucfg;
  std::cout << "[DiseaseLoader] Compartmental uptake mode '" << tmode.name
            << "' registered at index " << mode_idx << std::endl;
}

void loadModeCompartmentalDeposition(
    const YAML::Node& mode_node, TransmissionMode& tmode, int mode_idx,
    const std::vector<SymptomTag>& symptom_tags, bool verbose) {
  tmode.type = TransmissionModeType::CompartmentalDeposition;
  rejectUnknownKeys(
      mode_node,
      {"name", "mode_transmissibility_multiplier", "type", "deposition_stages"},
      "disease: ", " in compartmental_deposition mode '" + tmode.name + "'");
  CompartmentalDepositionConfig dcfg;
  dcfg.mode_index = mode_idx;
  parseDepositionStages(mode_node, dcfg.deposition_by_symptom,
                        "compartmental_deposition", tmode.name, symptom_tags,
                        verbose);
  tmode.config = std::move(dcfg);
  std::cout << "[DiseaseLoader] Compartmental deposition mode '" << tmode.name
            << "' registered at index " << mode_idx << std::endl;
}

void loadTransmissionTrajectoryDriven(const YAML::Node& trans_node,
                                      TransmissionParams& transmission) {
  if (trans_node["type"]) {
    transmission.type = trans_node["type"].as<std::string>();
  } else {
    transmission.type = "gamma";
  }

  // Load infectiousness distribution parameters
  if (trans_node["max_infectiousness"]) {
    transmission.max_infectiousness = parseDistribution(
        trans_node["max_infectiousness"], "transmission max_infectiousness");
  }
  if (trans_node["shape"]) {
    transmission.shape =
        parseDistribution(trans_node["shape"], "transmission shape");
  }
  if (trans_node["rate"]) {
    transmission.rate =
        parseDistribution(trans_node["rate"], "transmission rate");
  }
  if (trans_node["shift"]) {
    transmission.shift =
        parseDistribution(trans_node["shift"], "transmission shift");
  }

  // Note: symptom-specific infectiousness scaling is now handled via
  // infectiousness_factor on each TrajectoryDefinition, applied once at
  // infection creation (matching original JUNE behaviour).

  // Multi-mode support for TRAJECTORY_DRIVEN: parse optional modes list
  // for mode names and susceptibility multipliers. No per-mode curves
  // are loaded; getInfectiousness(m, t) returns the same scalar for all
  // modes; the multipliers scale susceptibility per-mode at the FOI step.
  if (trans_node["modes"]) {
    for (const auto& mode_node : trans_node["modes"]) {
      rejectUnknownKeys(mode_node, {"name", "mode_transmissibility_multiplier"},
                        "disease: ", " in a mode of 'transmission'");
      TransmissionMode tmode;
      tmode.name =
          mode_node["name"] ? mode_node["name"].as<std::string>() : "default";
      tmode.mode_transmissibility_multiplier =
          mode_node["mode_transmissibility_multiplier"]
              ? mode_node["mode_transmissibility_multiplier"].as<double>()
              : 1.0;
      transmission.modes.push_back(std::move(tmode));
    }
  } else {
    TransmissionMode tmode;
    tmode.name = "default";
    transmission.modes.push_back(std::move(tmode));
  }
}

void parseDepositionStages(
    const YAML::Node& mode_node,
    std::vector<std::shared_ptr<InfectiousnessCurve>>& deposition_by_symptom,
    const std::string& mode_type_prefix, const std::string& mode_name,
    const std::vector<SymptomTag>& symptom_tags, bool verbose) {
  deposition_by_symptom.resize(symptom_tags.size(), nullptr);
  if (!mode_node["deposition_stages"]) return;
  for (auto it = mode_node["deposition_stages"].begin();
       it != mode_node["deposition_stages"].end(); ++it) {
    std::string stage_name = it->first.as<std::string>();
    requireSymptomTag(
        stage_name,
        mode_type_prefix + " mode '" + mode_name + "' deposition_stages",
        symptom_tags);
    auto curve = parseCurve(
        it->second, mode_type_prefix + " / " + mode_name + " / " + stage_name,
        verbose);
    for (const auto& tag : symptom_tags) {
      if (tag.name == stage_name) {
        deposition_by_symptom[tag.id] = curve;
        break;
      }
    }
  }
}

void loadNaturalImmunity(const YAML::Node& config,
                         TransmissionParams& transmission) {
  if (!config["immunity"]) return;
  auto immunity_node = config["immunity"];
  rejectUnknownKeys(immunity_node, {"level", "waning_rate"},
                    "disease: ", " in 'immunity'");
  transmission.natural_immunity.level =
      immunity_node["level"] ? immunity_node["level"].as<double>() : 0.95;
  transmission.natural_immunity.waning_rate =
      immunity_node["waning_rate"] ? immunity_node["waning_rate"].as<double>()
                                   : 0.001;
}

void validateOutcomeRowSums(const OutcomeRates& outcome_rates,
                            const YAML::Node& config) {
  std::vector<std::string> populations;
  std::vector<std::string> sexes;
  std::vector<std::string> outcomes;

  const YAML::Node& csv_node = config["outcome_rates_csv"];
  if (csv_node && csv_node.IsMap()) {
    if (csv_node["populations"] && csv_node["populations"].IsMap()) {
      for (const auto& population : csv_node["populations"]) {
        populations.push_back(population.first.as<std::string>());
      }
    }
    if (csv_node["sexes"] && csv_node["sexes"].IsSequence()) {
      for (const auto& sex : csv_node["sexes"]) {
        sexes.push_back(sex.as<std::string>());
      }
    }
    if (csv_node["outcomes"] && csv_node["outcomes"].IsMap()) {
      for (const auto& outcome : csv_node["outcomes"]) {
        outcomes.push_back(outcome.first.as<std::string>());
      }
    }
  }

  // The declared population/outcome/sex names describe the wide column
  // convention: <population>_<outcome>_<sex>. Only select this path when at
  // least one declared wide column is present, so configs whose tables use
  // filter.* columns retain the legacy per-row validation.
  bool wide_layout = false;
  if (!populations.empty() && !sexes.empty() && !outcomes.empty()) {
    for (const auto& population : populations) {
      for (const auto& sex : sexes) {
        for (const auto& outcome : outcomes) {
          const std::string column = population + "_" + outcome + "_" + sex;
          for (const auto& row : outcome_rates.rows) {
            if (row.probabilities.count(column) > 0) {
              wide_layout = true;
              break;
            }
          }
          if (wide_layout) break;
        }
        if (wide_layout) break;
      }
      if (wide_layout) break;
    }
  }

  for (size_t row_i = 0; row_i < outcome_rates.rows.size(); ++row_i) {
    const auto& probabilities = outcome_rates.rows[row_i].probabilities;
    if (wide_layout) {
      for (const auto& population : populations) {
        for (const auto& sex : sexes) {
          double group_sum = 0.0;
          bool has_value = false;
          for (const auto& outcome : outcomes) {
            const std::string column = population + "_" + outcome + "_" + sex;
            auto it = probabilities.find(column);
            if (it != probabilities.end()) {
              group_sum += it->second;
              has_value = true;
            }
          }
          if (!has_value || std::abs(group_sum - 1.0) > 0.01) {
            throw std::runtime_error(
                "Outcome rates row " + std::to_string(row_i) + " group '" +
                population + "/" + sex + "' sums to " +
                std::to_string(group_sum) +
                ", expected 1.0. Outcome columns are a probability "
                "distribution over trajectories and must sum to 1.");
          }
        }
      }
    } else {
      double row_sum = 0.0;
      for (const auto& [key, prob] : probabilities) row_sum += prob;
      if (!probabilities.empty() && std::abs(row_sum - 1.0) > 0.01) {
        // Reject rows whose probability columns do not form a distribution.
        // The trajectory sampler cannot use incomplete or oversized
        // probability mass.
        throw std::runtime_error(
            "Outcome rates row " + std::to_string(row_i) + " sums to " +
            std::to_string(row_sum) +
            ", expected 1.0. Outcome columns are a "
            "probability distribution over trajectories and must sum to 1.");
      }
    }
    for (const auto& [key, prob] : probabilities) {
      // Second net, not the first one: a negative rate can hide inside a row
      // that still sums to 1.
      if (prob < 0.0) {
        throw std::runtime_error("Outcome rates row " + std::to_string(row_i) +
                                 " has negative '" + key +
                                 "' = " + std::to_string(prob) + ".");
      }
    }
  }
}

OutcomeRates loadOutcomeRatesFromConfig(
    const YAML::Node& config, const std::string& yaml_path,
    std::vector<std::string>* referenced_paths) {
  OutcomeRates outcome_rates;
  if (!config["outcome_rates_csv"]) return outcome_rates;
  const auto& csv_node = config["outcome_rates_csv"];

  // Accept either a scalar path ("outcome_rates_csv: path/to/file.csv")
  // or a mapping with a "file:" key (legacy format).
  std::string csv_rel_path;
  if (csv_node.IsScalar()) {
    csv_rel_path = csv_node.as<std::string>();
  } else if (csv_node["file"]) {
    rejectUnknownKeys(csv_node, {"file", "populations", "sexes", "outcomes"},
                      "disease: ", " in 'outcome_rates_csv'");
    csv_rel_path = csv_node["file"].as<std::string>();
  } else {
    throw std::runtime_error(
        "outcome_rates_csv must be a file path or a mapping with a 'file:' "
        "key");
  }

  std::string yaml_dir = yaml_path.substr(0, yaml_path.find_last_of("/\\") + 1);
  const std::string csv_path = yaml_dir + csv_rel_path;
  if (referenced_paths) referenced_paths->push_back(csv_path);
  return DiseaseLoader::loadOutcomeRatesFromCSV(csv_path);
}

void validateTrajectoryStageRefs(
    const std::vector<TrajectoryDefinition>& trajectories,
    const std::vector<SymptomTag>& symptom_tags) {
  std::unordered_set<std::string> valid_tags;
  for (const auto& tag : symptom_tags) {
    valid_tags.insert(tag.name);
  }
  for (size_t ti = 0; ti < trajectories.size(); ++ti) {
    const auto& traj = trajectories[ti];
    const std::string traj_label =
        traj.description.empty() ? ("trajectory[" + std::to_string(ti) + "]")
                                 : ("trajectory \"" + traj.description + "\"");
    // A configured start_stage must name a defined symptom tag.
    if (traj.start_stage.has_value() &&
        valid_tags.find(*traj.start_stage) == valid_tags.end()) {
      throw std::runtime_error("Disease config error: " + traj_label +
                               " has start_stage \"" + *traj.start_stage +
                               "\" which is not a defined symptom tag.");
    }
    // Every trajectory stage must name a defined symptom tag.
    for (size_t si = 0; si < traj.stages.size(); ++si) {
      const auto& stage = traj.stages[si];
      if (valid_tags.find(stage.symptom_tag) == valid_tags.end()) {
        throw std::runtime_error("Disease config error: " + traj_label +
                                 " stage[" + std::to_string(si) +
                                 "] has symptom_tag \"" + stage.symptom_tag +
                                 "\" which is not a defined symptom tag.");
      }
    }
  }
}

DiseaseStageSettings loadStageSettings(const YAML::Node& config) {
  DiseaseStageSettings stage_settings;
  if (!config["settings"]) return stage_settings;
  auto settings = config["settings"];
  rejectUnknownKeys(settings,
                    {"hospitalised_stage", "intensive_care_stage",
                     "fatality_stage", "recovered_stage"},
                    "disease: ", " in 'settings'");

  auto loadStageList = [](const YAML::Node& node,
                          std::vector<std::string>& out) {
    if (!node) return;
    for (const auto& item : node) {
      rejectUnknownKeys(item, {"name"}, "disease: ", " in 'settings'");
      if (item["name"]) {
        out.push_back(item["name"].as<std::string>());
      }
    }
  };

  loadStageList(settings["hospitalised_stage"],
                stage_settings.hospitalised_stages);
  loadStageList(settings["intensive_care_stage"],
                stage_settings.intensive_care_stages);
  loadStageList(settings["fatality_stage"], stage_settings.fatality_stages);
  loadStageList(settings["recovered_stage"], stage_settings.recovered_stages);
  return stage_settings;
}

std::vector<SymptomTag> loadSymptomTags(const YAML::Node& config) {
  std::vector<SymptomTag> symptom_tags;
  if (!config["symptom_tags"]) return symptom_tags;
  // Symptom ids are narrowed to uint8_t for event logging and cross-rank
  // transmission (PendingInfection, InfectionEvent); kNoSymptomId reserves
  // 255 as the "absent" sentinel, so a real id can't reach it.
  if (config["symptom_tags"].size() >= kNoSymptomId) {
    throw std::runtime_error(
        "DiseaseLoader::loadSymptomTags: too many symptom_tags (" +
        std::to_string(config["symptom_tags"].size()) +
        "); ids are stored as uint8_t and must stay below " +
        std::to_string(static_cast<int>(kNoSymptomId)));
  }
  uint16_t current_id = 0;
  for (const auto& tag_node : config["symptom_tags"]) {
    rejectUnknownKeys(tag_node, {"name", "value"},
                      "disease: ", " in a symptom tag");
    SymptomTag tag;
    tag.name = tag_node["name"].as<std::string>();
    tag.value = tag_node["value"].as<int>();
    tag.id = current_id++;
    symptom_tags.push_back(tag);
  }
  return symptom_tags;
}

}  // namespace

// =============================================================================
// Main Loading Function
// =============================================================================

Disease DiseaseLoader::loadFromYAML(
    const std::string& yaml_path, bool verbose,
    std::vector<std::string>* referenced_paths) {
  try {
    YAML::Node root = YAML::LoadFile(yaml_path);

    // Get disease node
    if (!root["disease"]) {
      throw std::runtime_error("YAML must contain 'disease' section");
    }

    YAML::Node config = root["disease"];
    rejectUnknownKeys(root, {"disease"}, "disease: ");
    rejectUnknownKeys(config,
                      {"name", "symptom_tags", "settings", "trajectories",
                       "outcome_rates_csv", "transmission", "immunity"},
                      "disease: ", " in 'disease'");

    // Get disease name
    std::string disease_name = config["name"].as<std::string>("covid19");

    std::vector<SymptomTag> symptom_tags = loadSymptomTags(config);

    DiseaseStageSettings stage_settings = loadStageSettings(config);

    // === Parse Trajectories ===
    std::vector<TrajectoryDefinition> trajectories;
    if (config["trajectories"]) {
      trajectories = parseTrajectories(config["trajectories"]);
    }

    validateTrajectoryStageRefs(trajectories, symptom_tags);

    OutcomeRates outcome_rates =
        loadOutcomeRatesFromConfig(config, yaml_path, referenced_paths);

    TransmissionParams transmission;
    loadTransmission(config, transmission, symptom_tags, verbose);

    loadNaturalImmunity(config, transmission);
    validateOutcomeRowSums(outcome_rates, config);

    // Create and return Disease object
    return Disease(disease_name, symptom_tags, stage_settings, trajectories,
                   outcome_rates, transmission);

  } catch (const YAML::Exception& e) {
    throw std::runtime_error("Failed to load disease config from " + yaml_path +
                             ": " + e.what());
  } catch (const std::exception& e) {
    throw std::runtime_error("Error loading disease config " + yaml_path +
                             ": " + std::string(e.what()));
  }
}

// =============================================================================
// Parse an infectiousness curve from YAML
// =============================================================================

namespace {

std::shared_ptr<InfectiousnessCurve> parseCurve(
    const YAML::Node& curve_node, const std::string& context_label,
    bool verbose) {
  if (!curve_node || !curve_node["type"]) {
    return std::make_shared<ConstantCurve>(0.0);
  }

  std::string type = curve_node["type"].as<std::string>();

  if (type == "constant") {
    rejectUnknownKeys(curve_node, {"type", "value"}, "disease: ",
                      " in curve '" + context_label + "' (constant)");
  } else if (type == "gamma") {
    rejectUnknownKeys(curve_node,
                      {"type", "max_infectiousness", "shape", "rate", "shift"},
                      "disease: ", " in curve '" + context_label + "' (gamma)");
  } else if (type == "exponential_decay") {
    rejectUnknownKeys(
        curve_node, {"type", "initial_value", "decay_rate", "delay"},
        "disease: ", " in curve '" + context_label + "' (exponential_decay)");
  } else if (type == "linear_ramp") {
    rejectUnknownKeys(
        curve_node, {"type", "start_value", "end_value", "ramp_duration"},
        "disease: ", " in curve '" + context_label + "' (linear_ramp)");
  } else if (type == "lognormal") {
    rejectUnknownKeys(
        curve_node, {"type", "max_infectiousness", "mu", "sigma"},
        "disease: ", " in curve '" + context_label + "' (lognormal)");
  } else if (type == "beta") {
    rejectUnknownKeys(
        curve_node, {"type", "max_infectiousness", "alpha", "beta", "duration"},
        "disease: ", " in curve '" + context_label + "' (beta)");
  }

  static constexpr double kTableMaxDays = 90.0;
  static constexpr int kTableNPoints = 2700;

  const auto valueOr = [&curve_node](const char* key, double default_value) {
    return curve_node[key] ? curve_node[key].as<double>() : default_value;
  };

  std::shared_ptr<InfectiousnessCurve> curve;
  double max_inf = 1.0;
  double peak_scaling_factor = 0.0;
  const char* rescaled_type = nullptr;
  if (type == "constant") {
    curve = std::make_shared<ConstantCurve>(valueOr("value", 1.0));
  } else if (type == "gamma") {
    max_inf = valueOr("max_infectiousness", 1.0);
    auto gamma = std::make_shared<GammaCurve>(max_inf, valueOr("shape", 1.56),
                                              valueOr("rate", 0.53),
                                              valueOr("shift", 0.0));
    peak_scaling_factor = gamma->peakScalingFactor();
    curve = std::move(gamma);
    rescaled_type = "gamma";
  } else if (type == "exponential_decay") {
    curve = std::make_shared<ExponentialDecayCurve>(
        valueOr("initial_value", 1.0), valueOr("decay_rate", 0.5),
        valueOr("delay", 0.0));
  } else if (type == "linear_ramp") {
    curve = std::make_shared<LinearRampCurve>(valueOr("start_value", 0.0),
                                              valueOr("end_value", 1.0),
                                              valueOr("ramp_duration", 1.0));
  } else if (type == "lognormal") {
    max_inf = valueOr("max_infectiousness", 1.0);
    auto lognormal = std::make_shared<LognormalCurve>(
        max_inf, valueOr("mu", 0.0), valueOr("sigma", 1.0));
    peak_scaling_factor = lognormal->peakScalingFactor();
    curve = std::move(lognormal);
    rescaled_type = "lognormal";
  } else if (type == "beta") {
    max_inf = valueOr("max_infectiousness", 1.0);
    auto beta = std::make_shared<BetaCurve>(max_inf, valueOr("alpha", 2.0),
                                            valueOr("beta", 2.0),
                                            valueOr("duration", 1.0));
    peak_scaling_factor = beta->peakScalingFactor();
    curve = std::move(beta);
    rescaled_type = "beta";
  } else {
    throw std::runtime_error("Unknown infectiousness curve type: " + type);
  }

  if (verbose && rescaled_type) {
    logCurveRescale(context_label, rescaled_type, max_inf, peak_scaling_factor);
  }

  curve->buildIntegralTable(kTableMaxDays, kTableNPoints);
  return curve;
}

}  // namespace

}  // namespace june
