#include <yaml-cpp/yaml.h>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "loaders/config_loader.h"
#include "loaders/config_loader_detail.h"
#include "utils/mpi_utils.h"

namespace june {

namespace {

using ::june::config_detail::parseSelectionCriteria;
using ::june::config_detail::rejectUnknownKeys;

// Campaign and tier ids are uint16_t in the campaign_days dataset
constexpr size_t kMaxEventIds = std::numeric_limits<uint16_t>::max();

// An efficacy or a probability, so it has to be from 0 to 1.
double readUnitInterval(const YAML::Node& node, const std::string& where) {
  double value = std::nan("");
  try {
    value = node.as<double>();
  } catch (const YAML::Exception&) {
  }
  // Written as a negation so a NaN fails the check too
  if (!(value >= 0.0 && value <= 1.0))
    throw std::runtime_error(where + " must be a number from 0 to 1, got '" +
                             YAML::Dump(node) + "'");
  return value;
}

// Parse a vaccine efficacy node: either a flat `disease -> scalar` map
// (taken as 0-100 age-uniform efficacy) or a nested `disease -> {age-range
// -> value}` map (each age-range encoded as "min-max"). Used for both
// infection_efficacy and symptom_efficacy on a DoseConfig. `where` names the
// dose and the node for error messages.
std::unordered_map<std::string, std::vector<AgeEfficacy>> parseDoseEfficacies(
    const YAML::Node& node, const std::string& where) {
  std::unordered_map<std::string, std::vector<AgeEfficacy>> res;
  if (!node) return res;
  for (const auto& eff_kv : node) {
    const std::string disease = eff_kv.first.as<std::string>();
    const std::string here = where + " " + disease;
    if (eff_kv.second.IsMap()) {
      for (const auto& age_kv : eff_kv.second) {
        const std::string band = age_kv.first.as<std::string>();
        const std::string band_where = here + " band '" + band + "'";
        // A band that is not min-max would otherwise be dropped without a word
        const size_t dash = band.find('-');
        AgeEfficacy ae{};
        size_t min_used = 0, max_used = 0;
        try {
          if (dash != std::string::npos) {
            ae.min_age = std::stoi(band.substr(0, dash), &min_used);
            ae.max_age = std::stoi(band.substr(dash + 1), &max_used);
          }
        } catch (const std::exception&) {
          min_used = max_used = 0;
        }
        if (dash == std::string::npos || min_used != dash ||
            max_used != band.size() - dash - 1 || ae.min_age < 0 ||
            ae.min_age > ae.max_age)
          throw std::runtime_error(
              band_where +
              " must be written min-max in whole years, such as "
              "0-17");
        ae.efficacy = readUnitInterval(age_kv.second, band_where + " efficacy");
        res[disease].push_back(ae);
      }
    } else {
      res[disease].push_back(
          {0, 100, readUnitInterval(eff_kv.second, here + " efficacy")});
    }
  }
  return res;
}

// Read a timing field every dose has to state. A dose with no schedule has no
// meaning, so a missing one throws instead of taking a value.
double requireDoseTiming(const YAML::Node& d_node, const std::string& key,
                         const std::string& vaccine, int number) {
  if (!d_node[key])
    throw std::runtime_error("vaccine '" + vaccine + "' dose " +
                             std::to_string(number) +
                             ": missing required field '" + key + "'");
  return d_node[key].as<double>();
}

// Parse one entry from a vaccine's `doses:` sequence. The position in the
// sequence becomes DoseConfig::number.
DoseConfig parseDose(const std::string& vaccine, int number,
                     const YAML::Node& d_node) {
  DoseConfig dc;
  dc.number = number;
  const std::string where =
      "vaccine '" + vaccine + "' dose " + std::to_string(number);
  // A misspelt key would otherwise leave its field at the default, so
  // days_to_onsett would quietly turn a step into a ramp
  rejectUnknownKeys(d_node,
                    {"days_to_onset", "days_to_effective", "days_to_waning",
                     "days_to_finished", "waning_factor", "infection_efficacy",
                     "symptom_efficacy"},
                    where + ": ");
  dc.days_to_effective =
      requireDoseTiming(d_node, "days_to_effective", vaccine, number);
  dc.days_to_waning =
      requireDoseTiming(d_node, "days_to_waning", vaccine, number);
  dc.days_to_finished =
      requireDoseTiming(d_node, "days_to_finished", vaccine, number);
  if (d_node["days_to_onset"])
    dc.days_to_onset = d_node["days_to_onset"].as<double>();
  // Written as a negation so a NaN fails the check too
  if (!(dc.days_to_onset >= 0.0 && dc.days_to_onset <= dc.days_to_effective)) {
    std::ostringstream msg;
    msg << "vaccine '" << vaccine << "' dose " << number
        << ": days_to_onset must be between 0 and days_to_effective, got "
        << dc.days_to_onset << " with days_to_effective "
        << dc.days_to_effective;
    throw std::runtime_error(msg.str());
  }
  if (d_node["waning_factor"])
    dc.waning_factor = d_node["waning_factor"].as<double>();
  dc.infection_efficacy = parseDoseEfficacies(d_node["infection_efficacy"],
                                              where + " infection_efficacy");
  dc.symptom_efficacy = parseDoseEfficacies(d_node["symptom_efficacy"],
                                            where + " symptom_efficacy");
  return dc;
}

// Parse one entry from the `vaccines:` map (a named vaccine and its dose
// schedule).
VaccineConfig parseVaccine(const std::string& name,
                           const YAML::Node& vacc_node) {
  VaccineConfig vc;
  vc.name = name;
  rejectUnknownKeys(vacc_node, {"doses"}, "vaccine '" + name + "': ");
  if (vacc_node["doses"]) {
    const auto& doses_node = vacc_node["doses"];
    for (size_t i = 0; i < doses_node.size(); ++i) {
      vc.doses.push_back(parseDose(name, static_cast<int>(i), doses_node[i]));
    }
  }
  return vc;
}

// A YAML number as a whole number, or an error naming where it was read.
// yaml-cpp refuses "1.5" and text for an integer type, which is the check
// wanted.
template <typename T>
T readWholeNumber(const YAML::Node& node, const std::string& where) {
  try {
    return node.as<T>();
  } catch (const YAML::Exception&) {
    throw std::runtime_error(where + " must be a whole number, got '" +
                             YAML::Dump(node) + "'");
  }
}

// The text of a key that must be one of a few words
template <typename T>
T readChoice(const YAML::Node& node, const std::string& where,
             std::initializer_list<std::pair<const char*, T>> choices) {
  std::string known;
  for (const auto& [word, value] : choices) {
    if (node.IsScalar() && node.Scalar() == word) return value;
    known += (known.empty() ? "" : ", ") + std::string(word);
  }
  throw std::runtime_error(where + " must be one of " + known + ", got '" +
                           YAML::Dump(node) + "'");
}

// The uptake of a tier: a number for everyone, a list of rows, or a CSV file.
// A relative file is taken as given, else relative to the YAML file's folder.
UptakeTable parseUptake(const YAML::Node& node, const std::string& where,
                        const std::filesystem::path& yaml_dir) {
  UptakeTable table;
  if (node.IsScalar()) {
    table.rows.push_back({{}, readUnitInterval(node, where)});
  } else if (node.IsSequence()) {
    if (node.size() == 0)
      throw std::runtime_error(where +
                               " is an empty list; leave 'uptake' out "
                               "for everyone to accept");
    int n = 0;
    for (const auto& row_node : node) {
      const std::string row_where = where + " row " + std::to_string(++n);
      if (!row_node.IsMap())
        throw std::runtime_error(row_where +
                                 " must be a mapping with 'probability' and "
                                 "optionally 'selection'");
      rejectUnknownKeys(row_node, {"selection", "probability"},
                        row_where + ": ");
      if (!row_node["probability"])
        throw std::runtime_error(row_where + " needs 'probability'");
      UptakeRow row;
      row.probability =
          readUnitInterval(row_node["probability"], row_where + " probability");
      parseSelectionCriteria(row_node["selection"], row.criteria);
      table.rows.push_back(std::move(row));
    }
  } else if (node.IsMap()) {
    rejectUnknownKeys(node, {"file"}, where + ": ");
    if (!node["file"] || !node["file"].IsScalar() ||
        node["file"].Scalar().empty())
      throw std::runtime_error(where + " needs 'file' with a path");
    std::filesystem::path path(node["file"].Scalar());
    if (path.is_relative() && !std::filesystem::exists(path))
      path = yaml_dir / path;
    table.csv_path = path.string();
  } else {
    throw std::runtime_error(where +
                             " must be a number, a list of rows or "
                             "{file: path}, got '" +
                             YAML::Dump(node) + "'");
  }
  return table;
}

// Either the word `person` or a map naming a shared venue, who is inside it and
// who decides for it.
void parseDecisionUnit(const YAML::Node& node, const std::string& where,
                       VaccinationTier& tier) {
  if (node.IsScalar() && node.Scalar() == "person") return;
  if (!node.IsMap())
    throw std::runtime_error(where + " must be 'person' or a mapping with " +
                             "shared_venue, scope and deciders, got '" +
                             YAML::Dump(node) + "'");
  rejectUnknownKeys(
      node, {"shared_venue", "scope", "deciders", "combine", "if_no_decider"},
      where + ": ");
  for (const char* key : {"shared_venue", "scope", "deciders"}) {
    if (!node[key]) throw std::runtime_error(where + " needs '" + key + "'");
  }
  if (!node["shared_venue"].IsScalar() || node["shared_venue"].Scalar().empty())
    throw std::runtime_error(where + " shared_venue needs an activity name");
  tier.decision_unit.shared_activity = node["shared_venue"].Scalar();
  DecisionRules& rules = tier.decision_rules;
  parseSelectionCriteria(node["scope"], rules.scope);
  parseSelectionCriteria(node["deciders"], rules.deciders);
  if (node["combine"])
    rules.combine =
        readChoice<DeciderCombine>(node["combine"], where + " combine",
                                   {{"min", DeciderCombine::Min},
                                    {"mean", DeciderCombine::Mean},
                                    {"max", DeciderCombine::Max}});
  if (node["if_no_decider"])
    rules.if_no_decider = readChoice<NoDecider>(
        node["if_no_decider"], where + " if_no_decider",
        {{"error", NoDecider::Error}, {"individual", NoDecider::Individual}});
}

VaccinationTier parseTier(const YAML::Node& node, const std::string& campaign,
                          size_t position,
                          const std::filesystem::path& yaml_dir) {
  std::string where = campaign + " tier " + std::to_string(position);
  if (!node.IsMap())
    throw std::runtime_error(where + " must be a mapping, got '" +
                             YAML::Dump(node) + "'");
  // Name the tier by its name in every message below once it has one
  const bool named =
      node["name"] && node["name"].IsScalar() && !node["name"].Scalar().empty();
  if (named) where = campaign + " tier '" + node["name"].Scalar() + "'";
  rejectUnknownKeys(node,
                    {"name", "dose", "selection", "uptake", "decision_unit"},
                    where + ": ");
  if (!named) throw std::runtime_error(where + " needs a 'name'");
  VaccinationTier tier;
  tier.name = node["name"].Scalar();
  if (!node["dose"]) throw std::runtime_error(where + " needs a 'dose'");
  tier.dose = readWholeNumber<int>(node["dose"], where + " dose");
  if (tier.dose < 0)
    throw std::runtime_error(where + " dose must be 0 or more, got " +
                             std::to_string(tier.dose));
  parseSelectionCriteria(node["selection"], tier.selection_criteria);
  if (node["uptake"])
    tier.uptake = parseUptake(node["uptake"], where + " uptake", yaml_dir);
  if (node["decision_unit"])
    parseDecisionUnit(node["decision_unit"], where + " decision_unit", tier);
  return tier;
}

DoseBudget parseDoseBudget(const YAML::Node& node, const std::string& where) {
  if (!node.IsMap())
    throw std::runtime_error(where + " must be a mapping with 'doses' and " +
                             "'per_days', got '" + YAML::Dump(node) + "'");
  rejectUnknownKeys(node, {"doses", "per_days"}, where + ": ");
  const auto field = [&](const char* key) {
    if (!node[key]) throw std::runtime_error(where + " needs '" + key + "'");
    const auto value = readWholeNumber<int64_t>(node[key], where + " " + key);
    if (value < 1)
      throw std::runtime_error(where + " " + key +
                               " must be greater than zero, got " +
                               std::to_string(value));
    return value;
  };
  return {field("doses"), field("per_days")};
}

// Parse one entry from the `campaigns:` map (a named vaccination campaign).
VaccinationCampaignConfig parseVaccinationCampaign(
    const std::string& name, const YAML::Node& camp_node,
    const std::filesystem::path& yaml_dir) {
  VaccinationCampaignConfig camp;
  camp.name = name;
  const std::string where = "vaccination campaign '" + name + "'";
  rejectUnknownKeys(
      camp_node,
      {"start_date", "end_date", "vaccine_type", "days_to_next_dose",
       "courses_begun_with", "dose_budget", "coverage_target", "tiers"},
      where + ": ");
  if (camp_node["start_date"])
    camp.start_date = camp_node["start_date"].as<std::string>();
  if (camp_node["end_date"])
    camp.end_date = camp_node["end_date"].as<std::string>();
  if (!camp_node["vaccine_type"])
    throw std::runtime_error(where + " needs 'vaccine_type'");
  camp.vaccine_type = camp_node["vaccine_type"].as<std::string>();
  if (camp_node["days_to_next_dose"])
    camp.days_to_next_dose =
        camp_node["days_to_next_dose"].as<std::vector<double>>();
  if (camp_node["courses_begun_with"]) {
    camp.courses_begun_with =
        camp_node["courses_begun_with"].as<std::vector<std::string>>();
  }

  if (camp_node["dose_budget"])
    camp.dose_budget =
        parseDoseBudget(camp_node["dose_budget"], where + " dose_budget");
  if (camp_node["coverage_target"]) {
    double target = std::nan("");
    try {
      target = camp_node["coverage_target"].as<double>();
    } catch (const YAML::Exception&) {
    }
    // Written as a negation so a NaN fails the check too
    if (!(target > 0.0 && target <= 1.0))
      throw std::runtime_error(
          where +
          " coverage_target must be a number above 0 and at most 1, "
          "got '" +
          YAML::Dump(camp_node["coverage_target"]) + "'");
    camp.coverage_target = target;
  }
  if (!camp_node["tiers"]) throw std::runtime_error(where + " needs 'tiers'");
  const YAML::Node& tiers = camp_node["tiers"];
  if (!tiers.IsSequence() || tiers.size() == 0)
    throw std::runtime_error(where + " tiers must be a non-empty list");
  if (tiers.size() > kMaxEventIds)
    throw std::runtime_error(where + " has " + std::to_string(tiers.size()) +
                             " tiers; at most " + std::to_string(kMaxEventIds) +
                             " are allowed");
  std::unordered_set<std::string> names;
  for (size_t i = 0; i < tiers.size(); ++i) {
    camp.tiers.push_back(parseTier(tiers[i], where, i + 1, yaml_dir));
    if (!names.insert(camp.tiers.back().name).second)
      throw std::runtime_error(where + " has two tiers named '" +
                               camp.tiers.back().name + "'");
  }
  return camp;
}

// A tier gives a dose of its campaign's vaccine, so both must exist. Checked
// once every vaccine is loaded, since the file may list them after the
// campaigns.
void checkTierDoses(const VaccinationConfig& config) {
  for (const auto& camp : config.campaigns) {
    const std::string where = "vaccination campaign '" + camp.name + "'";
    const auto vaccine = config.vaccines.find(camp.vaccine_type);
    if (vaccine == config.vaccines.end())
      throw std::runtime_error(where + " vaccine_type '" + camp.vaccine_type +
                               "' is not one of the vaccines in this file");
    const auto doses = static_cast<int>(vaccine->second.doses.size());
    for (const auto& tier : camp.tiers) {
      if (tier.dose >= doses)
        throw std::runtime_error(where + " tier '" + tier.name + "' dose " +
                                 std::to_string(tier.dose) + " but vaccine '" +
                                 camp.vaccine_type + "' has " +
                                 std::to_string(doses) + " doses (0 to " +
                                 std::to_string(doses - 1) + ")");
    }
  }
}

}  // namespace

VaccinationConfig ConfigLoader::loadVaccination(const std::string& filename) {
  VaccinationConfig config;
  // No vaccines_file in config_paths means no vaccination. A path that was
  // written down but leads nowhere is a mistake, not a choice.
  if (filename.empty()) {
    config.enabled = false;
    return config;
  }
  if (!std::filesystem::exists(filename)) {
    throw std::runtime_error("Vaccination file not found: " + filename +
                             ". simulation.yaml's `vaccines_file` must point "
                             "at an existing YAML, or be left out to run "
                             "without vaccination.");
  }
  try {
    YAML::Node root = YAML::LoadFile(filename);
    if (!root) return config;
    rejectUnknownKeys(root, {"enabled", "vaccines", "campaigns"}, "");

    if (root["enabled"]) {
      config.enabled = root["enabled"].as<bool>();
    } else {
      config.enabled = true;  // Fallback to true if file exists but no flag
    }

    if (root["vaccines"]) {
      for (const auto& vacc_kv : root["vaccines"]) {
        std::string name = vacc_kv.first.as<std::string>();
        config.vaccines[name] = parseVaccine(name, vacc_kv.second);
      }
    }

    if (root["campaigns"]) {
      const auto yaml_dir = std::filesystem::path(filename).parent_path();
      for (const auto& camp_kv : root["campaigns"]) {
        config.campaigns.push_back(parseVaccinationCampaign(
            camp_kv.first.as<std::string>(), camp_kv.second, yaml_dir));
      }
      if (config.campaigns.size() > kMaxEventIds)
        throw std::runtime_error("the file has " +
                                 std::to_string(config.campaigns.size()) +
                                 " campaigns; at most " +
                                 std::to_string(kMaxEventIds) + " are allowed");
    }
    checkTierDoses(config);

  } catch (const config_detail::MalformedSelectionError&) {
    throw;
  } catch (const std::exception& e) {
    throw std::runtime_error("Error loading vaccination file '" + filename +
                             "': " + e.what());
  }
  return config;
}

}  // namespace june
