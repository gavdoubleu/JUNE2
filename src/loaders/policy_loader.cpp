#include "loaders/policy_loader.h"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "loaders/config_loader_detail.h"
#include "utils/time_utils.h"

namespace june {

namespace PolicyLoader {

namespace {

using ::june::config_detail::rejectUnknownKeys;

}  // namespace

void loadPolicies(PolicyManager& policy_manager, const std::string& filename,
                  const std::string& simulation_start_date) {
  try {
    YAML::Node root = YAML::LoadFile(filename);
    rejectUnknownKeys(root, {"observed_metrics", "policies"}, filename + ": ");

    if (root["observed_metrics"]) {
      loadObservedMetrics(policy_manager, root["observed_metrics"], filename);
    }

    if (!root["policies"]) {
      if (!root["observed_metrics"]) {
        std::cout << "No policies section found in " << filename << std::endl;
      }
      return;
    }

    const YAML::Node& policies = root["policies"];
    rejectUnknownKeys(policies, {"symptom_policies", "temporal_policies"},
                      filename + ": ", " in 'policies'");

    if (policies["symptom_policies"]) {
      loadSymptomPolicies(policy_manager, policies["symptom_policies"],
                          simulation_start_date, filename);
    }

    if (policies["temporal_policies"]) {
      loadTemporalPolicies(policy_manager, policies["temporal_policies"],
                           simulation_start_date, filename);
    }
  } catch (const YAML::Exception& e) {
    throw std::runtime_error("Error loading policy file '" + filename +
                             "': " + e.what());
  }
}

void loadObservedMetrics(PolicyManager& policy_manager, const YAML::Node& node,
                         const std::string& policies_filename) {
  if (!node.IsMap()) {
    throw std::runtime_error("observed_metrics in '" + policies_filename +
                             "' must be a mapping of metric name to fields.");
  }
  for (const auto& kv : node) {
    const std::string name = kv.first.as<std::string>();
    const YAML::Node& body = kv.second;
    if (!body.IsMap()) {
      throw std::runtime_error("observed_metrics '" + name + "' in '" +
                               policies_filename + "' must be a mapping.");
    }
    rejectUnknownKeys(
        body, {"source", "scope", "activity", "stage"},
        "observed_metrics '" + name + "' in '" + policies_filename + "': ");
    ObservedMetricConfig cfg;
    for (const auto& field : body) {
      const std::string key = field.first.as<std::string>();
      const std::string value = field.second.as<std::string>();
      if (key == "source") cfg.source = value;
      if (key == "scope") cfg.scope = value;
      if (key == "activity") cfg.activity = value;
      if (key == "stage") cfg.stage = value;
    }
    policy_manager.observedMetrics().declare(name, cfg);
  }
}

PolicyAction loadPolicyAction(const YAML::Node& node) {
  PolicyAction action;
  const auto parseStringSet = [](const YAML::Node& value) {
    std::unordered_set<std::string> result;
    if (value.IsSequence()) {
      for (const auto& item : value) {
        result.insert(item.as<std::string>());
      }
    } else {
      result.insert(value.as<std::string>());
    }
    return result;
  };

  if (node["override_activities"]) {
    action.override_activities = parseStringSet(node["override_activities"]);
  }

  if (node["override_venue_types"]) {
    action.override_venue_types = parseStringSet(node["override_venue_types"]);
  }

  if (node["exempt_venue_types"]) {
    action.exempt_venue_types = parseStringSet(node["exempt_venue_types"]);
  }

  if (node["replacement"]) {
    action.replacement_activity = node["replacement"].as<std::string>();
  } else {
    action.replacement_activity = "residence";
  }

  if (node["replacement_schedule"]) {
    action.replacement_schedule =
        node["replacement_schedule"].as<std::string>();
  }

  if (node["exempt"]) {
    if (node["exempt"].IsSequence()) {
      for (const auto& exempt_node : node["exempt"]) {
        rejectUnknownKeys(exempt_node, {"activity", "selection"},
                          "policy: ", " in an 'exempt' entry");
        ActivityExemption exemption;
        exemption.activity_name = exempt_node["activity"].as<std::string>();
        if (exempt_node["selection"].IsDefined()) {
          config_detail::parseSelectionCriteria(exempt_node["selection"],
                                                exemption.criteria);
        }
        action.exemptions.push_back(exemption);
      }
    }
  }

  if (node["compliance_rate"]) {
    action.compliance_rate = node["compliance_rate"].as<double>();
  }

  return action;
}

bool hasExplicitLocationAction(const YAML::Node& node) {
  return node["override_activities"] || node["override_venue_types"] ||
         node["exempt_venue_types"] || node["replacement"] ||
         node["replacement_schedule"] || node["exempt"];
}

ActiveWindow loadActiveWindow(const YAML::Node& node,
                              const std::string& simulation_start_date) {
  ActiveWindow window;
  const std::string policy_name =
      node["name"] ? node["name"].as<std::string>() : std::string("<unnamed>");

  if (node["start_date"] && node["start_time"]) {
    throw std::runtime_error("policy '" + policy_name +
                             "' declares both start_date and start_time; "
                             "give one");
  }
  if (node["end_date"] && node["end_time"]) {
    throw std::runtime_error("policy '" + policy_name +
                             "' declares both end_date and end_time; give one");
  }

  if (node["start_date"]) {
    std::tm sim_start_tm = parseDate(simulation_start_date);
    std::tm policy_start_tm = parseDate(node["start_date"].as<std::string>());
    window.start_time =
        static_cast<double>(daysBetween(sim_start_tm, policy_start_tm));
  } else if (node["start_time"]) {
    window.start_time = node["start_time"].as<double>();
  }

  if (node["end_date"]) {
    std::tm sim_start_tm = parseDate(simulation_start_date);
    std::tm policy_end_tm = parseDate(node["end_date"].as<std::string>());
    window.end_time =
        static_cast<double>(daysBetween(sim_start_tm, policy_end_tm));
  } else if (node["end_time"]) {
    window.end_time = node["end_time"].as<double>();
  }

  const bool start_stated = node["start_date"] || node["start_time"];
  if (start_stated && window.end_time &&
      *window.end_time <= window.start_time) {
    throw std::runtime_error("policy '" + policy_name + "' ends on day " +
                             std::to_string(*window.end_time) +
                             ", at or before its start on day " +
                             std::to_string(window.start_time));
  }

  return window;
}

namespace {

void loadPolicyName(std::string& name, const YAML::Node& node,
                    const char* policy_kind) {
  if (!node["name"]) {
    throw std::runtime_error(std::string(policy_kind) +
                             " must have a 'name' field");
  }
  name = node["name"].as<std::string>();
}

// Every policy entry takes these fields; a symptom policy adds a few of its
// own. A misspelt key would otherwise leave its field at the default, so
// complience_rate would quietly mean full compliance.
void rejectUnknownPolicyKeys(const YAML::Node& node, const std::string& kind,
                             const std::string& name,
                             const std::string& filename) {
  const std::string where = filename + ": ";
  const std::string within = " in " + kind + " '" + name + "'";
  if (kind == "symptom_policy") {
    rejectUnknownKeys(
        node,
        {"name", "start_date", "start_time", "end_date", "end_time",
         "override_activities", "override_venue_types", "exempt_venue_types",
         "replacement", "replacement_schedule", "exempt", "compliance_rate",
         "applies_to", "trigger", "transmission_effects_file", "symptoms",
         "follow_up_policy", "inherit_compliance", "inherit_refusal"},
        where, within);
  } else {
    rejectUnknownKeys(
        node,
        {"name", "start_date", "start_time", "end_date", "end_time",
         "override_activities", "override_venue_types", "exempt_venue_types",
         "replacement", "replacement_schedule", "exempt", "compliance_rate",
         "applies_to", "trigger", "transmission_effects_file"},
        where, within);
  }
}

template <typename Policy>
void loadCommonPolicyFields(Policy& policy, const YAML::Node& node,
                            const std::string& simulation_start_date) {
  policy.window = loadActiveWindow(node, simulation_start_date);
  policy.action = loadPolicyAction(node);
  if (node["applies_to"].IsDefined()) {
    config_detail::parseSelectionCriteria(node["applies_to"],
                                          policy.applies_to);
  }
}

int64_t readWholeNumber(const YAML::Node& value, const std::string& where,
                        const std::string& key) {
  try {
    return value.as<int64_t>();
  } catch (const YAML::Exception&) {
    throw std::runtime_error(where + "'" + key + "' must be a whole number.");
  }
}

double readNumber(const YAML::Node& value, const std::string& where,
                  const std::string& key) {
  try {
    const double number = value.as<double>();
    if (std::isfinite(number)) return number;
  } catch (const YAML::Exception&) {
  }
  throw std::runtime_error(where + "'" + key + "' must be a finite number.");
}

// The response form of a trigger: a base level, terms of metric and weight, and
// the decision unit. The level is base + sum(weight * count), read at the day
// boundary and clamped to 1 where it is used. Base is held to 0 to 1 and every
// weight above 0, so the level never needs a lower clamp.
PolicyTrigger loadResponse(const YAML::Node& node,
                           const ObservedMetrics& metrics,
                           const std::string& where) {
  for (const char* key : {"metric", "activate_at", "deactivate_at", "calm_days",
                          "tolerance", "tolerance_pct"}) {
    if (node[key]) {
      throw std::runtime_error(where +
                               "takes either the threshold fields or a "
                               "'response' list, never both; found '" +
                               key + "' next to 'response'.");
    }
  }
  const YAML::Node response = node["response"];
  if (!response.IsMap()) {
    throw std::runtime_error(where +
                             "'response' must be a mapping with 'base' and "
                             "'terms'.");
  }
  rejectUnknownKeys(response, {"base", "terms"}, where, " in 'response'");

  PolicyTrigger trigger;
  if (response["base"]) {
    trigger.base = readNumber(response["base"], where, "base");
    if (trigger.base < 0.0 || trigger.base > 1.0) {
      throw std::runtime_error(where + "'base' must be from 0 to 1, got " +
                               std::to_string(trigger.base) + ".");
    }
  }
  trigger.world_level = trigger.base;

  if (!response["terms"] || !response["terms"].IsSequence() ||
      response["terms"].size() == 0) {
    throw std::runtime_error(where +
                             "'response' needs a non empty 'terms' list of "
                             "{metric, weight}.");
  }
  for (const auto& term_node : response["terms"]) {
    if (!term_node.IsMap()) {
      throw std::runtime_error(where +
                               "each term must be a mapping with 'metric' and "
                               "'weight'.");
    }
    rejectUnknownKeys(term_node, {"metric", "weight"}, where, " in a term");
    if (!term_node["metric"] || !term_node["weight"]) {
      throw std::runtime_error(where +
                               "each term needs both 'metric' and 'weight'.");
    }
    ResponseTerm term;
    term.metric = term_node["metric"].as<std::string>();
    metrics.requireDeclared(term.metric, where);
    term.weight = readNumber(term_node["weight"], where, "weight");
    if (term.weight <= 0.0) {
      throw std::runtime_error(where + "the weight of '" + term.metric +
                               "' must be above 0, got " +
                               std::to_string(term.weight) + ".");
    }
    for (const auto& earlier : trigger.terms) {
      if (earlier.metric == term.metric) {
        throw std::runtime_error(where + "metric '" + term.metric +
                                 "' appears in two terms; give it one weight.");
      }
    }
    trigger.terms.push_back(term);
  }

  // One scope and one activity for all terms. Mixing them in one response is
  // not supported yet.
  const std::string& first = trigger.terms[0].metric;
  trigger.scope = metrics.scope(first);
  for (const auto& term : trigger.terms) {
    if (metrics.scope(term.metric) != trigger.scope ||
        metrics.activityOf(term.metric) != metrics.activityOf(first)) {
      throw std::runtime_error(
          where + "metrics '" + first + "' and '" + term.metric +
          "' differ in scope or activity. A response reads one scope and one "
          "activity; mixing them is not built yet.");
    }
  }

  if (!node["decision_unit"]) {
    throw std::runtime_error(where +
                             "a 'response' needs 'decision_unit': 'person' or "
                             "'{shared_venue: <activity>}'.");
  }
  const YAML::Node unit = node["decision_unit"];
  if (unit.IsScalar() && unit.as<std::string>() == "person") {
    // each person decides alone
  } else if (unit.IsMap() && unit.size() == 1 && unit["shared_venue"] &&
             unit["shared_venue"].IsScalar()) {
    trigger.decision_unit.shared_activity =
        unit["shared_venue"].as<std::string>();
    if (trigger.decision_unit.shared_activity.empty()) {
      throw std::runtime_error(where +
                               "'shared_venue' needs an activity name.");
    }
  } else {
    throw std::runtime_error(where +
                             "'decision_unit' must be 'person' or "
                             "'{shared_venue: <activity>}'.");
  }
  return trigger;
}

// Parse and check one policy's trigger block. Everything that needs only the
// YAML and the declared metrics is checked here, so a bad trigger fails at
// load with the policy and metric named.
PolicyTrigger loadTrigger(const YAML::Node& node,
                          const ObservedMetrics& metrics,
                          const std::string& label) {
  const std::string where = label + " trigger: ";
  if (!node.IsMap()) throw std::runtime_error(where + "must be a mapping.");

  rejectUnknownKeys(node,
                    {"metric", "activate_at", "deactivate_at", "calm_days",
                     "tolerance", "tolerance_pct", "response", "decision_unit"},
                    where);
  if (node["response"]) return loadResponse(node, metrics, where);
  if (node["decision_unit"]) {
    throw std::runtime_error(where +
                             "'decision_unit' belongs to the 'response' form; "
                             "a threshold trigger has no decision unit.");
  }
  if (!node["metric"])
    throw std::runtime_error(where + "'metric' is required.");
  const std::string metric = node["metric"].as<std::string>();
  metrics.requireDeclared(metric, where);

  const auto number = [&](const char* key) -> std::optional<int64_t> {
    if (!node[key]) return std::nullopt;
    return readWholeNumber(node[key], where, key);
  };
  const auto activate_at = number("activate_at");
  const auto deactivate_at = number("deactivate_at");
  const auto calm_days = number("calm_days");
  const auto tolerance = number("tolerance");
  const auto tolerance_pct = number("tolerance_pct");

  const bool has_threshold =
      activate_at || deactivate_at || calm_days || tolerance || tolerance_pct;
  if (!has_threshold) {
    throw std::runtime_error(where +
                             "needs the threshold fields ('activate_at' and "
                             "the fields that go with it) or a 'response' "
                             "list.");
  }

  PolicyTrigger trigger;
  trigger.metric = metric;
  trigger.cumulative = metrics.isCumulative(metric);
  trigger.scope = metrics.scope(metric);
  if (!activate_at) {
    throw std::runtime_error(where + "'activate_at' is required.");
  }
  if (*activate_at < 1) {
    throw std::runtime_error(where + "'activate_at' must be at least 1, got " +
                             std::to_string(*activate_at) + ".");
  }
  trigger.activate_at = *activate_at;

  if (trigger.cumulative) {
    if (deactivate_at) {
      throw std::runtime_error(where + "'deactivate_at' is not valid for the " +
                               "cumulative metric '" + metric +
                               "'; it reopens after 'calm_days' without " +
                               "growth.");
    }
    if (!calm_days) {
      throw std::runtime_error(where + "'calm_days' is required for the " +
                               "cumulative metric '" + metric + "'.");
    }
    if (*calm_days < 1) {
      throw std::runtime_error(where + "'calm_days' must be at least 1, got " +
                               std::to_string(*calm_days) + ".");
    }
    trigger.calm_days = *calm_days;
    if (tolerance) {
      if (*tolerance < 0) {
        throw std::runtime_error(where + "'tolerance' must not be negative.");
      }
      trigger.tolerance = *tolerance;
    }
    if (tolerance_pct) {
      if (*tolerance_pct < 0 || *tolerance_pct > 100) {
        throw std::runtime_error(
            where + "'tolerance_pct' must be a whole percent from 0 to 100.");
      }
      trigger.tolerance_pct = *tolerance_pct;
    }
  } else {
    if (!deactivate_at) {
      throw std::runtime_error(where + "'deactivate_at' is required for the " +
                               "current metric '" + metric + "'.");
    }
    if (calm_days || tolerance || tolerance_pct) {
      throw std::runtime_error(where +
                               "'calm_days', 'tolerance' and 'tolerance_pct' " +
                               "are only valid for a cumulative metric, and " +
                               "'" + metric + "' is a current one.");
    }
    if (*deactivate_at < 0) {
      throw std::runtime_error(where + "'deactivate_at' must not be negative.");
    }
    if (*deactivate_at >= trigger.activate_at) {
      throw std::runtime_error(where + "'deactivate_at' (" +
                               std::to_string(*deactivate_at) +
                               ") must be below 'activate_at' (" +
                               std::to_string(trigger.activate_at) + ").");
    }
    trigger.deactivate_at = *deactivate_at;
  }
  return trigger;
}

// A Policy whose trigger counts venues of one activity must be allowed to act
// on that activity: an empty list or '*' means all.
void checkActivityListed(const std::string& where, const PolicyAction& action,
                         const ObservedMetrics& metrics,
                         const std::string& metric) {
  const auto& listed = action.override_activities;
  const std::string& activity = metrics.activityOf(metric);
  if (!listed.empty() && !listed.count("*") && !listed.count(activity)) {
    throw std::runtime_error(
        where + "'override_activities' must include '" + activity +
        "', the activity of metric '" + metric +
        "'; otherwise the policy never acts on the venues the metric counts.");
  }
}

// A voluntary Policy derives who avoids from the metric every time it is asked.
// The compliance rate it replaces must be absent, and whatever assumes a
// latched decision is refused rather than half honoured.
void checkVoluntaryPolicy(const std::string& kind, const std::string& name,
                          const YAML::Node& node, const PolicyTrigger& trigger,
                          const ObservedMetrics& metrics,
                          const PolicyAction& action) {
  const std::string where = kind + "_policy '" + name + "' trigger: ";
  if (node["compliance_rate"]) {
    throw std::runtime_error(
        where +
        "'compliance_rate' must be absent on a policy with a 'response', "
        "because the avoidance level replaces it. Remove the key, even a "
        "value of 1.");
  }
  if (node["transmission_effects_file"]) {
    throw std::runtime_error(
        where +
        "transmission effects on a policy with a 'response' are not built "
        "yet; the policy can only send people elsewhere.");
  }
  if (kind == "symptom") {
    if (trigger.decision_unit.shared()) {
      throw std::runtime_error(
          where +
          "'shared_venue' is not built yet for a symptom policy; use "
          "'decision_unit: person'.");
    }
    for (const char* key : {"replacement_schedule", "follow_up_policy"}) {
      if (node[key]) {
        throw std::runtime_error(
            where + "'" + key +
            "' is not built yet on a policy with a 'response': it assumes a "
            "decision that stays latched, and a voluntary policy keeps none.");
      }
    }
  } else if (trigger.perUnit()) {
    checkActivityListed(where, action, metrics, trigger.terms[0].metric);
  }
}

// A Policy triggered per venue closes a venue for everyone scheduled there.
// Whatever would make that partial, or would need venue by venue rebuilding of
// the transmission modifiers, is refused here rather than half honoured.
void checkVenueTriggeredPolicy(const std::string& where, const YAML::Node& node,
                               const PolicyTrigger& trigger,
                               const ObservedMetrics& metrics,
                               const PolicyAction& action) {
  if (action.compliance_rate != 1.0) {
    throw std::runtime_error(
        where +
        "'compliance_rate' must be 1 on a policy triggered per venue, "
        "because a closed venue closes for everyone scheduled there; "
        "got " +
        std::to_string(action.compliance_rate) + ".");
  }
  if (node["transmission_effects_file"]) {
    throw std::runtime_error(
        where +
        "transmission effects on a policy triggered per venue are not "
        "built yet; the policy can only send people elsewhere.");
  }
  checkActivityListed(where, action, metrics, trigger.metric);
}

// Load the trigger of one policy, if it has one, and check what the form of
// the trigger rules out. `kind` is "symptom" or "temporal".
template <typename Policy>
void loadPolicyTrigger(Policy& policy, const YAML::Node& node,
                       const ObservedMetrics& metrics,
                       const std::string& kind) {
  if (!node["trigger"]) return;
  const std::string label = kind + "_policy '" + policy.name + "'";
  policy.trigger = loadTrigger(node["trigger"], metrics, label);
  const PolicyTrigger& trigger = *policy.trigger;
  // A threshold trigger on a world scope metric has nothing more to check.
  if (!trigger.decidesPerSlot()) return;

  const std::string where = label + " trigger: ";
  if (kind == "symptom" && trigger.perUnit()) {
    throw std::runtime_error(
        where + "metric '" +
        (trigger.voluntary() ? trigger.terms[0].metric : trigger.metric) +
        "' is counted per venue, which is not built yet for a symptom "
        "policy. Use a world scope metric, or write a temporal policy.");
  }
  if (trigger.voluntary()) {
    checkVoluntaryPolicy(kind, policy.name, node, trigger, metrics,
                         policy.action);
  } else {
    checkVenueTriggeredPolicy(where, node, trigger, metrics, policy.action);
  }
}

}  // namespace

void loadSymptomPolicies(PolicyManager& policy_manager, const YAML::Node& node,
                         const std::string& simulation_start_date,
                         const std::string& policies_filename) {
  if (!node.IsSequence()) {
    throw std::runtime_error("symptom_policies must be a list");
  }

  for (const auto& policy_node : node) {
    SymptomPolicy policy;

    loadPolicyName(policy.name, policy_node, "symptom_policy");
    rejectUnknownPolicyKeys(policy_node, "symptom_policy", policy.name,
                            policies_filename);

    if (!policy_node["symptoms"]) {
      throw std::runtime_error("symptom_policy '" + policy.name +
                               "' must have 'symptoms' field");
    }
    policy.trigger_symptoms =
        policy_node["symptoms"].as<std::vector<std::string>>();

    loadCommonPolicyFields(policy, policy_node, simulation_start_date);
    loadPolicyTrigger(policy, policy_node, policy_manager.observedMetrics(),
                      "symptom");

    if (policy_node["transmission_effects_file"]) {
      if (!hasExplicitLocationAction(policy_node))
        policy.action.has_location_override = false;
      loadTransmissionEffects(
          policy_manager,
          policy_node["transmission_effects_file"].as<std::string>(),
          TransmissionPolicyKind::Symptom,
          static_cast<uint16_t>(policy_manager.getSymptomPolicyCount()),
          policy.action.compliance_rate, policy.name, policies_filename);
    }

    if (policy_node["follow_up_policy"]) {
      policy.follow_up_policy_name =
          policy_node["follow_up_policy"].as<std::string>();
    }
    if (policy_node["inherit_compliance"]) {
      policy.inherit_compliance = policy_node["inherit_compliance"].as<bool>();
    }
    if (policy_node["inherit_refusal"]) {
      policy.inherit_refusal = policy_node["inherit_refusal"].as<bool>();
    }

    policy_manager.addSymptomPolicy(policy);
  }
}

void loadTemporalPolicies(PolicyManager& policy_manager, const YAML::Node& node,
                          const std::string& simulation_start_date,
                          const std::string& policies_filename) {
  if (!node.IsSequence()) {
    throw std::runtime_error("temporal_policies must be a list");
  }

  for (const auto& policy_node : node) {
    TemporalPolicy policy;

    loadPolicyName(policy.name, policy_node, "temporal_policy");
    rejectUnknownPolicyKeys(policy_node, "temporal_policy", policy.name,
                            policies_filename);
    loadCommonPolicyFields(policy, policy_node, simulation_start_date);
    loadPolicyTrigger(policy, policy_node, policy_manager.observedMetrics(),
                      "temporal");

    if (policy_node["transmission_effects_file"]) {
      if (!hasExplicitLocationAction(policy_node))
        policy.action.has_location_override = false;
      loadTransmissionEffects(
          policy_manager,
          policy_node["transmission_effects_file"].as<std::string>(),
          TransmissionPolicyKind::Temporal,
          static_cast<uint16_t>(policy_manager.getTemporalPolicyCount()),
          policy.action.compliance_rate, policy.name, policies_filename);
    }

    policy_manager.addTemporalPolicy(policy);
  }
}

}  // namespace PolicyLoader

}  // namespace june
