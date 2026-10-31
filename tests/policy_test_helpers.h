#pragma once

// Helpers the observed metric and Policy trigger tests share.

#include <filesystem>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

#include "epidemiology/disease.h"
#include "epidemiology/policy.h"
#include "loaders/policy_loader.h"

namespace june::policy_test {

// Message of the exception a call throws, or "" if it does not throw.
inline std::string thrown(const std::function<void()>& f) {
  try {
    f();
  } catch (const std::exception& e) {
    return e.what();
  }
  return "";
}

inline bool has(const std::string& text, const std::string& part) {
  return text.find(part) != std::string::npos;
}

inline PolicyMask maskOf(std::initializer_list<size_t> global_indices) {
  PolicyMask mask;
  for (size_t index : global_indices) mask.set(index);
  return mask;
}

// Set the policies whose criteria `person` meets, as global-index masks, with
// no decisions held. Normally set by precomputePolicyApplicability.
inline void importApplicable(ApplicablePolicies& applicable_policies,
                             Person& person, const PolicyMask& symptom_mask,
                             const PolicyMask& temporal_mask) {
  ApplicablePolicyMasks masks;
  masks.symptom.applicable = symptom_mask;
  masks.temporal.applicable = temporal_mask;
  applicable_policies.importState(person, masks);
}

inline void importApplicable(PolicyManager& pm, Person& person,
                             const PolicyMask& symptom_mask,
                             const PolicyMask& temporal_mask) {
  importApplicable(pm.applicablePolicies(), person, symptom_mask,
                   temporal_mask);
}

// `person`'s Applicable Policy State as global-index masks.
inline ApplicablePolicyMasks exportedState(const PolicyManager& pm,
                                           const Person& person) {
  return pm.applicablePolicies().exportState(person);
}

// Temporal policy 0 applies to every Person in `world`. Normally set by
// precomputePolicyApplicability.
inline void applyTemporalPolicyZero(PolicyManager& pm, WorldState& world) {
  for (Person& person : world.people)
    importApplicable(pm, person, PolicyMask{}, 1);
}

// Every policy `pm` holds applies to every Person in `world`. Normally set by
// precomputePolicyApplicability.
inline void applyEveryPolicy(PolicyManager& pm, WorldState& world) {
  PolicyMask symptom_mask;
  PolicyMask temporal_mask;
  for (size_t i = 0; i < pm.getSymptomPolicyCount(); ++i) symptom_mask.set(i);
  for (size_t i = 0; i < pm.getTemporalPolicyCount(); ++i) temporal_mask.set(i);
  for (Person& person : world.people)
    importApplicable(pm, person, symptom_mask, temporal_mask);
}

// healthy(0) -> mild(1) -> dead(2), with the days each stage lasts. Only here
// so the metrics can resolve.
inline Disease makeMildDeadDisease(const std::string& name, double mild_days,
                                   double dead_days) {
  TransmissionParams trans;
  trans.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(1.0);
  trans.symptom_id_curves = {nullptr, curve, curve};
  std::vector<SymptomTag> tags = {
      {"healthy", -1, 0}, {"mild", 1, 1}, {"dead", 2, 2}};
  DiseaseStageSettings settings;
  settings.fatality_stages = {"dead"};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", mild_days}}}});
  td.stages.push_back({"dead", {"constant", {{"value", dead_days}}}});
  return Disease(name, tags, settings, {td}, {}, trans);
}

// Loads `text` into `pm`; returns the thrown message or "".
inline std::string loadPolicyText(PolicyManager& pm, const std::string& text) {
  const auto path =
      std::filesystem::temp_directory_path() / "june_policy_test.yaml";
  {
    std::ofstream out(path);
    out << text;
  }
  const std::string msg = thrown(
      [&] { PolicyLoader::loadPolicies(pm, path.string(), "2020-02-03"); });
  std::filesystem::remove(path);
  return msg;
}

}  // namespace june::policy_test
