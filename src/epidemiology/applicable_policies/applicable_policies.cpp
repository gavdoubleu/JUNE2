#include "epidemiology/applicable_policies/applicable_policies.h"

#include <stdexcept>
#include <string>

namespace june {

namespace {

// Rethrow a set table refusal naming the Person it was for.
template <typename Body>
void namingPerson(const Person& person, Body body) {
  try {
    body();
  } catch (const std::runtime_error& error) {
    throw std::runtime_error("person " + std::to_string(person.id) + ": " +
                             error.what());
  }
}

}  // namespace

ApplicablePolicyMasks ApplicablePolicies::exportState(
    const Person& person) const {
  const ApplicablePolicyState& state = person.applicable_policy_state;
  ApplicablePolicyMasks masks;
  for (PolicyKind kind : {PolicyKind::Symptom, PolicyKind::Temporal}) {
    const size_t slot = ApplicablePolicyState::slot(kind);
    masks.of(kind) = {
        sets_.globalMask(state.set_id_, kind),
        sets_.expand(state.set_id_, kind, state.decided_[slot]),
        sets_.expand(state.set_id_, kind, state.participating_[slot])};
  }
  return masks;
}

void ApplicablePolicies::importState(Person& person,
                                     const ApplicablePolicyMasks& masks) {
  // Build the whole state before writing, so a refused import leaves the
  // Person as they were.
  ApplicablePolicyState imported;
  namingPerson(person, [&] {
    for (PolicyKind kind : {PolicyKind::Symptom, PolicyKind::Temporal}) {
      const char* kind_name = policyKindName(kind);
      // Participation is the Sticky Compliance of a decided policy.
      const PolicyMask undecided_participation =
          masks.of(kind).participating & ~masks.of(kind).decided;
      for (size_t global_index = 0; global_index < kMaxPoliciesPerKind;
           ++global_index)
        if (undecided_participation.test(global_index))
          throw std::runtime_error(std::string(kind_name) + " policy " +
                                   std::to_string(global_index) +
                                   " participates without a decision");
      // Refused before interning, so a refusal leaves no orphan set in the
      // table. Participation is within the decisions, so this covers it too.
      const PolicyMask decided_outside =
          masks.of(kind).decided & ~masks.of(kind).applicable;
      for (size_t global_index = 0; global_index < kMaxPoliciesPerKind;
           ++global_index)
        if (decided_outside.test(global_index))
          throw std::runtime_error(
              std::string(kind_name) + " policy " +
              std::to_string(global_index) +
              " is outside the person's applicable policy set");
      const size_t loaded_count =
          loaded_counts_[ApplicablePolicyState::slot(kind)];
      if ((masks.of(kind).applicable >> loaded_count).any())
        throw std::runtime_error(std::string(kind_name) +
                                 " policy mask names a policy past the " +
                                 std::to_string(loaded_count) + " loaded");
    }
    imported.set_id_ =
        sets_.intern(masks.symptom.applicable, masks.temporal.applicable);
    for (PolicyKind kind : {PolicyKind::Symptom, PolicyKind::Temporal}) {
      const size_t slot = ApplicablePolicyState::slot(kind);
      imported.decided_[slot] =
          sets_.compact(imported.set_id_, kind, masks.of(kind).decided);
      imported.participating_[slot] =
          sets_.compact(imported.set_id_, kind, masks.of(kind).participating);
    }
  });
  person.applicable_policy_state = imported;
}

}  // namespace june
