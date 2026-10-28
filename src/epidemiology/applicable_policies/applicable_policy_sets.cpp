#include "epidemiology/applicable_policies/applicable_policy_sets.h"

#include <stdexcept>
#include <string>

namespace june {

namespace {

std::vector<uint8_t> setBits(const PolicyMask& mask) {
  std::vector<uint8_t> indices;
  for (size_t index = 0; index < kMaxPoliciesPerKind; ++index)
    if (mask.test(index)) indices.push_back(static_cast<uint8_t>(index));
  return indices;
}

}  // namespace

ApplicablePolicySets::ApplicablePolicySets() {
  intern(PolicyMask{}, PolicyMask{});
}

uint32_t ApplicablePolicySets::intern(const PolicyMask& symptom_mask,
                                      const PolicyMask& temporal_mask) {
  Key key;
  for (size_t word = 0; word < kPolicyMaskWords; ++word) {
    key[word] = policyMaskWord(symptom_mask, word);
    key[kPolicyMaskWords + word] = policyMaskWord(temporal_mask, word);
  }
  auto it = ids_.find(key);
  if (it != ids_.end()) return it->second;
  const PolicyMask* masks[kPolicyKindCount] = {&symptom_mask, &temporal_mask};
  const char* kind_names[kPolicyKindCount] = {"symptom", "temporal"};
  for (size_t kind = 0; kind < kPolicyKindCount; ++kind)
    if (masks[kind]->count() > kMaxApplicablePoliciesPerKind)
      throw std::runtime_error("applicable policy set holds " +
                               std::to_string(masks[kind]->count()) + " " +
                               kind_names[kind] + " policies; at most " +
                               std::to_string(kMaxApplicablePoliciesPerKind) +
                               " of each kind may apply to one person");
  Entry entry;
  for (size_t kind = 0; kind < kPolicyKindCount; ++kind) {
    entry.indices[kind] = setBits(*masks[kind]);
    entry.positions[kind].fill(kAbsent);
    for (size_t position = 0; position < entry.indices[kind].size(); ++position)
      entry.positions[kind][entry.indices[kind][position]] =
          static_cast<uint8_t>(position);
  }
  const uint32_t set_id = static_cast<uint32_t>(entries_.size());
  entries_.push_back(std::move(entry));
  ids_.emplace(key, set_id);
  return set_id;
}

PolicyMask ApplicablePolicySets::globalMask(uint32_t set_id,
                                            PolicyKind kind) const {
  PolicyMask mask;
  for (uint8_t global_index : indices(set_id, kind)) mask.set(global_index);
  return mask;
}

}  // namespace june
