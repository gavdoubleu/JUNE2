#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "core/types.h"
#include "epidemiology/applicable_policies/applicable_policy_sets.h"
#include "epidemiology/policy.h"

namespace june {

// Temporal policy i keys its compliance draw on i + this offset, so no
// temporal key collides with a symptom policy's key (its index).
inline constexpr uint32_t kTemporalComplianceIndexOffset =
    static_cast<uint32_t>(kMaxPoliciesPerKind);

template <typename Policy>
inline constexpr PolicyKind kPolicyKindOf = PolicyKind::Symptom;
template <>
inline constexpr PolicyKind kPolicyKindOf<TemporalPolicy> =
    PolicyKind::Temporal;

// One Applicable Policy: a policy of a Person's Applicable Policy Set, seen
// through that Person. PersonT is const Person for a query walk and Person for
// the relocating walk.
template <typename Policy, typename PersonT>
class ApplicablePolicy {
 public:
  ApplicablePolicy(const Policy& policy, PersonT& person, size_t index,
                   size_t position)
      : policy_(&policy),
        person_(&person),
        index_(index),
        position_(position) {}

  const Policy& policy() const { return *policy_; }
  // Global index into the PolicyManager's policies of this kind.
  size_t index() const { return index_; }
  // Where the policy sits in the Person's set: its bit in the decision masks.
  size_t position() const { return position_; }
  // Key of this Person's compliance draw on the policy, unique across kinds.
  uint32_t complianceKey() const {
    const auto key = static_cast<uint32_t>(index_);
    return kPolicyKindOf<Policy> == PolicyKind::Temporal
               ? key + kTemporalComplianceIndexOffset
               : key;
  }

 private:
  const Policy* policy_;
  PersonT* person_;
  size_t index_;
  size_t position_;
};

// The Applicable Policies of one kind for one Person, in ascending global
// index. Every index the set holds must address `policies`;
// PolicyManager::setApplicablePolicies refuses any that does not.
template <typename Policy, typename PersonT>
class ApplicablePolicyWalk {
 public:
  using Entry = ApplicablePolicy<Policy, PersonT>;
  static constexpr PolicyKind kKind = kPolicyKindOf<Policy>;

  ApplicablePolicyWalk(const ApplicablePolicySets& sets, PersonT& person,
                       std::span<const Policy> policies)
      : sets_(&sets),
        person_(&person),
        policies_(policies),
        indices_(sets.indices(person.applicable_policy_set_id, kKind)) {}

  class iterator {
   public:
    iterator(const ApplicablePolicyWalk* walk, size_t position)
        : walk_(walk), position_(position) {}
    Entry operator*() const {
      const size_t index = walk_->indices_[position_];
      return Entry(walk_->policies_[index], *walk_->person_, index, position_);
    }
    iterator& operator++() {
      ++position_;
      return *this;
    }
    bool operator!=(const iterator& other) const {
      return position_ != other.position_;
    }

   private:
    const ApplicablePolicyWalk* walk_;
    size_t position_;
  };

  iterator begin() const { return iterator(this, 0); }
  iterator end() const { return iterator(this, indices_.size()); }

  // The policy at a global index, if it is in the Person's set.
  std::optional<Entry> find(size_t global_index) const {
    const uint8_t position =
        sets_->position(person_->applicable_policy_set_id, kKind, global_index);
    if (position == ApplicablePolicySets::kAbsent) return std::nullopt;
    return Entry(policies_[global_index], *person_, global_index, position);
  }

 private:
  const ApplicablePolicySets* sets_;
  PersonT* person_;
  std::span<const Policy> policies_;
  std::span<const uint8_t> indices_;
};

// Query walk for a const Person, relocating walk for a mutable one.
template <typename Policy, typename PersonT>
ApplicablePolicyWalk<Policy, PersonT> walkApplicablePolicies(
    const ApplicablePolicySets& sets, PersonT& person,
    const std::vector<Policy>& policies) {
  return ApplicablePolicyWalk<Policy, PersonT>(sets, person, policies);
}

}  // namespace june
