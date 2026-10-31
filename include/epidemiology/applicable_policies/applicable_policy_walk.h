#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <type_traits>
#include <vector>

#include "core/types.h"
#include "epidemiology/applicable_policies/applicable_policy_sets.h"
#include "utils/deterministic_rng.h"

namespace june {

struct SymptomPolicy;
struct TemporalPolicy;

// Temporal policy i keys its compliance draw on i + this offset, so no
// temporal key collides with a symptom policy's key (its index).
inline constexpr uint32_t kTemporalComplianceIndexOffset =
    static_cast<uint32_t>(kMaxPoliciesPerKind);

// Does this Person comply with a policy of this compliance rate? A pure draw:
// SplitMix64 over (base_seed, person_id, compliance_key), no stored RNG state,
// so the same triple always answers the same way. That is what lets a query
// reach a verdict for an undecided Person without latching the decision.
inline bool drawsCompliance(uint64_t base_seed, PersonId person_id,
                            uint32_t compliance_key, double compliance_rate) {
  SplitMix64 rng(mix_seed(base_seed, person_id, compliance_key, 0xC0E1A9ULL));
  return uniform01(rng) < compliance_rate;
}

// A voluntary policy keeps no latched decision: whether the Person avoids is
// derived from the draw and the level every time it is asked.
template <typename Policy>
bool isVoluntary(const Policy& policy) {
  return policy.trigger && policy.trigger->voluntary();
}

// Left undefined for any other type, so a walk over one fails to compile
// rather than reading another kind's decisions.
template <typename Policy>
struct PolicyKindOf;
template <>
struct PolicyKindOf<SymptomPolicy> {
  static constexpr PolicyKind value = PolicyKind::Symptom;
};
template <>
struct PolicyKindOf<TemporalPolicy> {
  static constexpr PolicyKind value = PolicyKind::Temporal;
};
template <typename Policy>
inline constexpr PolicyKind kPolicyKindOf = PolicyKindOf<Policy>::value;

// One Applicable Policy: a policy of a Person's Applicable Policy Set, seen
// through that Person. PersonT is const Person for a query walk, which only
// reads its Applicable Policy State, and Person for the relocating walk, which
// also latches, clears and hands it on.
template <typename Policy, typename PersonT>
class ApplicablePolicy {
  static constexpr bool kRelocating = !std::is_const_v<PersonT>;
  static constexpr size_t kSlot = static_cast<size_t>(kPolicyKindOf<Policy>);

 public:
  ApplicablePolicy(const Policy& policy, PersonT& person,
                   const ApplicablePolicySets& sets, uint64_t base_seed,
                   size_t index, size_t position)
      : policy_(&policy),
        person_(&person),
        sets_(&sets),
        base_seed_(base_seed),
        index_(index),
        position_(position) {}

  const Policy& policy() const { return *policy_; }
  // Global index into the PolicyManager's policies of this kind.
  size_t index() const { return index_; }
  // Where the policy sits in the Person's set.
  size_t position() const { return position_; }
  // Key of this Person's compliance draw on the policy, unique across kinds.
  uint32_t complianceKey() const {
    const auto key = static_cast<uint32_t>(index_);
    return kPolicyKindOf<Policy> == PolicyKind::Temporal
               ? key + kTemporalComplianceIndexOffset
               : key;
  }

  // Does the Person comply? The latched answer if decided, else a redraw that
  // stores nothing. The redraw is exactly what latch() would store, so a query
  // that beats the relocating walk to a Person still agrees with it.
  bool complies() const {
    const ApplicablePolicyState& state = person_->applicable_policy_state;
    if (state.decided_[kSlot] & bit())
      return state.participating_[kSlot] & bit();
    return drawsCompliance(base_seed_, person_->id, complianceKey(),
                           policy_->action.compliance_rate);
  }
  // Is the policy acting on the Person? Always for a voluntary policy, which
  // asks per slot instead; otherwise whether they comply.
  bool participates() const { return isVoluntary(*policy_) || complies(); }

  // participates(), latching the compliance draw the first time it is asked.
  bool latch() const
    requires kRelocating
  {
    if (isVoluntary(*policy_)) return true;
    ApplicablePolicyState& state = person_->applicable_policy_state;
    if (!(state.decided_[kSlot] & bit())) {
      if (drawsCompliance(base_seed_, person_->id, complianceKey(),
                          policy_->action.compliance_rate))
        state.participating_[kSlot] |= bit();
      state.decided_[kSlot] |= bit();
    }
    return state.participating_[kSlot] & bit();
  }
  // Forget the Person's decision on this policy.
  void clear() const
    requires kRelocating
  {
    ApplicablePolicyState& state = person_->applicable_policy_state;
    state.decided_[kSlot] &= ~bit();
    state.participating_[kSlot] &= ~bit();
  }
  // Hand a latched decision on to the follow-up policy at a global index:
  // compliance if inherit_compliance, refusal if inherit_refusal. Inheriting
  // compliance marks the follow-up decided and complying; inheriting refusal
  // only marks it decided. A follow-up outside the Person's set never applies
  // to them, so there is nothing to inherit.
  void inheritTo(size_t follow_up_index, bool inherit_compliance,
                 bool inherit_refusal) const
    requires kRelocating
  {
    ApplicablePolicyState& state = person_->applicable_policy_state;
    if (!(state.decided_[kSlot] & bit())) return;
    const uint8_t follow_up_position =
        sets_->position(person_->applicable_policy_state.set_id_,
                        kPolicyKindOf<Policy>, follow_up_index);
    if (follow_up_position == ApplicablePolicySets::kAbsent) return;
    const PolicyPositionMask follow_up_bit = PolicyPositionMask{1}
                                             << follow_up_position;
    const bool complying = state.participating_[kSlot] & bit();
    if (complying && inherit_compliance) {
      state.participating_[kSlot] |= follow_up_bit;
      state.decided_[kSlot] |= follow_up_bit;
    } else if (!complying && inherit_refusal) {
      state.decided_[kSlot] |= follow_up_bit;
    }
  }

 private:
  PolicyPositionMask bit() const { return PolicyPositionMask{1} << position_; }

  const Policy* policy_;
  PersonT* person_;
  const ApplicablePolicySets* sets_;
  uint64_t base_seed_;
  size_t index_;
  size_t position_;
};

// The Applicable Policies of one kind for one Person, in ascending global
// index. Every index the set holds must address `policies`;
// ApplicablePolicies::importState refuses any that does not, and loaded
// policies are never removed. Built by ApplicablePolicies::walk.
template <typename Policy, typename PersonT>
class ApplicablePolicyWalk {
 public:
  using Entry = ApplicablePolicy<Policy, PersonT>;
  static constexpr PolicyKind kKind = kPolicyKindOf<Policy>;

  ApplicablePolicyWalk(const ApplicablePolicySets& sets, PersonT& person,
                       std::span<const Policy> policies, uint64_t base_seed)
      : sets_(&sets),
        person_(&person),
        policies_(policies),
        base_seed_(base_seed),
        indices_(sets.indices(person.applicable_policy_state.set_id_, kKind)) {}

  class iterator {
   public:
    iterator(const ApplicablePolicyWalk* walk, size_t position)
        : walk_(walk), position_(position) {}
    Entry operator*() const {
      const size_t index = walk_->indices_[position_];
      return Entry(walk_->policies_[index], *walk_->person_, *walk_->sets_,
                   walk_->base_seed_, index, position_);
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
    const uint8_t position = sets_->position(
        person_->applicable_policy_state.set_id_, kKind, global_index);
    if (position == ApplicablePolicySets::kAbsent) return std::nullopt;
    return Entry(policies_[global_index], *person_, *sets_, base_seed_,
                 global_index, position);
  }

 private:
  const ApplicablePolicySets* sets_;
  PersonT* person_;
  std::span<const Policy> policies_;
  uint64_t base_seed_;
  std::span<const uint8_t> indices_;
};

}  // namespace june
