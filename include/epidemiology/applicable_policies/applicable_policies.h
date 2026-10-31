#pragma once

#include <cstdint>
#include <vector>

#include "core/types.h"
#include "epidemiology/applicable_policies/applicable_policy_sets.h"
#include "epidemiology/applicable_policies/applicable_policy_walk.h"

namespace june {

// One kind's Applicable Policy State as global-index masks.
struct KindPolicyMasks {
  PolicyMask applicable;
  PolicyMask decided;
  PolicyMask participating;
  bool operator==(const KindPolicyMasks&) const = default;
};

// A Person's whole Applicable Policy State as global-index masks: the
// checkpoint layout, independent of the per-rank set table.
struct ApplicablePolicyMasks {
  KindPolicyMasks symptom;
  KindPolicyMasks temporal;
  KindPolicyMasks& of(PolicyKind kind) {
    return kind == PolicyKind::Symptom ? symptom : temporal;
  }
  const KindPolicyMasks& of(PolicyKind kind) const {
    return kind == PolicyKind::Symptom ? symptom : temporal;
  }
  bool operator==(const ApplicablePolicyMasks&) const = default;
};

// Every read or write of a Person's Applicable Policy State. Owns the per-rank
// Applicable Policy Set table and the run's base seed, so a walk can never pair
// a table or seed with the wrong Person. The seed's only owner: PolicyManager's
// own draws read it here.
class ApplicablePolicies {
 public:
  void setBaseSeed(uint64_t seed) { base_seed_ = seed; }
  uint64_t baseSeed() const { return base_seed_; }
  // How many policies of each kind are loaded. A walk indexes the policies by
  // every position in the set unchecked, so import refuses any index past
  // these.
  void setLoadedCounts(size_t symptom_count, size_t temporal_count) {
    loaded_counts_[ApplicablePolicyState::slot(PolicyKind::Symptom)] =
        symptom_count;
    loaded_counts_[ApplicablePolicyState::slot(PolicyKind::Temporal)] =
        temporal_count;
  }

  // A Person's Applicable Policies of one kind, in ascending global index.
  // Query walk from a const module or for a const Person; the relocating walk
  // needs both mutable, so read-only holders cannot change decisions.
  template <typename Policy>
  ApplicablePolicyWalk<Policy, const Person> walk(
      const Person& person, const std::vector<Policy>& policies) const {
    return ApplicablePolicyWalk<Policy, const Person>(sets_, person, policies,
                                                      base_seed_);
  }
  template <typename Policy>
  ApplicablePolicyWalk<Policy, Person> walk(
      Person& person, const std::vector<Policy>& policies) {
    return ApplicablePolicyWalk<Policy, Person>(sets_, person, policies,
                                                base_seed_);
  }
  // The walk keeps a span into the policies, so a temporary would dangle.
  template <typename Policy, typename PersonT>
  void walk(PersonT& person, const std::vector<Policy>&& policies) = delete;
  template <typename Policy, typename PersonT>
  void walk(PersonT& person,
            const std::vector<Policy>&& policies) const = delete;

  // Forget a Person's decisions of one kind (e.g. symptom ones on recovery).
  // Stores only when something is held: getOverride clears every uninfected
  // Person each slot, and must not dirty their cache lines for nothing.
  void clear(Person& person, PolicyKind kind) {
    ApplicablePolicyState& state = person.applicable_policy_state;
    const size_t slot = ApplicablePolicyState::slot(kind);
    if ((state.decided_[slot] | state.participating_[slot]) == 0) return;
    state.decided_[slot] = 0;
    state.participating_[slot] = 0;
  }

  ApplicablePolicyMasks exportState(const Person& person) const;
  // Replace a Person's whole Applicable Policy State. Throws, naming the
  // Person, and changes nothing, if either kind exceeds
  // kMaxApplicablePoliciesPerKind or names a policy past those loaded, or a
  // decision falls outside the set, or a policy participates undecided.
  // Checkpoint restore remaps stored masks to the current policy list by name
  // and intersects them with the fresh set, so it never trips these; they
  // guard the invariant.
  void importState(Person& person, const ApplicablePolicyMasks& masks);

 private:
  ApplicablePolicySets sets_;
  uint64_t base_seed_ = 0;
  size_t loaded_counts_[kPolicyKindCount] = {};
};

}  // namespace june
