#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "core/applicable_policy_state.h"
#include "core/types.h"

namespace june {

// Identity of a policy list: each kind's policy names in global-index order.
// A policy is its name within its kind, so a checkpoint's roster matched with
// the current one says where each stored policy now sits.
struct PolicyRoster {
  std::vector<std::string> symptom;
  std::vector<std::string> temporal;
  std::vector<std::string>& of(PolicyKind kind) {
    return kind == PolicyKind::Symptom ? symptom : temporal;
  }
  const std::vector<std::string>& of(PolicyKind kind) const {
    return kind == PolicyKind::Symptom ? symptom : temporal;
  }
  bool operator==(const PolicyRoster&) const = default;
};

// Names of each kind's policies, in order.
template <typename SymptomPolicies, typename TemporalPolicies>
PolicyRoster makePolicyRoster(const SymptomPolicies& symptom_policies,
                              const TemporalPolicies& temporal_policies) {
  PolicyRoster roster;
  for (const auto& policy : symptom_policies)
    roster.symptom.push_back(policy.name);
  for (const auto& policy : temporal_policies)
    roster.temporal.push_back(policy.name);
  return roster;
}

// Throws, naming the kind, the name and both indices, if a name repeats within
// a kind. The same name in both kinds is two different policies.
void rejectDuplicateNames(const PolicyRoster& roster);

// One kind's stored roster matched with the current one.
struct KindRosterMatch {
  struct Move {
    std::string name;
    size_t stored_index;
    size_t current_index;
    bool operator==(const Move&) const = default;
  };
  // Indexed by stored index: the policy's current index, or nullopt if gone.
  std::vector<std::optional<size_t>> current_index;
  // In stored order.
  std::vector<std::string> dropped;
  // In current order.
  std::vector<std::string> added;
  // Kept policies whose index changed, in stored order.
  std::vector<Move> moved;
};

// Where each stored policy sits in the current roster, by kind.
class PolicyRosterMatch {
 public:
  PolicyRosterMatch(KindRosterMatch symptom, KindRosterMatch temporal)
      : symptom_(std::move(symptom)), temporal_(std::move(temporal)) {}

  const KindRosterMatch& of(PolicyKind kind) const {
    return kind == PolicyKind::Symptom ? symptom_ : temporal_;
  }
  // True when neither kind dropped, added or moved a policy.
  bool identical() const;

  // A stored policy's current index, or nullopt if it is gone. Throws on an
  // index past the stored roster.
  std::optional<size_t> currentIndex(PolicyKind kind,
                                     size_t stored_index) const;
  // A global-index mask in stored order to current order; gone policies drop
  // out. Throws on a bit past the stored roster.
  PolicyMask remap(const PolicyMask& stored_mask, PolicyKind kind) const;

  // One line per changed kind; empty when identical.
  std::string describeChanges() const;

 private:
  KindRosterMatch symptom_;
  KindRosterMatch temporal_;
};

// Matches by name within each kind. Throws if either roster repeats a name.
PolicyRosterMatch match(const PolicyRoster& stored,
                        const PolicyRoster& current);

}  // namespace june
