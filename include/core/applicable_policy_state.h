#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace june {

enum class PolicyKind : uint8_t { Symptom, Temporal };
inline constexpr size_t kPolicyKindCount = 2;
inline constexpr const char* policyKindName(PolicyKind kind) {
  return kind == PolicyKind::Symptom ? "symptom" : "temporal";
}

// Sticky Compliance word of one kind: bit k is position k of the Person's
// Applicable Policy Set.
using PolicyPositionMask = uint32_t;
// One set holds at most this many policies of each kind, so every position
// has a bit.
inline constexpr size_t kMaxApplicablePoliciesPerKind = 32;
static_assert(kMaxApplicablePoliciesPerKind <=
                  std::numeric_limits<PolicyPositionMask>::digits,
              "raise PolicyPositionMask's width with the per-Person cap");

class ApplicablePolicies;
template <typename Policy, typename PersonT>
class ApplicablePolicy;
template <typename Policy, typename PersonT>
class ApplicablePolicyWalk;

// A Person's Applicable Policy State: their Applicable Policy Set id and their
// Sticky Compliance on it. Decision words are indexed by position in that set,
// so the two are held together and only the ApplicablePolicies module reads or
// writes them.
class ApplicablePolicyState {
 public:
  bool operator==(const ApplicablePolicyState&) const = default;

 private:
  friend class ApplicablePolicies;
  template <typename Policy, typename PersonT>
  friend class ApplicablePolicy;
  template <typename Policy, typename PersonT>
  friend class ApplicablePolicyWalk;

  static size_t slot(PolicyKind kind) { return static_cast<size_t>(kind); }

  uint32_t set_id_ = 0;
  std::array<PolicyPositionMask, kPolicyKindCount> decided_{};
  std::array<PolicyPositionMask, kPolicyKindCount> participating_{};
};

}  // namespace june
