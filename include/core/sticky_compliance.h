#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

namespace june {

enum class PolicyKind : uint8_t { Symptom, Temporal };
inline constexpr size_t kPolicyKindCount = 2;

template <typename Policy, typename PersonT>
class ApplicablePolicy;
class PolicyManager;

// A Person's Sticky Compliance: for each policy of their Applicable Policy Set,
// whether they have decided to comply and, if so, what they decided. Words are
// indexed by set position, which only the Applicable Policy walk and
// PolicyManager's global-index helpers (decidedPolicies, setPolicyState, ...)
// translate.
class StickyCompliance {
 public:
  // Forget every decision of one kind (e.g. symptom ones on recovery).
  void clear(PolicyKind kind) {
    decided_[slot(kind)] = 0;
    participating_[slot(kind)] = 0;
  }
  bool empty() const {
    return (decided_[0] | decided_[1] | participating_[0] |
            participating_[1]) == 0;
  }
  bool operator==(const StickyCompliance&) const = default;

 private:
  template <typename Policy, typename PersonT>
  friend class ApplicablePolicy;
  friend class PolicyManager;

  static size_t slot(PolicyKind kind) { return static_cast<size_t>(kind); }

  std::array<uint32_t, kPolicyKindCount> decided_{};
  std::array<uint32_t, kPolicyKindCount> participating_{};
};

}  // namespace june
