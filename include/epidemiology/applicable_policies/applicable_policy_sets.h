#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <span>
#include <vector>

#include "core/types.h"

namespace june {

enum class PolicyKind : uint8_t { Symptom, Temporal };
inline constexpr size_t kPolicyKindCount = 2;

// Per-Person decision and participation masks are uint32_t indexed by set
// position, so one set holds at most this many policies of each kind.
inline constexpr size_t kMaxApplicablePoliciesPerKind = 32;

// Per-rank table of the distinct sets of policies a Person's selection criteria
// match. Each Person keeps only a set id; few distinct sets exist because they
// follow demographic criteria. Id 0 is always the empty set.
class ApplicablePolicySets {
 public:
  static constexpr uint32_t kEmptySet = 0;
  // position() result for a policy outside the set.
  static constexpr uint8_t kAbsent = 0xFF;
  static_assert(kMaxPoliciesPerKind <= kAbsent,
                "global policy indices are stored as uint8_t");

  ApplicablePolicySets();

  // Throws if either mask holds more than kMaxApplicablePoliciesPerKind bits.
  uint32_t intern(const PolicyMask& symptom_mask,
                  const PolicyMask& temporal_mask);

  // Global indices of the set's policies of one kind, ascending.
  std::span<const uint8_t> indices(uint32_t set_id, PolicyKind kind) const {
    return entries_[set_id].indices[static_cast<size_t>(kind)];
  }

  // Where a global policy index sits in indices(), or kAbsent.
  uint8_t position(uint32_t set_id, PolicyKind kind,
                   size_t global_index) const {
    return entries_[set_id].positions[static_cast<size_t>(kind)][global_index];
  }

  // The set's policies of one kind as a global-index mask (checkpoint layout).
  PolicyMask globalMask(uint32_t set_id, PolicyKind kind) const;

  // Position mask (bit k = indices()[k]) to global-index mask.
  PolicyMask expand(uint32_t set_id, PolicyKind kind,
                    uint32_t position_mask) const;
  // Global-index mask to position mask. Throws on a policy outside the set.
  uint32_t compact(uint32_t set_id, PolicyKind kind,
                   const PolicyMask& global_mask) const;

 private:
  struct Entry {
    std::vector<uint8_t> indices[kPolicyKindCount];
    std::array<uint8_t, kMaxPoliciesPerKind> positions[kPolicyKindCount];
  };
  using Key = std::array<uint64_t, kPolicyKindCount * kPolicyMaskWords>;

  std::vector<Entry> entries_;
  std::map<Key, uint32_t> ids_;
};

}  // namespace june
