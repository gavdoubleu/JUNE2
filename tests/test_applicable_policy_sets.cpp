#include <doctest.h>

#include "core/types.h"
#include "epidemiology/applicable_policies/applicable_policy_sets.h"

using namespace june;

namespace {

PolicyMask maskOf(std::initializer_list<size_t> global_indices) {
  PolicyMask mask;
  for (size_t index : global_indices) mask.set(index);
  return mask;
}

}  // namespace

TEST_CASE("applicable policy sets start with the empty set at id 0") {
  ApplicablePolicySets sets;
  CHECK(sets.intern(PolicyMask{}, PolicyMask{}) ==
        ApplicablePolicySets::kEmptySet);
  CHECK(sets.indices(ApplicablePolicySets::kEmptySet, PolicyKind::Symptom)
            .empty());
  CHECK(sets.indices(ApplicablePolicySets::kEmptySet, PolicyKind::Temporal)
            .empty());
}

TEST_CASE("interning the same masks twice returns the same id") {
  ApplicablePolicySets sets;
  const auto first = sets.intern(maskOf({3}), maskOf({0, 64}));
  CHECK(first != ApplicablePolicySets::kEmptySet);
  CHECK(sets.intern(maskOf({3}), maskOf({0, 64})) == first);
  CHECK(sets.intern(maskOf({0, 64}), maskOf({3})) != first);
}

TEST_CASE("a set's indices ascend and position lookup agrees with them") {
  ApplicablePolicySets sets;
  const auto set_id = sets.intern(maskOf({}), maskOf({127, 5, 64, 31, 32}));
  const auto indices = sets.indices(set_id, PolicyKind::Temporal);
  REQUIRE(indices.size() == 5);
  const std::vector<uint8_t> expected{5, 31, 32, 64, 127};
  CHECK(std::vector<uint8_t>(indices.begin(), indices.end()) == expected);
  for (size_t position = 0; position < indices.size(); ++position)
    CHECK(sets.position(set_id, PolicyKind::Temporal, indices[position]) ==
          position);
  CHECK(sets.position(set_id, PolicyKind::Temporal, 6) ==
        ApplicablePolicySets::kAbsent);
  CHECK(sets.position(set_id, PolicyKind::Symptom, 5) ==
        ApplicablePolicySets::kAbsent);
}

TEST_CASE("a set's global masks round trip through interning") {
  ApplicablePolicySets sets;
  const auto symptom_mask = maskOf({0, 63, 100});
  const auto temporal_mask = maskOf({1, 64, 127});
  const auto set_id = sets.intern(symptom_mask, temporal_mask);
  CHECK(sets.globalMask(set_id, PolicyKind::Symptom) == symptom_mask);
  CHECK(sets.globalMask(set_id, PolicyKind::Temporal) == temporal_mask);
  CHECK(sets.globalMask(ApplicablePolicySets::kEmptySet, PolicyKind::Symptom)
            .none());
}

TEST_CASE("a set holds at most 32 policies of each kind") {
  ApplicablePolicySets sets;
  PolicyMask thirty_two;
  for (size_t index = 0; index < 32; ++index) thirty_two.set(index * 4);
  CHECK_NOTHROW(sets.intern(thirty_two, thirty_two));

  PolicyMask thirty_three = thirty_two;
  thirty_three.set(1);
  CHECK_THROWS_WITH_AS(sets.intern(PolicyMask{}, thirty_three),
                       doctest::Contains("33 temporal policies"),
                       std::runtime_error);
  CHECK_THROWS_WITH_AS(sets.intern(thirty_three, PolicyMask{}),
                       doctest::Contains("33 symptom policies"),
                       std::runtime_error);
}

TEST_CASE("position masks expand to global masks and compact back") {
  ApplicablePolicySets sets;
  const auto set_id = sets.intern(maskOf({2}), maskOf({0, 63, 64, 127}));
  // Positions 1 and 3 of {0, 63, 64, 127}.
  const uint32_t position_mask = 0b1010;
  const auto global_mask =
      sets.expand(set_id, PolicyKind::Temporal, position_mask);
  CHECK(global_mask == maskOf({63, 127}));
  CHECK(sets.compact(set_id, PolicyKind::Temporal, global_mask) ==
        position_mask);
  CHECK(sets.compact(set_id, PolicyKind::Symptom, maskOf({2})) == 1u);
  CHECK(sets.expand(ApplicablePolicySets::kEmptySet, PolicyKind::Symptom, 0)
            .none());
  CHECK(sets.compact(ApplicablePolicySets::kEmptySet, PolicyKind::Symptom,
                     PolicyMask{}) == 0u);
}

TEST_CASE("compacting a global mask throws on a policy outside the set") {
  ApplicablePolicySets sets;
  const auto set_id = sets.intern(maskOf({2}), maskOf({0, 64}));
  CHECK_THROWS_WITH_AS(
      sets.compact(set_id, PolicyKind::Temporal, maskOf({0, 65})),
      doctest::Contains("temporal policy 65"), std::runtime_error);
  CHECK_THROWS_WITH_AS(sets.compact(ApplicablePolicySets::kEmptySet,
                                    PolicyKind::Symptom, maskOf({2})),
                       doctest::Contains("symptom policy 2"),
                       std::runtime_error);

  PolicyMask thirty_two;
  for (size_t index = 0; index < 32; ++index) thirty_two.set(index * 4);
  const auto full_set_id = sets.intern(thirty_two, PolicyMask{});
  PolicyMask thirty_three = thirty_two;
  thirty_three.set(1);
  CHECK(sets.compact(full_set_id, PolicyKind::Symptom, thirty_two) ==
        0xFFFFFFFFu);
  CHECK_THROWS_AS(sets.compact(full_set_id, PolicyKind::Symptom, thirty_three),
                  std::runtime_error);
}
