#include <doctest.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/applicable_policies/applicable_policy_walk.h"

using namespace june;

namespace {

PolicyMask maskOf(std::initializer_list<size_t> global_indices) {
  PolicyMask mask;
  for (size_t index : global_indices) mask.set(index);
  return mask;
}

template <typename Policy>
std::vector<Policy> namedPolicies(size_t count) {
  std::vector<Policy> policies(count);
  for (size_t index = 0; index < count; ++index)
    policies[index].name = "policy " + std::to_string(index);
  return policies;
}

}  // namespace

TEST_CASE("the walk yields a Person's applicable policies in ascending order") {
  ApplicablePolicySets sets;
  const auto temporal_policies = namedPolicies<TemporalPolicy>(100);
  Person person;
  person.applicable_policy_set_id =
      sets.intern(PolicyMask{}, maskOf({64, 5, 99, 31}));

  std::vector<size_t> indices;
  std::vector<size_t> positions;
  for (const auto& applicable :
       walkApplicablePolicies(sets, person, temporal_policies)) {
    indices.push_back(applicable.index());
    positions.push_back(applicable.position());
    CHECK(&applicable.policy() == &temporal_policies[applicable.index()]);
  }
  CHECK(indices == std::vector<size_t>{5, 31, 64, 99});
  CHECK(positions == std::vector<size_t>{0, 1, 2, 3});
}

TEST_CASE("a compliance key never collides across policy kinds") {
  ApplicablePolicySets sets;
  const auto symptom_policies = namedPolicies<SymptomPolicy>(10);
  const auto temporal_policies = namedPolicies<TemporalPolicy>(10);
  Person person;
  person.applicable_policy_set_id = sets.intern(maskOf({3}), maskOf({3}));

  for (const auto& applicable :
       walkApplicablePolicies(sets, person, symptom_policies))
    CHECK(applicable.complianceKey() == 3);
  for (const auto& applicable :
       walkApplicablePolicies(sets, person, temporal_policies))
    CHECK(applicable.complianceKey() == 3 + kTemporalComplianceIndexOffset);
}

TEST_CASE("looking up one policy finds it only inside the Person's set") {
  ApplicablePolicySets sets;
  const auto temporal_policies = namedPolicies<TemporalPolicy>(10);
  const Person person = [&] {
    Person person;
    person.applicable_policy_set_id = sets.intern(PolicyMask{}, maskOf({2, 7}));
    return person;
  }();
  const auto walk = walkApplicablePolicies(sets, person, temporal_policies);

  const auto found = walk.find(7);
  REQUIRE(found.has_value());
  CHECK(&found->policy() == &temporal_policies[7]);
  CHECK(found->index() == 7);
  CHECK(found->position() == 1);
  CHECK(found->complianceKey() == 7 + kTemporalComplianceIndexOffset);

  CHECK_FALSE(walk.find(3).has_value());
}

TEST_CASE("a policy index past the loaded policies is refused at set time") {
  WorldState world;
  PolicyManager policy_manager(world);
  policy_manager.addTemporalPolicy(TemporalPolicy{});
  Person person;
  person.id = 42;

  CHECK_NOTHROW(
      policy_manager.setApplicablePolicies(person, PolicyMask{}, maskOf({0})));
  std::string message;
  try {
    policy_manager.setApplicablePolicies(person, PolicyMask{}, maskOf({0, 1}));
  } catch (const std::runtime_error& error) {
    message = error.what();
  }
  CHECK(message.find("person 42") != std::string::npos);
  CHECK(message.find("temporal") != std::string::npos);
  CHECK(policy_manager.isApplicable(person, PolicyKind::Temporal, 0));
  CHECK_FALSE(policy_manager.isApplicable(person, PolicyKind::Temporal, 1));

  CHECK_THROWS_AS(
      policy_manager.setApplicablePolicies(person, maskOf({0}), PolicyMask{}),
      std::runtime_error);
}
