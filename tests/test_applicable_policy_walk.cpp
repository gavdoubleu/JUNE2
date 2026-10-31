#include <doctest.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/applicable_policies/applicable_policies.h"
#include "epidemiology/policy.h"
#include "policy_test_helpers.h"

using namespace june;
using june::policy_test::importApplicable;
using june::policy_test::maskOf;

namespace {

template <typename Policy>
std::vector<Policy> namedPolicies(size_t count) {
  std::vector<Policy> policies(count);
  for (size_t index = 0; index < count; ++index)
    policies[index].name = "policy " + std::to_string(index);
  return policies;
}

template <typename Policy>
std::vector<Policy> policiesComplyingAt(std::initializer_list<double> rates) {
  auto policies = namedPolicies<Policy>(rates.size());
  size_t index = 0;
  for (double rate : rates) policies[index++].action.compliance_rate = rate;
  return policies;
}

// The Applicable Policy at a global index, which the test knows is in the set.
template <typename Policy, typename PersonT>
auto applicableAt(ApplicablePolicies& applicable_policies, PersonT& person,
                  const std::vector<Policy>& policies, size_t global_index) {
  auto applicable =
      applicable_policies.walk(person, policies).find(global_index);
  REQUIRE(applicable.has_value());
  return *applicable;
}

// A module accepting every global index, for tests not about loading.
ApplicablePolicies allLoaded() {
  ApplicablePolicies applicable_policies;
  applicable_policies.setLoadedCounts(kMaxPoliciesPerKind, kMaxPoliciesPerKind);
  return applicable_policies;
}

bool holdsNoDecision(const ApplicablePolicies& applicable_policies,
                     const Person& person) {
  const ApplicablePolicyMasks masks = applicable_policies.exportState(person);
  for (PolicyKind kind : {PolicyKind::Symptom, PolicyKind::Temporal})
    if (masks.of(kind).decided.any() || masks.of(kind).participating.any())
      return false;
  return true;
}

}  // namespace

TEST_CASE("the walk yields a Person's applicable policies in ascending order") {
  ApplicablePolicies applicable_policies = allLoaded();
  const auto temporal_policies = namedPolicies<TemporalPolicy>(100);
  Person person;
  importApplicable(applicable_policies, person, PolicyMask{},
                   maskOf({64, 5, 99, 31}));

  std::vector<size_t> indices;
  for (const auto& applicable :
       applicable_policies.walk(person, temporal_policies)) {
    indices.push_back(applicable.index());
    CHECK(&applicable.policy() == &temporal_policies[applicable.index()]);
  }
  CHECK(indices == std::vector<size_t>{5, 31, 64, 99});
}

TEST_CASE("a compliance key never collides across policy kinds") {
  ApplicablePolicies applicable_policies = allLoaded();
  const auto symptom_policies = namedPolicies<SymptomPolicy>(10);
  const auto temporal_policies = namedPolicies<TemporalPolicy>(10);
  Person person;
  importApplicable(applicable_policies, person, maskOf({3}), maskOf({3}));

  for (const auto& applicable :
       applicable_policies.walk(person, symptom_policies))
    CHECK(applicable.complianceKey() == 3);
  for (const auto& applicable :
       applicable_policies.walk(person, temporal_policies))
    CHECK(applicable.complianceKey() == 3 + kTemporalComplianceIndexOffset);
}

TEST_CASE("looking up one policy finds it only inside the Person's set") {
  ApplicablePolicies applicable_policies = allLoaded();
  const auto temporal_policies = namedPolicies<TemporalPolicy>(10);
  const Person person = [&] {
    Person person;
    importApplicable(applicable_policies, person, PolicyMask{}, maskOf({2, 7}));
    return person;
  }();
  const auto walk = applicable_policies.walk(person, temporal_policies);

  const auto found = walk.find(7);
  REQUIRE(found.has_value());
  CHECK(&found->policy() == &temporal_policies[7]);
  CHECK(found->index() == 7);
  CHECK(found->complianceKey() == 7 + kTemporalComplianceIndexOffset);

  CHECK_FALSE(walk.find(3).has_value());
}

TEST_CASE("a latched decision sticks when the compliance rate changes") {
  ApplicablePolicies applicable_policies = allLoaded();
  auto temporal_policies = policiesComplyingAt<TemporalPolicy>({1.0, 0.0});
  Person person;
  importApplicable(applicable_policies, person, PolicyMask{}, maskOf({0, 1}));

  CHECK(
      applicableAt(applicable_policies, person, temporal_policies, 0).latch());
  CHECK_FALSE(
      applicableAt(applicable_policies, person, temporal_policies, 1).latch());

  temporal_policies[0].action.compliance_rate = 0.0;
  temporal_policies[1].action.compliance_rate = 1.0;
  CHECK(
      applicableAt(applicable_policies, person, temporal_policies, 0).latch());
  CHECK(applicableAt(applicable_policies, person, temporal_policies, 0)
            .participates());
  CHECK_FALSE(
      applicableAt(applicable_policies, person, temporal_policies, 1).latch());
  CHECK_FALSE(applicableAt(applicable_policies, person, temporal_policies, 1)
                  .participates());

  SUBCASE("clearing forgets the decision, so the next latch redraws") {
    applicableAt(applicable_policies, person, temporal_policies, 0).clear();
    CHECK_FALSE(applicableAt(applicable_policies, person, temporal_policies, 0)
                    .latch());
    CHECK_FALSE(applicableAt(applicable_policies, person, temporal_policies, 1)
                    .latch());
  }
}

TEST_CASE("a voluntary policy participates without latching a decision") {
  ApplicablePolicies applicable_policies = allLoaded();
  auto temporal_policies = policiesComplyingAt<TemporalPolicy>({0.0});
  temporal_policies[0].trigger = PolicyTrigger{};
  temporal_policies[0].trigger->terms.emplace_back();
  Person person;
  importApplicable(applicable_policies, person, PolicyMask{}, maskOf({0}));

  const auto applicable =
      applicableAt(applicable_policies, person, temporal_policies, 0);
  CHECK(applicable.latch());
  CHECK(applicable.participates());
  CHECK_FALSE(applicable.complies());
  CHECK(holdsNoDecision(applicable_policies, person));
}

TEST_CASE("the query walk never changes a Person's Sticky Compliance") {
  ApplicablePolicies applicable_policies = allLoaded();
  const auto symptom_policies = policiesComplyingAt<SymptomPolicy>({1.0, 0.0});
  const auto temporal_policies =
      policiesComplyingAt<TemporalPolicy>({1.0, 0.0, 1.0});
  Person mutable_person;
  importApplicable(applicable_policies, mutable_person, maskOf({0, 1}),
                   maskOf({0, 1, 2}));
  applicableAt(applicable_policies, mutable_person, temporal_policies, 2)
      .latch();
  const Person& person = mutable_person;
  const ApplicablePolicyMasks before = applicable_policies.exportState(person);

  std::vector<bool> symptom_answers;
  for (const auto& applicable :
       applicable_policies.walk(person, symptom_policies))
    symptom_answers.push_back(applicable.participates());
  std::vector<bool> temporal_answers;
  for (const auto& applicable :
       applicable_policies.walk(person, temporal_policies))
    temporal_answers.push_back(applicable.participates());

  CHECK(symptom_answers == std::vector<bool>{true, false});
  CHECK(temporal_answers == std::vector<bool>{true, false, true});
  CHECK(applicable_policies.exportState(person) == before);
}

TEST_CASE("clearing symptom decisions leaves temporal ones") {
  ApplicablePolicies applicable_policies = allLoaded();
  const auto symptom_policies = policiesComplyingAt<SymptomPolicy>({1.0});
  auto temporal_policies = policiesComplyingAt<TemporalPolicy>({1.0});
  Person person;
  importApplicable(applicable_policies, person, maskOf({0}), maskOf({0}));
  applicableAt(applicable_policies, person, symptom_policies, 0).latch();
  applicableAt(applicable_policies, person, temporal_policies, 0).latch();

  applicable_policies.clear(person, PolicyKind::Symptom);

  temporal_policies[0].action.compliance_rate = 0.0;
  CHECK(
      applicableAt(applicable_policies, person, temporal_policies, 0).latch());
  applicableAt(applicable_policies, person, temporal_policies, 0).clear();
  CHECK(holdsNoDecision(applicable_policies, person));
}
