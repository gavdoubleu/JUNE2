#include <doctest.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/applicable_policies/applicable_policy_walk.h"

using namespace june;

namespace {

constexpr uint64_t kSeed = 12345;

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

template <typename Policy>
std::vector<Policy> policiesComplyingAt(std::initializer_list<double> rates) {
  auto policies = namedPolicies<Policy>(rates.size());
  size_t index = 0;
  for (double rate : rates) policies[index++].action.compliance_rate = rate;
  return policies;
}

// The Applicable Policy at a global index, which the test knows is in the set.
template <typename Policy, typename PersonT>
auto applicableAt(const ApplicablePolicySets& sets, PersonT& person,
                  const std::vector<Policy>& policies, size_t global_index) {
  auto applicable =
      walkApplicablePolicies(sets, person, policies, kSeed).find(global_index);
  REQUIRE(applicable.has_value());
  return *applicable;
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
       walkApplicablePolicies(sets, person, temporal_policies, kSeed)) {
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
       walkApplicablePolicies(sets, person, symptom_policies, kSeed))
    CHECK(applicable.complianceKey() == 3);
  for (const auto& applicable :
       walkApplicablePolicies(sets, person, temporal_policies, kSeed))
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
  const auto walk =
      walkApplicablePolicies(sets, person, temporal_policies, kSeed);

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

TEST_CASE("a latched decision sticks when the compliance rate changes") {
  ApplicablePolicySets sets;
  auto temporal_policies = policiesComplyingAt<TemporalPolicy>({1.0, 0.0});
  Person person;
  person.applicable_policy_set_id = sets.intern(PolicyMask{}, maskOf({0, 1}));

  CHECK(applicableAt(sets, person, temporal_policies, 0).latch());
  CHECK_FALSE(applicableAt(sets, person, temporal_policies, 1).latch());

  temporal_policies[0].action.compliance_rate = 0.0;
  temporal_policies[1].action.compliance_rate = 1.0;
  CHECK(applicableAt(sets, person, temporal_policies, 0).latch());
  CHECK(applicableAt(sets, person, temporal_policies, 0).participates());
  CHECK_FALSE(applicableAt(sets, person, temporal_policies, 1).latch());
  CHECK_FALSE(applicableAt(sets, person, temporal_policies, 1).participates());

  SUBCASE("clearing forgets the decision, so the next latch redraws") {
    applicableAt(sets, person, temporal_policies, 0).clear();
    CHECK_FALSE(applicableAt(sets, person, temporal_policies, 0).latch());
    CHECK_FALSE(applicableAt(sets, person, temporal_policies, 1).latch());
  }
}

TEST_CASE("a follow-up inherits compliance, refusal, or nothing") {
  ApplicablePolicySets sets;
  auto symptom_policies = policiesComplyingAt<SymptomPolicy>({1.0, 0.0});
  Person person;
  person.applicable_policy_set_id = sets.intern(maskOf({0, 1}), PolicyMask{});

  SUBCASE("compliance") {
    const auto current = applicableAt(sets, person, symptom_policies, 0);
    REQUIRE(current.latch());
    current.inheritTo(1, true, false);
    CHECK(applicableAt(sets, person, symptom_policies, 1).latch());
  }
  SUBCASE("refusal") {
    symptom_policies[0].action.compliance_rate = 0.0;
    symptom_policies[1].action.compliance_rate = 1.0;
    const auto current = applicableAt(sets, person, symptom_policies, 0);
    REQUIRE_FALSE(current.latch());
    current.inheritTo(1, false, true);
    CHECK_FALSE(applicableAt(sets, person, symptom_policies, 1).latch());
  }
  SUBCASE("neither: compliance not inherited leaves the follow-up undecided") {
    const auto current = applicableAt(sets, person, symptom_policies, 0);
    REQUIRE(current.latch());
    current.inheritTo(1, false, true);
    symptom_policies[1].action.compliance_rate = 1.0;
    CHECK(applicableAt(sets, person, symptom_policies, 1).latch());
  }
  SUBCASE("an undecided policy hands nothing on") {
    const StickyCompliance before = person.sticky_compliance;
    applicableAt(sets, person, symptom_policies, 0).inheritTo(1, true, true);
    CHECK(person.sticky_compliance == before);
  }
}

TEST_CASE("a follow-up outside the Person's set inherits nothing") {
  ApplicablePolicySets sets;
  const auto symptom_policies = policiesComplyingAt<SymptomPolicy>({1.0, 1.0});
  Person person;
  person.applicable_policy_set_id = sets.intern(maskOf({0}), PolicyMask{});

  const auto current = applicableAt(sets, person, symptom_policies, 0);
  REQUIRE(current.latch());
  const StickyCompliance before = person.sticky_compliance;
  current.inheritTo(1, true, true);
  CHECK(person.sticky_compliance == before);
}

TEST_CASE("a voluntary policy participates without latching a decision") {
  ApplicablePolicySets sets;
  auto temporal_policies = policiesComplyingAt<TemporalPolicy>({0.0});
  temporal_policies[0].trigger = PolicyTrigger{};
  temporal_policies[0].trigger->terms.emplace_back();
  Person person;
  person.applicable_policy_set_id = sets.intern(PolicyMask{}, maskOf({0}));

  const auto applicable = applicableAt(sets, person, temporal_policies, 0);
  CHECK(applicable.latch());
  CHECK(applicable.participates());
  CHECK_FALSE(applicable.complies());
  CHECK(person.sticky_compliance.empty());
}

TEST_CASE("the query walk never changes a Person's Sticky Compliance") {
  ApplicablePolicySets sets;
  const auto symptom_policies = policiesComplyingAt<SymptomPolicy>({1.0, 0.0});
  const auto temporal_policies =
      policiesComplyingAt<TemporalPolicy>({1.0, 0.0, 1.0});
  Person mutable_person;
  mutable_person.applicable_policy_set_id =
      sets.intern(maskOf({0, 1}), maskOf({0, 1, 2}));
  applicableAt(sets, mutable_person, temporal_policies, 2).latch();
  const Person& person = mutable_person;
  const StickyCompliance before = person.sticky_compliance;

  std::vector<bool> symptom_answers;
  for (const auto& applicable :
       walkApplicablePolicies(sets, person, symptom_policies, kSeed))
    symptom_answers.push_back(applicable.participates());
  std::vector<bool> temporal_answers;
  for (const auto& applicable :
       walkApplicablePolicies(sets, person, temporal_policies, kSeed))
    temporal_answers.push_back(applicable.participates());

  CHECK(symptom_answers == std::vector<bool>{true, false});
  CHECK(temporal_answers == std::vector<bool>{true, false, true});
  CHECK(person.sticky_compliance == before);
}

TEST_CASE("clearing symptom decisions leaves temporal ones") {
  ApplicablePolicySets sets;
  const auto symptom_policies = policiesComplyingAt<SymptomPolicy>({1.0});
  auto temporal_policies = policiesComplyingAt<TemporalPolicy>({1.0});
  Person person;
  person.applicable_policy_set_id = sets.intern(maskOf({0}), maskOf({0}));
  applicableAt(sets, person, symptom_policies, 0).latch();
  applicableAt(sets, person, temporal_policies, 0).latch();

  person.sticky_compliance.clear(PolicyKind::Symptom);

  temporal_policies[0].action.compliance_rate = 0.0;
  CHECK(applicableAt(sets, person, temporal_policies, 0).latch());
  applicableAt(sets, person, temporal_policies, 0).clear();
  CHECK(person.sticky_compliance.empty());
}
