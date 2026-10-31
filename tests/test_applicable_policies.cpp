#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/applicable_policies/applicable_policies.h"
#include "epidemiology/policy.h"
#include "policy_test_helpers.h"

using namespace june;

using june::policy_test::has;
using june::policy_test::maskOf;

namespace {

void addTemporalPolicies(PolicyManager& manager, size_t count) {
  for (size_t index = 0; index < count; ++index) {
    TemporalPolicy policy;
    policy.name = "temporal " + std::to_string(index);
    manager.addTemporalPolicy(policy);
  }
}

void addSymptomPolicies(PolicyManager& manager, size_t count) {
  for (size_t index = 0; index < count; ++index) {
    SymptomPolicy policy;
    policy.name = "symptom " + std::to_string(index);
    manager.addSymptomPolicy(policy);
  }
}

// A const ApplicablePolicies only reads: no clear, and only the query walk even
// for a mutable Person. The relocating walk needs a mutable module.
template <typename Module>
concept ClearsDecisions =
    requires(Module& applicable_policies, Person& person) {
      applicable_policies.clear(person, PolicyKind::Symptom);
    };
static_assert(!ClearsDecisions<const ApplicablePolicies>);
static_assert(ClearsDecisions<ApplicablePolicies>);
static_assert(
    std::is_same_v<decltype(std::declval<const ApplicablePolicies&>().walk(
                       std::declval<Person&>(),
                       std::declval<const std::vector<TemporalPolicy>&>())),
                   ApplicablePolicyWalk<TemporalPolicy, const Person>>);
static_assert(
    std::is_same_v<decltype(std::declval<ApplicablePolicies&>().walk(
                       std::declval<Person&>(),
                       std::declval<const std::vector<TemporalPolicy>&>())),
                   ApplicablePolicyWalk<TemporalPolicy, Person>>);

// A walk keeps a span into the policies, so a temporary vector is refused.
template <typename Module, typename PersonT>
concept WalksTemporary =
    requires(Module& applicable_policies, PersonT& person) {
      applicable_policies.walk(person, std::vector<TemporalPolicy>{});
    };
static_assert(!WalksTemporary<ApplicablePolicies, Person>);
static_assert(!WalksTemporary<const ApplicablePolicies, Person>);
static_assert(!WalksTemporary<const ApplicablePolicies, const Person>);

// Only the two policy types have a kind; anything else fails to compile.
template <typename Policy>
concept HasPolicyKind = requires { PolicyKindOf<Policy>::value; };
static_assert(HasPolicyKind<SymptomPolicy>);
static_assert(HasPolicyKind<TemporalPolicy>);
static_assert(!HasPolicyKind<const TemporalPolicy>);

}  // namespace

TEST_CASE(
    "a query walk yields a Person's applicable policies in ascending "
    "global index") {
  WorldState world;
  PolicyManager manager(world);
  addTemporalPolicies(manager, 100);
  Person person;
  policy_test::importApplicable(manager, person, PolicyMask{},
                                maskOf({64, 5, 99, 31}));

  const Person& query_person = person;
  const auto& policies = manager.getTemporalPolicies();
  std::vector<size_t> indices;
  for (const auto& applicable :
       manager.applicablePolicies().walk(query_person, policies)) {
    indices.push_back(applicable.index());
    CHECK(&applicable.policy() == &policies[applicable.index()]);
  }
  CHECK(indices == std::vector<size_t>{5, 31, 64, 99});
}

TEST_CASE("the query walk draws compliance on the PolicyManager's base seed") {
  WorldState world;
  PolicyManager manager(world);
  TemporalPolicy policy;
  policy.action.compliance_rate = 0.5;
  manager.addTemporalPolicy(policy);
  constexpr uint64_t kSeed = 12345;
  manager.setBaseSeed(kSeed);

  for (PersonId person_id = 0; person_id < 64; ++person_id) {
    Person person;
    person.id = person_id;
    policy_test::importApplicable(manager, person, PolicyMask{}, maskOf({0}));
    const Person& query_person = person;
    for (const auto& applicable : manager.applicablePolicies().walk(
             query_person, manager.getTemporalPolicies()))
      CHECK(applicable.complies() ==
            drawsCompliance(kSeed, person_id, applicable.complianceKey(), 0.5));
  }
}

TEST_CASE("an added policy draws on its kind and name") {
  WorldState world;
  PolicyManager manager(world);
  SymptomPolicy symptom_policy;
  symptom_policy.name = "isolate";
  manager.addSymptomPolicy(symptom_policy);
  TemporalPolicy temporal_policy;
  temporal_policy.name = "isolate";
  manager.addTemporalPolicy(temporal_policy);
  CHECK(manager.getSymptomPolicies()[0].draw_key ==
        policyDrawKey(PolicyKind::Symptom, "isolate"));
  CHECK(manager.getTemporalPolicies()[0].draw_key ==
        policyDrawKey(PolicyKind::Temporal, "isolate"));
}

TEST_CASE("a relocating walk latches a decision once") {
  WorldState world;
  PolicyManager manager(world);
  TemporalPolicy policy;
  policy.action.compliance_rate = 1.0;
  manager.addTemporalPolicy(policy);
  Person person;
  policy_test::importApplicable(manager, person, PolicyMask{}, maskOf({0}));
  const auto latchPolicyZero = [&] {
    return manager.applicablePolicies()
        .walk(person, manager.getTemporalPolicies())
        .find(0)
        ->latch();
  };

  CHECK(latchPolicyZero());
  TemporalPolicy refused;
  refused.action.compliance_rate = 0.0;
  const std::vector<TemporalPolicy> rate_changed{refused};
  CHECK(manager.applicablePolicies().walk(person, rate_changed).find(0)->latch());
  CHECK(latchPolicyZero());
}

TEST_CASE("find answers only for a policy in the Person's set") {
  WorldState world;
  PolicyManager manager(world);
  addTemporalPolicies(manager, 10);
  Person person;
  policy_test::importApplicable(manager, person, PolicyMask{}, maskOf({2, 7}));
  const Person& query_person = person;
  const auto walk = manager.applicablePolicies().walk(
      query_person, manager.getTemporalPolicies());

  const auto found = walk.find(7);
  REQUIRE(found.has_value());
  CHECK(found->index() == 7);
  CHECK(&found->policy() == &manager.getTemporalPolicies()[7]);
  CHECK_FALSE(walk.find(3).has_value());
}

TEST_CASE("clearing symptom decisions forgets them and keeps temporal ones") {
  WorldState world;
  PolicyManager manager(world);
  SymptomPolicy symptom_policy;
  symptom_policy.action.compliance_rate = 1.0;
  manager.addSymptomPolicy(symptom_policy);
  TemporalPolicy temporal_policy;
  temporal_policy.action.compliance_rate = 1.0;
  manager.addTemporalPolicy(temporal_policy);
  Person person;
  policy_test::importApplicable(manager, person, maskOf({0}), maskOf({0}));
  auto& applicable_policies = manager.applicablePolicies();
  applicable_policies.walk(person, manager.getSymptomPolicies())
      .find(0)
      ->latch();
  applicable_policies.walk(person, manager.getTemporalPolicies())
      .find(0)
      ->latch();

  applicable_policies.clear(person, PolicyKind::Symptom);

  SymptomPolicy refused_symptom;
  refused_symptom.action.compliance_rate = 0.0;
  TemporalPolicy refused_temporal;
  refused_temporal.action.compliance_rate = 0.0;
  const std::vector<SymptomPolicy> symptom_rate_changed{refused_symptom};
  const std::vector<TemporalPolicy> temporal_rate_changed{refused_temporal};
  CHECK_FALSE(applicable_policies.walk(person, symptom_rate_changed)
                  .find(0)
                  ->latch());
  CHECK(applicable_policies.walk(person, temporal_rate_changed)
            .find(0)
            ->latch());
}

TEST_CASE("a Person with no applicable policies exports and imports empty") {
  WorldState world;
  PolicyManager manager(world);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;

  CHECK(applicable_policies.exportState(person) == ApplicablePolicyMasks{});
  applicable_policies.importState(person, ApplicablePolicyMasks{});
  CHECK(applicable_policies.exportState(person) == ApplicablePolicyMasks{});
}

TEST_CASE("applicable masks of both kinds round trip through import") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, 10);
  addTemporalPolicies(manager, 10);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  ApplicablePolicyMasks masks;
  masks.symptom.applicable = maskOf({1, 4});
  masks.temporal.applicable = maskOf({0, 9});

  applicable_policies.importState(person, masks);

  CHECK(applicable_policies.exportState(person) == masks);
}

TEST_CASE("imported decisions round trip and a later latch keeps them") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, 10);
  addTemporalPolicies(manager, 10);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  ApplicablePolicyMasks masks;
  masks.symptom = {maskOf({1, 4, 7}), maskOf({1, 4}), maskOf({1})};
  masks.temporal = {maskOf({2, 9}), maskOf({9}), PolicyMask{}};

  applicable_policies.importState(person, masks);

  CHECK(applicable_policies.exportState(person) == masks);
  SymptomPolicy complies;
  complies.action.compliance_rate = 1.0;
  const std::vector<SymptomPolicy> all_comply(10, complies);
  auto walk = applicable_policies.walk(person, all_comply);
  CHECK(walk.find(1)->latch());
  CHECK_FALSE(walk.find(4)->latch());
  CHECK(walk.find(7)->latch());
}

TEST_CASE("decisions on policies above global index 64 round trip") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, kMaxPoliciesPerKind);
  addTemporalPolicies(manager, kMaxPoliciesPerKind);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  ApplicablePolicyMasks masks;
  masks.symptom = {maskOf({3, 64, 100, 127}), maskOf({64, 127}), maskOf({127})};
  masks.temporal = {maskOf({63, 65, 126}), maskOf({63, 65}), maskOf({65})};

  applicable_policies.importState(person, masks);

  CHECK(applicable_policies.exportState(person) == masks);
}

TEST_CASE("a Person at the cap of both kinds round trips every decision") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, kMaxPoliciesPerKind);
  addTemporalPolicies(manager, kMaxPoliciesPerKind);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  ApplicablePolicyMasks masks;
  for (size_t position = 0; position < kMaxApplicablePoliciesPerKind;
       ++position) {
    const size_t symptom_index = position * 4;
    const size_t temporal_index = kMaxPoliciesPerKind - 1 - position * 3;
    masks.symptom.applicable.set(symptom_index);
    masks.temporal.applicable.set(temporal_index);
    masks.symptom.decided.set(symptom_index);
    masks.temporal.decided.set(temporal_index);
    if (position % 2 == 0) masks.symptom.participating.set(symptom_index);
    if (position % 3 == 0) masks.temporal.participating.set(temporal_index);
  }

  applicable_policies.importState(person, masks);

  CHECK(applicable_policies.exportState(person) == masks);
}

TEST_CASE("an import replaces decisions rather than carrying them") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, 10);
  addTemporalPolicies(manager, 10);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  ApplicablePolicyMasks decided_first;
  decided_first.symptom = {maskOf({1, 4}), maskOf({1, 4}), maskOf({1})};
  decided_first.temporal = {maskOf({2}), maskOf({2}), maskOf({2})};
  applicable_policies.importState(person, decided_first);
  ApplicablePolicyMasks undecided_after;
  undecided_after.symptom.applicable = maskOf({1, 4, 5});
  undecided_after.temporal.applicable = maskOf({2});

  applicable_policies.importState(person, undecided_after);

  CHECK(applicable_policies.exportState(person) == undecided_after);
}

namespace {

// The message of the exception an import throws, or empty.
std::string importRefusal(ApplicablePolicies& applicable_policies,
                          Person& person, const ApplicablePolicyMasks& masks) {
  return policy_test::thrown(
      [&] { applicable_policies.importState(person, masks); });
}

}  // namespace

TEST_CASE("an import naming a policy past those loaded is refused by name") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, 2);
  addTemporalPolicies(manager, 3);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  person.id = 42;
  ApplicablePolicyMasks loaded;
  loaded.symptom.applicable = maskOf({1});
  loaded.temporal.applicable = maskOf({2});
  CHECK(importRefusal(applicable_policies, person, loaded).empty());

  ApplicablePolicyMasks symptom_unloaded = loaded;
  symptom_unloaded.symptom.applicable.set(2);
  const std::string symptom_message =
      importRefusal(applicable_policies, person, symptom_unloaded);
  CHECK(has(symptom_message, "person 42"));
  CHECK(has(symptom_message, "symptom"));

  ApplicablePolicyMasks temporal_unloaded = loaded;
  temporal_unloaded.temporal.applicable.set(3);
  const std::string temporal_message =
      importRefusal(applicable_policies, person, temporal_unloaded);
  CHECK(has(temporal_message, "person 42"));
  CHECK(has(temporal_message, "temporal"));
}

TEST_CASE("an import over the cap of either kind is refused by name") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, kMaxPoliciesPerKind);
  addTemporalPolicies(manager, kMaxPoliciesPerKind);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  person.id = 7;
  PolicyMask over_cap;
  for (size_t index = 0; index <= kMaxApplicablePoliciesPerKind; ++index)
    over_cap.set(index);
  const std::string over_cap_count =
      std::to_string(kMaxApplicablePoliciesPerKind + 1);

  ApplicablePolicyMasks symptom_over;
  symptom_over.symptom.applicable = over_cap;
  const std::string symptom_message =
      importRefusal(applicable_policies, person, symptom_over);
  CHECK(has(symptom_message, "person 7"));
  CHECK(has(symptom_message, over_cap_count + " symptom"));

  ApplicablePolicyMasks temporal_over;
  temporal_over.temporal.applicable = over_cap;
  const std::string temporal_message =
      importRefusal(applicable_policies, person, temporal_over);
  CHECK(has(temporal_message, "person 7"));
  CHECK(has(temporal_message, over_cap_count + " temporal"));
}

TEST_CASE("an import deciding a policy outside the set is refused by name") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, 10);
  addTemporalPolicies(manager, 10);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  person.id = 9;
  ApplicablePolicyMasks within;
  within.symptom = {maskOf({1, 3}), maskOf({1}), maskOf({1})};
  within.temporal = {maskOf({2}), maskOf({2}), PolicyMask{}};
  CHECK(importRefusal(applicable_policies, person, within).empty());

  SUBCASE("decided") {
    ApplicablePolicyMasks outside = within;
    outside.temporal.decided.set(5);
    const std::string message =
        importRefusal(applicable_policies, person, outside);
    CHECK(has(message, "person 9"));
    CHECK(has(message, "temporal policy 5"));
  }
  SUBCASE("participating") {
    ApplicablePolicyMasks outside = within;
    outside.symptom.decided.set(4);
    outside.symptom.participating.set(4);
    const std::string message =
        importRefusal(applicable_policies, person, outside);
    CHECK(has(message, "person 9"));
    CHECK(has(message, "symptom policy 4"));
  }
}

TEST_CASE("an import participating in an undecided policy is refused by name") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, 10);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  person.id = 11;
  ApplicablePolicyMasks undecided_participation;
  undecided_participation.symptom = {maskOf({1, 3}), maskOf({1}),
                                     maskOf({1, 3})};

  const std::string message =
      importRefusal(applicable_policies, person, undecided_participation);

  CHECK(has(message, "person 11"));
  CHECK(has(message, "symptom policy 3"));
}

TEST_CASE("a refused import leaves the Person's state as it was") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, 10);
  addTemporalPolicies(manager, 10);
  auto& applicable_policies = manager.applicablePolicies();
  Person person;
  ApplicablePolicyMasks before;
  before.symptom = {maskOf({1, 3}), maskOf({1}), maskOf({1})};
  before.temporal = {maskOf({2}), maskOf({2}), PolicyMask{}};
  applicable_policies.importState(person, before);
  ApplicablePolicyMasks refused;
  refused.symptom = {maskOf({4}), maskOf({4}), PolicyMask{}};
  refused.temporal = {maskOf({5}), maskOf({6}), PolicyMask{}};

  CHECK_FALSE(importRefusal(applicable_policies, person, refused).empty());

  CHECK(applicable_policies.exportState(person) == before);
}

TEST_CASE("a refused import leaves the set table as it was") {
  ApplicablePolicyMasks unseen;
  unseen.symptom = {maskOf({4}), PolicyMask{}, PolicyMask{}};
  unseen.temporal = {maskOf({5}), PolicyMask{}, PolicyMask{}};

  WorldState refusing_world;
  PolicyManager refusing_manager(refusing_world);
  addSymptomPolicies(refusing_manager, 10);
  addTemporalPolicies(refusing_manager, 10);
  ApplicablePolicyMasks refused;
  refused.symptom = {maskOf({1}), PolicyMask{}, PolicyMask{}};
  refused.temporal = {maskOf({2}), maskOf({6}), PolicyMask{}};
  Person refused_person;
  CHECK_FALSE(importRefusal(refusing_manager.applicablePolicies(),
                            refused_person, refused)
                  .empty());
  Person after_refusal;
  refusing_manager.applicablePolicies().importState(after_refusal, unseen);

  WorldState fresh_world;
  PolicyManager fresh_manager(fresh_world);
  addSymptomPolicies(fresh_manager, 10);
  addTemporalPolicies(fresh_manager, 10);
  Person without_refusal;
  fresh_manager.applicablePolicies().importState(without_refusal, unseen);

  // Set ids are handed out in turn, so an orphan set shifts the next one.
  CHECK(after_refusal.applicable_policy_state ==
        without_refusal.applicable_policy_state);
}

TEST_CASE("walk rules read and write the state import and export carry") {
  WorldState world;
  PolicyManager manager(world);
  addSymptomPolicies(manager, 10);
  auto& applicable_policies = manager.applicablePolicies();
  SymptomPolicy complies;
  complies.action.compliance_rate = 1.0;
  const std::vector<SymptomPolicy> all_comply(10, complies);
  Person person;
  ApplicablePolicyMasks imported;
  imported.symptom.applicable = maskOf({0, 1, 2});

  SUBCASE("an imported refusal sticks through a latch") {
    imported.symptom.decided = maskOf({0});
    applicable_policies.importState(person, imported);
    CHECK_FALSE(applicable_policies.walk(person, all_comply).find(0)->latch());
    CHECK(applicable_policies.exportState(person) == imported);
  }
  SUBCASE("a follow-up inherits imported compliance") {
    imported.symptom.decided = maskOf({0});
    imported.symptom.participating = maskOf({0});
    applicable_policies.importState(person, imported);
    applicable_policies.walk(person, all_comply)
        .find(0)
        ->inheritTo(1, true, false);
    imported.symptom.decided.set(1);
    imported.symptom.participating.set(1);
    CHECK(applicable_policies.exportState(person) == imported);
  }
  SUBCASE("a follow-up inherits imported refusal") {
    imported.symptom.decided = maskOf({0});
    applicable_policies.importState(person, imported);
    applicable_policies.walk(person, all_comply)
        .find(0)
        ->inheritTo(1, false, true);
    imported.symptom.decided.set(1);
    CHECK(applicable_policies.exportState(person) == imported);
  }
  SUBCASE("compliance not inherited leaves the follow-up undecided") {
    imported.symptom.decided = maskOf({0});
    imported.symptom.participating = maskOf({0});
    applicable_policies.importState(person, imported);
    applicable_policies.walk(person, all_comply)
        .find(0)
        ->inheritTo(1, false, true);
    CHECK(applicable_policies.exportState(person) == imported);
  }
  SUBCASE("an undecided policy hands nothing on") {
    applicable_policies.importState(person, imported);
    applicable_policies.walk(person, all_comply)
        .find(0)
        ->inheritTo(1, true, true);
    CHECK(applicable_policies.exportState(person) == imported);
  }
  SUBCASE("a follow-up outside the set inherits nothing") {
    imported.symptom.decided = maskOf({0});
    imported.symptom.participating = maskOf({0});
    applicable_policies.importState(person, imported);
    applicable_policies.walk(person, all_comply)
        .find(0)
        ->inheritTo(5, true, true);
    CHECK(applicable_policies.exportState(person) == imported);
  }
  SUBCASE("a query walk agrees with a later latch and changes nothing") {
    SymptomPolicy halves;
    halves.action.compliance_rate = 0.5;
    const std::vector<SymptomPolicy> half_comply(10, halves);
    manager.setBaseSeed(2024);
    imported.symptom.decided = maskOf({1});
    applicable_policies.importState(person, imported);
    for (PersonId person_id = 0; person_id < 32; ++person_id) {
      person.id = person_id;
      const Person& query_person = person;
      std::vector<bool> queried;
      for (const auto& applicable :
           applicable_policies.walk(query_person, half_comply))
        queried.push_back(applicable.participates());
      CHECK(applicable_policies.exportState(person) == imported);
      std::vector<bool> latched;
      for (const auto& applicable :
           applicable_policies.walk(person, half_comply))
        latched.push_back(applicable.latch());
      CHECK(queried == latched);
      applicable_policies.importState(person, imported);
    }
  }
}
