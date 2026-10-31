// A checkpoint restored under a changed policy list carries each person's
// decisions, participation and freeze over by policy name.

#include <H5Cpp.h>

#include <filesystem>
#include <initializer_list>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/policy.h"
#include "simulation/simulator.h"
#include "simulator_fixtures.h"
#include "test_utils.h"

using namespace june;
using namespace june::simulator_fixture;

using Names = std::set<std::string>;

namespace doctest {
template <>
struct StringMaker<Names> {
  static String convert(const Names& names) {
    std::string text = "{";
    for (const std::string& name : names)
      text += (text.size() > 1 ? ", " : "") + name;
    return (text + "}").c_str();
  }
};
}  // namespace doctest

namespace {

// The checkpoint is written after this day; the run never steps, so any day
// before the end will do.
constexpr int kCompletedDay = 2;
constexpr const char* kCompletedDate = "2020-01-03";

// Schedule ids of a frozen person's hop. The run never steps, so they need not
// name real schedules.
constexpr int16_t kTravelScheduleId = 1;
constexpr int16_t kReturnScheduleId = 0;
constexpr int16_t kReplacementScheduleId = 2;

// One policy entry of a policies.yaml. A symptom policy (triggered by mild)
// applies to people younger than `younger_than`, or to everyone if it is 0.
std::string symptomPolicy(const std::string& name, int younger_than = 0) {
  std::string entry = "    - name: \"" + name +
                      "\"\n"
                      "      symptoms: [\"mild\"]\n"
                      "      override_activities: [\"primary_activity\"]\n"
                      "      replacement: \"residence\"\n";
  if (younger_than > 0)
    entry +=
        "      applies_to:\n"
        "        - property: \"age\"\n"
        "          operator: \"<\"\n"
        "          value: " +
        std::to_string(younger_than) + "\n";
  return entry;
}

std::string temporalPolicy(const std::string& name) {
  return "    - name: \"" + name +
         "\"\n"
         "      override_activities: [\"primary_activity\"]\n"
         "      replacement: \"residence\"\n";
}

std::string policiesYaml(std::initializer_list<std::string> symptom_entries,
                         std::initializer_list<std::string> temporal_entries) {
  std::string text = "policies:\n";
  if (symptom_entries.size() > 0) {
    text += "  symptom_policies:\n";
    for (const std::string& entry : symptom_entries) text += entry;
  }
  if (temporal_entries.size() > 0) {
    text += "  temporal_policies:\n";
    for (const std::string& entry : temporal_entries) text += entry;
  }
  return text;
}

// The policy list checkpoints are written under. symptom_young applies to the
// people aged 20 to 29 (ids 0 to 9).
std::string checkpointedPolicies() {
  return policiesYaml(
      {symptomPolicy("symptom_all"), symptomPolicy("symptom_also"),
       symptomPolicy("symptom_young", 30)},
      {temporalPolicy("temporal_first"), temporalPolicy("temporal_second")});
}

// A Simulator over its own small world (people aged 20 + id).
// `policies_file` empty means no policies.
struct SmallWorldRun {
  WorldState world;
  Config config;
  std::unique_ptr<Simulator> simulator;

  SmallWorldRun(const ScopedTestFiles& files, const std::string& run_name,
                const std::string& policies_file) {
    buildSmallWorld(world);
    config = makeSmallWorldConfig(world, "2020-01-10");
    config.simulation.policies_file = policies_file;
    simulator = std::make_unique<Simulator>(
        world, config, nullptr, "tests/configs/simulator_infection_seeds.yaml",
        (files.directory() / (run_name + ".h5")).string());
  }

  PolicyManager& policies() { return simulator->getPolicyManager(); }

  // Names of the policies of `kind` set in `mask`.
  Names namesIn(const PolicyMask& mask, PolicyKind kind) {
    Names names;
    const PolicyManager& manager = policies();
    if (kind == PolicyKind::Symptom) {
      for (size_t i = 0; i < manager.getSymptomPolicyCount(); ++i)
        if (mask.test(i)) names.insert(manager.getSymptomPolicies()[i].name);
    } else {
      for (size_t i = 0; i < manager.getTemporalPolicyCount(); ++i)
        if (mask.test(i)) names.insert(manager.getTemporalPolicies()[i].name);
    }
    return names;
  }

  KindPolicyMasks masks(PersonId person_id, PolicyKind kind) {
    return policies()
        .applicablePolicies()
        .exportState(world.people[person_id])
        .of(kind);
  }
  Names decided(PersonId person_id, PolicyKind kind) {
    return namesIn(masks(person_id, kind).decided, kind);
  }
  Names participating(PersonId person_id, PolicyKind kind) {
    return namesIn(masks(person_id, kind).participating, kind);
  }

  size_t symptomIndex(const std::string& name) {
    const auto& symptom_policies = policies().getSymptomPolicies();
    for (size_t i = 0; i < symptom_policies.size(); ++i)
      if (symptom_policies[i].name == name) return i;
    FAIL("no symptom policy " << name);
    return 0;
  }
};

// Every person has decided on every policy that applies to them except
// symptom_all and temporal_first, the first of each kind: left undecided, so a
// stored mask lined up by position rather than matched by name shows. The even
// ids participate in what they decided, the odd ids refused it.
void decideApplicablePolicies(SmallWorldRun& run) {
  ApplicablePolicies& applicable_policies = run.policies().applicablePolicies();
  for (Person& person : run.world.people) {
    ApplicablePolicyMasks masks = applicable_policies.exportState(person);
    for (PolicyKind kind : {PolicyKind::Symptom, PolicyKind::Temporal}) {
      KindPolicyMasks& kind_masks = masks.of(kind);
      kind_masks.decided = kind_masks.applicable;
      kind_masks.decided.reset(0);
      kind_masks.participating =
          person.id % 2 == 0 ? kind_masks.decided : PolicyMask{};
    }
    applicable_policies.importState(person, masks);
  }
}

// `person_id` is mid-hop and frozen by the symptom policy `policy_name`.
void freezeBy(SmallWorldRun& run, PersonId person_id,
              const std::string& policy_name) {
  Person& person = run.world.people[person_id];
  person.schedule_hop.restoreTargets(kReplacementScheduleId, kReturnScheduleId);
  run.policies().setFrozenStates(
      {{person_id,
        FrozenPersonState{static_cast<uint8_t>(run.symptomIndex(policy_name)),
                          kTravelScheduleId, kReturnScheduleId, 1, -1}}});
}

// Writes a checkpoint under checkpointedPolicies() with the decisions of
// decideApplicablePolicies, person 0 frozen by symptom_young, then restores it
// under `resumed_policies` (a policies.yaml text, or "" for none).
std::unique_ptr<SmallWorldRun> resumeUnder(
    const ScopedTestFiles& files, const std::string& resumed_policies) {
  {
    SmallWorldRun checkpointed(
        files, "checkpointed",
        files.write("checkpointed_policies.yaml", checkpointedPolicies())
            .string());
    decideApplicablePolicies(checkpointed);
    freezeBy(checkpointed, 0, "symptom_young");
    checkpointed.simulator->writeCheckpoint(kCompletedDay, kCompletedDate);
  }
  const std::string resumed_file =
      resumed_policies.empty()
          ? ""
          : files.write("resumed_policies.yaml", resumed_policies).string();
  auto resumed =
      std::make_unique<SmallWorldRun>(files, "resumed", resumed_file);
  resumed->simulator->restoreFromCheckpoint(
      (files.directory() / "checkpoints" / "checkpoint_20200103_day002")
          .string());
  return resumed;
}

// Ids 0 and 1 (aged 20, 21) stand for everyone young enough for
// symptom_young, 12 and 13 (aged 32, 33) for everyone too old. Even ids
// participate in what they decided, odd ids refused it.
void checkSymptom(SmallWorldRun& run, PersonId person_id,
                  const Names& expected_decided) {
  CHECK(run.decided(person_id, PolicyKind::Symptom) == expected_decided);
  CHECK(run.participating(person_id, PolicyKind::Symptom) ==
        (person_id % 2 == 0 ? expected_decided : Names{}));
}
void checkTemporal(SmallWorldRun& run, PersonId person_id,
                   const Names& expected_decided) {
  CHECK(run.decided(person_id, PolicyKind::Temporal) == expected_decided);
  CHECK(run.participating(person_id, PolicyKind::Temporal) ==
        (person_id % 2 == 0 ? expected_decided : Names{}));
}

bool isFrozen(SmallWorldRun& run, PersonId person_id) {
  return run.policies().getFrozenStates().count(person_id) > 0;
}

// Person 0's freeze is kept, now pointing at symptom_young's current index,
// and the hop still targets the replacement schedule.
void checkFreezeKept(SmallWorldRun& run) {
  REQUIRE(isFrozen(run, 0));
  const FrozenPersonState& state = run.policies().getFrozenStates().at(0);
  CHECK(state.triggering_policy_index == run.symptomIndex("symptom_young"));
  CHECK(state.paused_hopped_schedule_id == kTravelScheduleId);
  CHECK(run.world.people[0].schedule_hop.hopped_schedule_id ==
        kReplacementScheduleId);
}

// Person 0's freeze is lifted and the hop put back on the travel schedule, as
// on recovery.
void checkFreezeReleased(SmallWorldRun& run) {
  CHECK_FALSE(isFrozen(run, 0));
  CHECK(run.world.people[0].schedule_hop.hopped_schedule_id ==
        kTravelScheduleId);
  CHECK(run.world.people[0].schedule_hop.return_schedule_id ==
        kReturnScheduleId);
}

}  // namespace

TEST_CASE("Checkpoint restore - unchanged policy list keeps every decision") {
  ScopedTestFiles files{"june_test_policy_restore"};
  auto run = resumeUnder(files, checkpointedPolicies());

  for (PersonId person_id : {0, 1}) {
    checkSymptom(*run, person_id, {"symptom_also", "symptom_young"});
    checkTemporal(*run, person_id, {"temporal_second"});
  }
  for (PersonId person_id : {12, 13}) {
    checkSymptom(*run, person_id, {"symptom_also"});
    checkTemporal(*run, person_id, {"temporal_second"});
  }
  CHECK(run->world.people[1].schedule_hop.hopped_schedule_id == -1);
  checkFreezeKept(*run);
}

TEST_CASE("Checkpoint restore - fewer policies drop the gone ones") {
  ScopedTestFiles files{"june_test_policy_restore"};

  SUBCASE("last of each kind removed") {
    auto run =
        resumeUnder(files, policiesYaml({symptomPolicy("symptom_all"),
                                         symptomPolicy("symptom_also")},
                                        {temporalPolicy("temporal_first")}));
    for (PersonId person_id : {0, 1, 12, 13}) {
      checkSymptom(*run, person_id, {"symptom_also"});
      checkTemporal(*run, person_id, {});
    }
    checkFreezeReleased(*run);
  }

  SUBCASE("middle and first removed") {
    auto run =
        resumeUnder(files, policiesYaml({symptomPolicy("symptom_all"),
                                         symptomPolicy("symptom_young", 30)},
                                        {temporalPolicy("temporal_second")}));
    for (PersonId person_id : {0, 1}) {
      checkSymptom(*run, person_id, {"symptom_young"});
      checkTemporal(*run, person_id, {"temporal_second"});
    }
    for (PersonId person_id : {12, 13}) {
      checkSymptom(*run, person_id, {});
      checkTemporal(*run, person_id, {"temporal_second"});
    }
    checkFreezeKept(*run);
  }

  SUBCASE("none at all") {
    auto run = resumeUnder(files, "");
    for (PersonId person_id : {0, 1, 12, 13}) {
      checkSymptom(*run, person_id, {});
      checkTemporal(*run, person_id, {});
    }
    checkFreezeReleased(*run);
  }
}

TEST_CASE("Checkpoint restore - added policies start undecided") {
  ScopedTestFiles files{"june_test_policy_restore"};
  auto run = resumeUnder(
      files,
      policiesYaml(
          {symptomPolicy("symptom_new"), symptomPolicy("symptom_all"),
           symptomPolicy("symptom_also"), symptomPolicy("symptom_young", 30)},
          {temporalPolicy("temporal_first"), temporalPolicy("temporal_new"),
           temporalPolicy("temporal_second")}));

  for (PersonId person_id : {0, 1}) {
    checkSymptom(*run, person_id, {"symptom_also", "symptom_young"});
    checkTemporal(*run, person_id, {"temporal_second"});
  }
  for (PersonId person_id : {12, 13}) {
    checkSymptom(*run, person_id, {"symptom_also"});
    checkTemporal(*run, person_id, {"temporal_second"});
  }
  // The new policies apply, but nobody has decided on them yet.
  CHECK(run->namesIn(run->masks(0, PolicyKind::Symptom).applicable,
                     PolicyKind::Symptom)
            .count("symptom_new") == 1);
  CHECK(run->namesIn(run->masks(0, PolicyKind::Temporal).applicable,
                     PolicyKind::Temporal)
            .count("temporal_new") == 1);
  checkFreezeKept(*run);
}

TEST_CASE("Checkpoint restore - reordered policies keep decisions by name") {
  ScopedTestFiles files{"june_test_policy_restore"};
  auto run =
      resumeUnder(files, policiesYaml({symptomPolicy("symptom_young", 30),
                                       symptomPolicy("symptom_also"),
                                       symptomPolicy("symptom_all")},
                                      {temporalPolicy("temporal_second"),
                                       temporalPolicy("temporal_first")}));

  for (PersonId person_id : {0, 1}) {
    checkSymptom(*run, person_id, {"symptom_also", "symptom_young"});
    checkTemporal(*run, person_id, {"temporal_second"});
  }
  for (PersonId person_id : {12, 13}) {
    checkSymptom(*run, person_id, {"symptom_also"});
    checkTemporal(*run, person_id, {"temporal_second"});
  }
  REQUIRE(run->symptomIndex("symptom_young") == 0);
  checkFreezeKept(*run);
}

TEST_CASE("Checkpoint restore - edited criteria follow the current config") {
  ScopedTestFiles files{"june_test_policy_restore"};

  SUBCASE("narrowed: people no longer applicable lose the decision") {
    // symptom_young now applies to ids 0 to 4 only.
    auto run =
        resumeUnder(files, policiesYaml({symptomPolicy("symptom_all"),
                                         symptomPolicy("symptom_also"),
                                         symptomPolicy("symptom_young", 25)},
                                        {temporalPolicy("temporal_first"),
                                         temporalPolicy("temporal_second")}));
    for (PersonId person_id : {0, 1})
      checkSymptom(*run, person_id, {"symptom_also", "symptom_young"});
    for (PersonId person_id : {6, 7, 12})
      checkSymptom(*run, person_id, {"symptom_also"});
    checkFreezeKept(*run);
  }

  SUBCASE("narrowed past the frozen person: their freeze is released") {
    // symptom_young now applies to nobody aged 20 or over.
    auto run =
        resumeUnder(files, policiesYaml({symptomPolicy("symptom_all"),
                                         symptomPolicy("symptom_also"),
                                         symptomPolicy("symptom_young", 20)},
                                        {temporalPolicy("temporal_first"),
                                         temporalPolicy("temporal_second")}));
    checkSymptom(*run, 0, {"symptom_also"});
    checkFreezeReleased(*run);
  }

  SUBCASE("widened: newly applicable people start undecided") {
    // symptom_young now applies to ids 0 to 14.
    auto run =
        resumeUnder(files, policiesYaml({symptomPolicy("symptom_all"),
                                         symptomPolicy("symptom_also"),
                                         symptomPolicy("symptom_young", 35)},
                                        {temporalPolicy("temporal_first"),
                                         temporalPolicy("temporal_second")}));
    checkSymptom(*run, 2, {"symptom_also", "symptom_young"});
    checkSymptom(*run, 12, {"symptom_also"});
    CHECK(run->namesIn(run->masks(12, PolicyKind::Symptom).applicable,
                       PolicyKind::Symptom)
              .count("symptom_young") == 1);
    checkFreezeKept(*run);
  }
}

TEST_CASE("Checkpoint restore - a checkpoint without its policy roster is "
          "refused") {
  ScopedTestFiles files{"june_test_policy_restore"};
  const auto checkpoint =
      files.directory() / "checkpoints" / "checkpoint_20200103_day002";
  {
    SmallWorldRun checkpointed(
        files, "checkpointed",
        files.write("checkpointed_policies.yaml", checkpointedPolicies())
            .string());
    decideApplicablePolicies(checkpointed);
    checkpointed.simulator->writeCheckpoint(kCompletedDay, kCompletedDate);
  }
  {
    H5::H5File state((checkpoint / "state.h5").string(), H5F_ACC_RDWR);
    REQUIRE(H5Ldelete(state.getId(), "/policy_roster/temporal_names",
                      H5P_DEFAULT) >= 0);
  }
  SmallWorldRun resumed(
      files, "resumed",
      files.write("resumed_policies.yaml", checkpointedPolicies()).string());
  CHECK_THROWS_WITH_AS(
      resumed.simulator->restoreFromCheckpoint(checkpoint.string()),
      doctest::Contains("/policy_roster/temporal_names"), std::runtime_error);
}
