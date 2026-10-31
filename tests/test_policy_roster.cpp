#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

#include <optional>
#include <string>
#include <vector>

#include "epidemiology/policy_roster/policy_roster.h"
#include "policy_test_helpers.h"

using namespace june;

using june::policy_test::has;
using june::policy_test::maskOf;
using june::policy_test::thrown;

namespace {

PolicyRoster rosterOf(std::vector<std::string> symptom_names,
                      std::vector<std::string> temporal_names) {
  PolicyRoster roster;
  roster.symptom = std::move(symptom_names);
  roster.temporal = std::move(temporal_names);
  return roster;
}

using Indices = std::vector<std::optional<size_t>>;

}  // namespace

TEST_CASE("an identical roster maps every policy to itself, silently") {
  const PolicyRoster roster = rosterOf({"isolate", "test"}, {"lockdown"});
  const PolicyRosterMatch roster_match = match(roster, roster);

  CHECK(roster_match.identical());
  CHECK(roster_match.of(PolicyKind::Symptom).current_index == Indices{0, 1});
  CHECK(roster_match.of(PolicyKind::Temporal).current_index == Indices{0});
  CHECK(roster_match.remap(maskOf({0, 1}), PolicyKind::Symptom) ==
        maskOf({0, 1}));
  CHECK(roster_match.describeChanges().empty());
}

TEST_CASE("a policy absent from the current roster is dropped") {
  const PolicyRosterMatch roster_match =
      match(rosterOf({}, {"lockdown", "curfew"}), rosterOf({}, {"lockdown"}));
  const KindRosterMatch& temporal = roster_match.of(PolicyKind::Temporal);

  CHECK_FALSE(roster_match.identical());
  CHECK(temporal.current_index == Indices{0, std::nullopt});
  CHECK(temporal.dropped == std::vector<std::string>{"curfew"});
  CHECK(temporal.added.empty());
  CHECK(temporal.moved.empty());
  CHECK(roster_match.remap(maskOf({0, 1}), PolicyKind::Temporal) ==
        maskOf({0}));
  CHECK(has(roster_match.describeChanges(), "dropped 'curfew'"));
}

TEST_CASE("a policy new to the current roster is added, mapped from nothing") {
  const PolicyRosterMatch roster_match =
      match(rosterOf({}, {"lockdown"}), rosterOf({}, {"lockdown", "curfew"}));
  const KindRosterMatch& temporal = roster_match.of(PolicyKind::Temporal);

  CHECK_FALSE(roster_match.identical());
  CHECK(temporal.current_index == Indices{0});
  CHECK(temporal.added == std::vector<std::string>{"curfew"});
  CHECK(temporal.dropped.empty());
  CHECK(roster_match.remap(maskOf({0}), PolicyKind::Temporal) == maskOf({0}));
  CHECK(has(roster_match.describeChanges(), "added 'curfew'"));
}

TEST_CASE("a reordered roster moves decisions to each policy's new index") {
  const PolicyRosterMatch roster_match =
      match(rosterOf({"a", "b", "c"}, {}), rosterOf({"c", "a", "b"}, {}));
  const KindRosterMatch& symptom = roster_match.of(PolicyKind::Symptom);

  CHECK_FALSE(roster_match.identical());
  CHECK(symptom.current_index == Indices{1, 2, 0});
  CHECK(symptom.dropped.empty());
  CHECK(symptom.added.empty());
  REQUIRE(symptom.moved.size() == 3);
  CHECK(symptom.moved[0].name == "a");
  CHECK(symptom.moved[0].stored_index == 0);
  CHECK(symptom.moved[0].current_index == 1);
  CHECK(roster_match.remap(maskOf({0, 2}), PolicyKind::Symptom) ==
        maskOf({0, 1}));
  CHECK(has(roster_match.describeChanges(),
            "moved 'a' 0->1, 'b' 1->2, 'c' 2->0"));
}

TEST_CASE("removing a middle policy shifts the later ones down") {
  const PolicyRosterMatch roster_match =
      match(rosterOf({"a", "b", "c", "d"}, {}), rosterOf({"a", "c", "d"}, {}));
  const KindRosterMatch& symptom = roster_match.of(PolicyKind::Symptom);

  CHECK(symptom.current_index == Indices{0, std::nullopt, 1, 2});
  CHECK(symptom.dropped == std::vector<std::string>{"b"});
  REQUIRE(symptom.moved.size() == 2);
  CHECK(symptom.moved[0].name == "c");
  CHECK(symptom.moved[1].name == "d");
  CHECK(roster_match.remap(maskOf({0, 1, 3}), PolicyKind::Symptom) ==
        maskOf({0, 2}));
  CHECK(roster_match.currentIndex(PolicyKind::Symptom, 3) == 2);
  CHECK(roster_match.currentIndex(PolicyKind::Symptom, 1) == std::nullopt);
}

TEST_CASE("empty rosters match, and emptying a roster drops everything") {
  const PolicyRoster empty;
  CHECK(match(empty, empty).identical());
  CHECK(match(empty, empty).describeChanges().empty());

  const PolicyRosterMatch emptied = match(rosterOf({"a"}, {"b"}), empty);
  CHECK(emptied.of(PolicyKind::Symptom).dropped ==
        std::vector<std::string>{"a"});
  CHECK(emptied.of(PolicyKind::Temporal).dropped ==
        std::vector<std::string>{"b"});
  CHECK(emptied.remap(maskOf({0}), PolicyKind::Symptom).none());

  const PolicyRosterMatch filled = match(empty, rosterOf({"a"}, {}));
  CHECK(filled.of(PolicyKind::Symptom).added == std::vector<std::string>{"a"});
  CHECK(filled.remap(PolicyMask{}, PolicyKind::Symptom).none());
}

TEST_CASE("each kind is matched on its own names") {
  // The same name in both kinds names two different policies.
  const PolicyRosterMatch roster_match =
      match(rosterOf({"shared", "isolate"}, {"shared"}),
            rosterOf({"isolate", "shared"}, {"curfew", "shared"}));

  CHECK(roster_match.of(PolicyKind::Symptom).current_index == Indices{1, 0});
  CHECK(roster_match.of(PolicyKind::Symptom).added.empty());
  CHECK(roster_match.of(PolicyKind::Temporal).current_index == Indices{1});
  CHECK(roster_match.of(PolicyKind::Temporal).added ==
        std::vector<std::string>{"curfew"});
  CHECK(roster_match.remap(maskOf({0}), PolicyKind::Symptom) == maskOf({1}));
  CHECK(roster_match.remap(maskOf({0}), PolicyKind::Temporal) == maskOf({1}));

  const std::string changes = roster_match.describeChanges();
  CHECK(has(changes, "symptom policies"));
  CHECK(has(changes, "temporal policies"));
}

TEST_CASE("only the changed kind is described") {
  const std::string changes =
      match(rosterOf({"a"}, {"b"}), rosterOf({"a"}, {})).describeChanges();
  CHECK_FALSE(has(changes, "symptom"));
  CHECK(has(changes, "temporal policies: dropped 'b'"));
}

TEST_CASE("remap refuses a stored policy past the stored roster") {
  const PolicyRosterMatch roster_match =
      match(rosterOf({"a"}, {}), rosterOf({"a"}, {}));
  CHECK(
      has(thrown([&] { roster_match.remap(maskOf({1}), PolicyKind::Symptom); }),
          "symptom policy 1"));
  CHECK(has(thrown([&] { roster_match.currentIndex(PolicyKind::Symptom, 1); }),
            "symptom policy 1"));
}

TEST_CASE("a name repeated within a kind is refused, naming kind and name") {
  CHECK_NOTHROW(rejectDuplicateNames(rosterOf({"same"}, {"same"})));

  const std::string message =
      thrown([] { rejectDuplicateNames(rosterOf({}, {"a", "b", "a"})); });
  CHECK(has(message, "temporal"));
  CHECK(has(message, "'a'"));
  CHECK(has(message, "0 and 2"));

  CHECK(
      has(thrown([] { match(rosterOf({"a", "a"}, {}), rosterOf({"a"}, {})); }),
          "symptom"));
}

TEST_CASE("a roster is built from each kind's policy names in order") {
  struct NamedPolicy {
    std::string name;
  };
  const PolicyRoster roster =
      makePolicyRoster(std::vector<NamedPolicy>{{"isolate"}, {"test"}},
                       std::vector<NamedPolicy>{{"lockdown"}});
  CHECK(roster.symptom == std::vector<std::string>{"isolate", "test"});
  CHECK(roster.temporal == std::vector<std::string>{"lockdown"});
}
