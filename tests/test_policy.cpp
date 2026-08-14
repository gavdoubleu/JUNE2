#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "epidemiology/policy.h"

#include "doctest.h"
#include "epidemiology/disease.h"
#include "test_utils.h"

using namespace june;

// =============================================================================
// Scheduled-venue filter — both directions (override_venue_types,
// exempt_venue_types).
//
// A policy action targets Activities only, and one Activity reaches many Venue
// types (leisure → pub, grocery, care_home...). The gate narrows an override to
// a named set of venue types, ANDed with the activity mask, so pubs can close
// while groceries stay open.
// =============================================================================

namespace {

constexpr int16_t kResidence = 0;
constexpr int16_t kLeisure = 1;
constexpr int16_t kPrimaryActivity = 2;

constexpr VenueId kHouseholdVenue = 0;
constexpr VenueId kPubVenue = 1;
constexpr VenueId kGroceryVenue = 2;
constexpr VenueId kCareHomeVenue = 3;
constexpr VenueId kClassroomVenue = 4;
constexpr VenueId kOfficeVenue = 5;

constexpr int kVenueTypeCount = 6;

// venue types: household=0, pub=1, grocery=2, care_home=3, classroom=4,
//              office=5
// activities:  residence=0, leisure=1, primary_activity=2
// Every person's residence is the household venue, so an override that fires
// redirects them there.
WorldState buildVenueGateWorld(int num_people = 1) {
  WorldState world =
      TestWorldFactory::createMinimalWorld(num_people, kVenueTypeCount);
  world.activity_names = {"residence", "leisure", "primary_activity"};
  world.venue_type_names = {"household",  "pub",       "grocery",
                            "care_home",  "classroom", "office"};
  for (int i = 0; i < kVenueTypeCount; ++i) {
    world.venues[i].type_id = static_cast<uint8_t>(i);
  }
  world.buildIndices();

  for (auto& person : world.people) {
    person.activity_meta_start =
        static_cast<uint32_t>(world.activity_meta.size());
    person.activity_meta_count = 1;
    world.activity_meta.push_back(
        {kResidence, static_cast<uint32_t>(world.activity_venues.size()), 1});
    world.activity_venues.push_back({kHouseholdVenue, 0});
    // Normally set by precomputePolicyApplicability.
    person.applicable_temporal_policy_mask = 1;
  }
  return world;
}

TemporalPolicy makeClosurePolicy(
    const std::vector<std::string>& activities,
    const std::vector<std::string>& venue_types,
    VenueGateDirection direction = VenueGateDirection::RestrictTo) {
  TemporalPolicy policy;
  policy.name = "closure";
  policy.start_time = 0.0;
  policy.end_time = 10.0;
  policy.action.override_activities.insert(activities.begin(),
                                           activities.end());
  auto& gated_types = (direction == VenueGateDirection::ExemptFrom)
                          ? policy.action.exempt_venue_types
                          : policy.action.override_venue_types;
  gated_types.insert(venue_types.begin(), venue_types.end());
  policy.action.replacement_activity = "residence";
  policy.action.compliance_rate = 1.0;
  return policy;
}

std::optional<PersonLocation> runOverride(PolicyManager& policy_manager,
                                          Person& person,
                                          int16_t activity_index,
                                          VenueId venue_id) {
  return policy_manager.getOverride(person, activity_index, venue_id, 0, 5.0,
                                    0);
}

}  // namespace

TEST_CASE("Venue gate - restrict-to fires only at listed venue types") {
  WorldState world = buildVenueGateWorld();
  PolicyManager policy_manager(world);

  TemporalPolicy close_pubs = makeClosurePolicy({"leisure"}, {"pub"});
  close_pubs.resolve(world);
  policy_manager.addTemporalPolicy(close_pubs);

  Person& person = world.people[0];

  SUBCASE("fires at a listed venue type") {
    auto override_location =
        runOverride(policy_manager, person, kLeisure, kPubVenue);
    REQUIRE(override_location.has_value());
    CHECK(override_location->venue_id == kHouseholdVenue);
    CHECK(override_location->activity_index == kResidence);
  }

  SUBCASE("does not fire at an unlisted venue type, same activity") {
    CHECK_FALSE(
        runOverride(policy_manager, person, kLeisure, kGroceryVenue).has_value());
  }

  SUBCASE("absent venue (venue_id < 0) does not fire") {
    CHECK_FALSE(
        runOverride(policy_manager, person, kLeisure, kInvalidVenueId)
            .has_value());
  }

  SUBCASE("unknown venue id reads as absent and does not fire") {
    CHECK_FALSE(
        runOverride(policy_manager, person, kLeisure, 9999).has_value());
  }
}

TEST_CASE("Venue gate - ANDs with the activity mask") {
  WorldState world = buildVenueGateWorld();
  PolicyManager policy_manager(world);

  SUBCASE("named activity: venue matches but activity does not") {
    TemporalPolicy close_pubs = makeClosurePolicy({"leisure"}, {"pub"});
    close_pubs.resolve(world);
    policy_manager.addTemporalPolicy(close_pubs);

    CHECK_FALSE(runOverride(policy_manager, world.people[0], kPrimaryActivity,
                            kPubVenue)
                    .has_value());
  }

  SUBCASE("wildcard activity: gate still restricts by venue type") {
    TemporalPolicy close_pubs = makeClosurePolicy({"*"}, {"pub"});
    close_pubs.resolve(world);
    policy_manager.addTemporalPolicy(close_pubs);

    Person& person = world.people[0];
    CHECK(runOverride(policy_manager, person, kPrimaryActivity, kPubVenue)
              .has_value());
    CHECK_FALSE(runOverride(policy_manager, person, kLeisure, kGroceryVenue)
                    .has_value());
  }
}

TEST_CASE("Venue gate - multi-role venue type separates the three roles") {
  // care_home is reached by leisure (visitors), primary_activity (staff) and
  // residence (residents). Banning visits must not evict residents or stop
  // staff working.
  WorldState world = buildVenueGateWorld(3);
  PolicyManager policy_manager(world);

  TemporalPolicy ban_visits = makeClosurePolicy({"leisure"}, {"care_home"});
  ban_visits.resolve(world);
  policy_manager.addTemporalPolicy(ban_visits);

  Person& visitor = world.people[0];
  Person& worker = world.people[1];
  Person& resident = world.people[2];

  CHECK(runOverride(policy_manager, visitor, kLeisure, kCareHomeVenue)
            .has_value());
  CHECK_FALSE(runOverride(policy_manager, worker, kPrimaryActivity,
                          kCareHomeVenue)
                  .has_value());
  CHECK_FALSE(
      runOverride(policy_manager, resident, kResidence, kCareHomeVenue)
          .has_value());
}

TEST_CASE("Venue gate - absent means unchanged behaviour") {
  WorldState world = buildVenueGateWorld();
  PolicyManager policy_manager(world);

  TemporalPolicy close_leisure = makeClosurePolicy({"leisure"}, {});
  close_leisure.resolve(world);
  policy_manager.addTemporalPolicy(close_leisure);

  Person& person = world.people[0];
  CHECK(runOverride(policy_manager, person, kLeisure, kPubVenue).has_value());
  CHECK(runOverride(policy_manager, person, kLeisure, kGroceryVenue)
            .has_value());
  CHECK(runOverride(policy_manager, person, kLeisure, kInvalidVenueId)
            .has_value());
}

TEST_CASE("Venue gate - validation throws from resolve") {
  WorldState world = buildVenueGateWorld();

  SUBCASE("unknown venue type name") {
    TemporalPolicy policy = makeClosurePolicy({"leisure"}, {"tavern"});
    policy.name = "close_taverns";
    CHECK_THROWS_AS(policy.resolve(world), std::runtime_error);
  }

  SUBCASE("venue type index at or beyond the 64-bit mask width") {
    world.venue_type_names.resize(70);
    for (size_t i = 4; i < world.venue_type_names.size(); ++i) {
      world.venue_type_names[i] = "filler_" + std::to_string(i);
    }
    TemporalPolicy policy = makeClosurePolicy({"leisure"}, {"filler_64"});
    CHECK_THROWS_AS(policy.resolve(world), std::runtime_error);
  }
}

// =============================================================================
// Exempt-from direction
// =============================================================================

TEST_CASE("Venue gate - exempt-from is the inverse of restrict-to") {
  WorldState world = buildVenueGateWorld();
  PolicyManager policy_manager(world);

  // Close all leisure except the pub.
  TemporalPolicy close_all_but_pubs = makeClosurePolicy(
      {"leisure"}, {"pub"}, VenueGateDirection::ExemptFrom);
  close_all_but_pubs.resolve(world);
  policy_manager.addTemporalPolicy(close_all_but_pubs);

  Person& person = world.people[0];

  SUBCASE("does not fire at a listed venue type") {
    CHECK_FALSE(
        runOverride(policy_manager, person, kLeisure, kPubVenue).has_value());
  }

  SUBCASE("fires at an unlisted venue type, same activity") {
    auto override_location =
        runOverride(policy_manager, person, kLeisure, kGroceryVenue);
    REQUIRE(override_location.has_value());
    CHECK(override_location->venue_id == kHouseholdVenue);
    CHECK(override_location->activity_index == kResidence);
  }

  SUBCASE("absent venue (venue_id < 0) fires: absent is never listed") {
    CHECK(runOverride(policy_manager, person, kLeisure, kInvalidVenueId)
              .has_value());
  }

  SUBCASE("unknown venue id reads as absent and fires") {
    CHECK(runOverride(policy_manager, person, kLeisure, 9999).has_value());
  }

  SUBCASE("still ANDs with the activity mask") {
    CHECK_FALSE(runOverride(policy_manager, person, kPrimaryActivity,
                            kGroceryVenue)
                    .has_value());
  }
}

TEST_CASE("Venue gate - both directions on one action throws from resolve") {
  WorldState world = buildVenueGateWorld();

  TemporalPolicy policy = makeClosurePolicy({"leisure"}, {"pub"});
  policy.name = "contradictory_closure";
  policy.action.exempt_venue_types = {"grocery"};
  CHECK_THROWS_AS(policy.resolve(world), std::runtime_error);
}

// =============================================================================
// Composition with the conjuncts either side of the gate
// =============================================================================

TEST_CASE("Venue gate - compliance latch is unaffected by the gate") {
  // The gate sits after the compliance latch, so a slot the gate blocks must
  // still leave the person's compliance decision exactly as an ungated policy
  // would, and that decision must stick across slots.
  constexpr int kNumPeople = 16;

  WorldState gated_world = buildVenueGateWorld(kNumPeople);
  PolicyManager gated_manager(gated_world);
  TemporalPolicy close_pubs = makeClosurePolicy({"leisure"}, {"pub"});
  close_pubs.action.compliance_rate = 0.5;
  close_pubs.resolve(gated_world);
  gated_manager.addTemporalPolicy(close_pubs);

  WorldState ungated_world = buildVenueGateWorld(kNumPeople);
  PolicyManager ungated_manager(ungated_world);
  TemporalPolicy close_leisure = makeClosurePolicy({"leisure"}, {});
  close_leisure.action.compliance_rate = 0.5;
  close_leisure.resolve(ungated_world);
  ungated_manager.addTemporalPolicy(close_leisure);

  int compliers = 0;
  for (int i = 0; i < kNumPeople; ++i) {
    Person& gated_person = gated_world.people[i];

    // First slot at the gated type reveals the compliance decision.
    const bool complies =
        runOverride(gated_manager, gated_person, kLeisure, kPubVenue)
            .has_value();
    compliers += complies ? 1 : 0;

    // An intervening slot the gate blocks must not disturb the latch.
    CHECK_FALSE(runOverride(gated_manager, gated_person, kLeisure, kGroceryVenue)
                    .has_value());
    CHECK(runOverride(gated_manager, gated_person, kLeisure, kPubVenue)
              .has_value() == complies);

    // Same decision as the identical policy without a gate.
    CHECK(runOverride(ungated_manager, ungated_world.people[i], kLeisure,
                      kGroceryVenue)
              .has_value() == complies);
  }

  // Test has power only if both outcomes occur.
  CHECK(compliers > 0);
  CHECK(compliers < kNumPeople);
}

TEST_CASE("Venue gate - gate and exemption are independent conjuncts") {
  // Close schools (primary_activity at classroom), key workers exempt.
  // Key worker is proxied by age: people are aged 20 + index.
  WorldState world = buildVenueGateWorld(2);
  PolicyManager policy_manager(world);

  TemporalPolicy close_schools =
      makeClosurePolicy({"primary_activity"}, {"classroom"});
  close_schools.name = "close_schools";
  ActivityExemption key_worker_exemption;
  key_worker_exemption.activity_name = "primary_activity";
  SelectionCriterion age_criterion;
  age_criterion.property_path = "age";
  age_criterion.operator_type = ">";
  age_criterion.value = 20.5;
  key_worker_exemption.criteria.push_back(age_criterion);
  close_schools.action.exemptions.push_back(key_worker_exemption);
  close_schools.resolve(world);
  policy_manager.addTemporalPolicy(close_schools);

  Person& non_key_worker = world.people[0];  // age 20
  Person& key_worker = world.people[1];      // age 21

  // Blocked by the exemption, at a gated venue type.
  CHECK_FALSE(runOverride(policy_manager, key_worker, kPrimaryActivity,
                          kClassroomVenue)
                  .has_value());
  // Blocked by the gate, exemption does not apply.
  CHECK_FALSE(runOverride(policy_manager, non_key_worker, kPrimaryActivity,
                          kOfficeVenue)
                  .has_value());
  // Neither blocks: overridden.
  CHECK(runOverride(policy_manager, non_key_worker, kPrimaryActivity,
                    kClassroomVenue)
            .has_value());
}

TEST_CASE("Venue gate - symptom path is gated too") {
  WorldState world = buildVenueGateWorld();

  TransmissionParams transmission;
  transmission.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(1.0);
  transmission.stage_curves["sick"] = curve;
  transmission.symptom_id_curves = {nullptr, curve};
  std::vector<SymptomTag> symptom_tags = {{"healthy", -1, 0}, {"sick", 1, 1}};
  DiseaseStageSettings stage_settings;
  stage_settings.recovered_stages = {"healthy"};
  TrajectoryDefinition trajectory;
  trajectory.selection_key = "general";
  trajectory.severity = 1.0;
  trajectory.stages.push_back({"sick", {"constant", {{"value", 100.0}}}});
  Disease disease("TestDisease", symptom_tags, stage_settings, {trajectory}, {},
                  transmission);

  PolicyManager policy_manager(world);
  SymptomPolicy stay_out_of_pubs;
  stay_out_of_pubs.name = "sick_avoid_pubs";
  stay_out_of_pubs.trigger_symptoms = {"sick"};
  stay_out_of_pubs.action.override_activities = {"leisure"};
  stay_out_of_pubs.action.override_venue_types = {"pub"};
  stay_out_of_pubs.action.replacement_activity = "residence";
  stay_out_of_pubs.action.compliance_rate = 1.0;
  policy_manager.addSymptomPolicy(stay_out_of_pubs);
  policy_manager.resolveAll(disease);

  Person& person = world.people[0];
  person.infection = std::make_unique<Infection>(&disease, 0.0, &person, 42,
                                                 nullptr, "household", 0);
  person.applicable_symptom_policy_mask = 1;

  CHECK(runOverride(policy_manager, person, kLeisure, kPubVenue).has_value());
  CHECK_FALSE(runOverride(policy_manager, person, kLeisure, kGroceryVenue)
                  .has_value());
}
