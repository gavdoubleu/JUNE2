#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <memory>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/emission/emission.h"
#include "epidemiology/interaction_manager.h"
#include "test_utils.h"

using namespace june;

// =============================================================================
// Sibling mixing: susceptibles in one child venue (classroom) are exposed to
// infectious people in the other children of the same parent (school).
//
// World: school (venue 0) with classrooms A (venue 1) and B (venue 2).
// Person 0 sits in A; persons 1-5 sit in B. The only possible source for B
// is its sibling A, so any infection in B is a sibling infection.
// =============================================================================

namespace {

constexpr VenueId kSchoolId = 0;
constexpr VenueId kClassroomAId = 1;
constexpr VenueId kClassroomBId = 2;
constexpr int kNumPeople = 6;  // 0 in A, 1-5 in B
constexpr double kCurrentTime = 5.0;
constexpr double kDeltaHours = 8.0;

struct SiblingMixingWorld {
  WorldState world;
  std::unique_ptr<Disease> disease;
  ContactMatrixConfig contact_matrices;

  SiblingMixingWorld() {
    world.venue_type_names = {"school", "classroom"};
    world.geo_level_names = {"city"};
    auto add_venue = [&](VenueId id, uint8_t type_id, VenueId parent_id) {
      Venue venue;
      venue.id = id;
      venue.type_id = type_id;
      venue.geo_unit_id = -1;
      venue.parent_id = parent_id;
      world.venues.push_back(venue);
    };
    add_venue(kSchoolId, 0, -1);
    add_venue(kClassroomAId, 1, kSchoolId);
    add_venue(kClassroomBId, 1, kSchoolId);

    for (int i = 0; i < kNumPeople; ++i) {
      auto& person = world.people.emplace_back();
      person.id = i;
      person.age = 10.0f;
      person.sex = Sex::MALE;
      person.geo_unit_id = -1;
    }
    world.buildIndices();

    // "exposed" has no curve, so an exposed person is infected but emits
    // nothing.
    TransmissionParams transmission;
    transmission.mode = InfectiousnessMode::STAGE_DRIVEN;
    auto curve = std::make_shared<ConstantCurve>(1.0);
    transmission.stage_curves["infectious"] = curve;
    transmission.symptom_id_curves = {nullptr, nullptr, curve};

    TrajectoryDefinition trajectory;
    trajectory.selection_key = "general_population";
    trajectory.stages.push_back({"exposed", {"constant", {{"value", 1.0}}}});
    trajectory.stages.push_back(
        {"infectious", {"constant", {{"value", 100.0}}}});
    trajectory.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});

    SymptomTag healthy{.name = "healthy", .value = -1, .id = 0};
    SymptomTag exposed{.name = "exposed", .value = 0, .id = 1};
    SymptomTag infectious{.name = "infectious", .value = 1, .id = 2};
    DiseaseStageSettings stage_settings;
    disease = std::make_unique<Disease>(
        "Flu", std::vector<SymptomTag>{healthy, exposed, infectious},
        stage_settings, std::vector<TrajectoryDefinition>{trajectory},
        OutcomeRates{}, transmission);

    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{100.0}};
    contact_matrices.default_matrix = default_contact_matrix;
    contact_matrices.allow_default_matrix = true;
    finalizeContactMatrices(contact_matrices, world, *disease);
  }

  void infect(int person_index, double infection_time) {
    world.people[person_index].infection = std::make_unique<Infection>(
        disease.get(), infection_time, &world.people[person_index], 42, nullptr,
        "classroom", 0);
  }

  // Person 0 in classroom A, persons 1-5 in classroom B.
  std::vector<PersonLocation> locations() const {
    std::vector<PersonLocation> locations;
    locations.push_back({0, kClassroomAId, -1, 0, 255, 0});
    for (int i = 1; i < kNumPeople; ++i)
      locations.push_back(
          {i, kClassroomBId, -1, 0, 255, static_cast<size_t>(i)});
    return locations;
  }

  // Returns the number of new infections among persons 1-5 (classroom B).
  int runTickAndCountClassroomBInfections() {
    std::vector<bool> infected_before(kNumPeople);
    for (int i = 0; i < kNumPeople; ++i)
      infected_before[i] = world.people[i].infection != nullptr;

    SimulationConfig simulation_config;
    ParallelConfig parallel_config;
    InteractionManager interaction_manager(world, contact_matrices,
                                           simulation_config, parallel_config,
                                           disease.get(), nullptr);
    interaction_manager.processTransmissions(locations(), kCurrentTime,
                                             kDeltaHours, nullptr);

    int new_infections = 0;
    for (int i = 1; i < kNumPeople; ++i)
      if (!infected_before[i] && world.people[i].infection) new_infections++;
    return new_infections;
  }
};

}  // namespace

TEST_CASE("Sibling mixing: clean child is infected by an infectious sibling") {
  SiblingMixingWorld setup;
  setup.infect(0, 0.0);  // infectious by t=5

  CHECK(setup.runTickAndCountClassroomBInfections() > 0);
}

TEST_CASE(
    "Sibling mixing: child whose only case is not yet infectious is infected "
    "by an infectious sibling") {
  SiblingMixingWorld setup;
  setup.infect(0, 0.0);                 // infectious by t=5
  setup.infect(5, kCurrentTime - 0.1);  // still exposed at t=5

  CHECK(setup.runTickAndCountClassroomBInfections() > 0);
}

TEST_CASE("Sibling mixing: no infectious sibling means no infections") {
  SiblingMixingWorld setup;
  setup.infect(0, kCurrentTime - 0.1);  // still exposed at t=5

  CHECK(setup.runTickAndCountClassroomBInfections() == 0);
}

// Classroom bins take a local's Emission from EmissionCalculator::emit, so
// the school's aggregate must carry that same Emission for the infector.
TEST_CASE("Sibling mixing: parent aggregate carries the infector's Emission") {
  SiblingMixingWorld setup;
  setup.infect(0, 0.0);  // infectious by t=5

  SimulationConfig simulation_config;
  ParallelConfig parallel_config;
  InteractionManager interaction_manager(setup.world, setup.contact_matrices,
                                         simulation_config, parallel_config,
                                         setup.disease.get(), nullptr);
  interaction_manager.processTransmissions(setup.locations(), kCurrentTime,
                                           kDeltaHours, nullptr);

  const ParentAggregate* school =
      interaction_manager.getParentAggregate(kSchoolId);
  REQUIRE(school != nullptr);
  REQUIRE(school->infectors_by_bin.size() == 1);
  REQUIRE(school->infectors_by_bin[0].size() == 1);
  const ParentInfectorEntry& infector = school->infectors_by_bin[0][0];
  CHECK(infector.person_id == 0);

  Emission emission;
  EmissionCalculator(*setup.disease, kDeltaHours)
      .emit(setup.world.people[0], kCurrentTime, emission);
  REQUIRE_FALSE(emission.infectiousness_by_mode.empty());
  CHECK(infector.inf_by_mode == emission.infectiousness_by_mode);
}
