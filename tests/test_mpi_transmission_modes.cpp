// MPI tests for cross-rank transmission mode correctness.
//
// Must be run with exactly 2 MPI ranks:
//   mpirun -np 2 ./test_mpi_transmission_modes
//
// Tests cover stage-driven and trajectory-driven infectiousness propagation
// across rank boundaries, pending infection routing, multi-mode dispatch,
// immunity, bidirectional exchange, visitor/local fomite deposit parity, and
// mixed-state exchanges where records carry different tails.

#define DOCTEST_CONFIG_IMPLEMENT  // custom main so we can wrap MPI
                                  // init/finalize
#include "doctest.h"

#ifdef USE_MPI
#include <mpi.h>

#include <algorithm>
#include <cmath>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/disease.h"
#include "epidemiology/infectiousness_curves.h"
#include "epidemiology/interaction_manager.h"
#include "parallel/domain.h"
#include "parallel/domain_manager.h"
#include "test_utils.h"

using namespace june;

// ---------------------------------------------------------------------------
// TwoRankFixture: minimal 2-rank world (same pattern as
// test_mpi_communication.cpp)
//
// rank 0 owns: geo_unit 0, venue 0, person 0
// rank 1 owns: geo_unit 1, venue 1, person 1
// ---------------------------------------------------------------------------
struct TwoRankFixture {
  int rank;
  int size;

  WorldState world;
  Config config;
  std::unique_ptr<DomainManager> dm;

  static constexpr PersonId PERSON_R0 = 0;
  static constexpr PersonId PERSON_R1 = 1;
  static constexpr VenueId VENUE_R0 = 0;
  static constexpr VenueId VENUE_R1 = 1;

  TwoRankFixture() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    world.venue_type_names = {"household"};
    world.geo_level_names = {"MGU"};
    world.activity_names = {"residence", "work", "visiting", "none", "dead"};
    world.encounter_type_names = {};
    world.subset_type_names = {};

    GeographicalUnit gu;
    gu.id = rank;
    gu.name = "MGU_" + std::to_string(rank);
    gu.level_id = 0;
    gu.parent_id = -1;
    world.geo_units.push_back(gu);

    Venue v;
    v.id = rank;
    v.type_id = 0;
    v.geo_unit_id = rank;
    v.is_residence = true;
    world.venues.push_back(v);

    Person& p = world.people.emplace_back();
    p.id = rank;
    p.age = 30.0f;
    p.sex = Sex::MALE;
    p.geo_unit_id = rank;

    world.buildIndices();

    config.parallel.partition_level = "MGU";
    config.parallel.geo_unit_chunk_size = 1000;

    dm = std::make_unique<DomainManager>(world, config);
    dm->setMPI(rank, size);
    dm->setMaxPersonId(1);

    dm->setGeoUnitRank(0, 0);
    dm->setGeoUnitRank(1, 1);
    dm->setPersonRank(PERSON_R0, 0);
    dm->setPersonRank(PERSON_R1, 1);
    dm->setVenueRank(VENUE_R0, 0);
    dm->setVenueRank(VENUE_R1, 1);

    Domain& domain = dm->getDomain();
    domain.addGeoUnit(rank);
    domain.resident_ids.push_back(rank);
    domain.resident_set.insert(rank);
    domain.local_venue_ids.push_back(rank);
    domain.local_venue_set.insert(rank);
  }
};

// ---------------------------------------------------------------------------
// Helper: build a PersonLocation targeting the remote rank's venue
// ---------------------------------------------------------------------------
static PersonLocation makeRemoteLocation(int rank) {
  PersonLocation loc;
  loc.person_id = rank;
  loc.venue_id = 1 - rank;  // remote venue
  loc.subset_index = 0;
  loc.activity_index = 1;
  loc.encounter_type_id = 255;
  return loc;
}

// ---------------------------------------------------------------------------
// Helper: the receiving rank's VisitorInfo for an incoming visitor, as
// Simulator builds it, plus a home_array_index mapping back to the person
// ---------------------------------------------------------------------------
static VisitorInfo toVisitorInfo(const Domain::VisitorData& vis) {
  VisitorInfo vi;
  vi.person_id = vis.person_id;
  vi.is_infected = vis.is_infected;
  vi.is_infectious = vis.is_infectious;
  vi.immunity_level = vis.immunity_level;
  vi.home_array_index = vis.person_id;
  vi.symptom_id = vis.symptom_id;
  std::copy(std::begin(vis.integrated_infectiousness),
            std::end(vis.integrated_infectiousness),
            std::begin(vi.integrated_infectiousness));
  vi.fomite_deposition_sub = vis.fomite_deposition_sub;
  return vi;
}

// ---------------------------------------------------------------------------
// H1: Stage-driven visitor infects local susceptible
// ---------------------------------------------------------------------------
TEST_CASE("H1: Stage-driven visitor infects local susceptible") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(5.0);
  tp.stage_curves["mild"] = curve;
  tp.symptom_id_curves = {nullptr, curve};

  {
    TransmissionMode m;
    m.name = "default";
    m.symptom_curves = tp.symptom_id_curves;
    tp.modes.push_back(std::move(m));
  }

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("StageFlu", stags, {}, {td}, {}, tp);
  f.dm->setDisease(&disease);

  // Rank 0: infect person 0
  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    p->infection = std::make_unique<Infection>(
        &disease, -1.0, p, 42u, &f.world, "household", 0, 1.0f, 0, "general");
  }

  // Person visits remote venue
  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.0, 1.0);

  Domain& domain = f.dm->getDomain();
  REQUIRE(domain.incoming_visitors.size() == 1);

  // On rank 1: build visitor data and run processTransmissions
  if (f.rank == 1) {
    const auto& vis = domain.incoming_visitors[0];

    // Build visitor data map
    std::unordered_map<PersonId, VisitorInfo> visitor_data;
    visitor_data[vis.person_id] = toVisitorInfo(vis);

    std::unordered_set<PersonId> visitor_ids = {vis.person_id};

    // Locations: visitor (person 0) + local (person 1) at venue 1
    std::vector<PersonLocation> locs = {
        {vis.person_id, f.rank, -1, 0, 255, 0},  // visitor
        {f.rank, f.rank, -1, 0, 255, 0}          // local person (person 1)
    };

    ContactMatrixConfig cm;
    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{100.0}};
    cm.default_matrix = default_contact_matrix;
    SimulationConfig sim;
    ParallelConfig par;
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, f.world, disease);
    InteractionManager im(f.world, cm, sim, par, &disease, nullptr);

    std::vector<PendingInfection> pending;
    im.processTransmissions(locs, 0.0, 1.0, nullptr, &visitor_ids, &pending,
                            &visitor_data);

    // Person 1 (local) should be infected
    CHECK(f.world.getPerson(f.rank)->infection != nullptr);
  }
}

// ---------------------------------------------------------------------------
// H2: Trajectory-driven visitor infects local susceptible
// ---------------------------------------------------------------------------
TEST_CASE("H2: Trajectory-driven visitor infects local susceptible") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::TRAJECTORY_DRIVEN;
  tp.type = "gamma";
  tp.max_infectiousness = {"constant", {{"value", 10.0}}};
  tp.shape = {"constant", {{"value", 2.0}}};
  tp.rate = {"constant", {{"value", 1.0}}};
  tp.shift = {"constant", {{"value", 0.0}}};

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("TrajFlu", stags, {}, {td}, {}, tp);
  f.dm->setDisease(&disease);

  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    p->infection = std::make_unique<Infection>(
        &disease, -1.0, p, 42u, &f.world, "household", 0, 1.0f, 0, "general");
  }

  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.0, 1.0);

  Domain& domain = f.dm->getDomain();
  REQUIRE(domain.incoming_visitors.size() == 1);

  if (f.rank == 1) {
    const auto& vis = domain.incoming_visitors[0];

    std::unordered_map<PersonId, VisitorInfo> visitor_data;
    visitor_data[vis.person_id] = toVisitorInfo(vis);

    std::unordered_set<PersonId> visitor_ids = {vis.person_id};

    std::vector<PersonLocation> locs = {{vis.person_id, f.rank, -1, 0, 255, 0},
                                        {f.rank, f.rank, -1, 0, 255, 0}};

    ContactMatrixConfig cm;
    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{100.0}};
    cm.default_matrix = default_contact_matrix;
    SimulationConfig sim;
    ParallelConfig par;
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, f.world, disease);
    InteractionManager im(f.world, cm, sim, par, &disease, nullptr);

    std::vector<PendingInfection> pending;
    im.processTransmissions(locs, 0.0, 1.0, nullptr, &visitor_ids, &pending,
                            &visitor_data);

    CHECK(f.world.getPerson(f.rank)->infection != nullptr);
  }
}

// ---------------------------------------------------------------------------
// H3: Local infector infects visitor (pending infection routed back)
// ---------------------------------------------------------------------------
TEST_CASE("H3: Local infector infects visitor, pending routed back") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(5.0);
  tp.stage_curves["mild"] = curve;
  tp.symptom_id_curves = {nullptr, curve};

  {
    TransmissionMode m;
    m.name = "default";
    m.symptom_curves = tp.symptom_id_curves;
    tp.modes.push_back(std::move(m));
  }

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("StageFlu", stags, {}, {td}, {}, tp);
  f.dm->setDisease(&disease);

  // Rank 1: infect person 1 (local infector)
  if (f.rank == 1) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R1);
    p->infection = std::make_unique<Infection>(
        &disease, -1.0, p, 42u, &f.world, "household", 1, 1.0f, 0, "general");
  }

  // Person 0 (rank 0, susceptible) visits venue 1 (rank 1)
  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.0, 1.0);

  Domain& domain = f.dm->getDomain();

  std::vector<PendingInfection> pending;

  // On rank 1: visitor (person 0) + local infector (person 1) at venue 1
  if (f.rank == 1) {
    REQUIRE(domain.incoming_visitors.size() == 1);
    const auto& vis = domain.incoming_visitors[0];

    std::unordered_map<PersonId, VisitorInfo> visitor_data;
    visitor_data[vis.person_id] = toVisitorInfo(vis);

    std::unordered_set<PersonId> visitor_ids = {vis.person_id};

    std::vector<PersonLocation> locs = {
        {vis.person_id, f.rank, -1, 0, 255, 0},  // visitor (person 0)
        {f.rank, f.rank, -1, 0, 255, 0}          // local (person 1)
    };

    ContactMatrixConfig cm;
    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{100.0}};
    cm.default_matrix = default_contact_matrix;
    SimulationConfig sim;
    ParallelConfig par;
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, f.world, disease);
    InteractionManager im(f.world, cm, sim, par, &disease, nullptr);

    im.processTransmissions(locs, 0.0, 1.0, nullptr, &visitor_ids, &pending,
                            &visitor_data);
  }

  // Route pending infections back to home rank
  f.dm->receivePendingInfections(pending);

  // On rank 0: person 0 should now be infected
  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    CHECK(p->infection != nullptr);
  }
}

// ---------------------------------------------------------------------------
// H4: Stage-driven multi-mode: correct mode infectiousness across ranks
// ---------------------------------------------------------------------------
TEST_CASE("H4: Multi-mode stage-driven infectiousness across ranks") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve_resp = std::make_shared<ConstantCurve>(2.0);
  auto curve_bite = std::make_shared<ConstantCurve>(0.8);
  tp.stage_curves["mild"] = curve_resp;
  tp.symptom_id_curves = {nullptr, curve_resp};

  {
    TransmissionMode m;
    m.name = "respiratory";
    m.symptom_curves = {nullptr, curve_resp};
    tp.modes.push_back(std::move(m));
  }
  {
    TransmissionMode m;
    m.name = "animal_bite";
    m.symptom_curves = {nullptr, curve_bite};
    tp.modes.push_back(std::move(m));
  }

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("MultiModeFlu", stags, {}, {td}, {}, tp);
  f.dm->setDisease(&disease);

  // Rank 0: infect person 0
  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    p->infection = std::make_unique<Infection>(
        &disease, -1.0, p, 42u, &f.world, "household", 0, 1.0f, 0, "general");
  }

  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.0, 1.0);

  if (f.rank == 1) {
    const auto& vis = f.dm->getDomain().incoming_visitors[0];

    // Each mode's infectiousness arrives integrated over the 1 h slot:
    // 24 * rate * (1 / 24) d = rate.
    REQUIRE(vis.integrated_infectiousness.size() == 2);
    CHECK(vis.integrated_infectiousness[0] == doctest::Approx(2.0));
    CHECK(vis.integrated_infectiousness[1] == doctest::Approx(0.8));
  }
}

// ---------------------------------------------------------------------------
// H5: Transmission mode index preserved across ranks
// ---------------------------------------------------------------------------
TEST_CASE("H5: Transmission mode index preserved across ranks") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve_bite = std::make_shared<ConstantCurve>(1.0);
  auto curve_resp = std::make_shared<ConstantCurve>(2.0);
  tp.stage_curves["mild"] = curve_bite;
  tp.symptom_id_curves = {nullptr, curve_bite};

  {
    TransmissionMode m;
    m.name = "animal_bite";
    m.symptom_curves = {nullptr, curve_bite};
    tp.modes.push_back(std::move(m));
  }
  {
    TransmissionMode m;
    m.name = "respiratory";
    m.symptom_curves = {nullptr, curve_resp};
    tp.modes.push_back(std::move(m));
  }

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("Plague", stags, {}, {td}, {}, tp);
  f.dm->setDisease(&disease);

  // Person 0 visits rank 1's venue
  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.0, 1.0);

  // Rank 1 reports person 0 was infected via RESPIRATORY mode (index 1)
  std::vector<PendingInfection> pending;
  if (f.rank == 1) {
    const auto& vis = f.dm->getDomain().incoming_visitors[0];
    PendingInfection pi;
    pi.person_id = vis.person_id;
    pi.venue_type_id = 0;
    pi.venue_id = f.rank;
    pi.infector_symptom_id = 1;
    pi.transmission_mode_index = 1;  // Respiratory
    pending.push_back(pi);
  }

  f.dm->receivePendingInfections(pending);

  if (f.rank == 0) {
    Person* p = f.world.getPerson(0);
    REQUIRE(p->infection != nullptr);
    // Mode 1 (respiratory) has ConstantCurve(2.0)
    CHECK(p->infection->getInfectiousness(1, 1.0) == doctest::Approx(2.0));
    // Mode 0 (animal_bite) has ConstantCurve(1.0)
    CHECK(p->infection->getInfectiousness(0, 1.0) == doctest::Approx(1.0));
  }
}

// ---------------------------------------------------------------------------
// H6: Susceptible visitor with immunity resists cross-rank infection
// ---------------------------------------------------------------------------
TEST_CASE("H6: Immune visitor resists cross-rank infection") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(5.0);
  tp.stage_curves["mild"] = curve;
  tp.symptom_id_curves = {nullptr, curve};

  {
    TransmissionMode m;
    m.name = "default";
    m.symptom_curves = tp.symptom_id_curves;
    tp.modes.push_back(std::move(m));
  }

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("StageFlu", stags, {}, {td}, {}, tp);
  f.dm->setDisease(&disease);

  // Rank 1: infect person 1 (local infector)
  if (f.rank == 1) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R1);
    p->infection = std::make_unique<Infection>(
        &disease, -1.0, p, 42u, &f.world, "household", 1, 1.0f, 0, "general");
  }

  // Rank 0: person 0 has full immunity and visits rank 1
  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    p->immunity.natural_level = 1.0;
    p->immunity.natural_acquisition_time =
        0.0;  // Must be >= 0 (sentinel -1 = "never acquired")
    p->immunity.natural_waning_rate = 0.0;
  }

  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.0, 1.0);

  std::vector<PendingInfection> pending;

  if (f.rank == 1) {
    Domain& domain = f.dm->getDomain();
    REQUIRE(domain.incoming_visitors.size() == 1);
    const auto& vis = domain.incoming_visitors[0];

    // Visitor should have high immunity
    CHECK(vis.immunity_level == doctest::Approx(1.0).epsilon(0.01));

    std::unordered_map<PersonId, VisitorInfo> visitor_data;
    visitor_data[vis.person_id] = toVisitorInfo(vis);

    std::unordered_set<PersonId> visitor_ids = {vis.person_id};

    std::vector<PersonLocation> locs = {{vis.person_id, f.rank, -1, 0, 255, 0},
                                        {f.rank, f.rank, -1, 0, 255, 0}};

    ContactMatrixConfig cm;
    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{100.0}};
    cm.default_matrix = default_contact_matrix;
    SimulationConfig sim;
    ParallelConfig par;
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, f.world, disease);
    InteractionManager im(f.world, cm, sim, par, &disease, nullptr);

    im.processTransmissions(locs, 0.0, 1.0, nullptr, &visitor_ids, &pending,
                            &visitor_data);
  }

  // Route any pending infections back
  f.dm->receivePendingInfections(pending);

  // Person 0 should NOT be infected (immune)
  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    CHECK(p->infection == nullptr);
  }
}

// ---------------------------------------------------------------------------
// H7: Bidirectional cross-rank transmission
// ---------------------------------------------------------------------------
TEST_CASE("H7: Bidirectional cross-rank transmission") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(5.0);
  tp.stage_curves["mild"] = curve;
  tp.symptom_id_curves = {nullptr, curve};

  {
    TransmissionMode m;
    m.name = "default";
    m.symptom_curves = tp.symptom_id_curves;
    tp.modes.push_back(std::move(m));
  }

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("StageFlu", stags, {}, {td}, {}, tp);
  f.dm->setDisease(&disease);

  // Both persons are infectious
  Person* local_p = f.world.getPerson(f.rank);
  local_p->infection =
      std::make_unique<Infection>(&disease, -1.0, local_p, 42u, &f.world,
                                  "household", f.rank, 1.0f, 0, "general");

  // Each person visits the other rank's venue
  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.0, 1.0);

  Domain& domain = f.dm->getDomain();
  // Each rank should receive one visitor
  REQUIRE(domain.incoming_visitors.size() == 1);

  const auto& vis = domain.incoming_visitors[0];
  CHECK(vis.is_infectious == true);
  CHECK(vis.integrated_infectiousness[0] > 0.0);
}

// ---------------------------------------------------------------------------
// H8: a Visitor deposits exactly what an identical local would
// ---------------------------------------------------------------------------
// Symptom ids: 0 = healthy, 1 = exposed, 2 = mild. Healthy is infected but
// not infectious (incubation).
static constexpr uint16_t kHealthy = 0;
static constexpr uint16_t kExposed = 1;
static constexpr uint16_t kMild = 2;

// Direct mode 0 plus fomite mode 1, sub-binned at `sub_bin_time` hours. Each
// symptom deposits on a different curve; mild ramps, so a deposit depends on
// the exact time in stage.
static Disease makeFomiteDisease(double sub_bin_time) {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(5.0);
  tp.symptom_id_curves = {nullptr, curve, curve};

  TransmissionMode direct;
  direct.name = "direct";
  direct.symptom_curves = tp.symptom_id_curves;
  tp.modes.push_back(std::move(direct));

  TransmissionMode fomite;
  fomite.name = "fomite";
  fomite.type = TransmissionModeType::Fomite;
  fomite.symptom_curves = {nullptr, nullptr, nullptr};
  FomiteConfig fomite_config;
  fomite_config.mode_index = 1;
  fomite_config.max_age = 2.0;
  fomite_config.sub_bin_time = sub_bin_time;
  fomite_config.infectiousness_curve = std::make_shared<ConstantCurve>(1.0);
  fomite_config.deposition_by_symptom = {
      std::make_shared<ConstantCurve>(0.7),
      std::make_shared<ConstantCurve>(2.0),
      std::make_shared<LinearRampCurve>(1.0, 3.0, 1.0)};
  fomite.config = std::move(fomite_config);
  tp.modes.push_back(std::move(fomite));

  std::vector<SymptomTag> stags = {
      {"healthy", -1, kHealthy}, {"exposed", 0, kExposed}, {"mild", 1, kMild}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  return Disease("FomiteFlu", stags, {}, {td}, {}, tp);
}

// An Infection following exactly `transitions` (time, symptom id).
static std::unique_ptr<Infection> makeInfection(
    const Disease& disease, double infection_time,
    std::vector<std::pair<double, uint16_t>> transitions) {
  InfectionTrajectory trajectory;
  trajectory.infection_time = infection_time;
  trajectory.transitions = std::move(transitions);
  return Infection::fromCheckpoint(&disease, infection_time, trajectory, 1.0,
                                   1.0, 1.0, 0.0, /*last_checked_time=*/-1.0,
                                   kHealthy, infection_time);
}

// Fomite deposits made at `f`'s own venue by one slot over `locs`.
static std::deque<Venue::DepositEvent> depositsFromOneSlot(
    TwoRankFixture& f, Disease& disease,
    const std::vector<PersonLocation>& locs, double current_time,
    double delta_hours, const std::unordered_set<PersonId>* visitor_ids,
    std::vector<PendingInfection>* pending,
    const std::unordered_map<PersonId, VisitorInfo>* visitor_data) {
  Venue* venue = f.world.getVenue(f.rank);
  venue->fomite_history.assign(1, {});

  ContactMatrixConfig cm;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, f.world, disease);
  SimulationConfig sim;
  ParallelConfig par;
  InteractionManager im(f.world, cm, sim, par, &disease, nullptr);
  im.processTransmissions(locs, current_time, delta_hours, nullptr,
                          visitor_ids, pending, visitor_data);
  return venue->fomite_history[0];
}

// Person 0 visits rank 1 while person 1, in the same disease state, stays
// home at the same Venue. Their deposits there must be bitwise equal.
static void checkVisitorDepositsLikeLocal(
    double sub_bin_time, double infection_time,
    std::vector<std::pair<double, uint16_t>> transitions, double current_time,
    double delta_hours) {
  TwoRankFixture f;
  REQUIRE(f.size == 2);
  Disease disease = makeFomiteDisease(sub_bin_time);
  f.dm->setDisease(&disease);

  Person* local_person = f.world.getPerson(f.rank);
  local_person->infection =
      makeInfection(disease, infection_time, transitions);

  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, current_time,
                         delta_hours);
  if (f.rank != 1) return;

  const auto& incoming = f.dm->getDomain().incoming_visitors;
  REQUIRE(incoming.size() == 1);
  const PersonId visitor_id = incoming[0].person_id;
  std::unordered_map<PersonId, VisitorInfo> visitor_data = {
      {visitor_id, toVisitorInfo(incoming[0])}};
  std::unordered_set<PersonId> visitor_ids = {visitor_id};
  std::vector<PendingInfection> pending;

  auto visitor_deposits = depositsFromOneSlot(
      f, disease, {{visitor_id, f.rank, -1, 0, 255, 0}}, current_time,
      delta_hours, &visitor_ids, &pending, &visitor_data);
  auto local_deposits = depositsFromOneSlot(
      f, disease, {{f.rank, f.rank, -1, 0, 255, 0}}, current_time,
      delta_hours, nullptr, nullptr, nullptr);

  REQUIRE_FALSE(local_deposits.empty());
  REQUIRE(visitor_deposits.size() == local_deposits.size());
  for (size_t k = 0; k < local_deposits.size(); ++k) {
    CHECK(visitor_deposits[k].time == local_deposits[k].time);
    CHECK(visitor_deposits[k].amount == local_deposits[k].amount);
  }
}

TEST_CASE("H8: Visitor fomite deposit equals identical local's") {
  // Mild, part-way up the ramp; one sub-bin per slot.
  checkVisitorDepositsLikeLocal(/*sub_bin_time=*/0.0, -0.37,
                                {{-0.37, kExposed}, {-0.11, kMild}}, 0.3, 1.0);
}

TEST_CASE("H8b: Visitor fomite deposit equals local's across sub-bins") {
  // Three 2 h sub-bins; exposed -> mild inside the second.
  checkVisitorDepositsLikeLocal(/*sub_bin_time=*/2.0, 9.0,
                                {{9.0, kExposed}, {10.1, kMild}}, 10.0, 6.0);
}

TEST_CASE("H8c: Incubating Visitor alone still deposits") {
  // Infected, not infectious, and the only person at a fomite-free Venue.
  checkVisitorDepositsLikeLocal(/*sub_bin_time=*/0.0, 9.0,
                                {{9.0, kHealthy}, {11.0, kMild}}, 10.0, 6.0);
}

// ---------------------------------------------------------------------------
// H9: uninfected, incubating and infectious Visitors in one exchange
// ---------------------------------------------------------------------------
enum class VisitorState { Uninfected, Incubating, Infectious };
constexpr VisitorState kVisitorStates[] = {
    VisitorState::Uninfected, VisitorState::Incubating,
    VisitorState::Infectious};

// Person on `owner` in `state`: ids 2.. so the fixture's persons 0 and 1 stay.
static PersonId statePersonId(int owner, VisitorState state) {
  return 2 + 3 * owner + static_cast<int>(state);
}

static VisitorState stateOf(PersonId id) {
  return static_cast<VisitorState>((id - 2) % 3);
}

// Adds one person per state to each rank; this rank gets its own and an
// infection matching the state. Current time 10, 6 h slot.
static void addStatePeople(TwoRankFixture& f, const Disease& disease) {
  for (int owner = 0; owner < 2; ++owner) {
    for (VisitorState state : kVisitorStates) {
      const PersonId id = statePersonId(owner, state);
      f.dm->setPersonRank(id, owner);
      if (owner != f.rank) continue;
      Person& person = f.world.people.emplace_back();
      person.id = id;
      person.age = 30.0f;
      person.sex = Sex::MALE;
      person.geo_unit_id = f.rank;
      if (state == VisitorState::Incubating) {
        person.infection =
            makeInfection(disease, 9.0, {{9.0, kHealthy}, {11.0, kMild}});
      } else if (state == VisitorState::Infectious) {
        person.infection =
            makeInfection(disease, 9.0, {{9.0, kExposed}, {10.1, kMild}});
      }
      Domain& domain = f.dm->getDomain();
      domain.resident_ids.push_back(id);
      domain.resident_set.insert(id);
    }
  }
  f.dm->setMaxPersonId(statePersonId(1, VisitorState::Infectious));
  f.world.buildIndices();
}

// `senders` each send their three state people to the other rank's Venue.
// Every receiver checks each Visitor's tails against its state, and that the
// infected ones deposit exactly what their local twin on this rank does.
static void checkMixedStateExchange(const std::vector<int>& senders) {
  constexpr double kCurrentTime = 10.0;
  constexpr double kDeltaHours = 6.0;
  constexpr int kNumModes = 2;
  constexpr int kFomiteSubBins = 3;  // 6 h slot / 2 h sub-bins
  TwoRankFixture f;
  REQUIRE(f.size == 2);
  Disease disease = makeFomiteDisease(/*sub_bin_time=*/2.0);
  addStatePeople(f, disease);

  std::vector<PersonLocation> locations;
  const bool sends =
      std::find(senders.begin(), senders.end(), f.rank) != senders.end();
  if (sends) {
    for (VisitorState state : kVisitorStates) {
      PersonLocation location = makeRemoteLocation(f.rank);
      location.person_id = statePersonId(f.rank, state);
      locations.push_back(location);
    }
  }
  f.dm->exchangeVisitors(locations, disease, kCurrentTime, kDeltaHours);

  const auto& incoming = f.dm->getDomain().incoming_visitors;
  const bool receives =
      std::find(senders.begin(), senders.end(), 1 - f.rank) != senders.end();
  REQUIRE(incoming.size() == (receives ? 3u : 0u));

  for (const auto& visitor : incoming) {
    const VisitorState state = stateOf(visitor.person_id);
    const bool infected = state != VisitorState::Uninfected;
    const bool infectious = state == VisitorState::Infectious;
    CHECK(visitor.is_infected == infected);
    CHECK(visitor.is_infectious == infectious);
    // A tail arrives full length if its header gate sends it, else empty.
    CHECK(visitor.integrated_infectiousness.size() ==
          (infectious ? kNumModes : 0u));
    REQUIRE(visitor.fomite_deposition_sub.size() ==
            (infected ? kFomiteSubBins : 0u));

    const VisitorInfo info = toVisitorInfo(visitor);
    for (double integrated : info.integrated_infectiousness) {
      if (!infectious) CHECK(integrated == 0.0);
    }
    if (!infected) continue;

    std::unordered_map<PersonId, VisitorInfo> visitor_data = {
        {visitor.person_id, info}};
    std::unordered_set<PersonId> visitor_ids = {visitor.person_id};
    std::vector<PendingInfection> pending;
    auto visitor_deposits = depositsFromOneSlot(
        f, disease, {{visitor.person_id, f.rank, -1, 0, 255, 0}},
        kCurrentTime, kDeltaHours, &visitor_ids, &pending, &visitor_data);

    const PersonId twin_id = statePersonId(f.rank, state);
    const size_t twin_index = f.world.person_index.at(twin_id);
    auto local_deposits = depositsFromOneSlot(
        f, disease, {{twin_id, f.rank, -1, 0, 255, twin_index}}, kCurrentTime,
        kDeltaHours, nullptr, nullptr, nullptr);

    REQUIRE_FALSE(local_deposits.empty());
    REQUIRE(visitor_deposits.size() == local_deposits.size());
    for (size_t k = 0; k < local_deposits.size(); ++k) {
      CHECK(visitor_deposits[k].time == local_deposits[k].time);
      CHECK(visitor_deposits[k].amount == local_deposits[k].amount);
    }
  }
}

TEST_CASE("H9: mixed-state Visitors, one-way (P2P)") {
  checkMixedStateExchange({0});
}

TEST_CASE("H9b: mixed-state Visitors, two-way (all-to-all)") {
  checkMixedStateExchange({0, 1});
}

#endif  // USE_MPI

// ---------------------------------------------------------------------------
// Custom main: wrap MPI init/finalize around doctest
// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

  doctest::Context ctx;
  ctx.applyCommandLine(argc, argv);
  int result = ctx.run();

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return result;
}
