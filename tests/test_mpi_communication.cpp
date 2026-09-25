// MPI regression test for cross-rank visitor exchange and infection routing.
//
// Must be run with exactly 2 MPI ranks:
//   mpirun -np 2 ./test_mpi_communication
//
// Tests cover:
//   1. exchangeVisitors   – a person on rank 0 assigned to a venue on rank 1
//                           arrives in rank 1's incoming_visitors (and vice
//                           versa).
//   2. Infectious visitor – infectious status is correctly propagated across
//                           the rank boundary.
//   3. receivePendingInfections – an infection that occurred at a cross-rank
//                           venue is routed back and applied to the person on
//                           their home rank.

#define DOCTEST_CONFIG_IMPLEMENT  // custom main so we can wrap MPI
                                  // init/finalize
#include "doctest.h"

#ifdef USE_MPI
#include <mpi.h>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/disease.h"
#include "parallel/domain.h"
#include "parallel/domain_manager.h"

using namespace june;

// ---------------------------------------------------------------------------
// Build a minimal 2-rank world.
//
// rank 0 owns: geo_unit 0, venue 0, person 0
// rank 1 owns: geo_unit 1, venue 1, person 1
//
// Each rank loads only its own venue into world.venues (mimicking the real
// distributed HDF5 load), but knows about the remote venue via
// global_venue_rank_.
// ---------------------------------------------------------------------------
struct TwoRankFixture {
  int rank;
  int size;

  WorldState world;
  Config config;
  std::unique_ptr<DomainManager> dm;

  // Ids that are the same on every rank
  static constexpr PersonId PERSON_R0 = 0;
  static constexpr PersonId PERSON_R1 = 1;
  static constexpr VenueId VENUE_R0 = 0;
  static constexpr VenueId VENUE_R1 = 1;

  TwoRankFixture() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // --- Registries ---
    world.venue_type_names = {"household"};
    world.geo_level_names = {"MGU"};
    world.activity_names = {"residence", "work", "visiting", "none", "dead"};
    world.encounter_type_names = {};
    world.subset_type_names = {};

    // --- My geo_unit (each rank loads its own) ---
    GeographicalUnit gu;
    gu.id = rank;
    gu.name = "MGU_" + std::to_string(rank);
    gu.level_id = 0;
    gu.parent_id = -1;
    world.geo_units.push_back(gu);

    // --- My local venue only ---
    Venue v;
    v.id = rank;  // venue 0 on rank 0, venue 1 on rank 1
    v.type_id = 0;
    v.geo_unit_id = rank;
    v.is_residence = true;
    world.venues.push_back(v);

    // --- My local person ---
    Person& p = world.people.emplace_back();
    p.id = rank;
    p.age = 30.0f;
    p.sex = Sex::MALE;
    p.geo_unit_id = rank;

    world.buildIndices();

    // --- Config ---
    config.parallel.partition_level = "MGU";
    config.parallel.geo_unit_chunk_size = 1000;

    // --- DomainManager (test mode: bypass HDF5 initialize()) ---
    dm = std::make_unique<DomainManager>(world, config);
    dm->setMPI(rank, size);
    dm->setMaxPersonId(1);  // max person id is 1

    // Teach each rank about who owns what
    dm->setGeoUnitRank(0, 0);
    dm->setGeoUnitRank(1, 1);
    dm->setPersonRank(PERSON_R0, 0);
    dm->setPersonRank(PERSON_R1, 1);
    dm->setVenueRank(VENUE_R0, 0);
    dm->setVenueRank(VENUE_R1, 1);

    // Populate domain ownership directly (mirrors assignPeopleAndVenues)
    Domain& domain = dm->getDomain();
    domain.addGeoUnit(rank);
    domain.resident_ids.push_back(rank);
    domain.resident_set.insert(rank);
    domain.local_venue_ids.push_back(rank);
    domain.local_venue_set.insert(rank);
  }
};

// ---------------------------------------------------------------------------
// TEST 4: Stage-Driven Infectiousness Propagation
// ---------------------------------------------------------------------------
TEST_CASE(
    "exchangeVisitors: stage-driven infectiousness propagates correctly") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;

  // Define a curve that is NOT constant 1.0 to ensure we're actually
  // calculating it Mode 0 (default): linear ramp from 0 to 1 over 10 days
  auto curve = std::make_shared<LinearRampCurve>(0.0, 1.0, 10.0);
  tp.stage_curves["mild"] = curve;
  tp.symptom_id_curves = {nullptr, curve};

  TransmissionMode default_mode;
  default_mode.name = "default";
  default_mode.symptom_curves = tp.symptom_id_curves;
  tp.modes.push_back(std::move(default_mode));

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  std::vector<TrajectoryDefinition> trajectories = {td};
  DiseaseStageSettings stage_settings;
  OutcomeRates outcome_rates;
  Disease disease("StageFlu", stags, stage_settings, trajectories,
                  outcome_rates, tp);

  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    // Infected at time -5.0, so at time 0.0 five days into "mild"
    p->infection = std::make_unique<Infection>(
        &disease, -5.0, p, 42u, &f.world, "household", 0, 1.0f, 0, "general");
  }

  int remote_venue = 1 - f.rank;
  PersonLocation loc;
  loc.person_id = f.rank;
  loc.venue_id = remote_venue;
  loc.subset_index = 0;
  loc.activity_index = 1;
  loc.encounter_type_id = 255;

  // Pass a realistic delta_hours (6h slot) so the integration window is
  // non-zero.  exchangeVisitors pre-computes integrated_infectiousness over
  // [current_time, current_time + delta_hours/24].
  constexpr double delta_hours = 6.0;
  constexpr double current_time = 0.0;
  f.dm->exchangeVisitors({loc}, disease, current_time, delta_hours);

  Domain& domain = f.dm->getDomain();
  REQUIRE(domain.incoming_visitors.size() == 1);

  // Compute expected integrated infectiousness on rank 0 (where the person
  // lives) and broadcast to rank 1 for comparison.  This tests the core
  // invariant: the visitor's pre-computed value must exactly match what the
  // home rank's Infection object produces.
  double expected_ii = 0.0;
  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    expected_ii = p->infection->getIntegratedInfectiousness(
        0, current_time, current_time + delta_hours / 24.0);
  }
  MPI_Bcast(&expected_ii, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  CHECK(expected_ii > 0.0);  // sanity: non-trivial curve value

  if (f.rank == 1) {
    const auto& v = domain.incoming_visitors[0];
    CHECK(v.person_id == 0);
    CHECK(v.is_infectious == true);
    CHECK(v.symptom_id == 1);
    CHECK(v.integrated_infectiousness[0] == doctest::Approx(expected_ii));
  }
}

// ---------------------------------------------------------------------------
// TEST 5: Trajectory-Driven Infectiousness Propagation
// ---------------------------------------------------------------------------
TEST_CASE(
    "exchangeVisitors: trajectory-driven infectiousness propagates correctly") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::TRAJECTORY_DRIVEN;
  tp.type = "gamma";
  tp.max_infectiousness = {"constant", {{"value", 2.0}}};
  tp.shape = {"constant", {{"value", 2.0}}};
  tp.rate = {"constant", {{"value", 1.0}}};
  tp.shift = {"constant", {{"value", 0.0}}};

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("TrajFlu", stags, {}, {td}, {}, tp);

  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    // Infected at time -1.0.
    // Gamma PDF(shape=2, rate=1) at x=1 is 1.0 * exp(-1) = 0.3678...
    // Max infectiousness = 2.0, so 2.0 * 0.3678 = 0.7357...
    p->infection = std::make_unique<Infection>(
        &disease, -1.0, p, 42u, &f.world, "household", 0, 1.0f, 0, "general");
  }

  constexpr double delta_hours = 6.0;
  constexpr double current_time = 0.0;

  int remote_venue = 1 - f.rank;
  PersonLocation loc;
  loc.person_id = f.rank;
  loc.venue_id = remote_venue;
  loc.subset_index = 0;
  loc.activity_index = 1;
  f.dm->exchangeVisitors({loc}, disease, current_time, delta_hours);

  // Compute expected value on rank 0 and broadcast — verifying the
  // pre-computation invariant across ranks.
  double expected_ii = 0.0;
  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    expected_ii = p->infection->getIntegratedInfectiousness(
        0, current_time, current_time + delta_hours / 24.0);
  }
  MPI_Bcast(&expected_ii, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  CHECK(expected_ii > 0.0);

  if (f.rank == 1) {
    const auto& v = f.dm->getDomain().incoming_visitors[0];
    CHECK(v.is_infectious == true);
    CHECK(v.integrated_infectiousness[0] ==
          doctest::Approx(expected_ii).epsilon(0.001));
  }
}

// ---------------------------------------------------------------------------
// TEST 6: Transmission Mode Routing (Animal Bite / Respiratory)
// ---------------------------------------------------------------------------
TEST_CASE("receivePendingInfections: transmission mode index is preserved") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve_bite = std::make_shared<ConstantCurve>(1.0);
  auto curve_resp = std::make_shared<ConstantCurve>(2.0);

  TransmissionMode bite_mode;
  bite_mode.name = "animal_bite";
  bite_mode.symptom_curves = {nullptr, curve_bite};
  tp.modes.push_back(std::move(bite_mode));

  TransmissionMode resp_mode;
  resp_mode.name = "respiratory";
  resp_mode.symptom_curves = {nullptr, curve_resp};
  tp.modes.push_back(std::move(resp_mode));

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("Plague", stags, {}, {td}, {}, tp);

  // Step 1: Person 0 visits Rank 1's venue
  int remote_venue = 1 - f.rank;
  PersonLocation loc;
  loc.person_id = f.rank;
  loc.venue_id = remote_venue;
  loc.subset_index = 0;
  loc.activity_index = 1;
  f.dm->exchangeVisitors({loc}, disease, 0.0);

  // Step 2: Rank 1 reports that Person 0 was infected via RESPIRATORY mode
  // (index 1)
  std::vector<PendingInfection> pending;
  if (f.rank == 1) {
    const auto& v = f.dm->getDomain().incoming_visitors[0];
    PendingInfection pi;
    pi.person_id = v.person_id;
    pi.infection_time = 0.5;
    pi.venue_type_id = 0;
    pi.venue_id = f.rank;
    pi.infector_symptom_id = 1;
    pi.transmission_mode_index = 1;  // Respiratory!
    pending.push_back(pi);
  }

  f.dm->receivePendingInfections(pending, disease);

  // Step 3: Verify Person 0 on Rank 0 now has an infection with mode index 1
  if (f.rank == 0) {
    Person* p = f.world.getPerson(0);
    REQUIRE(p->infection != nullptr);
    // We can't directly query transmission_mode_index from Infection class
    // (it's private) but we can verify it via the resulting infectiousness if
    // we use stage-driven curves. At current_time 1.0, with infection_time 0.5,
    // time_in_stage is 0.5. Mode 1 (respiratory) has constant 2.0. Mode 0 has
    // constant 1.0.
    CHECK(p->infection->getInfectiousness(1, 1.0) == doctest::Approx(2.0));
    CHECK(p->infection->getInfectiousness(0, 1.0) == doctest::Approx(1.0));
  }
}

// ---------------------------------------------------------------------------
// TEST 7: Multi-record batch to the same destination rank (wire-format
// regression — see mission mpi_pending_infection_wire_format_desync)
// ---------------------------------------------------------------------------
TEST_CASE(
    "receivePendingInfections: two records to the same destination rank "
    "round-trip independently") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  constexpr PersonId PERSON_R0_SECOND = 2;

  // Rank 0 gains a second local resident so it can send two visitors to
  // rank 1's venue in the same round, forcing receivePendingInfections to
  // pack two PendingInfection records into one destination slot.
  if (f.rank == 0) {
    Person& p2 = f.world.people.emplace_back();
    p2.id = PERSON_R0_SECOND;
    p2.age = 30.0f;
    p2.sex = Sex::MALE;
    p2.geo_unit_id = f.rank;
    f.world.buildIndices();

    Domain& domain = f.dm->getDomain();
    domain.resident_ids.push_back(PERSON_R0_SECOND);
    domain.resident_set.insert(PERSON_R0_SECOND);
  }

  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general";
  td.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease("Plague2", stags, {}, {td}, {}, tp);

  // Step 1: rank 0's two people (0 and 2) visit rank 1's venue.
  std::vector<PersonLocation> locations;
  if (f.rank == 0) {
    PersonLocation loc0;
    loc0.person_id = 0;
    loc0.venue_id = TwoRankFixture::VENUE_R1;
    loc0.subset_index = 0;
    loc0.activity_index = 1;
    locations.push_back(loc0);

    PersonLocation loc2;
    loc2.person_id = PERSON_R0_SECOND;
    loc2.venue_id = TwoRankFixture::VENUE_R1;
    loc2.subset_index = 0;
    loc2.activity_index = 1;
    locations.push_back(loc2);
  }
  f.dm->exchangeVisitors(locations, disease, 0.0);

  if (f.rank == 1) {
    REQUIRE(f.dm->getDomain().incoming_visitors.size() == 2);
  }

  // Step 2: rank 1 reports two distinct infections, both destined for
  // rank 0, in one receivePendingInfections round.
  std::vector<PendingInfection> pending;
  if (f.rank == 1) {
    PendingInfection pi0;
    pi0.person_id = 0;
    pi0.infector_id = 100;
    pi0.infection_time = 0.5;
    pi0.venue_type_id = 0;
    pi0.encounter_type_id = 0;
    pi0.venue_id = TwoRankFixture::VENUE_R1;
    pi0.infector_symptom_id = 1;
    pi0.transmission_mode_index = 0;
    pending.push_back(pi0);

    PendingInfection pi2;
    pi2.person_id = PERSON_R0_SECOND;
    pi2.infector_id = 200;
    pi2.infection_time = 1.5;
    pi2.venue_type_id = 0;
    pi2.encounter_type_id = 1;
    pi2.venue_id = TwoRankFixture::VENUE_R1;
    pi2.infector_symptom_id = 2;
    pi2.transmission_mode_index = 1;
    pending.push_back(pi2);
  }

  auto applied = f.dm->receivePendingInfections(pending, disease);

  // Step 3: rank 0 must see both records with fields intact and
  // independent of each other — a wire-format desync corrupts the second
  // record with bytes borrowed from the first.
  if (f.rank == 0) {
    REQUIRE(applied.size() == 2);

    const PendingInfection* a0 = nullptr;
    const PendingInfection* a2 = nullptr;
    for (const auto& a : applied) {
      if (a.person_id == 0) a0 = &a;
      if (a.person_id == PERSON_R0_SECOND) a2 = &a;
    }
    REQUIRE(a0 != nullptr);
    REQUIRE(a2 != nullptr);

    CHECK(a0->infector_id == 100);
    CHECK(a0->infection_time == doctest::Approx(0.5));
    CHECK(a0->encounter_type_id == 0);
    CHECK(a0->infector_symptom_id == 1);
    CHECK(a0->transmission_mode_index == 0);

    CHECK(a2->infector_id == 200);
    CHECK(a2->infection_time == doctest::Approx(1.5));
    CHECK(a2->encounter_type_id == 1);
    CHECK(a2->infector_symptom_id == 2);
    CHECK(a2->transmission_mode_index == 1);
  }
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
