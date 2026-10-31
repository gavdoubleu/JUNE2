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
#include "transmission_fixtures.h"

#ifdef USE_MPI
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <optional>

#include "epidemiology/interaction_manager.h"
#include "epidemiology/policy.h"
#include "mpi_test_helpers.h"
#include "test_utils.h"

using namespace june;
using fomite_flu::kExposed;
using fomite_flu::kHealthy;
using fomite_flu::kMild;

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
  vi.immunity_level = vis.immunity_level;
  vi.home_array_index = vis.person_id;
  vi.age = vis.age;
  vi.sex = static_cast<Sex>(vis.sex_code);
  vi.symptom_id = vis.symptom_id;
  vi.time_in_stage = vis.time_in_stage;
  vi.emission = vis.emission;
  vi.target_susceptibility = vis.target_susceptibility;
  return vi;
}

// ---------------------------------------------------------------------------
// H1: Stage-driven visitor infects local susceptible
// ---------------------------------------------------------------------------
TEST_CASE("H1: Stage-driven visitor infects local susceptible") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  Disease disease = makeStageDisease(std::make_shared<ConstantCurve>(5.0));

  // The visitor wire carries effective values, not the sender's rank-local
  // modifier-set ID. Source and target channels deliberately differ here so
  // the test catches a collapsed scalar payload.
  PolicyManager policy_manager(f.world);
  TemporalPolicy policy;
  policy.name = "visitor_channels";
  policy.action.compliance_rate = 1.0;
  policy_manager.addTemporalPolicy(policy);
  PolicyTransmissionEffect source_effect;
  source_effect.policy_index = 0;
  source_effect.scope = TransmissionEffectScope::Person;
  source_effect.mode_name = "default";
  source_effect.channel = TransmissionEffectChannel::SourceInfectiousness;
  source_effect.multiplier = 0.50;
  policy_manager.addTransmissionEffect(source_effect);
  PolicyTransmissionEffect target_effect = source_effect;
  target_effect.channel = TransmissionEffectChannel::TargetSusceptibility;
  target_effect.multiplier = 0.25;
  policy_manager.addTransmissionEffect(target_effect);
  policy_manager.precomputePolicyApplicability(f.world.people);
  policy_manager.initializeTransmissionModifiers(disease, 0.0);
  f.dm->setPolicyManager(&policy_manager);

  // Rank 0: infect person 0
  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    p->infection = std::make_unique<Infection>(&disease, -1.0, p, 42u,
                                               kNoTransmissionContext, &f.world,
                                               "household", 0, "general");
  }

  double expected_integrated = 0.0;
  if (f.rank == 0) {
    expected_integrated =
        f.world.people[0].infection->getIntegratedInfectiousness(0, 0.0,
                                                                 1.0 / 24.0);
  }
  MPI_Bcast(&expected_integrated, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
    CHECK(vis.emission.infectiousness_by_mode[0] ==
          doctest::Approx(expected_integrated * 0.50));
    CHECK(vis.target_susceptibility[0] == doctest::Approx(0.25));
  }
}

// ---------------------------------------------------------------------------
// H2: Trajectory-driven visitor infects local susceptible
// ---------------------------------------------------------------------------
TEST_CASE("H2: Trajectory-driven visitor infects local susceptible") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  Disease disease = makeTrajectoryDisease(/*max_infectiousness=*/10.0);

  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    p->infection = std::make_unique<Infection>(&disease, -1.0, p, 42u,
                                               kNoTransmissionContext, &f.world,
                                               "household", 0, "general");
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

  Disease disease = makeStageDisease(std::make_shared<ConstantCurve>(5.0));

  // Rank 1: infect person 1 (local infector)
  if (f.rank == 1) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R1);
    p->infection = std::make_unique<Infection>(&disease, -1.0, p, 42u,
                                               kNoTransmissionContext, &f.world,
                                               "household", 1, "general");
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
  f.dm->receivePendingInfections(pending, disease);

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

  Disease disease =
      makeTwoModeStageDisease({"respiratory", 2.0}, {"animal_bite", 0.8});

  // Rank 0: infect person 0
  if (f.rank == 0) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R0);
    p->infection = std::make_unique<Infection>(&disease, -1.0, p, 42u,
                                               kNoTransmissionContext, &f.world,
                                               "household", 0, "general");
  }

  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.0, 1.0);

  if (f.rank == 1) {
    const auto& vis = f.dm->getDomain().incoming_visitors[0];

    // Each mode's infectiousness arrives integrated over the 1 h slot:
    // 24 * rate * (1 / 24) d = rate.
    REQUIRE(vis.emission.infectiousness_by_mode.size() == 2);
    CHECK(vis.emission.infectiousness_by_mode[0] == doctest::Approx(2.0));
    CHECK(vis.emission.infectiousness_by_mode[1] == doctest::Approx(0.8));
  }
}

// ---------------------------------------------------------------------------
// H5: Transmission mode index preserved across ranks
// ---------------------------------------------------------------------------
TEST_CASE("H5: Transmission mode index preserved across ranks") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  Disease disease =
      makeTwoModeStageDisease({"animal_bite", 1.0}, {"respiratory", 2.0});

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

  f.dm->receivePendingInfections(pending, disease);

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

  Disease disease = makeStageDisease(std::make_shared<ConstantCurve>(5.0));

  // Rank 1: infect person 1 (local infector)
  if (f.rank == 1) {
    Person* p = f.world.getPerson(TwoRankFixture::PERSON_R1);
    p->infection = std::make_unique<Infection>(&disease, -1.0, p, 42u,
                                               kNoTransmissionContext, &f.world,
                                               "household", 1, "general");
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
  f.dm->receivePendingInfections(pending, disease);

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

  Disease disease = makeStageDisease(std::make_shared<ConstantCurve>(5.0));

  // Both persons are infectious
  Person* local_p = f.world.getPerson(f.rank);
  local_p->infection = std::make_unique<Infection>(
      &disease, -1.0, local_p, 42u, kNoTransmissionContext, &f.world,
      "household", f.rank, "general");

  // Each person visits the other rank's venue
  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.0, 1.0);

  Domain& domain = f.dm->getDomain();
  // Each rank should receive one visitor
  REQUIRE(domain.incoming_visitors.size() == 1);

  const auto& vis = domain.incoming_visitors[0];
  CHECK(vis.is_infected == true);
  CHECK(vis.emission.infectiousness_by_mode[0] > 0.0);
}

// ---------------------------------------------------------------------------
// H8: a Visitor deposits exactly what an identical local would
// ---------------------------------------------------------------------------
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
    const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
    const PolicyManager* policy_manager = nullptr) {
  Venue* venue = f.world.getVenue(f.rank);
  venue->fomite_history.assign(1, {});

  ContactMatrixConfig cm;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, f.world, disease);
  SimulationConfig sim;
  ParallelConfig par;
  InteractionManager im(f.world, cm, sim, par, &disease, nullptr);
  im.setPolicyManager(policy_manager);
  im.processTransmissions(locs, current_time, delta_hours, nullptr, visitor_ids,
                          pending, visitor_data);
  return venue->fomite_history[0];
}

// Person 0 visits rank 1 while person 1, in the same disease state, stays
// home at the same Venue. Their deposits there must be bitwise equal. A
// fomite_source_multiplier other than 1 applies a person-scoped
// SourceInfectiousness policy on the fomite mode to both.
static void checkVisitorDepositsLikeLocal(
    double sub_bin_time, double infection_time,
    std::vector<std::pair<double, uint16_t>> transitions, double current_time,
    double delta_hours, double fomite_source_multiplier = 1.0) {
  TwoRankFixture f;
  REQUIRE(f.size == 2);
  Disease disease = makeFomiteDisease(sub_bin_time);

  std::optional<PolicyManager> policy_manager;
  if (fomite_source_multiplier != 1.0) {
    policy_manager.emplace(f.world);
    TemporalPolicy policy;
    policy.name = "fomite_source";
    policy.action.compliance_rate = 1.0;
    policy_manager->addTemporalPolicy(policy);
    PolicyTransmissionEffect source_effect;
    source_effect.policy_index = 0;
    source_effect.scope = TransmissionEffectScope::Person;
    source_effect.mode_name = "fomite";
    source_effect.channel = TransmissionEffectChannel::SourceInfectiousness;
    source_effect.multiplier = fomite_source_multiplier;
    policy_manager->addTransmissionEffect(source_effect);
    policy_manager->precomputePolicyApplicability(f.world.people);
    policy_manager->initializeTransmissionModifiers(disease, 0.0);
    f.dm->setPolicyManager(&*policy_manager);
  }
  const PolicyManager* policy = policy_manager ? &*policy_manager : nullptr;

  Person* local_person = f.world.getPerson(f.rank);
  local_person->infection = makeInfection(disease, infection_time, transitions);

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
      delta_hours, &visitor_ids, &pending, &visitor_data, policy);
  auto local_deposits = depositsFromOneSlot(
      f, disease, {{f.rank, f.rank, -1, 0, 255, 0}}, current_time, delta_hours,
      nullptr, nullptr, nullptr, policy);

  REQUIRE_FALSE(local_deposits.empty());
  REQUIRE(visitor_deposits.size() == local_deposits.size());
  for (size_t k = 0; k < local_deposits.size(); ++k) {
    CHECK(visitor_deposits[k].time == local_deposits[k].time);
    CHECK(visitor_deposits[k].amount == local_deposits[k].amount);
  }

  if (policy) {
    auto unmodified_deposits = depositsFromOneSlot(
        f, disease, {{f.rank, f.rank, -1, 0, 255, 0}}, current_time,
        delta_hours, nullptr, nullptr, nullptr, nullptr);
    REQUIRE(unmodified_deposits.size() == local_deposits.size());
    for (size_t k = 0; k < local_deposits.size(); ++k) {
      CHECK(local_deposits[k].amount ==
            doctest::Approx(unmodified_deposits[k].amount *
                            fomite_source_multiplier));
    }
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

TEST_CASE("H8d: Source policy scales Visitor and local deposits alike") {
  // Home rank applies the visitor's fomite source modifier, as for locals.
  checkVisitorDepositsLikeLocal(/*sub_bin_time=*/2.0, 9.0,
                                {{9.0, kExposed}, {10.1, kMild}}, 10.0, 6.0,
                                /*fomite_source_multiplier=*/0.3);
}

// ---------------------------------------------------------------------------
// H9: uninfected, incubating and infectious Visitors in one exchange
// ---------------------------------------------------------------------------
enum class VisitorState { Uninfected, Incubating, Infectious };
constexpr VisitorState kVisitorStates[] = {VisitorState::Uninfected,
                                           VisitorState::Incubating,
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
    CHECK(visitor.is_infected == infected);
    // Both tails arrive full length if the Visitor is infected, else empty.
    CHECK(visitor.emission.infectiousness_by_mode.size() ==
          (infected ? kNumModes : 0u));
    REQUIRE(visitor.emission.fomite_deposits.size() ==
            (infected ? kFomiteSubBins : 0u));
    // target_susceptibility always travels; the source multiplier with the
    // deposits.
    CHECK(visitor.target_susceptibility.size() == kNumModes);
    CHECK(visitor.deposition_source_multiplier.size() == (infected ? 1u : 0u));

    const VisitorInfo info = toVisitorInfo(visitor);
    CHECK(info.emission == visitor.emission);
    if (!infected) continue;

    std::unordered_map<PersonId, VisitorInfo> visitor_data = {
        {visitor.person_id, info}};
    std::unordered_set<PersonId> visitor_ids = {visitor.person_id};
    std::vector<PendingInfection> pending;
    auto visitor_deposits = depositsFromOneSlot(
        f, disease, {{visitor.person_id, f.rank, -1, 0, 255, 0}}, kCurrentTime,
        kDeltaHours, &visitor_ids, &pending, &visitor_data);

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

// ---------------------------------------------------------------------------
// H12: a Visitor goes into the contact-matrix bin its home rank would give it
// ---------------------------------------------------------------------------
// Person 0 visits rank 1's household Venue while person 1, identical in age
// and sex, stays there. An infectious local is also present. The matrix gives
// contacts only to observers in one bin, with counts so large that a person
// in that bin is infected for certain and anyone else is not, so the outcome
// says which bin the person was put in. Visitor and local must agree, and
// agree with the bin their age or sex selects.
static ContactMatrixConfig binnedMatrix(
    const WorldState& world, const Disease& disease,
    const std::vector<std::string>& bins,
    const std::vector<std::vector<double>>& contacts) {
  ContactMatrix matrix;
  matrix.bins = bins;
  matrix.contacts = contacts;
  ContactMatrixConfig cm;
  cm.matrices["household"] = matrix;
  cm.mode_matrices["household"]["default"] = matrix;
  cm.default_matrix = matrix;
  cm.mode_names = {"default"};
  finalizeContactMatrices(cm, world, disease);
  return cm;
}

// Whether the Visitor and the local twin were each infected in one slot.
struct BinOutcome {
  bool evaluated = false;  // only the rank holding the Venue evaluates
  bool visitor_infected = false;
  bool local_infected = false;
};

static BinOutcome infectionByBin(
    Sex sex, float age, const std::vector<std::string>& bins,
    const std::vector<std::vector<double>>& contacts) {
  constexpr PersonId kSource = 2;
  constexpr double kCurrentTime = 10.5;
  constexpr double kDeltaHours = 1.0;
  TwoRankFixture f;
  REQUIRE(f.size == 2);
  Disease disease = makeStageDisease(std::make_shared<ConstantCurve>(5.0));

  f.world.getPerson(f.rank)->sex = sex;
  f.world.getPerson(f.rank)->age = age;
  f.dm->setPersonRank(kSource, 1);
  if (f.rank == 1) {
    Person& source = f.world.people.emplace_back();
    source.id = kSource;
    source.age = 30.0f;
    source.sex = Sex::MALE;
    source.geo_unit_id = f.rank;
    source.infection = makeInfection(disease, 9.0, {{9.0, 1}});
    Domain& domain = f.dm->getDomain();
    domain.resident_ids.push_back(kSource);
    domain.resident_set.insert(kSource);
  }
  f.dm->setMaxPersonId(kSource);
  f.world.buildIndices();

  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, kCurrentTime,
                         kDeltaHours);
  BinOutcome outcome;
  if (f.rank != 1) return outcome;
  outcome.evaluated = true;

  const auto& incoming = f.dm->getDomain().incoming_visitors;
  REQUIRE(incoming.size() == 1);
  const PersonId visitor_id = incoming[0].person_id;
  std::unordered_map<PersonId, VisitorInfo> visitor_data = {
      {visitor_id, toVisitorInfo(incoming[0])}};
  std::unordered_set<PersonId> visitor_ids = {visitor_id};

  ContactMatrixConfig cm = binnedMatrix(f.world, disease, bins, contacts);
  SimulationConfig sim;
  ParallelConfig par;
  const size_t source_index = f.world.person_index.at(kSource);
  const PersonLocation source_loc{kSource, f.rank, -1, 0, 255, source_index};

  std::vector<PendingInfection> pending;
  {
    InteractionManager im(f.world, cm, sim, par, &disease, nullptr);
    im.processTransmissions({source_loc, {visitor_id, f.rank, -1, 0, 255, 0}},
                            kCurrentTime, kDeltaHours, nullptr, &visitor_ids,
                            &pending, &visitor_data);
  }
  for (const auto& infection : pending) {
    if (infection.person_id == visitor_id) outcome.visitor_infected = true;
  }

  Person* twin = f.world.getPerson(f.rank);
  const size_t twin_index = f.world.person_index.at(twin->id);
  {
    InteractionManager im(f.world, cm, sim, par, &disease, nullptr);
    im.processTransmissions(
        {source_loc, {twin->id, f.rank, -1, 0, 255, twin_index}}, kCurrentTime,
        kDeltaHours, nullptr, nullptr, nullptr, nullptr);
  }
  outcome.local_infected = twin->infection != nullptr;
  return outcome;
}

TEST_CASE("H12: Visitor is binned by sex like an identical local") {
  struct SexCase {
    std::vector<std::string> bins;
    std::vector<std::vector<double>> contacts;
    Sex infectable_sex;  // the sex whose bin has contacts
  };
  // Contacts sit on the male row in one scheme and the female row in the
  // other, so a Visitor dropped into the first bin is wrong for one sex in
  // each scheme.
  const SexCase cases[] = {
      {{"male", "female"}, {{1e4, 1e4}, {0, 0}}, Sex::MALE},
      {{"female", "male"}, {{0, 0}, {1e4, 1e4}}, Sex::MALE},
      {{"male", "female"}, {{0, 0}, {1e4, 1e4}}, Sex::FEMALE},
  };
  for (const SexCase& c : cases) {
    for (Sex sex : {Sex::MALE, Sex::FEMALE}) {
      const bool expect = sex == c.infectable_sex;
      const BinOutcome outcome = infectionByBin(sex, 30.0f, c.bins, c.contacts);
      if (!outcome.evaluated) continue;
      CHECK(outcome.local_infected == expect);
      CHECK(outcome.visitor_infected == expect);
    }
  }
}

TEST_CASE("H12b: Visitor is binned by age across band edges like a local") {
  const std::vector<std::string> bins = {"0-17", "18-64", "65+"};
  // Only the middle band has contacts, and the infectious local is in it.
  const std::vector<std::vector<double>> contacts = {
      {0, 0, 0}, {0, 1e4, 0}, {0, 0, 0}};
  struct AgeCase {
    float age;
    bool in_middle_band;
  };
  const AgeCase cases[] = {{0.0f, false},  {17.0f, false}, {17.99f, false},
                           {18.0f, true},  {64.0f, true},  {64.99f, true},
                           {65.0f, false}, {99.0f, false}, {100.0f, false},
                           {130.0f, false}};
  for (const AgeCase& c : cases) {
    CAPTURE(c.age);
    const BinOutcome outcome =
        infectionByBin(Sex::UNKNOWN, c.age, bins, contacts);
    if (!outcome.evaluated) continue;
    CHECK(outcome.local_infected == c.in_middle_band);
    CHECK(outcome.visitor_infected == c.in_middle_band);
  }
}

// One slot at venue 1, where nobody is infectious but a large fomite deposit
// lies. Person 0 visits it; person 1, its resident, is there too if
// `resident_present`. Returns the infections this rank applied from the
// cross-rank exchange.
static std::vector<PendingInfection> slotAtFomiteVenue(TwoRankFixture& f,
                                                       Disease& disease,
                                                       bool resident_present) {
  f.dm->exchangeVisitors({makeRemoteLocation(f.rank)}, disease, 0.3, 1.0);
  std::vector<PendingInfection> pending;
  const auto& incoming = f.dm->getDomain().incoming_visitors;
  // Not REQUIRE on rank 1: rank 0 would then wait forever on the exchange.
  if (f.rank == 1 && incoming.size() == 1) {
    const PersonId visitor_id = incoming[0].person_id;
    std::unordered_map<PersonId, VisitorInfo> visitor_data = {
        {visitor_id, toVisitorInfo(incoming[0])}};
    std::unordered_set<PersonId> visitor_ids = {visitor_id};
    std::vector<PersonLocation> locs = {{visitor_id, f.rank, -1, 0, 255, 0}};
    if (resident_present) locs.push_back({f.rank, f.rank, -1, 0, 255, 0});

    Venue* venue = f.world.getVenue(f.rank);
    venue->fomite_history.assign(1, {});
    venue->fomite_history[0].push_back({0.2, 1.0e6});
    ContactMatrixConfig cm;
    ContactMatrix default_contact_matrix;
    default_contact_matrix.bins = {"all"};
    default_contact_matrix.contacts = {{100.0}};
    cm.default_matrix = default_contact_matrix;
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, f.world, disease);
    SimulationConfig sim;
    ParallelConfig par;
    InteractionManager im(f.world, cm, sim, par, &disease, nullptr);
    im.processTransmissions(locs, 0.3, 1.0, nullptr, &visitor_ids, &pending,
                            &visitor_data);
    CHECK(pending.size() == 1);
  }
  if (f.rank == 1) CHECK(incoming.size() == 1);
  return f.dm->receivePendingInfections(pending, disease);
}

// ---------------------------------------------------------------------------
// H10: a fomite infection of a Visitor keeps its source on the home rank
// ---------------------------------------------------------------------------
TEST_CASE("H10: Cross-rank fomite infection is applied as Fomite") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);
  Disease disease = makeFomiteDisease(/*sub_bin_time=*/0.0);

  auto applied = slotAtFomiteVenue(f, disease, /*resident_present=*/false);
  if (f.rank == 0) {
    REQUIRE(applied.size() == 1);
    CHECK(applied[0].source == InfectionSource::Fomite);
    CHECK(applied[0].transmission_mode_index == 1);
  }
}

// ---------------------------------------------------------------------------
// H11: the outcome row a fomite infection takes does not depend on whether
// the infected person lives on the rank that holds the Venue
// ---------------------------------------------------------------------------
TEST_CASE("H11: Local and cross-rank fomite infections take the same row") {
  TwoRankFixture f;
  REQUIRE(f.size == 2);

  SelectionCriterion source_is_fomite;
  source_is_fomite.property_path = "infection_source";
  source_is_fomite.operator_type = "==";
  source_is_fomite.value = "fomite";
  OutcomeRates outcomes;
  outcomes.rows = {{{source_is_fomite}, {{"fomite_route", 1.0}}},
                   {{}, {{"default_route", 1.0}}}};
  TrajectoryDefinition fomite_route;
  fomite_route.selection_key = "fomite_route";
  fomite_route.stages.push_back({"exposed", {"constant", {{"value", 100.0}}}});
  TrajectoryDefinition default_route;
  default_route.selection_key = "default_route";
  default_route.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  Disease disease = makeFomiteDisease(/*sub_bin_time=*/0.0, outcomes,
                                      {fomite_route, default_route});
  disease.resolve(f.world);

  // Rank 0's person is infected as a Visitor, rank 1's at home.
  slotAtFomiteVenue(f, disease, /*resident_present=*/true);
  const Person* person = f.world.getPerson(f.rank);
  REQUIRE(person->infection != nullptr);
  CHECK(disease.getSymptomName(
            person->infection->getTrajectory().transitions.at(0).second) ==
        "exposed");
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
