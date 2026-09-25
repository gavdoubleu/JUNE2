// MPI tests for the three visitor exchange paths (all-to-all, P2P, very-sparse
// P2P). DomainCommunicator picks the path from the share of rank pairs that
// send; each test picks its send pattern to land on one path.
//
// Must be run with exactly 4 MPI ranks (12 ordered pairs, so one sending pair
// is under the very-sparse threshold):
//   mpirun -np 4 ./test_mpi_visitor_exchange_paths

#define DOCTEST_CONFIG_IMPLEMENT  // custom main so we can wrap MPI
                                  // init/finalize
#include "doctest.h"

#ifdef USE_MPI
#include <mpi.h>

#include <algorithm>
#include <map>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/disease.h"
#include "epidemiology/infectiousness_curves.h"
#include "parallel/domain.h"
#include "parallel/domain_manager.h"

using namespace june;

namespace {

constexpr int kNumRanks = 4;

// Symptom ids: 0 = healthy, 1 = mild.
constexpr uint16_t kHealthy = 0;
constexpr uint16_t kMild = 1;

// Rank r owns geo unit r, venue r and person r.
struct RankPerPersonFixture {
  int rank;
  int size;

  WorldState world;
  Config config;
  std::unique_ptr<DomainManager> dm;

  RankPerPersonFixture() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    world.venue_type_names = {"household"};
    world.geo_level_names = {"MGU"};
    world.activity_names = {"residence", "work", "visiting", "none", "dead"};

    GeographicalUnit geo_unit;
    geo_unit.id = rank;
    geo_unit.name = "MGU_" + std::to_string(rank);
    geo_unit.level_id = 0;
    geo_unit.parent_id = -1;
    world.geo_units.push_back(geo_unit);

    Venue venue;
    venue.id = rank;
    venue.type_id = 0;
    venue.geo_unit_id = rank;
    venue.is_residence = true;
    world.venues.push_back(venue);

    Person& person = world.people.emplace_back();
    person.id = rank;
    person.age = 30.0f;
    person.sex = Sex::MALE;
    person.geo_unit_id = rank;

    world.buildIndices();

    config.parallel.partition_level = "MGU";
    config.parallel.geo_unit_chunk_size = 1000;

    dm = std::make_unique<DomainManager>(world, config);
    dm->setMPI(rank, size);
    dm->setMaxPersonId(size - 1);
    for (int owner = 0; owner < size; ++owner) {
      dm->setGeoUnitRank(owner, owner);
      dm->setPersonRank(owner, owner);
      dm->setVenueRank(owner, owner);
    }

    Domain& domain = dm->getDomain();
    domain.addGeoUnit(rank);
    domain.resident_ids.push_back(rank);
    domain.resident_set.insert(rank);
    domain.local_venue_ids.push_back(rank);
    domain.local_venue_set.insert(rank);
  }
};

// Direct mode 0 plus fomite mode 1 in 2 h sub-bins.
Disease makeFomiteDisease() {
  TransmissionParams params;
  params.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(5.0);
  params.symptom_id_curves = {nullptr, curve};

  TransmissionMode direct;
  direct.name = "direct";
  direct.symptom_curves = params.symptom_id_curves;
  params.modes.push_back(std::move(direct));

  TransmissionMode fomite;
  fomite.name = "fomite";
  fomite.type = TransmissionModeType::Fomite;
  fomite.symptom_curves = {nullptr, nullptr};
  FomiteConfig fomite_config;
  fomite_config.mode_index = 1;
  fomite_config.max_age = 2.0;
  fomite_config.sub_bin_time = 2.0;
  fomite_config.infectiousness_curve = std::make_shared<ConstantCurve>(1.0);
  fomite_config.deposition_by_symptom = {std::make_shared<ConstantCurve>(0.7),
                                         std::make_shared<ConstantCurve>(2.0)};
  fomite.config = std::move(fomite_config);
  params.modes.push_back(std::move(fomite));

  std::vector<SymptomTag> symptom_tags = {{"healthy", -1, kHealthy},
                                          {"mild", 1, kMild}};
  TrajectoryDefinition trajectory_definition;
  trajectory_definition.selection_key = "general";
  trajectory_definition.severity = 1.0;
  trajectory_definition.stages.push_back(
      {"mild", {"constant", {{"value", 100.0}}}});
  return Disease("FomiteFlu", symptom_tags, {}, {trajectory_definition}, {},
                 params);
}

constexpr int kNumModes = 2;
constexpr double kCurrentTime = 10.0;
constexpr double kDeltaHours = 6.0;
constexpr int kFomiteSubBins = 3;  // 6 h slot / 2 h sub-bins

// Each rank's person is in one disease state, so every exchange mixes
// record lengths: infectious (both tails), incubating (deposits only) and
// uninfected (header only).
enum class HomeState { Uninfected, Incubating, Infectious };
HomeState homeState(int home_rank) {
  switch (home_rank % 3) {
    case 0:
      return HomeState::Infectious;
    case 1:
      return HomeState::Uninfected;
    default:
      return HomeState::Incubating;
  }
}

// Every (source, destination) pair in `pairs` sends person `source` to venue
// `destination` for one exchange. Each rank then checks it received exactly
// its expected visitors, intact, whatever mix of record lengths it gets.
void checkExchange(const std::vector<std::pair<int, int>>& pairs) {
  RankPerPersonFixture fixture;
  REQUIRE(fixture.size == kNumRanks);
  Disease disease = makeFomiteDisease();
  fixture.dm->setDisease(&disease);

  const HomeState state = homeState(fixture.rank);
  if (state != HomeState::Uninfected) {
    InfectionTrajectory trajectory;
    trajectory.infection_time = 9.0;
    // Healthy has no direct-contact curve, so incubating isn't infectious
    // but still deposits.
    trajectory.transitions =
        state == HomeState::Infectious
            ? std::vector<std::pair<double, uint16_t>>{{9.0, kMild}}
            : std::vector<std::pair<double, uint16_t>>{{9.0, kHealthy},
                                                       {11.0, kMild}};
    fixture.world.getPerson(fixture.rank)->infection =
        Infection::fromCheckpoint(&disease, 9.0, trajectory, 1.0, 1.0, 1.0,
                                  0.0, /*last_checked_time=*/-1.0, kHealthy,
                                  9.0);
  }

  std::vector<PersonLocation> locations;
  std::vector<int> expected_senders;
  for (const auto& [source, destination] : pairs) {
    if (source == fixture.rank) {
      PersonLocation location;
      location.person_id = source;
      location.venue_id = destination;
      location.subset_index = 0;
      location.activity_index = 1;
      location.encounter_type_id = 255;
      locations.push_back(location);
    }
    if (destination == fixture.rank) expected_senders.push_back(source);
  }

  fixture.dm->exchangeVisitors(locations, disease, kCurrentTime, kDeltaHours);

  const auto& incoming = fixture.dm->getDomain().incoming_visitors;
  std::vector<int> senders;
  for (const auto& visitor : incoming) {
    senders.push_back(visitor.home_rank);
    CHECK(visitor.person_id == visitor.home_rank);
    CHECK(visitor.venue_id == fixture.rank);
    const HomeState sender_state = homeState(visitor.home_rank);
    const bool infected = sender_state != HomeState::Uninfected;
    const bool infectious = sender_state == HomeState::Infectious;
    CHECK(visitor.is_infected == infected);
    CHECK(visitor.is_infectious == infectious);
    // A tail arrives full length if its header gate sends it, else empty.
    REQUIRE(visitor.integrated_infectiousness.size() ==
            (infectious ? kNumModes : 0));
    REQUIRE(visitor.fomite_deposition_sub.size() ==
            (infected ? kFomiteSubBins : 0));
    if (infectious) CHECK(visitor.integrated_infectiousness[0] > 0.0);
    for (double deposit : visitor.fomite_deposition_sub) CHECK(deposit > 0.0);
  }
  std::sort(senders.begin(), senders.end());
  std::sort(expected_senders.begin(), expected_senders.end());
  CHECK(senders == expected_senders);
}

}  // namespace

TEST_CASE("Very-sparse P2P: one sending pair of 12") {
  checkExchange({{0, 3}});
}

TEST_CASE("P2P: ring, 4 sending pairs of 12") {
  checkExchange({{0, 1}, {1, 2}, {2, 3}, {3, 0}});
}

TEST_CASE("All-to-all: every pair sends") {
  std::vector<std::pair<int, int>> pairs;
  for (int source = 0; source < kNumRanks; ++source)
    for (int destination = 0; destination < kNumRanks; ++destination)
      if (source != destination) pairs.emplace_back(source, destination);
  checkExchange(pairs);
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
