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

#include "mpi_test_helpers.h"

using namespace june;

namespace {

constexpr int kNumRanks = 4;

using fomite_flu::kHealthy;
using fomite_flu::kMild;

// Direct mode 0 plus fomite mode 1 in 2 h sub-bins.
constexpr double kSubBinTime = 2.0;
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
  Disease disease = makeFomiteDisease(kSubBinTime);

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
        Infection::fromCheckpoint(&disease, 9.0, trajectory, 1.0, 1.0, 1.0, 0.0,
                                  /*last_checked_time=*/-1.0, kHealthy, 9.0);
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
