#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#ifdef USE_MPI
#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <optional>
#include <vector>

#include "epidemiology/vaccination_allocation.h"
#include "utils/deterministic_rng.h"
#include "utils/mpi_utils.h"

using namespace june;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();
  MPI_Finalize();
  return result;
}

namespace {

// Every rank builds the same key set from a fixed seed, which stands for the
// whole world, and keeps only its own slice of it.
std::vector<uint64_t> worldKeys(size_t n) {
  SplitMix64 g(314);
  std::vector<uint64_t> keys(n);
  for (auto& k : keys) k = g();
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  return keys;
}

enum class Slicing { Interleaved, ByKeyRange };

std::vector<uint64_t> mySlice(const std::vector<uint64_t>& sorted_world,
                              Slicing how, int rank, int size) {
  std::vector<uint64_t> mine;
  for (size_t i = 0; i < sorted_world.size(); ++i) {
    const size_t owner =
        how == Slicing::Interleaved ? i % size : i * size / sorted_world.size();
    if (static_cast<int>(owner) == rank) mine.push_back(sorted_world[i]);
  }
  return mine;  // still ascending
}

}  // namespace

// Each rank holds a different slice of one key set. The cutoff found through
// real MPI reduces must equal the one a single process finds from the union,
// and the people served must be the same people.
TEST_CASE("lotteryCutoff - real reduces give the serial selection") {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const auto world = worldKeys(5000);
  const int64_t total = static_cast<int64_t>(world.size());

  for (Slicing how : {Slicing::Interleaved, Slicing::ByKeyRange}) {
    const auto mine = mySlice(world, how, rank, size);
    // The slices together are the world and are not all the same size.
    CHECK(mpi_runtime::globalSum(static_cast<int64_t>(mine.size())) == total);

    for (int64_t take : {int64_t{0}, int64_t{1}, int64_t{37}, total / 3,
                         total / 2, total - 1, total}) {
      CAPTURE(take);
      const auto cutoff = lotteryCutoff(
          mpi_runtime::globalSum(static_cast<int64_t>(mine.size())), take,
          [&](uint64_t t) {
            return mpi_runtime::globalSum(static_cast<int64_t>(
                std::upper_bound(mine.begin(), mine.end(), t) - mine.begin()));
          });

      // The serial answer, from the union.
      if (take == 0) {
        CHECK_FALSE(cutoff.has_value());
        continue;
      }
      REQUIRE(cutoff.has_value());
      CHECK(*cutoff == world[take - 1]);

      // This rank serves exactly its keys inside the serial selection.
      int64_t served = 0;
      for (uint64_t k : mine) {
        CHECK((k <= *cutoff) == (k <= world[take - 1]));
        if (k <= *cutoff) ++served;
      }
      CHECK(mpi_runtime::globalSum(served) == take);

      // Every rank agrees on the cutoff.
      uint64_t lowest = 0;
      uint64_t highest = 0;
      MPI_Allreduce(&*cutoff, &lowest, 1, MPI_UINT64_T, MPI_MIN,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&*cutoff, &highest, 1, MPI_UINT64_T, MPI_MAX,
                    MPI_COMM_WORLD);
      CHECK(lowest == highest);
    }
  }
}

// The reduce the engine will use sums plain integers over the ranks.
TEST_CASE("globalSum - adds one integer over all ranks") {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  CHECK(mpi_runtime::globalSum(rank + 1) ==
        static_cast<int64_t>(size) * (size + 1) / 2);
  CHECK(mpi_runtime::globalSum(int64_t{1} << 40) ==
        static_cast<int64_t>(size) * (int64_t{1} << 40));
}

#else
TEST_CASE("MPI vaccination allocation requires MPI") {}
#endif
