#pragma once

#include <cstdint>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace june::mpi_runtime {

struct State {
  int rank = 0;
  int size = 1;
  bool active = false;
};

// Return the current communicator state, or serial defaults before MPI_Init,
// after MPI_Finalize, or in a serial build.
inline State state() {
  State result;
#ifdef USE_MPI
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (!initialized) return result;

  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized) return result;

  MPI_Comm_rank(MPI_COMM_WORLD, &result.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &result.size);
  result.active = true;
#endif
  return result;
}

// Sum of one integer over all ranks. Integer only, so every rank gets the same
// value whatever the process count. Serial and single rank runs return the
// value as it is. Collective: every rank must call it.
inline int64_t globalSum(int64_t local) {
#ifdef USE_MPI
  const State s = state();
  if (s.active && s.size > 1)
    MPI_Allreduce(MPI_IN_PLACE, &local, 1, MPI_INT64_T, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
  return local;
}

// Smallest of one integer over all ranks, with the same rules as globalSum.
inline int64_t globalMin(int64_t local) {
#ifdef USE_MPI
  const State s = state();
  if (s.active && s.size > 1)
    MPI_Allreduce(MPI_IN_PLACE, &local, 1, MPI_INT64_T, MPI_MIN,
                  MPI_COMM_WORLD);
#endif
  return local;
}

}  // namespace june::mpi_runtime
