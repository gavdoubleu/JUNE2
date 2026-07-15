#pragma once

#ifdef USE_MPI

#include <mpi.h>

#include <map>
#include <optional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "activity/coordinated_encounter_types.h"
#include "activity/runtime_bin_allocator.h"
#include "core/config.h"
#include "core/world_state.h"
#include "domain.h"
#include "epidemiology/disease.h"

namespace june {
class DomainManager;

/**
 * DomainCommunicator - Handles adaptive MPI communication for visitors and
 * infections
 */
class DomainCommunicator {
 public:
  DomainCommunicator(WorldState& world, const Config& config, Domain& domain);

  void setDisease(const Disease* disease) { disease_ = disease; }

  // Exchange visitors across ranks
  void exchangeVisitors(const std::vector<PersonLocation>& locations,
                        const DomainManager& dm, double current_time,
                        double delta_hours = 0.0,
                        const RuntimeBinAllocator* alloc = nullptr);

  // Exchange pending infections across ranks. Returns the PendingInfection
  // records that this rank applied (i.e. for which a new local infection was
  // created on the home rank), so the caller can log InfectionEvents on the
  // owning rank.
  std::vector<PendingInfection> receivePendingInfections(
      const std::vector<PendingInfection>& pending_infections);

  // Cross-rank coordinated encounter exchange
  void exchangeEncounterProposals(
      const std::vector<EncounterProposal>& local_proposals,
      const DomainManager& dm,
      std::vector<EncounterProposal>& proposals_for_this_rank);
  void exchangeEncounterReplies(
      const std::vector<EncounterReply>& local_replies, const DomainManager& dm,
      std::vector<EncounterReply>& replies_for_this_rank);
  void exchangeFinalizedEncounters(
      const std::vector<CoordinatedEncounter>& local_finalized,
      const DomainManager& dm,
      std::vector<CoordinatedEncounter>& finalized_for_this_rank);

 private:
  void exchangeAllToAll(
      const std::vector<std::vector<Domain::VisitorData>>& outgoing,
      const std::vector<int>& send_counts);
  void exchangePointToPoint(
      const std::vector<std::vector<Domain::VisitorData>>& outgoing,
      const std::vector<int>& send_counts);
  void exchangePointToPointVerySparse(
      const std::vector<std::vector<Domain::VisitorData>>& outgoing,
      const std::vector<int>& send_counts);

  // Shared P2P logic for visitor exchange
  void performP2PVisitorExchange(
      const std::vector<std::vector<Domain::VisitorData>>& outgoing,
      const std::vector<int>& send_counts, const std::vector<int>& recv_counts);

  // Sparsity-driven dispatch: chooses among all-to-all, P2P, or very-sparse
  // P2P based on the cross-rank send-pair density (MPI_Allreduce), then
  // records the outgoing visitors locally.
  void dispatchVisitorExchange(
      const std::vector<std::vector<Domain::VisitorData>>& outgoing,
      const std::vector<int>& send_counts);

  // Join the local PendingInfections against incoming_visitors to find the
  // home rank that needs to be notified for each. Returns the per-rank
  // updates table and corresponding send_counts vector.
  void routePendingByHomeRank(
      const std::vector<PendingInfection>& pending,
      std::vector<std::vector<PendingInfection>>& updates,
      std::vector<int>& send_counts);

  // Walks the receive buffer of an Alltoallv exchange of PendingInfections,
  // unpacks each record, and delegates to applyOnePendingInfection. Returns
  // the vector of newly-applied records.
  std::vector<PendingInfection> unpackAndApplyIncoming(
      const std::vector<char>& rbuf, const std::vector<int>& rd,
      const std::vector<int>& recv_counts);

  // Construct the local Infection for one successfully-unpacked pending
  // record and return the resulting PendingInfection. Returns std::nullopt
  // if the record was skipped (person not owned, already infected, or no
  // disease loaded).
  std::optional<PendingInfection> applyOnePendingInfection(
      const PendingInfection& pending);

  WorldState& world_;
  const Config& config_;
  Domain& domain_;
  const Disease* disease_;
  int rank_, num_ranks_;
};

}  // namespace june

#endif  // USE_MPI
