#ifdef USE_MPI

#include "parallel/domain_communicator.h"

#include <cstddef>
#include <cstring>
#include <iostream>
#include <type_traits>
#include <vector>

#include "parallel/domain_communicator_detail.h"
#include "parallel/domain_manager.h"
#include "parallel/mpi_utils.h"
#include "utils/profiler.h"

namespace {

using june::domain_comm_detail::makeWireRecord;
using june::domain_comm_detail::packField;
using june::domain_comm_detail::unpackField;

// Fixed header of the visitor wire format (everything before the
// integrated_infectiousness payload); the variable-length per-mode payload
// is appended manually after this, outside WireRecord.
// Total wire size = VISITOR_WIRE_HEADER + num_modes * sizeof(double).
constexpr auto kVisitorWire = makeWireRecord(
    &june::Domain::VisitorData::person_id,
    &june::Domain::VisitorData::home_rank, &june::Domain::VisitorData::venue_id,
    &june::Domain::VisitorData::subset_idx,
    &june::Domain::VisitorData::is_infected,
    &june::Domain::VisitorData::is_infectious,
    &june::Domain::VisitorData::immunity_level,
    &june::Domain::VisitorData::encounter_type_id,
    &june::Domain::VisitorData::symptom_id,
    &june::Domain::VisitorData::time_in_stage);
constexpr int VISITOR_WIRE_HEADER = kVisitorWire.size();
// Tripwire: VisitorData's trailing integrated_infectiousness is a
// std::vector<double> packed manually as a count-known-elsewhere tail (not
// via WireRecord), and fields after it (newly_infected etc.) are pure
// return data never on the wire, so sizeof(VisitorData) isn't a useful
// proxy here. offsetof(integrated_infectiousness) instead marks where the
// fixed header covered by kVisitorWire ends - it moves if a field is
// added/removed/resized anywhere before the tail.
// offsetof is only standard-guaranteed for standard-layout types; guard that
// assumption explicitly so a future member (e.g. a std::set) that breaks it
// fails loudly here rather than degrading to a silent -Winvalid-offsetof.
static_assert(std::is_standard_layout_v<june::Domain::VisitorData>,
              "VisitorData must stay standard-layout for the offsetof check "
              "below to be well-defined");
static_assert(offsetof(june::Domain::VisitorData, integrated_infectiousness) ==
                  40,
              "VisitorData's fixed-header region changed - check kVisitorWire "
              "covers every field, then update this literal");
inline int visitorWireSize(int num_modes) {
  return VISITOR_WIRE_HEADER + num_modes * static_cast<int>(sizeof(double));
}

char* packVisitor(char* ptr, const june::Domain::VisitorData& v,
                  int num_modes) {
  ptr = kVisitorWire.pack(ptr, v);
  // Variable-length per-mode payload. The vector must be exactly num_modes
  // long; sender and receiver agree on this from disease->numModes(), which
  // is loaded identically on every rank. Mismatches indicate a config bug
  // and are caught loud at the build site rather than papered over here.
  if (static_cast<int>(v.integrated_infectiousness.size()) != num_modes) {
    throw std::runtime_error(
        "packVisitor: integrated_infectiousness size " +
        std::to_string(v.integrated_infectiousness.size()) + " != num_modes " +
        std::to_string(num_modes));
  }
  if (num_modes > 0) {
    std::memcpy(ptr, v.integrated_infectiousness.data(),
                num_modes * sizeof(double));
    ptr += num_modes * sizeof(double);
  }
  return ptr;
}

const char* unpackVisitor(const char* ptr, june::Domain::VisitorData& v,
                          int num_modes) {
  ptr = kVisitorWire.unpack(ptr, v);
  v.integrated_infectiousness.assign(num_modes, 0.0);
  if (num_modes > 0) {
    std::memcpy(v.integrated_infectiousness.data(), ptr,
                num_modes * sizeof(double));
    ptr += num_modes * sizeof(double);
  }
  return ptr;
}

// Builds a fully-populated VisitorData for a person attending a remote
// venue. Pre-computes integrated_infectiousness per mode using the SAME
// code path as local people (getIntegratedInfectiousness) so FP results
// are bit-identical regardless of where a person is processed. The
// integrated_infectiousness vector is always sized to num_modes
// (zero-filled for non-infectious people); the wire format expects
// exactly that many doubles.
june::Domain::VisitorData buildVisitorPayload(const june::PersonLocation& loc,
                                              const june::Person& person,
                                              int home_rank,
                                              double current_time,
                                              double delta_hours, int num_modes,
                                              const june::Disease* disease) {
  june::Domain::VisitorData visitor;
  visitor.person_id = loc.person_id;
  visitor.home_rank = home_rank;
  visitor.venue_id = loc.venue_id;
  visitor.subset_idx = loc.subset_index;
  visitor.is_infected = (person.infection != nullptr);
  visitor.is_infectious =
      visitor.is_infected && person.infection->isInfectious(current_time);

  double susceptibility = 1.0;
  if (disease) {
    susceptibility = person.getSusceptibility(current_time, disease->getName());
  } else {
    susceptibility = 1.0 - person.immunity.natural_level;
  }
  visitor.immunity_level = static_cast<float>(1.0 - susceptibility);

  visitor.encounter_type_id = loc.encounter_type_id;
  visitor.newly_infected = false;
  visitor.new_infection_time = -1.0;

  visitor.symptom_id = 0;
  visitor.time_in_stage = 0.0;
  visitor.integrated_infectiousness.assign(num_modes, 0.0);
  if (visitor.is_infected && person.infection) {
    const june::InfectionTrajectory& traj = person.infection->getTrajectory();
    double stage_start_time = traj.infection_time;
    uint16_t cur_symptom_id = 0;
    for (const auto& trans : traj.transitions) {
      if (current_time >= trans.first) {
        stage_start_time = trans.first;
        cur_symptom_id = trans.second;
      } else {
        break;
      }
    }
    visitor.symptom_id = cur_symptom_id;
    visitor.time_in_stage = current_time - stage_start_time;

    if (visitor.is_infectious && disease) {
      double t1 = current_time + delta_hours / 24.0;
      for (int m = 0; m < num_modes; ++m) {
        visitor.integrated_infectiousness[m] =
            person.infection->getIntegratedInfectiousness(m, current_time, t1);
      }
    }
  }
  return visitor;
}

}  // anonymous namespace

namespace june {

DomainCommunicator::DomainCommunicator(WorldState& world, const Config& config,
                                       Domain& domain)
    : world_(world), config_(config), domain_(domain), disease_(nullptr) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks_);
}

void DomainCommunicator::exchangeVisitors(
    const std::vector<PersonLocation>& locations, const DomainManager& dm,
    double current_time, double delta_hours, const RuntimeBinAllocator* alloc) {
  domain_.clearVisitors();
  std::vector<std::vector<Domain::VisitorData>> outgoing(num_ranks_);
  std::vector<int> send_counts(num_ranks_, 0);
  // Mode count is fixed for the run (Disease YAML defines it once).
  // Hoist here so the same value sizes every visitor's payload below.
  // Wire size for MPI buffers is computed inside the helpers (Pack/Unpack
  // and exchangeAllToAll/performP2PVisitorExchange) since they each own
  // their own buffer layout.
  const int num_modes =
      (disease_ && disease_->numModes() > 0) ? disease_->numModes() : 1;

  auto send = [&](const PersonLocation& loc, Person& person, int target_rank) {
    outgoing[target_rank].push_back(buildVisitorPayload(
        loc, person, rank_, current_time, delta_hours, num_modes, disease_));
    send_counts[target_rank]++;
  };

  for (const auto& loc : locations) {
    if (loc.venue_id == -1) continue;
    if (!domain_.ownsPerson(loc.person_id)) continue;
    // A rider's lines are shipped below, one per leg. Their location names
    // only one of those legs, so sending on it would leave the other legs
    // short of a passenger they are carrying.
    if (alloc && alloc->isPartialPresenceVenue(loc.venue_id)) continue;
    if (domain_.ownsVenue(loc.venue_id)) continue;

    int target_rank = dm.getVenueRank(loc.venue_id);
    if (target_rank == -1) continue;
    if (target_rank == rank_) {
      // Venue is on this rank (e.g., virtual venue hosted locally), so no
      // need to send as visitor; interaction will be computed locally.
      continue;
    }

    Person* person = world_.getPerson(loc.person_id);
    if (!person) continue;

    send(loc, *person, target_rank);
  }

  // One visitor record per (rider, leg): a commuter changing trains three
  // times is carried by three lines, and each line's owner needs their disease
  // state to work out who they infected on board.
  if (alloc) {
    for (const auto& [vid, riders] : alloc->ridersByVenue()) {
      if (domain_.ownsVenue(vid)) continue;
      const int target_rank = dm.getVenueRank(vid);
      if (target_rank == -1 || target_rank == rank_) continue;

      for (const auto& r : riders) {
        if (!domain_.ownsPerson(r.pid)) continue;
        Person* person = world_.getPerson(r.pid);
        if (!person) continue;
        auto it = world_.person_index.find(r.pid);
        if (it == world_.person_index.end()) continue;

        PersonLocation leg = locations[it->second];
        leg.venue_id = vid;
        leg.subset_index = r.subset;
        send(leg, *person, target_rank);
      }
    }
  }

  dispatchVisitorExchange(outgoing, send_counts);
}

void DomainCommunicator::dispatchVisitorExchange(
    const std::vector<std::vector<Domain::VisitorData>>& outgoing,
    const std::vector<int>& send_counts) {
  int local_pairs = 0;
  for (int c : send_counts)
    if (c > 0) local_pairs++;
  int global_pairs;

  MPI_Allreduce(&local_pairs, &global_pairs, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  double sparsity = 100.0 * global_pairs / (num_ranks_ * (num_ranks_ - 1));

  if (sparsity >= 70.0)
    exchangeAllToAll(outgoing, send_counts);
  else if (sparsity >= 15.0)
    exchangePointToPoint(outgoing, send_counts);
  else
    exchangePointToPointVerySparse(outgoing, send_counts);

  for (int r = 0; r < num_ranks_; ++r) {
    for (const auto& v : outgoing[r]) domain_.addOutgoingVisitor(v);
  }
}

void DomainCommunicator::exchangeAllToAll(
    const std::vector<std::vector<Domain::VisitorData>>& outgoing,
    const std::vector<int>& send_counts_in) {
  std::vector<int> send_counts = send_counts_in;
  std::vector<int> recv_counts(num_ranks_, 0);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
               MPI_COMM_WORLD);

  const int num_modes =
      (disease_ && disease_->numModes() > 0) ? disease_->numModes() : 1;
  const int wire_size = visitorWireSize(num_modes);

  std::vector<int> sdisp, rdisp;
  int stotal, rtotal;
  mpi_utils::computeByteDisplacements(send_counts, wire_size, sdisp, stotal);
  mpi_utils::computeByteDisplacements(recv_counts, wire_size, rdisp, rtotal);

  std::vector<char> sbuf(stotal);
  std::vector<char> rbuf(rtotal);

  for (int r = 0; r < num_ranks_; ++r) {
    char* ptr = sbuf.data() + sdisp[r];
    for (const auto& v : outgoing[r]) {
      ptr = packVisitor(ptr, v, num_modes);
    }
  }

  std::vector<int> sc(num_ranks_), rc(num_ranks_);
  for (int i = 0; i < num_ranks_; ++i) {
    sc[i] = send_counts[i] * wire_size;
    rc[i] = recv_counts[i] * wire_size;
  }

  MPI_Alltoallv(sbuf.data(), sc.data(), sdisp.data(), MPI_BYTE, rbuf.data(),
                rc.data(), rdisp.data(), MPI_BYTE, MPI_COMM_WORLD);

  for (int r = 0; r < num_ranks_; ++r) {
    const char* ptr = rbuf.data() + rdisp[r];
    for (int i = 0; i < recv_counts[r]; ++i) {
      Domain::VisitorData v;
      ptr = unpackVisitor(ptr, v, num_modes);
      if (domain_.ownsVenue(v.venue_id)) domain_.addIncomingVisitor(v);
    }
  }
}

void DomainCommunicator::exchangePointToPoint(
    const std::vector<std::vector<Domain::VisitorData>>& outgoing,
    const std::vector<int>& send_counts) {
  std::vector<int> recv_counts(num_ranks_, 0);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
               MPI_COMM_WORLD);

  performP2PVisitorExchange(outgoing, send_counts, recv_counts);
}

void DomainCommunicator::performP2PVisitorExchange(
    const std::vector<std::vector<Domain::VisitorData>>& outgoing,
    const std::vector<int>& send_counts, const std::vector<int>& recv_counts) {
  std::vector<std::vector<char>> sbufs(num_ranks_);
  std::vector<MPI_Request> sreqs, rreqs;
  std::vector<std::vector<char>> rbufs(num_ranks_);

  const int num_modes =
      (disease_ && disease_->numModes() > 0) ? disease_->numModes() : 1;
  const int wire_size = visitorWireSize(num_modes);

  for (int r = 0; r < num_ranks_; ++r) {
    if (r != rank_ && recv_counts[r] > 0) {
      try {
        rbufs[r].resize(recv_counts[r] * wire_size);
      } catch (const std::exception& e) {
        std::cerr << "[MPI] rbufs resize failed: recv_counts[r]="
                  << recv_counts[r] << " error: " << e.what() << std::endl;
        throw;
      }
      MPI_Request req;
      MPI_Irecv(rbufs[r].data(), rbufs[r].size(), MPI_BYTE, r, 101,
                MPI_COMM_WORLD, &req);
      rreqs.push_back(req);
    }
  }

  for (int r = 0; r < num_ranks_; ++r) {
    if (r != rank_ && send_counts[r] > 0) {
      try {
        sbufs[r].resize(send_counts[r] * wire_size);
      } catch (const std::exception& e) {
        std::cerr << "[MPI] sbufs resize failed: send_counts[r]="
                  << send_counts[r] << " error: " << e.what() << std::endl;
        throw;
      }
      char* ptr = sbufs[r].data();
      for (const auto& v : outgoing[r]) {
        ptr = packVisitor(ptr, v, num_modes);
      }
      MPI_Request req;
      MPI_Isend(sbufs[r].data(), sbufs[r].size(), MPI_BYTE, r, 101,
                MPI_COMM_WORLD, &req);
      sreqs.push_back(req);
    }
  }

  if (!rreqs.empty())
    MPI_Waitall(rreqs.size(), rreqs.data(), MPI_STATUSES_IGNORE);

  for (int r = 0; r < num_ranks_; ++r) {
    if (r != rank_ && recv_counts[r] > 0) {
      const char* ptr = rbufs[r].data();
      for (int i = 0; i < recv_counts[r]; ++i) {
        Domain::VisitorData v;
        ptr = unpackVisitor(ptr, v, num_modes);
        if (domain_.ownsVenue(v.venue_id)) domain_.addIncomingVisitor(v);
      }
    }
  }

  if (!sreqs.empty())
    MPI_Waitall(sreqs.size(), sreqs.data(), MPI_STATUSES_IGNORE);
}

void DomainCommunicator::exchangePointToPointVerySparse(
    const std::vector<std::vector<Domain::VisitorData>>& outgoing,
    const std::vector<int>& send_counts) {
  std::vector<int> pattern;
  for (int r = 0; r < num_ranks_; ++r) {
    if (send_counts[r] > 0) {
      pattern.push_back(r);
      pattern.push_back(send_counts[r]);
    }
  }

  int sz = static_cast<int>(pattern.size());
  std::vector<int> all_sz(num_ranks_);
  MPI_Allgather(&sz, 1, MPI_INT, all_sz.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> disp(num_ranks_, 0);
  int tot = 0;
  for (int r = 0; r < num_ranks_; ++r) {
    disp[r] = tot;
    tot += all_sz[r];
  }
  std::vector<int> all_pat(tot);
  MPI_Allgatherv(pattern.data(), sz, MPI_INT, all_pat.data(), all_sz.data(),
                 disp.data(), MPI_INT, MPI_COMM_WORLD);

  std::vector<int> recv_counts(num_ranks_, 0);
  for (int src = 0; src < num_ranks_; ++src) {
    for (int i = 0; i < all_sz[src] / 2; ++i) {
      if (all_pat[disp[src] + i * 2] == rank_) {
        recv_counts[src] = all_pat[disp[src] + i * 2 + 1];
      }
    }
  }

  // Call the shared point-to-point logic with reconstructed receive counts
  performP2PVisitorExchange(outgoing, send_counts, recv_counts);
}

}  // namespace june

#endif  // USE_MPI
