#ifdef USE_MPI

#include "parallel/visitor_wire.h"

#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "parallel/domain_communicator_detail.h"

namespace june::visitor_wire {

namespace {

using june::domain_comm_detail::makeWireRecord;

// Fixed header of the visitor wire format (everything before the
// integrated_infectiousness payload); the variable-length per-mode payload and
// the per-fomite-sub-bin deposits are appended manually after this, outside
// WireRecord.
constexpr auto kVisitorWire = makeWireRecord(
    &Domain::VisitorData::person_id, &Domain::VisitorData::home_rank,
    &Domain::VisitorData::venue_id, &Domain::VisitorData::subset_idx,
    &Domain::VisitorData::is_infected, &Domain::VisitorData::is_infectious,
    &Domain::VisitorData::immunity_level,
    &Domain::VisitorData::encounter_type_id, &Domain::VisitorData::symptom_id);
constexpr int VISITOR_WIRE_HEADER = kVisitorWire.size();
// Tripwire: VisitorData's trailing integrated_infectiousness and
// fomite_deposition_sub are std::vector<double>s packed manually as
// count-known-elsewhere tails (not via WireRecord), and fields after them
// (newly_infected etc.) are pure return data never on the wire, so
// sizeof(VisitorData) isn't a useful proxy here.
// offsetof(integrated_infectiousness) instead marks where the fixed header
// covered by kVisitorWire ends - it moves if a field is added/removed/resized
// anywhere before the tails.
// offsetof is only standard-guaranteed for standard-layout types; guard that
// assumption explicitly so a future member (e.g. a std::set) that breaks it
// fails loudly here rather than degrading to a silent -Winvalid-offsetof.
static_assert(std::is_standard_layout_v<Domain::VisitorData>,
              "VisitorData must stay standard-layout for the offsetof check "
              "below to be well-defined");
static_assert(offsetof(Domain::VisitorData, integrated_infectiousness) == 32,
              "VisitorData's fixed-header region changed - check kVisitorWire "
              "covers every field, then update this literal");

// Tails are exactly `count` long; sender and receiver agree on the counts from
// the Disease and timestep, which every rank loads identically. A mismatch is
// a config bug, caught loud rather than papered over.
char* packTail(char* ptr, const std::vector<double>& tail, int count,
               const char* name) {
  if (static_cast<int>(tail.size()) != count) {
    throw std::runtime_error(std::string("visitor_wire::pack: ") + name + " size " +
                             std::to_string(tail.size()) + " != " +
                             std::to_string(count));
  }
  if (count > 0) {
    std::memcpy(ptr, tail.data(), count * sizeof(double));
    ptr += count * sizeof(double);
  }
  return ptr;
}

const char* unpackTail(const char* ptr, std::vector<double>& tail,
                       int count) {
  tail.assign(count, 0.0);
  if (count > 0) {
    std::memcpy(tail.data(), ptr, count * sizeof(double));
    ptr += count * sizeof(double);
  }
  return ptr;
}

}  // namespace

int recordSize(const Domain::VisitorData& /*visitor*/,
               const TailCounts& tails) {
  return VISITOR_WIRE_HEADER + (tails.num_modes + tails.fomite_sub_bins) *
                                   static_cast<int>(sizeof(double));
}

char* pack(char* ptr, const Domain::VisitorData& visitor,
           const TailCounts& tails) {
  ptr = kVisitorWire.pack(ptr, visitor);
  ptr = packTail(ptr, visitor.integrated_infectiousness, tails.num_modes,
                 "integrated_infectiousness");
  return packTail(ptr, visitor.fomite_deposition_sub, tails.fomite_sub_bins,
                  "fomite_deposition_sub");
}

const char* unpack(const char* ptr, Domain::VisitorData& visitor,
                   const TailCounts& tails) {
  ptr = kVisitorWire.unpack(ptr, visitor);
  ptr = unpackTail(ptr, visitor.integrated_infectiousness, tails.num_modes);
  return unpackTail(ptr, visitor.fomite_deposition_sub, tails.fomite_sub_bins);
}

namespace detail {

const char* unpackWithin(const char* ptr, const char* end,
                         Domain::VisitorData& visitor,
                         const TailCounts& tails) {
  // Header first: the record's full size may depend on it.
  if (end - ptr < VISITOR_WIRE_HEADER) {
    throw std::runtime_error(
        "visitor_wire::unpackSlice: " + std::to_string(end - ptr) +
        " bytes left, too few for a record header");
  }
  const char* tails_begin = kVisitorWire.unpack(ptr, visitor);
  const int size = recordSize(visitor, tails);
  if (end - ptr < size) {
    throw std::runtime_error("visitor_wire::unpackSlice: record of " +
                             std::to_string(size) + " bytes runs past slice "
                             "end (" + std::to_string(end - ptr) + " left)");
  }
  const char* next =
      unpackTail(tails_begin, visitor.integrated_infectiousness,
                 tails.num_modes);
  return unpackTail(next, visitor.fomite_deposition_sub,
                    tails.fomite_sub_bins);
}

}  // namespace detail

int sliceSize(const std::vector<Domain::VisitorData>& visitors,
              const TailCounts& tails) {
  int size = 0;
  for (const auto& visitor : visitors) size += recordSize(visitor, tails);
  return size;
}

}  // namespace june::visitor_wire

#endif  // USE_MPI
