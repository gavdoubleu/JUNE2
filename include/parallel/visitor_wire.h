#pragma once

#ifdef USE_MPI

#include "parallel/domain.h"

// Wire format of one Visitor record: a fixed header (WireRecord over
// VisitorData's plain fields) followed by two count-known-elsewhere tails,
// integrated_infectiousness then fomite_deposition_sub.
namespace june::visitor_wire {

// Lengths of a visitor record's two tails. Derived from the Disease and
// timestep, so identical on every rank and fixed for one exchange.
struct TailCounts {
  int num_modes;        // integrated_infectiousness
  int fomite_sub_bins;  // fomite_deposition_sub
};

// Bytes `visitor` occupies on the wire.
int recordSize(const Domain::VisitorData& visitor, const TailCounts& tails);

// Writes `visitor` at `ptr`, returns the end of the record. Throws if a tail's
// length differs from its count in `tails`.
char* pack(char* ptr, const Domain::VisitorData& visitor,
           const TailCounts& tails);

// Reads one record at `ptr` into `visitor`, returns the end of the record.
const char* unpack(const char* ptr, Domain::VisitorData& visitor,
                   const TailCounts& tails);

}  // namespace june::visitor_wire

#endif  // USE_MPI
