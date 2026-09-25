#pragma once

#include <vector>

#include "core/types.h"
#include "epidemiology/disease.h"
#include "epidemiology/fomite/fomite_sub_bins.h"

namespace june {

// What one Person emits over one slot. Each list's emptiness is its gate, so
// a gate and its data cannot disagree.
struct Emission {
  // Integrated infectiousness per Transmission Mode (hour units, 24 * ∫I dt),
  // Disease::numModes() long. Empty unless infectious at slot start.
  std::vector<double> infectiousness_by_mode;
  // Deposit per (fomite mode, sub-bin), flat in FomiteSubBinSchedule order,
  // totalSubBins() long. Empty unless infected (or no fomite modes).
  std::vector<double> fomite_deposits;
};

// Answers "what does this Person emit this slot?" from the Person, the
// Disease, the slot start and the slot length alone: it knows nothing of
// Venues, contact matrices or beta, which the Emission Site applies.
//
// The slot length is fixed at construction and the fomite sub-bin schedule is
// derived from it, so integrals and deposits always cover the same slot.
class EmissionCalculator {
 public:
  EmissionCalculator(const Disease& disease, double slot_hours);

  // For sites that size their fomite bins.
  const FomiteSubBinSchedule& fomiteSchedule() const {
    return fomite_schedule_;
  }

  // Overwrite `out` (caller-owned and reused, so the binning loop does not
  // allocate per Person) with `person`'s Emission for the slot starting at
  // `slot_start` (days). Callers decide who is present; dead Persons are
  // theirs to skip.
  void emit(const Person& person, double slot_start, Emission& out) const;

 private:
  const Disease& disease_;
  double slot_hours_;
  FomiteSubBinSchedule fomite_schedule_;
};

}  // namespace june
