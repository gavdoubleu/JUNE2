#pragma once

#include <vector>

namespace june {

// What one Person emits over one slot. Each list's emptiness is its gate, so
// a gate and its data cannot disagree. Kept apart from EmissionCalculator so
// records (VisitorData, VisitorInfo) can carry it without pulling in Disease.
struct Emission {
  // Integrated infectiousness per Transmission Mode (hour units, 24 * ∫I dt),
  // Disease::numModes() long. Empty unless infectious at slot start.
  std::vector<double> infectiousness_by_mode;
  // Deposit per (fomite mode, sub-bin), flat in FomiteSubBinSchedule order,
  // totalSubBins() long. Empty unless infected (or no fomite modes).
  std::vector<double> fomite_deposits;

  bool operator==(const Emission&) const = default;
};

}  // namespace june
