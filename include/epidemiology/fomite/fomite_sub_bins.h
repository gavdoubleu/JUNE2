#pragma once

#include <vector>

#include "epidemiology/disease.h"

namespace june {

// Handle to one fomite transmission mode. mode_index points into
// TransmissionParams::modes; config points into the same mode's config variant.
struct FomiteModeRef {
  int mode_index;
  const FomiteConfig* config;
};

// How one slot of `delta_hours` splits into fomite sub-bins, and the deposit
// an infected Person makes into each.
//
// The schedule lists the fomite modes in mode-index order, each with its
// sub-bin count n_sub = floor(delta_hours / sub_bin_time), clamped to at least
// 1 (sub_bin_time unset means one sub-bin per slot). It depends only on the
// Disease and the timestep, so every rank derives the same one.
//
// Deposits are flat over (fomite mode, sub-bin) in schedule order, length
// totalSubBins(). Every emission site goes through integrateDeposits, so a
// Person deposits the same amount wherever it is computed.
class FomiteSubBinSchedule {
 public:
  FomiteSubBinSchedule(const TransmissionParams& transmission,
                       double delta_hours);

  const std::vector<FomiteModeRef>& modes() const { return modes_; }
  const std::vector<int>& subBinsPerMode() const { return sub_bins_per_mode_; }
  int numModes() const { return static_cast<int>(modes_.size()); }
  int totalSubBins() const { return total_sub_bins_; }

  // Fill `deposits_out` (resized to totalSubBins()) with the deposit per
  // sub-bin of the slot [current_time, current_time + delta_hours / 24).
  // All zero when `infection` is null. Deposition is keyed by symptom, so an
  // infected Person may deposit before becoming infectious; a sub-bin
  // straddling a transition deposits on each symptom's curve in turn.
  void integrateDeposits(const Infection* infection, double current_time,
                         std::vector<double>& deposits_out) const;

 private:
  double delta_hours_;
  std::vector<FomiteModeRef> modes_;
  std::vector<int> sub_bins_per_mode_;
  int total_sub_bins_ = 0;
};

}  // namespace june
