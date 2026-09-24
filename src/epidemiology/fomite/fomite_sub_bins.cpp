#include "epidemiology/fomite/fomite_sub_bins.h"

#include <algorithm>

namespace june {

FomiteSubBinSchedule::FomiteSubBinSchedule(
    const TransmissionParams& transmission, double delta_hours)
    : delta_hours_(delta_hours) {
  for (int mode_index = 0; mode_index < (int)transmission.modes.size();
       ++mode_index) {
    const auto& mode = transmission.modes[mode_index];
    if (mode.type == TransmissionModeType::Fomite) {
      modes_.push_back(
          FomiteModeRef{mode_index, &std::get<FomiteConfig>(mode.config)});
    }
  }
  sub_bins_per_mode_.assign(modes_.size(), 1);
  for (int local_fm = 0; local_fm < numModes(); ++local_fm) {
    double sub_bin_time = modes_[local_fm].cfg->sub_bin_time;
    sub_bins_per_mode_[local_fm] =
        (sub_bin_time > 0.0) ? std::max(1, (int)(delta_hours / sub_bin_time))
                             : 1;
    total_sub_bins_ += sub_bins_per_mode_[local_fm];
  }
}

void FomiteSubBinSchedule::integrateDeposits(
    const Infection* infection, double current_time,
    std::vector<double>& deposits_out) const {
  deposits_out.assign(total_sub_bins_, 0.0);
  if (!infection) return;
  const double t1 = current_time + delta_hours_ / 24.0;
  int offset = 0;
  for (int local_fm = 0; local_fm < numModes(); ++local_fm) {
    const int n_sub = sub_bins_per_mode_[local_fm];
    const double dt_sub = (t1 - current_time) / n_sub;
    for (int k = 0; k < n_sub; ++k) {
      const double t_sub_s = current_time + k * dt_sub;
      const double t_sub_e = current_time + (k + 1) * dt_sub;
      deposits_out[offset + k] = infection->getIntegratedFomiteDeposition(
          modes_[local_fm].mode_index, t_sub_s, t_sub_e);
    }
    offset += n_sub;
  }
}

}  // namespace june
