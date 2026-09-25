#include "epidemiology/emission/emission.h"

namespace june {

EmissionCalculator::EmissionCalculator(const Disease& disease,
                                       double slot_hours)
    : disease_(disease),
      slot_hours_(slot_hours),
      fomite_schedule_(disease.getTransmissionParams(), slot_hours) {}

void EmissionCalculator::emit(const Person& person, double slot_start,
                              Emission& out) const {
  out.infectiousness_by_mode.clear();
  out.fomite_deposits.clear();
  const Infection* infection = person.infection.get();
  if (!infection) return;
  fomite_schedule_.integrateDeposits(infection, slot_start,
                                     out.fomite_deposits);
  if (!infection->isInfectious(slot_start)) return;
  const int num_modes = disease_.numModes();
  const double slot_end = slot_start + slot_hours_ / 24.0;
  out.infectiousness_by_mode.resize(num_modes);
  for (int mode = 0; mode < num_modes; ++mode) {
    out.infectiousness_by_mode[mode] =
        infection->getIntegratedInfectiousness(mode, slot_start, slot_end);
  }
}

}  // namespace june
