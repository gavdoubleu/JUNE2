#ifdef USE_MPI

#include "parallel/visitor_payload.h"

#include "epidemiology/policy.h"

namespace june {

namespace {

double personModifier(const PolicyManager* policy_manager, const Person& person,
                      size_t mode, TransmissionEffectChannel channel) {
  return policy_manager ? policy_manager->personModifier(person, mode, channel)
                        : 1.0;
}

// Start of the stage `trajectory` is in at `time` by a forward walk; the
// infection time before the first transition.
double stageStartTime(const InfectionTrajectory& trajectory, double time) {
  double stage_start_time = trajectory.infection_time;
  for (const auto& transition : trajectory.transitions) {
    if (time < transition.first) break;
    stage_start_time = transition.first;
  }
  return stage_start_time;
}

}  // namespace

Domain::VisitorData buildVisitorPayload(const PersonLocation& location,
                                        const Person& person, int home_rank,
                                        double slot_start,
                                        const Disease& disease,
                                        const PolicyManager* policy_manager,
                                        const EmissionCalculator& calculator) {
  Domain::VisitorData visitor;
  visitor.person_id = location.person_id;
  visitor.home_rank = home_rank;
  visitor.venue_id = location.venue_id;
  visitor.subset_idx = location.subset_index;
  visitor.is_infected = (person.infection != nullptr);

  const double susceptibility =
      person.getSusceptibility(slot_start, disease.getName());
  visitor.immunity_level = static_cast<float>(1.0 - susceptibility);

  visitor.encounter_type_id = location.encounter_type_id;
  visitor.age = person.age;
  visitor.sex_code = static_cast<uint8_t>(person.sex);
  visitor.newly_infected = false;
  visitor.new_infection_time = -1.0;

  visitor.symptom_id =
      visitor.is_infected
          ? person.infection->getTrajectory().getCurrentSymptomId(slot_start)
          : 0;
  visitor.time_in_stage =
      visitor.is_infected
          ? slot_start -
                stageStartTime(person.infection->getTrajectory(), slot_start)
          : 0.0;

  const auto& modes = disease.getTransmissionParams().modes;
  const int num_modes = disease.numModes();
  visitor.target_susceptibility.resize(num_modes);
  for (int m = 0; m < num_modes; ++m) {
    visitor.target_susceptibility[m] =
        susceptibility *
        personModifier(policy_manager, person, static_cast<size_t>(m),
                       TransmissionEffectChannel::TargetSusceptibility);
  }

  calculator.emit(person, slot_start, visitor.emission);
  for (size_t m = 0; m < visitor.emission.infectiousness_by_mode.size(); ++m) {
    visitor.emission.infectiousness_by_mode[m] *=
        personModifier(policy_manager, person, m,
                       TransmissionEffectChannel::SourceInfectiousness);
  }
  if (visitor.is_infected) {
    for (size_t m = 0; m < modes.size(); ++m) {
      if (modes[m].type != TransmissionModeType::Fomite &&
          modes[m].type != TransmissionModeType::CompartmentalDeposition)
        continue;
      visitor.deposition_source_multiplier.push_back(
          personModifier(policy_manager, person, m,
                         TransmissionEffectChannel::SourceInfectiousness));
    }
  }
  // Same product as InteractionManager::addFomiteDeposits for locals.
  const FomiteSubBinSchedule& fomite_schedule = calculator.fomiteSchedule();
  const auto& sub_bins_per_mode = fomite_schedule.subBinsPerMode();
  size_t offset = 0;
  for (int local_fm = 0; local_fm < fomite_schedule.numModes() &&
                         offset < visitor.emission.fomite_deposits.size();
       ++local_fm) {
    const double source_modifier = personModifier(
        policy_manager, person, fomite_schedule.modes()[local_fm].mode_index,
        TransmissionEffectChannel::SourceInfectiousness);
    for (int k = 0; k < sub_bins_per_mode[local_fm]; ++k) {
      visitor.emission.fomite_deposits[offset + k] *= source_modifier;
    }
    offset += sub_bins_per_mode[local_fm];
  }
  return visitor;
}

}  // namespace june

#endif  // USE_MPI
