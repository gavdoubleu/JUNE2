#include "epidemiology/trajectory/stage_curve_integral.h"

#include "epidemiology/disease.h"

namespace june {

namespace {

// Integral of one stage's curve over [piece_start, piece_end], in hours.
double integrateStage(
    const std::vector<std::shared_ptr<InfectiousnessCurve>>& curves_by_symptom,
    uint16_t symptom_id, double stage_start_time, double piece_start,
    double piece_end) {
  if (symptom_id >= curves_by_symptom.size()) return 0.0;
  const auto& curve = curves_by_symptom[symptom_id];
  if (!curve) return 0.0;
  return curve->integrate(piece_start - stage_start_time,
                          piece_end - stage_start_time) *
         24.0;
}

}  // namespace

double integrateStageCurves(
    const std::vector<std::shared_ptr<InfectiousnessCurve>>& curves_by_symptom,
    const InfectionTrajectory& trajectory, double interval_start,
    double interval_end) {
  uint16_t symptom_id = 0;
  double stage_start_time = trajectory.infection_time;
  double piece_start = interval_start;
  double integral = 0.0;
  for (const auto& [transition_time, next_symptom_id] :
       trajectory.transitions) {
    if (transition_time >= interval_end) break;
    if (transition_time > interval_start) {
      integral +=
          integrateStage(curves_by_symptom, symptom_id, stage_start_time,
                         piece_start, transition_time);
      piece_start = transition_time;
    }
    symptom_id = next_symptom_id;
    stage_start_time = transition_time;
  }
  return integral + integrateStage(curves_by_symptom, symptom_id,
                                   stage_start_time, piece_start, interval_end);
}

}  // namespace june
