#pragma once

#include <memory>
#include <vector>

namespace june {

class InfectiousnessCurve;
struct InfectionTrajectory;

// Integral over [interval_start, interval_end] (days) of the curve of the
// symptom held at each instant, in hour units (24 × the day integral).
//
// The interval splits at each transition inside it; each piece is integrated
// on its symptom's curve, with time measured from that stage's start. Before
// the first transition the symptom is 0, starting at the infection time. A
// transition at or before interval_start sets the stage without splitting; one
// at or after interval_end is ignored. A symptom with no curve (index out of
// range, or null) contributes 0.
double integrateStageCurves(
    const std::vector<std::shared_ptr<InfectiousnessCurve>>& curves_by_symptom,
    const InfectionTrajectory& trajectory, double interval_start,
    double interval_end);

}  // namespace june
