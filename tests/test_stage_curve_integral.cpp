#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <memory>
#include <vector>

#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/infectiousness_curves.h"
#include "epidemiology/trajectory/stage_curve_integral.h"

using namespace june;

using Curves = std::vector<std::shared_ptr<InfectiousnessCurve>>;

// Trajectory infected at `infection_time` with the given (time, symptom_id)
// transitions.
static InfectionTrajectory makeTrajectory(
    double infection_time,
    std::vector<std::pair<double, uint16_t>> transitions) {
  InfectionTrajectory trajectory;
  trajectory.infection_time = infection_time;
  trajectory.transitions = std::move(transitions);
  return trajectory;
}

TEST_CASE("No transition inside the interval integrates one curve") {
  // Symptom 1 from day 0; ramp 0 → 10 over 10 days, so value = t.
  Curves curves = {nullptr, std::make_shared<LinearRampCurve>(0.0, 10.0, 10.0)};
  auto trajectory = makeTrajectory(0.0, {{0.0, 1}});

  // ∫_2^4 t dt = 6 days → 144 hours.
  CHECK(integrateStageCurves(curves, trajectory, 2.0, 4.0) ==
        doctest::Approx(144.0));
}

TEST_CASE("One transition inside the interval splits into two pieces") {
  // Symptom 1 from day 0: value = t. Symptom 2 from day 3: value = 2τ, τ from
  // its own stage start.
  Curves curves = {nullptr, std::make_shared<LinearRampCurve>(0.0, 10.0, 10.0),
                   std::make_shared<LinearRampCurve>(0.0, 20.0, 10.0)};
  auto trajectory = makeTrajectory(0.0, {{0.0, 1}, {3.0, 2}});

  // ∫_2^3 t dt + ∫_0^2 2τ dτ = 2.5 + 4 = 6.5 days → 156 hours.
  CHECK(integrateStageCurves(curves, trajectory, 2.0, 5.0) ==
        doctest::Approx(156.0));
}

TEST_CASE("Two transitions inside the interval split into three pieces") {
  // Symptom 1 from day 0: value = t. Symptom 2 from day 3: value = 2τ.
  // Symptom 3 from day 4: value = 5.
  Curves curves = {nullptr, std::make_shared<LinearRampCurve>(0.0, 10.0, 10.0),
                   std::make_shared<LinearRampCurve>(0.0, 20.0, 10.0),
                   std::make_shared<ConstantCurve>(5.0)};
  auto trajectory = makeTrajectory(0.0, {{0.0, 1}, {3.0, 2}, {4.0, 3}});

  // ∫_2^3 t dt + ∫_0^1 2τ dτ + 5 × 2 = 2.5 + 1 + 10 = 13.5 days → 324 hours.
  CHECK(integrateStageCurves(curves, trajectory, 2.0, 6.0) ==
        doctest::Approx(324.0));
}

TEST_CASE("A symptom with no curve contributes zero") {
  // Symptom 1 from day 0: value = 3. Symptom 0 (null curve) from day 2;
  // symptom 5 (beyond the curves) from day 3.
  Curves curves = {nullptr, std::make_shared<ConstantCurve>(3.0)};

  SUBCASE("Null curve") {
    auto trajectory = makeTrajectory(0.0, {{0.0, 1}, {2.0, 0}});
    // 3 × 1 day before the transition, nothing after → 72 hours.
    CHECK(integrateStageCurves(curves, trajectory, 1.0, 4.0) ==
          doctest::Approx(72.0));
  }

  SUBCASE("Symptom beyond the curves") {
    auto trajectory = makeTrajectory(0.0, {{0.0, 1}, {3.0, 5}});
    // 3 × 2 days before the transition, nothing after → 144 hours.
    CHECK(integrateStageCurves(curves, trajectory, 1.0, 4.0) ==
          doctest::Approx(144.0));
  }
}

TEST_CASE("A transition on an interval boundary belongs to the new stage") {
  // Symptom 1 from day 0: value = t. Symptom 2 from day 2: value = 2τ.
  Curves curves = {nullptr, std::make_shared<LinearRampCurve>(0.0, 10.0, 10.0),
                   std::make_shared<LinearRampCurve>(0.0, 20.0, 10.0)};
  auto trajectory = makeTrajectory(0.0, {{0.0, 1}, {2.0, 2}});

  SUBCASE("At the interval start the whole interval is the new stage") {
    // ∫_0^2 2τ dτ = 4 days → 96 hours.
    CHECK(integrateStageCurves(curves, trajectory, 2.0, 4.0) ==
          doctest::Approx(96.0));
  }

  SUBCASE("At the interval end the whole interval is the old stage") {
    // ∫_1^2 t dt = 1.5 days → 36 hours.
    CHECK(integrateStageCurves(curves, trajectory, 1.0, 2.0) ==
          doctest::Approx(36.0));
  }
}

TEST_CASE("Without transitions the symptom is 0 from the infection time") {
  // Symptom 0: value = τ from the infection time at day 1.
  Curves curves = {std::make_shared<LinearRampCurve>(0.0, 10.0, 10.0)};
  auto trajectory = makeTrajectory(1.0, {});

  // ∫_1^3 τ dτ = 4 days → 96 hours.
  CHECK(integrateStageCurves(curves, trajectory, 2.0, 4.0) ==
        doctest::Approx(96.0));
}
