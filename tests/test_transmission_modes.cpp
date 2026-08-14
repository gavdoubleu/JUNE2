#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <limits>

#include "core/config.h"
#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/infectiousness_curves.h"
#include "epidemiology/interaction_manager.h"
#include "test_utils.h"
#include "utils/random.h"

using namespace june;

// =============================================================================
// Helper: Minimal STAGE_DRIVEN disease (healthy + mild, one trajectory)
// =============================================================================
static Disease makeStageDrivenDisease(
    std::shared_ptr<InfectiousnessCurve> curve, double stage_duration = 10.0) {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  tp.stage_curves["mild"] = curve;
  tp.symptom_id_curves = {nullptr, curve};  // id 0=healthy(null), id 1=mild

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general_population";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", stage_duration}}}});
  td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});

  return Disease("StageFlu", stags, {}, {td}, {}, tp);
}

// =============================================================================
// Helper: Three-symptom STAGE_DRIVEN disease (healthy, exposed, mild)
// =============================================================================
static Disease makeThreeStageStageDrivenDisease(
    std::shared_ptr<InfectiousnessCurve> exposed_curve,
    std::shared_ptr<InfectiousnessCurve> mild_curve, double exposed_dur = 2.0,
    double mild_dur = 8.0) {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;
  tp.stage_curves["exposed"] = exposed_curve;
  tp.stage_curves["mild"] = mild_curve;
  // id 0=healthy(null), id 1=exposed, id 2=mild
  tp.symptom_id_curves = {nullptr, exposed_curve, mild_curve};

  std::vector<SymptomTag> stags = {
      {"healthy", -1, 0}, {"exposed", 0, 1}, {"mild", 1, 2}};

  TrajectoryDefinition td;
  td.selection_key = "general_population";
  td.severity = 1.0;
  td.stages.push_back({"exposed", {"constant", {{"value", exposed_dur}}}});
  td.stages.push_back({"mild", {"constant", {{"value", mild_dur}}}});
  td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});

  return Disease("ThreeStageFlu", stags, {}, {td}, {}, tp);
}

// =============================================================================
// Helper: TRAJECTORY_DRIVEN disease with constant (deterministic) params
// =============================================================================
static Disease makeTrajectoryDrivenDisease(double max_inf = 5.0,
                                           double shape = 2.0,
                                           double rate = 1.0,
                                           double shift = 1.0,
                                           double inf_factor = 1.0) {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::TRAJECTORY_DRIVEN;
  tp.type = "gamma";
  tp.max_infectiousness = {"constant", {{"value", max_inf}}};
  tp.shape = {"constant", {{"value", shape}}};
  tp.rate = {"constant", {{"value", rate}}};
  tp.shift = {"constant", {{"value", shift}}};

  std::vector<SymptomTag> stags = {
      {"healthy", -1, 0}, {"asymptomatic", 0, 1}, {"mild", 1, 2}};

  TrajectoryDefinition td;
  td.selection_key = "general_population";
  td.severity = 1.0;
  td.infectiousness_factor = inf_factor;
  td.stages.push_back({"asymptomatic", {"constant", {{"value", 2.0}}}});
  td.stages.push_back({"mild", {"constant", {{"value", 8.0}}}});
  td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});

  return Disease("TrajFlu", stags, {}, {td}, {}, tp);
}

// =============================================================================
// Helper: Multi-mode STAGE_DRIVEN disease (respiratory + animal_bite)
// =============================================================================
static Disease makeMultiModeStageDrivenDisease() {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;

  auto curve_resp = std::make_shared<ConstantCurve>(2.0);
  auto curve_bite = std::make_shared<ConstantCurve>(0.8);

  // Legacy single-mode fields (populated for fallback)
  tp.stage_curves["mild"] = curve_resp;
  tp.symptom_id_curves = {nullptr, curve_resp};

  // mode 0 = respiratory, mode 1 = animal_bite
  {
    TransmissionMode m;
    m.name = "respiratory";
    m.mode_transmissibility_multiplier = 1.0;
    m.symptom_curves = {nullptr, curve_resp};
    tp.modes.push_back(std::move(m));
  }
  {
    TransmissionMode m;
    m.name = "animal_bite";
    m.mode_transmissibility_multiplier = 0.5;
    m.symptom_curves = {nullptr, curve_bite};
    tp.modes.push_back(std::move(m));
  }

  std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
  TrajectoryDefinition td;
  td.selection_key = "general_population";
  td.severity = 1.0;
  td.stages.push_back({"mild", {"constant", {{"value", 10.0}}}});
  td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});

  return Disease("MultiModeFlu", stags, {}, {td}, {}, tp);
}

// =============================================================================
// A. InfectiousnessCurve Unit Tests (6 tests)
// =============================================================================

TEST_CASE("A1: ConstantCurve returns fixed value at all times") {
  ConstantCurve c(3.5);
  CHECK(c.evaluate(0.0) == doctest::Approx(3.5));
  CHECK(c.evaluate(1.0) == doctest::Approx(3.5));
  CHECK(c.evaluate(100.0) == doctest::Approx(3.5));
  CHECK(c.evaluate(-5.0) == doctest::Approx(3.5));
}

TEST_CASE("A2: GammaCurve returns zero before shift and positive after") {
  GammaCurve g(5.0, 2.0, 1.0, 3.0);  // shift=3
  CHECK(g.evaluate(0.0) == doctest::Approx(0.0));
  CHECK(g.evaluate(1.0) == doctest::Approx(0.0));
  CHECK(g.evaluate(3.0) == doctest::Approx(0.0));  // shifted_t == 0
  CHECK(g.evaluate(4.0) > 0.0);
  CHECK(g.evaluate(5.0) > 0.0);
}

TEST_CASE("A3: GammaCurve peaks and decays (gamma shape)") {
  GammaCurve g(10.0, 2.0, 1.0, 0.0);  // no shift, shape=2, rate=1
  // Gamma PDF shape=2, rate=1: peaks at (shape-1)/rate = 1.0
  double val_before = g.evaluate(0.5);
  double val_peak = g.evaluate(1.0);
  double val_after = g.evaluate(3.0);
  double val_late = g.evaluate(10.0);

  CHECK(val_peak > val_before);  // Rise
  CHECK(val_peak > val_after);   // Decay after peak
  CHECK(val_after > val_late);   // Continued decay
  CHECK(val_peak > 0.0);
}

TEST_CASE("A4: LinearRampCurve interpolates and clamps") {
  LinearRampCurve lr(0.0, 1.0, 10.0);  // 0→1 over 10 days
  CHECK(lr.evaluate(0.0) == doctest::Approx(0.0));
  CHECK(lr.evaluate(5.0) == doctest::Approx(0.5));
  CHECK(lr.evaluate(10.0) == doctest::Approx(1.0));
  CHECK(lr.evaluate(15.0) == doctest::Approx(1.0));  // Clamped at end
  CHECK(lr.evaluate(-1.0) == doctest::Approx(0.0));  // Clamped at start
}

TEST_CASE("A5: ExponentialDecayCurve decays correctly") {
  ExponentialDecayCurve ed(2.0, 0.5, 0.0);  // initial=2, decay=0.5, no delay
  CHECK(ed.evaluate(0.0) == doctest::Approx(2.0));
  CHECK(ed.evaluate(1.0) == doctest::Approx(2.0 * std::exp(-0.5)));
  CHECK(ed.evaluate(2.0) == doctest::Approx(2.0 * std::exp(-1.0)));
  CHECK(ed.evaluate(10.0) == doctest::Approx(2.0 * std::exp(-5.0)));
}

TEST_CASE("A6: ExponentialDecayCurve respects delay") {
  ExponentialDecayCurve ed(3.0, 1.0, 2.0);          // delay=2
  CHECK(ed.evaluate(0.0) == doctest::Approx(3.0));  // Before delay
  CHECK(ed.evaluate(1.0) == doctest::Approx(3.0));  // Still before delay
  CHECK(ed.evaluate(2.0) == doctest::Approx(3.0));  // At delay boundary
  CHECK(ed.evaluate(3.0) == doctest::Approx(3.0 * std::exp(-1.0)));
  CHECK(ed.evaluate(5.0) == doctest::Approx(3.0 * std::exp(-3.0)));
}

// =============================================================================
// B. Trajectory-Driven Infectiousness (5 tests)
// =============================================================================

TEST_CASE("B1: Zero infectiousness during pre-infectious period") {
  Disease disease = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 1.0);  // shift=1
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // At t=0, time_since_infection=0, shifted_time = 0 - 1 = -1 → 0
  CHECK(p.infection->getInfectiousness(0.0) == doctest::Approx(0.0));
  // At t=0.5, shifted_time = 0.5 - 1 = -0.5 → 0
  CHECK(p.infection->getInfectiousness(0.5) == doctest::Approx(0.0));
  // At t=1.0, shifted_time = 0 → 0 (boundary)
  CHECK(p.infection->getInfectiousness(1.0) == doctest::Approx(0.0));
}

TEST_CASE("B2: Positive infectiousness after shift") {
  Disease disease = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 1.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // At t=2.0, time_since_infection=2, shifted_time=1 → gamma PDF > 0
  CHECK(p.infection->getInfectiousness(2.0) > 0.0);
  CHECK(p.infection->getInfectiousness(3.0) > 0.0);
}

TEST_CASE("B3: infectiousness_factor scales entire profile uniformly") {
  // Create two diseases identical except for infectiousness_factor
  Disease disease_full = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 0.0, 1.0);
  Disease disease_half = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 0.0, 0.5);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  // Infection with factor=1.0
  p.infection = std::make_unique<Infection>(&disease_full, 0.0, &p, 42, nullptr,
                                            "office", 0);
  double inf_full_t1 = p.infection->getInfectiousness(1.0);
  double inf_full_t3 = p.infection->getInfectiousness(3.0);

  // Infection with factor=0.5 (same seed so same gamma params)
  p.infection = std::make_unique<Infection>(&disease_half, 0.0, &p, 42, nullptr,
                                            "office", 0);
  double inf_half_t1 = p.infection->getInfectiousness(1.0);
  double inf_half_t3 = p.infection->getInfectiousness(3.0);

  CHECK(inf_full_t1 > 0.0);
  CHECK(inf_full_t3 > 0.0);

  // The factor is applied once at creation (baked into max_infectiousness_),
  // so the halved infection should be exactly 0.5x at every time point.
  CHECK(inf_half_t1 == doctest::Approx(inf_full_t1 * 0.5).epsilon(0.001));
  CHECK(inf_half_t3 == doctest::Approx(inf_full_t3 * 0.5).epsilon(0.001));

  // Verify absolute value: at t=1 (shifted_time=1), gamma(2,1) at x=1
  // = 1*exp(-1) ≈ 0.368. Full: 5.0 * 0.368 ≈ 1.839
  CHECK(inf_full_t1 ==
        doctest::Approx(5.0 * 1.0 * std::exp(-1.0)).epsilon(0.01));
}

TEST_CASE("B4: Gamma temporal profile shape") {
  Disease disease = makeTrajectoryDrivenDisease(1.0, 2.0, 1.0, 0.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // Gamma(shape=2, rate=1) peaks at x=1
  double v_early = p.infection->getInfectiousness(0.5);  // shifted=0.5
  double v_peak = p.infection->getInfectiousness(1.0);   // shifted=1.0 (peak)
  double v_later = p.infection->getInfectiousness(1.5);  // shifted=1.5

  CHECK(v_peak > v_early);
  CHECK(v_peak > v_later);
}

TEST_CASE("B5: Fixed seed produces identical infectiousness") {
  Disease disease = makeTrajectoryDrivenDisease();
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  // Create infection with same seed twice
  p.infection = std::make_unique<Infection>(&disease, 0.0, &p, 12345, nullptr,
                                            "office", 0);
  double inf1 = p.infection->getInfectiousness(3.0);

  p.infection = std::make_unique<Infection>(&disease, 0.0, &p, 12345, nullptr,
                                            "office", 0);
  double inf2 = p.infection->getInfectiousness(3.0);

  CHECK(inf1 == doctest::Approx(inf2));
}

// =============================================================================
// B6-B10: infectiousness_factor and shift+incubation tests
// =============================================================================

// Helper: TRAJECTORY_DRIVEN disease with exposed stage
static Disease makeTrajectoryDrivenDiseaseWithExposed(
    double max_inf = 5.0, double shift = -2.0, double inf_factor = 1.0,
    double exposed_dur = 5.0) {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::TRAJECTORY_DRIVEN;
  tp.type = "gamma";
  tp.max_infectiousness = {"constant", {{"value", max_inf}}};
  tp.shape = {"constant", {{"value", 2.0}}};
  tp.rate = {"constant", {{"value", 1.0}}};
  tp.shift = {"constant", {{"value", shift}}};

  std::vector<SymptomTag> stags = {{"healthy", -1, 0},
                                   {"exposed", 0, 1},
                                   {"asymptomatic", 1, 2},
                                   {"mild", 2, 3}};

  TrajectoryDefinition td;
  td.selection_key = "general_population";
  td.severity = 1.0;
  td.infectiousness_factor = inf_factor;
  td.stages.push_back({"exposed", {"constant", {{"value", exposed_dur}}}});
  td.stages.push_back({"asymptomatic", {"constant", {{"value", 7.0}}}});
  td.stages.push_back({"mild", {"constant", {{"value", 8.0}}}});
  td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});

  return Disease("TrajFluExposed", stags, {}, {td}, {}, tp);
}

TEST_CASE("B6: infectiousness_factor defaults to 1.0 when omitted") {
  // makeTrajectoryDrivenDisease defaults inf_factor to 1.0
  Disease disease = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 0.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // At t=1 (shifted_time=1), gamma(2,1) at x=1 = exp(-1) ≈ 0.368
  // Full profile: 5.0 * 0.368 = 1.839 (no scaling)
  double inf = p.infection->getInfectiousness(1.0);
  CHECK(inf == doctest::Approx(5.0 * 1.0 * std::exp(-1.0)).epsilon(0.01));
}

TEST_CASE("B7: infectiousness_factor is applied once, not per-timestep") {
  // With factor=0.5, the profile should be uniformly halved at ALL time points,
  // regardless of which symptom stage the person is in.
  Disease disease_half = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 0.0, 0.5);
  Disease disease_full = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 0.0, 1.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  // Trajectory: asymptomatic 0-2, mild 2-10
  // Test at t=1 (asymptomatic) and t=3 (mild) — ratio should be the same
  p.infection = std::make_unique<Infection>(&disease_full, 0.0, &p, 42, nullptr,
                                            "office", 0);
  double full_asymp = p.infection->getInfectiousness(1.0);
  double full_mild = p.infection->getInfectiousness(3.0);

  p.infection = std::make_unique<Infection>(&disease_half, 0.0, &p, 42, nullptr,
                                            "office", 0);
  double half_asymp = p.infection->getInfectiousness(1.0);
  double half_mild = p.infection->getInfectiousness(3.0);

  // Both stages should show exactly 0.5x — the factor is baked in uniformly
  CHECK(half_asymp == doctest::Approx(full_asymp * 0.5).epsilon(0.001));
  CHECK(half_mild == doctest::Approx(full_mild * 0.5).epsilon(0.001));
}

TEST_CASE("B8: Shift is linked to exposed stage duration") {
  // With shift=-2.0 and exposed_dur=5.0, effective shift = -2.0 + 5.0 = 3.0
  // So shifted_time = time_since_infection - 3.0
  // Infectiousness should be zero before t=3.0
  Disease disease = makeTrajectoryDrivenDiseaseWithExposed(5.0, -2.0, 1.0, 5.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // Before effective shift: no infectiousness
  CHECK(p.infection->getInfectiousness(0.0) == doctest::Approx(0.0));
  CHECK(p.infection->getInfectiousness(2.0) == doctest::Approx(0.0));
  CHECK(p.infection->getInfectiousness(3.0) == doctest::Approx(0.0));

  // After effective shift: infectiousness rises
  // At t=4.0, shifted_time=1.0, gamma(2,1) at x=1 = exp(-1)
  CHECK(p.infection->getInfectiousness(4.0) ==
        doctest::Approx(5.0 * 1.0 * std::exp(-1.0)).epsilon(0.01));
}

TEST_CASE("B9: infectiousness_factor works with exposed+shift linkage") {
  // Both features combined: shift linked to incubation AND factor applied
  Disease disease = makeTrajectoryDrivenDiseaseWithExposed(5.0, -2.0, 0.5, 5.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // Effective shift = -2.0 + 5.0 = 3.0
  // At t=4.0, shifted_time=1.0, gamma peak region
  // Expected: 5.0 * 0.5 * gamma(1.0) = 2.5 * exp(-1) ≈ 0.920
  CHECK(p.infection->getInfectiousness(4.0) ==
        doctest::Approx(5.0 * 0.5 * 1.0 * std::exp(-1.0)).epsilon(0.01));
}

TEST_CASE("B10: Stage-Driven mode ignores infectiousness_factor") {
  // Stage-Driven uses per-stage curves, not the global gamma profile.
  // infectiousness_factor should have no effect.
  auto curve = std::make_shared<ConstantCurve>(3.0);
  Disease disease_full = makeStageDrivenDisease(curve, 10.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection = std::make_unique<Infection>(&disease_full, 0.0, &p, 42, nullptr,
                                            "office", 0);
  double inf = p.infection->getInfectiousness(1.0);

  // Stage-Driven constant curve = 3.0, unaffected by any trajectory factor
  CHECK(inf == doctest::Approx(3.0).epsilon(0.01));
}

// =============================================================================
// C. Stage-Driven Infectiousness (6 tests)
// =============================================================================

TEST_CASE("C1: Infectiousness uses time_in_stage") {
  auto ramp = std::make_shared<LinearRampCurve>(0.0, 1.0, 10.0);
  Disease disease = makeStageDrivenDisease(ramp, 10.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // At t=5.0, person is in "mild" stage, time_in_stage=5.0
  // LinearRamp(0,1,10) at t=5 → 0.5
  CHECK(p.infection->getInfectiousness(5.0) == doctest::Approx(0.5));
  // At t=0.0, time_in_stage=0 → 0.0
  CHECK(p.infection->getInfectiousness(0.0) == doctest::Approx(0.0));
}

TEST_CASE("C2: Infectiousness resets at stage transition") {
  auto exposed_curve = std::make_shared<LinearRampCurve>(0.0, 1.0, 2.0);
  auto mild_curve = std::make_shared<ConstantCurve>(5.0);
  Disease disease =
      makeThreeStageStageDrivenDisease(exposed_curve, mild_curve, 2.0, 8.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // At t=1.0: exposed stage, time_in_stage=1.0, ramp at 1.0 → 0.5
  double inf_exposed = p.infection->getInfectiousness(1.0);
  CHECK(inf_exposed == doctest::Approx(0.5));

  // At t=2.5: mild stage (transition at t=2.0), time_in_stage=0.5
  // ConstantCurve(5.0) → 5.0
  double inf_mild = p.infection->getInfectiousness(2.5);
  CHECK(inf_mild == doctest::Approx(5.0));
}

TEST_CASE("C3: Different curve types per stage") {
  auto exposed_curve = std::make_shared<ConstantCurve>(1.0);
  auto mild_curve = std::make_shared<LinearRampCurve>(0.0, 2.0, 4.0);
  Disease disease =
      makeThreeStageStageDrivenDisease(exposed_curve, mild_curve, 2.0, 8.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // Exposed at t=1: constant → 1.0
  CHECK(p.infection->getInfectiousness(1.0) == doctest::Approx(1.0));

  // Mild at t=4.0: transition at 2.0, time_in_stage=2.0
  // LinearRamp(0,2,4) at t=2.0 → 1.0
  CHECK(p.infection->getInfectiousness(4.0) == doctest::Approx(1.0));

  // Mild at t=6.0: time_in_stage=4.0 → ramp clamped at 2.0
  CHECK(p.infection->getInfectiousness(6.0) == doctest::Approx(2.0));
}

TEST_CASE("C4: Null curve returns zero infectiousness") {
  // Healthy symptom has no curve (nullptr) → should return 0
  auto curve = std::make_shared<ConstantCurve>(1.0);
  Disease disease = makeStageDrivenDisease(curve, 5.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // After recovery (t=6.0, healthy stage), curve is nullptr → 0
  CHECK(p.infection->getInfectiousness(6.0) == doctest::Approx(0.0));
}

TEST_CASE("C5: evaluateStageDrivenInfectiousness direct test") {
  auto ramp = std::make_shared<LinearRampCurve>(0.0, 1.0, 10.0);
  Disease disease = makeStageDrivenDisease(ramp, 10.0);

  // symptom_id=1 (mild), time_in_stage=5.0 → ramp at 5 = 0.5
  CHECK(disease.evaluateStageDrivenInfectiousness(0, 1, 5.0f) ==
        doctest::Approx(0.5));

  // symptom_id=0 (healthy), curve is nullptr → 0
  CHECK(disease.evaluateStageDrivenInfectiousness(0, 0, 5.0f) ==
        doctest::Approx(0.0));
}

TEST_CASE(
    "C6: evaluateStageDrivenInfectiousness returns 0 for trajectory-driven") {
  Disease disease = makeTrajectoryDrivenDisease();

  // Should return 0 and print a warning (once)
  CHECK(disease.evaluateStageDrivenInfectiousness(0, 1, 5.0f) ==
        doctest::Approx(0.0));
}

// =============================================================================
// D. Mode Comparison (2 tests)
// =============================================================================

TEST_CASE("D1: Trajectory-driven is continuous across stage transitions") {
  Disease disease = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 0.0);
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // Transition from asymptomatic to mild at t=2.0
  // With infectiousness_factor applied once at creation (not per-timestep),
  // there should be NO discontinuity at stage boundaries — only the smooth
  // gamma profile.
  double inf_before = p.infection->getInfectiousness(1.99);
  double inf_after = p.infection->getInfectiousness(2.01);

  CHECK(inf_before > 0.0);
  CHECK(inf_after > 0.0);

  // The ratio should be ~1.0 (smooth gamma, no multiplier jump)
  double ratio = inf_after / inf_before;
  CHECK(ratio == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("D2: Both modes return zero for recovered/dead stages") {
  // Stage-driven: healthy has nullptr curve → 0
  auto curve = std::make_shared<ConstantCurve>(1.0);
  Disease sd_disease = makeStageDrivenDisease(curve, 5.0);
  WorldState world1 = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p1 = world1.people[0];
  p1.infection = std::make_unique<Infection>(&sd_disease, 0.0, &p1, 42, nullptr,
                                             "office", 0);
  CHECK(p1.infection->getInfectiousness(6.0) == doctest::Approx(0.0));

  // Trajectory-driven: after all stages the gamma profile has decayed
  // to near-zero (no per-stage multiplier to force exact zero).
  Disease td_disease = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 0.0);
  WorldState world2 = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p2 = world2.people[0];
  p2.infection = std::make_unique<Infection>(&td_disease, 0.0, &p2, 42, nullptr,
                                             "office", 0);
  // After all stages (asymptomatic 2d + mild 8d = 10d), gamma(2,1) at x=11
  // = 11*exp(-11) ≈ 0.00016 — negligible infectiousness
  CHECK(p2.infection->getInfectiousness(11.0) < 0.01);
}

// =============================================================================
// E. Multi-Mode Stage-Driven (4 tests)
// =============================================================================

TEST_CASE("E1: Different mode indices yield different values") {
  Disease disease = makeMultiModeStageDrivenDisease();
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // Mode 0 (respiratory): ConstantCurve(2.0) → 2.0
  double inf_mode0 = p.infection->getInfectiousness(0, 1.0);
  // Mode 1 (animal_bite): ConstantCurve(0.8) → 0.8
  double inf_mode1 = p.infection->getInfectiousness(1, 1.0);

  CHECK(inf_mode0 == doctest::Approx(2.0));
  CHECK(inf_mode1 == doctest::Approx(0.8));
}

TEST_CASE("E2: Susceptibility multipliers stored correctly") {
  Disease disease = makeMultiModeStageDrivenDisease();
  const auto& tp = disease.getTransmissionParams();

  REQUIRE(tp.modes.size() == 2);
  CHECK(tp.modes[0].mode_transmissibility_multiplier == doctest::Approx(1.0));
  CHECK(tp.modes[1].mode_transmissibility_multiplier == doctest::Approx(0.5));
}

TEST_CASE("E3: Out-of-range mode_index falls back to mode 0") {
  Disease disease = makeMultiModeStageDrivenDisease();
  WorldState world = TestWorldFactory::createMinimalWorld(1, 0);
  Person& p = world.people[0];

  p.infection =
      std::make_unique<Infection>(&disease, 0.0, &p, 42, nullptr, "office", 0);

  // Mode 99 (out of range) should fall back to mode 0 (respiratory, 2.0)
  double inf_oor = p.infection->getInfectiousness(99, 1.0);
  double inf_mode0 = p.infection->getInfectiousness(0, 1.0);
  CHECK(inf_oor == doctest::Approx(inf_mode0));
}

TEST_CASE("E4: evaluateStageDrivenInfectiousness multi-mode dispatch") {
  Disease disease = makeMultiModeStageDrivenDisease();

  // Mode 0 (respiratory), symptom_id=1 (mild), time=1.0 → 2.0
  CHECK(disease.evaluateStageDrivenInfectiousness(0, 1, 1.0f) ==
        doctest::Approx(2.0));
  // Mode 1 (animal_bite), symptom_id=1 (mild), time=1.0 → 0.8
  CHECK(disease.evaluateStageDrivenInfectiousness(1, 1, 1.0f) ==
        doctest::Approx(0.8));
}

// =============================================================================
// F. Transmission Integration (5 tests)
// =============================================================================

TEST_CASE("F1: Trajectory-driven transmission occurs") {
  WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
  world.venues[0].type_id = 0;

  Disease disease = makeTrajectoryDrivenDisease(5.0, 2.0, 1.0, 0.0);

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix =
      default_contact_matrix;  // High contacts for near-certain transmission
  cm.default_beta = 1.0;       // Neutralize beta dampening for this test
  cm.default_characteristic_time = 1.0;  // Match delta_hours for unit scaling
  SimulationConfig sim;
  sim.random_seed = 123;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease, nullptr);

  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 123, nullptr, "office", 0);

  std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                      {1, 0, -1, 0, 255, 1}};

  im.processTransmissions(locs, 1.0, 1.0, nullptr);
  CHECK(world.people[1].infection != nullptr);
}

TEST_CASE("F2: Stage-driven transmission occurs") {
  WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
  world.venues[0].type_id = 0;

  auto curve = std::make_shared<ConstantCurve>(1.0);
  Disease disease = makeStageDrivenDisease(curve, 10.0);

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix =
      default_contact_matrix;  // High contacts for near-certain transmission
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease, nullptr);

  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 123, nullptr, "office", 0);

  std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                      {1, 0, -1, 0, 255, 1}};

  im.processTransmissions(locs, 1.0, 1.0, nullptr);
  CHECK(world.people[1].infection != nullptr);
}

TEST_CASE("F3: Zero infectiousness prevents transmission") {
  WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
  world.venues[0].type_id = 0;

  auto zero_curve = std::make_shared<ConstantCurve>(0.0);
  Disease disease = makeStageDrivenDisease(zero_curve, 10.0);

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease, nullptr);

  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 123, nullptr, "office", 0);

  std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                      {1, 0, -1, 0, 255, 1}};

  im.processTransmissions(locs, 1.0, 1.0, nullptr);
  CHECK(world.people[1].infection == nullptr);
}

TEST_CASE("F4: Immune person resists infection") {
  WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
  world.venues[0].type_id = 0;

  auto curve = std::make_shared<ConstantCurve>(1.0);
  Disease disease = makeStageDrivenDisease(curve, 10.0);

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease, nullptr);

  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 123, nullptr, "office", 0);

  // Give person 1 full natural immunity
  world.people[1].immunity.natural_level = 1.0;
  world.people[1].immunity.natural_acquisition_time = 0.0;
  world.people[1].immunity.natural_waning_rate = 0.0;

  std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                      {1, 0, -1, 0, 255, 1}};

  im.processTransmissions(locs, 1.0, 1.0, nullptr);
  CHECK(world.people[1].infection == nullptr);
}

TEST_CASE("F5: Higher infectiousness yields higher infection probability") {
  // Statistical test: high infectiousness should produce more infections
  // than low infectiousness over many trials
  int high_infections = 0;
  int low_infections = 0;
  const int N = 200;

  for (int trial = 0; trial < N; ++trial) {
    GlobalRNG::seed(trial * 1000);

    // High infectiousness trial
    {
      WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
      world.venues[0].type_id = 0;

      auto high_curve = std::make_shared<ConstantCurve>(10.0);
      Disease disease = makeStageDrivenDisease(high_curve, 10.0);

      ContactMatrixConfig cm;
      ContactMatrix default_contact_matrix;
      default_contact_matrix.bins = {"all"};
      default_contact_matrix.contacts = {{1.0}};
      cm.default_matrix = default_contact_matrix;
      SimulationConfig sim;
      ParallelConfig par;
      cm.allow_default_matrix = true;
      finalizeContactMatrices(cm, world, disease);
      InteractionManager im(world, cm, sim, par, &disease, nullptr);

      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);

      std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                          {1, 0, -1, 0, 255, 1}};

      im.processTransmissions(locs, 1.0, 1.0, nullptr);
      if (world.people[1].infection != nullptr) high_infections++;
    }

    // Low infectiousness trial
    {
      WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
      world.venues[0].type_id = 0;

      auto low_curve = std::make_shared<ConstantCurve>(0.01);
      Disease disease = makeStageDrivenDisease(low_curve, 10.0);

      ContactMatrixConfig cm;
      ContactMatrix default_contact_matrix;
      default_contact_matrix.bins = {"all"};
      default_contact_matrix.contacts = {{1.0}};
      cm.default_matrix = default_contact_matrix;
      SimulationConfig sim;
      ParallelConfig par;
      cm.allow_default_matrix = true;
      finalizeContactMatrices(cm, world, disease);
      InteractionManager im(world, cm, sim, par, &disease, nullptr);

      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);

      std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                          {1, 0, -1, 0, 255, 1}};

      im.processTransmissions(locs, 1.0, 1.0, nullptr);
      if (world.people[1].infection != nullptr) low_infections++;
    }
  }

  CHECK(high_infections > low_infections);
}

// =============================================================================
// G. Edge Cases (5 tests)
// =============================================================================

TEST_CASE("G1: Negative time_in_stage handled gracefully") {
  ConstantCurve cc(1.0);
  CHECK(std::isfinite(cc.evaluate(-1.0)));

  GammaCurve gc(5.0, 2.0, 1.0, 0.0);
  CHECK(gc.evaluate(-1.0) == doctest::Approx(0.0));  // shifted_t <= 0

  LinearRampCurve lr(0.0, 1.0, 10.0);
  CHECK(lr.evaluate(-1.0) == doctest::Approx(0.0));  // Clamped at start

  ExponentialDecayCurve ed(2.0, 0.5, 0.0);
  CHECK(std::isfinite(ed.evaluate(-1.0)));
}

TEST_CASE("G2: Very large time values: no overflow/NaN") {
  double big_t = 1e10;

  ConstantCurve cc(1.0);
  CHECK(std::isfinite(cc.evaluate(big_t)));

  GammaCurve gc(5.0, 2.0, 1.0, 0.0);
  CHECK(std::isfinite(gc.evaluate(big_t)));

  LinearRampCurve lr(0.0, 1.0, 10.0);
  CHECK(std::isfinite(lr.evaluate(big_t)));

  ExponentialDecayCurve ed(2.0, 0.5, 0.0);
  CHECK(std::isfinite(ed.evaluate(big_t)));
}

TEST_CASE(
    "G3: All people infectious, no susceptibles yields zero new infections") {
  WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
  world.venues[0].type_id = 0;

  auto curve = std::make_shared<ConstantCurve>(1.0);
  Disease disease = makeStageDrivenDisease(curve, 10.0);

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease, nullptr);

  // Both people infected
  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 123, nullptr, "office", 0);
  world.people[1].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[1], 456, nullptr, "office", 0);

  std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                      {1, 0, -1, 0, 255, 1}};

  int new_infections = im.processTransmissions(locs, 1.0, 1.0, nullptr);
  CHECK(new_infections == 0);
}

TEST_CASE("G4: Empty venue returns zero infections") {
  WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
  world.venues[0].type_id = 0;

  auto curve = std::make_shared<ConstantCurve>(1.0);
  Disease disease = makeStageDrivenDisease(curve, 10.0);

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease, nullptr);

  // Empty locations
  std::vector<PersonLocation> locs;
  int new_infections = im.processTransmissions(locs, 1.0, 1.0, nullptr);
  CHECK(new_infections == 0);
}

TEST_CASE("G5: Single person at venue, no self-infection") {
  WorldState world = TestWorldFactory::createMinimalWorld(1, 1);
  world.venues[0].type_id = 0;

  auto curve = std::make_shared<ConstantCurve>(1.0);
  Disease disease = makeStageDrivenDisease(curve, 10.0);

  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{100.0}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease, nullptr);

  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 123, nullptr, "office", 0);

  std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0}};

  int new_infections = im.processTransmissions(locs, 1.0, 1.0, nullptr);
  CHECK(new_infections == 0);
}

// =============================================================================
// H. Plague-like Multi-Mode Transmission Regression Tests
//    These tests verify the FoI formula matches the original
//    Bubonic_and_pneumonic branch:
//      lambda = delta_hours * C / N * I * mode_transmissibility_multiplier
//    and that mode_transmissibility_multipliers are correctly applied per mode.
// =============================================================================

// Helper to run processTransmissions with a fresh InteractionManager
static int runTransmission(WorldState& world, Disease& disease,
                           std::vector<PersonLocation>& locs, double time,
                           double delta, double contacts = 1.0,
                           unsigned int seed = 0) {
  ContactMatrixConfig cm;
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{contacts}};
  cm.default_matrix = default_contact_matrix;
  SimulationConfig sim;
  sim.random_seed = seed;
  ParallelConfig par;
  cm.allow_default_matrix = true;
  finalizeContactMatrices(cm, world, disease);
  InteractionManager im(world, cm, sim, par, &disease, nullptr);
  return im.processTransmissions(locs, time, delta, nullptr);
}

// Helper: Creates a plague-like multi-mode stage-driven disease with
// animal_bite and respiratory modes, matching Bubonic_and_pneumonic config.
static Disease makePlagueLikeDisease() {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;

  // Symptom tags: 0=recovered, 1=healthy, 2=exposed, 3=buboes_and_fever,
  //               4=bacteraemia, 5=secondary_pneumonic
  auto zero = std::make_shared<ConstantCurve>(0.0);
  auto bite_buboes = std::make_shared<LinearRampCurve>(0.0, 0.5, 4.0);
  auto bite_bacteraemia = std::make_shared<ConstantCurve>(1.0);
  auto bite_secondary = std::make_shared<ConstantCurve>(1.0);
  auto resp_secondary = std::make_shared<ConstantCurve>(1.5);

  std::vector<std::shared_ptr<InfectiousnessCurve>> bite_curves = {
      zero, zero, zero, bite_buboes, bite_bacteraemia, bite_secondary};
  std::vector<std::shared_ptr<InfectiousnessCurve>> resp_curves = {
      zero, zero, zero, zero, zero, resp_secondary};

  tp.symptom_id_curves = bite_curves;

  {
    TransmissionMode m;
    m.name = "animal_bite";
    m.mode_transmissibility_multiplier = 0.08;
    m.symptom_curves = bite_curves;
    tp.modes.push_back(std::move(m));
  }
  {
    TransmissionMode m;
    m.name = "respiratory";
    m.mode_transmissibility_multiplier = 0.12;
    m.symptom_curves = resp_curves;
    tp.modes.push_back(std::move(m));
  }

  std::vector<SymptomTag> stags = {
      {"recovered", -2, 0},  {"healthy", -1, 1},
      {"exposed", 0, 2},     {"buboes_and_fever", 2, 3},
      {"bacteraemia", 4, 4}, {"secondary_pneumonic", 5, 5}};

  DiseaseStageSettings settings;
  settings.recovered_stages = {"recovered"};

  // Trajectory: exposed -> buboes -> bacteraemia -> secondary_pneumonic ->
  // recovered
  TrajectoryDefinition td;
  td.selection_key = "bubonic_path";
  td.severity = 1.0;
  td.stages.push_back({"exposed", {"constant", {{"value", 2.0}}}});
  td.stages.push_back({"buboes_and_fever", {"constant", {{"value", 4.0}}}});
  td.stages.push_back({"bacteraemia", {"constant", {{"value", 1.0}}}});
  td.stages.push_back({"secondary_pneumonic", {"constant", {{"value", 4.0}}}});
  td.stages.push_back({"recovered", {"constant", {{"value", 100.0}}}});

  OutcomeRates rates;
  OutcomeRow row;
  row.probabilities["bubonic_path"] = 1.0;
  rates.rows.push_back(row);

  return Disease("plague", stags, settings, {td}, rates, tp);
}

TEST_CASE("H1: FoI scales with delta_hours (no beta/char_time division)") {
  // The ported code introduced beta * (delta/char_time) which reduced FoI by
  // ~480x with default parameters. This test verifies FoI scales linearly
  // with delta_hours, catching any formula regressions.
  int infections_short = 0;
  int infections_long = 0;
  const int N = 200;

  for (int trial = 0; trial < N; ++trial) {
    GlobalRNG::seed(trial * 100);
    {
      WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
      world.venues[0].type_id = 0;
      auto curve = std::make_shared<ConstantCurve>(1.0);
      Disease disease = makeStageDrivenDisease(curve, 10.0);
      ContactMatrixConfig cm;
      ContactMatrix default_contact_matrix;
      default_contact_matrix.bins = {"all"};
      default_contact_matrix.contacts = {{1.0}};
      cm.default_matrix = default_contact_matrix;
      SimulationConfig sim;
      sim.random_seed = trial * 100;
      ParallelConfig par;
      cm.allow_default_matrix = true;
      finalizeContactMatrices(cm, world, disease);
      InteractionManager im(world, cm, sim, par, &disease, nullptr);
      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);
      std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                          {1, 0, -1, 0, 255, 1}};
      im.processTransmissions(locs, 1.0, 1.0, nullptr);
      if (world.people[1].infection) infections_short++;
    }

    GlobalRNG::seed(trial * 100);
    {
      WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
      world.venues[0].type_id = 0;
      auto curve = std::make_shared<ConstantCurve>(1.0);
      Disease disease = makeStageDrivenDisease(curve, 10.0);
      ContactMatrixConfig cm;
      ContactMatrix default_contact_matrix;
      default_contact_matrix.bins = {"all"};
      default_contact_matrix.contacts = {{1.0}};
      cm.default_matrix = default_contact_matrix;
      SimulationConfig sim;
      sim.random_seed = trial * 100;
      ParallelConfig par;
      cm.allow_default_matrix = true;
      finalizeContactMatrices(cm, world, disease);
      InteractionManager im(world, cm, sim, par, &disease, nullptr);
      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);
      std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                          {1, 0, -1, 0, 255, 1}};
      im.processTransmissions(locs, 1.0, 8.0, nullptr);
      if (world.people[1].infection) infections_long++;
    }
  }

  // 1h: lambda = 1*2/1*1 = 2, prob ≈ 86.5%
  // 8h: lambda = 8*2/1*1 = 16, prob ≈ 100%
  CHECK(infections_long >= infections_short);
  CHECK(infections_short > 100);  // 86.5% → expect ~173/200
}

TEST_CASE("H2: Susceptibility multiplier dampens per-mode transmission") {
  // Verifies mode_transmissibility_multipliers are applied. The original
  // Bubonic_and_pneumonic used 0.08 for animal_bite. During the port these
  // were changed to 1.0, removing mode-specific dampening.
  int infections_low = 0;
  int infections_high = 0;
  const int N = 500;

  for (int trial = 0; trial < N; ++trial) {
    // Low mode_transmissibility_multiplier (original plague: 0.08)
    GlobalRNG::seed(trial * 200);
    {
      WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
      world.venues[0].type_id = 0;
      TransmissionParams tp;
      tp.mode = InfectiousnessMode::STAGE_DRIVEN;
      auto curve = std::make_shared<ConstantCurve>(1.0);
      tp.symptom_id_curves = {nullptr, curve};
      {
        TransmissionMode m;
        m.name = "animal_bite";
        m.mode_transmissibility_multiplier = 0.08;
        m.symptom_curves = {nullptr, curve};
        tp.modes.push_back(std::move(m));
      }
      std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
      TrajectoryDefinition td;
      td.selection_key = "general";
      td.severity = 1.0;
      td.stages.push_back({"mild", {"constant", {{"value", 10.0}}}});
      td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});
      Disease disease("Low", stags, {}, {td}, {}, tp);
      ContactMatrixConfig cm;
      ContactMatrix default_contact_matrix;
      default_contact_matrix.bins = {"all"};
      default_contact_matrix.contacts = {{1.0}};
      cm.default_matrix = default_contact_matrix;
      SimulationConfig sim;
      sim.random_seed = trial * 200;
      ParallelConfig par;
      cm.allow_default_matrix = true;
      finalizeContactMatrices(cm, world, disease);
      InteractionManager im(world, cm, sim, par, &disease, nullptr);
      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);
      std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                          {1, 0, -1, 0, 255, 1}};
      im.processTransmissions(locs, 1.0, 8.0, nullptr);
      if (world.people[1].infection) infections_low++;
    }

    // High mode_transmissibility_multiplier (incorrectly ported: 1.0)
    GlobalRNG::seed(trial * 200);
    {
      WorldState world = TestWorldFactory::createMinimalWorld(2, 1);
      world.venues[0].type_id = 0;
      TransmissionParams tp;
      tp.mode = InfectiousnessMode::STAGE_DRIVEN;
      auto curve = std::make_shared<ConstantCurve>(1.0);
      tp.symptom_id_curves = {nullptr, curve};
      {
        TransmissionMode m;
        m.name = "animal_bite";
        m.mode_transmissibility_multiplier = 1.0;
        m.symptom_curves = {nullptr, curve};
        tp.modes.push_back(std::move(m));
      }
      std::vector<SymptomTag> stags = {{"healthy", -1, 0}, {"mild", 1, 1}};
      TrajectoryDefinition td;
      td.selection_key = "general";
      td.severity = 1.0;
      td.stages.push_back({"mild", {"constant", {{"value", 10.0}}}});
      td.stages.push_back({"healthy", {"constant", {{"value", 100.0}}}});
      Disease disease("High", stags, {}, {td}, {}, tp);
      ContactMatrixConfig cm;
      ContactMatrix default_contact_matrix;
      default_contact_matrix.bins = {"all"};
      default_contact_matrix.contacts = {{1.0}};
      cm.default_matrix = default_contact_matrix;
      SimulationConfig sim;
      sim.random_seed = trial * 200;
      ParallelConfig par;
      cm.allow_default_matrix = true;
      finalizeContactMatrices(cm, world, disease);
      InteractionManager im(world, cm, sim, par, &disease, nullptr);
      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);
      std::vector<PersonLocation> locs = {{0, 0, -1, 0, 255, 0},
                                          {1, 0, -1, 0, 255, 1}};
      im.processTransmissions(locs, 1.0, 8.0, nullptr);
      if (world.people[1].infection) infections_high++;
    }
  }

  // susc=0.08, contacts=1.0: lambda = 8*1/1*1*0.08 = 0.64, prob ≈ 47%
  // susc=1.0, contacts=1.0:  lambda = 8*1/1*1*1.0  = 8,    prob ≈ 99.97%
  CHECK(infections_high > infections_low);
  CHECK(infections_low > 150);  // ~47% → expect ~235
  CHECK(infections_low < 400);  // but not close to 100%
}

TEST_CASE("H3: Plague multi-mode: no infection during exposed stage") {
  // All infectiousness curves are zero during exposed → no transmission
  Disease disease = makePlagueLikeDisease();

  WorldState world = TestWorldFactory::createMinimalWorld(6, 1);
  world.venues[0].type_id = 0;
  world.people[0].infection = std::make_unique<Infection>(
      &disease, 0.0, &world.people[0], 42, nullptr, "office", 0);

  std::vector<PersonLocation> locs;
  for (int i = 0; i < 6; ++i)
    locs.push_back(
        {static_cast<PersonId>(i), 0, -1, 0, 255, static_cast<uint32_t>(i)});

  // At t=1, person 0 is in exposed stage (lasts 2 days)
  int new_inf = runTransmission(world, disease, locs, 1.0, 8.0);
  CHECK(new_inf == 0);
}

TEST_CASE("H4: Plague multi-mode: animal_bite infects during bacteraemia") {
  // During bacteraemia, animal_bite infectiousness = 1.0
  Disease disease = makePlagueLikeDisease();
  int total_infections = 0;
  const int N = 200;

  for (int trial = 0; trial < N; ++trial) {
    GlobalRNG::seed(trial * 300);
    WorldState world = TestWorldFactory::createMinimalWorld(6, 1);
    world.venues[0].type_id = 0;
    world.people[0].infection = std::make_unique<Infection>(
        &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);

    std::vector<PersonLocation> locs;
    for (int i = 0; i < 6; ++i)
      locs.push_back(
          {static_cast<PersonId>(i), 0, -1, 0, 255, static_cast<uint32_t>(i)});

    // t=6.5: bacteraemia stage (exposed=2d + buboes=4d → bacteraemia at t=6)
    // animal_bite: I=1.0, susc=0.08; respiratory: I=0.0
    // FoI per susceptible = 8*1.0/4 * 1.0 * 0.08 = 0.16
    // prob = 1 - exp(-0.16) ≈ 14.8%
    total_infections +=
        runTransmission(world, disease, locs, 6.5, 8.0, 1.0, trial * 300);
  }

  // 200 trials * 5 susceptible * 14.8% ≈ 148 expected infections
  CHECK(total_infections > 50);
}

TEST_CASE(
    "H5: Plague multi-mode: respiratory boosts infection during "
    "secondary_pneumonic") {
  // During secondary_pneumonic, both modes contribute:
  //   animal_bite: I=1.0 * 0.08 = 0.08
  //   respiratory: I=1.5 * 0.12 = 0.18
  // Total should be higher than animal_bite alone
  Disease disease = makePlagueLikeDisease();
  int infections_bacteraemia = 0;
  int infections_pneumonic = 0;
  const int N = 300;

  for (int trial = 0; trial < N; ++trial) {
    GlobalRNG::seed(trial * 400);
    {
      WorldState world = TestWorldFactory::createMinimalWorld(6, 1);
      world.venues[0].type_id = 0;
      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);
      std::vector<PersonLocation> locs;
      for (int i = 0; i < 6; ++i)
        locs.push_back({static_cast<PersonId>(i), 0, -1, 0, 255,
                        static_cast<uint32_t>(i)});
      // t=6.5: bacteraemia (only animal_bite active)
      infections_bacteraemia +=
          runTransmission(world, disease, locs, 6.5, 8.0, 1.0, trial * 400);
    }

    GlobalRNG::seed(trial * 400);
    {
      WorldState world = TestWorldFactory::createMinimalWorld(6, 1);
      world.venues[0].type_id = 0;
      world.people[0].infection = std::make_unique<Infection>(
          &disease, 0.0, &world.people[0], trial, nullptr, "office", 0);
      std::vector<PersonLocation> locs;
      for (int i = 0; i < 6; ++i)
        locs.push_back({static_cast<PersonId>(i), 0, -1, 0, 255,
                        static_cast<uint32_t>(i)});
      // t=7.5: secondary_pneumonic (both modes active)
      infections_pneumonic +=
          runTransmission(world, disease, locs, 7.5, 8.0, 1.0, trial * 400);
    }
  }

  // Secondary pneumonic should produce more infections because respiratory
  // mode adds to the force of infection
  CHECK(infections_pneumonic > infections_bacteraemia);
}
