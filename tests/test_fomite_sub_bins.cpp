#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <memory>
#include <vector>

#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/fomite/fomite_sub_bins.h"

using namespace june;

// A fomite mode depositing at `deposition_by_symptom`, sub-binned at
// `sub_bin_time` hours (0 = one sub-bin per slot).
static TransmissionMode makeFomiteMode(
    const std::string& name, double sub_bin_time,
    std::vector<std::shared_ptr<InfectiousnessCurve>> deposition_by_symptom =
        {}) {
  TransmissionMode mode;
  mode.name = name;
  mode.type = TransmissionModeType::Fomite;
  FomiteConfig config;
  config.sub_bin_time = sub_bin_time;
  config.infectiousness_curve = std::make_shared<ConstantCurve>(1.0);
  config.deposition_by_symptom = std::move(deposition_by_symptom);
  mode.config = std::move(config);
  return mode;
}

static TransmissionMode makeDirectMode(const std::string& name) {
  TransmissionMode mode;
  mode.name = name;
  return mode;
}

static TransmissionMode makeCompartmentalUptakeMode(const std::string& name) {
  TransmissionMode mode;
  mode.name = name;
  mode.type = TransmissionModeType::CompartmentalUptake;
  mode.config = CompartmentalUptakeConfig{};
  return mode;
}

TEST_CASE("Schedule: sub-bin count per fomite mode from sub_bin_time") {
  TransmissionParams transmission;
  transmission.modes.push_back(makeFomiteMode("unset", 0.0));
  transmission.modes.push_back(makeFomiteMode("divisor", 2.0));
  transmission.modes.push_back(makeFomiteMode("wider_than_slot", 12.0));

  const double delta_hours = 8.0;
  FomiteSubBinSchedule schedule(transmission, delta_hours);

  REQUIRE(schedule.numModes() == 3);
  CHECK(schedule.subBinsPerMode() == std::vector<int>{1, 4, 1});
  CHECK(schedule.totalSubBins() == 6);
}

TEST_CASE("Schedule: fomite modes only, in mode-index order") {
  TransmissionParams transmission;
  transmission.modes.push_back(makeDirectMode("direct"));
  transmission.modes.push_back(makeFomiteMode("fomite_a", 0.0));
  transmission.modes.push_back(makeCompartmentalUptakeMode("uptake"));
  transmission.modes.push_back(makeFomiteMode("fomite_b", 0.0));

  FomiteSubBinSchedule schedule(transmission, 8.0);

  REQUIRE(schedule.numModes() == 2);
  CHECK(schedule.modes()[0].mode_index == 1);
  CHECK(schedule.modes()[1].mode_index == 3);
  CHECK(schedule.modes()[0].cfg ==
        &std::get<FomiteConfig>(transmission.modes[1].config));
  CHECK(schedule.modes()[1].cfg ==
        &std::get<FomiteConfig>(transmission.modes[3].config));
}

// Symptom ids: 0 = healthy, 1 = exposed, 2 = mild.
static constexpr uint16_t kHealthy = 0;
static constexpr uint16_t kExposed = 1;
static constexpr uint16_t kMild = 2;

static Disease makeDisease(std::vector<TransmissionMode> modes) {
  TransmissionParams transmission;
  transmission.mode = InfectiousnessMode::STAGE_DRIVEN;
  transmission.modes = std::move(modes);
  std::vector<SymptomTag> symptom_tags = {
      {"healthy", -1, kHealthy}, {"exposed", 0, kExposed}, {"mild", 1, kMild}};
  TrajectoryDefinition trajectory;
  trajectory.selection_key = "general";
  trajectory.severity = 1.0;
  trajectory.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  return Disease("TestFomiteSubBins", symptom_tags, {}, {trajectory}, {},
                 transmission);
}

// An Infection following exactly `transitions` (time, symptom id), with no
// sampled state and a cold symptom cache.
static std::unique_ptr<Infection> makeInfection(
    const Disease& disease, double infection_time,
    std::vector<std::pair<double, uint16_t>> transitions) {
  InfectionTrajectory trajectory;
  trajectory.infection_time = infection_time;
  trajectory.transitions = std::move(transitions);
  return Infection::fromCheckpoint(&disease, infection_time, trajectory, 1.0,
                                   1.0, 1.0, 0.0, /*last_checked_time=*/-1.0,
                                   kHealthy, infection_time);
}

TEST_CASE("Deposits: sub-bins against hand-computed values") {
  // Mode 0 (index 1): 3 sub-bins of 2 h over a 6 h slot. Exposed deposits a
  // constant 2 per hour; mild ramps 1 -> 3 over its first day in stage.
  // Mode 1 (index 2): one sub-bin, exposed only, constant 1.
  Disease disease = makeDisease({
      makeDirectMode("direct"),
      makeFomiteMode("fomite_sub_binned", 2.0,
                     {nullptr, std::make_shared<ConstantCurve>(2.0),
                      std::make_shared<LinearRampCurve>(1.0, 3.0, 1.0)}),
      makeFomiteMode("fomite_whole_slot", 0.0,
                     {nullptr, std::make_shared<ConstantCurve>(1.0), nullptr}),
  });
  FomiteSubBinSchedule schedule(disease.getTransmissionParams(), 6.0);
  REQUIRE(schedule.subBinsPerMode() == std::vector<int>{3, 1});

  // Slot [10, 10.25] d; sub-bins of 1/12 d start at 10, 10.0833, 10.1667.
  // Exposed -> mild at 10.1, inside sub-bin 1.
  const double transition_time = 10.1;
  auto infection =
      makeInfection(disease, 9.0, {{9.0, kExposed}, {transition_time, kMild}});

  std::vector<double> deposits;
  schedule.integrateDeposits(infection.get(), 10.0, deposits);

  // Deposit = 24 * integral over the sub-bin of the curve of the symptom held
  // at each instant; the untabulated curves integrate by midpoint.
  const double sub_bin_days = 0.25 / 3.0;
  // Sub-bin 0 is wholly exposed.
  const double exposed_sub_bin = 24.0 * 2.0 * sub_bin_days;  // 4
  // Sub-bin 1 straddles the transition: exposed until 10.1, then mild with
  // time in stage [0, 10.1667 - 10.1].
  const double straddle_exposed_days = transition_time - (10.0 + sub_bin_days);
  const double straddle_mild_days = 10.0 + 2 * sub_bin_days - transition_time;
  const double straddle_sub_bin =
      24.0 * 2.0 * straddle_exposed_days +
      24.0 * (1.0 + 2.0 * 0.5 * straddle_mild_days) * straddle_mild_days;
  // Sub-bin 2 is mild, time in stage [10.1667 - 10.1, 10.25 - 10.1].
  const double mild_midpoint =
      0.5 * ((10.0 + 2 * sub_bin_days - transition_time) +
             (10.25 - transition_time));
  const double mild_sub_bin =
      24.0 * (1.0 + 2.0 * mild_midpoint) * sub_bin_days;  // ~2.4333
  // Second mode: one sub-bin, deposits only while exposed (mild has no curve).
  const double whole_slot = 24.0 * 1.0 * (transition_time - 10.0);  // 2.4

  REQUIRE(deposits.size() == 4);
  CHECK(deposits[0] == doctest::Approx(exposed_sub_bin));
  CHECK(deposits[1] == doctest::Approx(straddle_sub_bin));
  CHECK(deposits[2] == doctest::Approx(mild_sub_bin));
  CHECK(deposits[3] == doctest::Approx(whole_slot));
}

TEST_CASE("Deposits: no infection gives an all-zero tail of full length") {
  Disease disease = makeDisease({
      makeFomiteMode("fomite_sub_binned", 2.0,
                     {std::make_shared<ConstantCurve>(1.0),
                      std::make_shared<ConstantCurve>(1.0),
                      std::make_shared<ConstantCurve>(1.0)}),
      makeFomiteMode("fomite_whole_slot", 0.0,
                     {std::make_shared<ConstantCurve>(1.0)}),
  });
  FomiteSubBinSchedule schedule(disease.getTransmissionParams(), 6.0);

  std::vector<double> deposits = {7.0};  // stale contents must not survive
  schedule.integrateDeposits(nullptr, 10.0, deposits);

  CHECK(deposits == std::vector<double>{0.0, 0.0, 0.0, 0.0});
}

TEST_CASE("Deposits: infected but not yet infectious still deposits") {
  // Deposition is keyed by symptom, so an incubating Person ("healthy" stage,
  // not infectious) deposits wherever that symptom's curve is nonzero.
  Disease disease = makeDisease({
      makeFomiteMode("fomite", 0.0,
                     {std::make_shared<ConstantCurve>(3.0), nullptr, nullptr}),
  });
  FomiteSubBinSchedule schedule(disease.getTransmissionParams(), 6.0);
  auto infection =
      makeInfection(disease, 9.0, {{9.0, kHealthy}, {11.0, kMild}});
  REQUIRE_FALSE(infection->isInfectious(10.0));

  std::vector<double> deposits;
  schedule.integrateDeposits(infection.get(), 10.0, deposits);

  REQUIRE(deposits.size() == 1);
  CHECK(deposits[0] == doctest::Approx(24.0 * 3.0 * 0.25));
}
