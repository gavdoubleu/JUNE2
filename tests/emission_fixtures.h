#pragma once

// Hand-built Disease and Person shared by the Emission tests and the Visitor
// payload tests.

#include <memory>
#include <utility>
#include <vector>

#include "core/types.h"
#include "epidemiology/disease.h"

namespace emission_fixtures {

using namespace june;

// Symptom ids: 0 = healthy (infected, not infectious), 1 = exposed, 2 = mild.
inline constexpr uint16_t kHealthy = 0;
inline constexpr uint16_t kExposed = 1;
inline constexpr uint16_t kMild = 2;

// Two direct modes emitting different constants while mild, then a fomite
// mode sub-binned at `fomite_sub_bin_hours` depositing on every symptom.
inline Disease makeDisease(double fomite_sub_bin_hours) {
  TransmissionParams transmission;
  transmission.mode = InfectiousnessMode::STAGE_DRIVEN;

  TransmissionMode first_direct;
  first_direct.name = "first_direct";
  first_direct.symptom_curves = {
      nullptr, nullptr, std::make_shared<LinearRampCurve>(1.0, 3.0, 1.0)};
  transmission.modes.push_back(std::move(first_direct));

  TransmissionMode second_direct;
  second_direct.name = "second_direct";
  second_direct.symptom_curves = {nullptr, nullptr,
                                  std::make_shared<ConstantCurve>(0.3)};
  transmission.modes.push_back(std::move(second_direct));

  TransmissionMode fomite;
  fomite.name = "fomite";
  fomite.type = TransmissionModeType::Fomite;
  fomite.symptom_curves = {nullptr, nullptr, nullptr};
  FomiteConfig fomite_config;
  fomite_config.sub_bin_time = fomite_sub_bin_hours;
  fomite_config.infectiousness_curve = std::make_shared<ConstantCurve>(1.0);
  fomite_config.deposition_by_symptom = {
      std::make_shared<ConstantCurve>(0.7),
      std::make_shared<ConstantCurve>(2.0),
      std::make_shared<LinearRampCurve>(1.0, 3.0, 1.0)};
  fomite.config = std::move(fomite_config);
  transmission.modes.push_back(std::move(fomite));

  std::vector<SymptomTag> symptom_tags = {
      {"healthy", -1, kHealthy}, {"exposed", 0, kExposed}, {"mild", 1, kMild}};
  TrajectoryDefinition trajectory;
  trajectory.selection_key = "general";
  trajectory.severity = 1.0;
  trajectory.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  return Disease("TestEmission", symptom_tags, {}, {trajectory}, {},
                 transmission);
}

// A Person infected at `infection_time` following exactly `transitions`
// (time, symptom id), with a cold symptom cache.
inline Person makeInfectedPerson(
    const Disease& disease, double infection_time,
    std::vector<std::pair<double, uint16_t>> transitions) {
  InfectionTrajectory trajectory;
  trajectory.infection_time = infection_time;
  trajectory.transitions = std::move(transitions);
  Person person{};
  person.infection = Infection::fromCheckpoint(
      &disease, infection_time, trajectory, 1.0, 1.0, 1.0, 0.0,
      /*last_checked_time=*/-1.0, kHealthy, infection_time);
  return person;
}

}  // namespace emission_fixtures
