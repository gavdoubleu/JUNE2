#pragma once

// Shared fixtures and Disease factories for the MPI tests. Include only
// under USE_MPI.

#include <mpi.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/disease.h"
#include "epidemiology/infectiousness_curves.h"
#include "parallel/domain.h"
#include "parallel/domain_manager.h"

namespace june {

// ---------------------------------------------------------------------------
// Rank r owns geo unit r, venue r and person r.
//
// Each rank loads only its own venue into world.venues (mimicking the real
// distributed HDF5 load), but knows every venue's owner via the
// DomainManager's rank tables. Bypasses DomainManager::initialize().
// ---------------------------------------------------------------------------
struct RankPerPersonFixture {
  int rank;
  int size;

  WorldState world;
  Config config;
  std::unique_ptr<DomainManager> dm;

  RankPerPersonFixture() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    world.venue_type_names = {"household"};
    world.geo_level_names = {"MGU"};
    world.activity_names = {"residence", "work", "visiting", "none", "dead"};

    GeographicalUnit geo_unit;
    geo_unit.id = rank;
    geo_unit.name = "MGU_" + std::to_string(rank);
    geo_unit.level_id = 0;
    geo_unit.parent_id = -1;
    world.geo_units.push_back(geo_unit);

    Venue venue;
    venue.id = rank;
    venue.type_id = 0;
    venue.geo_unit_id = rank;
    venue.is_residence = true;
    world.venues.push_back(venue);

    Person& person = world.people.emplace_back();
    person.id = rank;
    person.age = 30.0f;
    person.sex = Sex::MALE;
    person.geo_unit_id = rank;

    world.buildIndices();

    config.parallel.partition_level = "MGU";
    config.parallel.geo_unit_chunk_size = 1000;

    dm = std::make_unique<DomainManager>(world, config);
    dm->setMPI(rank, size);
    dm->setMaxPersonId(size - 1);
    for (int owner = 0; owner < size; ++owner) {
      dm->setGeoUnitRank(owner, owner);
      dm->setPersonRank(owner, owner);
      dm->setVenueRank(owner, owner);
    }

    // Mirrors assignPeopleAndVenues
    Domain& domain = dm->getDomain();
    domain.addGeoUnit(rank);
    domain.resident_ids.push_back(rank);
    domain.resident_set.insert(rank);
    domain.local_venue_ids.push_back(rank);
    domain.local_venue_set.insert(rank);
  }
};

// The 2-rank case, with its ids named.
struct TwoRankFixture : RankPerPersonFixture {
  static constexpr PersonId PERSON_R0 = 0;
  static constexpr PersonId PERSON_R1 = 1;
  static constexpr VenueId VENUE_R0 = 0;
  static constexpr VenueId VENUE_R1 = 1;
};

// ---------------------------------------------------------------------------
// Diseases with symptom ids 0 = healthy, 1 = mild, where every infection is
// mild for 100 days.
// ---------------------------------------------------------------------------
inline std::vector<SymptomTag> healthyMildSymptomTags() {
  return {{"healthy", -1, 0}, {"mild", 1, 1}};
}

inline TrajectoryDefinition longMildTrajectory() {
  TrajectoryDefinition trajectory;
  trajectory.selection_key = "general";
  trajectory.severity = 1.0;
  trajectory.stages.push_back({"mild", {"constant", {{"value", 100.0}}}});
  return trajectory;
}

// One stage-driven mode: mild emits `mild_curve`, healthy nothing.
inline Disease makeStageDisease(
    const std::shared_ptr<InfectiousnessCurve>& mild_curve) {
  TransmissionParams params;
  params.mode = InfectiousnessMode::STAGE_DRIVEN;
  params.stage_curves["mild"] = mild_curve;
  params.symptom_id_curves = {nullptr, mild_curve};

  TransmissionMode default_mode;
  default_mode.name = "default";
  default_mode.symptom_curves = params.symptom_id_curves;
  params.modes.push_back(std::move(default_mode));

  return Disease("StageFlu", healthyMildSymptomTags(), {},
                 {longMildTrajectory()}, {}, params);
}

// Gamma(shape 2, rate 1) infectiousness profile scaled by
// `max_infectiousness`.
inline Disease makeTrajectoryDisease(double max_infectiousness) {
  TransmissionParams params;
  params.mode = InfectiousnessMode::TRAJECTORY_DRIVEN;
  params.type = "gamma";
  params.max_infectiousness = {"constant", {{"value", max_infectiousness}}};
  params.shape = {"constant", {{"value", 2.0}}};
  params.rate = {"constant", {{"value", 1.0}}};
  params.shift = {"constant", {{"value", 0.0}}};

  return Disease("TrajFlu", healthyMildSymptomTags(), {},
                 {longMildTrajectory()}, {}, params);
}

// A stage-driven mode whose mild stage emits a constant.
struct ConstantMildMode {
  std::string name;
  double mild_infectiousness;
};

// Two stage-driven modes, indexed 0 and 1 in the order given.
inline Disease makeTwoModeStageDisease(const ConstantMildMode& first_mode,
                                       const ConstantMildMode& second_mode) {
  TransmissionParams params;
  params.mode = InfectiousnessMode::STAGE_DRIVEN;
  for (const ConstantMildMode& mode_spec : {first_mode, second_mode}) {
    TransmissionMode mode;
    mode.name = mode_spec.name;
    mode.symptom_curves = {nullptr, std::make_shared<ConstantCurve>(
                                        mode_spec.mild_infectiousness)};
    params.modes.push_back(std::move(mode));
  }
  params.stage_curves["mild"] = params.modes[0].symptom_curves[1];
  params.symptom_id_curves = params.modes[0].symptom_curves;

  return Disease("TwoModeFlu", healthyMildSymptomTags(), {},
                 {longMildTrajectory()}, {}, params);
}

// ---------------------------------------------------------------------------
// Fomite disease: symptom ids 0 = healthy, 1 = exposed, 2 = mild. Healthy is
// infected but not infectious (incubation).
// ---------------------------------------------------------------------------
namespace fomite_flu {
constexpr uint16_t kHealthy = 0;
constexpr uint16_t kExposed = 1;
constexpr uint16_t kMild = 2;
}  // namespace fomite_flu

// Direct mode 0 plus fomite mode 1, sub-binned at `sub_bin_time` hours. Each
// symptom deposits on a different curve; mild ramps, so a deposit depends on
// the exact time in stage.
inline Disease makeFomiteDisease(double sub_bin_time) {
  TransmissionParams params;
  params.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto curve = std::make_shared<ConstantCurve>(5.0);
  params.symptom_id_curves = {nullptr, curve, curve};

  TransmissionMode direct;
  direct.name = "direct";
  direct.symptom_curves = params.symptom_id_curves;
  params.modes.push_back(std::move(direct));

  TransmissionMode fomite;
  fomite.name = "fomite";
  fomite.type = TransmissionModeType::Fomite;
  fomite.symptom_curves = {nullptr, nullptr, nullptr};
  FomiteConfig fomite_config;
  fomite_config.mode_index = 1;
  fomite_config.max_age = 2.0;
  fomite_config.sub_bin_time = sub_bin_time;
  fomite_config.infectiousness_curve = std::make_shared<ConstantCurve>(1.0);
  fomite_config.deposition_by_symptom = {
      std::make_shared<ConstantCurve>(0.7),
      std::make_shared<ConstantCurve>(2.0),
      std::make_shared<LinearRampCurve>(1.0, 3.0, 1.0)};
  fomite.config = std::move(fomite_config);
  params.modes.push_back(std::move(fomite));

  std::vector<SymptomTag> symptom_tags = {{"healthy", -1, fomite_flu::kHealthy},
                                          {"exposed", 0, fomite_flu::kExposed},
                                          {"mild", 1, fomite_flu::kMild}};
  return Disease("FomiteFlu", symptom_tags, {}, {longMildTrajectory()}, {},
                 params);
}

}  // namespace june
