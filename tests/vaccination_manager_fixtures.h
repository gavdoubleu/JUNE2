#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/world_state.h"
#include "epidemiology/vaccine.h"
#include "vaccination_uptake_fixtures.h"

// Configs, campaigns and bookkeeping for the vaccination manager tests, on the
// hand-built worlds of the uptake fixtures.
namespace june::manager_fixture {

using namespace june::uptake_fixture;

constexpr uint64_t kSeed = 20261028;

inline VaccineConfig makeVaccine(const std::string& name, int doses = 3) {
  VaccineConfig vaccine;
  vaccine.name = name;
  for (int d = 0; d < doses; ++d) {
    DoseConfig dose;
    dose.number = d;
    dose.days_to_onset = 0.0;
    dose.days_to_effective = 14.0;
    dose.days_to_waning = 100.0;
    dose.days_to_finished = 200.0;
    vaccine.doses.push_back(dose);
  }
  return vaccine;
}

// A run that starts on a Monday, with vaccination on and two three-dose
// vaccines to choose from.
inline Config makeConfig(std::vector<VaccinationCampaignConfig> campaigns) {
  Config config;
  config.simulation.random_seed = kSeed;
  config.simulation.start_date = "2024-01-01";
  config.vaccination.enabled = true;
  for (const char* name : {"TestVax", "Other"})
    config.vaccination.vaccines[name] = makeVaccine(name);
  config.vaccination.campaigns = std::move(campaigns);
  return config;
}

inline VaccinationTier tierFor(std::string name, int dose,
                               std::vector<SelectionCriterion> group = {}) {
  VaccinationTier tier = makeTier(std::move(name), std::move(group));
  tier.dose = dose;
  return tier;
}

inline VaccinationCampaignConfig vaxCampaign(std::string name,
                                             std::vector<VaccinationTier> tiers,
                                             std::string vaccine = "TestVax") {
  VaccinationCampaignConfig campaign = makeCampaign(std::move(tiers));
  campaign.name = std::move(name);
  campaign.vaccine_type = std::move(vaccine);
  return campaign;
}

// Homes of two adults and a child: ids 3h and 3h + 1 are adults, 3h + 2 is the
// child.
inline std::vector<Home> homesOfThree(size_t count) {
  return std::vector<Home>(count, Home{{30, 40, 5}});
}

// A dose the manager did not give, to set up a world part way through a
// rollout.
inline void giveDose(Person& person, const std::string& vaccine, double day) {
  if (!person.vaccine_trajectory) {
    person.vaccine_trajectory = std::make_unique<VaccineTrajectory>();
    person.vaccine_trajectory->vaccine_name = vaccine;
  }
  Dose dose;
  dose.number = person.vaccine_trajectory->getNumDoses();
  dose.day_administered = day;
  dose.days_to_onset = 0.0;
  dose.days_to_effective = 14.0;
  dose.days_to_waning = 100.0;
  dose.days_to_finished = 200.0;
  person.vaccine_trajectory->addDose(dose);
}

// Who got a dose on `day` and which one (0 is the first). A person gets at most
// one dose a day, so the last dose being today's is the same as being dosed
// today.
using DosesGiven = std::map<PersonId, int>;

inline DosesGiven dosesGivenOn(const WorldState& world, double day) {
  DosesGiven given;
  for (const Person& person : world.people) {
    const VaccineTrajectory* trajectory = person.vaccine_trajectory.get();
    if (trajectory && trajectory->getNumDoses() > 0 &&
        trajectory->getLastDoseTime() == day)
      given[person.id] = trajectory->getNumDoses() - 1;
  }
  return given;
}

inline int64_t dosesOfIndex(const DosesGiven& given, int index) {
  int64_t count = 0;
  for (const auto& [id, dose_index] : given) count += dose_index == index;
  return count;
}

inline int64_t totalDoses(const WorldState& world) {
  int64_t total = 0;
  for (const Person& person : world.people)
    if (person.vaccine_trajectory)
      total += person.vaccine_trajectory->getNumDoses();
  return total;
}

}  // namespace june::manager_fixture
