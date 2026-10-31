#pragma once

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/vaccination_allocation.h"
#include "epidemiology/vaccination_uptake.h"
#include "utils/event_logging/event_logger.h"
#include "utils/mpi_utils.h"
#include "utils/time_utils.h"
#include "vaccine.h"

namespace june {

// The most first doses a campaign may give: floor(target * population).
//
// The product is a floating point number, and 0.29 is not exactly 29/100, so
// 0.29 of 100 comes out as 28.999999999999996 and a plain floor would drop a
// whole person. The target and the multiply each carry at most half an ulp of
// error, so lifting the product by 4 * DBL_EPSILON (relative) covers both. It
// cannot reach the next integer unless the target has more digits than a
// double holds at that size.
inline int64_t firstDoseCap(double coverage_target, int64_t population) {
  if (!(coverage_target > 0.0 && coverage_target <= 1.0))
    throw std::runtime_error("vaccination: coverage_target " +
                             std::to_string(coverage_target) +
                             " must be above 0 and at most 1");
  if (population < 0)
    throw std::runtime_error("vaccination: population " +
                             std::to_string(population) +
                             " must not be negative");
  const double product = coverage_target * static_cast<double>(population);
  return static_cast<int64_t>(std::floor(product * (1.0 + 4.0 * DBL_EPSILON)));
}

// Spends each campaign's dose budget, one day at a time.
//
// At construction, which is collective, it settles who is in each tier and who
// accepts (settleUptake) and fixes the cap on first doses. Nothing else is kept
// between days that a restored checkpoint would not give back: the number of
// people who have started a course is counted again every day from the
// trajectories. Every rank runs the same code, in the same order, with the same
// number of collective calls.
class VaccinationManager {
 public:
  VaccinationManager(WorldState& world, const Config& config,
                     EventLogger* event_logger = nullptr)
      : world_(world), config_(config), event_logger_(event_logger) {
    if (!config.vaccination.enabled) return;
    const auto& configs = config.vaccination.campaigns;
    const std::tm sim_start = parseDate(config.simulation.start_date);
    const int64_t population =
        mpi_runtime::globalSum(static_cast<int64_t>(world.people.size()));
    campaigns_.reserve(configs.size());
    for (size_t c = 0; c < configs.size(); ++c) {
      const VaccinationCampaignConfig& cfg = configs[c];
      Campaign& campaign = campaigns_.emplace_back();
      campaign.config = &cfg;
      campaign.index = c;

      const auto vaccine = config.vaccination.vaccines.find(cfg.vaccine_type);
      if (vaccine == config.vaccination.vaccines.end())
        throw std::runtime_error("vaccination campaign '" + cfg.name +
                                 "': vaccine_type '" + cfg.vaccine_type +
                                 "' is not one of the vaccines");
      campaign.vaccine = &vaccine->second;
      if (cfg.tiers.empty())
        throw std::runtime_error("vaccination campaign '" + cfg.name +
                                 "' has no tiers");
      for (const auto& tier : cfg.tiers) {
        if (tier.dose >= static_cast<int>(vaccine->second.doses.size()))
          throw std::runtime_error(
              "vaccination campaign '" + cfg.name + "' tier '" + tier.name +
              "' gives dose " + std::to_string(tier.dose) + " but vaccine '" +
              cfg.vaccine_type + "' has " +
              std::to_string(vaccine->second.doses.size()) + " doses");
        campaign.starts_course.push_back(tier.dose == 0);
      }

      campaign.uptake =
          settleUptake(world, cfg, c, config.simulation.random_seed);
      if (cfg.coverage_target)
        campaign.cap = firstDoseCap(*cfg.coverage_target, population);

      // An absent start means campaign day 0 is the run's day 0. The end day
      // is included.
      if (!cfg.start_date.empty())
        campaign.start_day = daysBetween(sim_start, parseDate(cfg.start_date));
      if (!cfg.end_date.empty())
        campaign.end_day = daysBetween(sim_start, parseDate(cfg.end_date));
      campaign.candidates.resize(cfg.tiers.size());
    }
  }

  // This is called every simulation day, with the day index, to hand out
  // vaccines
  void update(int day) {
    const auto current_time = static_cast<double>(day);
    for (Campaign& campaign : campaigns_) {
      if (day < campaign.start_day || day > campaign.end_day) continue;
      runCampaign(campaign, day, current_time);
    }
  }

  // How far each campaign got, for the end of a run. Collective: every rank
  // calls it, and only rank 0 writes. A campaign with a coverage target says
  // whether it was reached, since a target that uptake, eligibility or the
  // budget made unreachable shows nowhere else but in the daily rows.
  void printCoverage(std::ostream& out) const {
    for (const Campaign& campaign : campaigns_) {
      const VaccinationCampaignConfig& cfg = *campaign.config;
      int64_t started = 0;
      for (const Person& person : world_.people)
        started += hasStarted(person, cfg.vaccine_type);
      started = mpi_runtime::globalSum(started);
      if (mpi_runtime::state().rank != 0) continue;

      out << "vaccination campaign '" << cfg.name << "': ";
      if (!campaign.cap) {
        out << started << " people started a course (no coverage target)\n";
        continue;
      }
      out << started << " of " << *campaign.cap
          << " first doses the coverage target allows";
      if (*campaign.cap > 0)
        out << " (" << std::fixed << std::setprecision(1)
            << 100.0 * static_cast<double>(started) /
                   static_cast<double>(*campaign.cap)
            << "%)";
      out << ", target "
          << (started >= *campaign.cap ? "reached" : "not reached") << "\n";
    }
  }

 private:
  // A candidate is kept as its place in the tier's lottery and its position in
  // world.people.
  using Candidate = std::pair<uint64_t, uint32_t>;

  struct Campaign {
    const VaccinationCampaignConfig* config = nullptr;
    const VaccineConfig* vaccine = nullptr;
    size_t index = 0;
    CampaignUptake uptake;
    std::optional<int64_t> cap;  // most first doses over the whole run
    int64_t start_day = 0;       // run day of campaign day 0
    int64_t end_day = std::numeric_limits<int64_t>::max();
    std::vector<bool> starts_course;  // per tier: does it give dose 0
    // Reused between days; cleared, not reallocated
    std::vector<std::vector<Candidate>> candidates;
  };

  // A course counts from its first dose, under the vaccine that gave it, and
  // the dead included.
  static bool hasStarted(const Person& person, const std::string& vaccine) {
    const VaccineTrajectory* trajectory = person.vaccine_trajectory.get();
    return trajectory && trajectory->getNumDoses() > 0 &&
           trajectory->vaccine_name == vaccine;
  }

  // The brand and gap rules for a person who holds `held` doses. Nobody with no
  // dose is held back by them. A course goes on with the vaccine it began
  // with, because coverage counts people under the first dose's vaccine; a
  // campaign that lists the vaccines a course may have begun with (a
  // booster, say) takes those instead.
  static bool nextDoseDue(const Person& person, int held,
                          const VaccinationCampaignConfig& cfg,
                          double current_time) {
    if (held == 0) return true;
    const VaccineTrajectory& trajectory = *person.vaccine_trajectory;
    const auto& accepted = cfg.courses_begun_with;
    if (accepted.empty()) {
      if (trajectory.vaccine_name != cfg.vaccine_type) return false;
    } else if (std::find(accepted.begin(), accepted.end(),
                         trajectory.vaccine_name) == accepted.end()) {
      return false;
    }
    const auto gap = static_cast<size_t>(held - 1);
    const double wait =
        gap < cfg.days_to_next_dose.size() ? cfg.days_to_next_dose[gap] : 0.0;
    return current_time >= trajectory.getLastDoseTime() + wait;
  }

  void runCampaign(Campaign& campaign, int64_t day, double current_time) {
    const VaccinationCampaignConfig& cfg = *campaign.config;
    const size_t tiers = cfg.tiers.size();
    const std::optional<int64_t> allotment =
        cfg.dose_budget ? std::optional(doseAllotment(cfg.dose_budget->doses,
                                                      cfg.dose_budget->per_days,
                                                      day - campaign.start_day))
                        : std::nullopt;

    // One pass over the local people: who has started a course (the dead
    // included, a course counts from its first dose) and who could be dosed
    // today.
    for (auto& list : campaign.candidates) list.clear();
    std::vector<int64_t> declined(tiers, 0);
    int64_t started = 0;
    for (size_t i = 0; i < world_.people.size(); ++i) {
      const Person& person = world_.people[i];
      const VaccineTrajectory* trajectory = person.vaccine_trajectory.get();
      const int held = trajectory ? trajectory->getNumDoses() : 0;
      started += hasStarted(person, cfg.vaccine_type);
      if (person.is_dead) continue;
      const uint64_t selected = campaign.uptake.selected[i];
      if (selected == 0) continue;
      // One dose a day, from the first campaign that reaches the person
      if (held > 0 && trajectory->getLastDoseTime() == current_time) continue;
      if (!nextDoseDue(person, held, cfg, current_time)) continue;
      for (size_t t = 0; t < tiers; ++t) {
        if (!((selected >> t) & 1) || cfg.tiers[t].dose != held) continue;
        // The first tier that fits decides; a decline does not move the
        // person on to a later tier.
        if ((campaign.uptake.accepted[i] >> t) & 1)
          campaign.candidates[t].emplace_back(
              lotteryKey(config_.simulation.random_seed,
                         static_cast<uint64_t>(person.id), campaign.index, t),
              static_cast<uint32_t>(i));
        else
          ++declined[t];
        break;
      }
    }

    started = mpi_runtime::globalSum(started);
    std::vector<int64_t> counts(tiers);
    for (size_t t = 0; t < tiers; ++t) {
      counts[t] = mpi_runtime::globalSum(
          static_cast<int64_t>(campaign.candidates[t].size()));
      declined[t] = mpi_runtime::globalSum(declined[t]);
    }

    std::optional<int64_t> cap_remaining;
    if (campaign.cap)
      cap_remaining = std::max<int64_t>(0, *campaign.cap - started);
    const std::vector<int64_t> take =
        fillTiers(allotment, counts, campaign.starts_course, cap_remaining);

    int64_t started_after = started;
    for (size_t t = 0; t < tiers; ++t) {
      if (campaign.starts_course[t]) started_after += take[t];
      if (take[t] > 0) serveTier(campaign, t, counts[t], take[t], current_time);
    }

    // Every count is a whole-world total, so only rank 0 writes the row
    if (event_logger_ && mpi_runtime::state().rank == 0) {
      for (size_t t = 0; t < tiers; ++t)
        event_logger_->logCampaignDay(
            static_cast<int>(day), static_cast<uint16_t>(campaign.index),
            static_cast<uint16_t>(t), counts[t], take[t], declined[t],
            allotment.value_or(kNoAllotment), started_after);
    }
  }

  // Doses `take` of the tier's `count` candidates over all ranks. When the
  // tier is not served in full, the ones with the lowest lottery keys go first.
  void serveTier(Campaign& campaign, size_t t, int64_t count, int64_t take,
                 double current_time) {
    std::vector<Candidate>& list = campaign.candidates[t];
    if (take < count) {
      std::sort(list.begin(), list.end());
      const auto local_at_or_below = [&](uint64_t key) {
        return std::upper_bound(
                   list.begin(), list.end(), key,
                   [](uint64_t k, const Candidate& c) { return k < c.first; }) -
               list.begin();
      };
      const std::optional<uint64_t> cutoff =
          lotteryCutoff(count, take, [&](uint64_t key) {
            return mpi_runtime::globalSum(local_at_or_below(key));
          });
      list.resize(static_cast<size_t>(local_at_or_below(*cutoff)));
    }
    for (const Candidate& candidate : list)
      administerVaccine(world_.people[candidate.second], campaign,
                        current_time);
  }

  void administerVaccine(Person& person, const Campaign& campaign,
                         double current_time) {
    const std::string& vaccine_type = campaign.config->vaccine_type;
    if (!person.vaccine_trajectory) {
      person.vaccine_trajectory = std::make_unique<VaccineTrajectory>();
      person.vaccine_trajectory->vaccine_name = vaccine_type;
    }
    const int dose_index = person.vaccine_trajectory->getNumDoses();
    // The tier's dose was checked against the vaccine at construction
    const DoseConfig& d_cfg = campaign.vaccine->doses[dose_index];

    Dose dose;
    dose.number = d_cfg.number;
    dose.day_administered = current_time;
    dose.days_to_onset = d_cfg.days_to_onset;
    dose.days_to_effective = d_cfg.days_to_effective;
    dose.days_to_waning = d_cfg.days_to_waning;
    dose.days_to_finished = d_cfg.days_to_finished;
    dose.waning_factor = d_cfg.waning_factor;
    dose.infection_efficacy = d_cfg.infection_efficacy;
    dose.symptom_efficacy = d_cfg.symptom_efficacy;

    person.vaccine_trajectory->addDose(dose);

    if (event_logger_)
      event_logger_->logVaccination(person.id, vaccine_type, dose_index,
                                    current_time);
  }

  WorldState& world_;
  const Config& config_;
  EventLogger* event_logger_;
  std::vector<Campaign> campaigns_;
};

}  // namespace june
