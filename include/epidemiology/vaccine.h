#pragma once

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

namespace june {

// Simple struct to hold efficacy for an age range
struct AgeEfficacy {
  int min_age;
  int max_age;
  double efficacy;
};

// Represents one dose of a vaccine and how its protection changes over time
struct Dose {
  int number;
  double day_administered;

  // Milestones (relative to day_administered)
  // No protection at all until this many days after the dose. From here
  // protection rises linearly to its peak at days_to_effective. If the two
  // are equal, protection jumps straight to the peak.
  double days_to_onset = 0.0;
  double days_to_effective;
  double days_to_waning;
  double days_to_finished;  // When it's basically gone or at baseline

  // Efficacy values (0.0 to 1.0)
  // Map: Disease Name -> list of age-stratified efficacies
  std::unordered_map<std::string, std::vector<AgeEfficacy>> infection_efficacy;
  std::unordered_map<std::string, std::vector<AgeEfficacy>> symptom_efficacy;

  double waning_factor = 1.0;  // Level it wanes to relative to peak

  // Figures out the efficacy for a specific disease right now.
  // We use linear interpolation between the timing milestones.
  double getEfficacy(double current_time, const std::string& disease_name,
                     float age, bool for_symptoms = false) const {
    const auto& efficacies_map =
        for_symptoms ? symptom_efficacy : infection_efficacy;
    auto it = efficacies_map.find(disease_name);
    if (it == efficacies_map.end()) return 0.0;

    // Find applicable age bracket
    double peak_efficacy = 0.0;
    bool found_bracket = false;
    for (const auto& bracket : it->second) {
      if (age >= bracket.min_age && age <= bracket.max_age) {
        peak_efficacy = bracket.efficacy;
        found_bracket = true;
        break;
      }
    }

    if (!found_bracket) return 0.0;

    double time_since = current_time - day_administered;

    if (time_since < 0) return 0.0;

    // Nothing before onset. This also keeps the ramp below from dividing by
    // zero when onset equals days_to_effective.
    if (time_since < days_to_onset) return 0.0;

    // 1. Ramp up from onset to effective
    if (time_since < days_to_effective) {
      return ((time_since - days_to_onset) /
              (days_to_effective - days_to_onset)) *
             peak_efficacy;
    }

    // 2. Peak efficacy period
    if (time_since < days_to_waning) {
      return peak_efficacy;
    }

    // 3. Waning period
    if (time_since < days_to_finished) {
      double waning_start = days_to_waning;
      double waning_duration = days_to_finished - days_to_waning;
      double progress = (time_since - waning_start) / waning_duration;

      double final_efficacy = peak_efficacy * waning_factor;
      return peak_efficacy - (progress * (peak_efficacy - final_efficacy));
    }

    // 4. Post-waning baseline
    return peak_efficacy * waning_factor;
  }
};

// Keeps track of which vaccines and doses a person has actually received
class VaccineTrajectory {
 public:
  std::string vaccine_name;
  std::vector<Dose> doses;

  void addDose(const Dose& dose) {
    doses.push_back(dose);
    // Ensure doses are sorted by time
    std::sort(doses.begin(), doses.end(), [](const Dose& a, const Dose& b) {
      return a.day_administered < b.day_administered;
    });
  }

  // Works out the total protection from all doses combined.
  // When multiple doses are active, return the maximum efficacy across doses.
  double getEfficacy(double current_time, const std::string& disease_name,
                     float age, bool for_symptoms = false) const {
    double max_efficacy = 0.0;
    for (const auto& dose : doses) {
      max_efficacy = std::max(
          max_efficacy,
          dose.getEfficacy(current_time, disease_name, age, for_symptoms));
    }
    return max_efficacy;
  }

  int getNumDoses() const { return static_cast<int>(doses.size()); }

  double getLastDoseTime() const {
    if (doses.empty()) return -1.0;
    return doses.back().day_administered;
  }
};

}  // namespace june
