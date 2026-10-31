#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "event_types.h"

namespace june {

class EventLogger;
namespace event_writer {
void saveToHDF5WithLookups(
    const EventLogger& logger, const std::string& filename,
    const WorldState& world, const Config& config,
    const std::unordered_set<PersonId>& person_ids_filter);
}

// =============================================================================
// EventLogger - Collects and saves epidemic events to HDF5
// =============================================================================

class EventLogger {
 public:
  EventLogger();

  // Log events
  void logInfection(PersonId person_id, PersonId infector_id, VenueId venue_id,
                    double time, uint8_t encounter_type_id,
                    const TransmissionRecord& transmission);
  void logSymptomChange(PersonId person_id, VenueId venue_id, double time,
                        uint8_t old_symptom_id, uint8_t new_symptom_id);
  void logDeath(PersonId person_id, VenueId venue_id, double time);

  // Log hospitalization events
  void logHospitalAdmission(PersonId person_id, VenueId hospital_id,
                            double time, const std::string& reason);
  void logICUAdmission(PersonId person_id, VenueId hospital_id, double time);
  void logHospitalDischarge(PersonId person_id, VenueId hospital_id,
                            double time, const std::string& outcome);
  void logVaccination(PersonId person_id, const std::string& vaccine_type,
                      int dose_index, double time);
  void logRelationship(PersonId person_a, PersonId person_b, double time,
                       double dissolution_time, const std::string& tie_tag);
  void logCoordinatedEncounter(PersonId person_a, PersonId person_b,
                               double time, uint8_t encounter_type_id, int slot,
                               uint64_t group_id);
  void logFollow(PersonId host, PersonId follower, double time, uint8_t rule_id,
                 int slot);
  // The counts are whole-world totals that every rank already holds, and the
  // logger does not know its rank. Per-rank files are concatenated when they
  // are merged, so a caller that logged on every rank would write each row once
  // per rank: log from rank 0 only.
  void logCampaignDay(int day, uint16_t campaign_id, uint16_t tier_id,
                      int64_t candidates, int64_t doses_given, int64_t declined,
                      int64_t allotment, int64_t started_total);

  // Save all events + lookup tables to HDF5 file (with config-based selective
  // saving)
  // `person_ids_filter` is the people whose lookup rows this write adds; an
  // empty one adds none.
  void saveToHDF5WithLookups(
      const std::string& filename, const WorldState& world,
      const Config& config,
      const std::unordered_set<PersonId>& person_ids_filter);

  // Flush current buffers to a specific file (incremental)
  void flush(const std::string& filename, const Config& config,
             const WorldState& world,
             const std::unordered_set<PersonId>& person_ids_filter);

  // Clear all events (useful for starting fresh)
  void clear();

  // Encounter stats per day type (day_type_idx = index into day_type_names)
  void logEncounterStats(int day_type_idx, bool is_actual, size_t count = 1) {
    auto& vec = is_actual ? actual_encounters_ : scheduled_encounters_;
    if (day_type_idx >= static_cast<int>(vec.size()))
      vec.resize(day_type_idx + 1, 0);
    vec[day_type_idx] += count;
  }
  void printEncounterStats(const std::vector<std::string>& day_type_names,
                           const std::vector<int>& day_type_counts) const;

  size_t getActualWeekdayEncounters() const {
    return actual_encounters_.empty() ? 0 : actual_encounters_[0];
  }
  size_t getActualWeekendEncounters() const {
    return actual_encounters_.size() > 1 ? actual_encounters_[1] : 0;
  }

  // Total record count across all event buffers
  size_t getTotalRecordCount() const {
    return infections_.size() + symptom_changes_.size() + deaths_.size() +
           hospital_admissions_.size() + icu_admissions_.size() +
           hospital_discharges_.size() + vaccinations_.size() +
           relationships_.size() + coordinated_encounters_.size() +
           follows_.size() + campaign_days_.size();
  }

  size_t getInfectionCount() const { return infections_.size(); }
  const std::vector<InfectionEvent>& getInfectionEvents() const {
    return infections_;
  }
  size_t getSymptomChangeCount() const { return symptom_changes_.size(); }
  size_t getDeathCount() const { return deaths_.size(); }
  size_t getHospitalAdmissionCount() const {
    return hospital_admissions_.size();
  }
  size_t getICUAdmissionCount() const { return icu_admissions_.size(); }
  size_t getHospitalDischargeCount() const {
    return hospital_discharges_.size();
  }
  size_t getVaccinationCount() const { return vaccinations_.size(); }
  const std::vector<CampaignDayEvent>& getCampaignDayEvents() const {
    return campaign_days_;
  }
  size_t getRelationshipCount() const { return relationships_.size(); }
  size_t getCoordinatedEncounterCount() const {
    return coordinated_encounters_.size();
  }

  // Helper: Get infected person IDs for the internal HDF5 writer.
  std::unordered_set<PersonId> getInfectedPersonIds() const;

 private:
  friend void event_writer::saveToHDF5WithLookups(
      const EventLogger& logger, const std::string& filename,
      const WorldState& world, const Config& config,
      const std::unordered_set<PersonId>& person_ids_filter);

  std::vector<InfectionEvent> infections_;
  std::vector<SymptomChangeEvent> symptom_changes_;
  std::vector<DeathEvent> deaths_;
  std::vector<HospitalAdmissionEvent> hospital_admissions_;
  std::vector<ICUAdmissionEvent> icu_admissions_;
  std::vector<HospitalDischargeEvent> hospital_discharges_;
  std::vector<VaccinationEvent> vaccinations_;
  std::vector<RelationshipEvent> relationships_;
  std::vector<CoordinatedEncounterEvent> coordinated_encounters_;
  std::vector<FollowEvent> follows_;
  std::vector<CampaignDayEvent> campaign_days_;

  std::vector<size_t> scheduled_encounters_;  // per day type index
  std::vector<size_t> actual_encounters_;     // per day type index
};

}  // namespace june
