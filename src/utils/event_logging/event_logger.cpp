#include "utils/event_logging/event_logger.h"

#include <iostream>

namespace june {

EventLogger::EventLogger() {
  // Reserve some space to avoid frequent reallocations
  infections_.reserve(10000);
  symptom_changes_.reserve(50000);
  deaths_.reserve(1000);
  hospital_admissions_.reserve(5000);
  icu_admissions_.reserve(2000);
  hospital_discharges_.reserve(5000);
  vaccinations_.reserve(5000);
  relationships_.reserve(10000);
  coordinated_encounters_.reserve(100000);
}

void EventLogger::logInfection(PersonId person_id, PersonId infector_id,
                               VenueId venue_id, double time,
                               uint8_t encounter_type_id,
                               const TransmissionRecord& transmission) {
  infections_.push_back(
      {.person_id = person_id,
       .infector_id = infector_id,
       .venue_id = venue_id,
       .time = time,
       .encounter_type_id = encounter_type_id,
       .transmission_mode_index = transmission.transmission_mode_index,
       .infector_symptom_id = transmission.infector_symptom_id,
       .source = transmission.source});
}

void EventLogger::logSymptomChange(PersonId person_id, VenueId venue_id,
                                   double time, uint8_t old_symptom_id,
                                   uint8_t new_symptom_id) {
  symptom_changes_.push_back(
      {person_id, venue_id, time, old_symptom_id, new_symptom_id});
}

void EventLogger::logDeath(PersonId person_id, VenueId venue_id, double time) {
  deaths_.push_back({person_id, venue_id, time});
}

void EventLogger::logHospitalAdmission(PersonId person_id, VenueId hospital_id,
                                       double time, const std::string& reason) {
  hospital_admissions_.push_back({person_id, hospital_id, time, reason});
}

void EventLogger::logICUAdmission(PersonId person_id, VenueId hospital_id,
                                  double time) {
  icu_admissions_.push_back({person_id, hospital_id, time});
}

void EventLogger::logHospitalDischarge(PersonId person_id, VenueId hospital_id,
                                       double time,
                                       const std::string& outcome) {
  hospital_discharges_.push_back({person_id, hospital_id, time, outcome});
}

void EventLogger::logVaccination(PersonId person_id,
                                 const std::string& vaccine_type,
                                 int dose_index, double time) {
  VaccinationEvent event;
  event.person_id = person_id;
  copyFixedField(event.vaccine_type, vaccine_type);
  event.dose_index = dose_index;
  event.time = time;
  vaccinations_.push_back(event);
}

void EventLogger::logRelationship(PersonId person_a, PersonId person_b,
                                  double time, double dissolution_time,
                                  const std::string& tie_tag) {
  RelationshipEvent event;
  event.person_a = person_a;
  event.person_b = person_b;
  event.time = time;
  event.dissolution_time = dissolution_time;
  copyFixedField(event.tie_tag, tie_tag);
  relationships_.push_back(event);
}

void EventLogger::logCoordinatedEncounter(PersonId person_a, PersonId person_b,
                                          double time,
                                          uint8_t encounter_type_id, int slot,
                                          uint64_t group_id) {
  coordinated_encounters_.push_back({.person_a = person_a,
                                     .person_b = person_b,
                                     .time = time,
                                     .encounter_type_id = encounter_type_id,
                                     .slot = slot,
                                     .group_id = group_id});
}

void EventLogger::logFollow(PersonId host, PersonId follower, double time,
                            uint8_t rule_id, int slot) {
  follows_.push_back({.host = host,
                      .follower = follower,
                      .time = time,
                      .rule_id = rule_id,
                      .slot = slot});
}

void EventLogger::logCampaignDay(int day, uint16_t campaign_id,
                                 uint16_t tier_id, int64_t candidates,
                                 int64_t doses_given, int64_t declined,
                                 int64_t allotment, int64_t started_total) {
  campaign_days_.push_back({.day = day,
                            .campaign_id = campaign_id,
                            .tier_id = tier_id,
                            .candidates = candidates,
                            .doses_given = doses_given,
                            .declined = declined,
                            .allotment = allotment,
                            .started_total = started_total});
}

void EventLogger::clear() {
  infections_.clear();
  symptom_changes_.clear();
  deaths_.clear();
  hospital_admissions_.clear();
  icu_admissions_.clear();
  hospital_discharges_.clear();
  vaccinations_.clear();
  relationships_.clear();
  coordinated_encounters_.clear();
  follows_.clear();
  campaign_days_.clear();
}

void EventLogger::printEncounterStats(
    const std::vector<std::string>& day_type_names,
    const std::vector<int>& day_type_counts) const {
  std::cout << "\n=== Encounter Statistics ===" << std::endl;
  size_t num_types = day_type_names.size();
  for (size_t i = 0; i < num_types; ++i) {
    const std::string& name = day_type_names[i];
    int days = (i < day_type_counts.size()) ? day_type_counts[i] : 0;

    size_t sched =
        (i < scheduled_encounters_.size()) ? scheduled_encounters_[i] : 0;
    size_t actual = (i < actual_encounters_.size()) ? actual_encounters_[i] : 0;

    std::cout << "Scheduled " << name << " encounters: " << sched;
    if (days > 0) std::cout << " (Daily Avg: " << sched / days << ")";
    std::cout << std::endl;

    std::cout << "Actual    " << name << " encounters: " << actual;
    if (days > 0) std::cout << " (Daily Avg: " << actual / days << ")";
    std::cout << std::endl;
  }
  std::cout << "============================" << std::endl;
}

void EventLogger::saveToHDF5WithLookups(
    const std::string& filename, const WorldState& world, const Config& config,
    const std::unordered_set<PersonId>& person_ids_filter) {
  event_writer::saveToHDF5WithLookups(*this, filename, world, config,
                                      person_ids_filter);
}

void EventLogger::flush(const std::string& filename, const Config& config,
                        const WorldState& world,
                        const std::unordered_set<PersonId>& person_ids_filter) {
  if (getTotalRecordCount() == 0) return;

  // Save current buffers and update lookups incrementally.
  saveToHDF5WithLookups(filename, world, config, person_ids_filter);

  // Clear buffers after successful flush
  clear();
}

std::unordered_set<PersonId> EventLogger::getInfectedPersonIds() const {
  std::unordered_set<PersonId> infected_ids;
  for (const auto& event : infections_) {
    infected_ids.insert(event.person_id);
  }
  return infected_ids;
}

}  // namespace june
