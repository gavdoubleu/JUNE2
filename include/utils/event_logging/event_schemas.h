#pragma once

#include <H5Cpp.h>

#include <cstddef>

#include "utils/event_logging/event_types.h"

namespace june::event_schema {

inline H5::CompType infection() {
  H5::CompType type(sizeof(InfectionEvent));
  type.insertMember("person_id", HOFFSET(InfectionEvent, person_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("infector_id", HOFFSET(InfectionEvent, infector_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("venue_id", HOFFSET(InfectionEvent, venue_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("time", HOFFSET(InfectionEvent, time),
                    H5::PredType::NATIVE_DOUBLE);
  type.insertMember("encounter_type_id",
                    HOFFSET(InfectionEvent, encounter_type_id),
                    H5::PredType::NATIVE_UINT8);
  type.insertMember("transmission_mode_index",
                    HOFFSET(InfectionEvent, transmission_mode_index),
                    H5::PredType::NATIVE_UINT8);
  type.insertMember("infector_symptom_id",
                    HOFFSET(InfectionEvent, infector_symptom_id),
                    H5::PredType::NATIVE_UINT8);
  type.insertMember("source", HOFFSET(InfectionEvent, source),
                    H5::PredType::NATIVE_UINT8);
  return type;
}

inline H5::CompType symptomChange() {
  H5::CompType type(sizeof(SymptomChangeEvent));
  type.insertMember("person_id", HOFFSET(SymptomChangeEvent, person_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("venue_id", HOFFSET(SymptomChangeEvent, venue_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("time", HOFFSET(SymptomChangeEvent, time),
                    H5::PredType::NATIVE_DOUBLE);
  type.insertMember("old_symptom_id",
                    HOFFSET(SymptomChangeEvent, old_symptom_id),
                    H5::PredType::NATIVE_UINT8);
  type.insertMember("new_symptom_id",
                    HOFFSET(SymptomChangeEvent, new_symptom_id),
                    H5::PredType::NATIVE_UINT8);
  return type;
}

inline H5::CompType personVenueTime(std::size_t size, std::size_t person_offset,
                                    std::size_t venue_offset,
                                    std::size_t time_offset,
                                    const char* venue_name) {
  H5::CompType type(size);
  type.insertMember("person_id", person_offset, H5::PredType::NATIVE_INT);
  type.insertMember(venue_name, venue_offset, H5::PredType::NATIVE_INT);
  type.insertMember("time", time_offset, H5::PredType::NATIVE_DOUBLE);
  return type;
}

inline H5::CompType death() {
  return personVenueTime(sizeof(DeathEvent), HOFFSET(DeathEvent, person_id),
                         HOFFSET(DeathEvent, venue_id),
                         HOFFSET(DeathEvent, time), "venue_id");
}

inline H5::CompType hospitalAdmission() {
  H5::StrType string_type(H5::PredType::C_S1, 64);
  H5::CompType type(sizeof(detail::HospitalAdmissionRecord));
  type.insertMember("person_id",
                    HOFFSET(detail::HospitalAdmissionRecord, person_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("hospital_id",
                    HOFFSET(detail::HospitalAdmissionRecord, hospital_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("time", HOFFSET(detail::HospitalAdmissionRecord, time),
                    H5::PredType::NATIVE_DOUBLE);
  type.insertMember("reason", HOFFSET(detail::HospitalAdmissionRecord, reason),
                    string_type);
  return type;
}

inline H5::CompType icuAdmission() {
  return personVenueTime(sizeof(ICUAdmissionEvent),
                         HOFFSET(ICUAdmissionEvent, person_id),
                         HOFFSET(ICUAdmissionEvent, hospital_id),
                         HOFFSET(ICUAdmissionEvent, time), "hospital_id");
}

inline H5::CompType hospitalDischarge() {
  H5::StrType string_type(H5::PredType::C_S1, 64);
  H5::CompType type(sizeof(detail::HospitalDischargeRecord));
  type.insertMember("person_id",
                    HOFFSET(detail::HospitalDischargeRecord, person_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("hospital_id",
                    HOFFSET(detail::HospitalDischargeRecord, hospital_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("time", HOFFSET(detail::HospitalDischargeRecord, time),
                    H5::PredType::NATIVE_DOUBLE);
  type.insertMember("outcome",
                    HOFFSET(detail::HospitalDischargeRecord, outcome),
                    string_type);
  return type;
}

inline H5::CompType vaccination() {
  H5::StrType string_type(H5::PredType::C_S1, 64);
  H5::CompType type(sizeof(VaccinationEvent));
  type.insertMember("person_id", HOFFSET(VaccinationEvent, person_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("vaccine_type", HOFFSET(VaccinationEvent, vaccine_type),
                    string_type);
  type.insertMember("dose_index", HOFFSET(VaccinationEvent, dose_index),
                    H5::PredType::NATIVE_INT);
  type.insertMember("time", HOFFSET(VaccinationEvent, time),
                    H5::PredType::NATIVE_DOUBLE);
  return type;
}

inline H5::CompType relationship() {
  H5::StrType string_type(H5::PredType::C_S1, 32);
  H5::CompType type(sizeof(RelationshipEvent));
  type.insertMember("person_a", HOFFSET(RelationshipEvent, person_a),
                    H5::PredType::NATIVE_INT);
  type.insertMember("person_b", HOFFSET(RelationshipEvent, person_b),
                    H5::PredType::NATIVE_INT);
  type.insertMember("time", HOFFSET(RelationshipEvent, time),
                    H5::PredType::NATIVE_DOUBLE);
  type.insertMember("dissolution_time",
                    HOFFSET(RelationshipEvent, dissolution_time),
                    H5::PredType::NATIVE_DOUBLE);
  type.insertMember("tie_tag", HOFFSET(RelationshipEvent, tie_tag),
                    string_type);
  return type;
}

inline H5::CompType coordinatedEncounter() {
  H5::CompType type(sizeof(CoordinatedEncounterEvent));
  type.insertMember("person_a", HOFFSET(CoordinatedEncounterEvent, person_a),
                    H5::PredType::NATIVE_INT);
  type.insertMember("person_b", HOFFSET(CoordinatedEncounterEvent, person_b),
                    H5::PredType::NATIVE_INT);
  type.insertMember("time", HOFFSET(CoordinatedEncounterEvent, time),
                    H5::PredType::NATIVE_DOUBLE);
  type.insertMember("encounter_type_id",
                    HOFFSET(CoordinatedEncounterEvent, encounter_type_id),
                    H5::PredType::NATIVE_UINT8);
  type.insertMember("slot", HOFFSET(CoordinatedEncounterEvent, slot),
                    H5::PredType::NATIVE_INT);
  type.insertMember("group_id", HOFFSET(CoordinatedEncounterEvent, group_id),
                    H5::PredType::NATIVE_UINT64);
  return type;
}

inline H5::CompType follow() {
  H5::CompType type(sizeof(FollowEvent));
  type.insertMember("host", HOFFSET(FollowEvent, host),
                    H5::PredType::NATIVE_INT);
  type.insertMember("follower", HOFFSET(FollowEvent, follower),
                    H5::PredType::NATIVE_INT);
  type.insertMember("time", HOFFSET(FollowEvent, time),
                    H5::PredType::NATIVE_DOUBLE);
  type.insertMember("rule_id", HOFFSET(FollowEvent, rule_id),
                    H5::PredType::NATIVE_UINT8);
  type.insertMember("slot", HOFFSET(FollowEvent, slot),
                    H5::PredType::NATIVE_INT);
  return type;
}

inline H5::CompType campaignDay() {
  H5::CompType type(sizeof(CampaignDayEvent));
  type.insertMember("day", HOFFSET(CampaignDayEvent, day),
                    H5::PredType::NATIVE_INT);
  type.insertMember("campaign_id", HOFFSET(CampaignDayEvent, campaign_id),
                    H5::PredType::NATIVE_UINT16);
  type.insertMember("tier_id", HOFFSET(CampaignDayEvent, tier_id),
                    H5::PredType::NATIVE_UINT16);
  type.insertMember("candidates", HOFFSET(CampaignDayEvent, candidates),
                    H5::PredType::NATIVE_INT64);
  type.insertMember("doses_given", HOFFSET(CampaignDayEvent, doses_given),
                    H5::PredType::NATIVE_INT64);
  type.insertMember("declined", HOFFSET(CampaignDayEvent, declined),
                    H5::PredType::NATIVE_INT64);
  type.insertMember("allotment", HOFFSET(CampaignDayEvent, allotment),
                    H5::PredType::NATIVE_INT64);
  type.insertMember("started_total", HOFFSET(CampaignDayEvent, started_total),
                    H5::PredType::NATIVE_INT64);
  return type;
}

}  // namespace june::event_schema
