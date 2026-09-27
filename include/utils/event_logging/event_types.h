#pragma once

#include <string>
#include <vector>

#include "core/types.h"
#include "epidemiology/transmission/transmission_record.h"

namespace june {

// Shared fallback strings for unresolvable/inapplicable HDF5 lookup-table
// values. Single point of definition so every writer/merger call site stays
// in sync.
inline constexpr const char* kCouldNotResolve = "unresolved";
inline constexpr const char* kNotApplicable = "not applicable";

// =============================================================================
// Lightweight Population Record for Denominators
// =============================================================================

struct PopulationSummaryRecord {
  int person_id;               // 4 bytes
  uint8_t age_group;           // 1 byte (0-17 for 5-year groups)
  uint8_t sex_code;            // 1 byte (0=M, 1=F, 2=other)
  uint8_t schedule_type_code;  // 1 byte
  uint8_t reserved;            // 1 byte (padding)
  int geo_unit_id;             // 4 bytes
  uint8_t extra_codes[4];      // 4 bytes for configurable properties
                               // Total: 16 bytes
};

// =============================================================================
// Event Types (In-memory buffers)
// =============================================================================

struct InfectionEvent {
  PersonId person_id;
  PersonId infector_id;
  VenueId venue_id;
  double time;
  uint8_t encounter_type_id;
  uint8_t transmission_mode_index;
  uint8_t infector_symptom_id;
  InfectionSource source;
  uint8_t reserved = 0;
};

struct SymptomChangeEvent {
  PersonId person_id;
  VenueId venue_id;
  double time;
  uint8_t old_symptom_id;
  uint8_t new_symptom_id;
};

struct DeathEvent {
  PersonId person_id;
  VenueId venue_id;
  double time;
};

struct HospitalAdmissionEvent {
  PersonId person_id;
  VenueId hospital_id;
  double time;
  std::string reason;
};

struct ICUAdmissionEvent {
  PersonId person_id;
  VenueId hospital_id;
  double time;
};

struct HospitalDischargeEvent {
  PersonId person_id;
  VenueId hospital_id;
  double time;
  std::string outcome;
};

struct VaccinationEvent {
  PersonId person_id;
  char vaccine_type[64];
  int dose_index;
  double time;
};

struct RelationshipEvent {
  PersonId person_a;
  PersonId person_b;
  double time;
  double dissolution_time;
  char tie_tag[32];  // Generic tag: "ooe", etc.
};

struct CoordinatedEncounterEvent {
  PersonId person_a;  // Always the host of the encounter.
  PersonId person_b;  // A guest. The host fans one pair-row out per guest.
  double time;
  uint8_t encounter_type_id;
  int slot;
  // Stable id shared by every pair-row belonging to the same real group
  // encounter. Group rows on this to recover participant sets without
  // depending on the (person_a, time, slot, encounter_type_id) heuristic.
  uint64_t group_id;
};

// A follower bound to a host by a follow rule. Not a coordinated encounter:
// nothing was proposed, negotiated or accepted, so these live in their own
// dataset rather than sharing the encounter one.
struct FollowEvent {
  PersonId host;
  PersonId follower;
  double time;
  // Index into the follows list, named by the follow_rules registry.
  uint8_t rule_id;
  int slot;
};

// =============================================================================
// HDF5 Record Structures (C-compatible for H5::CompType)
// =============================================================================

namespace detail {

struct SymptomChangeRecord {
  int person_id;
  int venue_id;
  double time;
  uint8_t old_symptom_id;
  uint8_t new_symptom_id;
};

struct HospitalAdmissionRecord {
  int person_id;
  int hospital_id;
  double time;
  char reason[64];
};

struct HospitalDischargeRecord {
  int person_id;
  int hospital_id;
  double time;
  char outcome[64];
};

struct VaccinationRecord {
  int person_id;
  char vaccine_type[64];
  int dose_index;
  double time;
};

struct RelationshipRecord {
  int person_a;
  int person_b;
  double time;
  double dissolution_time;
  char tie_tag[32];
};

struct CoordinatedEncounterRecord {
  int person_a;
  int person_b;
  double time;
  uint8_t encounter_type_id;
  int slot;
  uint64_t group_id;
};

struct FollowRecord {
  int host;
  int follower;
  double time;
  uint8_t rule_id;
  int slot;
};

struct PersonRecord {
  int person_id;
  double age;
  char sex[16];
  int geo_unit_id;
  int is_dead;
  double death_time;
  char schedule_type[64];
  int num_activities;
  int num_residence_venues;
  int num_primary_activities;
  int num_leisure_venues;
  int num_medical_facilities;
};

struct PersonActivityRecord {
  int person_id;
  char activity_name[64];
  int venue_id;
  int subset_index;
  int activity_index;
};

struct VenueRecord {
  int venue_id;
  char name[128];
  char type[64];
  int geo_unit_id;
  int n_subsets;
};

}  // namespace detail

}  // namespace june
