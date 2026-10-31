#pragma once

#include <algorithm>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

#include "core/types.h"
#include "epidemiology/transmission/transmission_record.h"

namespace june {

template <size_t N>
inline void copyFixedField(char (&destination)[N], std::string_view source) {
  static_assert(N > 0, "fixed-width fields must have room for a terminator");
  const size_t length = std::min(source.size(), N - 1);
  std::memset(destination, 0, N);
  std::memcpy(destination, source.data(), length);
}

// Records written to HDF5 as raw memory must not carry padding: a padding byte
// is whatever the allocator left behind, so the same events would be written
// differently from run to run. Where alignment leaves a gap the struct names it
// as a zero-initialised pad member at the same offset, so the schema does not
// move, and `packed` below fails the build if a gap is left unnamed.
// std::has_unique_object_representations is no use here: it is false for any
// struct holding a double.
template <typename T, typename... Members>
inline constexpr bool packed = sizeof(T) == (sizeof(Members) + ...);

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
  uint8_t pad_before_time_[4] = {};
  double time;
  uint8_t encounter_type_id;
  uint8_t transmission_mode_index;
  uint8_t infector_symptom_id;
  InfectionSource source;
  uint8_t reserved = 0;
  uint8_t pad_tail_[3] = {};
};
static_assert(
    packed<InfectionEvent, PersonId, PersonId, VenueId, uint8_t[4], double,
           uint8_t, uint8_t, uint8_t, InfectionSource, uint8_t, uint8_t[3]>);

struct SymptomChangeEvent {
  PersonId person_id;
  VenueId venue_id;
  double time;
  uint8_t old_symptom_id;
  uint8_t new_symptom_id;
  uint8_t pad_tail_[6] = {};
};
static_assert(packed<SymptomChangeEvent, PersonId, VenueId, double, uint8_t,
                     uint8_t, uint8_t[6]>);

struct DeathEvent {
  PersonId person_id;
  VenueId venue_id;
  double time;
};
static_assert(packed<DeathEvent, PersonId, VenueId, double>);

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
static_assert(packed<ICUAdmissionEvent, PersonId, VenueId, double>);

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
static_assert(packed<VaccinationEvent, PersonId, char[64], int, double>);

struct RelationshipEvent {
  PersonId person_a;
  PersonId person_b;
  double time;
  double dissolution_time;
  char tie_tag[32];  // Generic tag: "ooe", etc.
};
static_assert(
    packed<RelationshipEvent, PersonId, PersonId, double, double, char[32]>);

struct CoordinatedEncounterEvent {
  PersonId person_a;  // Always the host of the encounter.
  PersonId person_b;  // A guest. The host fans one pair-row out per guest.
  double time;
  uint8_t encounter_type_id;
  uint8_t pad_before_slot_[3] = {};
  int slot;
  // Stable id shared by every pair-row belonging to the same real group
  // encounter. Group rows on this to recover participant sets without
  // depending on the (person_a, time, slot, encounter_type_id) heuristic.
  uint64_t group_id;
};
static_assert(packed<CoordinatedEncounterEvent, PersonId, PersonId, double,
                     uint8_t, uint8_t[3], int, uint64_t>);

// A follower bound to a host by a follow rule. Not a coordinated encounter:
// nothing was proposed, negotiated or accepted, so these live in their own
// dataset rather than sharing the encounter one.
struct FollowEvent {
  PersonId host;
  PersonId follower;
  double time;
  // Index into the follows list, named by the follow_rules registry.
  uint8_t rule_id;
  uint8_t pad_before_slot_[3] = {};
  int slot;
};
static_assert(
    packed<FollowEvent, PersonId, PersonId, double, uint8_t, uint8_t[3], int>);

// What a vaccination campaign's dose budget did in one tier on one day. Every
// count is over the whole world, so every rank would write the same row.
struct CampaignDayEvent {
  int day;
  uint16_t campaign_id;  // index into the campaigns registry
  uint16_t tier_id;      // index within that campaign's tier list
  int64_t candidates;    // people in the tier who could be dosed today
  int64_t doses_given;   // doses given in the tier today
  int64_t declined;      // people in the tier who declined
  // The campaign's dose budget for the day, kNoAllotment when it has none.
  int64_t allotment;
  int64_t started_total;  // people who have started a course so far
};
static_assert(packed<CampaignDayEvent, int, uint16_t, uint16_t, int64_t,
                     int64_t, int64_t, int64_t, int64_t>);

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
  uint8_t pad_before_age_[4] = {};
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
  uint8_t pad_tail_[4] = {};
};
static_assert(packed<PersonRecord, int, uint8_t[4], double, char[16], int, int,
                     double, char[64], int, int, int, int, int, uint8_t[4]>);

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
