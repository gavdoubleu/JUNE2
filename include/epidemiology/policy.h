#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "utils/deterministic_rng.h"
#include "utils/random.h"

namespace june {

// Forward declarations
class Disease;

// =============================================================================
// Activity Exemption - Specific conditions to exempt an activity
// =============================================================================

struct ActivityExemption {
  std::string activity_name;
  int16_t activity_index = -1;
  std::vector<SelectionCriterion> criteria;

  bool appliesTo(const Person& person, const WorldState* world = nullptr,
                 const Person* partner = nullptr) const {
    for (const auto& criterion : criteria) {
      if (!criterion.evaluate(person, world, partner)) {
        return false;
      }
    }
    return true;
  }

  void resolve(const WorldState& world) {
    activity_index =
        static_cast<int16_t>(world.getActivityIndex(activity_name));
    for (auto& criterion : criteria) {
      criterion.resolve(world);
    }
  }
};

// Direction of the scheduled-venue gate on a PolicyAction. The two directions
// are mutually exclusive, so one mask plus this flag covers both.
enum class VenueGateDirection : uint8_t {
  None = 0,        // no venue filter: override at any venue type
  RestrictTo = 1,  // override only at the listed venue types
  ExemptFrom = 2   // override everywhere except the listed venue types
};

struct PolicyAction {
  // Activities to override (empty = override all activities with "*")
  std::unordered_set<std::string> override_activities;
  uint64_t override_activity_mask = 0;  // BITMASK: support up to 64 activities
  bool override_all = false;

  // Venue types the override is restricted to (empty = any venue type).
  // ANDed with the activity mask: one Activity reaches many Venue types, so
  // this is what lets pubs close while groceries stay open.
  std::unordered_set<std::string> override_venue_types;

  // Venue types the override is exempt from (empty = no exemption). The
  // inverse direction: close everything except the listed types. Mutually
  // exclusive with override_venue_types, so both share one mask.
  std::unordered_set<std::string> exempt_venue_types;

  uint64_t venue_gate_mask = 0;  // BITMASK: support up to 64 venue types
  VenueGateDirection venue_gate_direction = VenueGateDirection::None;

  // Generic exemptions
  std::vector<ActivityExemption> exemptions;

  // What activity to do instead
  std::string replacement_activity;  // e.g., "residence", "medical_facility"
  int16_t replacement_activity_index = -1;

  // Optional: hop to a schedule instead of replacing the activity.
  // When set, getOverride triggers a schedule hop rather than calling
  // getReplacementLocation. Only effective for persons already on a hop.
  std::string replacement_schedule;
  int16_t replacement_schedule_idx = -1;

  // Compliance rate (0.0 = no one complies, 1.0 = everyone complies)
  double compliance_rate = 1.0;

  // Check if this action should override a given activity
  bool shouldOverride(const std::string& activity_name) const {
    if (override_all) return true;
    return override_activities.count(activity_name) > 0;
  }

  // Check by index
  bool shouldOverride(int16_t activity_index) const {
    if (override_all) return true;
    if (activity_index < 0 || activity_index >= 64) return false;
    return (override_activity_mask & (1ULL << activity_index));
  }

  // True when a venue filter is configured at all. Cheap guard so ungated
  // actions never pay for a venue-type lookup.
  bool hasVenueGate() const {
    return venue_gate_direction != VenueGateDirection::None;
  }

  // Check whether the venue the person actually ends up in passes the gate.
  // kUnknownVenueTypeId (255, also used for "no venue") is never in the mask,
  // so an absent venue uniformly reads as "type not listed": no override under
  // restrict-to, override under exempt-from.
  bool passesVenueGate(uint8_t venue_type_id) const {
    if (venue_gate_direction == VenueGateDirection::None) return true;
    const bool is_listed =
        venue_type_id < 64 && (venue_gate_mask & (1ULL << venue_type_id));
    return venue_gate_direction == VenueGateDirection::RestrictTo ? is_listed
                                                                  : !is_listed;
  }

  // Check if this action has an exemption for a given activity
  bool isExempt(const Person& person, int16_t activity_index,
                const WorldState* world = nullptr,
                const Person* partner = nullptr) const {
    for (const auto& exemption : exemptions) {
      if (exemption.activity_index == activity_index) {
        if (exemption.appliesTo(person, world, partner)) {
          return true;
        }
      }
    }
    return false;
  }

  // Resolve venue type names to a bitmask. Unlike SimulationConfig::resolve,
  // which silently ignores unknown venue types, an unknown name here throws: a
  // misspelt type would otherwise mean "nobody qualifies" and the policy would
  // quietly do nothing. Consequence: a policies.yaml naming 'pub' hard-fails on
  // a pub-less world.
  static uint64_t resolveVenueTypeMask(
      const WorldState& world,
      const std::unordered_set<std::string>& venue_type_names,
      const std::string& field_name, const std::string& policy_name) {
    uint64_t mask = 0;
    for (const auto& venue_type : venue_type_names) {
      int index = world.getVenueTypeIndex(venue_type);
      if (index < 0) {
        throw std::runtime_error("PolicyAction::resolve: policy '" +
                                 policy_name + "' lists unknown venue type '" +
                                 venue_type + "' in " + field_name + ".");
      }
      if (index >= 64) {
        throw std::runtime_error(
            "PolicyAction::resolve: policy '" + policy_name + "' venue type id " +
            std::to_string(index) + " ('" + venue_type + "') in " + field_name +
            " exceeds 64-bit mask width; promote venue_gate_mask to a wider "
            "bitset.");
      }
      mask |= (1ULL << index);
    }
    return mask;
  }

  void resolve(const WorldState& world, const std::string& policy_name = "") {
    venue_gate_mask = 0;
    venue_gate_direction = VenueGateDirection::None;
    if (!override_venue_types.empty() && !exempt_venue_types.empty()) {
      throw std::runtime_error(
          "PolicyAction::resolve: policy '" + policy_name +
          "' sets both override_venue_types and exempt_venue_types; the two "
          "directions are mutually exclusive.");
    }
    if (!override_venue_types.empty()) {
      venue_gate_mask = resolveVenueTypeMask(world, override_venue_types,
                                             "override_venue_types",
                                             policy_name);
      venue_gate_direction = VenueGateDirection::RestrictTo;
    } else if (!exempt_venue_types.empty()) {
      venue_gate_mask = resolveVenueTypeMask(world, exempt_venue_types,
                                             "exempt_venue_types", policy_name);
      venue_gate_direction = VenueGateDirection::ExemptFrom;
    }

    if (override_activities.empty() || override_activities.count("*") > 0) {
      override_all = true;
    } else {
      override_activity_mask = 0;
      for (const auto& act : override_activities) {
        int index = world.getActivityIndex(act);
        if (index >= 0 && index < 64) {
          override_activity_mask |= (1ULL << index);
        } else {
          std::cerr << "  [Policy Warning] Activity '" << act
                    << "' not found or index out of range for bitmask."
                    << std::endl;
        }
      }
    }

    // Resolve exemptions
    for (auto& exemption : exemptions) {
      exemption.resolve(world);
    }

    replacement_activity_index =
        static_cast<int16_t>(world.getActivityIndex(replacement_activity));

    if (!replacement_schedule.empty()) {
      replacement_schedule_idx = static_cast<int16_t>(
          world.getScheduleTypeIndex(replacement_schedule));
      if (replacement_schedule_idx < 0) {
        std::cerr << "  [Policy Warning] replacement_schedule '"
                  << replacement_schedule << "' not found." << std::endl;
      }
    }
  }
};

// =============================================================================
// Symptom-Based Policy - Override behavior based on disease symptoms
// =============================================================================

struct SymptomPolicy {
  std::string name;

  // Symptoms that trigger this policy
  std::vector<std::string> trigger_symptoms;
  uint32_t trigger_symptom_mask = 0;  // BITMASK: support up to 32 symptoms

  // Link to another policy to follow up if this one ends
  std::string follow_up_policy_name;
  int16_t follow_up_policy_index = -1;

  // Behavioral inheritance
  bool inherit_compliance = true;
  bool inherit_refusal = false;

  // What to do
  PolicyAction action;

  // Optional: selection criteria (only apply to certain people)
  std::vector<SelectionCriterion> applies_to;

  // Check if this policy applies to a person's current symptom
  bool triggeredBy(const std::string& symptom) const {
    return std::find(trigger_symptoms.begin(), trigger_symptoms.end(),
                     symptom) != trigger_symptoms.end();
  }

  // Check by symptom ID
  bool triggeredBy(uint16_t symptom_id) const {
    if (symptom_id >= 32) return false;
    return (trigger_symptom_mask & (1u << symptom_id));
  }

  // Check if policy applies to this person (based on selection criteria)
  bool appliesTo(const Person& person,
                 const WorldState* world = nullptr) const {
    // Empty criteria = applies to everyone
    if (applies_to.empty()) {
      return true;
    }

    // All criteria must match
    for (const auto& criterion : applies_to) {
      if (!criterion.evaluate(person, world)) {
        return false;
      }
    }
    return true;
  }

  void resolve(const WorldState& world, const Disease& disease) {
    for (auto& crit : applies_to) {
      crit.resolve(world);
    }

    // Intern action
    action.resolve(world, name);

    // Intern symptoms
    trigger_symptom_mask = 0;
    for (const auto& sym : trigger_symptoms) {
      uint16_t id = disease.getSymptomId(sym);
      if (id < 32) {
        trigger_symptom_mask |= (1u << id);
      }
    }
  }
};

// =============================================================================
// Temporal Policy - Override behavior during a time period (lockdowns, etc.)
// =============================================================================

struct TemporalPolicy {
  std::string name;

  // Time range (in simulation time, days from start)
  double start_time = 0.0;
  double end_time = -1.0;  // -1 = no end

  // What to do
  PolicyAction action;

  // Optional: selection criteria (only apply to certain people)
  std::vector<SelectionCriterion> applies_to;

  // Check if policy is active at given time
  bool isActive(double current_time) const {
    // Policy hasn't started yet
    if (current_time < start_time) {
      return false;
    }
    // Policy has ended (check if end_time is set, regardless of sign)
    // -1.0 means no end time (policy runs indefinitely)
    if (end_time != -1.0 && current_time > end_time) {
      return false;
    }
    // Policy is active
    return true;
  }

  // Check if policy applies to this person (based on selection criteria)
  bool appliesTo(const Person& person,
                 const WorldState* world = nullptr) const {
    // Empty criteria = applies to everyone
    if (applies_to.empty()) {
      return true;
    }

    // All criteria must match
    for (const auto& criterion : applies_to) {
      if (!criterion.evaluate(person, world)) {
        return false;
      }
    }
    return true;
  }

  void resolve(const WorldState& world) {
    for (auto& crit : applies_to) {
      crit.resolve(world);
    }
    action.resolve(world, name);
  }
};

// =============================================================================
// FrozenPersonState - Sparse storage for persons frozen by a policy hop
// =============================================================================

struct FrozenPersonState {
  uint8_t triggering_policy_index;
  int16_t paused_hopped_schedule_id;  // travel schedule to resume on recovery
  int16_t paused_return_schedule_id;  // saved schedule_hop.return_schedule_id
  VenueId pin_venue_id;
  SubsetIndex pin_subset_index;
  // schedule_hop.temp_slot_progress NOT saved: preserved automatically
  // (non-temporary hops never touch it, so it holds the correct
  // travel-schedule resume position)
};

// =============================================================================
// PolicyManager - Manages all policies and determines activity overrides
// =============================================================================

class PolicyManager {
 public:
  PolicyManager(WorldState& world);

  // Set base seed for deterministic RNG (MPI reproducibility)
  void setBaseSeed(uint64_t seed) { base_seed_ = seed; }

  // --- Checkpoint serialization ---
  // frozen_states_ pins the small set of persons mid policy-hop. It must be
  // saved/restored so an interrupted hop resumes correctly on restart.
  const std::unordered_map<PersonId, FrozenPersonState>& getFrozenStates()
      const {
    return frozen_states_;
  }
  void setFrozenStates(
      const std::unordered_map<PersonId, FrozenPersonState>& s) {
    frozen_states_ = s;
  }

  // Register policies
  void addSymptomPolicy(const SymptomPolicy& policy);
  void addTemporalPolicy(const TemporalPolicy& policy);

  // Get all policies (for inspection/debugging)
  const std::vector<SymptomPolicy>& getSymptomPolicies() const {
    return symptom_policies_;
  }
  const std::vector<TemporalPolicy>& getTemporalPolicies() const {
    return temporal_policies_;
  }

  // Main function: Check if a person's scheduled activity should be overridden
  // Returns std::nullopt if no override applies, otherwise returns the override
  // location
  std::optional<PersonLocation> getOverride(Person& person,
                                            int16_t scheduled_activity_index,
                                            VenueId scheduled_venue_id,
                                            SubsetIndex scheduled_subset_index,
                                            double current_time,
                                            int time_slot_index,
                                            const Person* partner = nullptr);

  // Clear all policies
  void clear() {
    symptom_policies_.clear();
    temporal_policies_.clear();
  }

  // Statistics
  size_t getSymptomPolicyCount() const { return symptom_policies_.size(); }
  size_t getTemporalPolicyCount() const { return temporal_policies_.size(); }

  // Precompute which policies can apply to each person (based on selection
  // criteria) This caches the results in person.applicable_*_policy_mask
  void precomputePolicyApplicability(std::vector<Person>& people);

  // Resolve all policy criteria and intern activities/symptoms
  void resolveAll(const Disease& disease) {
    for (auto& p : symptom_policies_) p.resolve(world_, disease);
    for (auto& p : temporal_policies_) p.resolve(world_);

    // Resolve follow-up policy indices
    for (auto& p : symptom_policies_) {
      if (!p.follow_up_policy_name.empty()) {
        for (size_t i = 0; i < symptom_policies_.size(); ++i) {
          if (symptom_policies_[i].name == p.follow_up_policy_name) {
            p.follow_up_policy_index = static_cast<int16_t>(i);
            break;
          }
        }
      }
    }
  }

 private:
  WorldState& world_;
  uint64_t base_seed_ = 0;

  std::vector<SymptomPolicy> symptom_policies_;
  std::vector<TemporalPolicy> temporal_policies_;

  // Sparse map: persons currently frozen by a policy-triggered schedule hop.
  // Only populated for the small minority of persons who are both travelling
  // and sick simultaneously. Avoids touching the Person struct.
  std::unordered_map<PersonId, FrozenPersonState> frozen_states_;

  // Cached residence activity index (resolved on first use)
  int16_t residence_act_idx_ = -1;
  void ensureResidenceIndexCached() {
    if (residence_act_idx_ < 0) {
      residence_act_idx_ =
          static_cast<int16_t>(world_.getActivityIndex("residence"));
    }
  }

  // Helper: Apply compliance rate (returns true if person complies)
  bool checkCompliance(double compliance_rate, PersonId person_id,
                       uint32_t policy_index);

  // Helper: Get replacement location for a given activity name
  std::optional<PersonLocation> getReplacementLocation(
      const Person& person, const std::string& replacement_activity,
      int16_t replacement_activity_index = -1);
};

}  // namespace june
