#pragma once

#include <algorithm>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/applicable_policies/applicable_policies.h"
#include "epidemiology/observed_metrics.h"
#include "epidemiology/policy_roster/policy_roster.h"
#include "epidemiology/policy_trigger.h"
#include "epidemiology/transmission_modifiers.h"

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
    return matchesAllCriteria(person, criteria, world, partner);
  }

  void resolve(const WorldState& world) {
    activity_index =
        static_cast<int16_t>(world.getActivityIndex(activity_name));
    for (auto& criterion : criteria) {
      criterion.resolveOrThrow(
          world, "policy exemption for activity '" + activity_name + "'");
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

// The type of the Venue a Person physically occupies for one slot: the key a
// venue-gated PolicyAction is decided on. Deliberately not the Pin Venue,
// which is only where a freeze_in_place override anchors the Person.
//
// Three states, not two. Absent means the Person occupies no Venue this slot
// (in transit, at a virtual encounter) and never passes the filter, in either
// direction. Known carries a resolved type. Deferred holds an unresolved
// VenueId so the world lookup happens only if a gated action is actually
// consulted — getOverride runs once per Person per slot, so an ungated run
// must pay nothing.
//
// It also carries the VenueId whenever the caller has one (known() and
// fromVenue() take it, venueId() returns it), because a policy triggered per
// venue is decided on the venue itself and not only on its type.
//
// A Venue that exists but cannot be typed on this rank is a defect, not a
// state: resolveAgainst throws rather than admitting kUnknownVenueTypeId,
// which is exactly the silent fail-open this type exists to remove.
enum class SlotVenueTypeState : uint8_t { Absent, Known, Deferred };

class SlotVenueType {
 public:
  // Occupies no Venue this slot.
  static SlotVenueType absent() {
    return SlotVenueType(SlotVenueTypeState::Absent, kUnknownVenueTypeId,
                         kInvalidVenueId);
  }

  // The caller already resolved the type. Handing over kUnknownVenueTypeId
  // while claiming to know it is the defect this mission exists to catch, so
  // it throws here, at the producer. A caller that also knows the venue passes
  // it, because a policy triggered per venue asks about that venue and not
  // only about its type.
  static SlotVenueType known(uint8_t type_id,
                             VenueId venue_id = kInvalidVenueId) {
    if (type_id == kUnknownVenueTypeId) {
      throw std::runtime_error(
          "SlotVenueType::known: caller claims a known venue type but passed "
          "kUnknownVenueTypeId (255); the slot venue type is unresolvable.");
    }
    return SlotVenueType(SlotVenueTypeState::Known, type_id, venue_id);
  }

  // Look the type up later, if anything asks. A negative id is absent.
  static SlotVenueType fromVenue(VenueId venue_id) {
    if (venue_id < 0) return absent();
    return SlotVenueType(SlotVenueTypeState::Deferred, kUnknownVenueTypeId,
                         venue_id);
  }

  // Collapses Deferred to Absent or Known. Throws naming the venue if the
  // world cannot type it. Absent and Known pass through unchanged.
  SlotVenueType resolveAgainst(const WorldState& world) const {
    if (state_ != SlotVenueTypeState::Deferred) return *this;
    const uint8_t type_id = world.getVenueTypeId(venue_id_);
    if (type_id == kUnknownVenueTypeId) {
      throw std::runtime_error(
          "SlotVenueType::resolveAgainst: venue " + std::to_string(venue_id_) +
          " exists but this rank cannot name its type; the slot venue type is "
          "unresolvable.");
    }
    return SlotVenueType(SlotVenueTypeState::Known, type_id, venue_id_);
  }

  bool isAbsent() const { return state_ == SlotVenueTypeState::Absent; }

  // The venue the person occupies this slot, or kInvalidVenueId when there is
  // none or the caller did not say. Unlike the Pin Venue this is where the
  // person is scheduled, so it is what a closed venue is checked against.
  VenueId venueId() const { return venue_id_; }

  // Precondition: Known. A sentinel must never survive to the gate, so an
  // unresolved value is a programming error rather than a return value.
  uint8_t typeId() const {
    if (state_ != SlotVenueTypeState::Known) {
      throw std::runtime_error(
          "SlotVenueType::typeId: value is not Known; resolveAgainst() must "
          "run before the venue gate is consulted.");
    }
    return type_id_;
  }

 private:
  SlotVenueType(SlotVenueTypeState state, uint8_t type_id, VenueId venue_id)
      : venue_id_(venue_id), state_(state), type_id_(type_id) {}

  VenueId venue_id_;
  SlotVenueTypeState state_;
  uint8_t type_id_;
};

struct PolicyAction {
  // False for a transmission-only policy: no implicit residence replacement.
  bool has_location_override = true;
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

  // Non-empty when a venue gate is configured but has not resolved. Set before
  // resolution is attempted and cleared only on success, so a swallowed
  // resolve() failure leaves the gate poisoned rather than silently absent.
  // Holds the message the gate throws when consulted.
  std::string venue_gate_error;

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
  bool shouldOverride(int16_t activity_index) const {
    if (override_all) return true;
    if (activity_index < 0 || activity_index >= 64) return false;
    return (override_activity_mask & (1ULL << activity_index));
  }

  // True when a venue filter is configured at all. Cheap guard so ungated
  // actions never pay for a venue-type lookup. A poisoned gate counts as
  // configured, so an unresolved one is still consulted and still throws.
  bool hasVenueGate() const {
    return venue_gate_direction != VenueGateDirection::None ||
           !venue_gate_error.empty();
  }

  // A configured-but-unresolved gate must never read as "this policy is not
  // here". Checked before the activity test, not just at the gate: a failed
  // resolve() also leaves override_activity_mask unbuilt, so the action would
  // otherwise fall out on shouldOverride and never reach the gate at all.
  void throwIfVenueGateUnresolved() const {
    if (!venue_gate_error.empty()) throw std::runtime_error(venue_gate_error);
  }

  // Check whether the type of the Venue the person physically occupies this
  // slot passes the gate. A Person occupying no Venue never passes, in either
  // direction: "not at a listed venue type" and "at no venue at all" are
  // different questions, and only the first is what exempt-from asks. An
  // unresolvable type never reaches here — SlotVenueType throws at the
  // producer instead.
  bool passesVenueGate(const SlotVenueType& slot_venue_type) const {
    if (!venue_gate_error.empty()) throw std::runtime_error(venue_gate_error);
    if (venue_gate_direction == VenueGateDirection::None) return true;
    if (slot_venue_type.isAbsent()) return false;
    const uint8_t type_id = slot_venue_type.typeId();  // throws if Deferred
    const bool is_listed =
        type_id < 64 && (venue_gate_mask & (1ULL << type_id));
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
            "PolicyAction::resolve: policy '" + policy_name +
            "' venue type id " + std::to_string(index) + " ('" + venue_type +
            "') in " + field_name +
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
    // Poison the gate before attempting to resolve it. resolveVenueTypeMask
    // throws on a misspelt venue type precisely so a policy cannot look
    // configured yet never fire, and that guarantee must not rest on no caller
    // anywhere ever wrapping resolve() in a catch. A swallowed failure now
    // leaves the gate configured-but-unresolved, which throws when consulted
    // rather than degrading to the silently ungated state.
    venue_gate_error =
        (override_venue_types.empty() && exempt_venue_types.empty())
            ? ""
            : "PolicyAction: policy '" + policy_name +
                  "' has a venue gate that never resolved; a resolve() failure "
                  "was swallowed by a caller.";
    if (!override_venue_types.empty() && !exempt_venue_types.empty()) {
      throw std::runtime_error(
          "PolicyAction::resolve: policy '" + policy_name +
          "' sets both override_venue_types and exempt_venue_types; the two "
          "directions are mutually exclusive.");
    }
    if (!override_venue_types.empty()) {
      venue_gate_mask = resolveVenueTypeMask(
          world, override_venue_types, "override_venue_types", policy_name);
      venue_gate_direction = VenueGateDirection::RestrictTo;
    } else if (!exempt_venue_types.empty()) {
      venue_gate_mask = resolveVenueTypeMask(world, exempt_venue_types,
                                             "exempt_venue_types", policy_name);
      venue_gate_direction = VenueGateDirection::ExemptFrom;
    }
    venue_gate_error.clear();  // resolved: only reachable if nothing threw

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
// ActiveWindow - the span of simulation time a policy is in force for
// =============================================================================

// Half-open [start_time, end_time): end_time is the first moment the policy is
// NOT in force. This mirrors the simulation's own window (`day < total_days_`,
// simulator.cpp), so an end_date of 2020-03-12 means "through 11 March" and
// adjacent windows abut exactly, with no slot where both are active.
// No end is an empty end_time, not a magic value: any number, -1 included, is
// a real day (an end_date one day before the simulation starts is day -1).
//
// A trigger can also hold a policy off: trigger_state is false while it does.
// Only a threshold trigger on a world scope metric sets it, from its closed
// state. A policy with no trigger, one triggered per venue, or a voluntary one
// leaves it true, so contains() is the one place that answers whether a policy
// is in force at all. A per venue or voluntary trigger then decides slot by
// slot where it acts.
struct ActiveWindow {
  double start_time = 0.0;         // days from simulation start
  std::optional<double> end_time;  // empty = no end
  bool trigger_state = true;

  bool contains(double current_time) const {
    if (!trigger_state) return false;
    if (current_time < start_time) return false;
    if (end_time && current_time >= *end_time) return false;
    return true;
  }
};

// =============================================================================
// Symptom-Based Policy - Override behavior based on disease symptoms
// =============================================================================

struct SymptomPolicy {
  std::string name;
  // Keys this policy's compliance and avoidance draws; policyDrawKey of its
  // kind and name, set by PolicyManager when the policy is added.
  uint64_t draw_key = 0;

  // Span of simulation time this policy is in force for
  ActiveWindow window;

  // Optional: ties activation to an observed metric
  std::optional<PolicyTrigger> trigger;

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
    return matchesAllCriteria(person, applies_to, world);
  }

  void resolve(const WorldState& world, const Disease& disease) {
    for (auto& crit : applies_to) {
      crit.resolveOrThrow(world, "symptom policy '" + name + "' applies_to");
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
  // Keys this policy's compliance and avoidance draws; policyDrawKey of its
  // kind and name, set by PolicyManager when the policy is added.
  uint64_t draw_key = 0;

  // Span of simulation time this policy is in force for
  ActiveWindow window;

  // Optional: ties activation to an observed metric
  std::optional<PolicyTrigger> trigger;

  // What to do
  PolicyAction action;

  // Optional: selection criteria (only apply to certain people)
  std::vector<SelectionCriterion> applies_to;

  // Check if policy applies to this person (based on selection criteria)
  bool appliesTo(const Person& person,
                 const WorldState* world = nullptr) const {
    return matchesAllCriteria(person, applies_to, world);
  }

  void resolve(const WorldState& world) {
    for (auto& crit : applies_to) {
      crit.resolveOrThrow(world, "temporal policy '" + name + "' applies_to");
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
  void setBaseSeed(uint64_t seed) { applicable_policies_.setBaseSeed(seed); }

  // Every read or write of a Person's Applicable Policy State.
  ApplicablePolicies& applicablePolicies() { return applicable_policies_; }
  const ApplicablePolicies& applicablePolicies() const {
    return applicable_policies_;
  }

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

  // Metrics declared in the policy file. A Policy trigger names the ones it
  // reads.
  ObservedMetrics& observedMetrics() { return observed_metrics_; }
  const ObservedMetrics& observedMetrics() const { return observed_metrics_; }

  // Mark the metrics the triggers read as active. Run once, after the metrics
  // resolve and before anything asks whether a metric is active.
  void requireTriggerMetrics();

  // Evaluate every trigger once. `day` is the day that just ended, and the
  // metrics must already be reduced for it. A change holds from the first slot
  // of the next day.
  void advanceTriggers(int64_t day);

  // Recompute the avoidance level of every voluntary trigger from the reduced
  // metrics. advanceTriggers does it each day; a checkpoint restore calls it
  // once after the metrics are put back.
  void rebuildResponseLevels();

  // For every voluntary Policy that shares a decision among the members of a
  // venue, check at startup that each person it can affect has exactly one
  // such venue. Collective: every rank must call it. Throws the same message
  // on every rank for the lowest failing person id. Needs the policies
  // resolved and their applicability computed.
  void checkDecisionUnits();

  // Trigger state for the checkpoint, symptom policies first, keyed by kind and
  // policy name. Restoring needs exactly the same keys.
  struct TriggerState {
    std::string key;
    TriggerCell cell;  // the world scope state
    MetricScope scope = MetricScope::World;
    // Scopes counted per unit: one cell per unit the trigger has touched, in
    // unit order. The cell above is unused then.
    std::vector<std::pair<VenueId, TriggerCell>> units;
    // A voluntary trigger keeps no state of its own; the flag only lets a
    // restore refuse a checkpoint that read the same Policy in the other form.
    bool voluntary = false;
  };
  std::vector<TriggerState> triggerStates() const;
  void restoreTriggerStates(const std::vector<TriggerState>& stored);

  // Register policies
  void addSymptomPolicy(const SymptomPolicy& policy);
  void addTemporalPolicy(const TemporalPolicy& policy);
  void addTransmissionEffect(PolicyTransmissionEffect effect) {
    transmission_effects_.push_back(std::move(effect));
  }

  void initializeTransmissionModifiers(const Disease& disease,
                                       double current_time);
  void refreshTransmissionModifiers(
      double current_time,
      const std::unordered_set<PersonId>& active_infections);
  double personModifier(const Person& person, size_t mode,
                        TransmissionEffectChannel channel) const {
    if (!transmission_modifiers_initialized_) return 1.0;
    return transmission_modifier_table_.get(person.transmission_modifier_set_id,
                                            mode, channel);
  }
  double venueModifier(const Venue* venue, size_t mode,
                       TransmissionEffectChannel channel) const {
    if (!venue || !transmission_modifiers_initialized_) return 1.0;
    return transmission_modifier_table_.get(venue->transmission_modifier_set_id,
                                            mode, channel);
  }
  const TransmissionModifierTable& transmissionModifierTable() const {
    return transmission_modifier_table_;
  }

  // Get all policies (for inspection/debugging)
  const std::vector<SymptomPolicy>& getSymptomPolicies() const {
    return symptom_policies_;
  }
  const std::vector<TemporalPolicy>& getTemporalPolicies() const {
    return temporal_policies_;
  }

  // Main function: Check if a person's scheduled activity should be overridden
  // Returns std::nullopt if no override applies, otherwise returns the override
  // location.
  //
  // pin_venue_id/pin_subset_index are the Pin Venue: where a freeze_in_place
  // override anchors the person. slot_venue_type is the Slot Venue Type: the
  // key a venue-gated action is decided on, and through its venue id the venue
  // a per venue trigger is checked against. The two are separate arguments
  // because callers legitimately differ on them — a traveller in transit is
  // pinned at their last overnight venue while occupying no venue at all.
  std::optional<PersonLocation> getOverride(
      Person& person, int16_t scheduled_activity_index, VenueId pin_venue_id,
      SubsetIndex pin_subset_index, SlotVenueType slot_venue_type,
      double current_time, int time_slot_index,
      const Person* partner = nullptr);

  // Would a policy remove this Person from `activity_index` this slot? A
  // question, not an instruction: no freeze is established or released, no
  // schedule hop swapped, no compliance decision latched. There is deliberately
  // no Pin Venue argument, because nothing is pinned.
  //
  // A per venue trigger never reaches a slot whose slot_venue_type carries no
  // venue, so pass the venue (SlotVenueType::known(type, venue)) when it is
  // known.
  //
  // Asked by the subsystems that only need the verdict — coordinated-encounter
  // eligibility and follow mirroring — so that reading the answer cannot commit
  // the consequence of it.
  bool suppressesParticipation(const Person& person, int16_t activity_index,
                               const SlotVenueType& slot_venue_type,
                               double current_time,
                               const Person* partner = nullptr) const;

  // End a freeze, but only if this policy is the one that established it —
  // another policy's freeze is not this policy's to lift.
  void releaseFreeze(Person& person, size_t policy_index) {
    auto frozen_it = frozen_states_.find(person.id);
    if (frozen_it == frozen_states_.end()) return;
    if (frozen_it->second.triggering_policy_index !=
        static_cast<uint8_t>(policy_index)) {
      return;
    }
    releaseAnyFreeze(person);
  }

  // End a freeze whichever policy established it, putting the person back on
  // the travel schedule they were pinned off. A person nobody froze is left
  // alone. The only place frozen_states_ entries are erased.
  void releaseAnyFreeze(Person& person) {
    auto frozen_it = frozen_states_.find(person.id);
    if (frozen_it == frozen_states_.end()) return;
    person.schedule_hop.restoreTargets(
        frozen_it->second.paused_hopped_schedule_id,
        frozen_it->second.paused_return_schedule_id);
    frozen_states_.erase(frozen_it);
  }

  // Statistics
  size_t getSymptomPolicyCount() const { return symptom_policies_.size(); }
  size_t getTemporalPolicyCount() const { return temporal_policies_.size(); }

  // Precompute which policies can apply to each person (based on selection
  // criteria) This caches the results in person.applicable_policy_state
  void precomputePolicyApplicability(std::vector<Person>& people);

  // Resolve all policy criteria and intern activities/symptoms
  void resolveAll(const Disease& disease) {
    rejectPoliciesBeyondMaskWidth(symptom_policies_, "symptom");
    rejectPoliciesBeyondMaskWidth(temporal_policies_, "temporal");
    // A policy is its name within its kind (checkpoints match policies by it).
    rejectDuplicateNames(
        makePolicyRoster(symptom_policies_, temporal_policies_));
    rejectDrawKeyCollisions();

    for (auto& p : symptom_policies_) p.resolve(world_, disease);
    for (auto& p : temporal_policies_) p.resolve(world_);
    resolveResponses();
    resolveTransmissionEffects(disease);

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

    // A follow-up inherits a latched decision, and a voluntary policy keeps
    // none.
    for (const auto& p : symptom_policies_) {
      if (p.follow_up_policy_index < 0) continue;
      const auto& target = symptom_policies_[p.follow_up_policy_index];
      if (target.trigger && target.trigger->voluntary()) {
        throw std::runtime_error(
            "symptom_policy '" + p.name + "' names the voluntary policy '" +
            target.name +
            "' as its follow_up_policy, which is not built yet: a voluntary "
            "policy keeps no decision for a follow-up to inherit.");
      }
    }
  }

 private:
  WorldState& world_;

  std::vector<SymptomPolicy> symptom_policies_;
  std::vector<TemporalPolicy> temporal_policies_;
  std::vector<PolicyTransmissionEffect> transmission_effects_;
  ApplicablePolicies applicable_policies_;
  ObservedMetrics observed_metrics_;
  TransmissionModifierTable transmission_modifier_table_;
  bool transmission_modifiers_initialized_ = false;
  size_t transmission_mode_count_ = 0;
  // Policies of each kind whose window holds the current time; a change marks
  // a window edge.
  struct PolicyWindowBits {
    PolicyMask temporal;
    PolicyMask symptom;
    bool operator==(const PolicyWindowBits&) const = default;
  };
  PolicyWindowBits last_policy_window_bits_;
  std::unordered_map<PersonId, uint16_t> last_symptom_ids_;

  void resolveTransmissionEffects(const Disease& disease);
  // Look up the activities a voluntary trigger names and start its level at the
  // base. Needs the metrics resolved.
  void resolveResponses();
  PolicyWindowBits policyWindowBits(double current_time) const;
  void rebuildTransmissionModifiers(double current_time);
  void rebuildPersonTransmissionModifier(Person& person, double current_time);
  void rebuildVenueTransmissionModifier(Venue& venue, double current_time);

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

  // Global policy indices are bit positions in a kMaxPoliciesPerKind-bit
  // PolicyMask, so a policy past the cap could never be marked applicable. Name
  // the first such policy rather than drop it. One Person's set is capped
  // lower, at kMaxApplicablePoliciesPerKind (see ApplicablePolicySets::intern).
  template <typename PolicyVector>
  static void rejectPoliciesBeyondMaskWidth(const PolicyVector& policies,
                                            const std::string& kind) {
    if (policies.size() <= kMaxPoliciesPerKind) return;
    throw std::runtime_error(
        "PolicyManager: " + std::to_string(policies.size()) + " " + kind +
        " policies exceeds the limit of " +
        std::to_string(kMaxPoliciesPerKind) + "; '" +
        policies[kMaxPoliciesPerKind].name + "' would never take effect");
  }

  // Two policies sharing a draw key would draw alike for every Person. Names are
  // unique per kind and the kind is hashed in, so only a hash collision does
  // this; refuse it, naming both, rather than run with correlated draws.
  void rejectDrawKeyCollisions() const;

  // The part of a policy decision that is a pure question: does this action
  // reach this Person's activity, at this slot's venue type, unexempted?
  // Shared by getOverride and suppressesParticipation so the gate ordering and
  // the unresolved-gate throw exist in exactly one place.
  //
  // Takes a resolver rather than a resolved SlotVenueType so an ungated action
  // never pays for the venue-type lookup — and never throws on a venue this
  // rank cannot type when no policy was asking about venues in the first place.
  //
  // A Policy triggered per venue (slot_trigger non-null, threshold form)
  // reaches only the venues its metric has closed. That is one bit lookup, made
  // after the activity gate and before the venue type gate and the exemptions,
  // and it applies to everyone scheduled at a closed venue who is not exempt.
  // A voluntary trigger (draw_key names the Policy in the draw) asks instead
  // whether this person's decision unit avoids the venue.
  template <typename SlotVenueTypeResolver>
  bool actionApplies(const PolicyAction& action, const Person& person,
                     int16_t activity_index,
                     SlotVenueTypeResolver&& resolve_slot_venue_type,
                     const Person* partner, const PolicyTrigger* slot_trigger,
                     VenueId slot_venue_id, uint64_t draw_key) const {
    if (!action.has_location_override) return false;
    action.throwIfVenueGateUnresolved();

    if (!action.shouldOverride(activity_index)) return false;

    if (slot_trigger) {
      if (slot_trigger->voluntary()) {
        if (!avoids(*slot_trigger, person, slot_venue_id, draw_key))
          return false;
      } else if (!unitClosedAt(*slot_trigger, slot_venue_id)) {
        return false;
      }
    }

    if (action.hasVenueGate() &&
        !action.passesVenueGate(resolve_slot_venue_type())) {
      return false;
    }

    if (action.isExempt(person, activity_index, &world_, partner)) return false;

    return true;
  }

  // Is the unit this venue counts toward closed by this trigger?
  bool unitClosedAt(const PolicyTrigger& trigger, VenueId venue) const;

  // Does this person's decision unit avoid the venue they are scheduled at?
  // Pure: the level is the one derived at the last day boundary and the draw is
  // fixed by the base seed, the decision unit and the Policy, so asking changes
  // nothing and the answer never depends on the process count. A trigger
  // counted per unit only reaches people the metric counts at that venue's
  // unit. Throws when a shared decision unit finds the person with other than
  // exactly one such venue; checkDecisionUnits rules that out at startup.
  bool avoids(const PolicyTrigger& trigger, const Person& person,
              VenueId slot_venue_id, uint64_t draw_key) const;

  // Helper: Get replacement location for a given activity name
  std::optional<PersonLocation> getReplacementLocation(
      const Person& person, const std::string& replacement_activity,
      int16_t replacement_activity_index = -1);
};

}  // namespace june
