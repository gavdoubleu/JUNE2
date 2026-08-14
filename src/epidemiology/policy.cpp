#include "epidemiology/policy.h"

#include "epidemiology/disease.h"

namespace june {

PolicyManager::PolicyManager(WorldState& world)
    : world_(world), base_seed_(0) {}

void PolicyManager::addSymptomPolicy(const SymptomPolicy& policy) {
  symptom_policies_.push_back(policy);
}

void PolicyManager::addTemporalPolicy(const TemporalPolicy& policy) {
  temporal_policies_.push_back(policy);
}

bool PolicyManager::checkCompliance(double compliance_rate, PersonId person_id,
                                    uint32_t policy_index) {
  SplitMix64 rng(mix_seed(base_seed_, person_id, policy_index, 0xC0E1A9ULL));
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(rng) < compliance_rate;
}

std::optional<PersonLocation> PolicyManager::getReplacementLocation(
    const Person& person, const std::string& replacement_activity,
    int16_t replacement_activity_index) {
  int16_t act_idx = replacement_activity_index;
  if (act_idx < 0) {
    act_idx =
        static_cast<int16_t>(world_.getActivityIndex(replacement_activity));
  }

  auto venues = world_.getActivityVenues(person, act_idx);
  if (venues.empty()) {
    int16_t residence_idx =
        static_cast<int16_t>(world_.getActivityIndex("residence"));
    if (act_idx != residence_idx) {
      return getReplacementLocation(person, "residence", residence_idx);
    }
    return std::nullopt;
  }

  auto [venue_id, subset_idx] = venues[0];

  PersonLocation loc;
  loc.person_id = person.id;
  loc.venue_id = venue_id;
  loc.subset_index = subset_idx;
  loc.activity_index = act_idx;
  loc.encounter_type_id = 255;  // None

  return loc;
}

std::optional<PersonLocation> PolicyManager::getOverride(
    Person& person, int16_t scheduled_activity_index,
    VenueId scheduled_venue_id, SubsetIndex scheduled_subset_index,
    double current_time, int time_slot_index, const Person* partner) {
  // Venue type of the slot the person actually ends up in, resolved lazily so
  // ungated policies pay nothing. kUnknownVenueTypeId covers both "no venue"
  // and a cross-rank lookup miss; both read as absent at the gate.
  uint8_t scheduled_venue_type_id = kUnknownVenueTypeId;
  bool venue_type_resolved = false;
  auto scheduledVenueTypeId = [&]() {
    if (!venue_type_resolved) {
      if (scheduled_venue_id >= 0) {
        scheduled_venue_type_id = world_.getVenueTypeId(scheduled_venue_id);
      }
      venue_type_resolved = true;
    }
    return scheduled_venue_type_id;
  };

  // Priority 1: symptom-based policies
  if (person.infection != nullptr) {
    uint16_t current_symptom_id =
        person.infection->getTrajectory().getCurrentSymptomId(current_time);

    uint32_t symptom_mask = person.applicable_symptom_policy_mask;
    for (size_t i = 0; symptom_mask && i < symptom_policies_.size(); ++i) {
      if (!(symptom_mask & (1u << i))) continue;

      const auto& policy = symptom_policies_[i];

      bool is_triggered = policy.triggeredBy(current_symptom_id);

      // STICKY COMPLIANCE: if not triggered, clear participation/decision and
      // continue (with follow-up policy inheritance check)
      if (!is_triggered) {
        if (policy.follow_up_policy_index >= 0) {
          bool decision_made = (person.symptom_policy_decisions & (1u << i));
          bool participating =
              (person.active_symptom_policy_participation & (1u << i));

          if (decision_made) {
            bool should_inherit = false;
            if (participating && policy.inherit_compliance) {
              person.active_symptom_policy_participation |=
                  (1u << policy.follow_up_policy_index);
              should_inherit = true;
            }
            if (!participating && policy.inherit_refusal) {
              should_inherit = true;
            }

            if (should_inherit) {
              person.symptom_policy_decisions |=
                  (1u << policy.follow_up_policy_index);
            }
          }
        }

        // Restore paused hop state if this policy was responsible for freezing
        auto frozen_it = frozen_states_.find(person.id);
        if (frozen_it != frozen_states_.end() &&
            frozen_it->second.triggering_policy_index ==
                static_cast<uint8_t>(i)) {
          person.schedule_hop.restoreTargets(
              frozen_it->second.paused_hopped_schedule_id,
              frozen_it->second.paused_return_schedule_id);
          frozen_states_.erase(frozen_it);
        }

        person.active_symptom_policy_participation &= ~(1u << i);
        person.symptom_policy_decisions &= ~(1u << i);
        continue;
      }

      bool has_made_decision = (person.symptom_policy_decisions & (1u << i));
      if (!has_made_decision) {
        if (checkCompliance(policy.action.compliance_rate, person.id,
                            static_cast<uint32_t>(i))) {
          person.active_symptom_policy_participation |= (1u << i);
        }
        person.symptom_policy_decisions |= (1u << i);
      }

      bool is_participating =
          (person.active_symptom_policy_participation & (1u << i));
      if (!is_participating) continue;

      if (!policy.action.shouldOverride(scheduled_activity_index)) {
        continue;
      }

      if (policy.action.hasVenueGate() &&
          !policy.action.passesVenueGate(scheduledVenueTypeId())) {
        continue;
      }

      if (policy.action.isExempt(person, scheduled_activity_index, &world_,
                                 partner)) {
        continue;
      }

      // Freeze-hop path: pin person at current/last overnight venue. Only
      // applies to hopped persons.
      if (policy.action.replacement_schedule_idx >= 0 &&
          person.schedule_hop.isActive()) {
        ensureResidenceIndexCached();

        auto frozen_it = frozen_states_.find(person.id);
        if (frozen_it != frozen_states_.end() &&
            frozen_it->second.triggering_policy_index ==
                static_cast<uint8_t>(i)) {
          PersonLocation loc;
          loc.person_id = person.id;
          loc.venue_id = frozen_it->second.pin_venue_id;
          loc.subset_index = frozen_it->second.pin_subset_index;
          loc.activity_index = residence_act_idx_;
          loc.encounter_type_id = 255;
          return loc;
        }

        // First freeze: ActivityManager pre-resolves no_venue transit slots
        // to last overnight venue, so scheduled_venue_id is the best candidate.
        VenueId pin_venue = scheduled_venue_id;
        SubsetIndex pin_subset = scheduled_subset_index;
        if (pin_venue < 0) {
          auto home = world_.getActivityVenues(person, residence_act_idx_);
          if (!home.empty()) {
            pin_venue = home[0].first;
            pin_subset = home[0].second;
          }
        }

        int16_t saved_hop = person.schedule_hop.hopped_schedule_id;
        int16_t saved_return = person.schedule_hop.return_schedule_id;
        frozen_states_[person.id] =
            FrozenPersonState{static_cast<uint8_t>(i), saved_hop, saved_return,
                              pin_venue, pin_subset};
        person.schedule_hop.swapTarget(policy.action.replacement_schedule_idx);

        PersonLocation loc;
        loc.person_id = person.id;
        loc.venue_id = pin_venue;
        loc.subset_index = pin_subset;
        loc.activity_index = residence_act_idx_;
        loc.encounter_type_id = 255;
        return loc;
      }

      return getReplacementLocation(person, policy.action.replacement_activity,
                                    policy.action.replacement_activity_index);
    }
  } else {
    if (person.symptom_policy_decisions != 0) {
      person.resetPolicyState();
    }
  }

  // Priority 2: temporal policies (lockdowns, etc.)
  uint32_t temporal_mask = person.applicable_temporal_policy_mask;
  for (size_t i = 0; temporal_mask && i < temporal_policies_.size(); ++i) {
    if (!(temporal_mask & (1u << i))) continue;

    const auto& policy = temporal_policies_[i];

    bool is_active = policy.isActive(current_time);

    if (!is_active) {
      person.active_temporal_policy_participation &= ~(1u << i);
      person.temporal_policy_decisions &= ~(1u << i);
      continue;
    }

    bool has_made_decision = (person.temporal_policy_decisions & (1u << i));
    if (!has_made_decision) {
      if (checkCompliance(policy.action.compliance_rate, person.id,
                          static_cast<uint32_t>(i + 100))) {
        person.active_temporal_policy_participation |= (1u << i);
      }
      person.temporal_policy_decisions |= (1u << i);
    }

    bool is_participating =
        (person.active_temporal_policy_participation & (1u << i));

    if (!is_participating) continue;

    if (!policy.action.shouldOverride(scheduled_activity_index)) {
      continue;
    }

    if (policy.action.hasVenueGate() &&
        !policy.action.passesVenueGate(scheduledVenueTypeId())) {
      continue;
    }

    if (policy.action.isExempt(person, scheduled_activity_index, &world_,
                               partner)) {
      continue;
    }

    return getReplacementLocation(person, policy.action.replacement_activity,
                                  policy.action.replacement_activity_index);
  }

  return std::nullopt;
}

void PolicyManager::precomputePolicyApplicability(std::vector<Person>& people) {
  for (auto& person : people) {
    if (person.is_dead) continue;

    person.applicable_symptom_policy_mask = 0;
    person.applicable_temporal_policy_mask = 0;

    for (size_t i = 0; i < std::min(symptom_policies_.size(), size_t(32));
         ++i) {
      const auto& policy = symptom_policies_[i];
      if (policy.appliesTo(person, &world_)) {
        person.applicable_symptom_policy_mask |= (1u << i);
      }
    }

    for (size_t i = 0; i < std::min(temporal_policies_.size(), size_t(32));
         ++i) {
      const auto& policy = temporal_policies_[i];
      if (policy.appliesTo(person, &world_)) {
        person.applicable_temporal_policy_mask |= (1u << i);
      }
    }
  }
}

}  // namespace june
