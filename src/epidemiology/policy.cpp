#include "epidemiology/policy.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_map>

#include "epidemiology/applicable_policies/applicable_policy_walk.h"
#include "epidemiology/disease.h"
#include "utils/deterministic_rng.h"
#include "utils/filtering.h"
#include "utils/mpi_utils.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace june {

namespace {

class LazySlotVenueTypeResolver {
 public:
  LazySlotVenueTypeResolver(const SlotVenueType& slot_venue_type,
                            const WorldState& world)
      : slot_venue_type_(slot_venue_type), world_(world) {}

  const SlotVenueType& operator()() {
    if (!resolved_) {
      resolved_slot_venue_type_ = slot_venue_type_.resolveAgainst(world_);
      resolved_ = true;
    }
    return resolved_slot_venue_type_;
  }

 private:
  const SlotVenueType& slot_venue_type_;
  const WorldState& world_;
  SlotVenueType resolved_slot_venue_type_ = SlotVenueType::absent();
  bool resolved_ = false;
};

// Tag for the avoidance draw, so it never repeats a compliance draw.
constexpr uint64_t kAvoidanceDrawTag = 0xA701D4ULL;

// The trigger a slot has to ask about, or null: a policy that closes venues one
// at a time, or a voluntary one.
template <typename Policy>
const PolicyTrigger* slotTrigger(const Policy& policy) {
  return policy.trigger && policy.trigger->decidesPerSlot() ? &*policy.trigger
                                                            : nullptr;
}

// Visit every triggered policy, symptom policies first. Works for the const
// and non-const vectors alike.
template <typename Symptom, typename Temporal, typename Visit>
void eachTriggered(Symptom& symptom_policies, Temporal& temporal_policies,
                   Visit&& visit) {
  for (auto& policy : symptom_policies)
    if (policy.trigger) visit(std::string("symptom"), policy);
  for (auto& policy : temporal_policies)
    if (policy.trigger) visit(std::string("temporal"), policy);
}

}  // namespace

PolicyManager::PolicyManager(WorldState& world)
    : world_(world) {}

void PolicyManager::addSymptomPolicy(const SymptomPolicy& policy) {
  symptom_policies_.push_back(policy);
  applicable_policies_.setLoadedCounts(symptom_policies_.size(),
                                       temporal_policies_.size());
  // A policy with a world scope threshold trigger starts off until the trigger
  // closes. One triggered per venue, or a voluntary one, is in force from the
  // start: the closed venues, or the avoidance level, decide where it acts.
  if (policy.trigger)
    symptom_policies_.back().window.trigger_state =
        policy.trigger->holdsPolicy();
}

void PolicyManager::addTemporalPolicy(const TemporalPolicy& policy) {
  temporal_policies_.push_back(policy);
  applicable_policies_.setLoadedCounts(symptom_policies_.size(),
                                       temporal_policies_.size());
  if (policy.trigger)
    temporal_policies_.back().window.trigger_state =
        policy.trigger->holdsPolicy();
}

void PolicyManager::requireTriggerMetrics() {
  eachTriggered(symptom_policies_, temporal_policies_,
                [&](const std::string&, const auto& policy) {
                  if (!policy.trigger->voluntary()) {
                    observed_metrics_.require(policy.trigger->metric);
                    return;
                  }
                  for (const auto& term : policy.trigger->terms)
                    observed_metrics_.require(term.metric);
                });
}

void PolicyManager::advanceTriggers(int64_t day) {
  eachTriggered(
      symptom_policies_, temporal_policies_,
      [&](const std::string&, auto& policy) {
        PolicyTrigger& trigger = *policy.trigger;
        if (trigger.voluntary()) return;  // derived below, nothing to step
        if (!trigger.perUnit()) {
          trigger.step(day, observed_metrics_.value(trigger.metric));
          policy.window.trigger_state = trigger.closed;
          return;
        }
        // The units the metric counted today, plus the closed ones, which can
        // reopen on a day with no count at all. Every other unit's answer
        // cannot change.
        std::vector<VenueId> visit =
            observed_metrics_.changedUnits(trigger.metric);
        for (VenueId unit : trigger.closedUnits()) visit.push_back(unit);
        std::sort(visit.begin(), visit.end());
        visit.erase(std::unique(visit.begin(), visit.end()), visit.end());
        for (VenueId unit : visit) {
          trigger.stepUnit(unit, day,
                           observed_metrics_.unitValue(trigger.metric, unit));
        }
      });
  rebuildResponseLevels();
}

void PolicyManager::resolveResponses() {
  eachTriggered(
      symptom_policies_, temporal_policies_,
      [&](const std::string& kind, auto& policy) {
        PolicyTrigger& trigger = *policy.trigger;
        if (!trigger.voluntary()) return;
        const std::string label =
            kind + "_policy '" + policy.name + "' trigger";
        if (trigger.decision_unit.shared()) {
          const int activity =
              world_.getActivityIndex(trigger.decision_unit.shared_activity);
          if (activity < 0) {
            throw std::runtime_error(
                label + ": decision_unit names the unknown activity '" +
                trigger.decision_unit.shared_activity + "'.");
          }
          trigger.decision_unit.activity_index = static_cast<int16_t>(activity);
        }
        if (trigger.perUnit()) {
          trigger.metric_activity_index =
              observed_metrics_.activityIndexOf(trigger.terms[0].metric);
          if (trigger.metric_activity_index < 0) {
            throw std::logic_error(label +
                                   ": the metrics were not resolved before "
                                   "the policies.");
          }
        }
        trigger.world_level = trigger.base;
      });
}

void PolicyManager::rebuildResponseLevels() {
  eachTriggered(
      symptom_policies_, temporal_policies_,
      [&](const std::string&, auto& policy) {
        PolicyTrigger& trigger = *policy.trigger;
        if (!trigger.voluntary()) return;
        // The same integers and the same order of operations on every rank, so
        // the level is the same double at any process count. A unit a term
        // does not count adds nothing, exactly as a term of zero would.
        if (!trigger.perUnit()) {
          double level = trigger.base;
          for (const auto& term : trigger.terms)
            level += term.weight *
                     static_cast<double>(observed_metrics_.value(term.metric));
          trigger.world_level = level;
          return;
        }
        trigger.unit_levels.clear();
        for (const auto& term : trigger.terms) {
          for (const auto& [unit, count] :
               observed_metrics_.unitValues(term.metric)) {
            const auto it =
                trigger.unit_levels.try_emplace(unit, trigger.base).first;
            it->second += term.weight * static_cast<double>(count);
          }
        }
      });
}

bool PolicyManager::avoids(const PolicyTrigger& trigger, const Person& person,
                           VenueId slot_venue_id, uint32_t draw_index) const {
  double level = 0.0;
  VenueId unit = kInvalidVenueId;
  if (trigger.perUnit()) {
    if (slot_venue_id < 0) return false;
    unit = unitOf(world_, trigger.scope, slot_venue_id);
    level = trigger.levelAtUnit(unit);
  } else {
    level = trigger.levelAtWorld();
  }
  if (level <= 0.0) return false;

  if (trigger.perUnit()) {
    // Only people the metric counts at this unit are affected by its level.
    bool counted = false;
    for (const auto& registration :
         world_.getActivityVenues(person, trigger.metric_activity_index)) {
      if (unitOf(world_, trigger.scope, registration.first) == unit) {
        counted = true;
        break;
      }
    }
    if (!counted) return false;
  }

  int64_t key = person.id;
  if (trigger.decision_unit.shared()) {
    const auto venues =
        world_.getActivityVenues(person, trigger.decision_unit.activity_index);
    if (venues.size() != 1) {
      throw std::runtime_error(
          "policy trigger: person " + std::to_string(person.id) + " has " +
          std::to_string(venues.size()) + " venues for activity '" +
          trigger.decision_unit.shared_activity +
          "', but a shared decision unit needs exactly 1. The startup check "
          "should have caught this.");
    }
    key = venues[0].first;
  }
  SplitMix64 rng(mix_seed(applicable_policies_.baseSeed(),
                          static_cast<uint64_t>(key), draw_index,
                          kAvoidanceDrawTag));
  return level > uniform01(rng);
}

void PolicyManager::checkDecisionUnits() {
  eachTriggered(
      symptom_policies_, temporal_policies_,
      [&](const std::string& kind, const auto& policy) {
        const PolicyTrigger& trigger = *policy.trigger;
        if (!trigger.voluntary()) return;
        const std::string label = kind + "_policy '" + policy.name + "'";
        if (!trigger.decision_unit.shared()) return;

        // People the Policy applies to (exemptions are per activity and cannot
        // be settled here) whom the metric can count.
        int64_t counts[2] = {0, 0};  // checked, failed
        int64_t worst = INT64_MAX;   // lowest failing person id
        int64_t worst_venues = 0;
        for (const Person& person : world_.people) {
          if (person.is_dead || !policy.appliesTo(person, &world_)) continue;
          if (trigger.perUnit() &&
              world_.getActivityVenues(person, trigger.metric_activity_index)
                  .empty()) {
            continue;
          }
          ++counts[0];
          const auto n = world_
                             .getActivityVenues(
                                 person, trigger.decision_unit.activity_index)
                             .size();
          if (n == 1) continue;
          ++counts[1];
          if (person.id < worst) {
            worst = person.id;
            worst_venues = static_cast<int64_t>(n);
          }
        }
#ifdef USE_MPI
        const auto mpi = mpi_runtime::state();
        if (mpi.active && mpi.size > 1) {
          MPI_Allreduce(MPI_IN_PLACE, counts, 2, MPI_INT64_T, MPI_SUM,
                        MPI_COMM_WORLD);
          int64_t lowest = worst;
          MPI_Allreduce(MPI_IN_PLACE, &lowest, 1, MPI_INT64_T, MPI_MIN,
                        MPI_COMM_WORLD);
          // Only the rank that owns the lowest failing person knows their
          // venue count; the others add 0.
          int64_t venues = lowest == worst ? worst_venues : 0;
          MPI_Allreduce(MPI_IN_PLACE, &venues, 1, MPI_INT64_T, MPI_SUM,
                        MPI_COMM_WORLD);
          worst = lowest;
          worst_venues = venues;
        }
#endif
        const std::string unit_text = "decision unit shared_venue '" +
                                      trigger.decision_unit.shared_activity +
                                      "'";
        if (counts[1] > 0) {
          throw std::runtime_error(
              "policy trigger: " + label + " " + unit_text + ": person " +
              std::to_string(worst) + " has " + std::to_string(worst_venues) +
              " venues for activity '" + trigger.decision_unit.shared_activity +
              "', exactly 1 is required. " + std::to_string(counts[1]) +
              " of " + std::to_string(counts[0]) +
              " people checked failed. Fix the world data or use "
              "'decision_unit: person'.");
        }
      });
}

bool PolicyManager::unitClosedAt(const PolicyTrigger& trigger,
                                 VenueId venue) const {
  if (venue < 0) return false;
  return trigger.unitClosed(unitOf(world_, trigger.scope, venue));
}

std::vector<PolicyManager::TriggerState> PolicyManager::triggerStates() const {
  std::vector<TriggerState> out;
  eachTriggered(symptom_policies_, temporal_policies_,
                [&](const std::string& kind, const auto& policy) {
                  const PolicyTrigger& trigger = *policy.trigger;
                  TriggerState state{kind + "/" + policy.name,
                                     trigger,
                                     trigger.scope,
                                     {},
                                     trigger.voluntary()};
                  for (const auto& [unit, cell] : trigger.unitCells())
                    state.units.emplace_back(unit, cell);
                  out.push_back(std::move(state));
                });
  return out;
}

void PolicyManager::restoreTriggerStates(
    const std::vector<TriggerState>& stored) {
  const auto current = triggerStates();
  const auto refuse = [&](const std::string& why) {
    throw std::runtime_error("policy triggers: " + why);
  };
  const auto wrongSet = [&] {
    const auto join = [](const std::vector<TriggerState>& states) {
      std::string out;
      for (const auto& state : states)
        out += (out.empty() ? "" : ", ") + state.key;
      return out;
    };
    refuse("the checkpoint holds triggers [" + join(stored) +
           "] but this run has [" + join(current) +
           "]. A trigger's state cannot be rebuilt; resume with the same "
           "triggered Policies.");
  };
  // Matched by key, so triggered policies may be reordered between runs.
  std::unordered_map<std::string, const TriggerState*> stored_by_key;
  for (const auto& state : stored)
    if (!stored_by_key.emplace(state.key, &state).second) wrongSet();
  if (stored_by_key.size() != current.size()) wrongSet();
  for (const auto& now : current) {
    const auto found = stored_by_key.find(now.key);
    if (found == stored_by_key.end()) wrongSet();
    const TriggerState& was = *found->second;
    // The same name counted per unit in one run and as a world total in the
    // other is not the same trigger, whatever the key says. Nor is a Policy
    // read as a threshold in one run and as a response in the other.
    if (was.scope != now.scope) {
      refuse("'" + was.key +
             "' reads a metric of another scope than the checkpoint's. A "
             "trigger's state cannot be rebuilt under another scope; resume "
             "with the same metrics.");
    }
    if (was.voluntary != now.voluntary) {
      refuse("'" + was.key + "' is " +
             (was.voluntary ? "a voluntary response" : "a threshold") +
             " trigger in the checkpoint and " +
             (now.voluntary ? "a voluntary response" : "a threshold") +
             " trigger in this run. Resume with the same triggered Policies.");
    }
  }
  eachTriggered(symptom_policies_, temporal_policies_,
                [&](const std::string& kind, auto& policy) {
                  const TriggerState& state =
                      *stored_by_key.at(kind + "/" + policy.name);
                  static_cast<TriggerCell&>(*policy.trigger) = state.cell;
                  policy.trigger->restoreUnits(state.units);
                  policy.window.trigger_state = policy.trigger->holdsPolicy();
                });
}

void PolicyManager::resolveTransmissionEffects(const Disease& disease) {
  transmission_mode_count_ =
      static_cast<size_t>(std::max(1, disease.numModes()));
  for (auto& effect : transmission_effects_) {
    const bool person_scope = effect.scope == TransmissionEffectScope::Person;
    const bool person_channel =
        effect.channel == TransmissionEffectChannel::SourceInfectiousness ||
        effect.channel == TransmissionEffectChannel::TargetSusceptibility;
    if (person_scope != person_channel)
      throw std::runtime_error(
          "transmission policy effect scope/channel mismatch");
    if (effect.policy_kind == TransmissionPolicyKind::Symptom) {
      if (effect.policy_index >= symptom_policies_.size())
        throw std::runtime_error(
            "transmission effect names an unknown symptom policy");
      if (!person_scope)
        throw std::runtime_error("venue effects require a temporal policy");
    } else {
      if (effect.policy_index >= temporal_policies_.size())
        throw std::runtime_error(
            "transmission effect names an unknown temporal policy");
      if (!person_scope &&
          temporal_policies_[effect.policy_index].action.compliance_rate != 1.0)
        throw std::runtime_error("venue effects require compliance_rate: 1");
      if (!person_scope &&
          !temporal_policies_[effect.policy_index].applies_to.empty())
        throw std::runtime_error(
            "venue effects cannot use person-level applies_to criteria");
    }
    int mode = -1;
    for (int i = 0; i < disease.numModes(); ++i) {
      if (disease.getModeName(static_cast<uint8_t>(i)) == effect.mode_name) {
        mode = i;
        break;
      }
    }
    if (mode < 0)
      throw std::runtime_error("unknown transmission mode '" +
                               effect.mode_name + "'");
    const auto mode_type = disease.getTransmissionParams().modes[mode].type;
    if (mode_type == TransmissionModeType::CompartmentalDeposition &&
        effect.channel != TransmissionEffectChannel::SourceInfectiousness)
      throw std::runtime_error(
          "compartmental deposition only supports source_infectiousness "
          "effects");
    if (mode_type == TransmissionModeType::CompartmentalUptake &&
        effect.channel == TransmissionEffectChannel::SourceInfectiousness)
      throw std::runtime_error(
          "compartmental uptake has no person source infectiousness");
    effect.mode_index = static_cast<uint8_t>(mode);
    for (auto& criterion : effect.criteria) {
      const std::string context =
          "transmission effect for mode '" + effect.mode_name + "'";
      if (person_scope)
        criterion.resolveOrThrow(world_, context);
      else
        criterion.resolveVenueOrThrow(world_, context);
    }
  }
}

PolicyManager::PolicyWindowBits PolicyManager::policyWindowBits(
    double current_time) const {
  PolicyWindowBits bits;
  for (size_t i = 0; i < temporal_policies_.size(); ++i)
    if (temporal_policies_[i].window.contains(current_time))
      bits.temporal.set(i);
  for (size_t i = 0; i < symptom_policies_.size(); ++i)
    if (symptom_policies_[i].window.contains(current_time)) bits.symptom.set(i);
  return bits;
}

void PolicyManager::rebuildPersonTransmissionModifier(Person& person,
                                                      double current_time) {
  auto values = transmission_modifier_table_.get(0);
  // complies(), not participates(): an effect on a voluntary policy still
  // takes only the share of its applicable Persons its compliance rate draws.
  for (const auto& effect : transmission_effects_) {
    if (effect.scope != TransmissionEffectScope::Person) continue;
    if (effect.policy_kind == TransmissionPolicyKind::Temporal) {
      const auto applicable =
          applicable_policies_.walk(person, temporal_policies_)
              .find(effect.policy_index);
      if (!applicable || !applicable->policy().window.contains(current_time) ||
          !applicable->complies())
        continue;
    } else {
      const auto applicable =
          applicable_policies_.walk(person, symptom_policies_)
              .find(effect.policy_index);
      if (!applicable || !applicable->policy().window.contains(current_time) ||
          !person.infection ||
          !applicable->policy().triggeredBy(
              person.infection->getTrajectory().getCurrentSymptomId(
                  current_time)) ||
          !applicable->complies())
        continue;
    }
    if (!filtering::matchesCriteria(person, &world_, effect.criteria)) continue;
    double& value =
        values[effect.mode_index][static_cast<size_t>(effect.channel)];
    value *= effect.multiplier;
    if (!std::isfinite(value))
      throw std::runtime_error("transmission policy multipliers overflow");
  }
  person.transmission_modifier_set_id =
      transmission_modifier_table_.intern(values);
}

void PolicyManager::rebuildVenueTransmissionModifier(Venue& venue,
                                                     double current_time) {
  auto values = transmission_modifier_table_.get(0);
  for (const auto& effect : transmission_effects_) {
    if (effect.scope != TransmissionEffectScope::Venue) continue;
    if (!temporal_policies_[effect.policy_index].window.contains(current_time))
      continue;
    bool matches = true;
    for (const auto& criterion : effect.criteria) {
      if (!criterion.evaluate(venue, &world_)) {
        matches = false;
        break;
      }
    }
    if (matches) {
      double& value =
          values[effect.mode_index][static_cast<size_t>(effect.channel)];
      value *= effect.multiplier;
      if (!std::isfinite(value))
        throw std::runtime_error("transmission policy multipliers overflow");
    }
  }
  venue.transmission_modifier_set_id =
      transmission_modifier_table_.intern(values);
}

void PolicyManager::rebuildTransmissionModifiers(double current_time) {
  transmission_modifier_table_.reset(transmission_mode_count_);
  for (auto& person : world_.people)
    rebuildPersonTransmissionModifier(person, current_time);
  for (auto& venue : world_.venues)
    rebuildVenueTransmissionModifier(venue, current_time);
  last_policy_window_bits_ = policyWindowBits(current_time);
  last_symptom_ids_.clear();
  for (const auto& person : world_.people) {
    if (person.infection)
      last_symptom_ids_[person.id] =
          person.infection->getTrajectory().getCurrentSymptomId(current_time);
  }
}

void PolicyManager::initializeTransmissionModifiers(const Disease& disease,
                                                    double current_time) {
  if (transmission_effects_.empty()) return;
  if (transmission_mode_count_ == 0) resolveTransmissionEffects(disease);
  rebuildTransmissionModifiers(current_time);
  transmission_modifiers_initialized_ = true;
}

void PolicyManager::refreshTransmissionModifiers(
    double current_time,
    const std::unordered_set<PersonId>& active_infections) {
  if (!transmission_modifiers_initialized_) return;
  if (policyWindowBits(current_time) != last_policy_window_bits_) {
    rebuildTransmissionModifiers(current_time);
    return;
  }
  bool has_symptom_effect = false;
  for (const auto& effect : transmission_effects_)
    if (effect.policy_kind == TransmissionPolicyKind::Symptom) {
      has_symptom_effect = true;
      break;
    }
  if (!has_symptom_effect) return;

  for (auto it = last_symptom_ids_.begin(); it != last_symptom_ids_.end();) {
    if (active_infections.count(it->first)) {
      ++it;
      continue;
    }
    if (Person* person = world_.getPerson(it->first))
      rebuildPersonTransmissionModifier(*person, current_time);
    it = last_symptom_ids_.erase(it);
  }
  for (PersonId id : active_infections) {
    Person* person = world_.getPerson(id);
    if (!person || !person->infection) continue;
    uint16_t symptom =
        person->infection->getTrajectory().getCurrentSymptomId(current_time);
    auto [it, inserted] = last_symptom_ids_.emplace(id, symptom);
    if (inserted || it->second != symptom) {
      rebuildPersonTransmissionModifier(*person, current_time);
      it->second = symptom;
    }
  }
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
    Person& person, int16_t scheduled_activity_index, VenueId pin_venue_id,
    SubsetIndex pin_subset_index, SlotVenueType slot_venue_type,
    double current_time, int time_slot_index, const Person* partner) {
  // Resolved lazily, and cached across both policy loops: this runs once per
  // person per slot across the whole population, so an ungated run must never
  // pay for the lookup. A Deferred value naming a venue this rank cannot type
  // throws here rather than degrading to kUnknownVenueTypeId.
  LazySlotVenueTypeResolver slotVenueType(slot_venue_type, world_);

  // Priority 1: symptom-based policies
  if (person.infection != nullptr) {
    uint16_t current_symptom_id =
        person.infection->getTrajectory().getCurrentSymptomId(current_time);

    for (const auto& applicable :
         applicable_policies_.walk(person, symptom_policies_)) {
      const size_t i = applicable.index();
      const auto& policy = applicable.policy();

      // Out of window: the policy is not in force, so it holds nobody. A
      // calendar edge is a change of government instruction, not a change in
      // the person, so no follow-up policy inherits the decision.
      if (!policy.window.contains(current_time)) {
        releaseFreeze(person, i);
        applicable.clear();
        continue;
      }

      bool is_triggered = policy.triggeredBy(current_symptom_id);

      // STICKY COMPLIANCE: if not triggered, clear participation/decision and
      // continue (with follow-up policy inheritance check)
      if (!is_triggered) {
        if (policy.follow_up_policy_index >= 0)
          applicable.inheritTo(
              static_cast<size_t>(policy.follow_up_policy_index),
              policy.inherit_compliance, policy.inherit_refusal);

        // Restore paused hop state if this policy was responsible for freezing
        releaseFreeze(person, i);
        applicable.clear();
        continue;
      }

      if (!applicable.latch()) continue;

      if (!actionApplies(policy.action, person, scheduled_activity_index,
                         slotVenueType, partner, slotTrigger(policy),
                         slot_venue_type.venueId(),
                         applicable.complianceKey())) {
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
        // to the last overnight venue, so pin_venue_id is the best candidate.
        VenueId pin_venue = pin_venue_id;
        SubsetIndex pin_subset = pin_subset_index;
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
    applicable_policies_.clear(person, PolicyKind::Symptom);
  }

  // Priority 2: temporal policies (lockdowns, etc.)
  for (const auto& applicable :
       applicable_policies_.walk(person, temporal_policies_)) {
    const auto& policy = applicable.policy();

    bool is_active = policy.window.contains(current_time);

    if (!is_active) {
      applicable.clear();
      continue;
    }

    if (!applicable.latch()) continue;

    // Where the person is scheduled, not where a freeze would pin them: a
    // traveller in transit is pinned at their last overnight venue while
    // occupying none.
    if (!actionApplies(policy.action, person, scheduled_activity_index,
                       slotVenueType, partner, slotTrigger(policy),
                       slot_venue_type.venueId(), applicable.complianceKey())) {
      continue;
    }

    return getReplacementLocation(person, policy.action.replacement_activity,
                                  policy.action.replacement_activity_index);
  }

  return std::nullopt;
}

bool PolicyManager::suppressesParticipation(
    const Person& person, int16_t activity_index,
    const SlotVenueType& slot_venue_type, double current_time,
    const Person* partner) const {
  // A frozen Person is pinned at the venue the freeze anchored them to, so they
  // are not at the encounter's venue or the host's venue whatever those are.
  // Answered ahead of the loops: there is no meaningful venue-gate question to
  // ask about a venue they cannot be in.
  if (frozen_states_.count(person.id) > 0) return true;

  // Same lazy resolution as getOverride, and for the same reason: this runs
  // once per person per slot across the population, so an ungated run must
  // never pay for the lookup.
  LazySlotVenueTypeResolver slotVenueType(slot_venue_type, world_);

  // Priority 1: symptom-based policies. Untriggered policies are skipped
  // outright — getOverride uses that branch to release a freeze and propagate
  // follow-up inheritance, and neither is this function's business.
  if (person.infection != nullptr) {
    uint16_t current_symptom_id =
        person.infection->getTrajectory().getCurrentSymptomId(current_time);

    for (const auto& applicable :
         applicable_policies_.walk(person, symptom_policies_)) {
      const auto& policy = applicable.policy();
      if (!policy.window.contains(current_time)) continue;
      if (!policy.triggeredBy(current_symptom_id)) continue;

      if (!applicable.participates()) continue;

      if (actionApplies(policy.action, person, activity_index, slotVenueType,
                        partner, slotTrigger(policy), slot_venue_type.venueId(),
                        applicable.complianceKey())) {
        return true;
      }
    }
  }

  // Priority 2: temporal policies (lockdowns, etc.)
  for (const auto& applicable :
       applicable_policies_.walk(person, temporal_policies_)) {
    const auto& policy = applicable.policy();
    if (!policy.window.contains(current_time)) continue;

    if (!applicable.participates()) continue;

    if (actionApplies(policy.action, person, activity_index, slotVenueType,
                      partner, slotTrigger(policy), slot_venue_type.venueId(),
                      applicable.complianceKey())) {
      return true;
    }
  }

  return false;
}

void PolicyManager::precomputePolicyApplicability(std::vector<Person>& people) {
  for (auto& person : people) {
    if (person.is_dead) continue;

    PolicyMask symptom_mask;
    for (size_t i = 0;
         i < std::min(symptom_policies_.size(), kMaxPoliciesPerKind); ++i) {
      const auto& policy = symptom_policies_[i];
      if (policy.appliesTo(person, &world_)) {
        symptom_mask.set(i);
      }
    }

    PolicyMask temporal_mask;
    for (size_t i = 0;
         i < std::min(temporal_policies_.size(), kMaxPoliciesPerKind); ++i) {
      const auto& policy = temporal_policies_[i];
      if (policy.appliesTo(person, &world_)) {
        temporal_mask.set(i);
      }
    }

    ApplicablePolicyMasks masks;
    masks.symptom.applicable = symptom_mask;
    masks.temporal.applicable = temporal_mask;
    applicable_policies_.importState(person, masks);
  }
}

}  // namespace june
