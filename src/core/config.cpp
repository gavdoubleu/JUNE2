#include "core/config.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <unordered_set>

#include "core/world_state.h"
#include "utils/config_checks.h"
#include "utils/filtered_csv.h"
#include "utils/filtering.h"

namespace june {

DistributionType parseDistributionType(const std::string& s) {
  if (s == "poisson") return DistributionType::POISSON;
  if (s == "binomial") return DistributionType::BINOMIAL;
  if (s == "fixed") return DistributionType::FIXED;
  throw std::runtime_error("Unknown distribution type: '" + s +
                           "'. Must be 'poisson', 'binomial', or 'fixed'.");
}

const char* distributionTypeToString(DistributionType t) {
  switch (t) {
    case DistributionType::POISSON:
      return "poisson";
    case DistributionType::BINOMIAL:
      return "binomial";
    case DistributionType::FIXED:
      return "fixed";
    default:
      return "unknown";
  }
}

void SelectionCriterion::resolve(const WorldState& world) {
  // 1. First ensure type is cached
  if (cached_type == PropertyType::UNKNOWN) {
    if (property_path == "age")
      cached_type = PropertyType::AGE;
    else if (property_path == "sex")
      cached_type = PropertyType::SEX;
    else if (property_path == "geo_unit_id")
      cached_type = PropertyType::GEO_ID;
    else if (property_path == "id")
      cached_type = PropertyType::PERSON_ID;
    else if (property_path.compare(0, 11, "activities.") == 0) {
      size_t dot1 = property_path.find('.');
      size_t dot2 = property_path.find('.', dot1 + 1);
      if (dot2 != std::string::npos) {
        cached_activity_name = property_path.substr(dot1 + 1, dot2 - dot1 - 1);
        cached_sub_property = property_path.substr(dot2 + 1);
        if (cached_sub_property == "length")
          cached_type = PropertyType::ACTIVITY_LENGTH;
        else if (cached_sub_property == "venue_type")
          cached_type = PropertyType::ACTIVITY_VENUE_TYPE;
      }
    } else if (property_path.compare(0, 11, "properties.") == 0) {
      cached_type = PropertyType::CUSTOM_PROPERTY;
      cached_sub_property = property_path.substr(11);
      cached_prop_idx = world.getPersonPropertyIndex(cached_sub_property);
    } else if (property_path.compare(0, 9, "networks.") == 0) {
      size_t dot1 = property_path.find('.');
      size_t dot2 = property_path.find('.', dot1 + 1);
      if (dot2 != std::string::npos) {
        cached_activity_name = property_path.substr(
            dot1 + 1, dot2 - dot1 - 1);  // Use for network name
        cached_sub_property = property_path.substr(dot2 + 1);
        if (cached_sub_property == "length")
          cached_type = PropertyType::NETWORK_SIZE;
      }
    } else if (property_path == "is_alive") {
      cached_type = PropertyType::IS_ALIVE;
    } else if (property_path.compare(0, 19, "partner_in_network(") == 0) {
      size_t open_paren = property_path.find('(');
      size_t close_paren = property_path.find(')', open_paren);
      if (open_paren != std::string::npos && close_paren != std::string::npos) {
        cached_activity_name =
            property_path.substr(open_paren + 1, close_paren - open_paren - 1);
        cached_type = PropertyType::PARTNER_IN_NETWORK;
      }
    }
  }

  // 2. Resolve target_code for equality comparisons
  if (operator_type == "==" || operator_type == "!=") {
    if (std::holds_alternative<std::string>(value)) {
      const std::string& target_val = std::get<std::string>(value);

      if (cached_type == PropertyType::SEX) {
        if (target_val == "male" || target_val == "M" || target_val == "Male")
          target_code = 0;
        else if (target_val == "female" || target_val == "F" ||
                 target_val == "Female")
          target_code = 1;
        else
          target_code = 2;  // unknown
      } else if (cached_type == PropertyType::CUSTOM_PROPERTY &&
                 cached_prop_idx >= 0) {
        const auto& prop_name = world.person_property_names[cached_prop_idx];
        auto it_reg = world.person_property_value_registries.find(prop_name);
        if (it_reg != world.person_property_value_registries.end()) {
          const auto& registry = it_reg->second;
          auto it = std::find(registry.begin(), registry.end(), target_val);
          if (it != registry.end()) {
            target_code =
                static_cast<int32_t>(std::distance(registry.begin(), it));
          }
        }
      } else if (cached_type == PropertyType::ACTIVITY_VENUE_TYPE) {
        auto it = std::find(world.venue_type_names.begin(),
                            world.venue_type_names.end(), target_val);
        if (it != world.venue_type_names.end()) {
          target_code = static_cast<int32_t>(
              std::distance(world.venue_type_names.begin(), it));
        }
      }
    }
  }
}

bool SelectionCriterion::evaluate(const Person& person, const WorldState* world,
                                  const Person* partner) const {
  // 1. Resolve property type and path if not cached
  if (cached_type == PropertyType::UNKNOWN) {
    const_cast<SelectionCriterion*>(this)->resolve(*world);
    if (cached_type == PropertyType::UNKNOWN) return false;
  }

  // Boolean predicates: bypass the target-code / fallback machinery entirely.
  // The criterion's `value` is bool; the operator must be == or !=.
  auto eval_bool = [this](bool actual) -> bool {
    if (!std::holds_alternative<bool>(value)) return false;
    bool target = std::get<bool>(value);
    if (operator_type == "==") return actual == target;
    if (operator_type == "!=") return actual != target;
    return false;
  };
  if (cached_type == PropertyType::IS_ALIVE) {
    return eval_bool(!person.is_dead);
  }

  // 2. Integer comparison for interned properties
  if (target_code != -1) {
    int32_t p_val_code = -1;
    if (cached_type == PropertyType::SEX) {
      p_val_code = static_cast<int32_t>(person.sex);
    } else if (cached_type == PropertyType::CUSTOM_PROPERTY) {
      if (world) {
        auto props = world->getPersonProperties(person);
        if (cached_prop_idx >= 0 && cached_prop_idx < (int)props.size()) {
          p_val_code = props[cached_prop_idx];
        }
      }
    } else if (cached_type == PropertyType::ACTIVITY_VENUE_TYPE) {
      if (!world) return false;
      auto activity_venues =
          world->getActivityVenues(person, cached_activity_name);
      for (const auto& av : activity_venues) {
        if (av.first >= 0 && av.first < (int)world->venues.size()) {
          const auto& v = world->venues[av.first];
          int32_t v_type_code = static_cast<int32_t>(v.type_id);
          if (operator_type == "==") {
            if (v_type_code == target_code) return true;
          } else if (operator_type == "!=") {
            if (v_type_code != target_code) return true;
          }
        }
      }
      return (operator_type ==
              "!=");  // True if none matched and looking for !=
    }

    if (p_val_code != -1) {
      if (operator_type == "==") return p_val_code == target_code;
      if (operator_type == "!=") return p_val_code != target_code;
    }
  }

  // 3. Fallback: Fetch person value(s) based on type
  std::vector<PropertyValue> person_values;

  switch (cached_type) {
    case PropertyType::AGE:
      // Age is float on Person; preserve the fractional part. The compare
      // lambda below handles int-threshold ⇆ double-value comparisons by
      // promoting both sides to double, so YAML expressions like
      // `filter.age<=59` correctly reject an age of 59.5 (matching the
      // legacy hardcoded check against rd.max_active_age).
      person_values.push_back(PropertyValue(static_cast<double>(person.age)));
      break;
    case PropertyType::SEX: {
      std::string sex_str = (person.sex == Sex::MALE)     ? "male"
                            : (person.sex == Sex::FEMALE) ? "female"
                                                          : "unknown";
      person_values.push_back(PropertyValue(sex_str));
      break;
    }
    case PropertyType::GEO_ID:
      person_values.push_back(
          PropertyValue(static_cast<int>(person.geo_unit_id)));
      break;
    case PropertyType::PERSON_ID:
      person_values.push_back(PropertyValue(static_cast<int32_t>(person.id)));
      break;
    case PropertyType::ACTIVITY_LENGTH:
      if (!world) return false;
      person_values.push_back(PropertyValue(static_cast<int>(
          world->getActivityVenues(person, cached_activity_name).size())));
      break;
    case PropertyType::ACTIVITY_VENUE_TYPE: {
      if (!world) return false;
      auto activity_venues =
          world->getActivityVenues(person, cached_activity_name);
      for (const auto& av : activity_venues) {
        if (av.first >= 0 && av.first < (int)world->venues.size()) {
          const auto& v = world->venues[av.first];
          if (v.type_id < world->venue_type_names.size()) {
            person_values.push_back(
                PropertyValue(world->venue_type_names[v.type_id]));
          }
        }
      }
      break;
    }
    case PropertyType::CUSTOM_PROPERTY: {
      if (!world) return false;
      auto prop_opt = world->getPersonProperty(
          person, world->person_property_names[cached_prop_idx]);
      if (prop_opt.has_value()) {
        person_values.push_back(*prop_opt);
      }
      break;
    }
    case PropertyType::NETWORK_SIZE: {
      if (!world) return false;
      size_t n = world->getNetworkPartners(person, cached_activity_name).size();
      person_values.push_back(PropertyValue(static_cast<int>(n)));
      break;
    }
    case PropertyType::PARTNER_IN_NETWORK: {
      if (!world || !partner) return false;
      auto partners = world->getNetworkPartners(person, cached_activity_name);
      bool found = std::find(partners.begin(), partners.end(), partner->id) !=
                   partners.end();
      person_values.push_back(PropertyValue(found ? 1 : 0));
      break;
    }
    default:
      return false;
  }

  if (person_values.empty()) return false;

  // 4. Perform comparison for each fetched value (multi-venue support)
  auto compare = [this](const PropertyValue& p_val) -> bool {
    if (operator_type == "==") return p_val == value;
    if (operator_type == "!=") return p_val != value;

    // Numeric helper: extract a double from PropertyValue if it holds an
    // int or double. Returns false in `out_set` if the value is non-numeric
    // (e.g. a string), in which case the caller falls through to the
    // type-mismatch return below. This lets `filter.age<=59` work even
    // when age is stored as double and the threshold parses as int.
    auto as_number = [](const PropertyValue& pv, double& out) -> bool {
      if (std::holds_alternative<int>(pv)) {
        out = static_cast<double>(std::get<int>(pv));
        return true;
      }
      if (std::holds_alternative<double>(pv)) {
        out = std::get<double>(pv);
        return true;
      }
      return false;
    };

    double lhs = 0.0, rhs = 0.0;
    if (operator_type == ">" || operator_type == "<" || operator_type == ">=" ||
        operator_type == "<=") {
      if (!as_number(p_val, lhs) || !as_number(value, rhs)) return false;
      if (operator_type == ">") return lhs > rhs;
      if (operator_type == "<") return lhs < rhs;
      if (operator_type == ">=") return lhs >= rhs;
      if (operator_type == "<=") return lhs <= rhs;
    } else if (operator_type == "in") {
      if (std::holds_alternative<std::vector<int32_t>>(value) &&
          std::holds_alternative<int>(p_val)) {
        const auto& list = std::get<std::vector<int32_t>>(value);
        int val = std::get<int>(p_val);
        return std::find(list.begin(), list.end(), val) != list.end();
      }
    } else if (operator_type == "contains") {
      if (std::holds_alternative<std::string>(p_val) &&
          std::holds_alternative<std::string>(value)) {
        return std::get<std::string>(p_val).find(
                   std::get<std::string>(value)) != std::string::npos;
      }
    }
    return false;
  };

  for (const auto& pv : person_values) {
    if (compare(pv)) return true;
  }
  return false;
}

void SimulationConfig::resolve(const WorldState& world) {
  // Resolve the partial-presence venue type names into a bitmask of
  // venue_type_ids + a per-id target_group_size lookup. Unknown names are
  // silently ignored (the world may not contain every declared type, e.g.
  // tube_line in a Durham-only world has no successful tube routes and is
  // absent from the registry). Throws if a declared type id exceeds the
  // bitmask width or if target_group_size is non-positive.
  partial_presence.enabled_venue_type_mask = 0;
  partial_presence.target_group_size_by_type_id.assign(
      world.venue_type_names.size(), 0);
  for (const auto& [name, tgs] : partial_presence.target_group_size_by_name) {
    int idx = world.getVenueTypeIndex(name);
    if (idx < 0) continue;
    if (idx >= 64) {
      throw std::runtime_error(
          "SimulationConfig::resolve: partial_presence venue type id " +
          std::to_string(idx) + " ('" + name +
          "') exceeds 64-bit mask width; promote enabled_venue_type_mask to "
          "a wider bitset.");
    }
    if (tgs <= 0) {
      throw std::runtime_error(
          "SimulationConfig::resolve: partial_presence target_group_size for "
          "venue type '" +
          name + "' must be > 0 (got " + std::to_string(tgs) + ").");
    }
    partial_presence.enabled_venue_type_mask |= (uint64_t(1) << idx);
    partial_presence.target_group_size_by_type_id[idx] = tgs;
  }
}

void ContactMatrixConfig::resolve(const WorldState& world) {
  const auto& venue_names = world.venue_type_names;
  betas_by_id.assign(venue_names.size(), default_beta);
  matrices_by_id.assign(venue_names.size(), nullptr);

  for (const auto& [name, beta] : betas) {
    int idx = world.getVenueTypeIndex(name);
    if (idx >= 0 && idx < (int)betas_by_id.size()) {
      betas_by_id[idx] = beta;
    }
  }

  // Helper lambda to resolve age_to_bin for a ContactMatrix
  auto resolveAgeToBin = [](ContactMatrix& matrix) {
    std::fill(std::begin(matrix.age_to_bin), std::end(matrix.age_to_bin), -1);
    matrix.has_age_bins = false;
    for (size_t b = 0; b < matrix.bins.size(); ++b) {
      const std::string& bin_name = matrix.bins[b];
      int min_age = -1, max_age = -1;
      if (!bin_name.empty()) {
        if (bin_name.back() == '+') {
          try {
            min_age = std::stoi(bin_name.substr(0, bin_name.size() - 1));
            max_age = 99;
          } catch (...) {
          }
        } else {
          size_t sep_pos = bin_name.find_first_of("-,");
          if (sep_pos != std::string::npos) {
            try {
              size_t start_pos = (bin_name[0] == '[') ? 1 : 0;
              min_age =
                  std::stoi(bin_name.substr(start_pos, sep_pos - start_pos));
              size_t end_pos = bin_name.find(']', sep_pos + 1);
              if (end_pos == std::string::npos) end_pos = bin_name.size();
              max_age = std::stoi(
                  bin_name.substr(sep_pos + 1, end_pos - sep_pos - 1));
            } catch (...) {
              min_age = -1;
              max_age = -1;
            }
          }
        }
      }
      if (min_age >= 0 && max_age >= min_age) {
        matrix.has_age_bins = true;
        for (int a = std::max(0, min_age); a <= std::min(99, max_age); ++a) {
          if (matrix.age_to_bin[a] < 0) {
            matrix.age_to_bin[a] = static_cast<int>(b);
          }
        }
      }
    }
  };

  for (auto& [name, matrix] : matrices) {
    int idx = world.getVenueTypeIndex(name);
    if (idx >= 0 && idx < (int)matrices_by_id.size()) {
      matrices_by_id[idx] = &matrix;
    }

    // Pre-resolve bin lookups for this matrix
    matrix.male_bin = matrix.findBinIndex("male");
    matrix.female_bin = matrix.findBinIndex("female");

    // Pre-resolve subset_type -> bin mapping
    matrix.bin_by_subset_type.assign(world.subset_type_names.size(), -1);
    for (size_t st = 0; st < world.subset_type_names.size(); ++st) {
      matrix.bin_by_subset_type[st] =
          matrix.findBinIndex(world.subset_type_names[st]);
    }

    // Pre-resolve age -> bin index lookup table
    resolveAgeToBin(matrix);
  }

  // Populate mode_matrices_by_id for per-mode contact matrix lookups.
  // mode_matrices[venue_type][mode_name] → ContactMatrix
  if (!mode_matrices.empty() && !mode_names.empty()) {
    int n_modes = static_cast<int>(mode_names.size());
    mode_matrices_by_id.assign(
        venue_names.size(),
        std::vector<const ContactMatrix*>(n_modes, nullptr));

    for (auto& [venue_name, mode_map] : mode_matrices) {
      int venue_idx = world.getVenueTypeIndex(venue_name);
      if (venue_idx < 0 || venue_idx >= (int)venue_names.size()) continue;

      for (int m = 0; m < n_modes; ++m) {
        auto it = mode_map.find(mode_names[m]);
        if (it != mode_map.end()) {
          mode_matrices_by_id[venue_idx][m] = &it->second;
          // Pre-resolve bin lookups for mode matrices too
          auto& cm = it->second;
          cm.male_bin = cm.findBinIndex("male");
          cm.female_bin = cm.findBinIndex("female");
          cm.bin_by_subset_type.assign(world.subset_type_names.size(), -1);
          for (size_t st = 0; st < world.subset_type_names.size(); ++st) {
            cm.bin_by_subset_type[st] =
                cm.findBinIndex(world.subset_type_names[st]);
          }
          resolveAgeToBin(cm);
        }
      }
    }
  }

  // Populate encounter-id-indexed virtual-encounter matrix arrays.
  // Virtual encounter matrices (e.g. "group_sex", "romantic_encounter")
  // are stored under string keys in `matrices` / `mode_matrices` but don't
  // correspond to any venue type, so they never reach `matrices_by_id` /
  // `mode_matrices_by_id` above. Build a parallel integer-indexed lookup
  // keyed by encounter_type_id so the transmission hot path in
  // InteractionManager::processVenueTransmissions stays branch-and-index
  // without string work. Uses `virtual_matrix_names`, which was built by
  // CoordinatedEncounterConfig::resolve from each is_virtual encounter's
  // virtual_contact_matrix field.
  const size_t n_enc_types = world.encounter_type_names.size();
  virtual_matrices_by_encounter_id.assign(n_enc_types, nullptr);
  const int n_modes_virtual = static_cast<int>(mode_names.size());
  if (n_modes_virtual > 0) {
    virtual_mode_matrices_by_encounter_id.assign(
        n_enc_types,
        std::vector<const ContactMatrix*>(n_modes_virtual, nullptr));
  } else {
    virtual_mode_matrices_by_encounter_id.clear();
  }

  // Track which matrices we've already bin-resolved so we don't redo the
  // work when several virtual encounter types alias onto the same matrix
  // (e.g. ooe_encounter, romantic_encounters, cohabiting_encounters all
  // point at "romantic_encounter").
  std::unordered_set<const ContactMatrix*> bin_resolved;

  auto resolve_matrix_bins = [&](ContactMatrix& cm) {
    if (!bin_resolved.insert(&cm).second) return;
    cm.male_bin = cm.findBinIndex("male");
    cm.female_bin = cm.findBinIndex("female");
    cm.bin_by_subset_type.assign(world.subset_type_names.size(), -1);
    for (size_t st = 0; st < world.subset_type_names.size(); ++st) {
      cm.bin_by_subset_type[st] = cm.findBinIndex(world.subset_type_names[st]);
    }
    resolveAgeToBin(cm);
  };

  for (const auto& [eid, matrix_name] : virtual_matrix_names) {
    if (static_cast<size_t>(eid) >= n_enc_types) continue;

    auto flat_it = matrices.find(matrix_name);
    if (flat_it != matrices.end()) {
      virtual_matrices_by_encounter_id[eid] = &flat_it->second;
      resolve_matrix_bins(flat_it->second);
    }

    if (n_modes_virtual <= 0) continue;
    auto mode_it = mode_matrices.find(matrix_name);
    if (mode_it == mode_matrices.end()) continue;
    auto& per_mode = mode_it->second;
    for (int m = 0; m < n_modes_virtual; ++m) {
      auto it = per_mode.find(mode_names[m]);
      if (it == per_mode.end()) continue;
      virtual_mode_matrices_by_encounter_id[eid][m] = &it->second;
      resolve_matrix_bins(it->second);
    }
  }
}

void PreferenceProfile::resolve(const WorldState& world) {
  for (auto& crit : selection_criteria) {
    crit.resolve(world);
  }

  activity_id = world.getActivityIndex(activity);

  weights_by_id.assign(world.venue_type_names.size(), 1.0);
  for (const auto& [name, weight] : preference_weights) {
    int idx = world.getVenueTypeIndex(name);
    if (idx >= 0 && idx < (int)weights_by_id.size()) {
      weights_by_id[idx] = weight;
    }
  }
}

// Helper: Compute a bitmask of activity indices from a list of activity names
static ActivityMask computeActivityMaskFromNames(
    const std::vector<std::string>& activities,
    const std::vector<std::string>& activity_names) {
  ActivityMask mask = 0;
  for (const auto& act : activities) {
    for (size_t i = 0; i < activity_names.size(); ++i) {
      if (activity_names[i] == act) {
        mask |= (ActivityMask(1) << i);
        break;
      }
    }
  }
  return mask;
}

void ScheduleConfig::resolveSlots(const WorldState& world) {
  // Build cycle_to_type_idx: cycle position -> index into day_type_names
  cycle_to_type_idx.resize(day_type_cycle.size());
  for (size_t i = 0; i < day_type_cycle.size(); ++i) {
    auto it = std::find(day_type_names.begin(), day_type_names.end(),
                        day_type_cycle[i]);
    cycle_to_type_idx[i] =
        (it != day_type_names.end())
            ? static_cast<int>(std::distance(day_type_names.begin(), it))
            : 0;
  }

  int num_dt = static_cast<int>(day_type_names.size());
  size_t num_acts = world.activity_names.size();

  // Helper: resolve slot vector caches
  auto resolveSlotVec = [&](std::vector<TimeSlot>& slots) {
    for (auto& slot : slots) {
      slot.allowed_activity_mask = computeActivityMaskFromNames(
          slot.allowed_activities, world.activity_names);
      slot.coordinated_only_activity_mask = computeActivityMaskFromNames(
          slot.coordinated_only_activities, world.activity_names);
      slot.allowed_activity_indices.clear();
      for (const auto& act : slot.allowed_activities) {
        int idx = world.getActivityIndex(act);
        if (idx >= 0)
          slot.allowed_activity_indices.push_back(static_cast<int16_t>(idx));
      }
      if (slot.specified_activity.has_value()) {
        auto& spec = slot.specified_activity.value();
        spec.cached_activity_idx =
            static_cast<int16_t>(world.getActivityIndex(spec.type));
        if (spec.venue_type.has_value()) {
          spec.cached_venue_type_idx =
              world.getVenueTypeIndex(spec.venue_type.value());
        }
      }
      // Resolve hop_on_activity -> hop_schedule_by_activity_idx
      slot.hop_schedule_by_activity_idx.assign(world.activity_names.size(), -1);
      for (const auto& [act_name, sched_name] : slot.hop_on_activity) {
        int act_idx = world.getActivityIndex(act_name);
        auto s_it = std::find_if(
            schedule_types.begin(), schedule_types.end(),
            [&](const ScheduleType& s) { return s.name == sched_name; });
        if (act_idx >= 0 && s_it != schedule_types.end()) {
          slot.hop_schedule_by_activity_idx[act_idx] =
              static_cast<int16_t>(s_it - schedule_types.begin());
        }
      }
      // Resolve property_hop_dispatch -> property_hop_dispatch_by_activity_idx
      for (const auto& [act_name, dispatch] : slot.property_hop_dispatch) {
        int act_idx = world.getActivityIndex(act_name);
        if (act_idx >= 0) {
          slot.property_hop_dispatch_by_activity_idx[static_cast<int16_t>(
              act_idx)] = dispatch;
        }
      }
    }
  };

  for (auto& sched_type : schedule_types) {
    // Resolve force_hybrid_mask and linked_activities_mask. Activities listed
    // in linked_activities are implicitly force-hybrid (must be re-rolled at
    // runtime to honour the daily cached decision).
    sched_type.force_hybrid_mask = 0;
    for (const auto& act_name : sched_type.force_hybrid_activities) {
      int idx = world.getActivityIndex(act_name);
      if (idx >= 0) {
        sched_type.force_hybrid_mask |= (ActivityMask(1) << idx);
      }
    }
    sched_type.linked_activities_mask = 0;
    for (const auto& act_name : sched_type.linked_activities) {
      int idx = world.getActivityIndex(act_name);
      if (idx >= 0) {
        ActivityMask bit = (ActivityMask(1) << idx);
        sched_type.linked_activities_mask |= bit;
        sched_type.force_hybrid_mask |= bit;  // implies force_hybrid
      }
    }

    // Build participation_by_day_type_id[dt_idx][act_idx]
    sched_type.participation_by_day_type_id.assign(
        num_dt, std::vector<double>(num_acts, 0.0));
    for (int dt_idx = 0; dt_idx < num_dt; ++dt_idx) {
      const std::string& dt_name = day_type_names[dt_idx];
      auto it = sched_type.participation_by_day_type.find(dt_name);
      if (it == sched_type.participation_by_day_type.end()) continue;
      for (const auto& [act_name, rate] : it->second) {
        int act_idx = world.getActivityIndex(act_name);
        if (act_idx >= 0 && act_idx < static_cast<int>(num_acts)) {
          sched_type.participation_by_day_type_id[dt_idx][act_idx] = rate;
        }
      }
    }

    // Build slots_by_day_type_idx[dt_idx] and resolve slot caches
    sched_type.slots_by_day_type_idx.assign(num_dt, nullptr);
    for (int dt_idx = 0; dt_idx < num_dt; ++dt_idx) {
      const std::string& dt_name = day_type_names[dt_idx];
      auto it = sched_type.slots_by_day_type.find(dt_name);
      if (it == sched_type.slots_by_day_type.end()) continue;
      resolveSlotVec(it->second);
      sched_type.slots_by_day_type_idx[dt_idx] = &it->second;
    }

    // Resolve flat_slots for temporary schedules
    resolveSlotVec(sched_type.flat_slots);

    // A temporary schedule with no flat_slots is malformed: every hop consumer
    // divides by flat_slots.size() (hopStartDay, advanceAndCheckComplete), so
    // an empty one would divide by zero mid-simulation. Fail fast at load
    // instead.
    if (sched_type.is_temporary && sched_type.flat_slots.empty()) {
      throw std::runtime_error("temporary schedule '" + sched_type.name +
                               "' has no flat_slots");
    }
  }

  // Resolve return_schedule_idx for temporary schedules
  for (auto& sched_type : schedule_types) {
    if (!sched_type.return_schedule.empty()) {
      auto it = std::find_if(schedule_types.begin(), schedule_types.end(),
                             [&](const ScheduleType& s) {
                               return s.name == sched_type.return_schedule;
                             });
      if (it != schedule_types.end())
        sched_type.return_schedule_idx =
            static_cast<int16_t>(it - schedule_types.begin());
    }
  }
}

void PerformanceConfig::resolve(const WorldState& world) {
  deterministic_mask = computeActivityMaskFromNames(deterministic_activities,
                                                    world.activity_names);
  hybrid_mask =
      computeActivityMaskFromNames(hybrid_activities, world.activity_names);
  stochastic_mask =
      computeActivityMaskFromNames(stochastic_activities, world.activity_names);
  masks_resolved = true;
}

void CoordinatedEncounterConfig::resolve(
    WorldState& world, ContactMatrixConfig& contact_matrices) {
  if (!enabled) return;

  // Build deterministic name→id mapping from sorted keys
  {
    std::vector<std::string> sorted_names;
    for (auto& [name, _] : contact_matrices.matrices)
      sorted_names.push_back(name);
    std::sort(sorted_names.begin(), sorted_names.end());
    for (int i = 0; i < (int)sorted_names.size(); ++i)
      contact_matrices.matrix_name_to_id[sorted_names[i]] = i;
  }

  // Populate encounter_type_names in WorldState if not already there
  for (auto& enc : encounters) {
    if (std::find(world.encounter_type_names.begin(),
                  world.encounter_type_names.end(),
                  enc.name) == world.encounter_type_names.end()) {
      world.encounter_type_names.push_back(enc.name);
    }

    // Build encounter_type_id -> virtual_contact_matrix name mapping
    if (enc.is_virtual && !enc.virtual_contact_matrix.empty()) {
      int type_id = world.getEncounterTypeIndex(enc.name);
      if (type_id >= 0) {
        contact_matrices.virtual_matrix_names[static_cast<uint8_t>(type_id)] =
            enc.virtual_contact_matrix;
      }
    }

    // Pre-resolve caches for each encounter def
    enc.cached_encounter_type_id =
        static_cast<uint8_t>(world.getEncounterTypeIndex(enc.name));
    enc.trigger_mask =
        computeActivityMaskFromNames(enc.trigger_slots, world.activity_names);

    // allowed_venue_mask: bit i set if venue_type_names[i] is in allowed_venues
    enc.allowed_venue_mask = 0;
    for (const auto& vname : enc.allowed_venues) {
      int idx = world.getVenueTypeIndex(vname);
      if (idx >= 0) {
        enc.allowed_venue_mask |= (ActivityMask(1) << idx);
      }
    }

    // Network resolution: looked up from the static network registry.
    enc.cached_network_idx = world.getNetworkTypeIndex(enc.network);

    // Hop-schedule resolution: -1 (no hop) if hop_schedule is empty or
    // names an unknown schedule type.
    if (!enc.hop_schedule.empty()) {
      enc.cached_hop_schedule_idx =
          static_cast<int16_t>(world.getScheduleTypeIndex(enc.hop_schedule));
    }

    // Trigger-activity resolution: used to stamp activity_index on
    // injection so it reflects the trigger activity even when that
    // activity is coordinated_only (never chosen by ordinary resolution).
    if (!enc.trigger_slots.empty()) {
      enc.cached_trigger_activity_idx =
          static_cast<int16_t>(world.getActivityIndex(enc.trigger_slots[0]));
    }

    // Virtual venue type ID: use deterministic sorted registry
    if (enc.is_virtual && !enc.virtual_contact_matrix.empty()) {
      auto id_it =
          contact_matrices.matrix_name_to_id.find(enc.virtual_contact_matrix);
      if (id_it != contact_matrices.matrix_name_to_id.end()) {
        enc.cached_virtual_venue_type_id = id_it->second;
      }
    }
  }
}

const ScheduleType* ScheduleConfig::tryCSVAssignment(const Person& person,
                                                     const WorldState& world,
                                                     std::mt19937& rng) const {
  if (csv_rows.empty()) return nullptr;

  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for (const auto& row : csv_rows) {
    if (filtering::matchesCriteria(person, &world, row.criteria)) {
      double rand_val = dist(rng);
      for (const auto& [type_idx, cumulative] : row.schedule_probs) {
        if (rand_val < cumulative) {
          return &schedule_types[type_idx];
        }
      }
      // rand_val >= last cumulative bound → fallback to YAML
      return nullptr;
    }
  }
  return nullptr;  // No row matched
}

void ScheduleConfig::resolveCSV(const WorldState& world) {
  if (csv_path.empty()) return;
  csv_rows.clear();

  csv::FilteredTable table = csv::loadFilteredCSV(csv_path);

  // Classify value columns: schedule.<name> → schedule_type index, plus the
  // geo_level / geo_unit pair used for world-graph BFS.
  std::vector<std::pair<std::string, int>>
      schedule_cols;  // (value-col name, schedule_type_idx)
  for (const auto& col : table.value_columns) {
    if (col.find("schedule.") == 0) {
      std::string sched_name = col.substr(9);
      int idx = -1;
      for (int j = 0; j < (int)schedule_types.size(); ++j) {
        if (schedule_types[j].name == sched_name) {
          idx = j;
          break;
        }
      }
      if (idx < 0) {
        throw std::runtime_error("Schedule '" + sched_name +
                                 "' referenced in column 'schedule." +
                                 sched_name + "' not found in schedule_types");
      }
      schedule_cols.push_back({col, idx});
    }
  }

  auto get = [](const csv::FilteredRow& r,
                const std::string& name) -> std::string {
    auto it = r.values.find(name);
    return it == r.values.end() ? "" : it->second;
  };

  int row_num = 0;
  for (const auto& r : table.rows) {
    ++row_num;
    ScheduleAssignmentRow row;
    row.criteria = r.criteria;

    // Geo filter: resolve (geo_level, geo_unit) → "in" criterion over SGU ids
    std::string geo_level_val = get(r, "geo_level");
    std::string geo_unit_val = get(r, "geo_unit");
    if (geo_level_val.empty() != geo_unit_val.empty()) {
      throw std::runtime_error(
          "Schedule assignment CSV '" + csv_path + "' row " +
          std::to_string(row_num) +
          " has only one of geo_level/geo_unit set (both or neither required); "
          "got geo_level='" +
          geo_level_val + "', geo_unit='" + geo_unit_val + "'");
    }
    if (!geo_level_val.empty() && !geo_unit_val.empty()) {
      GeoUnitId found_id = -1;
      for (const auto& unit : world.geo_units) {
        if (unit.level_id < (int)world.geo_level_names.size() &&
            world.geo_level_names[unit.level_id] == geo_level_val &&
            unit.name == geo_unit_val) {
          found_id = unit.id;
          break;
        }
      }
      if (found_id == -1) {
        throw std::runtime_error(
            "Schedule assignment CSV '" + csv_path + "' row " +
            std::to_string(row_num) + ": geo_level='" + geo_level_val +
            "', geo_unit='" + geo_unit_val + "' not found in world");
      }
      // BFS to collect all descendant geo_unit_ids (including self)
      std::vector<int32_t> descendant_ids;
      std::vector<GeoUnitId> to_visit = {found_id};
      while (!to_visit.empty()) {
        GeoUnitId current = to_visit.back();
        to_visit.pop_back();
        descendant_ids.push_back(static_cast<int32_t>(current));
        for (const auto& unit : world.geo_units) {
          if (unit.parent_id == current) to_visit.push_back(unit.id);
        }
      }
      SelectionCriterion geo_crit;
      geo_crit.property_path = "geo_unit_id";
      geo_crit.operator_type = "in";
      geo_crit.value = descendant_ids;
      row.criteria.push_back(geo_crit);
    }

    // Parse schedule probabilities
    double sum = 0.0;
    std::vector<std::pair<int, double>> raw_probs;
    for (const auto& [col_name, type_idx] : schedule_cols) {
      std::string v = get(r, col_name);
      if (v.empty()) continue;
      double prob;
      try {
        prob = std::stod(v);
      } catch (const std::exception&) {
        throw std::runtime_error("Schedule assignment CSV '" + csv_path +
                                 "' row " + std::to_string(row_num) +
                                 " column '" + col_name +
                                 "' has non-numeric value '" + v + "'");
      }
      if (prob > 0.0) {
        raw_probs.push_back({type_idx, prob});
        sum += prob;
      }
    }

    if (sum > 1.02) {
      throw std::runtime_error("Schedule probabilities in CSV row " +
                               std::to_string(row_num) + " sum to " +
                               std::to_string(sum) + ", must not exceed 1.02");
    }

    if (sum >= 0.98) {
      if (std::abs(sum - 1.0) > 1e-6) {
        std::cout << "Warning: schedule probabilities in row " << row_num
                  << " sum to " << sum << ", normalizing to 1.0" << std::endl;
      }
      row.fallback_prob = 0.0;
      for (auto& [type_idx, prob] : raw_probs) prob /= sum;
    } else {
      row.fallback_prob = 1.0 - sum;
    }

    double cumulative = 0.0;
    for (const auto& [type_idx, prob] : raw_probs) {
      cumulative += prob;
      row.schedule_probs.push_back({type_idx, cumulative});
    }

    csv_rows.push_back(row);
  }

  std::cout << "Loaded " << csv_rows.size()
            << " schedule assignment rows from: " << csv_path << std::endl;
}

void Config::resolve(WorldState& world) {
  checkConfigConsistency(*this, world);
  simulation.resolve(world);
  schedule.resolve(world);
  schedule.resolveCSV(world);
  activity_preferences.resolve(world);
  vaccination.resolve(world);
  contact_matrices.resolve(world);
  coordinated_encounters.resolve(world, contact_matrices);
  performance.resolve(world);
}

}  // namespace june
