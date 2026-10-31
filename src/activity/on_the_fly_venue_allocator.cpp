#include "activity/on_the_fly_venue_allocator.h"

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <stdexcept>

#include "core/world_state.h"
#include "loaders/config_loader_detail.h"

namespace june {

const std::vector<VenueId> OnTheFlyVenueAllocator::empty_pool_{};

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

OnTheFlyVenueAllocator::OnTheFlyVenueAllocator(const std::string& config_path)
    : OnTheFlyVenueAllocator(YAML::LoadFile(config_path)) {}

OnTheFlyVenueAllocator::OnTheFlyVenueAllocator(const YAML::Node& root) {
  std::unordered_map<std::string, RuleConfig> parsed_rules;
  config_detail::rejectUnknownKeys(root, {"rules", "activity_rules"},
                                   "on_the_fly_venues: ");
  const auto& rules_node = root["rules"];
  for (const auto& entry : rules_node) {
    const std::string rule_name = entry.first.as<std::string>();
    const YAML::Node& cfg = entry.second;

    config_detail::rejectUnknownKeys(
        cfg, {"strategy", "venue_type", "venue_stability", "geo_unit_level"},
        "on_the_fly_venues: ", " in rule '" + rule_name + "'");

    RuleConfig rule;

    const std::string strategy_str = cfg["strategy"].as<std::string>();
    if (strategy_str == "hosting_geo_unit") {
      rule.strategy = Strategy::hosting_geo_unit;
    } else if (strategy_str == "resident_geo_unit") {
      rule.strategy = Strategy::resident_geo_unit;
    } else {
      throw std::runtime_error("Unknown OTF strategy: " + strategy_str);
    }

    rule.venue_type = cfg["venue_type"].as<std::string>();

    if (cfg["venue_stability"]) {
      const std::string stability = cfg["venue_stability"].as<std::string>();
      if (stability == "fixed") {
        rule.venue_stability = VenueStability::fixed;
      } else if (stability == "daily") {
        rule.venue_stability = VenueStability::daily;
      } else {
        throw std::runtime_error("Unknown venue_stability: " + stability);
      }
    }

    if (cfg["geo_unit_level"]) {
      rule.geo_unit_level = cfg["geo_unit_level"].as<std::string>();
    }

    parsed_rules[rule_name] = std::move(rule);
  }

  const auto& activity_rules_node = root["activity_rules"];
  for (const auto& entry : activity_rules_node) {
    const std::string activity_name = entry.first.as<std::string>();
    const std::string rule_name = entry.second.as<std::string>();
    auto rule_it = parsed_rules.find(rule_name);
    if (rule_it == parsed_rules.end()) {
      throw std::runtime_error("Unknown OTF rule: " + rule_name);
    }
    rules_[activity_name] = rule_it->second;
  }
}

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

void OnTheFlyVenueAllocator::checkConsistency(const WorldState& world) const {
  for (const auto& [activity_name, rule] : rules_) {
    if (rule.geo_unit_level.empty()) continue;
    bool found = false;
    for (const auto& level : world.geo_level_names) {
      if (level == rule.geo_unit_level) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::string known;
      for (const auto& level : world.geo_level_names) {
        if (!known.empty()) known += ", ";
        known += '"' + level + '"';
      }
      throw std::runtime_error(
          "OTF activity '" + activity_name + "' specifies geo_unit_level '" +
          rule.geo_unit_level +
          "' which is not a known geo level. Known levels: " + known + ".");
    }
  }
}

bool OnTheFlyVenueAllocator::hasRule(std::string_view activity_name) const {
  return rules_.count(activity_name) > 0;
}

bool OnTheFlyVenueAllocator::isFixed(std::string_view activity_name) const {
  auto rule_it = rules_.find(activity_name);
  return rule_it != rules_.end() &&
         rule_it->second.venue_stability == VenueStability::fixed;
}

const std::vector<VenueId>& OnTheFlyVenueAllocator::resolve(
    std::string_view activity_name, const VenueResolveContext& context,
    const WorldState& world) {
  auto rule_it = rules_.find(activity_name);
  if (rule_it == rules_.end()) return empty_pool_;
  const std::string& activity_key = rule_it->first;
  const RuleConfig& rule = rule_it->second;

  GeoUnitId geo_unit_id;
  if (rule.strategy == Strategy::hosting_geo_unit) {
    if (!context.hosting_geo_unit_id) return empty_pool_;
    geo_unit_id = *context.hosting_geo_unit_id;
  } else {
    geo_unit_id = context.resident_geo_unit_id;
    if (!rule.geo_unit_level.empty()) {
      GeoUnitId ancestor =
          world.ancestorAtLevel(geo_unit_id, rule.geo_unit_level);
      if (ancestor == -1) {
        std::cerr << "[OTF] Warning: activity '" << activity_name
                  << "': no ancestor at level '" << rule.geo_unit_level
                  << "' for geo_unit_id " << geo_unit_id
                  << ". Returning empty pool.\n";
        return empty_pool_;
      }
      geo_unit_id = ancestor;
    }
  }

  auto rule_cache_it = cache_.find(activity_key);
  if (rule_cache_it != cache_.end()) {
    auto cache_it = rule_cache_it->second.find(geo_unit_id);
    if (cache_it != rule_cache_it->second.end()) return cache_it->second;
  }

  // After sealing, the global venue maps have been freed, so every real query
  // must already be cached. A miss means precomputeAllPools() did not cover
  // this (rule, geo); throw instead of returning an empty venue pool.
  if (sealed_)
    throw std::runtime_error(
        "OnTheFlyVenueAllocator: pool not precomputed for activity '" +
        std::string(activity_name) + "' geo_unit " +
        std::to_string(geo_unit_id) + " (global venue maps already freed)");

  auto pool = world.getVenuesInGeoUnit(geo_unit_id, rule.venue_type);
  auto& rule_cache = cache_[activity_key];
  auto [inserted_it, _] = rule_cache.emplace(geo_unit_id, std::move(pool));
  return inserted_it->second;
}

void OnTheFlyVenueAllocator::precomputeAllPools(
    const WorldState& world, const std::vector<GeoUnitId>& hosting_geo_units) {
  for (const auto& [activity, rule] : rules_) {
    if (rule.strategy == Strategy::hosting_geo_unit) {
      for (GeoUnitId g : hosting_geo_units) {
        VenueResolveContext ctx;
        ctx.hosting_geo_unit_id = g;
        resolve(activity, ctx, world);  // fills cache_
      }
    } else {  // resident_geo_unit: every resident geo is a possible query
      for (const Person& p : world.people) {
        VenueResolveContext ctx;
        ctx.resident_geo_unit_id = p.geo_unit_id;
        resolve(activity, ctx, world);  // cache_ dedups by (rule, geo)
      }
    }
  }
  sealed_ = true;  // from now on a cache miss in resolve() throws
}

}  // namespace june
