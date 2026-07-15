#include "epidemiology/infection_seed.h"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "epidemiology/disease.h"
#include "utils/deterministic_rng.h"
#include "utils/filtered_csv.h"
#include "utils/random.h"

namespace june {

// =============================================================================
// Configuration Loader Implementation
// =============================================================================

std::vector<SelectionCriterion> InfectionSeedConfigLoader::parseCriterion(
    const std::string& key, const std::string& val) {
  std::string property = key;
  // Map age_groups to age for internal property evaluation
  if (property == "age_groups") {
    property = "age";
    // Silent mapping: internal config normalization
  }

  std::vector<SelectionCriterion> results;

  // Support ranges like "18-30" or "65-100"
  size_t dash = val.find('-');
  if (dash != std::string::npos && dash > 0 && dash < val.size() - 1) {
    try {
      SelectionCriterion c_min, c_max;
      c_min.property_path = property;
      c_min.operator_type = ">=";
      c_min.value = std::stoi(val.substr(0, dash));

      c_max.property_path = property;
      c_max.operator_type = "<=";
      c_max.value = std::stoi(val.substr(dash + 1));

      results.push_back(c_min);
      results.push_back(c_max);
      return results;
    } catch (...) {
      // If parsing failed, fall back to treats it as a single string
    }
  }

  SelectionCriterion c;
  c.property_path = property;
  c.operator_type = "==";

  // Try numeric/bool first, fallback to string
  try {
    if (val == "true" || val == "True")
      c.value = true;
    else if (val == "false" || val == "False")
      c.value = false;
    else if (val.find('.') != std::string::npos) {
      c.value = std::stod(val);
    } else {
      c.value = std::stoi(val);
    }
  } catch (...) {
    c.value = val;
  }

  results.push_back(c);
  return results;
}

InfectionSeedType InfectionSeedConfigLoader::parseSeedType(
    const std::string& type_str) {
  std::string lower_type = type_str;
  std::transform(lower_type.begin(), lower_type.end(), lower_type.begin(),
                 ::tolower);

  if (lower_type == "uniform") return InfectionSeedType::UNIFORM;
  if (lower_type == "exact") return InfectionSeedType::EXACT;
  if (lower_type == "clustered") return InfectionSeedType::CLUSTERED;

  throw std::runtime_error("Unknown infection seed type: " + type_str);
}

void InfectionSeedConfigLoader::loadBulkCsvSeeds(const std::string& csv_path,
                                                 InfectionSeedConfig& config) {
  config.bulk_csv_path = csv_path;
  csv::FilteredTable table = csv::loadFilteredCSV(csv_path);

  auto has_col = [&](const std::string& name) {
    for (const auto& c : table.value_columns)
      if (c == name) return true;
    return false;
  };
  if (!has_col("name") || !has_col("date") || !has_col("type")) {
    throw std::runtime_error(
        "Bulk CSV missing required columns: name, date, type");
  }

  auto get = [](const csv::FilteredRow& r,
                const std::string& name) -> std::string {
    auto it = r.values.find(name);
    return it == r.values.end() ? "" : it->second;
  };

  struct SeedKey {
    std::string name;
    std::string date;
    InfectionSeedType type;
    std::string trajectory_key;
    std::string start_symptom;
    bool operator<(const SeedKey& o) const {
      return std::tie(name, date, type, trajectory_key, start_symptom) <
             std::tie(o.name, o.date, o.type, o.trajectory_key,
                      o.start_symptom);
    }
  };

  struct SeedDraft {
    InfectionSeedEvent event;
    std::vector<std::pair<std::vector<SelectionCriterion>, size_t>> profiles;
    std::map<std::string, std::vector<std::pair<size_t, int>>> unit_pending;

    size_t getOrCreateGroup(const std::vector<SelectionCriterion>& profile) {
      for (const auto& [p, idx] : profiles) {
        if (p.size() != profile.size()) continue;
        bool match = true;
        for (size_t i = 0; i < p.size(); ++i) {
          const auto& a = p[i];
          const auto& b = profile[i];
          if (a.property_path != b.property_path ||
              a.operator_type != b.operator_type || a.value != b.value) {
            match = false;
            break;
          }
        }
        if (match) return idx;
      }
      size_t new_idx = event.structured_config.target_groups.size();
      profiles.push_back({profile, new_idx});
      SeedTargetGroup group;
      group.criteria = profile;
      event.structured_config.target_groups.push_back(group);
      return new_idx;
    }
  };

  std::map<SeedKey, SeedDraft> drafts;

  int row_num = 0;
  for (const auto& row : table.rows) {
    ++row_num;
    std::string name_val = get(row, "name");
    std::string date_val = get(row, "date");
    std::string type_val = get(row, "type");
    if (name_val.empty() || date_val.empty() || type_val.empty()) {
      throw std::runtime_error(
          "Bulk seed CSV '" + csv_path + "' row " + std::to_string(row_num) +
          " missing required value in name/date/type (got name='" + name_val +
          "', date='" + date_val + "', type='" + type_val + "')");
    }

    SeedKey key = {name_val, date_val, parseSeedType(type_val),
                   get(row, "trajectory_key"), get(row, "start_symptom")};
    auto& draft = drafts[key];
    if (draft.event.name.empty()) {
      draft.event.name = key.name;
      draft.event.date_time = key.date;
      draft.event.type = key.type;
      draft.event.trajectory_key = key.trajectory_key;
      draft.event.start_symptom = key.start_symptom;
    }

    if (key.type == InfectionSeedType::UNIFORM) {
      std::string pc = get(row, "cases_per_capita");
      if (!pc.empty()) {
        try {
          draft.event.uniform_config.cases_per_capita = std::stod(pc);
        } catch (const std::exception&) {
          throw std::runtime_error("Bulk seed CSV '" + csv_path + "' row " +
                                   std::to_string(row_num) +
                                   " has non-numeric cases_per_capita='" + pc +
                                   "'");
        }
      }
      draft.event.attribute_filters.insert(draft.event.attribute_filters.end(),
                                           row.criteria.begin(),
                                           row.criteria.end());
    } else {
      std::string geo_level = get(row, "geo_level");
      if (!geo_level.empty()) {
        draft.event.structured_config.geo_level = geo_level;
      }

      std::string geo_unit = get(row, "geo_unit");
      int cases = 0;
      std::string cases_val = get(row, "cases");
      if (!cases_val.empty()) {
        try {
          cases = std::stoi(cases_val);
        } catch (const std::exception&) {
          throw std::runtime_error("Bulk seed CSV '" + csv_path + "' row " +
                                   std::to_string(row_num) +
                                   " has non-integer cases='" + cases_val +
                                   "'");
        }
      }

      size_t profile_idx = draft.getOrCreateGroup(row.criteria);

      if (!geo_unit.empty() && cases > 0) {
        draft.unit_pending[geo_unit].push_back({profile_idx, cases});
      }
    }
  }

  // Finalize structured seeds
  for (auto& [key, draft] : drafts) {
    if (draft.event.type != InfectionSeedType::UNIFORM) {
      for (const auto& [unit_id, groups] : draft.unit_pending) {
        UnitCases uc;
        uc.unit_id = unit_id;
        uc.cases_per_target_group.resize(
            draft.event.structured_config.target_groups.size(), 0);
        for (const auto& [g_idx, count] : groups) {
          uc.cases_per_target_group[g_idx] += count;
        }
        draft.event.structured_config.unit_cases.push_back(uc);
      }
    }
    config.seeds.push_back(draft.event);
  }
  // CSV parse summary removed (not actionable)
}

InfectionSeedConfig InfectionSeedConfigLoader::loadFromFile(
    const std::string& filename) {
  InfectionSeedConfig config;

  try {
    YAML::Node root = YAML::LoadFile(filename);

    if (root["global_parameters"]) {
      auto global = root["global_parameters"];
      if (global["base_cases_per_capita"]) {
        config.global_params.base_cases_per_capita =
            global["base_cases_per_capita"].as<double>();
      }
    }

    if (root["bulk_csv"]) {
      loadBulkCsvSeeds(root["bulk_csv"].as<std::string>(), config);
    }

    if (root["infection_seeds"]) {
      for (const auto& seed_node : root["infection_seeds"]) {
        InfectionSeedEvent seed;
        seed.name = seed_node["name"].as<std::string>();
        seed.type = parseSeedType(seed_node["type"].as<std::string>());
        seed.date_time = seed_node["date"].as<std::string>();

        if (seed_node["trajectory_key"])
          seed.trajectory_key = seed_node["trajectory_key"].as<std::string>();
        if (seed_node["start_symptom"])
          seed.start_symptom = seed_node["start_symptom"].as<std::string>();

        if (seed_node["parameters"]) {
          auto params = seed_node["parameters"];
          seed.seed_strength = params["seed_strength"].as<double>(
              config.global_params.default_seed_strength);

          if (params["attribute_filters"]) {
            auto filters = params["attribute_filters"];
            for (auto it = filters.begin(); it != filters.end(); ++it) {
              auto cs = parseCriterion(it->first.as<std::string>(),
                                       it->second.as<std::string>());
              seed.attribute_filters.insert(seed.attribute_filters.end(),
                                            cs.begin(), cs.end());
            }
          }
        }

        if (seed.type == InfectionSeedType::UNIFORM) {
          if (seed_node["parameters"]) {
            auto params = seed_node["parameters"];
            if (params["cases_per_capita_multiplier"]) {
              seed.uniform_config.cases_per_capita =
                  config.global_params.base_cases_per_capita *
                  params["cases_per_capita_multiplier"].as<double>();
            }
          }
        } else {
          seed.structured_config.geo_level =
              seed_node["geo_level"].as<std::string>("MGU");

          if (seed_node["parameters"]) {
            auto params = seed_node["parameters"];

            // Age groups from YAML are converted to generic TargetGroups
            if (params["age_groups"]) {
              for (const auto& age_str : params["age_groups"]) {
                std::string s = age_str.as<std::string>();
                SeedTargetGroup g;
                auto cs = parseCriterion("age", s);
                g.criteria.insert(g.criteria.end(), cs.begin(), cs.end());
                seed.structured_config.target_groups.push_back(g);
              }
            }

            if (params["units"]) {
              for (const auto& entry : params["units"]) {
                UnitCases uc;
                uc.unit_id = entry.first.as<std::string>();
                if (entry.second.IsSequence()) {
                  for (const auto& cases : entry.second)
                    uc.cases_per_target_group.push_back(
                        static_cast<int>(cases.as<double>()));
                } else {
                  int c = static_cast<int>(entry.second.as<double>());
                  uc.cases_per_target_group.assign(
                      seed.structured_config.target_groups.size(), c);
                }
                seed.structured_config.unit_cases.push_back(uc);
              }
            }
          }
        }
        config.seeds.push_back(seed);
      }
    }
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to load infection seed config: " +
                             std::string(e.what()));
  }
  return config;
}

// =============================================================================
// Infection Seeder Implementation
// =============================================================================

InfectionSeeder::InfectionSeeder(WorldState& world, const Disease* disease,
                                 const InfectionSeedConfig& config,
                                 EventLogger* event_logger, uint64_t base_seed)
    : world_(world),
      disease_(disease),
      config_(config),
      event_logger_(event_logger),
      current_simulation_time_(0.0),
      base_seed_(base_seed) {}

std::vector<PersonId> InfectionSeeder::seedInfections(
    const std::string& current_datetime, double simulation_time) {
  current_simulation_time_ = simulation_time;
  std::vector<PersonId> all_infected;

  for (const auto& seed : config_.seeds) {
    // Standardized comparison: skip whitespace/case if needed,
    // though currently matching exact string.
    if (seed.date_time == current_datetime) {
      std::string seed_key =
          seed.name + "|" + seed.trajectory_key + "|" + seed.start_symptom;
      if (applied_seeds_.count(seed_key) > 0) {
        continue;
      }
      std::vector<PersonId> infected = applySeed(seed);
      applied_seeds_.insert(seed_key);
      all_infected.insert(all_infected.end(), infected.begin(), infected.end());

      // Per-seed message removed; global count reported by Simulator
    }
  }

  return all_infected;
}

std::vector<PersonId> InfectionSeeder::applySeed(
    const InfectionSeedEvent& seed) {
  switch (seed.type) {
    case InfectionSeedType::UNIFORM:
      return applyUniformSeed(seed);
    case InfectionSeedType::EXACT:
      return applyExactSeed(seed);
    case InfectionSeedType::CLUSTERED:
      return applyClusteredSeed(seed);
    default:
      throw std::runtime_error("Unknown seed type");
  }
}

std::vector<PersonId> InfectionSeeder::applyUniformSeed(
    const InfectionSeedEvent& seed) {
  std::vector<PersonId> infected_ids;

  double cases_per_capita =
      seed.uniform_config.cases_per_capita * seed.seed_strength;

  // Target rate message removed (not actionable)

  // MPI-reproducible seeding: each person gets a per-person deterministic
  // decision based on their ID. This ensures the same person is always
  // seeded regardless of which rank owns them or the local population size.
  uint64_t seed_name_hash = std::hash<std::string>{}(seed.name);
  uint64_t time_bits = static_cast<uint64_t>(current_simulation_time_ * 1000);

  for (auto& person : world_.people) {
    if (person.infection != nullptr) continue;
    if (person.getSusceptibility(current_simulation_time_,
                                 disease_->getName()) < 0.01)
      continue;
    if (!matchesAttributes(&person, seed.attribute_filters)) continue;

    // Per-person deterministic draw keyed to person ID
    SplitMix64 prng(mix_seed(base_seed_, person.id, seed_name_hash, time_bits));
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double rng_val = dist(prng);
    bool seeded = rng_val < cases_per_capita;
    if (seeded) {
      infectPerson(&person, seed.trajectory_key, seed.start_symptom);
      if (person.infection != nullptr) {
        infected_ids.push_back(person.id);
      }
    }
  }

  // Seeded count removed; global count reported by Simulator

  return infected_ids;
}

std::vector<PersonId> InfectionSeeder::applyExactSeed(
    const InfectionSeedEvent& seed) {
  std::vector<PersonId> infected_ids;

  for (const auto& unit_case : seed.structured_config.unit_cases) {
    // Step 1: Get people in this unit
    std::vector<Person*> candidates_in_unit = world_.getPeopleInUnit(
        seed.structured_config.geo_level, unit_case.unit_id);

    if (candidates_in_unit.empty()) {
      std::cerr << "[ExactSeed] Warning: No people found in unit '"
                << unit_case.unit_id << "' at level '"
                << seed.structured_config.geo_level << "'" << std::endl;
      continue;
    }

    // Step 2: Apply seeds for each target group
    for (size_t g_idx = 0; g_idx < seed.structured_config.target_groups.size();
         ++g_idx) {
      const auto& group = seed.structured_config.target_groups[g_idx];
      int num_cases = static_cast<int>(unit_case.cases_per_target_group[g_idx] *
                                       seed.seed_strength);
      if (num_cases <= 0) continue;

      // Filter people in this unit by group and global filters
      std::vector<Person*> candidates;
      for (auto* person : candidates_in_unit) {
        if (person->infection == nullptr &&
            person->getSusceptibility(current_simulation_time_,
                                      disease_->getName()) >= 0.01 &&
            matchesAttributes(person, seed.attribute_filters) &&
            group.matches(*person, &world_)) {
          candidates.push_back(person);
        }
      }

      if (candidates.empty()) {
        continue;
      }

      size_t start_idx = infected_ids.size();
      uint64_t unit_hash = std::hash<std::string>{}(unit_case.unit_id);
      SplitMix64 exact_rng(mix_seed(base_seed_, unit_hash, g_idx, 0xE4AC7));
      std::sort(candidates.begin(), candidates.end(),
                [](const Person* a, const Person* b) { return a->id < b->id; });
      std::shuffle(candidates.begin(), candidates.end(), exact_rng);
      for (int i = 0; i < num_cases && i < (int)candidates.size(); ++i) {
        infectPerson(candidates[i], seed.trajectory_key, seed.start_symptom);
        if (candidates[i]->infection != nullptr) {
          infected_ids.push_back(candidates[i]->id);
        }
      }
      // Per-unit success message removed; global count reported by Simulator
    }
  }
  return infected_ids;
}

std::vector<PersonId> InfectionSeeder::applyClusteredSeed(
    const InfectionSeedEvent& seed) {
  std::vector<PersonId> infected_ids;

  for (const auto& unit_case : seed.structured_config.unit_cases) {
    std::vector<Person*> unit_people_raw = world_.getPeopleInUnit(
        seed.structured_config.geo_level, unit_case.unit_id);
    std::vector<Person*> unit_candidates;
    for (Person* person : unit_people_raw) {
      if (person->infection == nullptr &&
          person->getSusceptibility(current_simulation_time_,
                                    disease_->getName()) >= 0.01 &&
          matchesAttributes(person, seed.attribute_filters)) {
        unit_candidates.push_back(person);
      }
    }

    if (unit_candidates.empty()) continue;

    std::vector<int> target_infections(
        seed.structured_config.target_groups.size());
    int total_target = 0;
    for (size_t i = 0; i < target_infections.size(); ++i) {
      target_infections[i] = static_cast<int>(
          unit_case.cases_per_target_group[i] * seed.seed_strength);
      total_target += target_infections[i];
    }

    if (total_target <= 0) continue;

    // Group people by household
    std::map<VenueId, std::vector<Person*>> households;
    for (Person* person : unit_candidates) {
      auto residence = world_.getActivityVenues(*person, "residence");
      if (!residence.empty()) {
        households[residence[0].first].push_back(person);
      }
    }

    if (households.empty()) continue;

    // Scoring and selection logic
    std::vector<std::pair<VenueId, double>> pool;
    for (const auto& [hh_id, members] : households) {
      double score = 0.0;
      for (Person* p : members) {
        for (size_t g_idx = 0;
             g_idx < seed.structured_config.target_groups.size(); ++g_idx) {
          if (seed.structured_config.target_groups[g_idx].matches(*p,
                                                                  &world_)) {
            score += 1.0;  // Basic score for matching a target group
            break;
          }
        }
      }
      if (score > 0)
        pool.push_back({hh_id, score / std::sqrt((double)members.size())});
    }

    uint64_t cluster_unit_hash = std::hash<std::string>{}(unit_case.unit_id);
    SplitMix64 cluster_rng(mix_seed(base_seed_, cluster_unit_hash, 0xC1057E8));
    std::shuffle(pool.begin(), pool.end(), cluster_rng);
    std::sort(pool.begin(), pool.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    // Per-unit clustered seed message removed; global count reported by
    // Simulator

    int unit_infected = 0;
    std::vector<int> infected_per_group(
        seed.structured_config.target_groups.size(), 0);

    for (const auto& [hh_id, score] : pool) {
      if (unit_infected >= total_target) break;

      for (Person* person : households[hh_id]) {
        for (size_t g_idx = 0;
             g_idx < seed.structured_config.target_groups.size(); ++g_idx) {
          if (infected_per_group[g_idx] < target_infections[g_idx] &&
              seed.structured_config.target_groups[g_idx].matches(*person,
                                                                  &world_)) {
            infectPerson(person, seed.trajectory_key, seed.start_symptom);
            if (person->infection != nullptr) {
              infected_per_group[g_idx]++;
              unit_infected++;
              infected_ids.push_back(person->id);
            }
            break;
          }
        }
      }
    }
  }
  return infected_ids;
}

bool InfectionSeeder::matchesAttributes(
    const Person* person, const std::vector<SelectionCriterion>& filters) {
  if (filters.empty()) return true;
  for (const auto& filter : filters) {
    if (!filter.evaluate(*person, &world_)) return false;
  }
  return true;
}

void InfectionSeeder::infectPerson(Person* person,
                                   const std::string& trajectory_key,
                                   const std::string& start_symptom) {
  if (person->infection != nullptr) return;
  if (person->getSusceptibility(current_simulation_time_, disease_->getName()) <
      0.01)
    return;

  uint64_t infection_seed =
      mix_seed(base_seed_, person->id,
               static_cast<uint64_t>(current_simulation_time_ * 1000), 0x5EED);

  float severity_factor = 1.0f;
  auto* gu = world_.getGeoUnit(person->geo_unit_id);
  if (gu) severity_factor = gu->severity_factor;

  person->infection = std::make_unique<Infection>(
      disease_, current_simulation_time_, person,
      static_cast<unsigned int>(infection_seed), &world_,
      "seed",  // venue type
      INFECTION_SEED_VENUE_ID, severity_factor,
      0,  // infector_symptom_id -- no infector for seeds
      trajectory_key, start_symptom);

  if (event_logger_ != nullptr) {
    event_logger_->logInfection(
        person->id, kInvalidPersonId, INFECTION_SEED_VENUE_ID,
        current_simulation_time_, kDefaultEncounterTypeId,
        kNoSymptomId);  // no infector for seeds
  }

  // Per-seed audit: id, sex, age at t=seed.
  std::cout << "[SeedAudit] id=" << person->id
            << " sex=" << static_cast<int>(person->sex)
            << " age=" << person->age << " t=" << current_simulation_time_
            << std::endl;
}

}  // namespace june
