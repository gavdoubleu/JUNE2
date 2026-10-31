#pragma once

#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "activity/on_the_fly_venue_allocator.h"
#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/disease.h"

namespace june {

inline Disease makeAlwaysSickDisease() {
  TransmissionParams transmission;
  transmission.mode = InfectiousnessMode::STAGE_DRIVEN;
  auto constant_curve = std::make_shared<ConstantCurve>(1.0);
  transmission.symptom_id_curves = {nullptr, constant_curve};

  std::vector<SymptomTag> symptom_tags = {{"healthy", -1, 0}, {"sick", 1, 1}};
  DiseaseStageSettings stage_settings;
  stage_settings.recovered_stages = {"healthy"};

  TrajectoryDefinition trajectory;
  trajectory.selection_key = "general";
  trajectory.severity = 1.0;
  trajectory.stages.push_back({"sick", {"constant", {{"value", 100.0}}}});

  return Disease("TestDisease", symptom_tags, stage_settings, {trajectory}, {},
                 transmission);
}

inline WorldState makeCatchmentWorld(int num_people = 1) {
  WorldState world;
  world.geo_level_names = {"sgu"};
  GeographicalUnit gu;
  gu.id = 0;
  gu.parent_id = -1;
  gu.level_id = 0;
  world.geo_units.push_back(gu);
  world.venue_type_names = {"fair"};
  Venue v;
  v.id = 0;
  v.type_id = 0;
  v.geo_unit_id = 0;
  world.venues.push_back(v);
  world.activity_names = {"Fair_accommodation"};
  world.schedule_type_names = {"regular", "Fair_day_trip"};
  for (int i = 0; i < num_people; ++i) {
    Person& person = world.people.emplace_back();
    person.id = i;
    person.geo_unit_id = 0;
  }
  world.buildIndices();
  return world;
}

inline void resolveCalendarEventSlotIndices(TimeSlot& slot,
                                            const WorldState& world) {
  slot.allowed_activity_indices.clear();
  for (const auto& activity : slot.allowed_activities) {
    int index = world.getActivityIndex(activity);
    if (index >= 0)
      slot.allowed_activity_indices.push_back(static_cast<int16_t>(index));
  }
}

class ScopedTestFiles {
 public:
  explicit ScopedTestFiles(std::string prefix = "june2_test") {
    static std::atomic<unsigned> counter{0};
    // The process id keeps two runs of a test binary apart, and a directory
    // left under this exact name by a crashed run is cleared, not tripped over.
    directory_ = std::filesystem::temp_directory_path() /
                 (std::move(prefix) + "_" + std::to_string(::getpid()) + "_" +
                  std::to_string(counter.fetch_add(1)));
    std::filesystem::remove_all(directory_);
    if (!std::filesystem::create_directories(directory_))
      throw std::runtime_error("could not create test temp directory: " +
                               directory_.string());
  }

  ScopedTestFiles(const ScopedTestFiles&) = delete;
  ScopedTestFiles& operator=(const ScopedTestFiles&) = delete;

  ~ScopedTestFiles() {
    std::error_code ec;
    std::filesystem::remove_all(directory_, ec);
  }

  const std::filesystem::path& directory() const { return directory_; }

  std::filesystem::path path(std::string_view filename) const {
    return directory_ / std::string(filename);
  }

  std::filesystem::path write(std::string_view filename,
                              std::string_view contents) const {
    const auto path = directory_ / std::string(filename);
    std::ofstream file(path);
    if (!file.is_open())
      throw std::runtime_error("could not create test temp file: " +
                               path.string());
    file << contents;
    if (!file)
      throw std::runtime_error("could not write test temp file: " +
                               path.string());
    return path;
  }

 private:
  std::filesystem::path directory_;
};

inline OnTheFlyVenueAllocator allocatorFromYaml(std::string_view yaml,
                                                std::string_view file_prefix) {
  ScopedTestFiles files{std::string(file_prefix)};
  const auto path = files.write("config.yaml", yaml);
  auto allocator = OnTheFlyVenueAllocator(path.string());
  return allocator;
}

class TestWorldFactory {
 public:
  static WorldState createMinimalWorld(int num_people = 2, int num_venues = 1) {
    WorldState world;

    // 1. Setup Venues
    world.venue_type_names = {"office", "home", "school"};
    for (int i = 0; i < num_venues; ++i) {
      Venue v;
      v.id = i;
      v.type_id = 0;  // office
      v.geo_unit_id = 0;
      world.venues.push_back(v);
    }

    // 2. Setup Geography
    world.geo_level_names = {"city"};
    GeographicalUnit gu;
    gu.id = 0;
    gu.name = "TestCity";
    gu.level_id = 0;
    gu.parent_id = -1;
    world.geo_units.push_back(gu);

    // 3. Setup People
    world.person_property_names = {"age", "sex"};
    for (int i = 0; i < num_people; ++i) {
      Person& p = world.people.emplace_back();
      p.id = i;
      p.age = 20.0f + i;
      p.sex = (i % 2 == 0) ? Sex::MALE : Sex::FEMALE;
      p.geo_unit_id = 0;
    }

    world.buildIndices();
    return world;
  }

  /**
   * Adds a network to the world. Supports multiple networks per person:
   * the first call sets network_meta_start; subsequent calls append
   * contiguously and increment network_meta_count.
   *
   * IMPORTANT: All addNetwork calls for a given person must be made
   * consecutively (no interleaving with other people's networks) because
   * the flat array layout requires contiguous NetworkMeta entries per person.
   * In practice, call addNetwork for network A (all people), then for
   * network B (all people) — this works because each call appends to the
   * end of network_meta.
   *
   * NOTE: If you need to add networks for a person who already has some,
   * the entries MUST be contiguous. This utility handles that correctly
   * as long as you don't manually modify network_meta between calls.
   */
  static void addNetwork(
      WorldState& world, const std::string& name,
      const std::vector<std::pair<PersonId, PersonId>>& pairs) {
    int network_type_id = world.getNetworkTypeIndex(name);
    if (network_type_id < 0) {
      world.network_names.push_back(name);
      network_type_id = static_cast<int>(world.network_names.size() - 1);
    }

    for (auto& person : world.people) {
      std::vector<PersonId> partners;
      for (const auto& pair : pairs) {
        if (pair.first == person.id)
          partners.push_back(pair.second);
        else if (pair.second == person.id)
          partners.push_back(pair.first);
      }

      if (!partners.empty()) {
        uint32_t new_meta_idx =
            static_cast<uint32_t>(world.network_meta.size());

        if (person.network_meta_count == 0) {
          // First network for this person
          person.network_meta_start = new_meta_idx;
          person.network_meta_count = 1;
        } else {
          // Additional network: verify contiguity and increment count
          uint32_t expected =
              person.network_meta_start + person.network_meta_count;
          if (new_meta_idx == expected) {
            person.network_meta_count++;
          } else {
            // Non-contiguous entries cannot share the current contiguous range,
            // so this entry starts a new range.
            person.network_meta_start = new_meta_idx;
            person.network_meta_count = 1;
          }
        }

        Person::NetworkMeta meta;
        meta.network_type_id = static_cast<uint16_t>(network_type_id);
        meta.partner_start =
            static_cast<uint32_t>(world.network_partners.size());
        meta.partner_count = static_cast<uint16_t>(partners.size());

        world.network_meta.push_back(meta);
        world.network_partners.insert(world.network_partners.end(),
                                      partners.begin(), partners.end());
      }
    }
  }
};

// =============================================================================
// Partial-presence test helpers
//
// Build synthetic worlds for partial-presence FOI tests. Tests assert physics
// properties (overlap-only exposure, bin isolation, gate non-regression) on
// anonymous IDs — see feedback_tests_name_physics_not_scenarios memory.
// =============================================================================

// One (rider, leg) entry. Used by addPartialPresenceLegs to populate a
// person's ActivityMeta + activity_venues + membership_metadata in one call.
struct PartialPresenceLegSpec {
  PersonId person_id;
  VenueId venue_id;
  float t_board_min;
  float t_alight_min;
};

// Register a partial-presence venue type. Appends to world.venue_type_names,
// sets the type's target_group_size in sim_cfg.partial_presence, and flips
// the corresponding bit in enabled_venue_type_mask. Returns the new type_id.
inline uint8_t addPartialPresenceVenueType(WorldState& world,
                                           SimulationConfig& sim_cfg,
                                           const std::string& name,
                                           int target_group_size) {
  world.venue_type_names.push_back(name);
  const uint8_t type_id =
      static_cast<uint8_t>(world.venue_type_names.size() - 1);
  sim_cfg.partial_presence.target_group_size_by_name[name] = target_group_size;
  if (sim_cfg.partial_presence.target_group_size_by_type_id.size() <= type_id) {
    sim_cfg.partial_presence.target_group_size_by_type_id.resize(type_id + 1,
                                                                 0);
  }
  sim_cfg.partial_presence.target_group_size_by_type_id[type_id] =
      target_group_size;
  sim_cfg.partial_presence.enabled_venue_type_mask |= (1ULL << type_id);
  return type_id;
}

// Create a venue of the given type and append to world. Returns its VenueId.
inline VenueId addPartialPresenceVenue(WorldState& world, uint8_t type_id) {
  Venue v;
  v.id = static_cast<VenueId>(world.venues.size());
  v.type_id = type_id;
  v.parent_id = -1;
  v.geo_unit_id = 0;
  world.venues.push_back(v);
  return v.id;
}

// Attach a single multi-leg presence activity to people. Each person referenced
// in `legs` gets one ActivityMeta with venue_count = number of their legs.
// All legs become entries in activity_venues; per-leg t_board_min /
// t_alight_min are written into membership_field_values (the side-table
// the RuntimeGroupAllocator consults).
//
// Pre-condition: each referenced person has activity_meta_count == 0 on entry
// (helper is for fresh persons; reusing it on the same person would break
// the contiguous ActivityMeta range assumption).
inline void addPartialPresenceLegs(
    WorldState& world, int16_t activity_index,
    const std::vector<PartialPresenceLegSpec>& legs) {
  if (world.membership_field_names.empty()) {
    world.membership_field_names = {"t_board_min", "t_alight_min"};
    world.membership_field_values.resize(2);
  }
  const int tb_idx = world.getMembershipFieldIndex("t_board_min");
  const int ta_idx = world.getMembershipFieldIndex("t_alight_min");

  // Preserve insertion order per person.
  std::unordered_map<PersonId, std::vector<const PartialPresenceLegSpec*>>
      legs_by_pid;
  std::vector<PersonId> person_order;
  for (const auto& leg : legs) {
    auto& bucket = legs_by_pid[leg.person_id];
    if (bucket.empty()) person_order.push_back(leg.person_id);
    bucket.push_back(&leg);
  }

  for (PersonId pid : person_order) {
    auto pit = std::find_if(world.people.begin(), world.people.end(),
                            [pid](const Person& p) { return p.id == pid; });
    if (pit == world.people.end()) continue;
    Person& person = *pit;

    person.activity_meta_start =
        static_cast<uint32_t>(world.activity_meta.size());
    person.activity_meta_count = 1;

    Person::ActivityMeta meta;
    meta.activity_index = activity_index;
    meta.venue_start = static_cast<uint32_t>(world.activity_venues.size());
    meta.venue_count = static_cast<uint16_t>(legs_by_pid[pid].size());

    for (const auto* leg : legs_by_pid[pid]) {
      const uint32_t flat_idx =
          static_cast<uint32_t>(world.activity_venues.size());
      world.activity_venues.push_back({leg->venue_id, /*subset=*/0});
      world.membership_field_values[tb_idx][flat_idx] = leg->t_board_min;
      world.membership_field_values[ta_idx][flat_idx] = leg->t_alight_min;
    }
    world.activity_meta.push_back(meta);
  }
}

// Build the PersonLocation list that drives processTransmissions for a
// partial-presence scenario. One entry per (rider, leg): venue_id is the
// leg's venue, activity_index marks the rider as active this slot.
// (Multi-leg riders appear N times — once per leg — matching the engine's
// multi-presence behavior.)
inline std::vector<PersonLocation> makePartialPresenceLocations(
    const WorldState& world, int16_t activity_index) {
  std::vector<PersonLocation> locs;
  for (size_t i = 0; i < world.people.size(); ++i) {
    const Person& p = world.people[i];
    for (const auto& meta : world.getActivityMetas(p)) {
      if (meta.activity_index != activity_index) continue;
      for (const auto& [vid, subset] : world.getActivityVenues(meta)) {
        PersonLocation loc;
        loc.person_id = p.id;
        loc.venue_id = vid;
        loc.subset_index = subset;
        loc.activity_index = activity_index;
        loc.encounter_type_id = 255;
        loc.person_array_index = i;
        locs.push_back(loc);
      }
    }
  }
  return locs;
}

// Run the same contact-matrix finalize chain Simulator runs at startup, for
// tests that build an InteractionManager directly instead of going through
// Simulator. Transmission reads a table resolved here, so a test that skips
// this gets an unresolved-lookup throw rather than silently different physics.
// `mode_names` defaults to a single unnamed mode, matching a disease that
// declares none.
inline void finalizeContactMatrices(
    ContactMatrixConfig& cm, const WorldState& world,
    const std::vector<std::string>& mode_names = {}) {
  // A test that declares no contact structure at all is testing something
  // other than direct contact (fomite deposition, compartmental uptake, the
  // partial-presence gate). Give it a matrix that transmits nothing, so
  // direct contact is off unless a test asks for it: a test that meant to
  // exercise contact transmission and forgot to declare a matrix then sees
  // zero infections and fails, rather than picking up a rate from nowhere.
  if (!cm.default_matrix.has_value() && cm.matrices.empty() &&
      cm.mode_matrices.empty()) {
    ContactMatrix silent;
    silent.bins = {"all"};
    silent.contacts = {{0.0}};
    cm.default_matrix = silent;
    cm.allow_default_matrix = true;
  }
  cm.resolve(world);
  cm.finalizeDefaultModeMatrices(world, mode_names);
  cm.finalizeDiseaseModeAlignment(mode_names);
  cm.finalizeResolvedMatrices(world, mode_names);
}

// Same, taking the mode list from the disease that will drive transmission.
// The resolved table is indexed by disease mode, so a multi-mode disease needs
// its modes named here or the force-of-infection loop asks for a mode the
// table does not have.
inline void finalizeContactMatrices(ContactMatrixConfig& cm,
                                    const WorldState& world,
                                    const Disease& disease) {
  finalizeContactMatrices(cm, world, disease.getModeNames());
}

}  // namespace june
