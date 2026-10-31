#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "core/world_state.h"

// Hand-built worlds for the vaccination uptake tests: a list of homes, each a
// residence venue with its people. A home is its own geo unit, so a world can
// be cut into rank slices by home and every slice holds its homes whole.
namespace june::uptake_fixture {

constexpr int16_t kResidence = 0;
constexpr uint8_t kPlainType = 0;
constexpr uint8_t kCareType = 1;

struct Home {
  std::vector<float> ages;
  bool care = false;  // a care home (venue type care_home) rather than a home
  // The residents' geo unit, when it should differ from the venue's
  std::optional<GeoUnitId> person_geo_unit;
  bool registered = true;    // false: nobody has a residence registration
  bool venue_loaded = true;  // false: the venue is missing from the world
  std::vector<VenueId> extra_venues;  // more residence registrations
};

// Person ids run on in home order over every home, kept or not, so a person
// has the same id in a slice as in the whole world. `keep` picks the homes this
// slice holds. Venue id and geo unit are the home's position in the list.
inline WorldState makeWorld(
    const std::vector<Home>& homes,
    const std::function<bool(size_t)>& keep = [](size_t) { return true; }) {
  WorldState world;
  world.activity_names = {"residence"};
  world.venue_type_names = {"home", "care_home"};
  PersonId next_id = 0;
  for (size_t h = 0; h < homes.size(); ++h) {
    const Home& home = homes[h];
    const bool mine = keep(h);
    const auto venue_id = static_cast<VenueId>(h);
    // A venue this slice does not hold still has a known type, as in a real
    // distributed load
    world.setGlobalVenueType(venue_id, home.care ? kCareType : kPlainType);
    if (mine && home.venue_loaded) {
      Venue venue{};
      venue.id = venue_id;
      venue.type_id = home.care ? kCareType : kPlainType;
      venue.geo_unit_id = static_cast<GeoUnitId>(h);
      venue.is_residence = true;
      world.venues.push_back(venue);
    }
    for (float age : home.ages) {
      const PersonId id = next_id++;
      if (!mine) continue;
      Person& person = world.people.emplace_back();
      person.id = id;
      person.age = age;
      person.sex = Sex::UNKNOWN;
      person.geo_unit_id =
          home.person_geo_unit.value_or(static_cast<GeoUnitId>(h));
      person.activity_meta_start =
          static_cast<uint32_t>(world.activity_meta.size());
      if (!home.registered) continue;
      std::vector<VenueId> venues = {venue_id};
      venues.insert(venues.end(), home.extra_venues.begin(),
                    home.extra_venues.end());
      world.activity_meta.push_back(
          {kResidence, static_cast<uint32_t>(world.activity_venues.size()),
           static_cast<uint16_t>(venues.size())});
      for (VenueId v : venues) world.activity_venues.push_back({v, 0});
      person.activity_meta_count = 1;
    }
  }
  world.buildIndices();
  return world;
}

// Every seventh home is a care home of three old residents. The others have two
// adults and none, one or two children. Adults get the lower ids.
inline std::vector<Home> mixedHomes(size_t count) {
  std::vector<Home> homes(count);
  for (size_t h = 0; h < count; ++h) {
    if (h % 7 == 6) {
      homes[h].care = true;
      homes[h].ages = {82, 85, 90};
      continue;
    }
    homes[h].ages = {static_cast<float>(30 + h % 5),
                     static_cast<float>(40 + h % 7)};
    for (size_t k = 0; k < h % 3; ++k)
      homes[h].ages.push_back(static_cast<float>(4 + k + h % 8));
  }
  return homes;
}

// Ordinary homes of two adults and a child with a bad home at position `before`
// and another at the end, so the lowest failing id is not the only failure.
inline std::vector<Home> homesWith(const Home& bad, size_t before,
                                   size_t after) {
  std::vector<Home> homes(before, Home{{30, 40, 5}});
  homes.push_back(bad);
  homes.insert(homes.end(), after, Home{{30, 40, 5}});
  homes.push_back(bad);
  return homes;
}

inline SelectionCriterion criterion(std::string path, std::string op,
                                    PropertyValue value) {
  SelectionCriterion c;
  c.property_path = std::move(path);
  c.operator_type = std::move(op);
  c.value = std::move(value);
  return c;
}

inline VaccinationTier makeTier(std::string name,
                                std::vector<SelectionCriterion> group) {
  VaccinationTier tier;
  tier.name = std::move(name);
  tier.selection_criteria = std::move(group);
  return tier;
}

inline VaccinationTier children() {
  return makeTier("children", {criterion("age", "<", 18.0)});
}

inline VaccinationTier adults() {
  return makeTier("adults", {criterion("age", ">=", 18.0)});
}

// Everyone in a home decides together, adults decide, care homes are outside
// because only the venue type "home" is in scope.
inline void shareResidence(VaccinationTier& tier, DeciderCombine combine,
                           NoDecider none = NoDecider::Error) {
  tier.decision_unit.shared_activity = "residence";
  tier.decision_unit.activity_index = kResidence;
  tier.decision_rules.scope = {
      criterion("activities.residence.venue_type", "==", std::string("home"))};
  tier.decision_rules.deciders = {criterion("age", ">=", 18.0)};
  tier.decision_rules.combine = combine;
  tier.decision_rules.if_no_decider = none;
}

inline VaccinationCampaignConfig makeCampaign(
    std::vector<VaccinationTier> tiers) {
  VaccinationCampaignConfig campaign;
  campaign.name = "rollout";
  campaign.tiers = std::move(tiers);
  return campaign;
}

inline UptakeRow row(std::vector<SelectionCriterion> criteria,
                     double probability) {
  return {std::move(criteria), probability};
}

}  // namespace june::uptake_fixture
