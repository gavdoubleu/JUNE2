#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "core/world_state.h"
#include "doctest.h"
#include "loaders/domain_loader_internals.h"

using namespace june;

// The seam behind buildGlobalVenueMaps(). A rank must be able to name the type
// of every Venue in the world, not just the ones decomposed onto it — otherwise
// kUnknownVenueTypeId means both "no such Venue" and "not mine", and
// venue-gated policy becomes rank-dependent.
TEST_CASE("fillGlobalVenueMaps types venues this rank does not own") {
  WorldState world;
  world.venue_type_names = {"household", "school"};

  GeographicalUnit geo_unit;
  geo_unit.id = 10;
  geo_unit.parent_id = -1;
  geo_unit.level_id = 0;
  world.geo_units.push_back(geo_unit);

  // One local venue; venue 200 belongs to another rank. No activity_venues, so
  // no local person can reach 200 — it is outside any halo.
  Venue local_venue;
  local_venue.id = 100;
  local_venue.type_id = 0;
  local_venue.geo_unit_id = 10;
  world.venues.push_back(local_venue);
  world.buildIndices();

  detail::fillGlobalVenueMaps(world, {100, 200}, {0, 1}, {10, 10});

  SUBCASE("both venues land in the type map") {
    REQUIRE(world.global_venue_type_map.count(100) == 1);
    REQUIRE(world.global_venue_type_map.count(200) == 1);
    CHECK(world.global_venue_type_map.at(100) == 0);
    CHECK(world.global_venue_type_map.at(200) == 1);
  }

  SUBCASE("a foreign venue types via getVenueTypeId") {
    CHECK(world.getVenueTypeId(200) == 1);
  }

  SUBCASE("an id naming no Venue is still unresolvable") {
    CHECK(world.getVenueTypeId(300) == kUnknownVenueTypeId);
  }
}
