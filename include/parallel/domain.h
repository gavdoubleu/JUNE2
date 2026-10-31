#pragma once

#ifdef USE_MPI

#include <mpi.h>

#include <unordered_set>
#include <utility>
#include <vector>

#include "core/types.h"
#include "core/world_state.h"
#include "epidemiology/emission/emission_record.h"

namespace june {

/**
 * Domain - Represents a geographic partition of the world for parallel
 * execution
 *
 * Each MPI rank owns one Domain, which contains:
 * - A subset of geographic units (e.g., MGUs)
 * - All people who RESIDE in those geo units
 * - All venues LOCATED in those geo units
 * - Tracking of visitors from other domains during simulation
 */
class Domain {
 public:
  // MPI identifiers
  int rank;       // This process's MPI rank
  int num_ranks;  // Total number of MPI processes

  // Geographic coverage
  std::vector<GeoUnitId> geo_unit_ids;  // Geo units in this domain (e.g., MGUs)
  std::unordered_set<GeoUnitId> geo_unit_set;

  // Population partitioning
  std::vector<PersonId> resident_ids;  // People who live in this domain
  std::unordered_set<PersonId> resident_set;

  // Venue partitioning
  std::vector<VenueId> local_venue_ids;  // Venues in this domain
  std::unordered_set<VenueId> local_venue_set;

  // Cross-domain visitor data (minimal transfer for performance)
  struct VisitorData {
    PersonId person_id;
    int home_rank;           // Which rank owns this person
    VenueId venue_id;        // Which venue they're visiting
    SubsetIndex subset_idx;  // Which subset within the venue

    // Infection state (for transmission calculations)
    bool is_infected;
    float immunity_level;
    uint8_t encounter_type_id;  // Coordinated encounter type (ID in registry)
    float age = 0.0f;           // Contact-matrix binning: the visitor is binned
    uint8_t sex_code = 2;       // as its home rank would bin the same person

    // Stage at packing time: infector attribution and compartmental
    // deposition writeback
    uint16_t symptom_id = 0;     // Current symptom ID at packing time
    double time_in_stage = 0.0;  // Time in current stage at packing time

    // Computed on the sending rank by EmissionCalculator::emit, the same
    // code path as locals (bit-identical FP), then scaled by the source
    // modifier as locals are. See parallel/visitor_wire.h.
    Emission emission;

    // Effective target susceptibility per transmission mode, computed on the
    // home rank. Visitor modifier-set IDs are rank-local, so the values—not
    // the compact IDs—cross MPI.
    std::vector<double> target_susceptibility;

    // Source multipliers for deposition modes. Compartmental deposition is
    // reconstructed on the venue rank, so its source effect must be sent
    // separately from the emission. Fomite entries are unused:
    // emission.fomite_deposits already carries the source modifier.
    std::vector<double> deposition_source_multiplier;

    // Return data: infection status changes
    bool newly_infected;
    double new_infection_time;
  };

  std::vector<VisitorData>
      incoming_visitors;  // Visitors at our venues (from other ranks)

  // Domain-local state
  WorldState* world;  // Pointer to shared read-only world state

  // Statistics
  size_t local_population;  // Number of residents
  size_t local_venues;      // Number of venues

  // Constructor
  Domain(int rank, int num_ranks, WorldState* world);

  // Query methods
  bool ownsGeoUnit(GeoUnitId id) const { return geo_unit_set.count(id) > 0; }

  bool ownsPerson(PersonId id) const { return resident_set.count(id) > 0; }

  // Dynamic virtual venue ownership set (populated by encounter injection)
  std::unordered_set<VenueId> virtual_venue_set;

  void registerVirtualVenue(VenueId id) { virtual_venue_set.insert(id); }
  void clearVirtualVenues() { virtual_venue_set.clear(); }

  bool ownsVenue(VenueId id) const {
    if (isVirtualVenue(id)) {
      // Check dynamic registry (set during encounter injection)
      return virtual_venue_set.count(id) > 0;
    }
    return local_venue_set.count(id) > 0;
  }

  // Add a geo unit to this domain
  void addGeoUnit(GeoUnitId id) {
    geo_unit_ids.push_back(id);
    geo_unit_set.insert(id);
  }

  // Assign people and venues based on geo unit ownership. Throws if the world
  // breaks the geo unit locality of residence: see checkResidencesAreLocal.
  void assignPeopleAndVenues();

  // Visitor management
  void clearVisitors() { incoming_visitors.clear(); }

  void addIncomingVisitor(const VisitorData& visitor) {
    incoming_visitors.push_back(visitor);
  }

  void addIncomingVisitor(VisitorData&& visitor) {
    incoming_visitors.push_back(std::move(visitor));
  }

  // Print domain statistics
  void printStatistics() const;

 private:
  // A residence venue and its occupants share a geo unit, so a household is
  // held whole by one rank. The partition assumes it, and clustered seeding's
  // local density truncation is exact only under it. Throws
  // naming the offenders rather than running on a world that breaks it.
  void checkResidencesAreLocal() const;
};

}  // namespace june

#endif  // USE_MPI
