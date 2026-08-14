// Simulator per-timeslot pipeline: visitor exchange, transmission, applying
// inbound results, post-transmission epidemiology update. Split from
// simulator.cpp (declared in simulation/simulator.h).
#include <algorithm>
#include <iostream>
#include <vector>

#include "simulation/simulator.h"

namespace june {

namespace {

// Per-slot transmission + epidemiology summary print: MPI_Reduce the
// per-rank totals onto rank 0 and print one line of transmissions + one
// line of epi transitions/recoveries/deaths if any are non-zero.
void printSlotEpiSummary(int local_new_infections,
                         const EpiSlotStats& epi_stats, double delta_hours,
                         DomainManager* domain_mgr, int rank) {
  int global_new_infections = local_new_infections;
  int local_epi[4] = {epi_stats.transitions, epi_stats.recoveries,
                      epi_stats.deaths, epi_stats.active_remaining};
  int global_epi[4];
#ifdef USE_MPI
  if (domain_mgr) {
    MPI_Reduce(&local_new_infections, &global_new_infections, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_epi, global_epi, 4, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    std::copy(std::begin(local_epi), std::end(local_epi),
              std::begin(global_epi));
  }
#else
  (void)domain_mgr;
  std::copy(std::begin(local_epi), std::end(local_epi), std::begin(global_epi));
#endif
  if (rank != 0) return;
  if (global_new_infections > 0) {
    std::cout << "      [Transmission] " << global_new_infections
              << " new infections (duration=" << delta_hours << "h)"
              << std::endl;
  }
  if (global_epi[0] > 0 || global_epi[1] > 0 || global_epi[2] > 0) {
    std::cout << "      [Epidemiology] Processed " << global_epi[0]
              << " symptom transitions. " << global_epi[1] << " recoveries, "
              << global_epi[2] << " deaths. "
              << "Active infections remaining: " << global_epi[3] << std::endl;
  }
}

// Per-slot venue distribution print: bin locations by activity index and
// MPI_Reduce onto rank 0 for printing. Indexed by activity_index over
// world.activity_names so the Reduce buffer size is identical on every
// rank. Collapsing by name via std::map would break with
// MPI_ERR_TRUNCATE as soon as ranks diverge in which activities they hold.
void printSlotVenueDistribution(const WorldState& world,
                                const std::vector<PersonLocation>& locations,
                                DomainManager* domain_mgr, int rank) {
  const size_t num_activities = world.activity_names.size();
  std::vector<int> local_counts(num_activities, 0);
  for (const auto& loc : locations) {
    if (loc.activity_index >= 0 &&
        loc.activity_index < static_cast<int>(num_activities)) {
      local_counts[loc.activity_index]++;
    }
  }
  std::vector<int> global_counts(num_activities, 0);
#ifdef USE_MPI
  if (domain_mgr) {
    // Rank-0 print only: Reduce instead of Allreduce; non-root ranks
    // don't need the aggregated result.
    MPI_Reduce(local_counts.data(), global_counts.data(),
               static_cast<int>(num_activities), MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);
  } else {
    global_counts = local_counts;
  }
#else
  (void)domain_mgr;
  global_counts = local_counts;
#endif
  if (rank == 0) {
    std::cout << "      → ";
    for (size_t i = 0; i < num_activities; ++i) {
      if (global_counts[i] > 0) {
        std::cout << world.activity_names[i] << ": " << global_counts[i]
                  << "  ";
      }
    }
    std::cout << std::endl;
  }
}

}  // namespace

EpiSlotStats Simulator::updateEpidemiologyAfterTransmission(
    double delta_hours) {
  EpiSlotStats epi_stats;
  try {
    ScopedTimer timer("04_InfectionStateUpdates");
    epi_stats = epidemiology_->updateInfectionStates(current_simulation_time_,
                                                     locations_);
  } catch (const std::exception& e) {
    std::cerr << "[Step 5 Infection Updates] Fatal error: " << e.what()
              << std::endl;
    throw;
  }
  try {
    epidemiology_->updateVenueFomites(current_simulation_time_, delta_hours);
  } catch (const std::exception& e) {
    std::cerr << "[Step 6 Fomites] Fatal error: " << e.what() << std::endl;
    throw;
  }
  return epi_stats;
}

int Simulator::runSlotTransmission(
    std::vector<PersonLocation>& transmission_locations, double delta_hours,
    int day_type_idx, std::unordered_set<PersonId>* visitor_ids,
    std::vector<PendingInfection>* pending_infections,
    std::unordered_map<PersonId, VisitorInfo>* visitor_data_map) {
  int local_new_infections = 0;
  try {
    ScopedTimer timer("03_TransmissionProcessing");
    if (interaction_manager_) {
      interaction_manager_->setCurrentDayTypeIdx(day_type_idx);
    }
    local_new_infections = interaction_manager_->processTransmissions(
        transmission_locations, current_simulation_time_, delta_hours,
        &epidemiology_->getActiveInfectionsMutable(), visitor_ids,
        pending_infections, visitor_data_map,
        compartmental_model_manager_.get());

    // Transport lines are handled apart from the venue loop above, because a
    // rider on several legs cannot be found from the location table. Every
    // rank walks the lines it owns, then all of them settle the results
    // together so a rider infected on two legs at once gets one infection.
    if (runtime_group_allocator_ && runtime_group_allocator_->isActive()) {
      std::vector<VenueId> owned_lines;
      owned_lines.reserve(runtime_group_allocator_->ridersByVenue().size());
      for (const auto& [vid, riders] :
           runtime_group_allocator_->ridersByVenue()) {
#ifdef USE_MPI
        if (domain_mgr_ && !domain_mgr_->getDomain().ownsVenue(vid)) continue;
#endif
        owned_lines.push_back(vid);
      }
      std::sort(owned_lines.begin(), owned_lines.end());

      interaction_manager_->processPartialPresenceLines(
          owned_lines, current_simulation_time_, delta_hours,
          &epidemiology_->getActiveInfectionsMutable(), visitor_data_map);
      local_new_infections +=
          interaction_manager_->resolvePartialPresenceInfections(
              current_simulation_time_,
              &epidemiology_->getActiveInfectionsMutable());
    }
  } catch (const std::exception& e) {
    std::cerr << "[Step 3 Transmission] Fatal error: " << e.what() << std::endl;
    throw;
  }
  return local_new_infections;
}

#ifdef USE_MPI
void Simulator::receivePendingAndApply(
    const std::vector<PendingInfection>& pending_infections) {
  if (domain_mgr_ == nullptr) return;
  try {
    auto mpi_infected =
        domain_mgr_->receivePendingInfections(pending_infections);
    for (const auto& applied : mpi_infected) {
      epidemiology_->trackInfection(applied.person_id);
      event_logger_.logInfection(
          applied.person_id, applied.infector_id, applied.venue_id,
          applied.infection_time, applied.encounter_type_id,
          applied.infector_symptom_id, applied.transmission_mode_index);
    }
  } catch (const std::exception& e) {
    std::cerr << "[Step 4 Receive Pending] Fatal error: " << e.what()
              << std::endl;
    throw;
  }
}

void Simulator::exchangeVisitorsAndBuildAugmented(
    double delta_hours, std::vector<PersonLocation>& augmented_locations,
    std::unordered_set<PersonId>& visitor_ids,
    std::unordered_map<PersonId, VisitorInfo>& visitor_data_map) {
  if (domain_mgr_ == nullptr) {
    augmented_locations = locations_;
    return;
  }
  try {
    ScopedTimer timer("02_MPI_VisitorExchange");
    domain_mgr_->exchangeVisitors(locations_, current_simulation_time_,
                                  delta_hours, runtime_group_allocator_.get());

    Domain& domain = domain_mgr_->getDomain();

    // Filter out outgoing visitors (local people at remote venues): we
    // only keep local people at LOCAL venues. People at remote venues are
    // handled exclusively as visitors on the owning rank.
    augmented_locations.reserve(locations_.size() +
                                domain.incoming_visitors.size());
    for (const auto& loc : locations_) {
      if (loc.venue_id == -1) {
        augmented_locations.push_back(loc);  // unallocated, keep
        continue;
      }
      if (domain.ownsVenue(loc.venue_id)) {
        augmented_locations.push_back(loc);  // local venue, keep
      }
      // remote venue: skip (handled as visitor on owning rank)
    }

    // Add incoming visitors to augmented locations
    size_t visitor_start = augmented_locations.size();
    for (const auto& visitor : domain.incoming_visitors) {
      PersonLocation visitor_loc;
      visitor_loc.person_id = visitor.person_id;
      visitor_loc.venue_id = visitor.venue_id;
      visitor_loc.subset_index = visitor.subset_idx;
      visitor_loc.activity_index =
          static_cast<int16_t>(world_.getActivityIndex("visiting"));
      visitor_loc.encounter_type_id = visitor.encounter_type_id;
      augmented_locations.push_back(visitor_loc);
    }

    // Sort visitor portion by person_id for deterministic processing order
    // (MPI message arrival order is non-deterministic)
    std::sort(augmented_locations.begin() + visitor_start,
              augmented_locations.end(),
              [](const PersonLocation& a, const PersonLocation& b) {
                return a.person_id < b.person_id;
              });

    // Get visitor IDs for InteractionManager
    visitor_ids = domain_mgr_->getVisitorIds();

    // Populate visitor data map for transmission calculations
    for (const auto& visitor : domain.incoming_visitors) {
      VisitorInfo info;
      info.person_id = visitor.person_id;
      info.is_infected = visitor.is_infected;
      info.is_infectious = visitor.is_infectious;
      info.immunity_level = visitor.immunity_level;
      info.symptom_id = visitor.symptom_id;
      info.time_in_stage = visitor.time_in_stage;
      std::copy(std::begin(visitor.integrated_infectiousness),
                std::end(visitor.integrated_infectiousness),
                std::begin(info.integrated_infectiousness));
      visitor_data_map[visitor.person_id] = info;
    }
  } catch (const std::exception& e) {
    std::cerr << "[Step 2 MPI] Fatal error: " << e.what() << std::endl;
    throw;
  }
}
#endif

void Simulator::simulateTimeSlot(const TimeSlot& slot, int time_slot_index,
                                 int day_type_idx, double delta_hours) {
  const int rank = getRank();

  printSimulationState(slot.name, delta_hours);

  // Compartmental coupling sequence (ORDER IS LOAD-BEARING):
  // 1. advance():                    plugin integrates ODE with previous
  //                                  slot's inputs
  // 2. processTransmissions():       humans exposed to FOI from plugin (reads
  //                                  buffer lazily)
  // 3. computeDepositionWriteback(): aggregate infections into plugin inputs
  // 4. maybeSnapshot():              record plugin state after full slot
  compartmental_model_manager_->advance(
      static_cast<float>(delta_hours / 24.0),
      static_cast<float>(current_simulation_time_));

  // Step 0: Apply infection seeds for this time
  std::string current_dt = formatDate(current_date_) + " " + slot.start;
  applyInfectionSeeds(current_dt);

  // Step 1: Assign people to activities using pre-computed schedules
  // Update current time for policy checks
  {
    ScopedTimer timer("01_ActivityAssignment");
    activity_manager_.setCurrentTime(current_simulation_time_);
    activity_manager_.assignActivitiesFromSchedule(time_slot_index,
                                                   day_type_idx, locations_);
    // Runtime group allocation for partial-presence venues (e.g. train
    // groups). One-test no-op when SimulationConfig::partial_presence is
    // empty, so non-commute scenarios pay nothing here.
    runtime_group_allocator_->allocateForSlot(time_slot_index, day_type_idx,
                                              slot, current_simulation_time_,
                                              delta_hours, locations_);
  }

#ifdef USE_MPI
  // Clear per-slot virtual venue registry before encounter injection.
  // Also clear stale virtual venue rank assignments so that hash-colliding
  // venue IDs from previous slots don't prevent correct re-registration.
  if (domain_mgr_) {
    domain_mgr_->getDomain().clearVirtualVenues();
    domain_mgr_->clearVirtualVenueRanks();
  }
#endif

  // Inject Coordinated Encounters for this timeslot into the locations array.
  injectCoordinatedEncountersIntoSlot(time_slot_index);

  // Place followers at their host's resolved location for this slot.
  injectFollowsIntoSlot(time_slot_index);

  // Everyone on a line must be a rider of it, and every rider must be on one.
  // Riding without a rider entry was the old ghost: aboard, but infecting
  // nobody and catching nothing. Holding an entry while standing somewhere
  // else is its mirror image, and would crowd a group with someone who
  // left. Both are bugs, so say so rather than quietly modelling a phantom.
  if (runtime_group_allocator_ && runtime_group_allocator_->isActive()) {
    for (const auto& loc : locations_) {
      if (loc.person_id < 0) continue;
      const bool on_line =
          runtime_group_allocator_->isPartialPresenceVenue(loc.venue_id);
      const bool rides =
          !runtime_group_allocator_->legsOf(loc.person_id).empty();
      if (on_line == rides) continue;

      const std::string who = "person " + std::to_string(loc.person_id);
      if (on_line)
        throw std::runtime_error(
            who + " was placed on partial-presence venue " +
            std::to_string(loc.venue_id) +
            " but rides no leg of it. Something put them on a line after the "
            "groups were dealt without telling the allocator.");
      throw std::runtime_error(
          who + " rides a partial-presence venue but was placed at venue " +
          std::to_string(loc.venue_id) +
          ". Something moved them off their commute after the groups were "
          "dealt, leaving them aboard a line they are no longer on.");
    }
  }

  // Per-slot venue distribution print (collective Reduce → rank 0 prints).
  printSlotVenueDistribution(world_, locations_, domain_mgr_, rank);

#ifdef USE_MPI
  // Step 2: Exchange visitors between domains (parallel mode only)
  std::vector<PersonLocation> augmented_locations;
  std::unordered_set<PersonId> visitor_ids;
  std::vector<PendingInfection> pending_infections;
  std::unordered_map<PersonId, VisitorInfo> visitor_data_map;
  exchangeVisitorsAndBuildAugmented(delta_hours, augmented_locations,
                                    visitor_ids, visitor_data_map);

  // Use augmented locations (locals + visitors) for transmission processing
  std::vector<PersonLocation>& transmission_locations =
      domain_mgr_ ? augmented_locations : locations_;
#else
  // Serial mode: use original locations
  std::vector<PersonLocation>& transmission_locations = locations_;
#endif

  // Step 3: Calculate contacts and transmission (pass active_infections)
  int local_new_infections;
#ifdef USE_MPI
  const bool have_mpi = (domain_mgr_ != nullptr);
  local_new_infections =
      runSlotTransmission(transmission_locations, delta_hours, day_type_idx,
                          have_mpi ? &visitor_ids : nullptr,
                          have_mpi ? &pending_infections : nullptr,
                          have_mpi ? &visitor_data_map : nullptr);
#else
  local_new_infections =
      runSlotTransmission(transmission_locations, delta_hours, day_type_idx,
                          nullptr, nullptr, nullptr);
#endif

#ifdef USE_MPI
  // Step 4: Send back pending infections to home ranks (parallel mode only)
  receivePendingAndApply(pending_infections);
#endif

  // Steps 5 + 6: infection state update + venue fomite decay (must run
  // AFTER transmission so newly-infected people are tracked and death
  // processing lands at end of slot).
  EpiSlotStats epi_stats = updateEpidemiologyAfterTransmission(delta_hours);

  // Per-slot transmission + epidemiology summary (rank-0 prints after
  // Reduce).
  printSlotEpiSummary(local_new_infections, epi_stats, delta_hours, domain_mgr_,
                      rank);

  // Deposition write-back: aggregate per-node contributions from infected
  // people at owned venues and forward to the plugin for the next advance()
  // call.
  compartmental_model_manager_->computeDepositionWriteback(
      locations_, world_, *disease_,
      current_simulation_time_ - delta_hours / 24.0, current_simulation_time_);

  compartmental_model_manager_->maybeSnapshot(
      static_cast<float>(current_simulation_time_));
}

}  // namespace june
