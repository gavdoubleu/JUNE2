// Pipeline for one simulation timeslot: exchange visitors, process
// transmission, apply inbound results, then update epidemiology.
// Declarations are in simulation/simulator.h.
#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

#include "simulation/simulator.h"

namespace june {

namespace {

// Reduce per-rank totals to rank 0, then print the transmission count and any
// non-zero epidemiology counts for this slot.
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

// Count locations by activity index, then reduce the counts to rank 0.
// world.activity_names gives every rank the same buffer size. Reducing by
// name with std::map would produce different buffer sizes when ranks hold
// different activities and can trigger MPI_ERR_TRUNCATE.
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
    // Only rank 0 prints, so Reduce is enough; other ranks do not need the
    // aggregated result.
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
  // Apply every transition up to the end of this slot, so someone who dies in
  // it is already dead when the next slot places people. The clock adds slot
  // lengths one by one and can stop a hair short of midnight; the day's last
  // slot ends exactly on it.
  double slot_end = current_simulation_time_ + delta_hours / 24.0;
  const double midnight = current_day_num_ + 1.0;
  if (midnight - slot_end < 1e-9) slot_end = midnight;
  try {
    epi_stats = epidemiology_->updateInfectionStates(slot_end, locations_);
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
    if (interaction_manager_) {
      interaction_manager_->setCurrentDayTypeIdx(day_type_idx);
    }
    local_new_infections = interaction_manager_->processTransmissions(
        transmission_locations, current_simulation_time_, delta_hours,
        &epidemiology_->getActiveInfectionsMutable(), visitor_ids,
        pending_infections, visitor_data_map,
        compartmental_model_manager_.get());

    // Process transport lines separately from venues because a rider on
    // multiple legs is not represented by one location entry. Each rank
    // processes its owned lines, then the ranks resolve infections together
    // so a rider infected on two legs in the same slot gets one infection.
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
        domain_mgr_->receivePendingInfections(pending_infections, *disease_);
    for (const auto& applied : mpi_infected) {
      const TransmissionRecord transmission{applied.source,
                                            applied.infector_symptom_id,
                                            applied.transmission_mode_index};
      epidemiology_->trackInfection(applied.person_id);
      if (interaction_manager_)
        interaction_manager_->countInfectorLookupGap(transmission);
      event_logger_.logInfection(applied.person_id, applied.infector_id,
                                 applied.venue_id, applied.infection_time,
                                 applied.encounter_type_id, transmission);
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
    domain_mgr_->exchangeVisitors(locations_, *disease_,
                                  current_simulation_time_, delta_hours,
                                  runtime_group_allocator_.get());

    Domain& domain = domain_mgr_->getDomain();

    // Keep unallocated locations and people at venues owned by this rank.
    // The owning rank processes people at remote venues as visitors.
    augmented_locations.reserve(locations_.size() +
                                domain.incoming_visitors.size());
    for (const auto& loc : locations_) {
      if (loc.venue_id == -1) {
        augmented_locations.push_back(loc);  // keep unallocated location
        continue;
      }
      if (domain.ownsVenue(loc.venue_id)) {
        augmented_locations.push_back(loc);  // keep local location
      }
      // Skip remote venues; the owning rank handles them as visitors.
    }

    // Add incoming visitors to the augmented locations.
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

    // Sort incoming visitors so processing order does not depend on MPI
    // message arrival order. A rider on a multi-leg journey arrives once per
    // leg under the same person_id, so person_id alone ties; the venue and
    // subset finish the key.
    std::stable_sort(augmented_locations.begin() + visitor_start,
                     augmented_locations.end(),
                     [](const PersonLocation& a, const PersonLocation& b) {
                       if (a.person_id != b.person_id)
                         return a.person_id < b.person_id;
                       if (a.venue_id != b.venue_id)
                         return a.venue_id < b.venue_id;
                       if (a.subset_index != b.subset_index)
                         return a.subset_index < b.subset_index;
                       return a.encounter_type_id < b.encounter_type_id;
                     });

    // Collect visitor IDs for InteractionManager.
    visitor_ids = domain_mgr_->getVisitorIds();

    // Build the visitor data used for transmission calculations. Move the
    // emission and susceptibility vectors into the map instead of copying
    // them.
    for (auto& visitor : domain.incoming_visitors) {
      VisitorInfo info;
      info.person_id = visitor.person_id;
      info.is_infected = visitor.is_infected;
      info.immunity_level = visitor.immunity_level;
      info.age = visitor.age;
      info.sex = static_cast<Sex>(visitor.sex_code);
      info.symptom_id = visitor.symptom_id;
      info.time_in_stage = visitor.time_in_stage;
      info.emission = std::move(visitor.emission);
      info.target_susceptibility = std::move(visitor.target_susceptibility);
      // The wire format stores one multiplier per deposition mode. Expand it
      // to an entry for each transmission mode.
      if (!visitor.deposition_source_multiplier.empty()) {
        const auto& modes = disease_->getTransmissionParams().modes;
        info.deposition_source_multiplier.assign(modes.size(), 1.0);
        size_t deposition_index = 0;
        for (size_t mode = 0; mode < modes.size(); ++mode) {
          if (modes[mode].type != TransmissionModeType::Fomite &&
              modes[mode].type != TransmissionModeType::CompartmentalDeposition)
            continue;
          if (deposition_index >= visitor.deposition_source_multiplier.size())
            break;
          info.deposition_source_multiplier[mode] =
              visitor.deposition_source_multiplier[deposition_index++];
        }
        if (deposition_index == 0) info.deposition_source_multiplier.clear();
      }
      visitor_data_map[visitor.person_id] = std::move(info);
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

  // Compartmental coupling sequence. Each step consumes the state produced by
  // the preceding one:
  // 1. advance(): integrates the ODE using the previous slot's inputs.
  // 2. processTransmissions(): exposes people to the plugin's FOI, reading the
  //    buffer lazily.
  // 3. computeDepositionWriteback(): aggregates infections into plugin inputs.
  // 4. maybeSnapshot(): records plugin state after the slot completes.
  compartmental_model_manager_->advance(
      static_cast<float>(delta_hours / 24.0),
      static_cast<float>(current_simulation_time_));

  // Step 0: Apply the infection seeds dated since the previous slot start
  applyInfectionSeeds(seedWindowForSlot(config_.schedule,
                                        config_.simulation.start_date,
                                        current_day_num_, time_slot_index));

  // Step 1: Update the time used by policy checks, then assign activities from
  // the precomputed schedule.
  {
    activity_manager_.setCurrentTime(current_simulation_time_);
    activity_manager_.assignActivitiesFromSchedule(time_slot_index,
                                                   day_type_idx, locations_);
    // Allocate runtime groups for partial-presence venues such as train
    // groups. This call is a no-op when SimulationConfig::partial_presence is
    // empty, so other scenarios do not do this work.
    runtime_group_allocator_->allocateForSlot(time_slot_index, day_type_idx,
                                              slot, current_simulation_time_,
                                              delta_hours, locations_);
  }

#ifdef USE_MPI
  // Clear the per-slot virtual venue registry before injecting encounters.
  // Remove stale rank assignments too, so a hash collision with a venue ID
  // from an earlier slot does not block re-registration.
  if (domain_mgr_) {
    domain_mgr_->getDomain().clearVirtualVenues();
    domain_mgr_->clearVirtualVenueRanks();
  }
#endif

  // Inject Coordinated Encounters for this timeslot into the locations array.
  injectCoordinatedEncountersIntoSlot(time_slot_index);

  // Place followers at their host's resolved location for this slot.
  injectFollowsIntoSlot(time_slot_index);

  // A person on a partial-presence venue must have at least one rider leg, and
  // every rider must be placed on a partial-presence venue. A mismatch would
  // either omit the person from transmission or keep them in a group after
  // they leave the line, so fail before processing the slot.
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

  // The collective Reduce gathers the per-slot venue distribution before rank
  // 0 prints it.
  if (policy_manager_) {
    policy_manager_->refreshTransmissionModifiers(
        current_simulation_time_, epidemiology_->getActiveInfectionsMutable());
  }
  printSlotVenueDistribution(world_, locations_, domain_mgr_, rank);

#ifdef USE_MPI
  // Step 2: Exchange visitors between domains in MPI mode.
  std::vector<PersonLocation> augmented_locations;
  std::unordered_set<PersonId> visitor_ids;
  std::vector<PendingInfection> pending_infections;
  std::unordered_map<PersonId, VisitorInfo> visitor_data_map;
  exchangeVisitorsAndBuildAugmented(delta_hours, augmented_locations,
                                    visitor_ids, visitor_data_map);

  // In MPI mode, use local locations plus incoming visitors. Otherwise use the
  // original locations.
  std::vector<PersonLocation>& transmission_locations =
      domain_mgr_ ? augmented_locations : locations_;
#else
  // Serial mode uses the original locations.
  std::vector<PersonLocation>& transmission_locations = locations_;
#endif

  // Step 3: Calculate contacts and transmission using active infections.
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
  // Step 4: Return pending infections to their home ranks in MPI mode.
  receivePendingAndApply(pending_infections);
#endif

  // Aggregate infected people's per-node deposition at owned venues and pass
  // it to the plugin before the next advance() call. This reads the infections
  // as they stood during the slot, so it runs before the state update clears
  // anyone who recovers or dies in it.
  const std::unordered_map<PersonId, VisitorInfo>* deposition_visitor_data =
      nullptr;
#ifdef USE_MPI
  if (have_mpi) deposition_visitor_data = &visitor_data_map;
#endif
  compartmental_model_manager_->computeDepositionWriteback(
      transmission_locations, world_, *disease_, current_simulation_time_,
      current_simulation_time_ + delta_hours / 24.0, policy_manager_.get(),
      deposition_visitor_data);

  // Steps 5 and 6: update infection states and decay venue fomites after
  // transmission. This ensures newly infected people are tracked and deaths
  // are processed at the end of the slot.
  EpiSlotStats epi_stats = updateEpidemiologyAfterTransmission(delta_hours);

  // Reduce and print the per-slot transmission and epidemiology summary.
  printSlotEpiSummary(local_new_infections, epi_stats, delta_hours, domain_mgr_,
                      rank);

  compartmental_model_manager_->maybeSnapshot(
      static_cast<float>(current_simulation_time_));
}

}  // namespace june
