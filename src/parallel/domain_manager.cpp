#ifdef USE_MPI

#include "parallel/domain_manager.h"

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>

#include "loaders/hdf5_loader.h"
#include "parallel/mpi_utils.h"

namespace june {

namespace {

// 2GB-safe chunked MPI_Bcast for large buffers.
void chunkedBroadcast(void* data, uint64_t count, size_t element_size,
                      MPI_Datatype type, int root) {
  const uint64_t MAX_MPI_INT =
      static_cast<uint64_t>(std::numeric_limits<int>::max());
  uint64_t sent = 0;
  while (sent < count) {
    uint64_t remaining = count - sent;
    int chunk = static_cast<int>(std::min(remaining, MAX_MPI_INT));
    MPI_Bcast(static_cast<char*>(data) + sent * element_size, chunk, type, root,
              MPI_COMM_WORLD);
    sent += chunk;
  }
}

// Pack a local registry of strings into a single null-separated buffer and
// MPI_Allgatherv-concatenate the same on every rank.
std::vector<char> gatherPackedValuesAcrossRanks(
    const std::vector<std::string>& local_registry, int num_ranks) {
  std::string local_packed;
  for (const auto& val : local_registry) {
    local_packed += val;
    local_packed += '\0';
  }

  int local_size = static_cast<int>(local_packed.size());
  std::vector<int> all_sizes(num_ranks);
  MPI_Allgather(&local_size, 1, MPI_INT, all_sizes.data(), 1, MPI_INT,
                MPI_COMM_WORLD);

  std::vector<int> displs;
  int total_size;
  mpi_utils::computeDisplacements(all_sizes, displs, total_size);

  std::vector<char> all_packed(total_size);
  MPI_Allgatherv(local_packed.data(), local_size, MPI_CHAR, all_packed.data(),
                 all_sizes.data(), displs.data(), MPI_CHAR, MPI_COMM_WORLD);
  return all_packed;
}

// Parse a null-separated packed buffer into a sorted, deduplicated registry.
std::vector<std::string> buildGlobalRegistryFromPacked(
    const std::vector<char>& all_packed) {
  std::set<std::string> unique_vals;
  size_t pos = 0;
  const size_t total_size = all_packed.size();
  while (pos < total_size) {
    size_t end = pos;
    while (end < total_size && all_packed[end] != '\0') ++end;
    if (end > pos)
      unique_vals.insert(std::string(all_packed.data() + pos, end - pos));
    pos = end + 1;
  }
  return std::vector<std::string>(unique_vals.begin(), unique_vals.end());
}

// Build a remap from a local registry's old indices to the corresponding
// indices in the global registry (-1 if not present).
std::vector<int32_t> buildRemapToGlobal(
    const std::vector<std::string>& local_registry,
    const std::vector<std::string>& global_registry) {
  std::unordered_map<std::string, int32_t> new_index;
  for (size_t i = 0; i < global_registry.size(); ++i) {
    new_index[global_registry[i]] = static_cast<int32_t>(i);
  }
  std::vector<int32_t> remap(local_registry.size(), -1);
  for (size_t i = 0; i < local_registry.size(); ++i) {
    auto it = new_index.find(local_registry[i]);
    if (it != new_index.end()) remap[i] = it->second;
  }
  return remap;
}

// Print per-rank and aggregate domain-load stats from already-gathered
// per-rank pop / venue / RSS vectors. Caller restricts to rank 0.
void printDomainStats(const std::string& label, int num_ranks,
                      const std::vector<uint64_t>& pops,
                      const std::vector<uint64_t>& venues,
                      const std::vector<double>& rss) {
  std::cout << "\n=== Domain Statistics [" << label << "] ===" << std::endl;
  uint64_t total_p = 0, max_p = 0;
  uint64_t total_v = 0;
  double total_rss = 0, max_rss = 0;

  for (int r = 0; r < num_ranks; ++r) {
    total_p += pops[r];
    max_p = std::max(max_p, pops[r]);
    total_v += venues[r];
    total_rss += rss[r];
    max_rss = std::max(max_rss, rss[r]);

    if (num_ranks <= 64 || r < 4 || r >= num_ranks - 4) {
      std::cout << "  Rank " << std::setw(3) << r << ": Pop=" << std::setw(10)
                << pops[r] << " | Venues=" << std::setw(8) << venues[r]
                << " | Mem=" << std::fixed << std::setprecision(2)
                << std::setw(6) << rss[r] << " GB" << std::endl;
    } else if (r == 4) {
      std::cout << "  ... (omitted) ..." << std::endl;
    }
  }

  double avg_p = static_cast<double>(total_p) / num_ranks;
  double avg_v = static_cast<double>(total_v) / num_ranks;
  double avg_rss = total_rss / num_ranks;

  std::cout << "  --------------------------------------------------"
            << std::endl;
  std::cout << "  TOTAL:    Pop=" << total_p << " | Venues=" << total_v
            << " | Mem=" << total_rss << " GB" << std::endl;
  std::cout << "  AVERAGE:  Pop=" << std::fixed << std::setprecision(0) << avg_p
            << " | Venues=" << avg_v << " | Mem=" << std::fixed
            << std::setprecision(2) << avg_rss << " GB" << std::endl;

  if (avg_p > 0)
    std::cout << "  Pop Imbalance: " << std::fixed << std::setprecision(2)
              << (max_p / avg_p - 1.0) * 100.0 << "%" << std::endl;
  if (avg_rss > 0)
    std::cout << "  Mem Imbalance: " << std::fixed << std::setprecision(2)
              << (max_rss / avg_rss - 1.0) * 100.0 << "%" << std::endl;
  std::cout << "====================================================\n"
            << std::endl;
}

// Apply a (old → new) index remap to every person's stored value for one
// property, in-place on the world's person_properties buffer.
void remapPersonPropertyValues(WorldState& world, const std::string& prop_name,
                               const std::vector<int32_t>& remap) {
  int prop_idx = world.getPersonPropertyIndex(prop_name);
  if (prop_idx < 0) return;
  for (auto& person : world.people) {
    size_t base = person.properties_start + prop_idx;
    if (base >= world.person_properties.size()) continue;
    int32_t old_val = world.person_properties[base];
    if (old_val >= 0 && old_val < static_cast<int32_t>(remap.size())) {
      world.person_properties[base] = remap[old_val];
    }
  }
}

}  // anonymous namespace

DomainManager::DomainManager(WorldState& world, const Config& config)
    : world_(world),
      config_(config),
      world_state_file_("world_state.h5"),
      rank_(0),
      num_ranks_(1),
      domain_(0, 1, &world, config.simulation.random_seed) {}

void DomainManager::setMPI(int rank, int num_ranks) {
  rank_ = rank;
  num_ranks_ = num_ranks;
  domain_.rank = rank;
  domain_.num_ranks = num_ranks;

  // Initialize required components for ownership lookups
  if (!partitioner_)
    partitioner_ = std::make_unique<GeographyPartitioner>(world_, config_);
  if (!communicator_)
    communicator_ =
        std::make_unique<DomainCommunicator>(world_, config_, domain_);
}

void DomainManager::setDisease(const Disease* disease) {
  if (communicator_) communicator_->setDisease(disease);
}

void DomainManager::loadGeographyOnNonZeroRanks() {
  WorldState temp = HDF5Loader::loadGeographyOnly(world_state_file_);
  world_.geo_units = std::move(temp.geo_units);
  // Copy registries so that findParentAtLevel() and buildVenueOwnershipMap()
  // can resolve geographic hierarchy levels during partitioning.
  world_.geo_level_names = std::move(temp.geo_level_names);
  world_.venue_type_names = std::move(temp.venue_type_names);
  world_.schedule_type_names = std::move(temp.schedule_type_names);
  world_.activity_names = std::move(temp.activity_names);
  world_.encounter_type_names = std::move(temp.encounter_type_names);
  world_.subset_type_names = std::move(temp.subset_type_names);
  world_.network_names = std::move(temp.network_names);
  world_.person_property_names = std::move(temp.person_property_names);
  world_.person_property_value_registries =
      std::move(temp.person_property_value_registries);
  world_.geo_unit_property_names = std::move(temp.geo_unit_property_names);
  world_.geo_unit_property_value_registries =
      std::move(temp.geo_unit_property_value_registries);
  world_.buildIndices();
}

void DomainManager::computeGlobalMaxPersonId() {
  PersonId local_max_p = -1;
  for (const auto& p : world_.people)
    if (p.id > local_max_p) local_max_p = p.id;
  MPI_Allreduce(&local_max_p, &max_person_id_, 1, MPI_INT32_T, MPI_MAX,
                MPI_COMM_WORLD);
}

void DomainManager::initialize() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks_);

  domain_.rank = rank_;
  domain_.num_ranks = num_ranks_;

  partitioner_ = std::make_unique<GeographyPartitioner>(world_, config_);
  communicator_ =
      std::make_unique<DomainCommunicator>(world_, config_, domain_);

  if (rank_ != 0) {
    loadGeographyOnNonZeroRanks();
  }

  partitioner_->loadGeographyMetadata();
  if (rank_ == 0) {
    partitioner_->loadPopulationWeights(world_state_file_);
  }
  broadcastPopulationCounts();

  std::unordered_set<GeoUnitId> owned_units;
  partitioner_->partition(num_ranks_, rank_, owned_units);
  for (auto id : owned_units) domain_.addGeoUnit(id);

  loadDomainData();
  domain_.assignPeopleAndVenues();
  buildVenueOwnershipMap();
  loadGlobalPersonMetadata();

  exchangeActivityMasks();
  buildGlobalVenueOwnershipMap();
  computeGlobalMaxPersonId();

  reportDomainStats("Initial Load");
}

void DomainManager::broadcastPopulationCounts() {
  auto& geo_data = partitioner_->getGeoData();
  size_t num_units_size = geo_data.size();
  uint64_t num_units = static_cast<uint64_t>(num_units_size);
  MPI_Bcast(&num_units, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  if (num_units <= 0) return;

  std::vector<int32_t> pops(num_units);
  std::vector<int32_t> name_lengths(num_units);
  uint64_t total_names_len = 0;

  if (rank_ == 0) {
    int i = 0;
    for (const auto& [name, data] : geo_data) {
      pops[i] = static_cast<int32_t>(data.population);
      int32_t len = static_cast<int32_t>(name.length());
      name_lengths[i] = len;
      total_names_len += len;
      i++;
    }
  }

  // 1. Broadcast populations and name lengths (chunked for > 2GB safety)
  chunkedBroadcast(pops.data(), num_units, sizeof(int32_t), MPI_INT32_T, 0);
  chunkedBroadcast(name_lengths.data(), num_units, sizeof(int32_t), MPI_INT32_T,
                   0);

  // 2. Broadcast total size of names buffer
  MPI_Bcast(&total_names_len, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

  // 3. Batch broadcast all names (CHUNKING for > 2GB safety)
  std::vector<char> names_buffer(total_names_len);
  if (rank_ == 0) {
    size_t offset = 0;
    int i = 0;
    for (const auto& [name, data] : geo_data) {
      memcpy(names_buffer.data() + offset, name.c_str(), name_lengths[i]);
      offset += name_lengths[i];
      i++;
    }
  }

  chunkedBroadcast(names_buffer.data(), total_names_len, sizeof(char), MPI_CHAR,
                   0);

  // 4. Update local geo_data on non-zero ranks
  if (rank_ != 0) {
    size_t offset = 0;
    for (uint64_t i = 0; i < num_units; ++i) {
      std::string name(names_buffer.data() + offset, name_lengths[i]);
      offset += name_lengths[i];

      if (geo_data.count(name)) {
        geo_data[name].population = pops[i];
      }
    }
  }
}

void DomainManager::loadDomainData() {
  world_ = HDF5Loader::loadDomainChunked(
      world_state_file_, domain_.geo_unit_set,
      config_.parallel.geo_unit_chunk_size, config_);
}

void DomainManager::buildVenueOwnershipMap() {
  // Gather which rank owns which partition-level GeoUnit
  // Get the base partition units (e.g. MGUs) owned by this rank
  std::vector<GeoUnitId> local_partition_units;
  const std::string& level = config_.parallel.partition_level;

  for (GeoUnitId id : domain_.geo_unit_ids) {
    const GeographicalUnit* gu = world_.getGeoUnit(id);
    if (!gu) continue;
    std::string current_level = (gu->level_id < world_.geo_level_names.size())
                                    ? world_.geo_level_names[gu->level_id]
                                    : "unknown";
    if (current_level == level) {
      local_partition_units.push_back(id);
    }
  }

  int local_count = local_partition_units.size();
  std::vector<int> counts(num_ranks_);
  MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD);

  std::vector<int> displs;
  int total;
  mpi_utils::computeDisplacements(counts, displs, total);

  std::vector<GeoUnitId> all_gu_ids(total);
  MPI_Allgatherv(local_partition_units.data(), local_count, MPI_INT,
                 all_gu_ids.data(), counts.data(), displs.data(), MPI_INT,
                 MPI_COMM_WORLD);

  geounit_to_rank_.clear();
  for (int r = 0; r < num_ranks_; ++r) {
    for (int i = 0; i < counts[r]; ++i) {
      geounit_to_rank_[all_gu_ids[displs[r] + i]] = r;
    }
  }
}

void DomainManager::buildGlobalVenueOwnershipMap() {
  // Each rank shares its local venue IDs with all other ranks so that
  // cross-rank venue references (from activity mappings) can be resolved.
  const auto& local_ids = domain_.local_venue_ids;
  int local_count = static_cast<int>(local_ids.size());

  std::vector<int> counts(num_ranks_);
  MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD);

  std::vector<int> displs;
  int total;
  mpi_utils::computeDisplacements(counts, displs, total);

  std::vector<VenueId> all_venue_ids(total);
  MPI_Allgatherv(local_ids.data(), local_count, MPI_INT, all_venue_ids.data(),
                 counts.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);

  global_venue_rank_.clear();
  global_venue_rank_.reserve(total);
  int duplicates = 0;
  for (int r = 0; r < num_ranks_; ++r) {
    for (int i = 0; i < counts[r]; ++i) {
      VenueId vid = all_venue_ids[displs[r] + i];
      auto it = global_venue_rank_.find(vid);
      if (it != global_venue_rank_.end() && it->second != r) {
        duplicates++;
      }
      global_venue_rank_[vid] = r;
    }
  }

  if (rank_ == 0 && duplicates > 0) {
    std::cerr << "  Warning: " << duplicates
              << " venues claimed by multiple ranks!" << std::endl;
  }
}

int DomainManager::getVenueRank(VenueId vid) const {
  // 1. Check dynamic registry first (set by setVenueRank for virtual venues)
  auto git = global_venue_rank_.find(vid);
  if (git != global_venue_rank_.end()) return git->second;

  // 2. Virtual venues not in the registry are unknown
  if (vid < 0) return -1;

  // 3. Handle Physical Venues (local only, legacy path)
  const Venue* venue = world_.getVenue(vid);
  if (!venue) return -1;

  // A venue is owned by the rank that owns its partition-level parent
  GeoUnitId partition_gu_id =
      (partitioner_) ? partitioner_->findParentAtLevel(
                           venue->geo_unit_id, config_.parallel.partition_level)
                     : -1;

  if (partition_gu_id == -1) {
    // If no partition parent found (e.g. venue at higher level than partition),
    // we need a deterministic assignment.
    // Try to find ANY descendant that is at partition level.
    partition_gu_id = venue->geo_unit_id;

    // Simple search for a child at partition_level
    for (const auto& gu : world_.geo_units) {
      if (partitioner_ &&
          partitioner_->findParentAtLevel(
              gu.id, config_.parallel.partition_level) == venue->geo_unit_id) {
        // gu is a child/descendant of the venue's LGU
        // We can use gu.id as the partition_gu_id
        partition_gu_id = gu.id;
        break;
      }
    }

    // If still not found in map, default to rank 0 to avoid -1
    auto it = geounit_to_rank_.find(partition_gu_id);
    if (it == geounit_to_rank_.end()) {
      return 0;  // Default to rank 0 for "orphaned" venues
    }
    return it->second;
  }

  auto it = geounit_to_rank_.find(partition_gu_id);
  return (it != geounit_to_rank_.end()) ? it->second : -1;
}

int DomainManager::getPersonRank(PersonId pid) const {
  if (pid < 0 || pid >= static_cast<PersonId>(global_person_rank_.size()))
    return -1;
  return global_person_rank_[pid];
}

void DomainManager::loadGlobalPersonMetadata() {
  // First, find max PersonId to pre-allocate
  size_t chunk_size = config_.parallel.person_metadata_chunk_size;
  if (chunk_size == 0) chunk_size = 100000;

  // First pass or use HDF5 info to find size
  HDF5Loader::loadPersonMetadataChunked(
      world_state_file_, chunk_size,
      [&](const std::vector<HDF5Loader::PersonMetadata>& chunk) {
        for (const auto& m : chunk) {
          if (m.person_id > max_person_id_) max_person_id_ = m.person_id;
        }
      });

  global_person_geounit_.assign(max_person_id_ + 1, -1);
  global_person_rank_.assign(max_person_id_ + 1, -1);

  // Second pass: fill the vector
  HDF5Loader::loadPersonMetadataChunked(
      world_state_file_, chunk_size,
      [&](const std::vector<HDF5Loader::PersonMetadata>& chunk) {
        for (const auto& m : chunk) {
          if (m.person_id >= 0) {
            global_person_geounit_[m.person_id] = m.geo_unit_id;

            // Pre-calculate and cache the owning rank for each person
            GeoUnitId partition_gu_id = partitioner_->findParentAtLevel(
                m.geo_unit_id, config_.parallel.partition_level);
            if (partition_gu_id == -1) partition_gu_id = m.geo_unit_id;

            auto it = geounit_to_rank_.find(partition_gu_id);
            global_person_rank_[m.person_id] =
                (it != geounit_to_rank_.end()) ? it->second : -1;
          }
        }
      });
}

void DomainManager::exchangeVisitors(
    const std::vector<PersonLocation>& locations, double current_time,
    double delta_hours, const RuntimeGroupAllocator* alloc) {
  communicator_->exchangeVisitors(locations, *this, current_time, delta_hours,
                                  alloc);
}

std::unordered_set<PersonId> DomainManager::getVisitorIds() const {
  std::unordered_set<PersonId> ids;
  for (const auto& v : domain_.incoming_visitors) ids.insert(v.person_id);
  return ids;
}

std::vector<PendingInfection> DomainManager::receivePendingInfections(
    const std::vector<PendingInfection>& pending) {
  return communicator_->receivePendingInfections(pending);
}

void DomainManager::exchangeEncounterProposals(
    const std::vector<EncounterProposal>& local_proposals,
    std::vector<EncounterProposal>& proposals_for_this_rank) {
  communicator_->exchangeEncounterProposals(local_proposals, *this,
                                            proposals_for_this_rank);
}

void DomainManager::exchangeEncounterReplies(
    const std::vector<EncounterReply>& local_replies,
    std::vector<EncounterReply>& replies_for_this_rank) {
  communicator_->exchangeEncounterReplies(local_replies, *this,
                                          replies_for_this_rank);
}

void DomainManager::exchangeFinalizedEncounters(
    const std::vector<CoordinatedEncounter>& local_finalized,
    std::vector<CoordinatedEncounter>& finalized_for_this_rank) {
  communicator_->exchangeFinalizedEncounters(local_finalized, *this,
                                             finalized_for_this_rank);
}

template <typename T>
void DomainManager::exchangeGlobalProperty(
    std::vector<T>& global_buf, T null_value, MPI_Datatype mpi_type,
    MPI_Op mpi_op, std::function<T(const Person&)> extract) {
  if (num_ranks_ <= 1) return;

  global_buf.assign(max_person_id_ + 1, null_value);
  std::vector<T> local_buffer(max_person_id_ + 1, null_value);

  for (const auto& person : world_.people) {
    if (domain_.ownsPerson(person.id)) {
      local_buffer[person.id] = extract(person);
    }
  }

  MPI_Allreduce(local_buffer.data(), global_buf.data(),
                static_cast<int>(global_buf.size()), mpi_type, mpi_op,
                MPI_COMM_WORLD);
}

void DomainManager::exchangeDeathFlags() {
  exchangeGlobalProperty<uint8_t>(
      global_death_flags_, 0, MPI_UINT8_T, MPI_MAX,
      [&](const Person& p) -> uint8_t { return p.is_dead ? 1 : 0; });
}

void DomainManager::exchangeScheduleTypes() {
  exchangeGlobalProperty<uint16_t>(
      global_person_schedule_type_, 65535, MPI_UNSIGNED_SHORT, MPI_MIN,
      [](const Person& p) -> uint16_t { return p.schedule_type_id; });
}

uint16_t DomainManager::getGlobalScheduleType(PersonId pid) const {
  if (pid < 0 ||
      pid >= static_cast<PersonId>(global_person_schedule_type_.size())) {
    return 65535;
  }
  return global_person_schedule_type_[pid];
}

void DomainManager::setGlobalScheduleType(PersonId pid, uint16_t type_id) {
  if (pid < 0) return;
  if (pid >= static_cast<PersonId>(global_person_schedule_type_.size())) {
    global_person_schedule_type_.resize(pid + 1, 65535);
  }
  global_person_schedule_type_[pid] = type_id;
}

void DomainManager::exchangeActivityMasks() {
  if (rank_ == 0)
    std::cout << "MPI: Exchanging global activity/venue availability masks..."
              << std::endl;

  size_t num_people = max_person_id_ + 1;
  global_person_activity_mask_.assign(num_people, 0);

  // MPI has no built-in type for __uint128_t, so we split each ActivityMask
  // into its low 64 bits and high 64 bits and do two separate Allreduce calls.
  std::vector<uint64_t> local_lo(num_people, 0), local_hi(num_people, 0);

  for (const auto& person : world_.people) {
    if (domain_.ownsPerson(person.id)) {
      ActivityMask mask = 0;
      auto metas = world_.getActivityMetas(person);
      for (const auto& meta : metas) {
        if (meta.venue_count > 0 && meta.activity_index >= 0) {
          mask |= (ActivityMask(1) << meta.activity_index);
        }
      }
      local_lo[person.id] = static_cast<uint64_t>(mask);
      local_hi[person.id] = static_cast<uint64_t>(mask >> 64);
    }
  }

  // Combine all ranks into the global array.
  std::vector<uint64_t> global_lo(num_people, 0), global_hi(num_people, 0);
  MPI_Allreduce(local_lo.data(), global_lo.data(), static_cast<int>(num_people),
                MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(local_hi.data(), global_hi.data(), static_cast<int>(num_people),
                MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  for (size_t i = 0; i < num_people; ++i) {
    global_person_activity_mask_[i] =
        (ActivityMask(global_hi[i]) << 64) | ActivityMask(global_lo[i]);
  }

  if (rank_ == 0) {
    std::cout << "  Synchronized activity masks for "
              << global_person_activity_mask_.size() << " people." << std::endl;
  }
}

ActivityMask DomainManager::getGlobalActivityMask(PersonId pid) const {
  if (pid < 0 ||
      pid >= static_cast<PersonId>(global_person_activity_mask_.size())) {
    return 0;
  }
  return global_person_activity_mask_[pid];
}

void DomainManager::setGlobalActivityMask(PersonId pid, ActivityMask mask) {
  if (pid < 0) return;
  if (pid >= static_cast<PersonId>(global_person_activity_mask_.size())) {
    global_person_activity_mask_.resize(pid + 1, 0);
  }
  global_person_activity_mask_[pid] = mask;
}

void DomainManager::setPersonRank(PersonId pid, int rank) {
  if (pid < 0) return;
  if (pid >= static_cast<PersonId>(global_person_geounit_.size())) {
    global_person_geounit_.resize(pid + 1, -1);
  }
  // Create a dummy geo unit that we map to the requested rank
  GeoUnitId dummy_gu = static_cast<GeoUnitId>(pid + 100000);
  global_person_geounit_[pid] = dummy_gu;
  geounit_to_rank_[dummy_gu] = rank;

  if (pid >= static_cast<PersonId>(global_person_rank_.size())) {
    global_person_rank_.resize(pid + 1, -1);
  }
  global_person_rank_[pid] = rank;
}

void DomainManager::setGeoUnitRank(GeoUnitId guid, int rank) {
  geounit_to_rank_[guid] = rank;
}

void DomainManager::synchronizeRegistries() {
  // Each rank may have discovered property values in a different order
  // during chunked loading.  Build a globally-consistent registry by
  // gathering all unique values and sorting them deterministically.
  // Then remap every person property integer code to the new indices.

  for (auto& [prop_name, local_registry] :
       world_.person_property_value_registries) {
    // 1. Gather all unique values from all ranks (null-separated packed)
    std::vector<char> all_packed =
        gatherPackedValuesAcrossRanks(local_registry, num_ranks_);

    // 2. Build a sorted, deduplicated global registry
    std::vector<std::string> global_registry =
        buildGlobalRegistryFromPacked(all_packed);

    // 3. Build remap from old local index → new global index
    std::vector<int32_t> remap =
        buildRemapToGlobal(local_registry, global_registry);

    // 4. Remap all person property values for this property
    remapPersonPropertyValues(world_, prop_name, remap);

    // 5. Replace local registry with global registry
    local_registry = global_registry;
  }
}

void DomainManager::reportDomainStats(const std::string& label) const {
  size_t local_pop = world_.people.size();
  size_t local_venues = world_.venues.size();
  double local_rss_gb = MemoryUtils::getRSS() / (1024.0 * 1024.0);

  // Exchange stats
  std::vector<uint64_t> pops(num_ranks_);
  std::vector<uint64_t> venues(num_ranks_);
  std::vector<double> rss(num_ranks_);

  uint64_t lp = local_pop;
  uint64_t lv = local_venues;

  MPI_Allgather(&lp, 1, MPI_UINT64_T, pops.data(), 1, MPI_UINT64_T,
                MPI_COMM_WORLD);
  MPI_Allgather(&lv, 1, MPI_UINT64_T, venues.data(), 1, MPI_UINT64_T,
                MPI_COMM_WORLD);
  MPI_Allgather(&local_rss_gb, 1, MPI_DOUBLE, rss.data(), 1, MPI_DOUBLE,
                MPI_COMM_WORLD);

  if (rank_ == 0) printDomainStats(label, num_ranks_, pops, venues, rss);
}

}  // namespace june

#endif  // USE_MPI
