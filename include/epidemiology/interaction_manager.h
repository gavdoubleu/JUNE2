#pragma once

#include <map>
#include <optional>
#include <random>
#include <unordered_set>
#include <vector>

#include "activity/activity_manager.h"
#include "activity/runtime_group_allocator.h"
#include "core/config.h"
#include "core/world_state.h"
#include "disease.h"
#include "policy.h"
#include "utils/age_utils.h"
#include "utils/event_logging/event_logger.h"
#include "utils/event_logging/event_types.h"
#include "utils/time_utils.h"

namespace june {

class CompartmentalModelManager;  // optional plugin; null when inactive

// =============================================================================
// InteractionMember - Lightweight person reference
// =============================================================================
struct InteractionMember {
  PersonId id;
  size_t array_index;
  SubsetIndex subset_index;
  uint8_t encounter_type_id;
};

// =============================================================================
// Structures used in transmission processing (moved to header for reuse)
// =============================================================================

struct SusceptibleMember {
  PersonId id;
  double susceptibility;
  const VisitorInfo* visitor;
  uint8_t encounter_type_id;
};

struct BinGroup {
  std::vector<SusceptibleMember> susceptible;
  std::vector<PersonId> infectious_ids;
  // [mode_index][person_idx_within_bin]: per-mode infectiousness
  std::vector<std::vector<double>> infectiousness_by_mode;
  std::vector<double> total_infectiousness_by_mode;
  int total_size = 0;
  // [mode_index] cumulative weights derived from I_{b,m} for that mode.
  // Empty if the mode has no positive weight. Built once per bin/mode in
  // STEP 2 of processVenueTransmissions and sampled many times in STEP 3b
  // via sampleFromCumulative, replacing a per-call std::discrete_distribution
  // construction (the 116M "smoking gun" from the 60M profile).
  std::vector<std::vector<double>> cumulative_by_mode;
  // Per-fomite-mode deposition totals for this bin
  // [local_fomite_mode_index][sub_bin_k]
  std::vector<std::vector<double>> total_fomite_deposition_sub;

  // Clear person-proportional fields post-transmission while cache-hot.
  // Does NOT touch total_fomite_deposition_sub; that is handled exclusively
  // by initFomiteSubBins, called pre-use with the correct slot parameters.
  void clearAfterUse(int num_modes) {
    susceptible.clear();
    infectious_ids.clear();
    infectiousness_by_mode.resize(num_modes);
    for (auto& v : infectiousness_by_mode) v.clear();
    total_infectiousness_by_mode.assign(num_modes, 0.0);
    total_size = 0;
    cumulative_by_mode.resize(num_modes);
    for (auto& c : cumulative_by_mode) c.clear();
  }

  // Size and zero fomite sub-bin accumulators for the current slot's
  // delta_hours. Called pre-use so total_fomite_deposition_sub[fm] has exactly
  // n_sub_per_mode[fm] elements before any deposition is accumulated into them.
  void initFomiteSubBins(int num_fomite_modes,
                         const std::vector<int>& n_sub_per_mode) {
    total_fomite_deposition_sub.resize(num_fomite_modes);
    for (int fm = 0; fm < num_fomite_modes; ++fm) {
      int n_sub = fm < (int)n_sub_per_mode.size() ? n_sub_per_mode[fm] : 1;
      total_fomite_deposition_sub[fm].assign(n_sub, 0.0);
    }
  }
};

struct SourceEntry {
  int mode;
  int inf_bin;
  // For sibling-mixing sources, inf_bin == SIBLING_INF_BIN_SENTINEL and
  // sibling_parent_inf_bin holds the bin under the PARENT matrix to use
  // when sampling the infector from the per-parent flat list.
  int sibling_parent_inf_bin = -1;
};
// Sentinel for SourceEntry::inf_bin meaning "this source is sibling-venue
// mixing; sample infector from ParentAggregate". Distinct from -1 (fomite)
// and -2 (compartmental).
constexpr int SIBLING_INF_BIN_SENTINEL = -3;

// =============================================================================
// ParentAggregate - per-tick aggregation of infectiousness across all child
// venues of a single parent venue (e.g. classrooms of one school). Built in
// processTransmissions::buildParentAggregates() before the main per-venue
// loop. Each child venue (classroom/office/uni_groups) reads from its
// ParentAggregate to add a "sibling-mixing" force-of-infection term.
//
// All child venues of a parent live on the same MPI rank (children share
// the parent's MGU), so this map is rank-local with no MPI communication.
// =============================================================================
struct ParentInfectorEntry {
  PersonId person_id;
  VenueId child_venue_id;  // origin child, used to exclude own-venue at sample
  // Per-mode integrated infectiousness (24*∫I dt), same units as
  // BinGroup::infectiousness_by_mode entries.
  std::vector<double> inf_by_mode;
};

struct ParentAggregate {
  // [parent_bin][mode] → sum of integrated infectiousness across all children
  std::vector<std::vector<double>> total_inf_by_bin_mode;
  // [parent_bin] → headcount across all children
  std::vector<int> size_by_bin;
  // [parent_bin] → flat list of infectors from any child, used for sibling
  // infector sampling (filter by child_venue_id at sample time to exclude own)
  std::vector<std::vector<ParentInfectorEntry>> infectors_by_bin;

  // Per-child contributions, so each child can subtract its own share when
  // computing the sibling FOI it sees. Keyed by child venue id.
  // child_size[child_venue_id][parent_bin]
  std::unordered_map<VenueId, std::vector<int>> child_size_by_bin;
  // child_inf[child_venue_id][parent_bin][mode]
  std::unordered_map<VenueId, std::vector<std::vector<double>>>
      child_inf_by_bin_mode;
  // Cached parent venue metadata (filled at first insertion)
  uint8_t parent_venue_type_id = kUnknownVenueTypeId;
};

// =============================================================================
// RuntimeGroupMember: per-member record in a runtime group of a
// partial-presence venue. Built by buildPartialPresenceGroups, consumed
// by the sub-interval accumulator inside computePartialPresenceLambda.
// =============================================================================
struct RuntimeGroupMember {
  PersonId pid;
  size_t array_index;
  SubsetIndex subset_index;
  uint8_t enc_type_id;
  Person* person;  // null for visitors
  const VisitorInfo* visitor;
  float eff_board;
  float eff_alight;
  // Per-rider presence cap f_p = min(1, slot / T_p). Scales this rider's FOI
  // contribution as BOTH source (f_I) and sink (f_S) so its total modeled
  // commute presence stays ≤ one slot-hour. 1.0 for journeys that fit a slot.
  float f_presence;
  int matrix_bin;
};

// =============================================================================
// FomiteModeRef: per-FomiteConfig handle pre-resolved before the binning
// loop in processVenueTransmissions. mode_index points into
// disease_->getTransmissionParams().modes; cfg points into the same
// underlying TransmissionMode::config variant.
// =============================================================================
struct FomiteModeRef {
  int mode_index;
  const FomiteConfig* cfg;
};

// =============================================================================
// PartialPresenceSubBin: per-(matrix-bin) scratch accumulator reused across
// sub-intervals in computePartialPresenceLambda. reset() clears storage to
// start a new sub-interval.
// =============================================================================
struct PartialPresenceSubBin {
  std::vector<double> total_inf_by_mode;  // size num_modes
  std::vector<PersonId> infectious_ids;
  // Per (mode) → flat list aligned with infectious_ids for sampling.
  std::vector<std::vector<double>> inf_per_person_by_mode;  // [mode][i]
  int total_size = 0;
  void reset(int num_modes_) {
    total_inf_by_mode.assign(num_modes_, 0.0);
    infectious_ids.clear();
    inf_per_person_by_mode.assign(num_modes_, {});
    total_size = 0;
  }
};

// =============================================================================
// InteractionManager - Handles disease transmission through interactions
// =============================================================================

class InteractionManager {
 public:
  struct PerformanceStats {
    size_t matrix_lookups = 0;
    size_t bin_lookups = 0;
    size_t grouping_ops = 0;

    void print() const {
      std::cout << "\n=== InteractionManager Performance Stats ==="
                << std::endl;
      std::cout << "  Matrix lookups:      " << matrix_lookups << std::endl;
      std::cout << "  Bin lookups:         " << bin_lookups << std::endl;
      std::cout << "  Grouping operations: " << grouping_ops << std::endl;
    }
  };

  InteractionManager(WorldState& world,
                     const ContactMatrixConfig& contact_matrices,
                     const SimulationConfig& simulation_config,
                     const ParallelConfig& parallel_config, Disease* disease,
                     EventLogger* event_logger = nullptr);

  // Calculate and apply transmissions for a time slot.
  // Returns the number of new infections.
  // Pass comp_model to enable compartmental uptake FOI; null = no plugin.
  int processTransmissions(
      const std::vector<PersonLocation>& locations, double current_time,
      double delta_hours,
      std::unordered_set<PersonId>* active_infections = nullptr,
      const std::unordered_set<PersonId>* visitor_ids = nullptr,
      std::vector<PendingInfection>* pending_infections = nullptr,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data = nullptr,
      const CompartmentalModelManager* comp_model = nullptr);

  PerformanceStats& getStats() { return stats_; }

  // Set the current day type index (used by EventLogger for per-type stats)
  void setCurrentDayTypeIdx(int idx) { current_day_type_idx_ = idx; }

  // Wire the runtime group allocator so partial-presence venues (commute lines,
  // etc.) can resolve per-rider group assignments and presence windows.
  // Caller (the Simulator) sets this once after both objects exist. Null is
  // a valid state; partial-presence venues degrade to full-slot FOI without
  // an allocator.
  void setRuntimeGroupAllocator(const RuntimeGroupAllocator* a) {
    runtime_group_allocator_ = a;
  }

  // The force-of-infection pass for transport lines, run after the ordinary
  // venue loop. Lines are not in that loop: a line's members come from the
  // allocator's rider table, since a rider on four legs holds only one
  // location. Each line records who it would infect, and nothing is applied.
  //
  // owned_lines must be the lines this rank owns, in ascending order.
  int processPartialPresenceLines(
      const std::vector<VenueId>& owned_lines, double current_time,
      double delta_hours, std::unordered_set<PersonId>* active_infections,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data);

  // Settle the candidates the pass above collected: one infection per rider,
  // taken from the lowest venue id that drew for them. Collective, and every
  // rank reaches the same verdict from the same gathered list, so the outcome
  // does not depend on how the lines were spread across ranks. Returns how many
  // infections this rank applied to its own residents.
  int resolvePartialPresenceInfections(
      double current_time, std::unordered_set<PersonId>* active_infections);

  // Per-susceptible accumulated λ + source attribution from a single
  // partial-presence venue. Output of computePartialPresenceLambda; consumed
  // by processPartialPresenceVenue's Bernoulli step. Exposed publicly so
  // engine tests can assert the FOI math without driving the stochastic
  // infection draw.
  struct PartialPresenceAccumSource {
    int mode;
    PersonId infector;
    double weighted;
  };
  struct PartialPresenceLambdaResult {
    std::unordered_map<PersonId, double> susc_lambda;
    std::unordered_map<PersonId, std::vector<PartialPresenceAccumSource>>
        susc_sources;
  };

  // Steps 1–2 of partial-presence transmission processing, extracted as a
  // pure accumulator: bucket members into runtime bins ("groups"), walk
  // sub-intervals delimited by effective presence windows, accumulate
  // per-susceptible λ and per-source attribution weights. No infection
  // side-effects, no RNG. Preconditions are enforced here (throws on
  // violation) so callers can't accidentally bypass them.
  //
  // The wrapping processPartialPresenceVenue calls this then performs the
  // single Bernoulli draw + infection write per susceptible. Tests call
  // this directly to assert λ deterministically.
  PartialPresenceLambdaResult computePartialPresenceLambda(
      const std::vector<InteractionMember>& members, Venue* venue,
      VenueId actual_venue_id, double current_time, double delta_hours,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      uint8_t encounter_type_id);

 private:
  WorldState& world_;
  const ContactMatrixConfig& contact_matrices_;
  const SimulationConfig& simulation_config_;
  const ParallelConfig& parallel_config_;
  Disease* disease_;
  EventLogger* event_logger_;
  int current_day_type_idx_ = 0;

  // Pre-pass over locations (already sorted by venue) building
  // parent_aggregates_ for the current tick. Called once per
  // processTransmissions invocation before the per-venue loop. Iterates each
  // venue group, computes per-bin per-mode integrated infectiousness under
  // the PARENT's contact matrix, and stores per-child contributions so each
  // child can subtract its own share later. Walks members in person_id order
  // (the same order the main loop uses) for deterministic FP accumulation.
  void buildParentAggregates(
      double current_time, double delta_hours,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data);

  // Emit the per-tick parent-aggregate summary lines under
  // JUNE_DEBUG_PARENT_MIXING. Iterates parent_aggregates_ in sorted-key order
  // (deterministic across runs) and prints headline counts plus details of
  // up to 3 parents with non-zero infectiousness.
  void dumpParentAggregatesDebug(double current_time, double delta_hours) const;

  // Update parent_aggregates_ for the venue group [group_start, group_end)
  // (already grouped by venue_id in active_locations_buffer_). Returns early
  // if the venue is virtual, has no parent, or the parent's contact matrix
  // is missing. Walks members in person_id order so sibling FP sums match
  // the main loop's STEP 1 ordering across rank counts.
  void aggregateOneVenueGroupForParent(
      size_t group_start, size_t group_end, double current_time,
      double delta_hours, int num_modes,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data);

  // Look up parent_aggregates_[parent_id], lazily initialising its per-bin
  // arrays + the per-child-venue contribution entries. Returns the
  // ParentAggregate reference so the caller can accumulate into it directly.
  ParentAggregate& ensureParentAggregateInitialised(VenueId parent_id,
                                                    VenueId child_venue_id,
                                                    uint8_t parent_type_id,
                                                    int parent_num_bins,
                                                    int num_modes);

  // Fill inf_by_mode with the per-mode integrated infectiousness (24*∫I dt)
  // contributed by this member over [current_time, current_time+delta_hours].
  // Visitor branch reads pre-computed values from the sending rank; local
  // branch calls Infection::getIntegratedInfectiousness. Returns true iff
  // the member contributes any positive infectiousness.
  bool gatherMemberInfectiousnessByMode(const Person* person,
                                        const VisitorInfo* visitor,
                                        double current_time, double delta_hours,
                                        int num_modes,
                                        std::vector<double>& inf_by_mode) const;

  // Copy active_locations_buffer_[group_start..group_end) into a fresh vector
  // and sort it by person_id. Used by the parent-aggregate pre-pass to walk
  // venue members in the same person_id order as the main loop's STEP 1,
  // keeping the sibling FP sum order bit-identical across rank counts.
  std::vector<PersonLocation> buildPersonIdSortedMembers(
      size_t group_start, size_t group_end) const;

  // Resolve `pid` to a local Person* (preferring the hint at array_index)
  // or to a VisitorInfo* from visitor_data. Returns true iff at least one
  // resolved. Either out param is set to nullptr if not found.
  bool resolvePersonAndVisitor(
      PersonId pid, size_t array_index,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      Person*& person_out, const VisitorInfo*& visitor_out) const;

  // Throws std::runtime_error if any v1 precondition of partial-presence FOI
  // is violated for this venue/encounter, checking every member rather than the
  // group's representative. Pure check, no side effects.
  void validatePartialPresencePreconditions(
      const std::vector<InteractionMember>& members, const Venue* venue,
      VenueId actual_venue_id, uint8_t encounter_type_id) const;

  // Step 1 of computePartialPresenceLambda: resolve each member's group
  // (via runtime_group_allocator_), matrix_bin, and effective presence window,
  // and group them into one group bucket each. Each bucket is sorted by
  // person_id at the end so per-group FP accumulation matches across
  // ranks. Returns one bucket per group; some may be empty.
  std::vector<std::vector<RuntimeGroupMember>> buildPartialPresenceGroups(
      const std::vector<InteractionMember>& members, Venue* venue,
      VenueId actual_venue_id, const ContactMatrix& bin_structure,
      int num_bins_needed, uint16_t num_groups,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data) const;

  // Collect the group's raw {eff_board, eff_alight} offsets as event times,
  // sort and de-duplicate within 1e-5 tolerance. Sub-intervals span the
  // venue's own [min board, max alight] range (the windows are raw line-local
  // offsets, not [0, slot]). The accumulator walks the sub-intervals delimited
  // by consecutive entries. Returns the deduped event times.
  std::vector<float> collectSubIntervalEventTimes(
      const std::vector<RuntimeGroupMember>& group,
      float slot_duration_min) const;

  // For sub-interval [t0, t1) of group `group`, walk every member present
  // throughout the sub-interval and classify it: (a) infectious, append to
  // sub_bins[bin].infectious_ids + push per-mode integrated infectiousness
  // scaled by `scale`; (b) susceptible, append &member to
  // susc_by_bin[bin]; (c) dead/no role, only update headcount. sub_bins is
  // pre-reset by the caller for this sub-interval.
  void classifyMembersInSubInterval(
      const std::vector<RuntimeGroupMember>& group, float t0, float t1,
      double scale, double current_time, double delta_hours, int num_modes,
      std::vector<PartialPresenceSubBin>& sub_bins,
      std::vector<std::vector<const RuntimeGroupMember*>>& susc_by_bin) const;

  // Return the contacts entry for (susc_bin, inf_bin), preferring
  // The per-(mode, bin-pair) contacts lookup used by both the main FOI loop
  // and partial-presence. The matrix is whatever getMatrix/getVirtualMatrix
  // resolved for this (type, mode), which is always a complete matrix, so an
  // out-of-range bin pair is a real bug (loops sized off the wrong thing) and
  // throws rather than inventing a contact rate.
  double lookupContactsForBinPair(const ContactMatrix& matrix, int susc_bin,
                                  int inf_bin) const;

  // For one (group, sub-interval), iterate (susc_bin, mode, inf_bin) over
  // pre-classified sub_bins + susc_by_bin and accumulate per-susceptible λ
  // contributions into susc_lambda and per-(susc, source) attribution
  // AccumSources into susc_sources. Pure FOI-math step; no infection writes,
  // no RNG. The mode-matrix lookup and the own-bin minus-one adjustment
  // match the main-loop semantics.
  void accumulatePartialLambdaContributions(
      const std::vector<PartialPresenceSubBin>& sub_bins,
      const std::vector<std::vector<const RuntimeGroupMember*>>& susc_by_bin,
      uint8_t venue_type_id, const ContactMatrix& bin_structure,
      int num_bins_needed, int num_modes,
      const TransmissionParams& trans_params,
      PartialPresenceLambdaResult& result) const;

  // Walk the sub-intervals of one group, classifying members and
  // accumulating per-susceptible λ + AccumSources into result. sub_bins is
  // the caller-owned scratch buffer reused per sub-interval.
  void accumulateOneGroup(const std::vector<RuntimeGroupMember>& group,
                          float slot_duration_min, double current_time,
                          double delta_hours, int num_modes,
                          int num_bins_needed, uint8_t venue_type_id,
                          const ContactMatrix& bin_structure,
                          const TransmissionParams& trans_params,
                          std::vector<PartialPresenceSubBin>& sub_bins,
                          PartialPresenceLambdaResult& result) const;

  // Per-member body of the parent-aggregate pre-pass: resolve person/visitor,
  // compute parent_bin under parent_matrix, bump headcount in agg/csize,
  // gather per-mode infectiousness, and (if positive) accumulate into
  // agg.total_inf_by_bin_mode + cinf + agg.infectors_by_bin.
  void accumulateOneMemberIntoParent(
      const PersonLocation& loc, Venue* venue,
      const ContactMatrix* parent_matrix, int parent_num_bins,
      ParentAggregate& agg, std::vector<int>& csize,
      std::vector<std::vector<double>>& cinf, VenueId child_venue_id,
      double current_time, double delta_hours, int num_modes,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      std::vector<double>& inf_by_mode_scratch) const;

  // Populate active_locations_buffer_ with non-unallocated entries from
  // `locations`, then sort by (venue_id, encounter_type_id-when-virtual).
  // Bumps stats_.grouping_ops.
  void filterAndSortActiveLocations(
      const std::vector<PersonLocation>& locations);

  // Count local (non-visitor) members of the current venue group whose
  // encounter_type_id names a coordinated encounter, and log via event_logger_
  // if any. group_start/group_end index into active_locations_buffer_.
  void logCoordinatedEncounterParticipants(
      size_t group_start, size_t group_end,
      const std::unordered_set<PersonId>* visitor_ids);

  // Fast pre-check: true iff processVenueTransmissions could possibly produce
  // transmission for this venue group, considering fomite history,
  // compartmental uptake, and the presence of any infectious member in
  // group_members_buffer_. Used to skip venues with zero infectious source
  // before paying the FOI cost.
  bool venueGroupHasTransmissionSource(
      const Venue* venue, VenueId venue_id,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      const CompartmentalModelManager* comp_model) const;

  // Starting at group_start (an index into active_locations_buffer_), advance
  // past every contiguous entry that shares the same venue_id (and, for virtual
  // venues, encounter_type_id) and copy each into group_members_buffer_. The
  // buffer is cleared first, then sorted by person_id at the end so that
  // FP accumulation in STEP 1 is deterministic across rank counts.
  // Returns the index one past the end of the group.
  size_t collectAndSortGroupMembers(size_t group_start);

  // Emit the per-tick parent-mixing summary line when JUNE_DEBUG_PARENT_MIXING
  // is enabled and any infections occurred this tick.
  void printTickParentMixingSummary(int total_new_infections,
                                    double current_time) const;

  // Process one contiguous venue group already loaded into
  // group_members_buffer_ (between group_start and group_end in
  // active_locations_buffer_). Resolves the venue, logs coordinated-encounter
  // participants, runs the fast pre-check, and dispatches to
  // processVenueTransmissions. Returns new infections (0 if group is skipped).
  int processOneVenueGroup(
      size_t group_start, size_t group_end, double current_time,
      double delta_hours, std::unordered_set<PersonId>* active_infections,
      const std::unordered_set<PersonId>* visitor_ids,
      std::vector<PendingInfection>* pending_infections,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      const CompartmentalModelManager* comp_model);

  // Compute the bin index for a person under a given contact matrix.
  // Mirrors the STEP 1 logic in processVenueTransmissions so the pre-pass
  // can bin people under the PARENT matrix (which may have a different bin
  // scheme than the child).
  int computeBinIndexForMatrix(const Person* person, const Venue* venue,
                               SubsetIndex subset_index,
                               const ContactMatrix* matrix,
                               int num_groups) const;

  // STEP 2c: pre-STEP-3 early exit check. Returns true iff every
  // transmission source is empty (no infectious bin, no positive fomite
  // lambda, no compartmental uptake potential) OR every bin has no
  // susceptibles. On true return the caller must still call clearAfterUse
  // on used_bins_. This helper only inspects state.
  bool venueHasNoTransmissionPossible(
      int num_bins_needed, const std::vector<int>& comp_uptake_modes,
      const std::vector<double>& lambda_fomite_by_mode, VenueId actual_venue_id,
      const CompartmentalModelManager* comp_model) const;

  // STEP 2d: look up the ParentAggregate (and its flat-bin contact matrix)
  // that this venue's susceptibles should read sibling-mixing FOI from. Returns
  // {nullptr, nullptr} if not applicable. Throws if the parent contact matrix
  // has more than one bin (V1 limitation).
  std::pair<const ParentAggregate*, const ContactMatrix*>
  getParentAggregateForVenue(const Venue* venue,
                             bool is_virtual_encounter) const;

  // STEP 2b: append per-sub-bin deposition entries to venue->fomite_history
  // for each fomite mode and compute the per-mode integrated lambda_fomite
  // contribution. Returns the lambda_fomite_by_mode vector
  // (size = num_fomite_modes; element = 0 when no infectiousness curve or
  // no deposition).
  std::vector<double> recordFomiteDepositionAndLambda(
      Venue* venue, int num_bins_needed, int num_fomite_modes,
      const std::vector<FomiteModeRef>& fomite_modes,
      const std::vector<int>& n_sub_per_mode, double current_time,
      double delta_hours);

  // STEP 1b: sort each bin's infectious_ids by person_id and apply the same
  // permutation to infectiousness_by_mode[m] for every mode. Ensures the
  // discrete_distribution index → person mapping is the same regardless of
  // which rank or in what order the venue was processed.
  void sortInfectiousByPersonId(int num_bins_needed, int num_modes);

  // STEP 1c: sort each bin's susceptible vector by person_id.
  void sortSusceptiblesByPersonId(int num_bins_needed);

  // STEP 2: pre-calculate per-mode cumulative weight arrays for infectious
  // bins. Sampling sites later read these via sampleFromCumulative. Bins with
  // all-zero weight per mode are left empty so callers can skip sampling.
  void buildCumulativeWeightsPerBin(int num_bins_needed, int num_modes);

  // STEP 3a: per-(mode, inf_bin) direct-contact source builder for one
  // susceptible bin. Walks every mode and every infectious bin, computes
  // omega * inf_total * susc_mult, and appends a SourceEntry{m, inf_bin}
  // for each positive weight. Mutates sources_buffer_, source_weights_buffer_,
  // and accumulates into total_lambda_eff. Routes through
  // lookupContactsForBinPair for the contacts lookup.
  void appendDirectContactSources(int susc_bin, int num_bins_needed,
                                  int num_modes, bool is_virtual_encounter,
                                  uint8_t encounter_type_id,
                                  uint8_t venue_type_id,
                                  const TransmissionParams& trans_params,
                                  double& total_lambda_eff);

  // STEP 3a.bis: append a single SIBLING-mixing source per mode for the
  // current susceptible bin. Reads from the parent ParentAggregate, subtracts
  // own-child contributions, and pushes one SourceEntry with the
  // SIBLING_INF_BIN_SENTINEL into sources_buffer_. Bumps
  // dbg_sample_susc_prints_ for verbose parent-mixing diagnostics. No-op if
  // parent_agg == nullptr.
  void appendSiblingMixingSources(int susc_bin, int num_modes,
                                  VenueId actual_venue_id, const Venue* venue,
                                  const ParentAggregate* parent_agg,
                                  const ContactMatrix* parent_flat_matrix,
                                  const TransmissionParams& trans_params,
                                  double& total_lambda_eff);

  // Verbose JUNE_DEBUG_PARENT_MIXING per-(susc_bin, mode) sibling-FOI dump.
  // Inlined out of appendSiblingMixingSources to keep that orchestrator
  // within the per-function LOC budget. dbg_sample_susc_prints_ is rate-
  // limited at 20 lines per tick.
  void logSiblingFOIDump(int m, int susc_bin, VenueId actual_venue_id,
                         const Venue* venue, const ParentAggregate& parent_agg,
                         double contacts, double parent_inf, double own_inf,
                         double sibling_inf, int parent_size, int own_size,
                         int sibling_size, double omega, double weighted);

  // STEP 3a: append per-fomite-mode sentinel sources (inf_bin = -1) to the
  // current susceptible bin. Multiplies lambda_fomite_by_mode[fm] by the
  // mode's susceptibility multiplier.
  void appendFomiteSources(int num_fomite_modes,
                           const std::vector<FomiteModeRef>& fomite_modes,
                           const std::vector<double>& lambda_fomite_by_mode,
                           const TransmissionParams& trans_params,
                           double& total_lambda_eff);

  // STEP 3a: append compartmental-uptake sentinel sources (inf_bin = -2) to
  // the current susceptible bin. Reads the plugin's coupling output buffer,
  // scales by the venue-type foi_scale, and applies each mode's susc_mult.
  void appendCompUptakeSources(VenueId actual_venue_id, uint8_t venue_type_id,
                               const std::vector<int>& comp_uptake_modes,
                               const CompartmentalModelManager* comp_model,
                               const TransmissionParams& trans_params,
                               double& total_lambda_eff);

  // Build cumulative weights over the parent's infector pool for one mode,
  // excluding entries whose origin child_venue_id matches actual_venue_id, then
  // weight-sample one entry with susc_rng. Returns the chosen infector's
  // PersonId or -1 if the pool is empty after exclusion. Uses the scratch
  // sibling_cum_buffer_ / sibling_pool_indices_buffer_.
  PersonId sampleSiblingInfector(int sampled_mode, int pbin,
                                 VenueId actual_venue_id,
                                 const ParentAggregate& parent_agg,
                                 SplitMix64& susc_rng);

  // STEP 3b infector dispatch: given a SourceEntry index into sources_buffer_,
  // resolve the InfectionSource enum + transmission mode index + infector
  // PersonId for one chosen susceptible. Handles all four sentinel cases:
  //   -2 → compartmental uptake (Person infector_id stays -1)
  //   -1 → fomite (infector_id stays -1)
  //   SIBLING_INF_BIN_SENTINEL → cross-child Person via parent_agg's per-bin
  //     pool, excluding the susceptible's own venue
  //   ≥0 → in-bin Person sampled from cumulative_by_mode[sampled_mode]
  // Reads + advances susc_rng. Mutates dbg_sibling_infections_ /
  // dbg_sample_infection_prints_ when applicable.
  PersonId sampleVenueInfector(int src_idx, int susc_bin,
                               VenueId actual_venue_id, const Venue* venue,
                               PersonId susceptible_id,
                               const ParentAggregate* parent_agg,
                               SplitMix64& susc_rng,
                               InfectionSource& infection_source_out,
                               uint8_t& transmission_mode_index_out);

  // Clear per-bin person-proportional state on every bin index touched this
  // call (per used_bins_), then reset used_bins_. Called at both early-exit
  // and end-of-call sites of processVenueTransmissions.
  void clearUsedBins(int num_modes);

  // Partial-presence dispatch. If runtime_group_allocator_ is wired AND the
  // venue's type is enabled in simulation_config_.partial_presence.enabled_
  // venue_type_mask, delegate to processPartialPresenceVenue and return its
  // new_infections. Otherwise returns std::nullopt to signal the caller
  // should proceed with the regular FOI path.
  std::optional<int> dispatchPartialPresenceIfApplicable(
      const std::vector<InteractionMember>& members, Venue* venue,
      VenueId actual_venue_id, double current_time, double delta_hours,
      std::unordered_set<PersonId>* active_infections,
      const std::unordered_set<PersonId>* visitor_ids,
      std::vector<PendingInfection>* pending_infections,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      uint8_t encounter_type_id, const CompartmentalModelManager* comp_model);

  // Bin every member, deterministically sort bins, build per-mode cumulative
  // weights, and record per-fomite-mode deposition + lambda. Bundles
  // STEP 1, 1b, 1c, 2, and 2b of processVenueTransmissions into one named
  // step. Returns the per-fomite-mode lambda vector.
  std::vector<double> binMembersAndPrepareBuffers(
      const std::vector<InteractionMember>& members, Venue* venue,
      const ContactMatrix& bin_structure, int num_bins_needed, int num_modes,
      int num_fomite_modes, const std::vector<FomiteModeRef>& fomite_modes,
      const std::vector<int>& n_sub_per_mode, double current_time,
      double delta_hours, uint8_t encounter_type_id,
      const std::string& venue_type, uint8_t venue_type_id,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data);

  // STEP 3 per-susceptible-bin orchestrator: build the four kinds of source
  // (direct contact, sibling mixing, fomite, comp uptake), apply regional-risk
  // multiplier, short-circuit if total_risk == 0, build the cumulative weights
  // once, and iterate susceptibles via processOneVenueSusceptible. Returns the
  // number of new infections produced.
  int processOneSuscBin(
      int susc_bin, int num_bins_needed, int num_modes, int num_fomite_modes,
      bool is_virtual_encounter, uint8_t encounter_type_id,
      uint8_t venue_type_id, Venue* venue, VenueId actual_venue_id,
      const std::vector<FomiteModeRef>& fomite_modes,
      const std::vector<int>& comp_uptake_modes,
      const std::vector<double>& lambda_fomite_by_mode,
      const TransmissionParams& trans_params, const ParentAggregate* parent_agg,
      const ContactMatrix* parent_flat_matrix,
      const CompartmentalModelManager* comp_model, double current_time,
      double delta_hours,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      std::unordered_set<PersonId>* active_infections,
      std::vector<PendingInfection>* pending_infections);

  // STEP 3b body for one susceptible: derive its deterministic RNG, Bernoulli-
  // draw against total_risk * susc.susceptibility, sample (mode, infector) +
  // symptom on success, and apply the infection. Returns true iff a new
  // infection was created.
  bool processOneVenueSusceptible(
      const SusceptibleMember& susc_mem, double total_risk, int susc_bin,
      bool have_source_dist, uint64_t time_bits, double current_time,
      VenueId actual_venue_id, Venue* venue, uint8_t venue_type_id,
      const ParentAggregate* parent_agg,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      std::unordered_set<PersonId>* active_infections,
      std::vector<PendingInfection>* pending_infections);

  // STEP 3b apply step for a single susceptible whose Bernoulli draw was
  // positive. Queues a PendingInfection for cross-rank visitors or creates
  // an Infection in place + logs it. Mirrors applyPartialPresenceInfection
  // but takes the susc_mem encounter_type_id and the venue-FOI InfectionSource.
  void applyVenueInfection(
      const SusceptibleMember& susc_mem, PersonId infector_id,
      InfectionSource infection_source, uint8_t transmission_mode_index,
      uint16_t infector_symptom_id, double current_time, uint8_t venue_type_id,
      VenueId actual_venue_id, uint64_t venue_key,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      std::unordered_set<PersonId>* active_infections,
      std::vector<PendingInfection>* pending_infections);

  // STEP 1 classification dispatch for one (member, person, visitor) tuple
  // already pinned to bin_index. Pushes susceptible into
  // bins_buffer_[bin_index] or routes to accumulate{Visitor,Local}* helpers
  // when infectious.
  void binMemberClassification(const InteractionMember& member, Person* person,
                               const VisitorInfo* visitor, int bin_index,
                               int num_modes, int num_fomite_modes,
                               const std::vector<FomiteModeRef>& fomite_modes,
                               const std::vector<int>& n_sub_per_mode,
                               double current_time, double delta_hours);

  // Resolve the matrix bin for a member: route through
  // computeBinIndexForMatrix, bump stats_.bin_lookups when matrix is non-null,
  // and clamp to [0, num_bins_needed). Emits a (rate-limited)
  // DEBUG_TRANSMISSION line if the clamp fired.
  int resolveMemberBinIndex(const InteractionMember& member,
                            const Person* person, Venue* venue,
                            const ContactMatrix& bin_structure,
                            int num_bins_needed, uint8_t encounter_type_id,
                            const std::string& venue_type,
                            uint8_t venue_type_id);

  // STEP 1 per-member body in processVenueTransmissions: resolve
  // person/visitor, compute bin_index (with DEBUG fallback log), update
  // used_bins_ + total_size, then dispatch infectiousness (visitor or local),
  // fomite deposition, or susceptible.
  void binOneMember(
      const InteractionMember& member, Venue* venue,
      const ContactMatrix& bin_structure, int num_bins_needed, int num_modes,
      int num_fomite_modes, const std::vector<FomiteModeRef>& fomite_modes,
      const std::vector<int>& n_sub_per_mode, double current_time,
      double delta_hours, uint8_t encounter_type_id,
      const std::string& venue_type, uint8_t venue_type_id,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data);

  // Append a cross-rank visitor's per-mode infectiousness (pre-computed on
  // the sending rank) into bins_buffer_[bin_index] and accumulate its
  // per-sub-bin fomite deposition. Caller has already confirmed
  // visitor->is_infectious. Uses im_scratch_buffer_ as scratch.
  void accumulateVisitorInfectiousnessAndFomite(
      const VisitorInfo* visitor, PersonId pid, int bin_index, int num_modes,
      int num_fomite_modes, const std::vector<FomiteModeRef>& fomite_modes,
      const std::vector<int>& n_sub_per_mode, double delta_hours);

  // Append a local infectious person's per-mode integrated infectiousness
  // into bins_buffer_[bin_index]. Caller has already confirmed
  // person->infection && person->infection->isInfectious(current_time). Uses
  // im_scratch_buffer_ as scratch.
  void accumulateLocalInfectiousness(const Person* person, PersonId pid,
                                     int bin_index, int num_modes,
                                     double current_time, double delta_hours);

  // Accumulate a local infected person's fomite deposition into
  // bins_buffer_[bin_index].total_fomite_deposition_sub. Only called when
  // person->infection != nullptr and num_fomite_modes > 0.
  void accumulateLocalFomiteDeposition(
      const Person* person, int bin_index, int num_fomite_modes,
      const std::vector<FomiteModeRef>& fomite_modes,
      const std::vector<int>& n_sub_per_mode, double current_time,
      double delta_hours);

  // Resize bins_buffer_ to at least num_bins_needed entries, then ensure
  // every active bin has correctly-sized per-mode vectors and pre-sized
  // fomite sub-bins for this slot's delta_hours / sub_bin_time. Re-initialises
  // bins whose infectiousness_by_mode size mismatches num_modes (which happens
  // for freshly-grown bins or when num_modes changed).
  void prepareBinsBuffer(int num_bins_needed, int num_modes,
                         int num_fomite_modes,
                         const std::vector<int>& n_sub_per_mode);

  // Walk the disease's transmission modes and split them into a flat list of
  // FomiteModeRefs and a separate list of CompartmentalUptake mode indices.
  // Also compute n_sub_per_mode from each fomite mode's sub_bin_time and
  // delta_hours.
  void collectFomiteAndCompUptakeModes(
      double delta_hours, std::vector<FomiteModeRef>& fomite_modes_out,
      std::vector<int>& comp_uptake_modes_out,
      std::vector<int>& n_sub_per_mode_out) const;

  // Resolve the venue type id, human-readable type label, and contact matrix
  // for a venue group. For VIRTUAL encounters (actual_venue_id < 0) the
  // encounter_type_id is used as the key (via virtual_matrix_names /
  // getVirtualMatrix). For PHYSICAL venues the venue's own type is used,
  // even if some members are encounter participants. Bumps
  // stats_.matrix_lookups for physical venues. Out params get sensible defaults
  // if no matrix found.
  void resolveVenueTypeAndMatrix(Venue* venue, VenueId actual_venue_id,
                                 uint8_t encounter_type_id,
                                 std::string& venue_type_out,
                                 uint8_t& venue_type_id_out,
                                 const ContactMatrix*& matrix_out);

  // Process transmissions within a single venue/subset
  int processVenueTransmissions(
      const std::vector<InteractionMember>& members, Venue* venue,
      VenueId actual_venue_id, double current_time, double delta_hours,
      std::unordered_set<PersonId>* active_infections = nullptr,
      const std::unordered_set<PersonId>* visitor_ids = nullptr,
      std::vector<PendingInfection>* pending_infections = nullptr,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data = nullptr,
      uint8_t encounter_type_id = kDefaultEncounterTypeId,
      const CompartmentalModelManager* comp_model = nullptr);

  // Return person_ids of susceptibles sorted ascending; the partial-presence
  // post-pass iterates susceptibles in this order so per-call work is
  // deterministic across MPI rank counts.
  std::vector<PersonId> orderSusceptibles(
      const std::unordered_map<PersonId, double>& susc_lambda) const;

  // Weight-sample one (mode, infector) from accumulated AccumSource entries.
  // Sorts in place by (mode, infector) for deterministic order, builds the
  // cumulative weights, and draws one sample with the given RNG. Returns
  // mode=0, infector=-1 when the source list is empty / all-zero-weight.
  std::pair<int, PersonId> sampleInfectorFromAccumSources(
      std::vector<PartialPresenceAccumSource>& srcs, SplitMix64& rng) const;

  // Look up the infector's current symptom id. For local persons reads from
  // Infection::getTrajectory(); for cross-rank visitors reads from
  // VisitorInfo::symptom_id. Returns 0 if infector_id is negative or neither
  // a local infection nor a visitor record exists.
  uint16_t resolveInfectorSymptomId(
      PersonId infector_id, double current_time,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data) const;

  // Compute susceptibility for either a local Person or a cross-rank
  // VisitorInfo. Returns 0.0 if both are null. Local persons go through
  // Person::getSusceptibility(current_time, disease_name); visitors use
  // (1.0 - immunity_level).
  double computeMemberSusceptibility(const Person* person,
                                     const VisitorInfo* visitor,
                                     double current_time) const;

  // Run the per-susceptible Bernoulli draw + (if infected) infector sampling
  // + apply step for one susceptible id in the partial-presence post-pass.
  // Returns true iff a new infection was created (counted toward
  // new_infections). All other susceptible filtering (zero λ, zero
  // susceptibility, missed roll) returns false.
  bool processOnePartialSusceptible(
      PersonId susc_id, const std::unordered_map<PersonId, double>& susc_lambda,
      std::unordered_map<PersonId, std::vector<PartialPresenceAccumSource>>&
          susc_sources,
      double current_time, Venue* venue, uint8_t venue_type_id,
      VenueId actual_venue_id,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      std::unordered_set<PersonId>* active_infections,
      std::vector<PendingInfection>* pending_infections, uint64_t time_bits,
      uint64_t venue_key);

  // Apply a single partial-presence infection: either queue a
  // Note down that this line would infect this rider, without acting on it.
  // The winner among a rider's legs is decided later, in
  // resolvePartialPresenceInfections.
  void recordPartialPresenceCandidate(PersonId susc_id, PersonId infector_id,
                                      uint8_t transmission_mode_index,
                                      uint16_t infector_symptom_id,
                                      double current_time,
                                      uint8_t venue_type_id,
                                      VenueId actual_venue_id);

  // What each line said it would do to its riders this slot, before any of it
  // is applied. Cleared at the start of every partial-presence pass.
  std::vector<PendingInfection> pp_candidates_;

  // Sibling of processVenueTransmissions for partial-presence venues
  // (transport_line, etc.; anything declared in
  // SimulationConfig::partial_presence). Per-rider group assignments come
  // from the runtime_group_allocator_; per-rider effective presence windows
  // come from the membership_metadata side-table + presence_window helper.
  //
  // v1 scope (assumed and enforced; throws on violation):
  //   - Physical venue (actual_venue_id >= 0); not a virtual encounter venue.
  //   - No parent venue (transport lines have none in current MAY output).
  //   - No coordinated encounter participants (encounter_type_id ==
  //   kDefaultEncounterTypeId).
  //   - Direct-contact FOI only; no fomite / compartmental uptake on
  //     partial-presence venue types in v1.
  // Violations throw with a descriptive error rather than silently
  // falling back, per the project's no-silent-fallbacks rule.
  int processPartialPresenceVenue(
      const std::vector<InteractionMember>& members, Venue* venue,
      VenueId actual_venue_id, double current_time, double delta_hours,
      std::unordered_set<PersonId>* active_infections,
      const std::unordered_set<PersonId>* visitor_ids,
      std::vector<PendingInfection>* pending_infections,
      const std::unordered_map<PersonId, VisitorInfo>* visitor_data,
      uint8_t encounter_type_id, const CompartmentalModelManager* comp_model);

  PerformanceStats stats_;

  // Optimization buffers (reused across calls to avoid allocation)
  std::vector<PersonLocation> active_locations_buffer_;
  std::vector<InteractionMember> group_members_buffer_;

  // BinGroup reuse buffer
  std::vector<BinGroup> bins_buffer_;
  // Tracks which bin indices were actually used (for selective clearing)
  std::vector<int> used_bins_;

  // Source/weight buffers for transmission sampling
  std::vector<SourceEntry> sources_buffer_;
  std::vector<double> source_weights_buffer_;
  // Cumulative-weight scratch for source sampling (rebuilt per susc_bin,
  // sampled once per susceptible). Avoids per-bin std::discrete_distribution
  // allocation; see sampleFromCumulative in utils/random.h.
  std::vector<double> source_cumulative_buffer_;

  // Scratch buffers for sibling-infector two-stage sampling. Reused across
  // infections so we don't reallocate on every sibling-attributed event.
  std::vector<double> sibling_cum_buffer_;
  std::vector<size_t> sibling_pool_indices_buffer_;

  // Per-mode infectiousness scratch buffer
  std::vector<double> im_scratch_buffer_;

  // Cached uniform distribution for transmission rolls
  std::uniform_real_distribution<double> uniform_dist_{0.0, 1.0};

  // Base seed for deterministic per-entity RNG (MPI reproducibility)
  uint64_t base_seed_ = 0;

  // Per-tick parent-venue aggregates for cross-child-venue ("sibling")
  // mixing. Cleared and rebuilt at the start of each processTransmissions
  // call. Rank-local by construction: all child venues of a parent share
  // its MGU, so no MPI sync is required.
  std::unordered_map<VenueId, ParentAggregate> parent_aggregates_;

  // Wired by Simulator after construction. Null = no partial-presence
  // bucketing (sub-interval FOI is skipped; processVenueTransmissions
  // gate stays a no-op).
  const RuntimeGroupAllocator* runtime_group_allocator_ = nullptr;

  // Cached env-flag enabling verbose parent-mixing debug prints
  // (JUNE_DEBUG_PARENT_MIXING=1). Read once at construction.
  bool debug_parent_mixing_ = false;
  // Per-tick debug counters so we can summarise at the end of each tick
  // instead of spamming once per encounter.
  mutable int dbg_sibling_infections_ = 0;
  mutable int dbg_sample_susc_prints_ = 0;
  mutable int dbg_sample_infection_prints_ = 0;
};

}  // namespace june
