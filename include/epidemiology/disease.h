#pragma once

#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <variant>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "infectiousness_curves.h"
#include "utils/deterministic_rng.h"
#include "utils/filtering.h"

namespace june {

// =============================================================================
// Symptom Tag - Represents a disease status/stage (loaded from YAML)
// =============================================================================

struct SymptomTag {
  std::string name;
  int value;        // Severity level (-2 = recovered, -1 = healthy, 0+ =
                    // progressively worse)
  uint16_t id = 0;  // Runtime unique ID for fast indexing

  bool operator<(const SymptomTag& other) const { return value < other.value; }

  bool operator==(const SymptomTag& other) const { return name == other.name; }
};

// =============================================================================
// Disease Stage Settings (loaded from YAML)
// =============================================================================

struct DiseaseStageSettings {
  std::string default_lowest_stage;
  std::string max_mild_symptom_tag;

  std::vector<std::string> stay_at_home_stages;
  std::vector<std::string> severe_symptoms_stay_at_home_stages;
  std::vector<std::string> hospitalised_stages;
  std::vector<std::string> intensive_care_stages;
  std::vector<std::string> fatality_stages;
  std::vector<std::string> recovered_stages;
};

// =============================================================================
// Distribution Parameters for Trajectory Timing
// =============================================================================

struct DistributionParams {
  std::string
      type;  // "constant", "normal", "lognormal", "beta", "exponweib", "gamma"
  std::map<std::string, double> params;  // Flexible parameter storage (mean,
                                         // std, a, b, loc, scale, s, c, etc.)
};

// =============================================================================
// Trajectory Stage - A single stage in a disease progression trajectory
// =============================================================================

struct TrajectoryStage {
  std::string symptom_tag;
  DistributionParams completion_time;
};

// =============================================================================
// Trajectory Definition - A complete disease progression path
// =============================================================================

struct TrajectoryDefinition {
  std::string description;
  std::string
      selection_key;  // Key to match in outcome rates (e.g., "hospital")
  std::optional<double>
      probability;  // Fixed selection weight; overrides CSV lookup when set
  double severity =
      0.0;  // 0.0 (mild) to 1.0 (dead); used by vaccine efficacy logic
  double infectiousness_factor =
      1.0;  // Scales the entire infectiousness profile for this trajectory
            // (e.g. 0.5 for asymptomatic cases). Applied once at infection

  std::vector<TrajectoryStage> stages;
  std::optional<std::string>
      start_stage;  // If set, trajectory enters at the first occurrence
                    // of this tag in stages rather than at index 0
};

// =============================================================================
// Filter-Based Outcome Rates
// =============================================================================

// A single demographic row: filter criteria + outcome probabilities.
struct OutcomeRow {
  std::vector<SelectionCriterion> criteria;     // empty = matches all persons
  std::map<std::string, double> probabilities;  // selection_key -> probability
};

struct OutcomeRates {
  std::vector<OutcomeRow> rows;

  // Return the probability for `outcome` for the first row that matches
  // `person`. `ctx` provides infection-specific context (infector symptom,
  // transmission mode). Returns 0.0 if no row matches or outcome is not present
  // in the matched row.
  double getRate(const Person& person, const WorldState* world,
                 const std::string& outcome,
                 const InfectionContext& ctx = {}) const;

  // Resolve property-name strings to integer codes after WorldState is built.
  void resolve(const WorldState& world);
};

// =============================================================================
// Transmission Parameters
// =============================================================================

enum class InfectiousnessMode {
  TRAJECTORY_DRIVEN,  // Default: based on time since infection
  STAGE_DRIVEN        // Based on time since current symptom started
};

// Configuration for a fomite transmission mode.
// The fomite_curve describes how the infectiousness of deposited material
// evolves as a function of time since deposition (same curve interface as
// person infectiousness). Deposition from infectious people depends on their
// current symptom stage via deposition_by_symptom.
struct FomiteConfig {
  int mode_index = -1;        // Index in TransmissionParams::modes
  double max_age = 14.0;      // Days: prune deposition history older than this
  double sub_bin_time = 0.0;  // Hours: target sub-bin width (0 = single bin)

  // Infectiousness of deposited material as a function of age (days since
  // deposit). Any InfectiousnessCurve type. Evaluated once per history entry
  // per time step.
  std::shared_ptr<InfectiousnessCurve> infectiousness_curve;

  // Deposition rate per symptom: [symptom_id] -> curve(time_in_stage) =
  // deposit/hour. Same InfectiousnessCurve format as stage_curves. nullptr = no
  // deposition.
  std::vector<std::shared_ptr<InfectiousnessCurve>> deposition_by_symptom;
};

// Configuration for a compartmental-uptake transmission mode.
struct CompartmentalUptakeConfig {
  int mode_index = -1;
};

// Configuration for a compartmental-deposition transmission mode.
// Per-symptom deposition rate curves, same layout as
// FomiteConfig::deposition_by_symptom.
struct CompartmentalDepositionConfig {
  int mode_index = -1;
  // [symptom_id] → deposition rate curve (nullptr = no contribution).
  std::vector<std::shared_ptr<InfectiousnessCurve>> deposition_by_symptom;
};

// =============================================================================
// Unified transmission mode. Replaces parallel is_fomite_mode /
// is_compartmental_uptake_mode / is_compartmental_deposition_mode vectors.
// =============================================================================

enum class TransmissionModeType {
  Standard,
  Fomite,
  CompartmentalUptake,
  CompartmentalDeposition,
};

struct TransmissionMode {
  std::string name;
  TransmissionModeType type = TransmissionModeType::Standard;
  double mode_transmissibility_multiplier = 1.0;
  // Per-symptom infectiousness curves for STAGE_DRIVEN Standard modes.
  // Empty for Fomite, CompartmentalUptake, CompartmentalDeposition modes.
  std::vector<std::shared_ptr<InfectiousnessCurve>> symptom_curves;
  // Type-specific configuration; monostate for Standard modes.
  std::variant<std::monostate, FomiteConfig, CompartmentalUptakeConfig,
               CompartmentalDepositionConfig>
      config;
};

struct NaturalImmunityParams {
  double level = 0.95;
  double waning_rate = 0.001;
};

struct TransmissionParams {
  InfectiousnessMode mode = InfectiousnessMode::TRAJECTORY_DRIVEN;
  std::string type;  // "gamma", "normal", etc. (used for TRAJECTORY_DRIVEN)

  // Infectiousness over time parameters (used for TRAJECTORY_DRIVEN)
  DistributionParams max_infectiousness;
  DistributionParams shape;
  DistributionParams rate;
  DistributionParams shift;

  // Stage-specific independent curves (used for STAGE_DRIVEN)
  std::map<std::string, std::shared_ptr<InfectiousnessCurve>> stage_curves;
  std::vector<std::shared_ptr<InfectiousnessCurve>> symptom_id_curves;

  // Natural immunity parameters
  NaturalImmunityParams natural_immunity;

  // Ordered transmission modes. Single-mode configs have one entry named
  // "default". Matched to ContactMatrixConfig::mode_names by name, not by
  // index, via ContactMatrixConfig::finalizeDiseaseModeAlignment.
  // Replaces the former parallel vectors: mode_names,
  // mode_transmissibility_multipliers, is_fomite_mode,
  // is_compartmental_uptake_mode, is_compartmental_deposition_mode,
  // fomite_configs, compartmental_uptake_configs,
  // compartmental_deposition_configs, and mode_symptom_curves.
  std::vector<TransmissionMode> modes;
};

// =============================================================================
// Disease Class - Fully configurable disease model
// =============================================================================

class Disease {
 public:
  Disease(const std::string& name, const std::vector<SymptomTag>& symptom_tags,
          const DiseaseStageSettings& stage_settings,
          const std::vector<TrajectoryDefinition>& trajectories,
          const OutcomeRates& outcome_rates,
          const TransmissionParams& transmission_params);

  // Getters
  const std::string& getName() const { return name_; }
  const std::vector<SymptomTag>& getSymptomTags() const {
    return symptom_tags_;
  }
  const DiseaseStageSettings& getStageSettings() const {
    return stage_settings_;
  }
  const TransmissionParams& getTransmissionParams() const {
    return transmission_params_;
  }
  const std::vector<TrajectoryDefinition>& getTrajectories() const {
    return trajectories_;
  }
  const OutcomeRates& getOutcomeRates() const { return outcome_rates_; }

  // Find symptom tag by name
  const SymptomTag* findSymptomTag(const std::string& name) const;

  // Check if a symptom tag is in a specific category
  bool isInCategory(const std::string& symptom_name,
                    const std::vector<std::string>& category) const;
  bool isFatalStage(const std::string& symptom_name) const;
  bool isRecoveredStage(const std::string& symptom_name) const;
  bool isHospitalisedStage(const std::string& symptom_name) const;
  bool isICUStage(const std::string& symptom_name) const;
  bool isInfectiousStage(const std::string& symptom_name) const;

  // Resolve outcome rate criteria after WorldState is built.
  void resolve(const WorldState& world) { outcome_rates_.resolve(world); }

  // Fast lookup
  uint16_t getSymptomId(const std::string& name) const;
  const std::string& getSymptomName(uint16_t id) const;
  const std::vector<std::string>& getSymptomNames() const {
    return id_to_name_;
  }

  /// Return the name of transmission mode at index. Empty string if out of
  /// range.
  const std::string& getModeName(uint8_t index) const;
  int numModes() const;

  /// Evaluate stage-driven infectiousness for a given mode, symptom, and
  /// time-in-stage. Returns 0.0 for TRAJECTORY_DRIVEN diseases.
  double evaluateStageDrivenInfectiousness(int mode_index, uint16_t symptom_id,
                                           float time_in_stage) const;

  /// Integrate stage-driven infectiousness over [t_in_stage_start,
  /// t_in_stage_end] (both in days). Returns 24 * integral (in hours).
  double integrateStageDrivenInfectiousness(int mode_index, uint16_t symptom_id,
                                            double t_in_stage_start,
                                            double t_in_stage_end) const;

  /// Evaluate fomite deposition rate for a given fomite mode, symptom, and
  /// time-in-stage. Used by the host rank to reconstruct per-mode deposition
  /// for cross-rank visitors.
  double evaluateFomiteDeposition(int fomite_mode_index, uint16_t symptom_id,
                                  double time_in_stage) const;

  /// Returns 24 * ∫_{t_start}^{t_end} dep_rate(τ) dτ (integral in hours).
  double integrateFomiteDeposition(int fomite_mode_index, uint16_t symptom_id,
                                   double t_in_stage_start,
                                   double t_in_stage_end) const;

 private:
  std::string name_;
  std::vector<SymptomTag> symptom_tags_;
  std::map<std::string, SymptomTag> symptom_tag_map_;  // For quick lookup
  DiseaseStageSettings stage_settings_;
  TransmissionParams transmission_params_;
  std::vector<TrajectoryDefinition> trajectories_;
  OutcomeRates outcome_rates_;

  std::vector<std::string> id_to_name_;  // Maps runtime ID to symptom name
};

// =============================================================================
// Infection Trajectory
// =============================================================================

struct InfectionTrajectory {
  double infection_time;
  double infectiousness_factor = 1.0;

  std::vector<std::pair<double, uint16_t>> transitions;  // (time, symptom_id)

  // Get current symptom ID at given time
  uint16_t getCurrentSymptomId(double current_time) const {
    for (auto it = transitions.rbegin(); it != transitions.rend(); ++it) {
      if (current_time >= it->first) {
        return it->second;
      }
    }
    return transitions.empty() ? 0 : transitions[0].second;
  }

  // Get time until next transition
  std::optional<double> getNextTransitionTime(double current_time) const {
    for (const auto& [time, symptom] : transitions) {
      if (time > current_time) {
        return time;
      }
    }
    return std::nullopt;
  }
};

// =============================================================================
// Infection Class
// =============================================================================

class Infection {
 public:
  Infection(const Disease* disease, double infection_time,
            const Person* person,  // Pass person for vaccine context
            unsigned int random_seed, const WorldState* world = nullptr,
            const std::string& venue_type = "", int venue_id = -1,
            float severity_factor = 1.0f, uint16_t infector_symptom_id = 0,
            const std::string& trajectory_key_override = "",
            const std::string& start_symptom_override = "",
            uint8_t transmission_mode_index = 0);

  // Getters
  const Disease* getDisease() const { return disease_; }
  std::string getCurrentSymptom(double current_time) const;
  double getInfectionTime() const { return infection_time_; }
  const InfectionTrajectory& getTrajectory() const { return trajectory_; }

  // --- Checkpoint serialization accessors (read-only) ---
  // Expose the privately-sampled transmission params + stage cache so a
  // checkpoint can reconstruct an Infection without re-sampling from RNG.
  // (Restore path added in P4.)
  double ckptMaxInfectiousness() const { return max_infectiousness_; }
  double ckptTransmissionShape() const { return transmission_shape_; }
  double ckptTransmissionRate() const { return transmission_rate_; }
  double ckptTransmissionShift() const { return transmission_shift_; }
  double ckptLastCheckedTime() const { return last_checked_time_; }
  uint16_t ckptCachedSymptomId() const { return cached_symptom_id_; }
  double ckptCachedSymptomStartTime() const {
    return cached_symptom_start_time_;
  }

  // Reconstruct an Infection from checkpoint state WITHOUT re-sampling any
  // RNG. Every behaviour-defining field is restored verbatim so a resumed
  // run is bit-identical.
  static std::unique_ptr<Infection> fromCheckpoint(
      const Disease* disease, double infection_time,
      const InfectionTrajectory& trajectory, double max_infectiousness,
      double transmission_shape, double transmission_rate,
      double transmission_shift, double last_checked_time,
      uint16_t cached_symptom_id, double cached_symptom_start_time);

  // Check status
  bool isInfectious(double current_time) const;
  bool isSymptomatic(double current_time) const;
  bool isRecovered(double current_time) const;
  bool isDead(double current_time) const;

  // Get infectiousness at a given time (based on time since infection)
  double getInfectiousness(double current_time) const;

  /// Per-mode infectiousness. Falls back to mode 0 behaviour if mode_index is
  /// out of range.
  double getInfectiousness(int mode_index, double current_time) const;

  /// Returns 24 * ∫_{t0}^{t1} I_mode(τ) dτ (integral in hours, for drop-in
  /// compatibility with the existing delta_hours * I_point accumulation).
  double getIntegratedInfectiousness(int mode_index, double t0,
                                     double t1) const;

  /// Fomite deposition rate for a given fomite mode at current_time.
  /// Evaluates modes[fomite_mode_index]
  /// FomiteConfig::deposition_by_symptom[symptom] at time_in_stage. Returns 0
  /// if no deposition curve is configured for the current symptom/mode.
  double getFomiteDepositRate(int fomite_mode_index, double current_time) const;

  /// Returns 24 * ∫_{t0}^{t1} dep_rate(τ) dτ (integral in hours).
  double getIntegratedFomiteDeposition(int fomite_mode_index, double t0,
                                       double t1) const;

  // Get time until next status change
  std::optional<double> getNextTransitionTime(double current_time) const;

 private:
  // Tag ctor for checkpoint restore: sets only disease_, skips all sampling.
  struct RestoreTag {};
  Infection(const Disease* disease, RestoreTag) : disease_(disease) {}

  // Resolve current symptom id + stage start time at `lookup_time`, using
  // the mutable cache when `lookup_time == last_checked_time_` and otherwise
  // scanning `trajectory_.transitions` and refreshing the cache.
  void cacheCurrentSymptom(double lookup_time, uint16_t& symptom_id,
                           double& stage_start_time) const;

  // Append transitions to `traj` by walking `traj_def.stages` from the first
  // stage whose symptom_tag matches `start_target` (or index 0 if empty / not
  // found), sampling each stage's completion_time from `rng` to advance time.
  void buildTransitionsFromStages(InfectionTrajectory& traj,
                                  const TrajectoryDefinition& traj_def,
                                  const std::string& start_target,
                                  SplitMix64& rng);

  // Build a trajectory directly from the trajectory whose selection_key
  // matches `key`, skipping rate-based selection. Returns empty if no
  // trajectory has that key.
  std::optional<InfectionTrajectory> tryBuildForcedTrajectory(
      const std::string& key, SplitMix64& rng);

  // Step 1 of rate-based selection: gather a raw rate per trajectory and
  // their sum. Trajectories with `probability` use that as a fixed weight;
  // otherwise the rate is looked up from outcome_rates via selection_key.
  // The returned rates are renormalised to sum to 1 if total > 1.0001.
  std::pair<std::vector<double>, double> gatherTrajectoryRates(
      const Person& person, const WorldState* world,
      const InfectionContext& ctx) const;

  // Step 2 of rate-based selection: shift probability mass from severe
  // trajectories to the lowest-severity ("safe") trajectory in proportion
  // to the person's symptom-reducing vaccine efficacy at `infection_time`.
  // No-op if the person is unvaccinated or efficacy is zero. Preserves the
  // sum (total_rate) since mass is redistributed.
  void applyVaccineEfficacyShift(std::vector<double>& rates,
                                 const Person& person, double infection_time,
                                 double total_rate) const;

  // Step 3 of rate-based selection: draw a uniform [0,1) sample and return
  // the index of the first trajectory whose cumulative mass crosses it.
  // Returns 0 if `total_rate <= 0` or if no cumulative bucket is hit (due
  // to normalisation / precision drift).
  static int sampleTrajectoryIndex(const std::vector<double>& rates,
                                   double total_rate, SplitMix64& rng);

  const Disease* disease_;
  double infection_time_;
  InfectionTrajectory trajectory_;

  // Transmission parameters (sampled at infection time)
  double max_infectiousness_ = 1.0;
  double transmission_shape_ = 1.0;
  double transmission_rate_ = 1.0;
  double transmission_shift_ = 0.0;

  // Cache for stage-driven lookups
  mutable double last_checked_time_ = -1.0;
  mutable uint16_t cached_symptom_id_ = 0;
  mutable double cached_symptom_start_time_ = -1.0;

  InfectionTrajectory generateTrajectoryFromRates(
      SplitMix64& rng, const Person* person, const WorldState* world,
      const std::string& venue_type = "", int venue_id = -1,
      float severity_factor = 1.0f, uint16_t infector_symptom_id = 0,
      const std::string& trajectory_key_override = "",
      const std::string& start_symptom_override = "",
      uint8_t transmission_mode_index = 0);

  // Sample transmission parameters from disease config
  void sampleTransmissionParameters(SplitMix64& rng);

  // Helper to sample from distribution
  double sampleFromDistribution(const DistributionParams& dist,
                                SplitMix64& rng);

  // Evaluate transmission profile (e.g. gamma PDF)
  double evaluateTransmissionProfile(double x) const;
};

}  // namespace june
