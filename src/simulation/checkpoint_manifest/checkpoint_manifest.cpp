#include "simulation/checkpoint_manifest/checkpoint_manifest.h"

#include <yaml-cpp/yaml.h>

#include <fstream>
#include <stdexcept>

namespace june {

namespace fs = std::filesystem;

// world_path / world_sha256 are populated in P4 (restore validation);
// intentionally empty here rather than a misleading placeholder.
void writeManifest(const fs::path& dir, int completed_day,
                   const std::string& date_iso, int nranks,
                   double current_simulation_time, unsigned int random_seed) {
  YAML::Emitter m;
  m << YAML::BeginMap;
  m << YAML::Key << "format_version" << YAML::Value << kCheckpointFormatVersion;
  m << YAML::Key << "completed_day" << YAML::Value << completed_day;
  m << YAML::Key << "date" << YAML::Value << date_iso;
  m << YAML::Key << "current_simulation_time" << YAML::Value
    << current_simulation_time;
  m << YAML::Key << "num_ranks" << YAML::Value << nranks;
  m << YAML::Key << "effective_random_seed" << YAML::Value
    << static_cast<unsigned long long>(random_seed);
  m << YAML::Key << "world_path" << YAML::Value << "";
  m << YAML::Key << "world_sha256" << YAML::Value << "";
  m << YAML::Key << "shard_index" << YAML::Value << "shard_index.yaml";
  m << YAML::Key << "state_file" << YAML::Value << "state.h5";
  m << YAML::EndMap;
  std::ofstream(dir / "manifest.yaml") << m.c_str() << "\n";
}

int readAndValidateManifest(const fs::path& cp, unsigned int config_seed) {
  if (!fs::exists(cp / "manifest.yaml"))
    throw std::runtime_error(
        "restoreFromCheckpoint: no manifest.yaml (incomplete/invalid "
        "checkpoint): " +
        cp.string());
  YAML::Node man = YAML::LoadFile((cp / "manifest.yaml").string());
  const int format_version = man["format_version"].as<int>();
  if (format_version != kCheckpointFormatVersion) {
    throw std::runtime_error(
        "restoreFromCheckpoint: unsupported checkpoint format_version " +
        std::to_string(format_version) + " (this build reads only " +
        std::to_string(kCheckpointFormatVersion) +
        "; restart from day 0): " + cp.string());
  }
  int completed_day = man["completed_day"].as<int>();
  unsigned int cp_seed = man["effective_random_seed"].as<unsigned int>();
  if (cp_seed != config_seed) {
    throw std::runtime_error(
        "restoreFromCheckpoint: effective_random_seed mismatch (checkpoint=" +
        std::to_string(cp_seed) + ", config=" + std::to_string(config_seed) +
        "). Resume with the checkpoint's seed (no silent override).");
  }
  return completed_day;
}

}  // namespace june
