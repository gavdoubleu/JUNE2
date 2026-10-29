#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "doctest.h"
#include "simulation/checkpoint_manifest/checkpoint_manifest.h"

using namespace june;

namespace {

const unsigned int kSeed = 42;

std::filesystem::path freshCheckpointDir() {
  const auto dir =
      std::filesystem::temp_directory_path() / "june_checkpoint_manifest";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  return dir;
}

void writeHandManifest(const std::filesystem::path& dir, int format_version) {
  std::ofstream(dir / "manifest.yaml")
      << "format_version: " << format_version << "\n"
      << "completed_day: 5\n"
      << "effective_random_seed: " << kSeed << "\n";
}

std::string restoreMessage(const std::filesystem::path& dir,
                           unsigned int config_seed) {
  try {
    readAndValidateManifest(dir, config_seed);
  } catch (const std::runtime_error& error) {
    return error.what();
  }
  FAIL("readAndValidateManifest accepted the manifest");
  return "";
}

}  // namespace

TEST_CASE("A written manifest reads back its completed day") {
  const auto dir = freshCheckpointDir();
  writeManifest(dir, 7, "2021-01-08", 2, 7.0, kSeed);

  CHECK(readAndValidateManifest(dir, kSeed) == 7);
}

TEST_CASE("A manifest of an earlier format version is refused") {
  const auto dir = freshCheckpointDir();
  writeHandManifest(dir, kCheckpointFormatVersion - 1);

  const std::string message = restoreMessage(dir, kSeed);
  CHECK(message.find("format_version " +
                     std::to_string(kCheckpointFormatVersion - 1)) !=
        std::string::npos);
  CHECK(message.find("reads only " +
                     std::to_string(kCheckpointFormatVersion)) !=
        std::string::npos);
}

TEST_CASE("A manifest whose seed differs from the config is refused") {
  const auto dir = freshCheckpointDir();
  writeHandManifest(dir, kCheckpointFormatVersion);

  CHECK(restoreMessage(dir, kSeed + 1).find("effective_random_seed mismatch") !=
        std::string::npos);
}
