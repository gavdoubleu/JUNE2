#include <H5Cpp.h>

#include <filesystem>
#include <string>
#include <vector>

#include "doctest.h"
#include "simulation/checkpoint_masks/policy_mask_io.h"

using namespace june;

namespace {

const std::string kMaskDataset = "/applicable_symptom_policy_mask";

std::filesystem::path maskFilePath() {
  return std::filesystem::temp_directory_path() / "june_policy_mask_io.h5";
}

}  // namespace

TEST_CASE("Reading masks whose row count differs from the population throws") {
  const std::vector<PolicyMask> masks(3);
  {
    H5::H5File file(maskFilePath().string(), H5F_ACC_TRUNC);
    writeMasks(file, kMaskDataset, masks, 0);
  }
  H5::H5File file(maskFilePath().string(), H5F_ACC_RDONLY);

  std::string message;
  try {
    readMasks(file, kMaskDataset, 4);
    FAIL("readMasks accepted 3 rows for a population of 4");
  } catch (const std::runtime_error& error) {
    message = error.what();
  }
  CHECK(message.find(kMaskDataset) != std::string::npos);
  CHECK(message.find("3") != std::string::npos);
  CHECK(message.find("4") != std::string::npos);
}

TEST_CASE("Masks survive a write and read through a checkpoint file") {
  std::vector<PolicyMask> masks(4);
  masks[1].set(0);
  masks[2].set(63).set(64);
  masks[3].set(kMaxPoliciesPerKind - 1);
  {
    H5::H5File file(maskFilePath().string(), H5F_ACC_TRUNC);
    writeMasks(file, kMaskDataset, masks, 4);
  }
  H5::H5File file(maskFilePath().string(), H5F_ACC_RDONLY);

  CHECK(readMasks(file, kMaskDataset, masks.size()) == masks);
}
