#include <H5Cpp.h>

#include <cstdint>
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

TEST_CASE("Reading a dataset that is not N x 2 words throws") {
  const std::vector<uint32_t> one_word_per_person(3, 0);
  {
    H5::H5File file(maskFilePath().string(), H5F_ACC_TRUNC);
    hsize_t dims[1] = {one_word_per_person.size()};
    H5::DataSet dataset = file.createDataSet(
        kMaskDataset, H5::PredType::NATIVE_UINT32, H5::DataSpace(1, dims));
    dataset.write(one_word_per_person.data(), H5::PredType::NATIVE_UINT32);
  }
  H5::H5File file(maskFilePath().string(), H5F_ACC_RDONLY);

  std::string message;
  try {
    readMasks(file, kMaskDataset, one_word_per_person.size());
    FAIL("readMasks accepted an N x 1 dataset");
  } catch (const std::runtime_error& error) {
    message = error.what();
  }
  CHECK(message.find(kMaskDataset) != std::string::npos);
  CHECK(message.find("corrupt checkpoint") != std::string::npos);
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

TEST_CASE("Reading an absent mask dataset throws, even for an empty population") {
  { H5::H5File file(maskFilePath().string(), H5F_ACC_TRUNC); }
  H5::H5File file(maskFilePath().string(), H5F_ACC_RDONLY);

  for (const size_t expected_rows : {size_t{4}, size_t{0}}) {
    CAPTURE(expected_rows);
    std::string message;
    try {
      readMasks(file, kMaskDataset, expected_rows);
      FAIL("readMasks accepted an absent dataset");
    } catch (const std::runtime_error& error) {
      message = error.what();
    }
    CHECK(message.find(kMaskDataset) != std::string::npos);
    CHECK(message.find("is missing") != std::string::npos);
  }
}
