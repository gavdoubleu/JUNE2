#include "simulation/checkpoint_masks/policy_mask_io.h"

#include <cstdint>
#include <stdexcept>

#include "simulation/checkpoint_datasets/chunked_dataset.h"

namespace june {

void writeMasks(H5::H5File& file, const std::string& name,
                const std::vector<PolicyMask>& masks, int compression) {
  std::vector<uint64_t> words;
  words.reserve(masks.size() * kPolicyMaskWords);
  for (const PolicyMask& mask : masks)
    for (size_t word = 0; word < kPolicyMaskWords; ++word)
      words.push_back(policyMaskWord(mask, word));
  const hsize_t dims[2] = {masks.size(), kPolicyMaskWords};
  writeChunkedDataSet(file, name, words.data(), H5::PredType::NATIVE_UINT64,
                      dims, compression);
}

std::vector<PolicyMask> readMasks(H5::H5File& file, const std::string& name,
                                  size_t expected_rows) {
  if (H5Lexists(file.getId(), name.c_str(), H5P_DEFAULT) <= 0) {
    throw std::runtime_error("checkpoint dataset '" + name +
                             "' is missing (corrupt or partially written "
                             "checkpoint)");
  }
  H5::DataSet dataset = file.openDataSet(name);
  H5::DataSpace space = dataset.getSpace();
  hsize_t dim[2] = {0, 0};
  if (space.getSimpleExtentNdims() == 2) space.getSimpleExtentDims(dim);
  if (dim[1] != kPolicyMaskWords) {
    throw std::runtime_error(
        "checkpoint dataset '" + name + "' is not an N x " +
        std::to_string(kPolicyMaskWords) + " policy mask (corrupt checkpoint)");
  }
  if (dim[0] != expected_rows) {
    throw std::runtime_error("checkpoint dataset '" + name + "' holds " +
                             std::to_string(dim[0]) + " rows; population has " +
                             std::to_string(expected_rows));
  }
  std::vector<uint64_t> words(dim[0] * kPolicyMaskWords);
  if (dim[0]) dataset.read(words.data(), H5::PredType::NATIVE_UINT64);
  std::vector<PolicyMask> masks(dim[0]);
  for (hsize_t i = 0; i < dim[0]; ++i)
    for (size_t w = 0; w < kPolicyMaskWords; ++w)
      setPolicyMaskWord(masks[i], w, words[i * kPolicyMaskWords + w]);
  return masks;
}

}  // namespace june
