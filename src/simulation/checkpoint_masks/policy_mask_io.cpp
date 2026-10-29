#include "simulation/checkpoint_masks/policy_mask_io.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>

namespace june {

void writeMasks(H5::H5File& file, const std::string& name,
                const std::vector<PolicyMask>& masks, int compression) {
  std::vector<uint64_t> words;
  words.reserve(masks.size() * kPolicyMaskWords);
  for (const PolicyMask& mask : masks)
    for (size_t w = 0; w < kPolicyMaskWords; ++w)
      words.push_back(policyMaskWord(mask, w));
  const auto& uint64_type = H5::PredType::NATIVE_UINT64;
  hsize_t dims[2] = {masks.size(), kPolicyMaskWords};
  H5::DataSpace space(2, dims);
  if (masks.empty()) {
    file.createDataSet(name, uint64_type, space);
    return;
  }
  H5::DSetCreatPropList plist;
  hsize_t chunk[2] = {std::min<hsize_t>(dims[0], 100000), kPolicyMaskWords};
  plist.setChunk(2, chunk);
  if (compression > 0) plist.setDeflate(compression);
  H5::DataSet dataset = file.createDataSet(name, uint64_type, space, plist);
  dataset.write(words.data(), uint64_type);
}

std::vector<PolicyMask> readMasks(H5::H5File& file, const std::string& name,
                                  size_t expected_rows) {
  hsize_t dim[2] = {0, kPolicyMaskWords};
  H5::DataSet dataset;
  const bool exists = H5Lexists(file.getId(), name.c_str(), H5P_DEFAULT) > 0;
  if (exists) {
    dataset = file.openDataSet(name);
    H5::DataSpace space = dataset.getSpace();
    dim[1] = 0;
    if (space.getSimpleExtentNdims() == 2) space.getSimpleExtentDims(dim);
    if (dim[1] != kPolicyMaskWords) {
      throw std::runtime_error(
          "checkpoint dataset '" + name + "' is not an N x " +
          std::to_string(kPolicyMaskWords) +
          " policy mask; checkpoint predates " +
          std::to_string(kMaxPoliciesPerKind) + "-policy masks");
    }
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
