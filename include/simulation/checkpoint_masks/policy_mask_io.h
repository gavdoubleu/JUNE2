#pragma once

#include <H5Cpp.h>

#include <cstddef>
#include <string>
#include <vector>

#include "core/types.h"

namespace june {

// Checkpoint I/O for per-Person PolicyMasks: one row per Person, stored as an
// N x kPolicyMaskWords U64 dataset.

void writeMasks(H5::H5File& file, const std::string& name,
                const std::vector<PolicyMask>& masks, int compression);

// Throws unless the dataset is N x kPolicyMaskWords with N == expected_rows.
// A missing dataset counts as zero rows.
std::vector<PolicyMask> readMasks(H5::H5File& file, const std::string& name,
                                  size_t expected_rows);

}  // namespace june
