#pragma once

#include <H5Cpp.h>

#include <span>
#include <string>

namespace june {

// Create and write a numeric checkpoint dataset of the given dims, chunked by
// up to 100000 rows and deflated at `compression` (0 for none). An empty one
// gets no chunk or filter, which HDF5 requires for zero extent.
void writeChunkedDataSet(H5::H5File& file, const std::string& name,
                         const void* data, const H5::PredType& type,
                         std::span<const hsize_t> dims, int compression);

}  // namespace june
