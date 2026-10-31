#include "simulation/checkpoint_datasets/chunked_dataset.h"

#include <algorithm>
#include <vector>

namespace june {

void writeChunkedDataSet(H5::H5File& file, const std::string& name,
                         const void* data, const H5::PredType& type,
                         std::span<const hsize_t> dims, int compression) {
  const int rank = static_cast<int>(dims.size());
  H5::DataSpace space(rank, dims.data());
  if (dims[0] == 0) {
    file.createDataSet(name, type, space);
    return;
  }
  std::vector<hsize_t> chunk(dims.begin(), dims.end());
  chunk[0] = std::min<hsize_t>(dims[0], 100000);
  H5::DSetCreatPropList plist;
  plist.setChunk(rank, chunk.data());
  if (compression > 0) plist.setDeflate(compression);
  H5::DataSet dataset = file.createDataSet(name, type, space, plist);
  dataset.write(data, type);
}

}  // namespace june
