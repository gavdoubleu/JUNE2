#pragma once

#include <vector>

#include "core/config.h"

// Forward-declare YAML::Node to keep yaml-cpp out of the public include
// surface. Only TUs that actually parse YAML need the full header.
namespace YAML {
class Node;
}  // namespace YAML

namespace june {
namespace config_detail {

// True on MPI rank 0 (and unconditionally true when MPI is not initialised
// or not compiled in). Shared by the rank-gated [CoordEnc] / vaccination
// load-time log lines that fire from multiple TUs.
bool logRank0();

// Parse a YAML sequence of `{property, operator, value}` entries into a
// vector of SelectionCriterion. The scalar `value` is dispatched
// int -> double -> string; a sequence `value` becomes vector<int32_t>, or
// vector<string> when it does not parse as ints.
// Shared by the loadSchedule / loadActivityPreferences code in
// config_loader.cpp and by loadVaccination in config_loader_vaccination.cpp.
void parseSelectionCriteria(const YAML::Node& selection_node,
                            std::vector<SelectionCriterion>& out);

}  // namespace config_detail
}  // namespace june
