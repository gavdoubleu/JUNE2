#pragma once

#include <string>
#include <vector>

#include "epidemiology/transmission/transmission_record.h"
#include "utils/filtering.h"

namespace june {

class Disease;

// Filter names of every InfectionSource, indexed by its uint8 value. The one
// list both context building and outcome-table validation use.
const std::vector<std::string>& infectionSourceNames();

// The Infection Context `disease` names for `record`. An absent fact (255)
// gives an empty name, so it fails every criterion on that fact.
InfectionContext buildInfectionContext(const TransmissionRecord& record,
                                       const Disease& disease);

}  // namespace june
