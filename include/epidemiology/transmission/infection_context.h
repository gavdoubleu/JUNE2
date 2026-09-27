#pragma once

#include "epidemiology/transmission/transmission_record.h"
#include "utils/filtering.h"

namespace june {

class Disease;

// The Infection Context `disease` names for `record`. An absent fact (255)
// gives an empty name, so it fails every criterion on that fact.
InfectionContext buildInfectionContext(const TransmissionRecord& record,
                                       const Disease& disease);

}  // namespace june
