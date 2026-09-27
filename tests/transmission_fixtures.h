#pragma once

#include "core/types.h"
#include "epidemiology/transmission/transmission_record.h"

namespace june {

// A Person-sourced transmission with no infector symptom or mode, for tests
// whose outcome rows don't filter on the Infection Context.
inline constexpr TransmissionRecord kNoTransmissionContext{
    InfectionSource::Person, kNoSymptomId, kNoModeIndex};

}  // namespace june
