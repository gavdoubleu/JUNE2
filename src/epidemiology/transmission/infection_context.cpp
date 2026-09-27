#include "epidemiology/transmission/infection_context.h"

#include "core/types.h"
#include "epidemiology/disease.h"

namespace june {

const std::vector<std::string>& infectionSourceNames() {
  static const std::vector<std::string> names = {"person", "fomite",
                                                 "compartmental", "seed"};
  return names;
}

InfectionContext buildInfectionContext(const TransmissionRecord& record,
                                       const Disease& disease) {
  InfectionContext context;
  context.infection_source =
      infectionSourceNames()[static_cast<uint8_t>(record.source)];
  if (record.infector_symptom_id != kNoSymptomId) {
    context.infector_symptom =
        disease.getSymptomName(record.infector_symptom_id);
  }
  if (record.transmission_mode_index != kNoModeIndex) {
    context.transmission_mode =
        disease.getModeName(record.transmission_mode_index);
  }
  return context;
}

}  // namespace june
