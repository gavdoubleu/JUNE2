#ifdef USE_MPI

#include "parallel/visitor_payload.h"

namespace june {

Domain::VisitorData buildVisitorPayload(const PersonLocation& location,
                                        const Person& person, int home_rank,
                                        double slot_start,
                                        const Disease& disease,
                                        const EmissionCalculator& calculator) {
  Domain::VisitorData visitor;
  visitor.person_id = location.person_id;
  visitor.home_rank = home_rank;
  visitor.venue_id = location.venue_id;
  visitor.subset_idx = location.subset_index;
  visitor.is_infected = (person.infection != nullptr);
  visitor.is_infectious =
      visitor.is_infected && person.infection->isInfectious(slot_start);

  const double susceptibility =
      person.getSusceptibility(slot_start, disease.getName());
  visitor.immunity_level = static_cast<float>(1.0 - susceptibility);

  visitor.encounter_type_id = location.encounter_type_id;
  visitor.newly_infected = false;
  visitor.new_infection_time = -1.0;

  visitor.symptom_id =
      visitor.is_infected ? person.infection->symptomIdAt(slot_start) : 0;

  calculator.emit(person, slot_start, visitor.emission);
  return visitor;
}

}  // namespace june

#endif  // USE_MPI
