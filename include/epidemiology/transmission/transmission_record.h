#pragma once

#include <cstdint>

namespace june {

// What infected a Person. Stored as uint8 in the event log; append new values
// only, so existing on-disk values keep their meaning.
enum class InfectionSource : uint8_t { Person, Fomite, Compartmental, Seed };

// The facts of one transmission, as ids: what infected the Person, the
// infector's symptom id and the Transmission Mode index. 255 (kNoSymptomId,
// kNoModeIndex) marks a fact as absent. The one value both an Infection and
// its event-log row are built from, so the two cannot disagree.
struct TransmissionRecord {
  InfectionSource source;
  uint8_t infector_symptom_id;
  uint8_t transmission_mode_index;
};

}  // namespace june
