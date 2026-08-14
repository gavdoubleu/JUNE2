# Mirrors the sentinel constants in include/core/types.h and
# include/utils/event_logging/event_types.h. Keep in sync when the engine
# adds/renames a sentinel (see the "Consistent sentinels" rollout, #16).

SEED_VENUE_ID = -999  # INFECTION_SEED_VENUE_ID
NO_INFECTOR_ID = -1  # kInvalidPersonId, as used for infector_id

# Distinct engine constants that happen to share the uint8 value 255 — kept as
# separate names (not one shared UNSET_REGISTRY_INDEX) because each marks a
# different registry/column and the engine itself no longer collapses them.
DEFAULT_ENCOUNTER_TYPE_ID = 255  # kDefaultEncounterTypeId — ordinary (non-coordinated) encounter, not an error
NO_SYMPTOM_ID = 255  # kNoSymptomId — "not applicable", e.g. infector_symptom_id with no infector
