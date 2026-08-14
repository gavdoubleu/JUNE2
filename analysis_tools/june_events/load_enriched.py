import logging

from .decode import (
    DEFAULT_ENCOUNTER_TYPE_ID,
    NO_SYMPTOM_ID,
    decode_registry_column,
    load_registry,
)
from .enrich import enrich_with_people, enrich_with_venues
from .io import load_people_lookup, load_raw_table, load_venues_lookup

logger = logging.getLogger(__name__)

DEFAULT_REGISTRY_COLUMNS = {
    "encounter_type_id": "encounter_types",
    "infector_symptom_id": "symptoms",
    "old_symptom_id": "symptoms",
    "new_symptom_id": "symptoms",
}

# encounter_type_id and infector_symptom_id each have a real engine-side
# sentinel (uint8_t 255, event_types.h/types.h — kDefaultEncounterTypeId,
# kNoSymptomId), written directly by the engine (no NaN-masking needed on
# this side). old_symptom_id/new_symptom_id are always a real transition, no
# sentinel, so they get no entry here — decode_registry_column then treats
# every index as ordinary registry data.
REGISTRY_SENTINELS = {
    "encounter_type_id": DEFAULT_ENCOUNTER_TYPE_ID,
    "infector_symptom_id": NO_SYMPTOM_ID,
}

# Each sentinel means something more specific than decode_registry_column's
# generic "unknown": encounter_type_id's is the majority "ordinary venue-level
# force of infection" case, not a missing value; infector_symptom_id's is "no
# infector to read a symptom from" (seed/fomite/compartmental infections).
UNSET_LABEL_OVERRIDES = {
    "encounter_type_id": "regular_non_coordinated_encounter",
    "infector_symptom_id": "no_infector",
}


def _decoded_column_name(id_column: str) -> str:
    if id_column.endswith("_id"):
        return id_column[: -len("_id")]
    return id_column


def load_enriched_events(
    path: str,
    dataset_path: str,
    *,
    with_people: bool = True,
    with_venues: bool = True,
    registry_columns: dict | None = None,
    people_prefix: str = "person_",
    venues_prefix: str = "venue_",
):
    events = load_raw_table(path, dataset_path)
    if events is None:
        return None

    if registry_columns is None:
        registry_columns = DEFAULT_REGISTRY_COLUMNS

    loaded_registries = {}
    for id_column, registry_name in registry_columns.items():
        if id_column not in events.columns:
            continue
        if registry_name not in loaded_registries:
            loaded_registries[registry_name] = load_registry(path, registry_name)
        registry = loaded_registries[registry_name]
        if registry is None:
            logger.warning(
                "registry %r not found in %r, skipping decode of %r",
                registry_name,
                path,
                id_column,
            )
            continue
        decode_kwargs = {}
        if id_column in REGISTRY_SENTINELS:
            decode_kwargs["unset_value"] = REGISTRY_SENTINELS[id_column]
        if id_column in UNSET_LABEL_OVERRIDES:
            decode_kwargs["unset_label"] = UNSET_LABEL_OVERRIDES[id_column]
        events[_decoded_column_name(id_column)] = decode_registry_column(
            events, id_column, registry, **decode_kwargs
        )

    if with_people and "person_id" in events.columns:
        people_lookup = load_people_lookup(path)
        if people_lookup is not None:
            events = enrich_with_people(events, people_lookup, prefix=people_prefix)

    if with_venues and "venue_id" in events.columns:
        venues_lookup = load_venues_lookup(path)
        if venues_lookup is not None:
            events = enrich_with_venues(events, venues_lookup, prefix=venues_prefix)

    return events
