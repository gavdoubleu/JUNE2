#include <algorithm>
#include <cstring>

#include "utils/event_logging/event_logger.h"
#include "utils/event_logging/event_writer.h"
#include "utils/event_logging/event_writer_detail.h"

namespace {

using june::event_writer_detail::openOrCreateGroup;

// Resolve which people should appear in a /lookups/* table for this write
// call. `mode` is the config switch ("none"/"all"/"infected_only"). Empty
// result means "skip the whole table". The mode=="all" + append branch
// returns empty intentionally (already written on first call).
std::unordered_set<june::PersonId> selectPeopleToSave(
    const june::WorldState& world, const std::string& mode,
    const june::EventLogger& logger, bool append,
    const std::unordered_set<june::PersonId>* person_ids_filter) {
  std::unordered_set<june::PersonId> people_to_save;
  if (person_ids_filter) {
    people_to_save = *person_ids_filter;
  } else if (mode == "all") {
    if (append) return people_to_save;
    for (const auto& person : world.people) people_to_save.insert(person.id);
  } else if (mode == "infected_only") {
    people_to_save = logger.getInfectedPersonIds();
  }
  return people_to_save;
}

// Project the selected people into the flat PersonRecord shape written to
// /lookups/people. Field order matches makePersonCompType().
std::vector<june::detail::PersonRecord> buildPersonRecords(
    const june::WorldState& world,
    const std::unordered_set<june::PersonId>& people_to_save) {
  std::vector<june::detail::PersonRecord> records;
  records.reserve(people_to_save.size());
  for (const auto& person : world.people) {
    if (!people_to_save.count(person.id)) continue;
    june::detail::PersonRecord record;
    record.person_id = person.id;
    record.age = person.age;
    std::string sex_str = (person.sex == june::Sex::MALE) ? "male"
                          : (person.sex == june::Sex::FEMALE)
                              ? "female"
                              : june::kCouldNotResolve;
    strncpy(record.sex, sex_str.c_str(), 15);
    record.sex[15] = '\0';
    record.geo_unit_id = person.geo_unit_id;
    record.is_dead = person.is_dead ? 1 : 0;
    record.death_time = person.death_time;
    std::string sched_type =
        person.schedule_type_id < world.schedule_type_names.size()
            ? world.schedule_type_names[person.schedule_type_id]
            : june::kCouldNotResolve;
    strncpy(record.schedule_type, sched_type.c_str(), 63);
    record.schedule_type[63] = '\0';
    record.num_activities =
        static_cast<int>(world.getActivityMetas(person).size());
    record.num_residence_venues =
        static_cast<int>(world.getActivityVenues(person, "residence").size());
    record.num_primary_activities = static_cast<int>(
        world.getActivityVenues(person, "primary_activity").size());
    record.num_leisure_venues =
        static_cast<int>(world.getActivityVenues(person, "leisure").size());
    record.num_medical_facilities = static_cast<int>(
        world.getActivityVenues(person, "medical_facility").size());
    records.push_back(record);
  }
  return records;
}

// Stringify one person-property column across the selected people, in
// world.people order. Falls back to network-partner serialisation for keys
// whose values live in Person::NetworkMeta rather than the flat property
// table; falls back to "not applicable" when neither source has a value.
std::vector<std::string> collectPropertyValues(
    const june::WorldState& world,
    const std::unordered_set<june::PersonId>& people_to_save,
    const std::string& key) {
  const int network_type_id = world.getNetworkTypeIndex(key);
  std::vector<std::string> values;
  values.reserve(people_to_save.size());
  for (const auto& person : world.people) {
    if (!people_to_save.count(person.id)) continue;
    auto prop = world.getPersonProperty(person, key);
    if (prop.has_value()) {
      const auto& val = *prop;
      if (std::holds_alternative<std::string>(val))
        values.push_back(std::get<std::string>(val));
      else if (std::holds_alternative<int32_t>(val)) {
        int32_t iv = std::get<int32_t>(val);
        auto rit = world.person_property_value_registries.find(key);
        if (rit != world.person_property_value_registries.end() && iv >= 0 &&
            (size_t)iv < rit->second.size())
          values.push_back(rit->second[iv]);
        else
          values.push_back(std::to_string(iv));
      } else if (std::holds_alternative<bool>(val))
        values.push_back(std::get<bool>(val) ? "true" : "false");
      else if (std::holds_alternative<double>(val))
        values.push_back(std::to_string(std::get<double>(val)));
      else
        values.push_back(june::kCouldNotResolve);
    } else if (network_type_id >= 0) {
      auto partners = world.getNetworkPartners(person, network_type_id);
      if (partners.empty()) {
        values.push_back(june::kNotApplicable);
      } else {
        std::string s = "[";
        for (size_t i = 0; i < partners.size(); ++i) {
          if (i > 0) s.push_back(' ');
          s += std::to_string(partners[i]);
        }
        s.push_back(']');
        values.push_back(std::move(s));
      }
    } else {
      values.push_back(june::kNotApplicable);
    }
  }
  return values;
}

// Write `values` to `name` under `prop_group` as a chunked var-length string
// dataset. Creates the dataset if absent, extends + appends if present.
void writeOrAppendStringDataset(H5::Group& prop_group, const std::string& name,
                                const std::vector<std::string>& values) {
  hsize_t out_dims[1] = {values.size()};
  hsize_t out_maxdims[1] = {H5S_UNLIMITED};
  H5::DataSpace out_space(1, out_dims, out_maxdims);
  H5::StrType out_type(H5::PredType::C_S1, H5T_VARIABLE);

  H5::DSetCreatPropList plist;
  hsize_t chunk_dims[1] = {std::min(out_dims[0], hsize_t(100000))};
  if (chunk_dims[0] == 0) chunk_dims[0] = 1;
  plist.setChunk(1, chunk_dims);

  std::vector<const char*> pc_strs;
  pc_strs.reserve(values.size());
  for (const auto& s : values) pc_strs.push_back(s.c_str());

  H5::DataSet pds;
  if (H5Lexists(prop_group.getId(), name.c_str(), H5P_DEFAULT)) {
    pds = prop_group.openDataSet(name);
    hsize_t current_dims[1];
    pds.getSpace().getSimpleExtentDims(current_dims);
    hsize_t new_dims[1] = {current_dims[0] + values.size()};
    pds.extend(new_dims);

    H5::DataSpace file_space = pds.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, out_dims, current_dims);
    pds.write(pc_strs.data(), out_type, out_space, file_space);
  } else {
    pds = prop_group.createDataSet(name, out_type, out_space, plist);
    pds.write(pc_strs.data(), out_type);
  }
}

// Dump (person_id, partner_id) pairs for one network type as two parallel
// int32 datasets under `networks_group/network_name`. No-op when no person
// has a partner in that network. Existing datasets are unlinked first so
// each call writes a fresh snapshot.
void writeOneNetworkLookup(H5::Group& networks_group,
                           const june::WorldState& world,
                           const std::string& network_name, int type_id,
                           int compression_level) {
  std::vector<int32_t> persons;
  std::vector<int32_t> partners;
  persons.reserve(world.people.size());
  partners.reserve(world.people.size());

  bool has_any = false;
  for (const auto& person : world.people) {
    auto partner_ids = world.getNetworkPartners(person, type_id);
    if (partner_ids.empty()) continue;
    has_any = true;
    for (const auto& pid : partner_ids) {
      persons.push_back(person.id);
      partners.push_back(static_cast<int32_t>(pid));
    }
  }
  if (!has_any) return;

  H5::Group net_group = openOrCreateGroup(networks_group, network_name);

  hsize_t dims[1] = {persons.size()};
  H5::DataSpace space(1, dims);
  H5::DSetCreatPropList plist;
  hsize_t chunk_dims[1] = {std::min(dims[0], hsize_t(100000))};
  if (chunk_dims[0] == 0) chunk_dims[0] = 1;
  plist.setChunk(1, chunk_dims);
  plist.setDeflate(compression_level);

  if (H5Lexists(net_group.getId(), "person_id", H5P_DEFAULT))
    net_group.unlink("person_id");
  if (H5Lexists(net_group.getId(), "partner_id", H5P_DEFAULT))
    net_group.unlink("partner_id");

  H5::DataSet p_ds = net_group.createDataSet(
      "person_id", H5::PredType::NATIVE_INT32, space, plist);
  p_ds.write(persons.data(), H5::PredType::NATIVE_INT32);

  H5::DataSet q_ds = net_group.createDataSet(
      "partner_id", H5::PredType::NATIVE_INT32, space, plist);
  q_ds.write(partners.data(), H5::PredType::NATIVE_INT32);
}

std::vector<june::detail::PersonActivityRecord> buildPersonActivityRecords(
    const june::WorldState& world,
    const std::unordered_set<june::PersonId>& people_to_save) {
  size_t total_entries = 0;
  for (const auto& person : world.people) {
    if (!people_to_save.count(person.id)) continue;
    for (const auto& meta : world.getActivityMetas(person))
      total_entries += world.getActivityVenues(meta).size();
  }
  std::vector<june::detail::PersonActivityRecord> records;
  records.reserve(total_entries);
  for (const auto& person : world.people) {
    if (!people_to_save.count(person.id)) continue;
    for (const auto& meta : world.getActivityMetas(person)) {
      if (meta.activity_index < 0 ||
          meta.activity_index >= (int16_t)world.activity_names.size())
        continue;
      const std::string& aname = world.activity_names[meta.activity_index];
      auto venues = world.getActivityVenues(meta);
      for (size_t idx = 0; idx < venues.size(); ++idx) {
        june::detail::PersonActivityRecord record;
        record.person_id = person.id;
        strncpy(record.activity_name, aname.c_str(), 63);
        record.activity_name[63] = '\0';
        record.venue_id = venues[idx].first;
        record.subset_index = venues[idx].second;
        record.activity_index = static_cast<int>(idx);
        records.push_back(record);
      }
    }
  }
  return records;
}

H5::CompType makePersonActivityCompType() {
  H5::StrType aname_type(H5::PredType::C_S1, 64);
  H5::CompType atype(sizeof(june::detail::PersonActivityRecord));
  atype.insertMember("person_id",
                     HOFFSET(june::detail::PersonActivityRecord, person_id),
                     H5::PredType::NATIVE_INT);
  atype.insertMember("activity_name",
                     HOFFSET(june::detail::PersonActivityRecord, activity_name),
                     aname_type);
  atype.insertMember("venue_id",
                     HOFFSET(june::detail::PersonActivityRecord, venue_id),
                     H5::PredType::NATIVE_INT);
  atype.insertMember("subset_index",
                     HOFFSET(june::detail::PersonActivityRecord, subset_index),
                     H5::PredType::NATIVE_INT);
  atype.insertMember(
      "activity_index",
      HOFFSET(june::detail::PersonActivityRecord, activity_index),
      H5::PredType::NATIVE_INT);
  return atype;
}

// Metadata for one /lookups/population_summary "extra" property slot.
// At most 4 are kept (matches PopulationSummaryRecord::extra_codes width).
struct SummaryPropInfo {
  std::string name;
  int index;
  const std::vector<std::string>* registry = nullptr;
};

std::vector<SummaryPropInfo> collectSummaryPropertyMetadata(
    const june::WorldState& world,
    const std::vector<std::string>& summary_props) {
  std::vector<SummaryPropInfo> props_metadata;
  for (const auto& name : summary_props) {
    int idx = world.getPersonPropertyIndex(name);
    if (idx < 0) continue;
    SummaryPropInfo pi;
    pi.name = name;
    pi.index = idx;
    auto it = world.person_property_value_registries.find(name);
    if (it != world.person_property_value_registries.end())
      pi.registry = &it->second;
    props_metadata.push_back(pi);
  }
  return props_metadata;
}

std::vector<june::PopulationSummaryRecord> buildPopulationSummaryRecords(
    const june::WorldState& world,
    const std::vector<SummaryPropInfo>& props_metadata) {
  const size_t n = world.people.size();
  std::vector<june::PopulationSummaryRecord> records(n);
  for (size_t i = 0; i < n; ++i) {
    const june::Person& person = world.people[i];
    records[i].person_id = person.id;
    records[i].age_group =
        std::min(static_cast<uint8_t>(person.age / 5), uint8_t(17));
    records[i].sex_code = static_cast<uint8_t>(person.sex);
    records[i].schedule_type_code =
        static_cast<uint8_t>(person.schedule_type_id % 256);
    records[i].reserved = 0;
    records[i].geo_unit_id = person.geo_unit_id;
    for (size_t k = 0; k < 4; ++k) {
      records[i].extra_codes[k] = 0;
      if (k < props_metadata.size()) {
        auto p_opt = world.getPersonProperty(person, props_metadata[k].name);
        if (p_opt) {
          if (std::holds_alternative<int32_t>(*p_opt))
            records[i].extra_codes[k] =
                (uint8_t)(std::get<int32_t>(*p_opt) % 256);
          else if (std::holds_alternative<bool>(*p_opt))
            records[i].extra_codes[k] = std::get<bool>(*p_opt) ? 1 : 0;
        }
      }
    }
  }
  return records;
}

H5::CompType makePopulationSummaryCompType() {
  H5::CompType ptype(sizeof(june::PopulationSummaryRecord));
  ptype.insertMember("person_id",
                     HOFFSET(june::PopulationSummaryRecord, person_id),
                     H5::PredType::NATIVE_INT);
  ptype.insertMember("age_group",
                     HOFFSET(june::PopulationSummaryRecord, age_group),
                     H5::PredType::NATIVE_UINT8);
  ptype.insertMember("sex_code",
                     HOFFSET(june::PopulationSummaryRecord, sex_code),
                     H5::PredType::NATIVE_UINT8);
  ptype.insertMember("schedule_type_code",
                     HOFFSET(june::PopulationSummaryRecord, schedule_type_code),
                     H5::PredType::NATIVE_UINT8);
  ptype.insertMember("reserved",
                     HOFFSET(june::PopulationSummaryRecord, reserved),
                     H5::PredType::NATIVE_UINT8);
  ptype.insertMember("geo_unit_id",
                     HOFFSET(june::PopulationSummaryRecord, geo_unit_id),
                     H5::PredType::NATIVE_INT);
  hsize_t adims[1] = {4};
  H5::ArrayType atype(H5::PredType::NATIVE_UINT8, 1, adims);
  ptype.insertMember("extra_codes",
                     HOFFSET(june::PopulationSummaryRecord, extra_codes),
                     atype);
  return ptype;
}

H5::CompType makePersonCompType() {
  H5::StrType sex_type(H5::PredType::C_S1, 16);
  H5::StrType schedule_type(H5::PredType::C_S1, 64);
  H5::CompType person_type(sizeof(june::detail::PersonRecord));
  person_type.insertMember("person_id",
                           HOFFSET(june::detail::PersonRecord, person_id),
                           H5::PredType::NATIVE_INT);
  person_type.insertMember("age", HOFFSET(june::detail::PersonRecord, age),
                           H5::PredType::NATIVE_DOUBLE);
  person_type.insertMember("sex", HOFFSET(june::detail::PersonRecord, sex),
                           sex_type);
  person_type.insertMember("geo_unit_id",
                           HOFFSET(june::detail::PersonRecord, geo_unit_id),
                           H5::PredType::NATIVE_INT);
  person_type.insertMember("is_dead",
                           HOFFSET(june::detail::PersonRecord, is_dead),
                           H5::PredType::NATIVE_INT);
  person_type.insertMember("death_time",
                           HOFFSET(june::detail::PersonRecord, death_time),
                           H5::PredType::NATIVE_DOUBLE);
  person_type.insertMember("schedule_type",
                           HOFFSET(june::detail::PersonRecord, schedule_type),
                           schedule_type);
  person_type.insertMember("num_activities",
                           HOFFSET(june::detail::PersonRecord, num_activities),
                           H5::PredType::NATIVE_INT);
  person_type.insertMember(
      "num_residence_venues",
      HOFFSET(june::detail::PersonRecord, num_residence_venues),
      H5::PredType::NATIVE_INT);
  person_type.insertMember(
      "num_primary_activities",
      HOFFSET(june::detail::PersonRecord, num_primary_activities),
      H5::PredType::NATIVE_INT);
  person_type.insertMember(
      "num_leisure_venues",
      HOFFSET(june::detail::PersonRecord, num_leisure_venues),
      H5::PredType::NATIVE_INT);
  person_type.insertMember(
      "num_medical_facilities",
      HOFFSET(june::detail::PersonRecord, num_medical_facilities),
      H5::PredType::NATIVE_INT);
  return person_type;
}

}  // namespace

namespace june {

void EventWriter::writePersonLookupTable(
    H5::H5File& file, const WorldState& world, const Config& config,
    const EventLogger& logger, bool append,
    const std::unordered_set<PersonId>* person_ids_filter) {
  auto people_to_save =
      selectPeopleToSave(world, config.simulation.save_full_person_details,
                         logger, append, person_ids_filter);
  if (people_to_save.empty()) return;

  auto records = buildPersonRecords(world, people_to_save);
  auto person_type = makePersonCompType();
  writeDatasetTemplate(file, "/lookups/people", records, person_type,
                       config.simulation.compression_level);

  if (!world.person_property_names.empty()) {
    H5::Group prop_group = event_writer_detail::openOrCreateGroup(
        file, "/lookups/people_properties");
    for (const auto& key : world.person_property_names) {
      auto values = collectPropertyValues(world, people_to_save, key);
      writeOrAppendStringDataset(prop_group, key, values);
    }
  }
}

void EventWriter::writeVenueLookupTable(H5::H5File& file,
                                        const WorldState& world,
                                        const Config& config) {
  size_t n = world.venues.size();
  std::vector<detail::VenueRecord> records(n + 1);
  records[0].venue_id = INFECTION_SEED_VENUE_ID;
  strncpy(records[0].name, "infection_seed", 127);
  records[0].name[127] = '\0';
  strncpy(records[0].type, "infection_seed", 63);
  records[0].type[63] = '\0';
  records[0].geo_unit_id = -1;
  records[0].n_subsets = 0;

  for (size_t i = 0; i < n; ++i) {
    const Venue& venue = world.venues[i];
    records[i + 1].venue_id = venue.id;
    std::string vname =
        (venue.id < 0 && venue.id != INFECTION_SEED_VENUE_ID)
            ? "Virtual Coordinated Site (Rank " +
                  std::to_string((-(static_cast<long long>(venue.id) + 1000)) %
                                 1000000) +
                  ")"
            : "Venue_" + std::to_string(venue.id);
    strncpy(records[i + 1].name, vname.c_str(), 127);
    records[i + 1].name[127] = '\0';
    std::string vtype = (venue.id < 0 && venue.id != INFECTION_SEED_VENUE_ID)
                            ? "coordinated_encounter"
                            : (venue.type_id < world.venue_type_names.size()
                                   ? world.venue_type_names[venue.type_id]
                                   : june::kCouldNotResolve);
    strncpy(records[i + 1].type, vtype.c_str(), 63);
    records[i + 1].type[63] = '\0';
    records[i + 1].geo_unit_id = venue.geo_unit_id;
    records[i + 1].n_subsets = static_cast<int>(venue.subset_count);
  }

  H5::StrType ntype(H5::PredType::C_S1, 128);
  H5::StrType ttype(H5::PredType::C_S1, 64);
  H5::CompType vtype(sizeof(detail::VenueRecord));
  vtype.insertMember("venue_id", HOFFSET(detail::VenueRecord, venue_id),
                     H5::PredType::NATIVE_INT);
  vtype.insertMember("name", HOFFSET(detail::VenueRecord, name), ntype);
  vtype.insertMember("type", HOFFSET(detail::VenueRecord, type), ttype);
  vtype.insertMember("geo_unit_id", HOFFSET(detail::VenueRecord, geo_unit_id),
                     H5::PredType::NATIVE_INT);
  vtype.insertMember("n_subsets", HOFFSET(detail::VenueRecord, n_subsets),
                     H5::PredType::NATIVE_INT);

  hsize_t vdims[1] = {n + 1};
  H5::DataSpace vspace(1, vdims);
  H5::DataSet vds =
      createCompressedDataset(file, "/lookups/venues", vtype, vspace,
                              config.simulation.compression_level);
  vds.write(records.data(), vtype);
}

void EventWriter::writePersonActivitiesTable(
    H5::H5File& file, const WorldState& world, const Config& config,
    const EventLogger& logger, bool append,
    const std::unordered_set<PersonId>* person_ids_filter) {
  auto people_to_save =
      selectPeopleToSave(world, config.simulation.save_person_activities,
                         logger, append, person_ids_filter);
  if (people_to_save.empty()) return;

  auto records = buildPersonActivityRecords(world, people_to_save);
  if (records.empty()) return;
  auto atype = makePersonActivityCompType();
  writeDatasetTemplate(file, "/lookups/person_activities", records, atype,
                       config.simulation.compression_level);
}

void EventWriter::writePopulationSummary(H5::H5File& file,
                                         const WorldState& world,
                                         const Config& config) {
  const size_t n = world.people.size();
  if (n == 0) return;

  std::vector<std::string> summary_props = config.simulation.summary_properties;
  if (summary_props.size() > 4) summary_props.resize(4);

  auto props_metadata = collectSummaryPropertyMetadata(world, summary_props);
  auto records = buildPopulationSummaryRecords(world, props_metadata);
  auto ptype = makePopulationSummaryCompType();

  hsize_t pdims[1] = {n};
  H5::DataSpace pspace(1, pdims);
  H5::DataSet pds =
      createCompressedDataset(file, "/lookups/population_summary", ptype,
                              pspace, config.simulation.compression_level);
  pds.write(records.data(), ptype);

  for (size_t k = 0; k < props_metadata.size(); ++k) {
    H5::StrType stype(H5::PredType::C_S1, H5T_VARIABLE);
    H5::Attribute attr = pds.createAttribute("extra_prop_" + std::to_string(k),
                                             stype, H5::DataSpace(H5S_SCALAR));
    attr.write(stype, props_metadata[k].name);
  }
}

void EventWriter::writePopulationNetworks(H5::H5File& file,
                                          const WorldState& world,
                                          const Config& config) {
  if (world.people.empty() || world.network_names.empty()) return;

  H5::Group lookups_group =
      event_writer_detail::openOrCreateGroup(file, "/lookups");
  H5::Group networks_group = event_writer_detail::openOrCreateGroup(
      lookups_group, "population_networks");

  for (const auto& network_name : world.network_names) {
    const int type_id = world.getNetworkTypeIndex(network_name);
    if (type_id < 0) continue;
    writeOneNetworkLookup(networks_group, world, network_name, type_id,
                          config.simulation.compression_level);
  }
}

}  // namespace june
