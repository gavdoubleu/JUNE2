#include <algorithm>

#include "utils/event_logging/event_lookup_schemas.h"
#include "utils/event_logging/event_writer_detail.h"

namespace {

using june::event_writer_detail::openOrCreateGroup;

template <typename T>
H5::DataSet writeFixedDataset(H5::H5File& file, const std::string& name,
                              const std::vector<T>& data,
                              const H5::CompType& type, int compression_level) {
  hsize_t dims[1] = {data.size()};
  H5::DataSpace space(1, dims);
  H5::DataSet dataset;
  if (compression_level <= 0) {
    dataset = file.createDataSet(name, type, space);
  } else {
    H5::DSetCreatPropList plist = june::event_writer_detail::chunkedProperties(
        dims[0], compression_level);
    dataset = file.createDataSet(name, type, space, plist);
  }
  dataset.write(data.data(), type);
  return dataset;
}

std::string jsonString(const std::string& value) {
  std::string escaped = "\"";
  for (char c : value) {
    if (c == '\\' || c == '\"') escaped.push_back('\\');
    escaped.push_back(c);
  }
  escaped.push_back('\"');
  return escaped;
}

// Resolve which people should appear in a /lookups/* table for this write
// call. `mode` is the config switch ("none"/"all"/"infected_only"). Empty
// result means "skip the whole table". A table in "all" mode is written on
// the first call only (the append branch returns empty intentionally, it is
// already written) and ignores the filter. The filter is the people an
// "infected_only" table still needs; an empty filter means nobody.
std::unordered_set<june::PersonId> selectPeopleToSave(
    const june::WorldState& world, const std::string& mode,
    const std::unordered_set<june::PersonId>& person_ids_filter, bool append) {
  std::unordered_set<june::PersonId> people_to_save;
  if (mode == "all") {
    if (append) return people_to_save;
    for (const auto& person : world.people) people_to_save.insert(person.id);
  } else if (mode == "infected_only") {
    people_to_save = person_ids_filter;
  }
  return people_to_save;
}

// Project the selected people into the flat PersonRecord shape written to
// /lookups/people. Field order matches event_lookup_schema::person().
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
    june::copyFixedField(record.sex, sex_str);
    record.geo_unit_id = person.geo_unit_id;
    record.is_dead = person.is_dead ? 1 : 0;
    record.death_time = person.death_time;
    std::string sched_type =
        person.schedule_type_id < world.schedule_type_names.size()
            ? world.schedule_type_names[person.schedule_type_id]
            : june::kCouldNotResolve;
    june::copyFixedField(record.schedule_type, sched_type);
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
      else if (std::holds_alternative<std::vector<std::string>>(val)) {
        const auto& list = std::get<std::vector<std::string>>(val);
        std::string encoded = "[";
        for (size_t i = 0; i < list.size(); ++i) {
          if (i > 0) encoded += ",";
          encoded += jsonString(list[i]);
        }
        encoded += "]";
        values.push_back(std::move(encoded));
      } else
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
  H5::DSetCreatPropList plist =
      june::event_writer_detail::chunkedProperties(dims[0], compression_level);

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
        june::copyFixedField(record.activity_name, aname);
        record.venue_id = venues[idx].first;
        record.subset_index = venues[idx].second;
        record.activity_index = static_cast<int>(idx);
        records.push_back(record);
      }
    }
  }
  return records;
}

std::vector<june::PopulationSummaryRecord> buildPopulationSummaryRecords(
    const june::WorldState& world,
    const std::vector<std::string>& property_names) {
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
      if (k < property_names.size()) {
        auto p_opt = world.getPersonProperty(person, property_names[k]);
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

}  // namespace

namespace june::event_writer {

void writePersonLookupTable(
    H5::H5File& file, const WorldState& world, const Config& config,
    const std::unordered_set<PersonId>& person_ids_filter, bool append) {
  auto people_to_save =
      selectPeopleToSave(world, config.simulation.save_full_person_details,
                         person_ids_filter, append);
  if (people_to_save.empty()) return;

  auto records = buildPersonRecords(world, people_to_save);
  auto person_type = event_lookup_schema::person();
  writeDatasetTemplate(file, "/lookups/people", records, person_type,
                       config.simulation.compression_level);

  if (!world.person_property_names.empty()) {
    H5::Group prop_group = event_writer_detail::openOrCreateGroup(
        file, "/lookups/people_properties");
    for (const auto& key : world.person_property_names) {
      auto values = collectPropertyValues(world, people_to_save, key);
      event_writer_detail::writeStringDataset(prop_group, key, values, true);
    }
  }
}

void writeVenueLookupTable(H5::H5File& file, const WorldState& world,
                           const Config& config) {
  size_t n = world.venues.size();
  std::vector<detail::VenueRecord> records(n + 1);
  records[0].venue_id = INFECTION_SEED_VENUE_ID;
  june::copyFixedField(records[0].name, "infection_seed");
  june::copyFixedField(records[0].type, "infection_seed");
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
    june::copyFixedField(records[i + 1].name, vname);
    std::string vtype = (venue.id < 0 && venue.id != INFECTION_SEED_VENUE_ID)
                            ? "coordinated_encounter"
                            : (venue.type_id < world.venue_type_names.size()
                                   ? world.venue_type_names[venue.type_id]
                                   : june::kCouldNotResolve);
    june::copyFixedField(records[i + 1].type, vtype);
    records[i + 1].geo_unit_id = venue.geo_unit_id;
    records[i + 1].n_subsets = static_cast<int>(venue.subset_count);
  }

  auto vtype = event_lookup_schema::venue();

  writeFixedDataset(file, "/lookups/venues", records, vtype,
                    config.simulation.compression_level);
}

void writePersonActivitiesTable(
    H5::H5File& file, const WorldState& world, const Config& config,
    const std::unordered_set<PersonId>& person_ids_filter, bool append) {
  auto people_to_save =
      selectPeopleToSave(world, config.simulation.save_person_activities,
                         person_ids_filter, append);
  if (people_to_save.empty()) return;

  auto records = buildPersonActivityRecords(world, people_to_save);
  if (records.empty()) return;
  auto atype = event_lookup_schema::personActivity();
  writeDatasetTemplate(file, "/lookups/person_activities", records, atype,
                       config.simulation.compression_level);
}

void writePopulationSummary(H5::H5File& file, const WorldState& world,
                            const Config& config) {
  const size_t n = world.people.size();
  if (n == 0) return;

  // At most 4 properties fit (PopulationSummaryRecord::extra_codes width).
  const std::vector<std::string> property_names =
      config.simulation.summaryPropertiesFor(world);
  auto records = buildPopulationSummaryRecords(world, property_names);
  auto ptype = event_lookup_schema::populationSummary();

  H5::DataSet pds =
      writeFixedDataset(file, "/lookups/population_summary", records, ptype,
                        config.simulation.compression_level);

  for (size_t k = 0; k < property_names.size(); ++k) {
    H5::StrType stype(H5::PredType::C_S1, H5T_VARIABLE);
    H5::Attribute attr = pds.createAttribute("extra_prop_" + std::to_string(k),
                                             stype, H5::DataSpace(H5S_SCALAR));
    attr.write(stype, property_names[k]);
  }
}

void writePopulationNetworks(H5::H5File& file, const WorldState& world,
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

}  // namespace june::event_writer
