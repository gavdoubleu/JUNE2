#include <H5Cpp.h>

#include <cstddef>
#include <cstring>
#include <new>
#include <string>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "test_utils.h"
#include "utils/event_logging/event_logger.h"
#include "utils/event_logging/event_merger.h"
#include "utils/event_logging/event_schemas.h"
#include "utils/event_logging/event_types.h"

using namespace june;

namespace {

// Constructs T over memory full of 0xAB and returns the bytes that no named
// data field covers. Every one of them must be zero, or the record leaks
// whatever the allocator left in it when it is written to HDF5.
template <typename T, typename Make>
size_t nonzeroGapBytes(Make make, const std::vector<std::pair<size_t, size_t>>&
                                      data_fields /* offset, size */) {
  alignas(T) unsigned char buffer[sizeof(T)];
  std::memset(buffer, 0xAB, sizeof buffer);
  make(buffer);
  std::vector<bool> is_data(sizeof(T), false);
  for (const auto& [offset, size] : data_fields)
    for (size_t i = 0; i < size; ++i) is_data[offset + i] = true;
  size_t bad = 0;
  for (size_t i = 0; i < sizeof(T); ++i)
    if (!is_data[i] && buffer[i] != 0) ++bad;
  return bad;
}

}  // namespace

TEST_CASE("event records written raw to HDF5 carry zero padding") {
  CHECK(nonzeroGapBytes<InfectionEvent>(
            [](void* p) {
              new (p) InfectionEvent{.person_id = 1,
                                     .infector_id = 2,
                                     .venue_id = 3,
                                     .time = 4.0,
                                     .encounter_type_id = 5,
                                     .transmission_mode_index = 6,
                                     .infector_symptom_id = 7,
                                     .source = InfectionSource::Person};
            },
            {{offsetof(InfectionEvent, person_id), 4},
             {offsetof(InfectionEvent, infector_id), 4},
             {offsetof(InfectionEvent, venue_id), 4},
             {offsetof(InfectionEvent, time), 8},
             {offsetof(InfectionEvent, encounter_type_id), 1},
             {offsetof(InfectionEvent, transmission_mode_index), 1},
             {offsetof(InfectionEvent, infector_symptom_id), 1},
             {offsetof(InfectionEvent, source), 1},
             {offsetof(InfectionEvent, reserved), 1}}) == 0);

  CHECK(nonzeroGapBytes<SymptomChangeEvent>(
            [](void* p) { new (p) SymptomChangeEvent{1, 2, 3.0, 4, 5}; },
            {{offsetof(SymptomChangeEvent, person_id), 4},
             {offsetof(SymptomChangeEvent, venue_id), 4},
             {offsetof(SymptomChangeEvent, time), 8},
             {offsetof(SymptomChangeEvent, old_symptom_id), 1},
             {offsetof(SymptomChangeEvent, new_symptom_id), 1}}) == 0);

  CHECK(nonzeroGapBytes<CoordinatedEncounterEvent>(
            [](void* p) {
              new (p) CoordinatedEncounterEvent{.person_a = 1,
                                                .person_b = 2,
                                                .time = 3.0,
                                                .encounter_type_id = 4,
                                                .slot = 5,
                                                .group_id = 6};
            },
            {{offsetof(CoordinatedEncounterEvent, person_a), 4},
             {offsetof(CoordinatedEncounterEvent, person_b), 4},
             {offsetof(CoordinatedEncounterEvent, time), 8},
             {offsetof(CoordinatedEncounterEvent, encounter_type_id), 1},
             {offsetof(CoordinatedEncounterEvent, slot), 4},
             {offsetof(CoordinatedEncounterEvent, group_id), 8}}) == 0);

  CHECK(
      nonzeroGapBytes<FollowEvent>(
          [](void* p) {
            new (p) FollowEvent{
                .host = 1, .follower = 2, .time = 3.0, .rule_id = 4, .slot = 5};
          },
          {{offsetof(FollowEvent, host), 4},
           {offsetof(FollowEvent, follower), 4},
           {offsetof(FollowEvent, time), 8},
           {offsetof(FollowEvent, rule_id), 1},
           {offsetof(FollowEvent, slot), 4}}) == 0);
}

// ---------------------------------------------------------------------------
// Campaign day summary rows
// ---------------------------------------------------------------------------

namespace {

CampaignDayEvent campaignDayRow(int day, uint16_t campaign, uint16_t tier,
                                int64_t allotment) {
  // Distinct values per field, so a swapped or dropped field shows.
  return {.day = day,
          .campaign_id = campaign,
          .tier_id = tier,
          .candidates = 1000 + day,
          .doses_given = 100 + day,
          .declined = 10 + day,
          .allotment = allotment,
          .started_total = 5000 + 7 * day};
}

void logRow(EventLogger& logger, const CampaignDayEvent& r) {
  logger.logCampaignDay(r.day, r.campaign_id, r.tier_id, r.candidates,
                        r.doses_given, r.declined, r.allotment,
                        r.started_total);
}

// The record has no padding, so equal bytes mean equal fields.
bool sameRows(const std::vector<CampaignDayEvent>& a,
              const std::vector<CampaignDayEvent>& b) {
  return a.size() == b.size() &&
         (a.empty() ||
          std::memcmp(a.data(), b.data(), a.size() * sizeof(a[0])) == 0);
}

Config configWithCampaigns(const std::vector<std::string>& names) {
  Config config;
  // Nothing but the events and registries is under test; skip the lookups that
  // would need a populated world.
  config.simulation.save_full_person_details = "none";
  config.simulation.save_population_summary = false;
  for (const auto& name : names) {
    VaccinationCampaignConfig campaign;
    campaign.name = name;
    config.vaccination.campaigns.push_back(campaign);
  }
  return config;
}

std::vector<CampaignDayEvent> readCampaignDays(const std::string& filename) {
  H5::H5File file(filename, H5F_ACC_RDONLY);
  H5::DataSet dataset = file.openDataSet("/events/campaign_days");
  hsize_t count = 0;
  dataset.getSpace().getSimpleExtentDims(&count);
  std::vector<CampaignDayEvent> rows(count);
  dataset.read(rows.data(), event_schema::campaignDay());
  return rows;
}

std::vector<std::string> readRegistry(const std::string& filename,
                                      const std::string& path) {
  H5::H5File file(filename, H5F_ACC_RDONLY);
  H5::DataSet dataset = file.openDataSet(path);
  hsize_t count = 0;
  dataset.getSpace().getSimpleExtentDims(&count);
  H5::StrType type(H5::PredType::C_S1, H5T_VARIABLE);
  std::vector<char*> raw(count);
  dataset.read(raw.data(), type);
  std::vector<std::string> names(raw.begin(), raw.end());
  H5Dvlen_reclaim(type.getId(), dataset.getSpace().getId(), H5P_DEFAULT,
                  raw.data());
  return names;
}

}  // namespace

TEST_CASE("campaign day rows carry no padding") {
  static_assert(sizeof(CampaignDayEvent) == 48);
  CHECK(nonzeroGapBytes<CampaignDayEvent>(
            [](void* p) { new (p) CampaignDayEvent{1, 2, 3, 4, 5, 6, 7, 8}; },
            {{offsetof(CampaignDayEvent, day), 4},
             {offsetof(CampaignDayEvent, campaign_id), 2},
             {offsetof(CampaignDayEvent, tier_id), 2},
             {offsetof(CampaignDayEvent, candidates), 8},
             {offsetof(CampaignDayEvent, doses_given), 8},
             {offsetof(CampaignDayEvent, declined), 8},
             {offsetof(CampaignDayEvent, allotment), 8},
             {offsetof(CampaignDayEvent, started_total), 8}}) == 0);
  // A real allotment is never negative, so "none" cannot be mistaken for one.
  CHECK(kNoAllotment < 0);
}

TEST_CASE("campaign days written through the real writer read back intact") {
  ScopedTestFiles dir("june_campaign_days");
  const std::string file = dir.path("events.h5").string();
  const Config config = configWithCampaigns({"Zeta", "Alpha", "Mid"});
  const WorldState world;

  EventLogger logger;
  const std::vector<CampaignDayEvent> first = {
      campaignDayRow(0, 0, 0, 40), campaignDayRow(0, 0, 1, 40),
      campaignDayRow(0, 2, 0, kNoAllotment)};
  for (const auto& r : first) logRow(logger, r);
  CHECK(logger.getTotalRecordCount() == first.size());
  logger.flush(file, config, world, {});
  CHECK(logger.getTotalRecordCount() == 0);  // a flush empties the buffer

  // The dataset has the field names and widths the schema promises.
  {
    H5::H5File f(file, H5F_ACC_RDONLY);
    H5::CompType stored = f.openDataSet("/events/campaign_days").getCompType();
    const std::vector<std::string> expected = {
        "day",         "campaign_id", "tier_id",   "candidates",
        "doses_given", "declined",    "allotment", "started_total"};
    REQUIRE(stored.getNmembers() == static_cast<int>(expected.size()));
    for (unsigned i = 0; i < expected.size(); ++i)
      CHECK(stored.getMemberName(i) == expected[i]);
    CHECK(stored.getSize() == sizeof(CampaignDayEvent));
  }
  CHECK(sameRows(readCampaignDays(file), first));

  // A later flush adds to the same dataset, after the earlier rows.
  const std::vector<CampaignDayEvent> second = {campaignDayRow(1, 0, 0, 41),
                                                campaignDayRow(1, 1, 3, 7)};
  for (const auto& r : second) logRow(logger, r);
  logger.flush(file, config, world, {});
  std::vector<CampaignDayEvent> all = first;
  all.insert(all.end(), second.begin(), second.end());
  CHECK(sameRows(readCampaignDays(file), all));
  CHECK(readCampaignDays(file)[2].allotment == kNoAllotment);
}

TEST_CASE("campaign registry lists the campaigns in config order") {
  ScopedTestFiles dir("june_campaign_days");
  const std::string file = dir.path("events.h5").string();
  const Config config = configWithCampaigns({"Zeta", "Alpha", "Mid"});

  EventLogger logger;
  logRow(logger, campaignDayRow(0, 0, 0, 5));
  logger.flush(file, config, WorldState{}, {});
  // Not alphabetical and not sorted by anything: the order of the config.
  CHECK(readRegistry(file, "/metadata/registries/campaigns") ==
        std::vector<std::string>{"Zeta", "Alpha", "Mid"});
}

TEST_CASE("campaign tier registry names each campaign's tiers by campaign_id") {
  ScopedTestFiles dir("june_campaign_days");
  const std::string file = dir.path("events.h5").string();
  Config config = configWithCampaigns({"Zeta", "Alpha"});
  for (const char* name : {"second", "first"}) {
    VaccinationTier tier;
    tier.name = name;
    config.vaccination.campaigns[0].tiers.push_back(tier);
  }
  VaccinationTier only;
  only.name = "only";
  config.vaccination.campaigns[1].tiers.push_back(only);

  EventLogger logger;
  logRow(logger, campaignDayRow(0, 0, 0, 5));
  logger.flush(file, config, WorldState{}, {});
  // tier_id indexes the list of the campaign named by campaign_id
  CHECK(readRegistry(file, "/metadata/registries/campaign_tiers/0") ==
        std::vector<std::string>{"second", "first"});
  CHECK(readRegistry(file, "/metadata/registries/campaign_tiers/1") ==
        std::vector<std::string>{"only"});
}

TEST_CASE("merging event files keeps campaign day rows from every file") {
  ScopedTestFiles dir("june_campaign_days");
  const std::string a = dir.path("a.h5").string();
  const std::string b = dir.path("b.h5").string();
  const Config config = configWithCampaigns({"First", "Second"});
  const WorldState world;

  const std::vector<CampaignDayEvent> rows_a = {campaignDayRow(0, 0, 0, 10),
                                                campaignDayRow(1, 0, 0, 11)};
  const std::vector<CampaignDayEvent> rows_b = {campaignDayRow(0, 1, 0, 20),
                                                campaignDayRow(0, 1, 1, 20),
                                                campaignDayRow(1, 1, 0, 21)};
  for (const auto& [file, rows] :
       {std::pair{a, rows_a}, std::pair{b, rows_b}}) {
    EventLogger logger;
    for (const auto& r : rows) logRow(logger, r);
    logger.flush(file, config, world, {});
  }

  const std::string merged = dir.path("merged.h5").string();
  mergeEventFiles({a, b}, merged);

  // Every row of every input, in input order.
  std::vector<CampaignDayEvent> expected = rows_a;
  expected.insert(expected.end(), rows_b.begin(), rows_b.end());
  CHECK(sameRows(readCampaignDays(merged), expected));
  // The registry comes through the merge.
  CHECK(readRegistry(merged, "/metadata/registries/campaigns") ==
        std::vector<std::string>{"First", "Second"});
}

TEST_CASE("clearing the logger drops buffered campaign days") {
  EventLogger logger;
  logRow(logger, campaignDayRow(0, 0, 0, 5));
  logRow(logger, campaignDayRow(1, 0, 0, 5));
  CHECK(logger.getTotalRecordCount() == 2);
  logger.clear();
  CHECK(logger.getTotalRecordCount() == 0);
}
