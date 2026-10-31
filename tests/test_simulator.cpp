#include <H5Cpp.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "simulation/simulator.h"
#include "simulator_fixtures.h"
#include "vaccination_manager_fixtures.h"

using namespace june;
using namespace june::simulator_fixture;

namespace {

// The person_id column of a dataset in an events file, in file order.
std::vector<int> readPersonIds(const std::string& file,
                               const std::string& dataset_path) {
  H5::H5File h5(file, H5F_ACC_RDONLY);
  H5::DataSet ds = h5.openDataSet(dataset_path);
  hsize_t rows = 0;
  ds.getSpace().getSimpleExtentDims(&rows);
  H5::CompType id_only(sizeof(int));
  id_only.insertMember("person_id", 0, H5::PredType::NATIVE_INT);
  std::vector<int> ids(rows);
  ds.read(ids.data(), id_only);
  return ids;
}

bool hasDataset(const std::string& file, const std::string& dataset_path) {
  H5::H5File h5(file, H5F_ACC_RDONLY);
  return H5Lexists(h5.getId(), dataset_path.c_str(), H5P_DEFAULT) > 0;
}

// Runs the small world to the end of January into `file`. Nobody keeps
// immunity, so people are infected again and again, and a one record buffer
// flushes after every slot. `lookup_mode` is save_full_person_details,
// `activities_mode` is save_person_activities and `contacts` is the daily
// contact rate (low means only the seeded people get infected).
void runFlushEverySlot(const std::string& file, const std::string& lookup_mode,
                       const std::string& activities_mode = "none",
                       double contacts = 5.0) {
  WorldState world;
  buildSmallWorld(world);
  Config config = makeSmallWorldConfig(world, "2020-01-31");
  config.contact_matrices.default_matrix->contacts = {{contacts}};
  config.simulation.max_event_buffer_size = 1;
  config.simulation.save_full_person_details = lookup_mode;
  config.simulation.save_person_activities = activities_mode;

  Simulator sim(world, config, nullptr,
                "tests/configs/simulator_infection_seeds.yaml", file);
  REQUIRE_NOTHROW(sim.run());
}

}  // namespace

TEST_CASE("Simulator Initialization") {
  WorldState world;
  Config config;
  config.simulation.start_date = "2020-01-01";
  config.simulation.end_date = "2020-01-02";  // 1 day simulation
  config.simulation.disease_file = "tests/configs/simulator_disease.yaml";

  // Minimal valid config setup for Simulator
  TimeSlot slot;
  slot.name = "Morning";
  slot.start = "08:00";
  slot.end = "12:00";

  config.schedule.day_type_cycle = {"workday"};
  config.schedule.day_type_names = {"workday"};

  ScheduleType st;
  st.slots_by_day_type["workday"].push_back(slot);
  config.schedule.schedule_types.push_back(st);

  // Flat default contact matrix so every disease mode has somewhere to
  // fall back to (finalizeDefaultModeMatrices requires this).
  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{0.2}};
  config.contact_matrices.default_matrix = default_contact_matrix;
  // This scenario declares no per-venue-type matrices; it means the default.
  config.contact_matrices.allow_default_matrix = true;

  // Initialise Simulator
  Simulator sim(world, config, nullptr,
                "tests/configs/simulator_infection_seeds.yaml", "test_sim.h5");

  CHECK(sim.getEventLogger() != nullptr);
}

TEST_CASE("Simulator run() multi-day smoke test") {
  WorldState world;
  buildSmallWorld(world);
  // 3 day simulation
  Config config = makeSmallWorldConfig(world, "2020-01-04");

  // Build and run simulator
  Simulator sim(world, config, nullptr,
                "tests/configs/simulator_infection_seeds.yaml",
                "test_sim_run.h5");

  // This is the critical test: run() should not crash
  REQUIRE_NOTHROW(sim.run());

  // Verify the event logger has some data (seeds should have been applied)
  EventLogger* logger = sim.getEventLogger();
  REQUIRE(logger != nullptr);

  // Clean up test output file
  std::remove("test_sim_run.h5");
}

TEST_CASE("Every infected person has exactly one /lookups/people row") {
  // Nobody keeps immunity, so people are infected again and again, and a
  // one record buffer flushes after every slot. A flush then often holds only
  // repeat infections of people who already have a row.
  const std::string out = "test_sim_lookup_rows.h5";
  runFlushEverySlot(out, "infected_only");

  std::vector<int> infections = readPersonIds(out, "/events/infections");
  std::vector<int> rows = readPersonIds(out, "/lookups/people");
  std::remove(out.c_str());

  const std::set<int> infected(infections.begin(), infections.end());
  // The run must contain repeat infections, or it checks nothing.
  REQUIRE(infections.size() > infected.size());
  CHECK(rows.size() == infected.size());
  CHECK(std::set<int>(rows.begin(), rows.end()) == infected);
}

TEST_CASE(
    "With save_full_person_details all, every person has exactly one "
    "/lookups/people row") {
  // The run flushes dozens of times, so this also covers the flushes after
  // the first write: they must not write anybody again.
  const std::string out = "test_sim_lookup_all.h5";
  runFlushEverySlot(out, "all");

  REQUIRE(hasDataset(out, "/lookups/people"));
  std::vector<int> rows = readPersonIds(out, "/lookups/people");
  std::remove(out.c_str());

  CHECK(rows.size() == 20);
  std::set<int> everyone;
  for (int id = 0; id < 20; ++id) everyone.insert(id);
  CHECK(std::set<int>(rows.begin(), rows.end()) == everyone);
}

TEST_CASE("With save_full_person_details none, no /lookups/people table") {
  const std::string out = "test_sim_lookup_none.h5";
  runFlushEverySlot(out, "none");

  CHECK_FALSE(hasDataset(out, "/lookups/people"));
  std::remove(out.c_str());
}

TEST_CASE(
    "With save_person_activities all, every person's activities are written "
    "once") {
  // Almost nobody is infected here, so "everyone" and "the infected" differ.
  const std::string out = "test_sim_activities_all.h5";
  runFlushEverySlot(out, "infected_only", "all", 0.0);

  REQUIRE(hasDataset(out, "/lookups/person_activities"));
  std::vector<int> rows = readPersonIds(out, "/lookups/person_activities");
  std::remove(out.c_str());

  // Everyone has two activities, primary_activity and residence.
  CHECK(rows.size() == 40);
  CHECK(std::set<int>(rows.begin(), rows.end()).size() == 20);
}

namespace {

struct DeathRow {
  int person_id;
  int venue_id;
  double time;
};

std::vector<DeathRow> readDeaths(const std::string& file) {
  H5::H5File h5(file, H5F_ACC_RDONLY);
  if (H5Lexists(h5.getId(), "/events/deaths", H5P_DEFAULT) <= 0) return {};
  H5::DataSet ds = h5.openDataSet("/events/deaths");
  hsize_t rows = 0;
  ds.getSpace().getSimpleExtentDims(&rows);
  H5::CompType type(sizeof(DeathRow));
  type.insertMember("person_id", HOFFSET(DeathRow, person_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("venue_id", HOFFSET(DeathRow, venue_id),
                    H5::PredType::NATIVE_INT);
  type.insertMember("time", HOFFSET(DeathRow, time),
                    H5::PredType::NATIVE_DOUBLE);
  std::vector<DeathRow> deaths(rows);
  ds.read(deaths.data(), type);
  return deaths;
}

// Everyone in the small world is infected at the start of day 0 and dies
// `days_to_death` days later, so the slot of the death is chosen by the caller.
void writeFatalDisease(const std::string& file, double days_to_death) {
  std::ofstream out(file);
  out << R"(disease:
  name: fatal_disease
  settings:
    fatality_stage:
      - name: dead_home
    recovered_stage:
      - name: recovered
  symptom_tags:
    - {name: recovered, value: -2}
    - {name: healthy, value: -1}
    - {name: exposed, value: 0}
    - {name: dead_home, value: 4}
  trajectories:
    - description: "exposed => dead at home"
      selection_key: home_ifr
      probability: 1.0
      severity: 1.0
      stages:
        - symptom_tag: exposed
          completion_time: {type: constant, value: )"
      << days_to_death << R"(}
        - symptom_tag: dead_home
          completion_time: {type: constant, value: 0.0}
  transmission:
    mode: Trajectory-Driven
    type: gamma
    max_infectiousness: {type: lognormal, s: 0.5, loc: 0.0, scale: 1}
    shape: {type: normal, loc: 1.56, scale: 0.08}
    rate: {type: normal, loc: 0.53, scale: 0.03}
    shift: {type: normal, loc: -2.12, scale: 0.1}
    modes:
      - {name: respiratory, mode_transmissibility_multiplier: 1.0}
  immunity: {level: 0.0, waning_rate: 0.001}
)";
}

void writeEveryoneSeed(const std::string& file) {
  std::ofstream out(file);
  out << R"(global_parameters:
  default_seed_strength: 1.0
  base_cases_per_capita: 1.0
infection_seeds:
  - name: "everyone"
    type: "uniform"
    date: "2020-01-01 08:00"
    parameters:
      cases_per_capita_multiplier: 1.0
)";
}

// The small world with three slots a day: two at the office and a night at
// home. Where someone is therefore says which slot they are in.
Config makeThreeSlotConfig(WorldState& world, const std::string& end_date,
                           const std::string& disease_file) {
  Config config;
  config.simulation.start_date = "2020-01-01";
  config.simulation.end_date = end_date;
  config.simulation.disease_file = disease_file;
  config.simulation.random_seed = 12345;

  ScheduleType st;
  st.name = "worker";
  for (const char* day_type : {"workday", "rest_day"}) {
    st.slots_by_day_type[day_type] = {
        {"Morning", "08:00", "12:00", {"primary_activity"}},
        {"Afternoon", "12:00", "20:00", {"primary_activity"}},
        {"Night", "20:00", "08:00", {"residence"}}};
    st.participation_by_day_type[day_type]["primary_activity"] = 1.0;
    st.participation_by_day_type[day_type]["residence"] = 1.0;
  }
  config.schedule.schedule_types.push_back(st);
  config.schedule.day_type_names = {"workday", "rest_day"};
  config.schedule.day_type_cycle = {"workday", "workday",  "workday", "workday",
                                    "workday", "rest_day", "rest_day"};

  ContactMatrix no_contact;
  no_contact.bins = {"all"};
  no_contact.contacts = {{0.0}};
  config.contact_matrices.default_matrix = no_contact;
  config.contact_matrices.allow_default_matrix = true;
  config.schedule.resolveSlots(world);
  return config;
}

// Runs the small world for `days` days in which everyone dies `days_to_death`
// days after the start. A campaign for everyone starts on day 2 when
// `vaccinate` is set. Returns the death rows and the people given a dose.
struct FatalRun {
  std::vector<DeathRow> deaths;
  std::vector<int> dosed;
  int dead_in_world = 0;
};

FatalRun runFatalWorld(double days_to_death, int days, bool vaccinate = false) {
  const std::string disease = "test_fatal_disease.yaml";
  const std::string seeds = "test_fatal_seeds.yaml";
  const std::string out = "test_fatal_run.h5";
  writeFatalDisease(disease, days_to_death);
  writeEveryoneSeed(seeds);

  WorldState world;
  buildSmallWorld(world);
  Config config = makeThreeSlotConfig(
      world, days == 2 ? "2020-01-03" : "2020-01-04", disease);
  if (vaccinate) {
    auto campaign = manager_fixture::vaxCampaign(
        "everyone", {manager_fixture::tierFor("all", 0)});
    campaign.start_date = "2020-01-03";  // day 2
    config.vaccination = manager_fixture::makeConfig({campaign}).vaccination;
  }

  FatalRun result;
  {
    Simulator sim(world, config, nullptr, seeds, out);
    REQUIRE_NOTHROW(sim.run());
  }
  result.deaths = readDeaths(out);
  if (hasDataset(out, "/events/vaccinations"))
    result.dosed = readPersonIds(out, "/events/vaccinations");
  for (const Person& p : world.people) result.dead_in_world += p.is_dead;
  std::remove(out.c_str());
  std::remove(disease.c_str());
  std::remove(seeds.c_str());
  return result;
}

}  // namespace

TEST_CASE("A death is applied in the slot it happens in") {
  // Death at 1.3: day 1, afternoon (the office). The death has to be applied
  // at the end of that slot, with the person where they died, not carried
  // into the night slot at home.
  const FatalRun run = runFatalWorld(1.3, 2);
  REQUIRE(run.deaths.size() == 20);
  for (const DeathRow& d : run.deaths) {
    CHECK(d.time == doctest::Approx(1.3));
    CHECK(d.venue_id == 0);
  }
}

TEST_CASE("A death in a day's last slot is applied before the next day") {
  // Death at 1.7: day 1, night. The run stops at the end of day 1, so a death
  // the next morning would still be pending.
  const FatalRun run = runFatalWorld(1.7, 2);
  CHECK(run.dead_in_world == 20);
  REQUIRE(run.deaths.size() == 20);
  for (const DeathRow& d : run.deaths) CHECK(d.time == doctest::Approx(1.7));
}

TEST_CASE("Nobody who died the night before is given a dose") {
  // Everyone dies on the night of day 1 and a campaign for everyone starts on
  // day 2, so there is nobody left to dose.
  const FatalRun run = runFatalWorld(1.7, 3, true);
  CHECK(run.deaths.size() == 20);
  CHECK(run.dosed.empty());
}
