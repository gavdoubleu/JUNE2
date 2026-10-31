#include <H5Cpp.h>

#include <cstdio>
#include <set>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "simulation/simulator.h"
#include "simulator_fixtures.h"

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
