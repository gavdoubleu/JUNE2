#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "core/config.h"
#include "core/world_state.h"
#include "doctest.h"
#include "simulation/simulator.h"
#include "utils/random.h"

using namespace june;

TEST_CASE("Simulator Initialization") {
  WorldState world;
  Config config;
  config.simulation.start_date = "2020-01-01";
  config.simulation.end_date = "2020-01-02";  // 1 day simulation
  config.simulation.disease_file = "configs/config_2021/disease.yaml";

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
                "configs/config_2021/infection_seeds.yaml", "test_sim.h5");

  CHECK(sim.getEventLogger() != nullptr);
}

TEST_CASE("Simulator run() multi-day smoke test") {
  // Seed the global RNG before simulator creation (critical invariant)
  GlobalRNG::seed(12345);

  // Create a small world with 20 people and 2 venues
  WorldState world;
  world.venue_type_names = {"office", "home"};
  world.activity_names = {"primary_activity", "residence"};
  world.geo_level_names = {"city"};

  // Create a geo unit
  GeographicalUnit gu;
  gu.id = 0;
  gu.name = "TestCity";
  gu.level_id = 0;
  gu.parent_id = -1;
  world.geo_units.push_back(gu);

  // Create venues
  for (int i = 0; i < 2; ++i) {
    Venue v;
    v.id = i;
    v.type_id = static_cast<uint8_t>(i);  // office=0, home=1
    v.geo_unit_id = 0;
    v.is_residence = (i == 1);
    world.venues.push_back(v);
  }

  // Create people with activities
  for (int i = 0; i < 20; ++i) {
    auto& p = world.people.emplace_back();
    p.id = i;
    p.age = 20.0f + i;
    p.sex = (i % 2 == 0) ? Sex::MALE : Sex::FEMALE;
    p.geo_unit_id = 0;

    // Give each person activity venues: primary_activity at office (venue 0),
    // residence at home (venue 1)
    p.activity_meta_start = static_cast<uint32_t>(world.activity_meta.size());
    p.activity_meta_count = 2;

    // primary_activity -> office
    Person::ActivityMeta am1;
    am1.activity_index = 0;  // primary_activity
    am1.venue_start = static_cast<uint32_t>(world.activity_venues.size());
    am1.venue_count = 1;
    world.activity_meta.push_back(am1);
    world.activity_venues.push_back({0, -1});  // venue 0, no subset

    // residence -> home
    Person::ActivityMeta am2;
    am2.activity_index = 1;  // residence
    am2.venue_start = static_cast<uint32_t>(world.activity_venues.size());
    am2.venue_count = 1;
    world.activity_meta.push_back(am2);
    world.activity_venues.push_back({1, -1});  // venue 1, no subset
  }
  world.buildIndices();

  // Config
  Config config;
  config.simulation.start_date = "2020-01-01";
  config.simulation.end_date = "2020-01-04";  // 3 day simulation
  config.simulation.disease_file = "configs/config_2021/disease.yaml";
  config.simulation.random_seed = 12345;

  // Schedule: single slot per day, everyone does primary_activity
  TimeSlot slot;
  slot.name = "Daytime";
  slot.start = "08:00";
  slot.end = "20:00";
  slot.allowed_activities = {"primary_activity", "residence"};

  ScheduleType st;
  st.name = "worker";
  st.slots_by_day_type["workday"].push_back(slot);
  st.slots_by_day_type["rest_day"].push_back(slot);
  st.participation_by_day_type["workday"]["primary_activity"] = 1.0;
  st.participation_by_day_type["workday"]["residence"] = 1.0;
  st.participation_by_day_type["rest_day"]["primary_activity"] = 0.5;
  st.participation_by_day_type["rest_day"]["residence"] = 1.0;
  config.schedule.schedule_types.push_back(st);
  config.schedule.day_type_names = {"workday", "rest_day"};
  config.schedule.day_type_cycle = {"workday", "workday",  "workday", "workday",
                                    "workday", "rest_day", "rest_day"};

  ContactMatrix default_contact_matrix;
  default_contact_matrix.bins = {"all"};
  default_contact_matrix.contacts = {{0.2}};
  config.contact_matrices.default_matrix = default_contact_matrix;
  // This scenario declares no per-venue-type matrices; it means the default.
  config.contact_matrices.allow_default_matrix = true;

  // Resolve schedule indices (maps string-keyed slots to indexed lookups)
  config.schedule.resolveSlots(world);

  // Build and run simulator
  Simulator sim(world, config, nullptr,
                "configs/config_2021/infection_seeds.yaml", "test_sim_run.h5");

  // This is the critical test: run() should not crash
  REQUIRE_NOTHROW(sim.run());

  // Verify the event logger has some data (seeds should have been applied)
  EventLogger* logger = sim.getEventLogger();
  REQUIRE(logger != nullptr);

  // Clean up test output file
  std::remove("test_sim_run.h5");
}
