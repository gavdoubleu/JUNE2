#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/world_state.h"
#include "doctest.h"
#include "epidemiology/vaccine.h"
#include "loaders/config_loader.h"
#include "test_utils.h"

using namespace june;

TEST_CASE("ConfigLoader - Vaccination Config from YAML") {
  VaccinationConfig config =
      ConfigLoader::loadVaccination("tests/configs/vaccines.yaml");

  SUBCASE("Pfizer vaccine is loaded correctly") {
    REQUIRE(config.vaccines.count("Pfizer") > 0);
    const auto& pfizer = config.vaccines.at("Pfizer");

    REQUIRE(pfizer.doses.size() == 3);

    // First dose: 52% infection efficacy for covid19
    REQUIRE(pfizer.doses[0].infection_efficacy.count("covid19") > 0);
    REQUIRE(!pfizer.doses[0].infection_efficacy.at("covid19").empty());
    CHECK(pfizer.doses[0].infection_efficacy.at("covid19")[0].efficacy ==
          doctest::Approx(0.52));
    CHECK(pfizer.doses[0].infection_efficacy.at("covid19")[0].min_age == 0);
    CHECK(pfizer.doses[0].infection_efficacy.at("covid19")[0].max_age == 100);

    // Second dose: 95% infection efficacy
    CHECK(pfizer.doses[1].infection_efficacy.at("covid19")[0].efficacy ==
          doctest::Approx(0.95));

    // Third dose (booster): 98% infection efficacy
    CHECK(pfizer.doses[2].infection_efficacy.at("covid19")[0].efficacy ==
          doctest::Approx(0.98));
  }

  SUBCASE("AstraZeneca vaccine is loaded correctly") {
    REQUIRE(config.vaccines.count("AstraZeneca") > 0);
    const auto& az = config.vaccines.at("AstraZeneca");

    REQUIRE(az.doses.size() == 3);
    CHECK(az.doses[0].infection_efficacy.at("covid19")[0].efficacy ==
          doctest::Approx(0.32));
    CHECK(az.doses[1].infection_efficacy.at("covid19")[0].efficacy ==
          doctest::Approx(0.75));
  }

  SUBCASE("Campaigns are loaded correctly") {
    REQUIRE(config.campaigns.size() == 3);

    const auto& pfizer = config.campaigns[0];
    CHECK(pfizer.vaccine_type == "Pfizer");
    REQUIRE(pfizer.dose_budget.has_value());
    CHECK(pfizer.dose_budget->doses == 300);
    CHECK(pfizer.dose_budget->per_days == 7);
    CHECK(*pfizer.coverage_target == doctest::Approx(0.5));
    REQUIRE(pfizer.tiers.size() == 2);
    CHECK(pfizer.tiers[0].dose == 1);
    CHECK(pfizer.tiers[1].dose == 0);
    REQUIRE(pfizer.tiers[1].selection_criteria.size() == 1);
    CHECK(pfizer.tiers[1].selection_criteria[0].property_path ==
          "activities.residence.venue_type");
    CHECK(pfizer.tiers[1].uptake.rows.size() == 2);

    CHECK(config.campaigns[1].vaccine_type == "AstraZeneca");
    CHECK_FALSE(config.campaigns[1].dose_budget.has_value());

    // Booster campaign
    const auto& booster = config.campaigns[2];
    REQUIRE(booster.tiers.size() == 1);
    CHECK(booster.tiers[0].dose == 2);
    CHECK(booster.courses_begun_with.size() == 2);
  }
}

namespace {

// A vaccines file whose first dose is complete and whose second dose has the
// given lines, so a failure on the second dose proves the index is reported.
std::string vaccineWithDose(const std::vector<std::string>& second_dose) {
  std::string yaml =
      "enabled: true\n"
      "vaccines:\n"
      "  TestVax:\n"
      "    doses:\n"
      "      - days_to_effective: 14\n"
      "        days_to_waning: 100\n"
      "        days_to_finished: 200\n"
      "      -\n";
  for (const auto& line : second_dose) yaml += "        " + line + "\n";
  return yaml;
}

}  // namespace

TEST_CASE("ConfigLoader - Vaccine dose onset") {
  ScopedTestFiles files("june_vaccine_onset");
  const auto load = [&](const std::vector<std::string>& second_dose) {
    const auto path =
        files.write("vaccines.yaml", vaccineWithDose(second_dose));
    return ConfigLoader::loadVaccination(path.string());
  };
  const std::vector<std::string> timing = {
      "days_to_effective: 14", "days_to_waning: 100", "days_to_finished: 200"};

  SUBCASE("days_to_onset is read from the dose") {
    auto fields = timing;
    fields.push_back("days_to_onset: 14");
    const auto config = load(fields);
    const auto& doses = config.vaccines.at("TestVax").doses;
    REQUIRE(doses.size() == 2);
    CHECK(doses[1].days_to_onset == 14.0);
    CHECK(doses[1].days_to_effective == 14.0);
  }

  SUBCASE("an absent days_to_onset means no delay") {
    const auto config = load(timing);
    const auto& doses = config.vaccines.at("TestVax").doses;
    REQUIRE(doses.size() == 2);
    CHECK(doses[0].days_to_onset == 0.0);
    CHECK(doses[1].days_to_onset == 0.0);
  }

  SUBCASE("onset between 0 and days_to_effective is accepted") {
    auto fields = timing;
    fields.push_back("days_to_onset: 7.5");
    CHECK(load(fields).vaccines.at("TestVax").doses[1].days_to_onset == 7.5);
  }

  SUBCASE("onset beyond days_to_effective throws and names the values") {
    auto fields = timing;
    fields.push_back("days_to_onset: 15");
    CHECK_THROWS_WITH(load(fields),
                      doctest::Contains("vaccine 'TestVax' dose 1"));
    CHECK_THROWS_WITH(load(fields), doctest::Contains("got 15"));
    CHECK_THROWS_WITH(load(fields), doctest::Contains("days_to_effective 14"));
  }

  SUBCASE("a negative onset throws and names the value") {
    auto fields = timing;
    fields.push_back("days_to_onset: -1");
    CHECK_THROWS_WITH(load(fields),
                      doctest::Contains("vaccine 'TestVax' dose 1"));
    CHECK_THROWS_WITH(load(fields), doctest::Contains("got -1"));
  }

  SUBCASE("each missing timing field throws, naming vaccine and dose") {
    for (size_t missing = 0; missing < timing.size(); ++missing) {
      auto fields = timing;
      const std::string key =
          fields[missing].substr(0, fields[missing].find(':'));
      fields.erase(fields.begin() + missing);
      CAPTURE(key);
      CHECK_THROWS_WITH(load(fields),
                        doctest::Contains("vaccine 'TestVax' dose 1"));
      const std::string quoted = "'" + key + "'";
      CHECK_THROWS_WITH(load(fields), doctest::Contains(quoted.c_str()));
    }
  }

  SUBCASE("waning_factor stays optional and defaults to 1") {
    CHECK(load(timing).vaccines.at("TestVax").doses[1].waning_factor == 1.0);
  }

  SUBCASE("a malformed file throws with its name instead of disabling") {
    const auto path =
        files.write("broken_vaccines.yaml", "enabled: true\nvaccines: [oops\n");
    CHECK_THROWS_WITH(ConfigLoader::loadVaccination(path.string()),
                      doctest::Contains("broken_vaccines.yaml"));
  }
}

TEST_CASE("ConfigLoader - Vaccine dose keys, age bands and efficacies") {
  ScopedTestFiles files("june_vaccine_dose_keys");
  const std::vector<std::string> timing = {
      "days_to_effective: 14", "days_to_waning: 100", "days_to_finished: 200"};
  // The second dose has the timing and then the given lines
  const auto load = [&](std::vector<std::string> extra) {
    auto lines = timing;
    lines.insert(lines.end(), extra.begin(), extra.end());
    return ConfigLoader::loadVaccination(
        files.write("vaccines.yaml", vaccineWithDose(lines)).string());
  };
  // What the loader says about those lines, or nothing if it accepts them
  const auto error = [&](std::vector<std::string> extra) -> std::string {
    try {
      load(std::move(extra));
    } catch (const std::runtime_error& e) {
      return e.what();
    }
    return "";
  };
  const auto names = [](const std::string& text,
                        std::initializer_list<const char*> fragments) {
    CHECK_MESSAGE(!text.empty(), "no error was thrown");
    for (const char* fragment : fragments)
      CHECK_MESSAGE(text.find(fragment) != std::string::npos,
                    "'" << fragment << "' is missing from: " << text);
  };

  SUBCASE("a misspelt dose key is an error, not a default") {
    names(error({"days_to_onsett: 14"}),
          {"vaccine 'TestVax' dose 1", "unknown field 'days_to_onsett'",
           "days_to_onset"});
  }

  SUBCASE("an age band must be min-max in whole years") {
    for (const char* band : {"0_17", "17-0", "a-b", "-5-10", "5-", "10"}) {
      CAPTURE(band);
      names(error({"infection_efficacy:",
                   std::string("  covid19: {\"") + band + "\": 0.5}"}),
            {"vaccine 'TestVax' dose 1 infection_efficacy covid19",
             (std::string("band '") + band + "'").c_str(), "min-max"});
    }
  }

  SUBCASE("an efficacy is a number from 0 to 1") {
    names(error({"symptom_efficacy:", "  covid19: 1.5"}),
          {"dose 1 symptom_efficacy covid19 efficacy",
           "must be a number from 0 to 1", "got '1.5'"});
    names(error({"infection_efficacy:", "  covid19: {\"0-17\": -0.1}"}),
          {"covid19 band '0-17' efficacy", "got '-0.1'"});
    names(error({"infection_efficacy:", "  covid19: high"}), {"got 'high'"});
  }

  SUBCASE("0, 1 and a set of bands are accepted") {
    const auto config =
        load({"infection_efficacy:", "  covid19: {\"0-17\": 0, \"18-100\": 1}",
              "symptom_efficacy:", "  covid19: 0.7"});
    const auto& dose = config.vaccines.at("TestVax").doses[1];
    const auto& bands = dose.infection_efficacy.at("covid19");
    REQUIRE(bands.size() == 2);
    CHECK(bands[0].efficacy == 0.0);
    CHECK(bands[1].min_age == 18);
    CHECK(bands[1].max_age == 100);
    CHECK(bands[1].efficacy == 1.0);
    CHECK(dose.symptom_efficacy.at("covid19")[0].efficacy == 0.7);
  }
}

TEST_CASE("ConfigLoader - Contact Matrix Config from YAML") {
  ContactMatrixConfig cm =
      ConfigLoader::loadContactMatrices("tests/configs/contact_matrices.yaml");

  SUBCASE("Venue matrices are loaded") {
    REQUIRE(cm.matrices.count("household") > 0);
    const auto& hh = cm.matrices.at("household");
    CHECK(hh.bins.size() == 1);
    CHECK(hh.bins[0] == "residents");
  }

  SUBCASE("resolution builds a per-venue-type lookup") {
    WorldState world;
    world.venue_type_names = {"office", "household", "pub"};
    world.buildIndices();
    // "pub" has no matrix in this fixture and can only be filled from the
    // default; accept that here so the lookup itself is what is under test.
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world);

    CHECK(cm.getBinStructure(world.getVenueTypeIndex("household")).bins[0] ==
          "residents");
    CHECK_FALSE(
        cm.getBinStructure(world.getVenueTypeIndex("office")).bins.empty());
  }

  SUBCASE("a venue type with no matrix of its own is refused at load") {
    WorldState world;
    world.venue_type_names = {"office", "household", "pub"};
    world.buildIndices();
    CHECK_THROWS_WITH_AS(finalizeContactMatrices(cm, world),
                         doctest::Contains("pub"), std::runtime_error);
  }
}

TEST_CASE("ContactMatrix beta scales contacts at load time") {
  ContactMatrixConfig cm =
      ConfigLoader::loadContactMatrices("tests/configs/contact_matrices.yaml");

  SUBCASE("single-mode venue with beta=3.0 triples contacts") {
    REQUIRE(cm.matrices.count("household") > 0);
    const auto& household_matrix = cm.matrices.at("household");
    // Raw YAML value is 2.5; beta=3.0 → effective = 7.5
    CHECK(household_matrix.contacts[0][0] == doctest::Approx(7.5));
  }

  SUBCASE("single-mode venue with no beta leaves contacts unchanged") {
    REQUIRE(cm.matrices.count("office") > 0);
    const auto& office_matrix = cm.matrices.at("office");
    CHECK(office_matrix.contacts[0][0] == doctest::Approx(5.0));
  }

  SUBCASE("multi-mode venue: respiratory beta=2.0 doubles contacts") {
    REQUIRE(cm.mode_matrices.count("hospital") > 0);
    REQUIRE(cm.mode_matrices.at("hospital").count("respiratory") > 0);
    const auto& respiratory_matrix =
        cm.mode_matrices.at("hospital").at("respiratory");
    // Raw YAML value is 1.5; beta=2.0 → effective = 3.0
    CHECK(respiratory_matrix.contacts[0][0] == doctest::Approx(3.0));
  }

  SUBCASE("multi-mode venue: mode with no beta leaves contacts unchanged") {
    REQUIRE(cm.mode_matrices.count("hospital") > 0);
    REQUIRE(cm.mode_matrices.at("hospital").count("physical_contact") > 0);
    const auto& physical_contact_matrix =
        cm.mode_matrices.at("hospital").at("physical_contact");
    CHECK(physical_contact_matrix.contacts[0][0] == doctest::Approx(0.5));
  }

  SUBCASE("default_contacts_matrix beta scales contacts") {
    REQUIRE(cm.default_matrix.has_value());
    // respiratory mode: raw=0.004, beta=1.5 → effective=0.006
    CHECK(cm.default_matrix->contacts[0][0] == doctest::Approx(0.006));
  }
}

TEST_CASE("ConfigLoader - default_contacts_matrix modes: form") {
  SUBCASE("modes-format default loads per-mode, covering every disease mode") {
    ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
        "tests/configs/contact_matrices_modes_default.yaml");

    REQUIRE(cm.default_mode_matrices.has_value());
    REQUIRE(cm.default_mode_matrices->count("respiratory") > 0);
    REQUIRE(cm.default_mode_matrices->count("physical_contact") > 0);
    // raw=0.004, beta=1.5 -> effective=0.006
    CHECK(cm.default_mode_matrices->at("respiratory").contacts[0][0] ==
          doctest::Approx(0.006));
    // no beta -> unchanged
    CHECK(cm.default_mode_matrices->at("physical_contact").contacts[0][0] ==
          doctest::Approx(0.006));

    WorldState world;
    world.venue_type_names = {"hospital", "gym"};
    world.buildIndices();
    // gym declares no matrix, so it can only come from the per-mode default.
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, {"respiratory", "physical_contact"});

    uint8_t gym_id = world.getVenueTypeIndex("gym");
    CHECK(cm.getMatrix(gym_id, 0).contacts[0][0] == doctest::Approx(0.006));
  }

  SUBCASE("missing default_contacts_matrix throws") {
    CHECK_THROWS_AS(ConfigLoader::loadContactMatrices(
                        "tests/configs/contact_matrices_missing_default.yaml"),
                    std::runtime_error);
  }

  SUBCASE("per-mode default missing a disease mode throws") {
    CHECK_THROWS_AS(
        ConfigLoader::loadContactMatrices(
            "tests/configs/contact_matrices_default_missing_mode.yaml"),
        std::runtime_error);
  }
}

TEST_CASE("ContactMatrixConfig::finalizeDefaultModeMatrices") {
  SUBCASE(
      "per-mode default reachable when contact_matrices declares no modes") {
    ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
        "tests/configs/contact_matrices_empty_mode_names_default.yaml");
    REQUIRE(cm.mode_names.empty());
    REQUIRE(cm.default_mode_matrices.has_value());

    WorldState world;
    world.venue_type_names = {"household", "gym"};
    world.buildIndices();
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, {"respiratory", "physical_contact"});

    uint8_t gym_id = world.getVenueTypeIndex("gym");
    // raw=0.004, beta=1.5 -> effective=0.006
    CHECK(cm.getMatrix(gym_id, 0).contacts[0][0] == doctest::Approx(0.006));
    CHECK(cm.getMatrix(gym_id, 1).contacts[0][0] == doctest::Approx(0.006));
  }

  SUBCASE(
      "disease mode absent from default matrix and no flat fallback throws") {
    ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
        "tests/configs/contact_matrices_empty_mode_names_default.yaml");

    WorldState world;
    world.venue_type_names = {"household", "gym"};
    world.buildIndices();
    cm.resolve(world);

    CHECK_THROWS_AS(cm.finalizeDefaultModeMatrices(
                        world, {"respiratory", "physical_contact", "fomite"}),
                    std::runtime_error);
  }
}

TEST_CASE("ContactMatrixConfig::finalizeDiseaseModeAlignment") {
  SUBCASE(
      "matrices matched by name even when mode order differs from disease") {
    // Fixture declares hospital's modes as physical_contact, respiratory (in
    // that order) -> cm.mode_names == ["physical_contact", "respiratory"].
    ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
        "tests/configs/contact_matrices_reversed_modes.yaml");
    REQUIRE(cm.mode_names.size() == 2);
    CHECK(cm.mode_names[0] == "physical_contact");
    CHECK(cm.mode_names[1] == "respiratory");

    WorldState world;
    world.venue_type_names = {"hospital"};
    world.buildIndices();

    // Disease declares the opposite order: respiratory=0, physical_contact=1.
    finalizeContactMatrices(cm, world, {"respiratory", "physical_contact"});

    uint8_t hospital_id = world.getVenueTypeIndex("hospital");
    // respiratory: raw 1.5, beta 2.0 -> 3.0
    CHECK(cm.getMatrix(hospital_id, 0).contacts[0][0] == doctest::Approx(3.0));
    // physical_contact: no beta -> unchanged
    CHECK(cm.getMatrix(hospital_id, 1).contacts[0][0] == doctest::Approx(0.5));
  }

  SUBCASE("a disease mode no venue type declares is reported by name") {
    // "school" has no entry at all in the fixture, and "fomite" has no
    // dedicated entry anywhere, so this pair can only come from the default.
    // The load names the pair rather than quietly substituting.
    WorldState world;
    world.venue_type_names = {"school"};
    world.buildIndices();

    ContactMatrixConfig strict = ConfigLoader::loadContactMatrices(
        "tests/configs/contact_matrices.yaml");
    CHECK_THROWS_WITH_AS(
        finalizeContactMatrices(strict, world, {"respiratory", "fomite"}),
        doctest::Contains("fomite"), std::runtime_error);

    ContactMatrixConfig opted_in = ConfigLoader::loadContactMatrices(
        "tests/configs/contact_matrices.yaml");
    opted_in.allow_default_matrix = true;
    finalizeContactMatrices(opted_in, world, {"respiratory", "fomite"});
    uint8_t school_id = world.getVenueTypeIndex("school");
    // Resolves to the flat default_contacts_matrix: raw 0.004, beta 1.5.
    CHECK(opted_in.getMatrix(school_id, 1).contacts[0][0] ==
          doctest::Approx(0.006));
  }

  SUBCASE("orphaned contact-matrix mode name warns but does not throw") {
    ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
        "tests/configs/contact_matrices.yaml");  // mode_names: [respiratory,
                                                 // physical_contact]

    WorldState world;
    world.venue_type_names = {"hospital"};
    world.buildIndices();
    cm.resolve(world);

    std::ostringstream captured;
    std::streambuf* old_cerr = std::cerr.rdbuf(captured.rdbuf());
    CHECK_NOTHROW(cm.finalizeDiseaseModeAlignment({"respiratory"}));
    std::cerr.rdbuf(old_cerr);

    CHECK(captured.str().find("physical_contact") != std::string::npos);
  }

  SUBCASE("duplicate disease mode names throw instead of aliasing") {
    ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
        "tests/configs/contact_matrices.yaml");  // mode_names: [respiratory,
                                                 // physical_contact]

    WorldState world;
    world.venue_type_names = {"hospital"};
    world.buildIndices();
    cm.resolve(world);

    CHECK_THROWS_AS(
        cm.finalizeDiseaseModeAlignment({"respiratory", "respiratory"}),
        std::runtime_error);
  }

  SUBCASE("a disease with no declared modes still resolves one channel") {
    // Nothing forces a disease to name its modes. The FOI loop still asks for
    // mode 0, so a nameless single channel has to resolve to something real.
    ContactMatrixConfig cm = ConfigLoader::loadContactMatrices(
        "tests/configs/contact_matrices.yaml");

    WorldState world;
    world.venue_type_names = {"hospital"};
    world.buildIndices();
    cm.allow_default_matrix = true;
    finalizeContactMatrices(cm, world, {});

    uint8_t hospital_id = world.getVenueTypeIndex("hospital");
    CHECK_FALSE(cm.getMatrix(hospital_id, 0).contacts.empty());
  }
}

TEST_CASE("SimulationConfig - calendar event paths parsed from config_paths") {
  SUBCASE("both paths present") {
    ScopedTestFiles files("june_config_calendar_test");
    const auto yaml_path =
        files.write("simulation.yaml",
                    "time:\n"
                    "  start_date: \"2020-01-01\"\n"
                    "  end_date: \"2020-01-10\"\n"
                    "config_paths:\n"
                    "  calendar_events_file: \"data/calendar_events.csv\"\n"
                    "  calendar_event_catchment_rules_file: \""
                    "data/catchment_rules.csv\"\n");
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path.string());
    CHECK(cfg.calendar_events_file == "data/calendar_events.csv");
    CHECK(cfg.calendar_event_catchment_rules_file ==
          "data/catchment_rules.csv");
  }

  SUBCASE("absent keys leave fields empty, no error") {
    ScopedTestFiles files("june_config_no_calendar_test");
    const auto yaml_path = files.write("simulation.yaml",
                                       "time:\n"
                                       "  start_date: \"2020-01-01\"\n"
                                       "  end_date: \"2020-01-10\"\n"
                                       "config_paths:\n"
                                       "  disease_file: \"disease.yaml\"\n");
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path.string());
    CHECK(cfg.calendar_events_file == "");
    CHECK(cfg.calendar_event_catchment_rules_file == "");
  }
}

TEST_CASE("SimulationConfig - save_coordinated_encounters parsed from output") {
  SUBCASE("absent defaults to false") {
    ScopedTestFiles files("june_config_output_test");
    const auto yaml_path = files.write("simulation.yaml",
                                       "time:\n"
                                       "  start_date: \"2020-01-01\"\n"
                                       "  end_date: \"2020-01-10\"\n"
                                       "config_paths: {}\n"
                                       "output:\n"
                                       "  stats_interval_days: 1\n");
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path.string());
    CHECK(cfg.save_coordinated_encounters == false);
  }

  SUBCASE("explicit true is honoured") {
    ScopedTestFiles files("june_config_output_true_test");
    const auto yaml_path = files.write("simulation.yaml",
                                       "time:\n"
                                       "  start_date: \"2020-01-01\"\n"
                                       "  end_date: \"2020-01-10\"\n"
                                       "config_paths: {}\n"
                                       "output:\n"
                                       "  stats_interval_days: 1\n"
                                       "  save_coordinated_encounters: true\n");
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path.string());
    CHECK(cfg.save_coordinated_encounters == true);
  }

  SUBCASE("explicit false is honoured") {
    ScopedTestFiles files("june_config_output_false_test");
    const auto yaml_path =
        files.write("simulation.yaml",
                    "time:\n"
                    "  start_date: \"2020-01-01\"\n"
                    "  end_date: \"2020-01-10\"\n"
                    "config_paths: {}\n"
                    "output:\n"
                    "  stats_interval_days: 1\n"
                    "  save_coordinated_encounters: false\n");
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path.string());
    CHECK(cfg.save_coordinated_encounters == false);
  }
}

TEST_CASE("SimulationConfig - lookup table modes are validated") {
  const std::string yaml_head =
      "time:\n"
      "  start_date: \"2020-01-01\"\n"
      "  end_date: \"2020-01-10\"\n"
      "config_paths: {}\n"
      "output:\n"
      "  stats_interval_days: 1\n";
  // Loads a simulation.yaml whose output block sets `key` to `value`.
  auto load = [&](const std::string& key, const std::string& value) {
    ScopedTestFiles files("june_config_lookup_mode_test");
    const auto yaml_path = files.write(
        "simulation.yaml", yaml_head + "  " + key + ": \"" + value + "\"\n");
    return ConfigLoader::loadSimulation(yaml_path.string());
  };

  SUBCASE("every valid value loads") {
    for (const char* value : {"none", "all", "infected_only"}) {
      CHECK(load("save_full_person_details", value).save_full_person_details ==
            value);
      CHECK(load("save_person_activities", value).save_person_activities ==
            value);
    }
  }

  SUBCASE("absent keys keep their defaults") {
    ScopedTestFiles files("june_config_lookup_mode_default_test");
    const auto yaml_path = files.write("simulation.yaml", yaml_head);
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path.string());
    CHECK(cfg.save_full_person_details == "infected_only");
    CHECK(cfg.save_person_activities == "none");
  }

  SUBCASE("an unknown value throws and says what is valid") {
    CHECK_THROWS_WITH_AS(
        load("save_full_person_details", "infected-only"),
        "output.save_full_person_details is 'infected-only'; valid values "
        "are 'none', 'all' and 'infected_only'.",
        std::runtime_error);
    CHECK_THROWS_WITH_AS(
        load("save_person_activities", "All"),
        "output.save_person_activities is 'All'; valid values are 'none', "
        "'all' and 'infected_only'.",
        std::runtime_error);
  }
}

TEST_CASE("ConfigLoader - Optional Files Handling") {
  SUBCASE("No vaccines file configured disables the module") {
    // Vaccines is the only truly-optional sub-config: leaving the path out of
    // config_paths means "don't run this module", which is fine, the rest of
    // the simulator can proceed.
    VaccinationConfig vc = ConfigLoader::loadVaccination("");
    CHECK(vc.enabled == false);
  }

  SUBCASE("A vaccines file that is named but missing is a misconfiguration") {
    CHECK_THROWS_WITH_AS(
        ConfigLoader::loadVaccination("non_existent_vaccines.yaml"),
        doctest::Contains("non_existent_vaccines.yaml"), std::runtime_error);
  }

  SUBCASE("Missing coordinated_encounters file is a misconfiguration") {
    // simulation.yaml's coordinated_encounters_file points at this YAML;
    // the required file must raise instead of silently disabling the feature.
    CHECK_THROWS_AS(
        ConfigLoader::loadCoordinatedEncounters("non_existent_ce.yaml"),
        std::runtime_error);
  }

  SUBCASE("Corrupted coordinated_encounters file must throw") {
    // A missing required field must throw. The loader must not catch the error
    // and disable the feature.
    ScopedTestFiles files("june_config_bad_ce_test");
    const auto bad_file = files.write(
        "bad.yaml",
        "coordinated_encounters:\n  encounters:\n    - name: bad\n      "
        "network: missing_fields\n");

    CHECK_THROWS_AS(ConfigLoader::loadCoordinatedEncounters(bad_file.string()),
                    std::runtime_error);
  }
}

namespace {

// Two doses, so a tier can ask for dose 0 or 1 and nothing above.
constexpr const char* kTwoDoseVaccine =
    "enabled: true\n"
    "vaccines:\n"
    "  TestVax:\n"
    "    doses:\n"
    "      - {days_to_effective: 14, days_to_waning: 100, days_to_finished: "
    "200}\n"
    "      - {days_to_effective: 14, days_to_waning: 100, days_to_finished: "
    "200}\n"
    "campaigns:\n";

// A campaign called Rollout whose body is the given lines, each already
// indented by four spaces.
std::string rolloutFile(const std::string& body) {
  return std::string(kTwoDoseVaccine) + "  Rollout:\n" + body;
}

// Runs the call, which must throw a runtime_error naming every fragment.
template <typename F>
void checkNames(F&& call, std::initializer_list<const char*> fragments) {
  std::string message;
  try {
    call();
  } catch (const std::runtime_error& e) {
    message = e.what();
  }
  CHECK_MESSAGE(!message.empty(), "no error was thrown");
  for (const char* fragment : fragments) {
    CHECK_MESSAGE(
        message.find(fragment) != std::string::npos,
        "'" << std::string(fragment) << "' is missing from: " << message);
  }
}

// Ages 5, 10, 17, 18, 30, 70, and a second activity so its index is not 0
WorldState campaignTestWorld() {
  WorldState world = TestWorldFactory::createMinimalWorld(6);
  const float ages[] = {5, 10, 17, 18, 30, 70};
  for (size_t i = 0; i < 6; ++i) world.people[i].age = ages[i];
  world.activity_names = {"work", "residence"};
  return world;
}

class CampaignFiles {
 public:
  VaccinationConfig load(const std::string& body) {
    return ConfigLoader::loadVaccination(
        files.write("vaccines.yaml", rolloutFile(body)).string());
  }
  // The first tier of a campaign whose tier body is the given lines
  VaccinationConfig loadTier(const std::string& tier_lines) {
    return load("    vaccine_type: TestVax\n    tiers:\n" + tier_lines);
  }
  // Loads, resolves against the test world and returns the first tier
  VaccinationTier resolvedTier(const std::string& tier_lines) {
    auto config = loadTier(tier_lines);
    const WorldState world = campaignTestWorld();
    config.resolve(world);
    return config.campaigns.at(0).tiers.at(0);
  }
  ScopedTestFiles files{"june_campaign_config"};
};

}  // namespace

TEST_CASE("ConfigLoader - a spec shaped campaign parses to its values") {
  CampaignFiles f;
  f.files.write("uptake_young.csv", "filter.age,probability\n<18,0.8\n");
  const auto config = f.load(R"(    vaccine_type: TestVax
    dose_budget: {doses: 10000000, per_days: 30}
    coverage_target: 0.7
    tiers:
      - name: young_second
        dose: 1
        selection:
          - {property: "age", operator: "<", value: 18}
        uptake: 0.9
        decision_unit: person
      - name: young_first
        dose: 0
        selection: [{property: "age", operator: "<", value: 18}]
        uptake: {file: "uptake_young.csv"}
        decision_unit:
          shared_venue: residence
          scope:    [{property: "activities.residence.venue_type", operator: "==", value: "family_home"}]
          deciders: [{property: "age", operator: ">=", value: 18}]
          combine: max
          if_no_decider: individual
      - name: everyone_else
        dose: 0
)");
  REQUIRE(config.campaigns.size() == 1);
  const auto& campaign = config.campaigns[0];
  CHECK(campaign.name == "Rollout");
  REQUIRE(campaign.dose_budget.has_value());
  CHECK(campaign.dose_budget->doses == 10000000);
  CHECK(campaign.dose_budget->per_days == 30);
  REQUIRE(campaign.coverage_target.has_value());
  CHECK(*campaign.coverage_target == doctest::Approx(0.7));

  REQUIRE(campaign.tiers.size() == 3);
  const auto& second = campaign.tiers[0];
  CHECK(second.name == "young_second");
  CHECK(second.dose == 1);
  REQUIRE(second.selection_criteria.size() == 1);
  CHECK(second.selection_criteria[0].property_path == "age");
  CHECK(second.selection_criteria[0].operator_type == "<");
  CHECK(std::get<int>(second.selection_criteria[0].value) == 18);
  REQUIRE(second.uptake.rows.size() == 1);
  CHECK(second.uptake.rows[0].probability == doctest::Approx(0.9));
  CHECK(second.uptake.rows[0].criteria.empty());
  CHECK(second.uptake.csv_path.empty());
  CHECK_FALSE(second.decision_unit.shared());

  const auto& first = campaign.tiers[1];
  CHECK(first.name == "young_first");
  CHECK(first.dose == 0);
  CHECK(first.uptake.rows.empty());
  CHECK(first.uptake.csv_path == f.files.path("uptake_young.csv").string());
  CHECK(first.decision_unit.shared_activity == "residence");
  CHECK(first.decision_rules.scope.size() == 1);
  CHECK(first.decision_rules.deciders.size() == 1);
  CHECK(first.decision_rules.combine == DeciderCombine::Max);
  CHECK(first.decision_rules.if_no_decider == NoDecider::Individual);

  // Absent keys mean everyone: no selection, uptake or shared decision
  const auto& rest = campaign.tiers[2];
  CHECK(rest.selection_criteria.empty());
  CHECK(rest.uptake.rows.empty());
  CHECK(rest.uptake.csv_path.empty());
  CHECK_FALSE(rest.decision_unit.shared());
}

TEST_CASE("ConfigLoader - a shared decision unit defaults to mean and error") {
  CampaignFiles f;
  const auto tier = f.loadTier(R"(      - name: t
        dose: 0
        decision_unit: {shared_venue: residence, scope: [], deciders: []}
)")
                        .campaigns[0]
                        .tiers[0];
  CHECK(tier.decision_rules.combine == DeciderCombine::Mean);
  CHECK(tier.decision_rules.if_no_decider == NoDecider::Error);

  const auto min = f.loadTier(R"(      - name: t
        dose: 0
        decision_unit: {shared_venue: residence, scope: [], deciders: [], combine: min}
)")
                       .campaigns[0]
                       .tiers[0];
  CHECK(min.decision_rules.combine == DeciderCombine::Min);
}

TEST_CASE("ConfigLoader - a campaign keeps its window, gaps and brand filter") {
  CampaignFiles f;
  const auto config = f.load(R"(    start_date: "2024-01-01"
    end_date: "2024-01-20"
    days_to_next_dose: [0, 14]
    vaccine_type: TestVax
    courses_begun_with: [TestVax]
    tiers: [{name: t, dose: 0}]
)");
  REQUIRE(config.campaigns.size() == 1);
  const auto& c = config.campaigns[0];
  CHECK(c.start_date == "2024-01-01");
  CHECK(c.end_date == "2024-01-20");
  CHECK(c.days_to_next_dose == std::vector<double>{0, 14});
  CHECK(c.vaccine_type == "TestVax");
  CHECK(c.courses_begun_with == std::vector<std::string>{"TestVax"});
  CHECK_FALSE(c.dose_budget.has_value());
  CHECK_FALSE(c.coverage_target.has_value());
  CHECK(c.tiers.size() == 1);
}

TEST_CASE("ConfigLoader - campaign validation names the campaign and value") {
  CampaignFiles f;
  const std::string head = "    vaccine_type: TestVax\n";
  const std::string one_tier = "    tiers: [{name: t, dose: 0}]\n";

  SUBCASE("dose_budget must give whole numbers above zero") {
    checkNames(
        [&] { f.load(head + "    dose_budget: {doses: 1.5, per_days: 30}\n"); },
        {"Rollout", "dose_budget", "doses", "whole number", "1.5"});
    checkNames(
        [&] {
          f.load(head + "    dose_budget: {doses: lots, per_days: 30}\n");
        },
        {"Rollout", "doses", "whole number", "lots"});
    checkNames(
        [&] { f.load(head + "    dose_budget: {doses: 10, per_days: 2.5}\n"); },
        {"Rollout", "per_days", "whole number", "2.5"});
    checkNames(
        [&] { f.load(head + "    dose_budget: {doses: 0, per_days: 30}\n"); },
        {"Rollout", "doses", "greater than zero", "got 0"});
    checkNames(
        [&] { f.load(head + "    dose_budget: {doses: 10, per_days: -3}\n"); },
        {"Rollout", "per_days", "greater than zero", "got -3"});
  }

  SUBCASE("dose_budget needs both numbers and no other key") {
    checkNames([&] { f.load(head + "    dose_budget: {doses: 10}\n"); },
               {"Rollout", "dose_budget", "needs 'per_days'"});
    checkNames([&] { f.load(head + "    dose_budget: {per_days: 10}\n"); },
               {"Rollout", "needs 'doses'"});
    checkNames(
        [&] {
          f.load(head +
                 "    dose_budget: {doses: 1, per_days: 1, per_week: 2}\n");
        },
        {"Rollout", "dose_budget", "unknown field 'per_week'"});
    checkNames([&] { f.load(head + "    dose_budget: 5\n"); },
               {"Rollout", "dose_budget", "mapping", "'5'"});
  }

  SUBCASE("a valid dose_budget is accepted") {
    CHECK(f.load(head + one_tier + "    dose_budget: {doses: 1, per_days: 1}\n")
              .campaigns[0]
              .dose_budget->doses == 1);
  }

  SUBCASE("coverage_target is above 0 and at most 1") {
    for (const char* bad : {"0", "-0.2", "1.5", "most"}) {
      checkNames([&] { f.load(head + "    coverage_target: " + bad + "\n"); },
                 {"Rollout", "coverage_target", bad});
    }
    CHECK(*f.load(head + one_tier + "    coverage_target: 1\n")
               .campaigns[0]
               .coverage_target == 1.0);
  }

  SUBCASE("a typo in a key cannot silently mean unlimited") {
    checkNames(
        [&] { f.load(head + "    dose_budjet: {doses: 1, per_days: 1}\n"); },
        {"Rollout", "unknown field 'dose_budjet'", "dose_budget"});
  }

  SUBCASE("a campaign without tiers is refused by name") {
    checkNames([&] { f.load(head); }, {"Rollout", "needs 'tiers'"});
    checkNames(
        [&] { f.load(head + "    dose_budget: {doses: 1, per_days: 1}\n"); },
        {"Rollout", "needs 'tiers'"});
  }

  SUBCASE("the keys of the old daily roll are unknown fields") {
    for (const char* old : {"daily_coverage: 0.5", "dose_sequence: [0]",
                            "selection: [{property: age, operator: '>', "
                            "value: 1}]"}) {
      checkNames([&] { f.load(head + one_tier + "    " + old + "\n"); },
                 {"Rollout", "unknown field", "tiers"});
    }
  }

  SUBCASE("more than 65535 campaigns are refused") {
    const auto file = [&](int count) {
      std::string yaml = kTwoDoseVaccine;
      for (int i = 0; i < count; ++i)
        yaml += "  c" + std::to_string(i) +
                ": {vaccine_type: TestVax, tiers: [{name: t, dose: 0}]}\n";
      return f.files.write("many.yaml", yaml).string();
    };
    CHECK(ConfigLoader::loadVaccination(file(65535)).campaigns.size() == 65535);
    checkNames([&] { ConfigLoader::loadVaccination(file(65536)); },
               {"65536 campaigns", "at most 65535"});
  }
}

TEST_CASE("ConfigLoader - tier validation names campaign, tier and value") {
  CampaignFiles f;

  SUBCASE("tiers must be a non-empty list") {
    checkNames([&] { f.load("    vaccine_type: TestVax\n    tiers: []\n"); },
               {"Rollout", "tiers", "non-empty"});
    checkNames(
        [&] { f.load("    vaccine_type: TestVax\n    tiers: {name: t}\n"); },
        {"Rollout", "tiers", "non-empty list"});
  }

  SUBCASE("a tier needs a name, unique in its campaign") {
    checkNames([&] { f.loadTier("      - {dose: 0}\n"); },
               {"Rollout", "tier 1", "needs a 'name'"});
    checkNames(
        [&] {
          f.loadTier(
              "      - {name: a, dose: 0}\n      - {name: a, dose: 1}\n");
        },
        {"Rollout", "two tiers named 'a'"});
  }

  SUBCASE("the same name in two campaigns is fine") {
    const auto config = f.load(
        "    vaccine_type: TestVax\n    tiers: [{name: a, dose: 0}]\n"
        "  Second:\n    vaccine_type: TestVax\n    tiers: [{name: a, dose: "
        "0}]\n");
    CHECK(config.campaigns.size() == 2);
  }

  SUBCASE("a tier needs a dose of 0 or more, as a whole number") {
    checkNames([&] { f.loadTier("      - {name: a}\n"); },
               {"Rollout", "tier 'a'", "needs a 'dose'"});
    checkNames([&] { f.loadTier("      - {name: a, dose: -1}\n"); },
               {"Rollout", "tier 'a'", "dose", "-1"});
    checkNames([&] { f.loadTier("      - {name: a, dose: 0.5}\n"); },
               {"Rollout", "tier 'a'", "dose", "whole number", "0.5"});
    checkNames([&] { f.loadTier("      - {name: a, dose: first}\n"); },
               {"Rollout", "tier 'a'", "dose", "whole number", "first"});
  }

  SUBCASE("the dose must exist in the campaign's vaccine") {
    CHECK(
        f.loadTier("      - {name: a, dose: 1}\n").campaigns[0].tiers[0].dose ==
        1);
    checkNames([&] { f.loadTier("      - {name: a, dose: 2}\n"); },
               {"Rollout", "tier 'a'", "dose 2", "'TestVax'", "2 doses"});
    checkNames(
        [&] {
          f.load(
              "    vaccine_type: Missing\n    tiers: [{name: a, dose: 0}]\n");
        },
        {"Rollout", "vaccine_type 'Missing'"});
    checkNames([&] { f.load("    tiers: [{name: a, dose: 0}]\n"); },
               {"Rollout", "needs 'vaccine_type'"});
  }

  SUBCASE("a vaccine listed after its campaign still counts") {
    const auto path =
        f.files.write("late.yaml",
                      "campaigns:\n  Rollout: {vaccine_type: TestVax, tiers: "
                      "[{name: a, dose: 0}]}\n"
                      "vaccines:\n  TestVax:\n    doses:\n"
                      "      - {days_to_effective: 1, days_to_waning: 2, "
                      "days_to_finished: 3}\n");
    CHECK_NOTHROW(ConfigLoader::loadVaccination(path.string()));
  }

  SUBCASE("more than 65535 tiers in one campaign are refused") {
    const auto file = [&](int count) {
      std::string tiers = "    vaccine_type: TestVax\n    tiers:\n";
      for (int i = 0; i < count; ++i)
        tiers += "      - {name: t" + std::to_string(i) + ", dose: 0}\n";
      return tiers;
    };
    CHECK(f.load(file(65535)).campaigns[0].tiers.size() == 65535);
    checkNames([&] { f.load(file(65536)); },
               {"Rollout", "65536 tiers", "at most 65535"});
  }

  SUBCASE("unknown keys in a tier are errors") {
    checkNames(
        [&] { f.loadTier("      - {name: a, dose: 0, uptake_rate: 0.5}\n"); },
        {"Rollout", "tier 'a'", "unknown field 'uptake_rate'"});
    checkNames([&] { f.loadTier("      - {name: a, dose: 0, group: []}\n"); },
               {"unknown field 'group'", "selection"});
  }

  SUBCASE("a selection that is not a list is refused") {
    CHECK_THROWS_WITH(
        f.loadTier("      - {name: a, dose: 0, selection: {property: age}}\n"),
        doctest::Contains("selection must be a sequence"));
  }

  SUBCASE("unknown keys in the campaign are errors") {
    checkNames([&] { f.load("    vaccine_type: TestVax\n    tier: []\n"); },
               {"Rollout", "unknown field 'tier'"});
  }
}

TEST_CASE("ConfigLoader - uptake and decision unit validation") {
  CampaignFiles f;
  const auto tier = [&](const std::string& extra) {
    return f.loadTier("      - name: a\n        dose: 0\n        " + extra +
                      "\n");
  };

  SUBCASE("a scalar uptake is a probability") {
    CHECK(tier("uptake: 0").campaigns[0].tiers[0].uptake.rows[0].probability ==
          0.0);
    CHECK(tier("uptake: 1").campaigns[0].tiers[0].uptake.rows[0].probability ==
          1.0);
    for (const char* bad : {"1.2", "-0.1", "most"}) {
      checkNames([&] { tier(std::string("uptake: ") + bad); },
                 {"Rollout", "tier 'a'", "uptake", bad});
    }
  }

  SUBCASE("a list of rows needs a probability in each and no other key") {
    checkNames([&] { tier("uptake: [{selection: []}]"); },
               {"Rollout", "tier 'a'", "uptake row 1", "needs 'probability'"});
    checkNames([&] { tier("uptake: [{probability: 0.5}, {probability: 2}]"); },
               {"Rollout", "tier 'a'", "uptake row 2", "probability", "2"});
    checkNames([&] { tier("uptake: [{probability: 0.5, odds: 3}]"); },
               {"Rollout", "tier 'a'", "uptake row 1", "unknown field 'odds'"});
    checkNames([&] { tier("uptake: []"); },
               {"Rollout", "tier 'a'", "empty list"});
    checkNames([&] { tier("uptake: [0.5]"); },
               {"Rollout", "tier 'a'", "uptake row 1", "mapping"});
  }

  SUBCASE("a file uptake takes only a path") {
    checkNames([&] { tier("uptake: {path: x.csv}"); },
               {"Rollout", "tier 'a'", "unknown field 'path'"});
    checkNames([&] { tier("uptake: {file: ''}"); },
               {"Rollout", "tier 'a'", "needs 'file'"});
    checkNames([&] { tier("uptake: null"); },
               {"Rollout", "tier 'a'", "uptake"});
  }

  SUBCASE(
      "a relative file is read as given when it exists, else beside the yaml") {
    const std::string here = "tests/configs/vaccines.yaml";
    CHECK(tier("uptake: {file: " + here + "}")
              .campaigns[0]
              .tiers[0]
              .uptake.csv_path == here);
    CHECK(tier("uptake: {file: nearby.csv}")
              .campaigns[0]
              .tiers[0]
              .uptake.csv_path == f.files.path("nearby.csv").string());
  }

  SUBCASE("decision_unit is person or a map with its three required keys") {
    CHECK_FALSE(tier("decision_unit: person")
                    .campaigns[0]
                    .tiers[0]
                    .decision_unit.shared());
    checkNames(
        [&] { tier("decision_unit: household"); },
        {"Rollout", "tier 'a'", "decision_unit", "'person'", "household"});
    checkNames([&] { tier("decision_unit: {scope: [], deciders: []}"); },
               {"Rollout", "tier 'a'", "needs 'shared_venue'"});
    checkNames([&] { tier("decision_unit: {shared_venue: r, deciders: []}"); },
               {"Rollout", "tier 'a'", "needs 'scope'"});
    checkNames([&] { tier("decision_unit: {shared_venue: r, scope: []}"); },
               {"Rollout", "tier 'a'", "needs 'deciders'"});
    checkNames(
        [&] {
          tier("decision_unit: {shared_venue: '', scope: [], deciders: []}");
        },
        {"Rollout", "tier 'a'", "shared_venue", "activity name"});
  }

  SUBCASE("combine and if_no_decider take only their words") {
    const std::string unit =
        "decision_unit: {shared_venue: r, scope: [], deciders: [], ";
    checkNames([&] { tier(unit + "combine: median}"); },
               {"Rollout", "tier 'a'", "combine", "min, mean, max", "median"});
    checkNames([&] { tier(unit + "if_no_decider: ignore}"); },
               {"Rollout", "tier 'a'", "if_no_decider", "error, individual",
                "ignore"});
    checkNames([&] { tier(unit + "extra: 1}"); },
               {"Rollout", "tier 'a'", "unknown field 'extra'"});
  }
}

TEST_CASE("VaccinationTier - uptake lookups") {
  CampaignFiles f;
  const WorldState world = campaignTestWorld();  // ages 5 10 17 18 30 70
  const auto probabilities = [&](const VaccinationTier& tier) {
    std::vector<std::optional<double>> out;
    for (const Person& p : world.people)
      out.push_back(tier.uptake.probability(p, world));
    return out;
  };
  using Probabilities = std::vector<std::optional<double>>;
  const std::string head = "      - name: a\n        dose: 0\n        uptake: ";

  SUBCASE("no uptake means everyone accepts") {
    const auto tier = f.resolvedTier("      - {name: a, dose: 0}\n");
    CHECK(probabilities(tier) == Probabilities(6, 1.0));
  }

  SUBCASE("a number applies to everyone") {
    const auto tier =
        f.resolvedTier("      - {name: a, dose: 0, uptake: 0.25}\n");
    CHECK(probabilities(tier) == Probabilities(6, 0.25));
  }

  SUBCASE("the first matching row wins") {
    const auto tier = f.resolvedTier(head + R"(
          - {selection: [{property: age, operator: "<", value: 18}], probability: 0.9}
          - {selection: [{property: age, operator: "<", value: 10}], probability: 0.1}
          - {selection: [{property: age, operator: ">=", value: 18}], probability: 0.5}
)");
    CHECK(probabilities(tier) == Probabilities{0.9, 0.9, 0.9, 0.5, 0.5, 0.5});
  }

  SUBCASE("a row without a selection is the catch-all") {
    const auto tier = f.resolvedTier(head + R"(
          - {selection: [{property: age, operator: "<", value: 18}], probability: 0.9}
          - {probability: 0.3}
)");
    CHECK(probabilities(tier) == Probabilities{0.9, 0.9, 0.9, 0.3, 0.3, 0.3});
  }

  SUBCASE("a person matching no row gets no answer, not a default") {
    const auto tier = f.resolvedTier(head + R"(
          - {selection: [{property: age, operator: "<", value: 18}], probability: 0.9}
)");
    CHECK(probabilities(tier) == Probabilities{0.9, 0.9, 0.9, std::nullopt,
                                               std::nullopt, std::nullopt});
  }

  SUBCASE(
      "rows from a file behave the same: first match, catch-all, no match") {
    f.files.write("uptake.csv",
                  "# a comment line\n"
                  "filter.age,probability\n"
                  "<18,0.9\n"
                  "<10,0.1\n"
                  "18-64,0.5\n");
    const auto tier = f.resolvedTier(head + "{file: uptake.csv}\n");
    CHECK(probabilities(tier) ==
          Probabilities{0.9, 0.9, 0.9, 0.5, 0.5, std::nullopt});

    f.files.write("uptake.csv", "filter.age,probability\n<18,0.9\n,0.3\n");
    const auto with_catch_all = f.resolvedTier(head + "{file: uptake.csv}\n");
    CHECK(probabilities(with_catch_all) ==
          Probabilities{0.9, 0.9, 0.9, 0.3, 0.3, 0.3});
  }

  SUBCASE("inline rows and file rows give the same lookup") {
    f.files.write("uptake.csv",
                  "filter.age,probability\n<10,0.2\n10-17,0.9\n>=65,0.6\n");
    const auto from_file = f.resolvedTier(head + "{file: uptake.csv}\n");
    const auto inline_rows = f.resolvedTier(head + R"(
          - {selection: [{property: age, operator: "<", value: 10}], probability: 0.2}
          - {selection: [{property: age, operator: ">=", value: 10}, {property: age, operator: "<=", value: 17}], probability: 0.9}
          - {selection: [{property: age, operator: ">=", value: 65}], probability: 0.6}
)");
    CHECK(probabilities(from_file) ==
          Probabilities{0.2, 0.9, 0.9, std::nullopt, std::nullopt, 0.6});
    CHECK(probabilities(inline_rows) == probabilities(from_file));
  }

  SUBCASE(
      "asking an unresolved file table for a probability is a bug, not 1.0") {
    f.files.write("uptake.csv", "filter.age,probability\n,0.5\n");
    const auto config = f.loadTier(head + "{file: uptake.csv}\n");
    const auto& tier = config.campaigns[0].tiers[0];
    CHECK_THROWS_AS(tier.uptake.probability(world.people[0], world),
                    std::logic_error);
  }
}

TEST_CASE("VaccinationTier - uptake file errors name the file and row") {
  CampaignFiles f;
  const auto resolve = [&](const std::string& csv) {
    f.files.write("uptake.csv", csv);
    f.resolvedTier("      - {name: a, dose: 0, uptake: {file: uptake.csv}}\n");
  };
  const std::string path = f.files.path("uptake.csv").string();

  SUBCASE("a file without a probability column") {
    checkNames(
        [&] { resolve("filter.age\n<18\n"); },
        {"Rollout", "tier 'a'", path.c_str(), "no 'probability' column"});
  }
  SUBCASE("a column that is neither a filter nor the probability") {
    checkNames(
        [&] { resolve("filter.age,probabilty\n<18,0.5\n"); },
        {"Rollout", "tier 'a'", path.c_str(), "unknown column 'probabilty'"});
    checkNames([&] { resolve("filer.age,probability\n<18,0.5\n"); },
               {path.c_str(), "unknown column 'filer.age'"});
  }
  SUBCASE("a file with no rows") {
    checkNames([&] { resolve("filter.age,probability\n"); },
               {"Rollout", "tier 'a'", path.c_str(), "no data rows"});
  }
  SUBCASE("a cell that is not a number") {
    checkNames([&] { resolve("filter.age,probability\n<18,0.5\n>=18,abc\n"); },
               {"Rollout", "tier 'a'", path.c_str(), "data row 2", "'abc'",
                "not a number"});
    checkNames([&] { resolve("filter.age,probability\n<18,0.5x\n"); },
               {path.c_str(), "data row 1", "'0.5x'"});
    checkNames([&] { resolve("filter.age,probability\n<18,\n"); },
               {path.c_str(), "data row 1", "not a number"});
  }
  SUBCASE("a probability outside 0 to 1") {
    checkNames([&] { resolve("filter.age,probability\n<18,0.5\n>=18,1.5\n"); },
               {"Rollout", "tier 'a'", path.c_str(), "data row 2", "1.5",
                "outside 0 to 1"});
    checkNames([&] { resolve("filter.age,probability\n<18,-0.1\n"); },
               {path.c_str(), "data row 1", "-0.1"});
    checkNames([&] { resolve("filter.age,probability\n<18,nan\n"); },
               {path.c_str(), "data row 1", "nan"});
  }
  SUBCASE("a file that does not exist") {
    auto config =
        f.loadTier("      - {name: a, dose: 0, uptake: {file: gone.csv}}\n");
    const WorldState world = campaignTestWorld();
    checkNames([&] { config.resolve(world); },
               {"Rollout", "tier 'a'", "gone.csv"});
  }
  SUBCASE("a filter on an infection fact describes no person") {
    checkNames(
        [&] {
          resolve("filter.transmission_mode,probability\nrespiratory,0.5\n");
        },
        {"Rollout", "tier 'a'", "row 1", "transmission_mode"});
  }
}

TEST_CASE("VaccinationConfig - resolve names campaign and tier") {
  CampaignFiles f;
  WorldState world = campaignTestWorld();

  SUBCASE("a shared venue is resolved to its activity index") {
    auto config = f.loadTier(R"(      - name: a
        dose: 0
        decision_unit: {shared_venue: residence, scope: [], deciders: []}
)");
    config.resolve(world);
    const auto& unit = config.campaigns[0].tiers[0].decision_unit;
    CHECK(unit.activity_index == world.getActivityIndex("residence"));
    CHECK(unit.activity_index == 1);
  }

  SUBCASE("an unknown activity throws") {
    auto config = f.loadTier(R"(      - name: a
        dose: 0
        decision_unit: {shared_venue: nowhere, scope: [], deciders: []}
)");
    checkNames([&] { config.resolve(world); },
               {"Rollout", "tier 'a'", "unknown activity 'nowhere'"});
  }

  SUBCASE("an unknown property in each criterion list throws") {
    const std::string bad =
        "{property: properties.no_such, operator: '==', value: 1}";
    const std::string unit = "decision_unit: {shared_venue: residence, ";
    for (const std::string& line :
         {"selection: [" + bad + "]",
          "uptake: [{selection: [" + bad + "], probability: 0.5}]",
          unit + "scope: [" + bad + "], deciders: []}",
          unit + "scope: [], deciders: [" + bad + "]}"}) {
      auto config = f.loadTier("      - {name: a, dose: 0, " + line + "}\n");
      checkNames([&] { config.resolve(world); },
                 {"Rollout", "tier 'a'", "no_such"});
    }
  }

  SUBCASE("a disabled vaccination config is not resolved") {
    auto config = f.loadTier(
        "      - {name: a, dose: 0, selection: [{property: properties.no_such, "
        "operator: '==', value: 1}]}\n");
    config.enabled = false;
    CHECK_NOTHROW(config.resolve(world));
  }
}
