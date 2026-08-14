#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <filesystem>
#include <iostream>
#include <sstream>

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

    CHECK(config.campaigns[0].vaccine_type == "Pfizer");
    CHECK(config.campaigns[0].daily_coverage == doctest::Approx(0.3));
    REQUIRE(config.campaigns[0].selection_criteria.size() == 1);
    CHECK(config.campaigns[0].selection_criteria[0].property_path ==
          "activities.residence.venue_type");

    CHECK(config.campaigns[1].vaccine_type == "AstraZeneca");

    // Booster campaign
    CHECK(config.campaigns[2].dose_sequence.size() == 1);
    CHECK(config.campaigns[2].dose_sequence[0] == 2);
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
    std::string yaml_path = "tmp_sim_cal_events.yaml";
    {
      std::ofstream f(yaml_path);
      f << "time:\n"
           "  start_date: \"2020-01-01\"\n"
           "  end_date: \"2020-01-10\"\n"
           "config_paths:\n"
           "  calendar_events_file: \"data/calendar_events.csv\"\n"
           "  calendar_event_catchment_rules_file: "
           "\"data/catchment_rules.csv\"\n";
    }
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path);
    CHECK(cfg.calendar_events_file == "data/calendar_events.csv");
    CHECK(cfg.calendar_event_catchment_rules_file ==
          "data/catchment_rules.csv");
    std::filesystem::remove(yaml_path);
  }

  SUBCASE("absent keys leave fields empty, no error") {
    std::string yaml_path = "tmp_sim_no_cal_events.yaml";
    {
      std::ofstream f(yaml_path);
      f << "time:\n"
           "  start_date: \"2020-01-01\"\n"
           "  end_date: \"2020-01-10\"\n"
           "config_paths:\n"
           "  disease_file: \"disease.yaml\"\n";
    }
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path);
    CHECK(cfg.calendar_events_file == "");
    CHECK(cfg.calendar_event_catchment_rules_file == "");
    std::filesystem::remove(yaml_path);
  }
}

TEST_CASE("SimulationConfig - save_coordinated_encounters parsed from output") {
  SUBCASE("absent defaults to false") {
    std::string yaml_path = "tmp_sim_output_no_sce.yaml";
    {
      std::ofstream f(yaml_path);
      f << "time:\n"
           "  start_date: \"2020-01-01\"\n"
           "  end_date: \"2020-01-10\"\n"
           "output:\n"
           "  stats_interval_days: 1\n";
    }
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path);
    CHECK(cfg.save_coordinated_encounters == false);
    std::filesystem::remove(yaml_path);
  }

  SUBCASE("explicit true is honoured") {
    std::string yaml_path = "tmp_sim_output_sce_true.yaml";
    {
      std::ofstream f(yaml_path);
      f << "time:\n"
           "  start_date: \"2020-01-01\"\n"
           "  end_date: \"2020-01-10\"\n"
           "output:\n"
           "  stats_interval_days: 1\n"
           "  save_coordinated_encounters: true\n";
    }
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path);
    CHECK(cfg.save_coordinated_encounters == true);
    std::filesystem::remove(yaml_path);
  }

  SUBCASE("explicit false is honoured") {
    std::string yaml_path = "tmp_sim_output_sce_false.yaml";
    {
      std::ofstream f(yaml_path);
      f << "time:\n"
           "  start_date: \"2020-01-01\"\n"
           "  end_date: \"2020-01-10\"\n"
           "output:\n"
           "  stats_interval_days: 1\n"
           "  save_coordinated_encounters: false\n";
    }
    SimulationConfig cfg = ConfigLoader::loadSimulation(yaml_path);
    CHECK(cfg.save_coordinated_encounters == false);
    std::filesystem::remove(yaml_path);
  }
}

TEST_CASE("ConfigLoader - Optional Files Handling") {
  SUBCASE("Missing vaccines file is optional and disables the module") {
    // Vaccines is the only truly-optional sub-config: a missing file means
    // "don't run this module", which is fine — the rest of the simulator
    // can proceed.
    VaccinationConfig vc =
        ConfigLoader::loadVaccination("non_existent_vaccines.yaml");
    CHECK(vc.enabled == false);
  }

  SUBCASE("Missing coordinated_encounters file is a misconfiguration") {
    // simulation.yaml's coordinated_encounters_file points at this YAML;
    // absence is a bug, not a silent opt-out.
    CHECK_THROWS_AS(
        ConfigLoader::loadCoordinatedEncounters("non_existent_ce.yaml"),
        std::runtime_error);
  }

  SUBCASE("Corrupted coordinated_encounters file must throw") {
    // Required field missing should fail loudly; previously the loader
    // wrapped this in a try/catch and disabled the feature with a cerr
    // warning, which masked real bugs in production runs.
    std::string bad_file = "tmp_bad_ce.yaml";
    {
      std::ofstream f(bad_file);
      f << "coordinated_encounters:\n  encounters:\n    - name: bad\n      "
           "network: missing_fields\n";
    }

    CHECK_THROWS_AS(ConfigLoader::loadCoordinatedEncounters(bad_file),
                    std::runtime_error);

    std::filesystem::remove(bad_file);
  }
}
