#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <fstream>
#include <memory>
#include <vector>

#include "doctest.h"
#include "epidemiology/disease.h"
#include "epidemiology/infectiousness_curves.h"
#include "loaders/disease_loader.h"

using namespace june;

// ---------------------------------------------------------------------------
// Minimal YAML with compartmental_uptake and compartmental_deposition modes.
// Written to /tmp so DiseaseLoader::loadFromYAML() can open it.
// ---------------------------------------------------------------------------
static std::string writeTempDiseaseYaml() {
  const char* path = "/tmp/test_compartmental_disease.yaml";
  std::ofstream f(path);
  f << R"(
disease:
  name: test_plague
  symptom_tags:
    - name: healthy
      value: -1
    - name: bacteraemia
      value: 1
    - name: buboes_and_fever
      value: 2
  transmission:
    mode: Stage-Driven
    modes:
      - name: direct_contact
        mode_transmissibility_multiplier: 1.0
      - name: compartmental_uptake
        type: compartmental_uptake
        mode_transmissibility_multiplier: 1.5
      - name: compartmental_deposition
        type: compartmental_deposition
        deposition_stages:
          bacteraemia:
            type: constant
            value: 0.03
          buboes_and_fever:
            type: linear_ramp
            start_value: 0.0
            end_value: 0.03
            ramp_duration: 4.0
    stage_curves:
      direct_contact:
        bacteraemia:
          type: constant
          value: 0.5
        buboes_and_fever:
          type: constant
          value: 1.0
  trajectories:
    - description: general
      selection_key: general_population
      probability: 1.0
      severity: 0.5
      stages:
        - symptom_tag: bacteraemia
          completion_time:
            type: constant
            value: 5.0
        - symptom_tag: buboes_and_fever
          completion_time:
            type: constant
            value: 5.0
        - symptom_tag: healthy
          completion_time:
            type: constant
            value: 100.0
)";
  return path;
}

// ---------------------------------------------------------------------------
// Direct-construction tests (no YAML, no disease_loader.cpp)
// ---------------------------------------------------------------------------

TEST_CASE("CompartmentalUptakeConfig default-constructs") {
  CompartmentalUptakeConfig cfg;
  CHECK(cfg.mode_index == -1);
}

TEST_CASE("CompartmentalDepositionConfig default-constructs empty") {
  CompartmentalDepositionConfig cfg;
  CHECK(cfg.mode_index == -1);
  CHECK(cfg.deposition_by_symptom.empty());
}

TEST_CASE("TransmissionParams carries compartmental mode fields") {
  TransmissionParams tp;
  tp.mode = InfectiousnessMode::STAGE_DRIVEN;

  TransmissionMode direct;
  direct.name = "direct";
  direct.mode_transmissibility_multiplier = 1.0;
  tp.modes.push_back(direct);

  TransmissionMode uptake;
  uptake.name = "uptake";
  uptake.type = TransmissionModeType::CompartmentalUptake;
  uptake.mode_transmissibility_multiplier = 1.5;
  CompartmentalUptakeConfig ucfg;
  ucfg.mode_index = 1;
  uptake.config = ucfg;
  tp.modes.push_back(uptake);

  TransmissionMode deposition;
  deposition.name = "deposition";
  deposition.type = TransmissionModeType::CompartmentalDeposition;
  CompartmentalDepositionConfig dcfg;
  dcfg.mode_index = 2;
  dcfg.deposition_by_symptom = {
      nullptr,
      std::make_shared<ConstantCurve>(0.03),
      std::make_shared<ConstantCurve>(0.06),
  };
  deposition.config = std::move(dcfg);
  tp.modes.push_back(std::move(deposition));

  CHECK(tp.modes[1].type == TransmissionModeType::CompartmentalUptake);
  CHECK(tp.modes[0].type == TransmissionModeType::Standard);
  CHECK(tp.modes[2].type == TransmissionModeType::CompartmentalDeposition);
  CHECK(std::get<CompartmentalUptakeConfig>(tp.modes[1].config).mode_index ==
        1);
  const auto& dep = std::get<CompartmentalDepositionConfig>(tp.modes[2].config);
  CHECK(dep.deposition_by_symptom.size() == 3);
  CHECK(dep.deposition_by_symptom[1] != nullptr);
}

// ---------------------------------------------------------------------------
// YAML-parsing tests — exercise disease_loader.cpp
// ---------------------------------------------------------------------------

TEST_CASE("DiseaseLoader: compartmental_uptake mode parsed") {
  std::string yaml_path = writeTempDiseaseYaml();
  Disease disease = DiseaseLoader::loadFromYAML(yaml_path);
  const auto& tp = disease.getTransmissionParams();

  // mode 0 = direct_contact, 1 = compartmental_uptake, 2 =
  // compartmental_deposition
  REQUIRE(tp.modes.size() == 3);
  CHECK(tp.modes[1].name == "compartmental_uptake");
  CHECK(tp.modes[1].type == TransmissionModeType::CompartmentalUptake);
  CHECK(tp.modes[0].type == TransmissionModeType::Standard);
  CHECK(tp.modes[2].type != TransmissionModeType::CompartmentalUptake);
  CHECK(std::get<CompartmentalUptakeConfig>(tp.modes[1].config).mode_index ==
        1);
  CHECK(tp.modes[1].mode_transmissibility_multiplier == doctest::Approx(1.5));
}

TEST_CASE("DiseaseLoader: compartmental_deposition mode parsed") {
  std::string yaml_path = writeTempDiseaseYaml();
  Disease disease = DiseaseLoader::loadFromYAML(yaml_path);
  const auto& tp = disease.getTransmissionParams();

  REQUIRE(tp.modes.size() == 3);
  CHECK(tp.modes[2].name == "compartmental_deposition");
  CHECK(tp.modes[2].type == TransmissionModeType::CompartmentalDeposition);
  CHECK(tp.modes[0].type != TransmissionModeType::CompartmentalDeposition);
  CHECK(tp.modes[1].type != TransmissionModeType::CompartmentalDeposition);

  const auto& dcfg =
      std::get<CompartmentalDepositionConfig>(tp.modes[2].config);
  CHECK(dcfg.mode_index == 2);
  // 3 symptom tags: healthy(0), bacteraemia(1), buboes_and_fever(2)
  REQUIRE(dcfg.deposition_by_symptom.size() == 3);
  CHECK(dcfg.deposition_by_symptom[0] == nullptr);
  CHECK(dcfg.deposition_by_symptom[1] != nullptr);
  CHECK(dcfg.deposition_by_symptom[2] != nullptr);
}

TEST_CASE("DiseaseLoader: non-compartmental modes unaffected") {
  std::string yaml_path = writeTempDiseaseYaml();
  Disease disease = DiseaseLoader::loadFromYAML(yaml_path);
  const auto& tp = disease.getTransmissionParams();

  CHECK(tp.modes[0].type == TransmissionModeType::Standard);
  CHECK(tp.modes[1].type != TransmissionModeType::Fomite);
  CHECK(tp.modes[2].type != TransmissionModeType::Fomite);
}
