cmake_minimum_required(VERSION 3.14)
project(disease_simulation VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ABI Compatibility Note:
# If you get "undefined reference" errors for HDF5 / yaml-cpp functions, check
# whether those libraries use the old C++11 ABI (pre-GCC 5.1). Modern builds
# use the new ABI by default; configure the compatibility definition explicitly
# when linking against libraries built with the old ABI.

# Find HDF5
if(USE_MPI)
  set(HDF5_PREFER_PARALLEL TRUE)
endif()

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(yaml-cpp CONFIG REQUIRED)

# Optional: Find MPI for parallel execution
option(USE_MPI "Enable MPI parallel execution" ON)
if(USE_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    message(STATUS "MPI found: Parallel execution enabled")
  else()
    message(WARNING "MPI not found: Building serial version only")
  endif()
endif()

# Optional: Find METIS for graph partitioning (better than RCB)
option(USE_METIS "Enable METIS graph partitioning" ON)
if(MPI_FOUND AND USE_METIS)
  # Try to find METIS using CMake variable, environment variable, or brew prefix on macOS
  if(METIS_DIR)
    set(METIS_HINTS HINTS ${METIS_DIR}/lib ${METIS_DIR}/lib64 ${METIS_DIR}/include)
  elseif(DEFINED ENV{METIS_DIR})
    set(METIS_HINTS HINTS $ENV{METIS_DIR}/lib $ENV{METIS_DIR}/lib64 $ENV{METIS_DIR}/include)
  elseif(APPLE)
    execute_process(COMMAND brew --prefix metis OUTPUT_VARIABLE METIS_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(METIS_PREFIX)
      set(METIS_HINTS HINTS ${METIS_PREFIX}/lib ${METIS_PREFIX}/include)
    endif()
  endif()

  find_library(METIS_LIBRARY metis ${METIS_HINTS})
  find_path(METIS_INCLUDE_DIR metis.h ${METIS_HINTS})

  if(METIS_LIBRARY AND METIS_INCLUDE_DIR)
    message(STATUS "METIS found: ${METIS_LIBRARY}")
    message(STATUS "METIS include: ${METIS_INCLUDE_DIR}")
    set(METIS_FOUND TRUE)
  else()
    message(WARNING "METIS not found: Will use RCB partitioning as fallback")
    set(METIS_FOUND FALSE)
  endif()
endif()

# ==============================================================================
# Profiling options - Use gperftools (Google CPU Profiler)
# ==============================================================================
# AUTO (the default) profiles when gperftools is found and warns when it is
# not. ON is a request: without gperftools the configure step fails, because a
# build that silently cannot profile is not what was asked for. OFF skips it.
set(ENABLE_PROFILING AUTO CACHE STRING
    "CPU profiling with gperftools: AUTO, ON or OFF")
set_property(CACHE ENABLE_PROFILING PROPERTY STRINGS AUTO ON OFF)

# Find gperftools
if(ENABLE_PROFILING)
  if(GPERFTOOLS_DIR)
    set(GPERFTOOLS_HINTS HINTS ${GPERFTOOLS_DIR}/lib ${GPERFTOOLS_DIR}/lib64)
  elseif(DEFINED ENV{GPERFTOOLS_DIR})
    set(GPERFTOOLS_HINTS HINTS $ENV{GPERFTOOLS_DIR}/lib $ENV{GPERFTOOLS_DIR}/lib64)
  elseif(APPLE)
    set(GPERFTOOLS_HINTS HINTS /opt/homebrew/lib)
  endif()

  find_library(PROFILER_LIBRARY profiler ${GPERFTOOLS_HINTS})

  # The header search needs its own hint. Only the library was hinted before,
  # which works on macOS because Homebrew's include prefix is searched by
  # default, but fails on any HPC install in a non-standard prefix: the build
  # links fine and then main.cpp cannot find gperftools/profiler.h.
  if(GPERFTOOLS_DIR)
    set(GPERFTOOLS_INC_HINTS HINTS ${GPERFTOOLS_DIR}/include)
  elseif(DEFINED ENV{GPERFTOOLS_DIR})
    set(GPERFTOOLS_INC_HINTS HINTS $ENV{GPERFTOOLS_DIR}/include)
  elseif(APPLE)
    set(GPERFTOOLS_INC_HINTS HINTS /opt/homebrew/include)
  endif()
  find_path(PROFILER_INCLUDE_DIR gperftools/profiler.h ${GPERFTOOLS_INC_HINTS})

  if(PROFILER_LIBRARY AND PROFILER_INCLUDE_DIR)
    message(STATUS "gperftools found: ${PROFILER_LIBRARY}")
    message(STATUS "gperftools headers: ${PROFILER_INCLUDE_DIR}")
    set(GPERFTOOLS_FOUND TRUE)
  else()
    if(PROFILER_LIBRARY)
      string(CONCAT GPERFTOOLS_MISSING
          "the gperftools library was found at ${PROFILER_LIBRARY} but "
          "gperftools/profiler.h was not; set GPERFTOOLS_DIR to a prefix "
          "containing both lib and include")
    else()
      string(CONCAT GPERFTOOLS_MISSING
          "gperftools was not found (brew install gperftools, or set "
          "GPERFTOOLS_DIR)")
    endif()
    if(ENABLE_PROFILING STREQUAL "AUTO")
      message(WARNING "${GPERFTOOLS_MISSING}. Building without profiling.")
    else()
      message(FATAL_ERROR
        "ENABLE_PROFILING=${ENABLE_PROFILING} asks for profiling, but "
        "${GPERFTOOLS_MISSING}. Fix that, or configure with "
        "-DENABLE_PROFILING=AUTO or -DENABLE_PROFILING=OFF.")
    endif()
    set(GPERFTOOLS_FOUND FALSE)
  endif()
endif()

# ==============================================================================
# Core Simulation Library
# ==============================================================================
set(JUNE_CORE_SOURCES
    src/epidemiology/infection_seed.cpp
    src/epidemiology/seeding/seed_selector.cpp
    src/epidemiology/seeding/seed_shortfall.cpp
    src/epidemiology/seeding/seed_cluster_planner.cpp
    src/epidemiology/seeding/seed_identity.cpp
    src/epidemiology/seeding/seed_window.cpp
    src/epidemiology/fomite/fomite_sub_bins.cpp
    src/epidemiology/emission/emission.cpp
    src/epidemiology/applicable_policies/applicable_policy_sets.cpp
    src/epidemiology/transmission/infection_context.cpp
    src/epidemiology/transmission/infector_symptom_lookup.cpp
    src/epidemiology/transmission/partial_presence_sources.cpp
    src/epidemiology/trajectory/stage_curve_integral.cpp
    src/epidemiology/disease.cpp
    src/loaders/disease_loader.cpp
    src/interaction/interaction_manager.cpp
    src/interaction/interaction_manager_venue.cpp
    src/interaction/interaction_manager_venue_bins.cpp
    src/interaction/interaction_manager_partial_presence.cpp
    src/activity/activity_manager.cpp
    src/activity/activity_manager_selection.cpp
    src/activity/on_the_fly_venue_allocator.cpp
    src/core/config.cpp
    src/core/config_checks.cpp
    src/core/world_state.cpp
    src/activity/coordinated_encounter_manager.cpp
    src/activity/coordinated_encounter_manager_stats.cpp
    src/activity/runtime_group_allocator.cpp
    src/activity/presence_window.cpp
    src/utils/event_logging/event_logger.cpp
    src/utils/filtering.cpp
    src/utils/filtered_csv.cpp
    src/utils/event_logging/event_writer.cpp
    src/utils/event_logging/event_writer_lookups.cpp
    src/utils/event_logging/event_merger.cpp
    src/utils/event_logging/event_merger_lookups.cpp
    src/epidemiology/epidemiology.cpp
    src/epidemiology/observed_metrics.cpp
    src/epidemiology/policy.cpp
    src/epidemiology/policy_trigger.cpp
    src/epidemiology/calendar_event.cpp
    src/loaders/calendar_event_loader.cpp
    src/loaders/catchment_rule_loader.cpp
    src/simulation/simulator.cpp
    src/simulation/encounters/eligibility.cpp
    src/simulation/simulator_checkpoint.cpp
    src/simulation/simulator_encounters.cpp
    src/simulation/simulator_reporting.cpp
    src/simulation/simulator_restore.cpp
    src/simulation/simulator_timeslot.cpp
    src/simulation/compartmental_model_manager.cpp
    src/loaders/hdf5_loader.cpp
    src/loaders/domain_loader.cpp
    src/loaders/domain_loader_internals.cpp
    src/loaders/config_loader.cpp
    src/loaders/config_loader_encounters.cpp
    src/loaders/config_loader_vaccination.cpp
    src/loaders/policy_loader.cpp
    src/loaders/policy_loader_transmission.cpp)

if(MPI_FOUND)
  list(APPEND JUNE_CORE_SOURCES
       src/parallel/domain.cpp
       src/parallel/domain_manager.cpp
       src/parallel/geography_partitioner.cpp
       src/parallel/domain_communicator.cpp
       src/parallel/domain_communicator_encounters.cpp
       src/parallel/visitor_wire.cpp
       src/parallel/visitor_payload.cpp)
endif()

# Create the library
add_library(june_core STATIC ${JUNE_CORE_SOURCES})

# Fused multiply-add changes the last bit of the random draws, and the default
# differs between compilers. Keep one rounding per operation everywhere.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(june_core PUBLIC -ffp-contract=off)
endif()

# GCC 15 false-positive: variant copy-assign with vector<int> triggers a
# spurious -Wfree-nonheap-object due to broken alias analysis in the compiler.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(june_core PRIVATE -Wno-free-nonheap-object)
endif()

# Define include directories for the library
# PUBLIC means anyone linking to june_core gets these headers too
target_include_directories(june_core PUBLIC ${CMAKE_SOURCE_DIR}/include)

if(MPI_FOUND)
  target_compile_definitions(june_core PUBLIC USE_MPI)
endif()

if(MPI_FOUND AND METIS_FOUND)
  target_include_directories(june_core PRIVATE ${METIS_INCLUDE_DIR})
  target_compile_definitions(june_core PRIVATE USE_METIS)
endif()

# Link mandatory libraries
target_link_libraries(june_core PUBLIC dl)

target_link_libraries(june_core PUBLIC yaml-cpp::yaml-cpp)
target_link_libraries(june_core PUBLIC hdf5::hdf5_cpp hdf5::hdf5)

if(MPI_FOUND)
  target_link_libraries(june_core PUBLIC MPI::MPI_CXX)
endif()

if(MPI_FOUND AND METIS_FOUND)
  target_link_libraries(june_core PRIVATE ${METIS_LIBRARY})
endif()

# ==============================================================================
# Main Executable
# ==============================================================================
add_executable(disease_sim src/main.cpp)
target_link_libraries(disease_sim PRIVATE june_core)
if(GPERFTOOLS_FOUND)
  target_link_libraries(disease_sim PRIVATE ${PROFILER_LIBRARY})
  target_compile_definitions(disease_sim PRIVATE USE_GPERFTOOLS)
  target_include_directories(disease_sim PRIVATE ${PROFILER_INCLUDE_DIR})
endif()

add_executable(merge_events src/merge_events.cpp)
target_link_libraries(merge_events PRIVATE june_core)

if(GPERFTOOLS_FOUND)
  message(STATUS "Profiling enabled: gperftools CPU profiler linked via disease_sim")
endif()

# ==============================================================================
# Unit Tests
# ==============================================================================
option(BUILD_TESTING "Build unit tests" ON)

if(BUILD_TESTING)
  enable_testing()

  # Add a standalone test target. Pass an MPI launcher as the optional third
  # argument, and an optional rank count as the fourth, for MPI tests.
  macro(add_june_test name sources)
    add_executable(${name} ${sources})
    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
    target_link_libraries(${name} PRIVATE june_core)
    if("${ARGC}" GREATER 2)
      set(_june_mpi_num_ranks 2)
      if("${ARGC}" GREATER 3)
        set(_june_mpi_num_ranks ${ARGV3})
      endif()
      add_test(NAME ${name}
               COMMAND ${ARGV2} -np ${_june_mpi_num_ranks}
                       --oversubscribe ${CMAKE_CURRENT_BINARY_DIR}/${name}
               WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    else()
      add_test(NAME ${name} COMMAND ${name} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
  endmacro()

  set(JUNE_STANDARD_TESTS
      test_domain_communicator_detail
      test_run_dir
      test_world_state
      test_selection_criterion_geography
      test_domain_loader_internals
      test_infection
      test_disease_loader
      test_transmission_lookups
      test_coordinated_encounters

      test_fomite_transmission
      test_vaccines
      test_vaccination_manager
      test_config_loading
      test_config_checks
      test_activity_manager

      test_transmission_modifiers
      test_applicable_policy_sets
      test_on_the_fly_venue_allocator
      test_schedule_csv
      test_day_types
      test_schedule_hopper
      test_schedule_hop_state
      test_calendar_event_trigger
      test_calendar_event_loader
      test_calendar_event_venue
      test_calendar_event_checkpoint
      test_catchment_rule_loader

      test_seed_selector
      test_seed_shortfall
      test_seed_cluster_planner
      test_cross_domain_venues
      test_membership_metadata_matching
      test_integration
      test_epidemiology
      test_observed_metrics
      test_policy_trigger
      test_policy_response
      test_simulator

      test_deterministic_rng
      test_event_record_layout
      test_encounter_group_id_rank_independent

      test_transmission_determinism

      test_presence_window
      test_partial_presence_foi

      test_follow
      test_policy
      test_infection_seed_loader
      test_infection_seeder
      test_transmission_modes
      test_foi_calculation
      test_visitor_wire
      test_visitor_payload
      test_filtering
      test_outcome_rate_context
      test_infector_symptom_lookup
      test_partial_presence_sources
      test_fomite_sub_bins
      test_emission
      test_stage_curve_integral
      test_seed_window)
  add_executable(june_tests tests/june_tests_main.cpp)
  set(JUNE_STANDARD_TEST_SOURCES)
  foreach(test_name IN LISTS JUNE_STANDARD_TESTS)
    list(APPEND JUNE_STANDARD_TEST_SOURCES "tests/${test_name}.cpp")
  endforeach()
  target_sources(june_tests PRIVATE ${JUNE_STANDARD_TEST_SOURCES})
  target_include_directories(june_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(june_tests PRIVATE june_core)
  add_test(NAME june_tests COMMAND june_tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  target_include_directories(june_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)
  add_june_test(test_sibling_mixing "tests/test_sibling_mixing.cpp")
  add_june_test(test_compartmental_disease_mode "tests/test_compartmental_disease_mode.cpp")
  add_june_test(test_compartmental_coupling_matrix "tests/test_compartmental_coupling_matrix.cpp")
  add_june_test(test_compartmental_transmission "tests/test_compartmental_transmission.cpp")
  add_june_test(test_compartmental_writeback "tests/test_compartmental_writeback.cpp")
  add_june_test(test_compartmental_integration "tests/test_compartmental_integration.cpp")
  # Mock compartmental model plugin .so used by compartmental model manager tests.
  add_library(mock_compartmental_model SHARED tests/mock_compartmental_model.cpp)
  target_include_directories(mock_compartmental_model PRIVATE ${CMAKE_SOURCE_DIR}/include)

  add_june_test(test_compartmental_model_manager "tests/test_compartmental_model_manager.cpp")
  # Pass the built .so path so the test can load it without hardcoding.
  target_compile_definitions(test_compartmental_model_manager PRIVATE
      MOCK_COMPARTMENTAL_MODEL_SO_PATH="$<TARGET_FILE:mock_compartmental_model>")
  add_dependencies(test_compartmental_model_manager mock_compartmental_model)
  # Slice 3: test writes synthetic HDF5 graphs using the C API directly.
  target_link_libraries(test_compartmental_model_manager PRIVATE hdf5::hdf5)

  if(MPI_FOUND)
    if(NOT MPIEXEC_EXECUTABLE)
      find_program(MPIEXEC_EXECUTABLE NAMES mpirun mpiexec)
    endif()
    if(NOT MPIEXEC_EXECUTABLE)
      message(FATAL_ERROR "MPI execution command (mpirun/mpiexec) not found, but MPI tests are enabled.")
    endif()

    add_june_test(test_mpi_communication "tests/test_mpi_communication.cpp" "${MPIEXEC_EXECUTABLE}")
    add_june_test(test_mpi_advanced "tests/test_mpi_advanced.cpp" "${MPIEXEC_EXECUTABLE}")
    add_june_test(test_mpi_registry_sync "tests/test_mpi_registry_sync.cpp" "${MPIEXEC_EXECUTABLE}")
    add_june_test(test_mpi_observed_metrics "tests/test_mpi_observed_metrics.cpp" "${MPIEXEC_EXECUTABLE}")
    add_june_test(test_mpi_policy_trigger "tests/test_mpi_policy_trigger.cpp" "${MPIEXEC_EXECUTABLE}")
    add_june_test(test_mpi_transmission_modes "tests/test_mpi_transmission_modes.cpp" "${MPIEXEC_EXECUTABLE}")
    add_june_test(test_mpi_visitor_exchange_paths "tests/test_mpi_visitor_exchange_paths.cpp" "${MPIEXEC_EXECUTABLE}" 4)

      # End-to-end disease_sim reproducibility: runs disease_sim at
      # np=1/2/3 against the real config + world and h5diff's the
      # /events/relationships, /events/coordinated_encounters (minus
      # partition-dependent group_id), and /events/infections datasets.
      # Covers R2 (OOE roster), E2 (coordinated encounters), and
      # infection counts bit-identical across np. Slow (~30s at 5 days);
      # opt out via ctest -E mpi_full_repro.
      add_test(NAME test_mpi_full_reproducibility
               COMMAND ${CMAKE_SOURCE_DIR}/tests/test_mpi_full_reproducibility.sh
               WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
      set_tests_properties(test_mpi_full_reproducibility PROPERTIES
                           ENVIRONMENT "BUILD_DIR=${CMAKE_CURRENT_BINARY_DIR};DAYS=5"
                           SKIP_RETURN_CODE 77 TIMEOUT 600)

      # The same end-to-end diff, seeded structurally above the partition
      # level (an exact seed in an LGU and in an XLGU, MGU being where
      # parallel.yaml partitions). A structured seed states an absolute
      # count, so this is what catches a unit being seeded once per rank
      # covering it.
      add_test(NAME test_mpi_structured_seed_reproducibility
               COMMAND ${CMAKE_SOURCE_DIR}/tests/test_mpi_full_reproducibility.sh
               WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
      set_tests_properties(test_mpi_structured_seed_reproducibility PROPERTIES
                           ENVIRONMENT "BUILD_DIR=${CMAKE_CURRENT_BINARY_DIR};DAYS=3;SEEDS=configs/config_2021/infection_seeds_structured.yaml"
                           SKIP_RETURN_CODE 77 TIMEOUT 600)

      # Checkpoint/restart determinism gate: a checkpoint leg + resume leg
      # must be set-identical to an uninterrupted baseline across a
      # (write_ranks -> resume_ranks) matrix, incl. asymmetric resume. Skips
      # cleanly (exit 77) if the world file or venv python is unavailable
      # (e.g. CI without data). Opt out via: ctest -E checkpoint_determinism.
      find_program(JUNE_VENV_PYTHON
                   NAMES python3
                   PATHS ${CMAKE_SOURCE_DIR}/.venv/bin NO_DEFAULT_PATH)
      if(JUNE_VENV_PYTHON)
        add_test(NAME test_checkpoint_determinism
                 COMMAND ${JUNE_VENV_PYTHON}
                         ${CMAKE_SOURCE_DIR}/tests/checkpoint_determinism_check.py
                         --bin ${CMAKE_CURRENT_BINARY_DIR}/disease_sim
                         --config configs/config_2021/simulation.yaml
                         --world worlds/world_2021.h5
                         --days 4 --ckpt-day 2
                         --matrix 1:1,2:2,2:1
                 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
        set_tests_properties(test_checkpoint_determinism PROPERTIES
                             SKIP_RETURN_CODE 77 TIMEOUT 900)
      endif()
  endif()
endif()
