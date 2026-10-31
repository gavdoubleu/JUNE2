#include "utils/event_logging/event_merger.h"

#include <iostream>
#include <stdexcept>

#include "utils/event_logging/event_schemas.h"

namespace june {

void mergeEventFiles(const std::vector<std::string>& input_files,
                     const std::string& output_file) {
  std::cout << "\n=== Merging Event Files ===" << std::endl;
  std::cout << "Input files: " << input_files.size() << std::endl;
  std::cout << "Output file: " << output_file << std::endl;

  try {
    H5::H5File out_file(output_file, H5F_ACC_TRUNC);
    out_file.createGroup("/events");
    out_file.createGroup("/lookups");

    event_merger_detail::mergeDatasetTemplate<InfectionEvent>(
        out_file, "/events/infections", input_files, event_schema::infection());
    event_merger_detail::mergeDatasetTemplate<SymptomChangeEvent>(
        out_file, "/events/symptom_changes", input_files,
        event_schema::symptomChange());
    event_merger_detail::mergeDatasetTemplate<DeathEvent>(
        out_file, "/events/deaths", input_files, event_schema::death());
    event_merger_detail::mergeDatasetTemplate<detail::HospitalAdmissionRecord>(
        out_file, "/events/hospital_admissions", input_files,
        event_schema::hospitalAdmission());
    event_merger_detail::mergeDatasetTemplate<ICUAdmissionEvent>(
        out_file, "/events/icu_admissions", input_files,
        event_schema::icuAdmission());
    event_merger_detail::mergeDatasetTemplate<detail::HospitalDischargeRecord>(
        out_file, "/events/hospital_discharges", input_files,
        event_schema::hospitalDischarge());
    event_merger_detail::mergeDatasetTemplate<VaccinationEvent>(
        out_file, "/events/vaccinations", input_files,
        event_schema::vaccination());
    event_merger_detail::mergeDatasetTemplate<RelationshipEvent>(
        out_file, "/events/relationships", input_files,
        event_schema::relationship());
    std::cout << "  Merged /events/relationships" << std::endl;
    event_merger_detail::mergeDatasetTemplate<CoordinatedEncounterEvent>(
        out_file, "/events/coordinated_encounters", input_files,
        event_schema::coordinatedEncounter());
    std::cout << "  Merged /events/coordinated_encounters" << std::endl;
    event_merger_detail::mergeDatasetTemplate<FollowEvent>(
        out_file, "/events/follows", input_files, event_schema::follow());
    std::cout << "  Merged /events/follows" << std::endl;
    event_merger_detail::mergeDatasetTemplate<CampaignDayEvent>(
        out_file, "/events/campaign_days", input_files,
        event_schema::campaignDay());
    std::cout << "  Merged /events/campaign_days" << std::endl;

    event_merger_detail::mergePeopleLookup(out_file, input_files);
    event_merger_detail::mergeVenueLookup(out_file, input_files);
    event_merger_detail::mergePersonActivityLookup(out_file, input_files);
    event_merger_detail::mergePopulationSummary(out_file, input_files);
    event_merger_detail::mergeProfileAssignments(out_file, input_files);
    event_merger_detail::mergePopulationNetworks(out_file, input_files);

    // Copy metadata (registries) from first input file
    if (!input_files.empty()) {
      try {
        H5::H5File first_file(input_files[0], H5F_ACC_RDONLY);
        if (H5Lexists(first_file.getId(), "/metadata", H5P_DEFAULT)) {
          H5Ocopy(first_file.getId(), "/metadata", out_file.getId(),
                  "/metadata", H5P_DEFAULT, H5P_DEFAULT);
          std::cout << "  Copied /metadata registries from " << input_files[0]
                    << std::endl;
        }
      } catch (const H5::Exception& e) {
        std::cerr << "Warning: Could not copy metadata: " << e.getDetailMsg()
                  << std::endl;
      }
    }

    std::cout << "\nMerge complete!" << std::endl;
  } catch (const H5::Exception& e) {
    std::cerr << "Error writing merged file: " << e.getDetailMsg() << std::endl;
    throw std::runtime_error("Failed to write merged event file");
  }
}

}  // namespace june
