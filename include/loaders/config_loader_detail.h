#pragma once

#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/config.h"

// Forward-declare YAML::Node to keep yaml-cpp out of the public include
// surface. Only TUs that actually parse YAML need the full header.
namespace YAML {
class Node;
}  // namespace YAML

namespace june {
namespace config_detail {

class MalformedSelectionError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Parse a YAML sequence of `{property, operator, value}` entries into a
// vector of SelectionCriterion. Scalar values use the shared bool -> int ->
// double -> string dispatch; sequence values are validated by the loader's
// local sequence parser.
// Shared by schedule, activity preferences, vaccination, encounters, and
// policy loading.
void parseSelectionCriteria(const YAML::Node& selection_node,
                            std::vector<SelectionCriterion>& out);

// Throws on the first key of `node` that is not one of `fields`, naming it and
// listing the ones that are. `within` says where the key was found when `where`
// does not already.
void rejectUnknownKeys(const YAML::Node& node,
                       std::initializer_list<const char*> fields,
                       const std::string& where,
                       const std::string& within = "");

// The same check for a list of fields only known at run time, such as the day
// type names a schedule file declares.
void rejectUnknownKeysAmong(const YAML::Node& node,
                            const std::vector<std::string>& fields,
                            const std::string& where,
                            const std::string& within = "");

}  // namespace config_detail
}  // namespace june
