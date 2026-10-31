#include "epidemiology/policy_roster/policy_roster.h"

#include <stdexcept>
#include <unordered_map>

namespace june {

namespace {

constexpr PolicyKind kKinds[] = {PolicyKind::Symptom, PolicyKind::Temporal};

// Each name's index; throws on a repeat.
std::unordered_map<std::string, size_t> indexByName(
    const std::vector<std::string>& names, PolicyKind kind) {
  std::unordered_map<std::string, size_t> index_by_name;
  for (size_t index = 0; index < names.size(); ++index) {
    const auto [existing, inserted] =
        index_by_name.emplace(names[index], index);
    if (inserted) continue;
    throw std::runtime_error(std::string(policyKindName(kind)) + " policies " +
                             std::to_string(existing->second) + " and " +
                             std::to_string(index) + " share the name '" +
                             names[index] +
                             "'; policy names must be unique within a kind");
  }
  return index_by_name;
}

KindRosterMatch matchKind(const std::vector<std::string>& stored,
                          const std::vector<std::string>& current,
                          PolicyKind kind) {
  indexByName(stored, kind);
  const auto current_index_by_name = indexByName(current, kind);

  KindRosterMatch kind_match;
  std::vector<bool> kept(current.size(), false);
  for (size_t stored_index = 0; stored_index < stored.size(); ++stored_index) {
    const auto found = current_index_by_name.find(stored[stored_index]);
    if (found == current_index_by_name.end()) {
      kind_match.current_index.push_back(std::nullopt);
      kind_match.dropped.push_back(stored[stored_index]);
      continue;
    }
    const size_t current_index = found->second;
    kind_match.current_index.push_back(current_index);
    kept[current_index] = true;
    if (current_index != stored_index)
      kind_match.moved.push_back(
          {stored[stored_index], stored_index, current_index});
  }
  for (size_t current_index = 0; current_index < current.size();
       ++current_index)
    if (!kept[current_index])
      kind_match.added.push_back(current[current_index]);
  return kind_match;
}

std::string quotedList(const std::vector<std::string>& names) {
  std::string text;
  for (const std::string& name : names)
    text += (text.empty() ? "'" : ", '") + name + "'";
  return text;
}

}  // namespace

void rejectDuplicateNames(const PolicyRoster& roster) {
  for (PolicyKind kind : kKinds) indexByName(roster.of(kind), kind);
}

PolicyRosterMatch match(const PolicyRoster& stored,
                        const PolicyRoster& current) {
  return PolicyRosterMatch(
      matchKind(stored.symptom, current.symptom, PolicyKind::Symptom),
      matchKind(stored.temporal, current.temporal, PolicyKind::Temporal));
}

bool PolicyRosterMatch::identical() const {
  for (PolicyKind kind : kKinds) {
    const KindRosterMatch& kind_match = of(kind);
    if (!kind_match.dropped.empty() || !kind_match.added.empty() ||
        !kind_match.moved.empty())
      return false;
  }
  return true;
}

std::optional<size_t> PolicyRosterMatch::currentIndex(
    PolicyKind kind, size_t stored_index) const {
  const auto& current_index = of(kind).current_index;
  if (stored_index >= current_index.size())
    throw std::out_of_range("stored " + std::string(policyKindName(kind)) +
                            " policy " + std::to_string(stored_index) +
                            " is past the stored roster of " +
                            std::to_string(current_index.size()));
  return current_index[stored_index];
}

PolicyMask PolicyRosterMatch::remap(const PolicyMask& stored_mask,
                                    PolicyKind kind) const {
  PolicyMask current_mask;
  for (size_t stored_index = 0; stored_index < kMaxPoliciesPerKind;
       ++stored_index) {
    if (!stored_mask.test(stored_index)) continue;
    if (const auto current_index = currentIndex(kind, stored_index))
      current_mask.set(*current_index);
  }
  return current_mask;
}

std::string PolicyRosterMatch::describeChanges() const {
  std::string text;
  for (PolicyKind kind : kKinds) {
    const KindRosterMatch& kind_match = of(kind);
    std::vector<std::string> parts;
    if (!kind_match.dropped.empty())
      parts.push_back("dropped " + quotedList(kind_match.dropped));
    if (!kind_match.added.empty())
      parts.push_back("added " + quotedList(kind_match.added));
    if (!kind_match.moved.empty()) {
      std::string moved = "moved";
      for (size_t index = 0; index < kind_match.moved.size(); ++index) {
        const KindRosterMatch::Move& move = kind_match.moved[index];
        moved += (index == 0 ? " '" : ", '") + move.name + "' " +
                 std::to_string(move.stored_index) + "->" +
                 std::to_string(move.current_index);
      }
      parts.push_back(moved);
    }
    if (parts.empty()) continue;
    if (!text.empty()) text += "\n";
    text += std::string(policyKindName(kind)) + " policies: ";
    for (size_t index = 0; index < parts.size(); ++index)
      text += (index == 0 ? "" : "; ") + parts[index];
  }
  return text;
}

}  // namespace june
