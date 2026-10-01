/// @file filtering.cpp
/// @brief Implementation of CSV-driven row filtering for outcome-rate tables
/// and
///        infection seeds.
///
/// See filtering.h for the full API documentation and a description of all
/// supported `filter.*` column types and value formats.
#include "utils/filtering.h"

#include <iostream>
#include <stdexcept>

namespace june {
namespace filtering {

namespace {

std::string trimSpaces(const std::string& s) {
  size_t a = s.find_first_not_of(" \t");
  if (a == std::string::npos) return "";
  size_t b = s.find_last_not_of(" \t");
  return s.substr(a, b - a + 1);
}

std::string stripFilterPrefix(const std::string& key) {
  static const std::string p = "filter.";
  if (key.size() > p.size() && key.compare(0, p.size(), p) == 0) {
    return key.substr(p.size());
  }
  return key;
}

// Find the first operator occurrence in `tok` (left-to-right scan). At each
// position multi-char ops are tried first so "==" isn't mis-split as "=".
// Returns (pos, op_length); pos==npos if none.
std::pair<size_t, size_t> findFirstOp(const std::string& tok) {
  static const std::vector<std::string> long_ops = {">=", "<=", "==", "!="};
  for (size_t i = 0; i < tok.size(); ++i) {
    for (const auto& op : long_ops) {
      if (i + op.size() <= tok.size() && tok.compare(i, op.size(), op) == 0) {
        return {i, op.size()};
      }
    }
    char c = tok[i];
    if (c == '>' || c == '<' || c == '=') return {i, 1};
  }
  return {std::string::npos, 0};
}

}  // namespace

std::vector<SelectionCriterion> parseConjunctiveExpression(
    const std::string& expr) {
  std::vector<SelectionCriterion> out;
  std::string s = trimSpaces(expr);
  if (s.empty()) return out;

  // Split on " AND ".
  std::vector<std::string> tokens;
  size_t pos = 0;
  while (pos < s.size()) {
    size_t and_pos = s.find(" AND ", pos);
    std::string tok = (and_pos == std::string::npos)
                          ? s.substr(pos)
                          : s.substr(pos, and_pos - pos);
    tok = trimSpaces(tok);
    if (!tok.empty()) tokens.push_back(tok);
    if (and_pos == std::string::npos) break;
    pos = and_pos + 5;
  }

  for (const auto& tok : tokens) {
    auto [op_pos, op_len] = findFirstOp(tok);
    if (op_pos == std::string::npos) {
      throw std::runtime_error(
          "filtering::parseConjunctiveExpression: no operator in token '" +
          tok + "'");
    }
    std::string key = trimSpaces(tok.substr(0, op_pos));
    std::string op = tok.substr(op_pos, op_len);
    std::string val = trimSpaces(tok.substr(op_pos + op_len));
    if (key.empty() || val.empty()) {
      throw std::runtime_error(
          "filtering::parseConjunctiveExpression: malformed token '" + tok +
          "'");
    }
    key = stripFilterPrefix(key);

    // Forward to parseCriterionFromKeyValue with op encoded in val for the
    // comparison ops (its existing convention). Plain == / = → bare val.
    if (op == "==" || op == "=") {
      auto cs = parseCriterionFromKeyValue(key, val);
      out.insert(out.end(), cs.begin(), cs.end());
    } else if (op == "!=") {
      // parseCriterionFromKeyValue does not produce "!=", so emit directly.
      auto cs = parseCriterionFromKeyValue(key, val);
      for (auto& c : cs) c.operator_type = "!=";
      out.insert(out.end(), cs.begin(), cs.end());
    } else {
      // >=, <=, >, <: parseCriterionFromKeyValue accepts these as a val
      // prefix.
      auto cs = parseCriterionFromKeyValue(key, op + val);
      out.insert(out.end(), cs.begin(), cs.end());
    }
  }
  return out;
}

bool matchesAnyGroup(
    const Person& person, const WorldState* world,
    const std::vector<std::vector<SelectionCriterion>>& criteria_groups) {
  if (criteria_groups.empty()) return true;
  for (const auto& group : criteria_groups) {
    if (matchesCriteria(person, world, group)) return true;
  }
  return false;
}

std::vector<SelectionCriterion> parseCriterionFromKeyValue(
    const std::string& key, const std::string& val) {
  std::string property = key;
  // Map legacy alias used in bulk seed CSV headers
  if (property == "age_groups") {
    property = "age";
    std::cout << "  [Config] Mapping property 'age_groups' -> 'age'"
              << std::endl;
  }

  std::vector<SelectionCriterion> results;

  // Support comparison prefixes: ">12", ">=12", "<64", "<=64"
  if (!val.empty() && (val[0] == '<' || val[0] == '>')) {
    std::string op;
    size_t num_start = 1;
    if (val.size() > 1 && val[1] == '=') {
      op = (val[0] == '>') ? ">=" : "<=";
      num_start = 2;
    } else {
      op = (val[0] == '>') ? ">" : "<";
    }
    try {
      SelectionCriterion c;
      c.property_path = property;
      c.operator_type = op;
      c.value = std::stoi(val.substr(num_start));
      return {c};
    } catch (...) {
      // Fall through to other parsers
    }
  }

  // Support ranges like "18-30" or "65-100"
  size_t dash = val.find('-');
  if (dash != std::string::npos && dash > 0 && dash < val.size() - 1) {
    try {
      SelectionCriterion c_min, c_max;
      c_min.property_path = property;
      c_min.operator_type = ">=";
      c_min.value = std::stoi(val.substr(0, dash));

      c_max.property_path = property;
      c_max.operator_type = "<=";
      c_max.value = std::stoi(val.substr(dash + 1));

      results.push_back(c_min);
      results.push_back(c_max);
      return results;
    } catch (...) {
      // Fall through to single-value parse
    }
  }

  SelectionCriterion c;
  c.property_path = property;
  c.operator_type = "==";

  // Try numeric/bool first, fallback to string
  try {
    if (val == "true" || val == "True")
      c.value = true;
    else if (val == "false" || val == "False")
      c.value = false;
    else if (val.find('.') != std::string::npos) {
      c.value = std::stod(val);
    } else {
      c.value = std::stoi(val);
    }
  } catch (...) {
    c.value = val;
  }

  results.push_back(c);
  return results;
}

bool matchesCriteria(const Person& person, const WorldState* world,
                     const std::vector<SelectionCriterion>& criteria,
                     const InfectionContext& ctx) {
  for (const SelectionCriterion& criterion : criteria) {
    if (!criterion.evaluate(person, world, nullptr, &ctx)) return false;
  }
  return true;
}

std::vector<std::pair<int, std::string>> findFilterColumns(
    const std::vector<std::string>& headers) {
  std::vector<std::pair<int, std::string>> result;
  for (int i = 0; i < (int)headers.size(); ++i) {
    if (headers[i].find("filter.") == 0) {
      result.push_back({i, headers[i].substr(7)});
    }
  }
  return result;
}

std::vector<SelectionCriterion> parseCriteriaFromRow(
    const std::vector<std::string>& fields,
    const std::vector<std::pair<int, std::string>>& filter_cols) {
  std::vector<SelectionCriterion> criteria;
  for (const auto& [col_idx, property_path] : filter_cols) {
    if (col_idx < (int)fields.size() && !fields[col_idx].empty()) {
      auto cs = parseCriterionFromKeyValue(property_path, fields[col_idx]);
      criteria.insert(criteria.end(), cs.begin(), cs.end());
    }
  }
  return criteria;
}

}  // namespace filtering
}  // namespace june
