#pragma once

#include <string>
#include <vector>

#include "core/config.h"  // SelectionCriterion
#include "core/types.h"   // Person

namespace june {

struct WorldState;

/// @namespace june::filtering
/// @brief CSV-driven row filtering used by outcome-rate tables and infection
/// seeds.
///
/// ## Overview
///
/// A CSV file can declare zero or more **filter columns**: columns whose
/// header begins with `filter.`. Each row then applies only to persons (and
/// infection events) that satisfy *all* non-empty filter cells on that row.
/// Rows are evaluated top-to-bottom; the first matching row wins.
///
/// An empty filter cell is silently skipped (no constraint imposed). A row with
/// all filter cells empty matches every person unconditionally and acts as a
/// catch-all default.
///
/// ## Convention
///
/// A `filter.X` header means "this cell is a direct comparison against person
/// (or infection-context) attribute `X`." Columns that need caller-side joins
/// or lookups (e.g. resolving a `geo_level` + `geo_unit` pair against the
/// world graph to emit an `in` criterion over descendant ids) MUST NOT use
/// the `filter.` prefix. They are value columns interpreted by the caller.
///
/// ## Available filter columns
///
/// ### Person attributes
///
/// | Header                        | Matches against       | Value format |
/// |-------------------------------|-----------------------|---------------------------------------|
/// | `filter.age`                  | `person.age`          | Exact int, range,
/// or comparison       | | `filter.sex`                  | `person.sex` |
/// `male` or `female`                    | | `filter.geo_unit_id`          |
/// `person.geo_unit_id`  | Integer ID                            | |
/// `filter.properties.<name>`    | named custom property | Integer, bool, or
/// string              | | `filter.activities.<act>.<p>` | activity
/// sub-property | Integer or string                     | |
/// `filter.networks.<net>.<p>`   | network sub-property  | Integer or string |
///
/// ### Infection-event context
///
/// | Header                        | Matches against                         |
/// Value format          |
/// |-------------------------------|-----------------------------------------|-----------------------|
/// | `filter.infector_symptom`     | symptom-tag name of the infector        |
/// String (exact match)  | | `filter.transmission_mode`    | transmission mode
/// that caused infection | String (exact match)  |
///
/// For `filter.infector_symptom` and `filter.transmission_mode`, a non-empty
/// filter cell **fails** when the corresponding context field is empty (i.e.
/// for seeded infections with no explicit infector). Rows that should match
/// seeds must leave those cells blank.
///
/// ## Value formats
///
/// | Format       | Example        | Meaning                          |
/// |--------------|----------------|----------------------------------|
/// | Exact int    | `30`           | `age == 30`                      |
/// | Range        | `18-64`        | `18 <= age <= 64`                |
/// | Comparison   | `>=65`         | `age >= 65`                      |
/// | Comparison   | `<18`          | `age < 18`                       |
/// | Boolean      | `true`         | property equals true             |
/// | String       | `male`         | string equality                  |
///
/// ## Example CSV
///
/// ```csv
/// filter.age, filter.sex, filter.transmission_mode, trajectory_a, trajectory_b
/// # Respiratory transmission → always trajectory_b
/// , , respiratory, 0.0, 1.0
/// # Animal-bite transmission, older males → higher mortality in trajectory_a
/// 50-99, male, animal_bite, 0.8, 0.2
/// # Animal-bite transmission, all others
/// , , animal_bite, 0.6, 0.4
/// # Default catch-all (seeds / no infector)
/// , , , 0.6, 0.4
/// ```
namespace filtering {

/// Parses a single filter key/value pair into one or two SelectionCriterion
/// objects.
///
/// A numeric range string such as `"18-60"` produces two criteria (`>= 18` AND
/// `<= 60`). Comparison prefixes (`>`, `>=`, `<`, `<=`) produce one criterion.
/// Plain integers, booleans, and strings produce an `==` criterion. The legacy
/// alias `"age_groups"` is silently mapped to `"age"`.
///
/// Args:
///   key: The property path (everything after `filter.` in the CSV header).
///   val: The cell value string from the CSV row.
///
/// Returns:
///   A vector of SelectionCriterion (one for plain values, two for ranges).
std::vector<SelectionCriterion> parseCriterionFromKeyValue(
    const std::string& key, const std::string& val);

/// Returns true if every criterion in `criteria` evaluates to true for
/// `person`.
///
/// Criteria are evaluated conjunctively (AND). An empty criteria list matches
/// all persons. Criteria with `property_path == "infector_symptom"` or
/// `"transmission_mode"` are matched against the corresponding field of `ctx`
/// rather than person attributes; a non-empty criterion on either field fails
/// when the context field is empty (e.g. seeded infections).
///
/// Args:
///   person:   The person being evaluated.
///   world:    Pointer to WorldState, used for property lookups. May be null
///             if no world-dependent properties are needed.
///   criteria: The list of criteria to evaluate (typically from
///   parseCriteriaFromRow). ctx:      Infection-event context (infector
///   symptom, transmission mode).
///             Defaults to empty strings (no context).
///
/// Returns:
///   True if all criteria pass; false if any criterion fails.
bool matchesCriteria(const Person& person, const WorldState* world,
                     const std::vector<SelectionCriterion>& criteria,
                     const InfectionContext& ctx = {});

/// Scans CSV headers and returns (column_index, property_path) for every
/// column whose header begins with `"filter."`.
///
/// Args:
///   headers: The full list of trimmed header strings from the CSV.
///
/// Returns:
///   A vector of (column_index, property_path) pairs where property_path is
///   the substring after `"filter."` (e.g. `"age"`, `"transmission_mode"`).
std::vector<std::pair<int, std::string>> findFilterColumns(
    const std::vector<std::string>& headers);

/// Parses an AND-conjunctive expression string into a list of criteria.
///
/// Format: `<key><op><value> [AND <key><op><value> ...]`
///   - Tokens are separated by literal " AND " (case-sensitive, surrounding
///     spaces required).
///   - Operators: `==`, `=`, `>=`, `<=`, `>`, `<`. (`!=` not supported.)
///   - Keys may carry a `filter.` prefix (stripped before parsing) or be a
///     bare path such as `is_alive` or a person property name.
///   - Values follow the same conventions as CSV cells: bool (`true`/`false`),
///     int, range (`16-59`), comparison (carried in val for `>=`/`<=`/`>`/`<`).
///
/// The criteria within one expression are evaluated conjunctively (AND).
/// To express OR, use a list of expressions and call matchesAnyGroup.
///
/// Throws std::runtime_error on malformed input.
std::vector<SelectionCriterion> parseConjunctiveExpression(
    const std::string& expr);

/// Returns true if `criteria_groups` is empty, or if any group fully matches.
/// Each group is itself AND-conjunctive (matchesCriteria semantics); the
/// outer combinator is OR. Used for eligibility-style filters where today's
/// logic depends on a partner-state-conditional choice between flag sets.
bool matchesAnyGroup(
    const Person& person, const WorldState* world,
    const std::vector<std::vector<SelectionCriterion>>& criteria_groups);

/// Builds a list of SelectionCriterion from the filter columns of a single CSV
/// row.
///
/// Empty cells are skipped; they impose no constraint. Only non-empty cells
/// contribute a criterion, meaning the row matches all values for that
/// property.
///
/// Args:
///   fields:      The full list of trimmed cell strings for one CSV row.
///   filter_cols: The (column_index, property_path) pairs from
///   findFilterColumns.
///
/// Returns:
///   A vector of SelectionCriterion representing the conjunction of all
///   non-empty filter cells in the row.
std::vector<SelectionCriterion> parseCriteriaFromRow(
    const std::vector<std::string>& fields,
    const std::vector<std::pair<int, std::string>>& filter_cols);

}  // namespace filtering
}  // namespace june
