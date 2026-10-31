#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/types.h"
#include "epidemiology/observed_metrics.h"

namespace june {

// The state of a trigger at one unit: whether it holds the Policy in force
// there, and one (day, value) pair that only a cumulative metric uses. While
// open the value is the baseline (0 at the start, then the value at the last
// reopening) and the day is not read. While closed the pair is the anchor: the
// last time the metric grew beyond the tolerance, or the moment of closing if
// it has not.
struct TriggerCell {
  bool closed = false;
  int64_t anchor_day = 0;
  int64_t anchor_value = 0;

  bool operator==(const TriggerCell&) const = default;
};

// One term of a response: a metric's reduced count times a weight.
struct ResponseTerm {
  std::string metric;
  double weight = 0.0;
};

// Ties a Policy's activation to an observed metric. "Closed" means the trigger
// holds the Policy in force, the way a venue closes; "open" means it does not.
//
// A trigger takes one of two forms. The threshold form closes and reopens as
// the metric crosses its thresholds. The response form (voluntary) keeps no
// open or closed state: a daily avoidance level is derived from the metric
// counts, and each decision unit avoids exactly when the level is above its own
// fixed draw.
//
// The state is read once per day, after the metric is reduced, so every rank
// reaches the same answer from the same integers. A world scope trigger keeps
// one cell, the inherited one. A trigger on a metric counted per unit keeps one
// cell per unit it has touched, and a bit per unit for the closed ones so the
// per slot question "is this unit closed" is a single lookup.
struct PolicyTrigger : TriggerCell {
  // Configuration. activate_at is a value of a current metric, or growth since
  // the baseline for a cumulative one. deactivate_at belongs to current
  // metrics, calm_days and the two tolerances to cumulative ones; the loader
  // rejects the other mix.
  std::string metric;
  int64_t activate_at = 0;
  std::optional<int64_t> deactivate_at;
  int64_t calm_days = 0;
  int64_t tolerance = 0;
  int64_t tolerance_pct = 0;  // whole percent of the anchor value
  bool cumulative = false;    // taken from the metric's source at load
  MetricScope scope = MetricScope::World;  // taken from the metric at load

  bool perUnit() const { return scope != MetricScope::World; }

  // Response form. level = base + sum(weight * count) over the terms, clamped
  // to 1 where it is read. The loader keeps base within 0 to 1 and every weight
  // above 0, so there is no lower clamp. All terms share one scope and one
  // activity.
  double base = 0.0;
  std::vector<ResponseTerm> terms;
  DecisionUnit decision_unit;
  int16_t metric_activity_index = -1;  // the terms' activity, set on resolve

  bool voluntary() const { return !terms.empty(); }

  // Does the trigger hold the Policy in force as a whole? A threshold trigger
  // on a world scope metric does while it is closed. A trigger counted per
  // unit and a voluntary one always do, and then decide slot by slot where the
  // Policy acts.
  bool holdsPolicy() const { return decidesPerSlot() || closed; }
  bool decidesPerSlot() const { return perUnit() || voluntary(); }

  // The level as of the last day boundary, before the clamp. A world scope
  // trigger has one; a trigger counted per unit has one for each unit some
  // term counted, and every other unit sits at the base.
  double world_level = 0.0;
  std::unordered_map<VenueId, double> unit_levels;

  double levelAtWorld() const { return std::min(1.0, world_level); }
  double levelAtUnit(VenueId unit) const {
    const auto it = unit_levels.find(unit);
    return std::min(1.0, it == unit_levels.end() ? base : it->second);
  }

  // One evaluation of the world scope state. `day` is the day that just ended
  // and `value` the reduced metric at its end.
  void step(int64_t day, int64_t value) { step(*this, day, value); }
  // The same rule applied to any cell.
  void step(TriggerCell& cell, int64_t day, int64_t value) const;

  // One evaluation at one unit. A cell is kept only while it differs from the
  // untouched state, so a unit that was never close to the threshold costs
  // nothing.
  void stepUnit(VenueId unit, int64_t day, int64_t value);
  bool unitClosed(VenueId unit) const;
  // Units whose cell says closed, in unit order.
  std::vector<VenueId> closedUnits() const;
  const std::map<VenueId, TriggerCell>& unitCells() const { return cells_; }
  // Replace every unit cell, as a checkpoint restore does.
  void restoreUnits(const std::vector<std::pair<VenueId, TriggerCell>>& cells);

 private:
  void setUnitBit(VenueId unit, bool closed);

  std::map<VenueId, TriggerCell> cells_;
  std::vector<uint64_t> closed_bits_;  // indexed by unit id
};

}  // namespace june
