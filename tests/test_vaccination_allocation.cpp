#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

#include "doctest.h"
#include "epidemiology/vaccination_allocation.h"
#include "utils/deterministic_rng.h"

using namespace june;

namespace {

constexpr int64_t kMax64 = std::numeric_limits<int64_t>::max();
constexpr uint64_t kMaxKey = std::numeric_limits<uint64_t>::max();

int64_t sumAllotments(int64_t doses, int64_t per_days, int64_t days) {
  int64_t total = 0;
  for (int64_t d = 0; d < days; ++d) total += doseAllotment(doses, per_days, d);
  return total;
}

std::vector<int64_t> fill(std::optional<int64_t> allotment,
                          std::vector<int64_t> counts,
                          std::vector<bool> starts_course,
                          std::optional<int64_t> cap) {
  return fillTiers(allotment, counts, starts_course, cap);
}

// Runs lotteryCutoff over a world held in `parts` (one sorted vector per
// rank). The count it is given adds up every part, which is what the global
// sum of a local upper_bound gives in an MPI run. Returns the cutoff and how
// many times the count was asked for.
std::pair<std::optional<uint64_t>, size_t> runCutoff(
    const std::vector<std::vector<uint64_t>>& parts, int64_t take) {
  int64_t total = 0;
  for (const auto& part : parts) total += static_cast<int64_t>(part.size());
  size_t calls = 0;
  const auto cutoff = lotteryCutoff(total, take, [&](uint64_t t) {
    ++calls;
    int64_t count = 0;
    for (const auto& part : parts)
      count += std::upper_bound(part.begin(), part.end(), t) - part.begin();
    return count;
  });
  return {cutoff, calls};
}

// Spreads the keys over `ranks` ranks by a rule that has nothing to do with
// their order, so each rank holds an arbitrary slice. Each part is sorted.
std::vector<std::vector<uint64_t>> splitKeys(const std::vector<uint64_t>& keys,
                                             size_t ranks, uint64_t seed) {
  std::vector<std::vector<uint64_t>> parts(ranks);
  for (size_t i = 0; i < keys.size(); ++i)
    parts[mix_seed(seed, i) % ranks].push_back(keys[i]);
  for (auto& p : parts) std::sort(p.begin(), p.end());
  return parts;
}

std::vector<uint64_t> randomUniqueKeys(size_t n, uint64_t seed) {
  SplitMix64 g(seed);
  std::vector<uint64_t> keys(n);
  for (auto& k : keys) k = g();
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  return keys;
}

}  // namespace

// ---------------------------------------------------------------------------
// doseAllotment
// ---------------------------------------------------------------------------

TEST_CASE("doseAllotment - whole periods add up to exactly the doses") {
  const std::vector<std::pair<int64_t, int64_t>> budgets = {
      {10, 3}, {7, 30}, {1, 5}, {30, 7}, {100, 100}, {1000003, 31}, {5, 1}};
  for (const auto& [doses, per_days] : budgets) {
    CAPTURE(doses);
    CAPTURE(per_days);
    for (int64_t periods : {1, 2, 3, 10})
      CHECK(sumAllotments(doses, per_days, periods * per_days) ==
            periods * doses);
  }
}

TEST_CASE("doseAllotment - no drift over many periods") {
  // 10 doses every 7 days for 1000 periods: the schedule must land on
  // exactly 10000 doses, not a few short or over from rounding.
  CHECK(sumAllotments(10, 7, 7000) == 10000);
  // 3 doses every 7 days is below one a day, so most days give none.
  CHECK(sumAllotments(3, 7, 7000) == 3000);
}

TEST_CASE("doseAllotment - days differ by at most one dose") {
  const std::vector<std::pair<int64_t, int64_t>> budgets = {
      {10, 3}, {7, 30}, {30, 7}, {1000003, 31}, {3, 7}, {6, 6}};
  for (const auto& [doses, per_days] : budgets) {
    CAPTURE(doses);
    CAPTURE(per_days);
    int64_t lowest = kMax64;
    int64_t highest = 0;
    for (int64_t d = 0; d < 3 * per_days; ++d) {
      const int64_t a = doseAllotment(doses, per_days, d);
      lowest = std::min(lowest, a);
      highest = std::max(highest, a);
    }
    CHECK(highest - lowest <= 1);
    // Never more than the rate rounded up, never less than it rounded down.
    CHECK(lowest >= doses / per_days);
    CHECK(highest <= (doses + per_days - 1) / per_days);
  }
}

TEST_CASE("doseAllotment - a budget as large as the arguments allow") {
  CHECK(doseAllotment(kMax64, kMax64, 0) == 1);
  CHECK(doseAllotment(kMax64, 1, 0) == kMax64);
}

TEST_CASE("doseAllotment - bad inputs throw") {
  CHECK_THROWS_WITH_AS(doseAllotment(0, 30, 0), doctest::Contains("doses"),
                       std::runtime_error);
  CHECK_THROWS_WITH_AS(doseAllotment(-5, 30, 0), doctest::Contains("doses"),
                       std::runtime_error);
  CHECK_THROWS_WITH_AS(doseAllotment(10, 0, 0), doctest::Contains("per_days"),
                       std::runtime_error);
  CHECK_THROWS_WITH_AS(doseAllotment(10, -1, 0), doctest::Contains("per_days"),
                       std::runtime_error);
  CHECK_THROWS_WITH_AS(doseAllotment(10, 30, -1), doctest::Contains("day"),
                       std::runtime_error);
  // The schedule needs (day + 1) * doses in 64 bits.
  CHECK_THROWS_WITH_AS(doseAllotment(kMax64, 30, 1),
                       doctest::Contains("overflows"), std::runtime_error);
  CHECK_THROWS_WITH_AS(doseAllotment(2, 30, kMax64 / 2),
                       doctest::Contains("overflows"), std::runtime_error);
  CHECK_THROWS_WITH_AS(doseAllotment(1, 30, kMax64),
                       doctest::Contains("overflows"), std::runtime_error);
  CHECK_NOTHROW(doseAllotment(2, 30, kMax64 / 2 - 1));
}

// ---------------------------------------------------------------------------
// fillTiers
// ---------------------------------------------------------------------------

TEST_CASE("fillTiers - tiers are served in strict order") {
  // Enough for the first tier and part of the second; the third gets nothing
  // even though it is the largest.
  CHECK(fill(12, {5, 10, 100}, {true, true, true}, std::nullopt) ==
        std::vector<int64_t>{5, 7, 0});
  CHECK(fill(0, {5, 10}, {true, true}, std::nullopt) ==
        std::vector<int64_t>{0, 0});
  CHECK(fill(1000, {5, 10}, {true, true}, std::nullopt) ==
        std::vector<int64_t>{5, 10});
  CHECK(fill(5, {}, {}, std::nullopt).empty());
}

TEST_CASE("fillTiers - a tier with no candidates passes its turn on") {
  CHECK(fill(6, {0, 4, 0, 4}, {true, true, true, true}, std::nullopt) ==
        std::vector<int64_t>{0, 4, 0, 2});
}

TEST_CASE("fillTiers - the coverage cap only cuts tiers that start a course") {
  // Tier 0 starts courses and is cut to the 4 the cap allows. Tier 1 is a
  // later dose of courses already begun, so the cap does not touch it and it
  // still takes budget. Tier 2 starts courses and the cap is now spent.
  CHECK(fill(100, {10, 10, 10}, {true, false, true}, 4) ==
        std::vector<int64_t>{4, 10, 0});
}

TEST_CASE(
    "fillTiers - the cap and the budget both apply to a first dose tier") {
  // Budget 6, cap 4: the first tier is held to 4 by the cap, the second tier
  // is held to the 2 doses left in the budget.
  CHECK(fill(6, {10, 10}, {true, false}, 4) == std::vector<int64_t>{4, 2});
  // Budget tighter than the cap.
  CHECK(fill(3, {10, 10}, {true, true}, 8) == std::vector<int64_t>{3, 0});
}

TEST_CASE("fillTiers - a cap of zero stops first doses but not later doses") {
  CHECK(fill(50, {10, 10}, {true, false}, 0) == std::vector<int64_t>{0, 10});
}

TEST_CASE("fillTiers - no budget serves everyone, within the cap") {
  CHECK(fill(std::nullopt, {5, 7, 9}, {true, false, true}, std::nullopt) ==
        std::vector<int64_t>{5, 7, 9});
  CHECK(fill(std::nullopt, {5, 7, 9}, {true, false, true}, 8) ==
        std::vector<int64_t>{5, 7, 3});
  // A count near the top of the range must not wrap around.
  CHECK(fill(std::nullopt, {kMax64, 1}, {false, false}, std::nullopt) ==
        std::vector<int64_t>{kMax64, 1});
}

TEST_CASE("fillTiers - invariants hold for random inputs") {
  SplitMix64 g(2024);
  for (int trial = 0; trial < 2000; ++trial) {
    const size_t tiers = 1 + bounded(g, 6);
    std::vector<int64_t> counts(tiers);
    std::vector<bool> starts(tiers);
    for (size_t t = 0; t < tiers; ++t) {
      counts[t] = static_cast<int64_t>(bounded(g, 50));
      starts[t] = bounded(g, 2) == 0;
    }
    const std::optional<int64_t> budget =
        bounded(g, 4) == 0 ? std::nullopt
                           : std::optional<int64_t>(bounded(g, 120));
    const std::optional<int64_t> cap =
        bounded(g, 3) == 0 ? std::nullopt
                           : std::optional<int64_t>(bounded(g, 60));
    CAPTURE(trial);
    const auto taken = fillTiers(budget, counts, starts, cap);
    REQUIRE(taken.size() == tiers);

    int64_t total = 0;
    int64_t courses = 0;
    for (size_t t = 0; t < tiers; ++t) {
      CHECK(taken[t] >= 0);
      CHECK(taken[t] <= counts[t]);
      total += taken[t];
      if (starts[t]) courses += taken[t];
    }
    if (budget) CHECK(total <= *budget);
    if (cap) CHECK(courses <= *cap);

    // Strict order, stated without reference to how the loop is written: a
    // tier that was cut short must have hit a limit right there, either no
    // budget left or, for a tier that starts courses, no cap left.
    int64_t spent = 0;
    int64_t started = 0;
    for (size_t t = 0; t < tiers; ++t) {
      spent += taken[t];
      if (starts[t]) started += taken[t];
      if (taken[t] < counts[t]) {
        const bool out_of_budget = budget && spent == *budget;
        const bool out_of_cap = starts[t] && cap && started == *cap;
        CHECK((out_of_budget || out_of_cap));
      }
    }
  }
}

TEST_CASE("fillTiers - bad inputs throw") {
  CHECK_THROWS_WITH_AS(fill(5, {1, 2}, {true}, std::nullopt),
                       doctest::Contains("starts_course"), std::runtime_error);
  CHECK_THROWS_WITH_AS(fill(5, {1, -2}, {true, true}, std::nullopt),
                       doctest::Contains("negative count"), std::runtime_error);
  CHECK_THROWS_WITH_AS(fill(-1, {1}, {true}, std::nullopt),
                       doctest::Contains("allotment"), std::runtime_error);
  CHECK_THROWS_WITH_AS(fill(5, {1}, {true}, -1),
                       doctest::Contains("cap_remaining"), std::runtime_error);
}

// ---------------------------------------------------------------------------
// lotteryKey
// ---------------------------------------------------------------------------

TEST_CASE("lotteryKey - unique per person for a fixed seed, campaign, tier") {
  constexpr uint64_t kPeople = 2'000'000;
  std::vector<uint64_t> keys(kPeople);
  for (uint64_t id = 0; id < kPeople; ++id) keys[id] = lotteryKey(7, id, 3, 2);
  std::sort(keys.begin(), keys.end());
  CHECK(std::adjacent_find(keys.begin(), keys.end()) == keys.end());
}

TEST_CASE("lotteryKey - a different seed, campaign or tier reorders people") {
  constexpr uint64_t kPeople = 2000;
  const auto order = [&](uint64_t seed, uint64_t campaign, uint64_t tier) {
    std::vector<uint64_t> ids(kPeople);
    std::iota(ids.begin(), ids.end(), uint64_t{0});
    std::sort(ids.begin(), ids.end(), [&](uint64_t a, uint64_t b) {
      return lotteryKey(seed, a, campaign, tier) <
             lotteryKey(seed, b, campaign, tier);
    });
    return ids;
  };
  const auto base = order(42, 1, 0);
  CHECK(order(42, 1, 0) == base);
  CHECK(order(43, 1, 0) != base);
  CHECK(order(42, 2, 0) != base);
  CHECK(order(42, 1, 1) != base);
  // The order is not the order of the ids.
  std::vector<uint64_t> by_id(kPeople);
  std::iota(by_id.begin(), by_id.end(), uint64_t{0});
  CHECK(base != by_id);
}

TEST_CASE("lotteryKey - one pinned value, so a change to the mix is noticed") {
  // Worked out outside the engine from the published SplitMix64 constants.
  CHECK(lotteryKey(42, 1000, 2, 1) == 0x6dcf52eedaec55d2ULL);
}

// ---------------------------------------------------------------------------
// lotteryCutoff
// ---------------------------------------------------------------------------

namespace {

// Everything the cutoff says should be true, checked against the plainest
// possible serial selection: sort all the keys and take the first `take`.
void checkCutoff(const std::vector<uint64_t>& keys, size_t ranks, int64_t take,
                 uint64_t split_seed) {
  CAPTURE(keys.size());
  CAPTURE(ranks);
  CAPTURE(take);
  const auto parts = splitKeys(keys, ranks, split_seed);
  const auto [cutoff, calls] = runCutoff(parts, take);
  // At most 64 rounds of the bisection.
  CHECK(calls <= 64);

  std::vector<uint64_t> sorted = keys;
  std::sort(sorted.begin(), sorted.end());
  if (take == 0) {
    CHECK_FALSE(cutoff.has_value());
    return;
  }
  REQUIRE(cutoff.has_value());
  std::vector<uint64_t> selected;
  for (const auto& part : parts)
    for (uint64_t k : part)
      if (k <= *cutoff) selected.push_back(k);
  std::sort(selected.begin(), selected.end());
  CHECK(static_cast<int64_t>(selected.size()) == take);
  CHECK(selected ==
        std::vector<uint64_t>(sorted.begin(), sorted.begin() + take));
}

}  // namespace

TEST_CASE(
    "lotteryCutoff - the selection equals the serial one at any rank count") {
  SplitMix64 g(11);
  for (int trial = 0; trial < 12; ++trial) {
    const auto keys = randomUniqueKeys(1 + bounded(g, 3000), g());
    const int64_t n = static_cast<int64_t>(keys.size());
    for (size_t ranks = 1; ranks <= 8; ++ranks) {
      const int64_t take = static_cast<int64_t>(bounded(g, n + 1));
      checkCutoff(keys, ranks, take, g());
    }
  }
}

TEST_CASE("lotteryCutoff - take of 0, 1 and everyone") {
  const auto keys = randomUniqueKeys(500, 5);
  const int64_t n = static_cast<int64_t>(keys.size());
  for (size_t ranks : {1, 3, 8})
    for (int64_t take : {int64_t{0}, int64_t{1}, n - 1, n})
      checkCutoff(keys, ranks, take, 99);
}

TEST_CASE("lotteryCutoff - ranks with no candidates, and a single candidate") {
  // Fewer keys than ranks: most ranks hold nothing.
  checkCutoff({40, 10, 30}, 8, 2, 1);
  checkCutoff({40}, 8, 1, 1);
}

TEST_CASE("lotteryCutoff - keys at both ends of the 64 bit range") {
  const std::vector<uint64_t> keys = {0, kMaxKey, kMaxKey / 2, 1, kMaxKey - 1};
  for (size_t ranks : {1, 2, 4})
    for (int64_t take = 0; take <= 5; ++take) checkCutoff(keys, ranks, take, 3);
  CHECK(runCutoff({{0}}, 1).first == 0);
  CHECK(runCutoff({{kMaxKey}}, 1).first == kMaxKey);
}

TEST_CASE("lotteryCutoff - the cutoff is the key of the last person served") {
  CHECK(runCutoff({{10, 30}, {20, 40}}, 3).first == 30);
}

TEST_CASE("lotteryCutoff - a plain worked example") {
  const std::vector<std::vector<uint64_t>> world = {{5, 12}, {9, 20}};
  CHECK(runCutoff(world, 3).first == 12);
  CHECK_FALSE(runCutoff(world, 0).first.has_value());
}

TEST_CASE("lotteryCutoff - duplicate keys at the cut throw") {
  // The cut falls between the two 5s, so no cutoff serves exactly one person.
  CHECK_THROWS_WITH_AS(runCutoff({{5, 5, 7}}, 1),
                       doctest::Contains("not unique"), std::runtime_error);
  // The same when the duplicate is split across two parts.
  CHECK_THROWS_WITH_AS(runCutoff({{5}, {5, 7}}, 1),
                       doctest::Contains("not unique"), std::runtime_error);
  // The cut at the very top of the key range.
  CHECK_THROWS_WITH_AS(runCutoff({{kMaxKey, kMaxKey}}, 1),
                       doctest::Contains("not unique"), std::runtime_error);
  // Duplicates wholly inside or wholly outside the cut do not matter.
  CHECK(runCutoff({{5, 5, 7}}, 2).first == 5);
  CHECK(runCutoff({{5, 5, 7}}, 3).first == 7);
  CHECK(runCutoff({{3, 5, 5}}, 1).first == 3);
}

TEST_CASE("lotteryCutoff - bad arguments throw") {
  CHECK_THROWS_WITH_AS(runCutoff({{1, 2, 3}}, -1),
                       doctest::Contains("cannot take"), std::runtime_error);
  CHECK_THROWS_WITH_AS(runCutoff({{1, 2, 3}}, 4),
                       doctest::Contains("cannot take"), std::runtime_error);
  CHECK_THROWS_WITH_AS(runCutoff({}, 1), doctest::Contains("cannot take"),
                       std::runtime_error);
  CHECK_FALSE(runCutoff({}, 0).first.has_value());
}
