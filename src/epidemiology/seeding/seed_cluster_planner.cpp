#include "epidemiology/seeding/seed_cluster_planner.h"

#include <algorithm>
#include <map>

namespace june {

namespace {

// One household as the pooled offers describe it: its candidate members, each
// with the budgets it matches. A household split across ranks arrives as
// several groups of offers sharing a key and merges here into one.
struct Household {
  uint64_t key = 0;
  std::map<PersonId, std::vector<uint32_t>> members;
  size_t matchedMembers() const {
    size_t matched = 0;
    for (const auto& [person_id, budgets] : members) {
      (void)person_id;
      if (!budgets.empty()) ++matched;
    }
    return matched;
  }
};

// Denser households first, by the seed's own score, matched / sqrt(size).
// Compared as the integer identity matched_a^2 * size_b > matched_b^2 * size_a
// so the order cannot turn on two ranks' libm agreeing bit for bit.
bool denserFirst(const Household& a, const Household& b) {
  const uint64_t matched_a = a.matchedMembers();
  const uint64_t matched_b = b.matchedMembers();
  const uint64_t left = matched_a * matched_a * b.members.size();
  const uint64_t right = matched_b * matched_b * a.members.size();
  if (left != right) return left > right;
  return a.key < b.key;
}

}  // namespace

ClusterPlan planClusteredSeed(
    const std::vector<std::vector<SeedOffer>>& offer_lists,
    const std::vector<int>& budgets) {
  std::map<uint64_t, Household> households_by_key;
  for (const auto& offers : offer_lists) {
    for (const auto& offer : offers) {
      Household& household = households_by_key[offer.key];
      household.key = offer.key;
      auto& matched_budgets = household.members[offer.person_id];
      if (offer.budget_slot > 0) {
        matched_budgets.push_back(offer.budget_slot - 1);
      }
    }
  }

  std::vector<Household> households;
  for (auto& [key, household] : households_by_key) {
    (void)key;
    if (household.matchedMembers() > 0) households.push_back(household);
  }
  std::sort(households.begin(), households.end(), denserFirst);

  int total_budget = 0;
  for (int budget : budgets) total_budget += budget;

  ClusterPlan plan;
  plan.filled_per_budget.assign(budgets.size(), 0);
  int filled = 0;
  for (const auto& household : households) {
    if (filled >= total_budget) break;
    for (const auto& [person_id, matched_budgets] : household.members) {
      // A person takes the first budget still open to them, in the order their
      // offers arrived, so a household whose members all match the same
      // exhausted budget yields nothing.
      for (uint32_t budget_index : matched_budgets) {
        if (plan.filled_per_budget[budget_index] >= budgets[budget_index]) {
          continue;
        }
        plan.assignments.push_back({person_id, budget_index});
        ++plan.filled_per_budget[budget_index];
        ++filled;
        break;
      }
    }
  }
  return plan;
}

}  // namespace june
