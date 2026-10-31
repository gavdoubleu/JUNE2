#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#ifdef USE_MPI
#include <mpi.h>

#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "epidemiology/vaccination_uptake.h"
#include "utils/mpi_utils.h"
#include "vaccination_uptake_fixtures.h"

using namespace june;
using namespace june::uptake_fixture;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  const int result = context.run();
  MPI_Finalize();
  return result;
}

namespace {

constexpr uint64_t kSeed = 20261028;
constexpr uint64_t kCampaignIndex = 2;

using Flags = std::map<PersonId, std::pair<uint64_t, uint64_t>>;

Flags flagsById(const WorldState& world, const CampaignUptake& result) {
  Flags flags;
  for (size_t i = 0; i < world.people.size(); ++i)
    flags[world.people[i].id] = {result.selected[i], result.accepted[i]};
  return flags;
}

VaccinationCampaignConfig childrenDecidedByParents() {
  VaccinationTier tier = children();
  shareResidence(tier, DeciderCombine::Mean);
  tier.uptake.rows = {row({criterion("age", "<", 35.0)}, 0.9), row({}, 0.3)};
  return makeCampaign({tier});
}

VaccinationCampaignConfig threeTiers() {
  VaccinationTier kids = children();
  shareResidence(kids, DeciderCombine::Mean);
  kids.uptake.rows = {row({criterion("age", "<", 33.0)}, 0.8), row({}, 0.4)};
  VaccinationTier grown = adults();
  grown.uptake.rows = {row({criterion("age", "<", 45.0)}, 0.6), row({}, 0.2)};
  VaccinationTier everyone = makeTier("everyone", {});
  shareResidence(everyone, DeciderCombine::Max);
  everyone.uptake.rows = {row({}, 0.5)};
  return makeCampaign({kids, grown, everyone});
}

int rankOf() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int sizeOf() {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

// This rank's homes only, as in the real distributed load
WorldState myWorld(const std::vector<Home>& homes) {
  const size_t rank = static_cast<size_t>(rankOf());
  const size_t size = static_cast<size_t>(sizeOf());
  return makeWorld(homes, [&](size_t h) { return h % size == rank; });
}

}  // namespace

// Every rank holds its own homes. The counts after the reduce and each
// person's bits must be what one process finds on the whole world.
TEST_CASE(
    "settleUptake - ranks holding their own homes give the serial result") {
  const auto homes = mixedHomes(1500);
  const WorldState whole = makeWorld(homes);
  const auto expected =
      settleUptakeLocal(whole, threeTiers(), kCampaignIndex, kSeed);
  const Flags expected_flags = flagsById(whole, expected);

  const WorldState mine = myWorld(homes);
  const auto got = settleUptake(mine, threeTiers(), kCampaignIndex, kSeed);

  CHECK(mpi_runtime::globalSum(static_cast<int64_t>(mine.people.size())) ==
        static_cast<int64_t>(whole.people.size()));
  REQUIRE(got.tiers.size() == expected.tiers.size());
  for (size_t t = 0; t < got.tiers.size(); ++t) {
    CAPTURE(t);
    CHECK(got.tiers[t].selected > 0);
    CHECK(got.tiers[t].selected == expected.tiers[t].selected);
    CHECK(got.tiers[t].accepted == expected.tiers[t].accepted);
    CHECK(got.tiers[t].units == expected.tiers[t].units);
    CHECK(got.tiers[t].unmatched.checked ==
          expected.tiers[t].unmatched.checked);
    CHECK(got.tiers[t].venue_count.checked ==
          expected.tiers[t].venue_count.checked);
    CHECK(got.tiers[t].split_unit.checked ==
          expected.tiers[t].split_unit.checked);
    CHECK(got.tiers[t].no_decider.checked ==
          expected.tiers[t].no_decider.checked);
  }
  const Flags mine_flags = flagsById(mine, got);
  CHECK_FALSE(mine_flags.empty());
  for (const auto& [id, bits] : mine_flags)
    CHECK(bits == expected_flags.at(id));
}

// A failure on any rank makes every rank throw, with the text one process
// gives for the whole world.
TEST_CASE("settleUptake - every rank throws the serial error text") {
  Home none{{30, 40, 5}};
  none.registered = false;
  Home elsewhere{{30, 40, 5}};
  elsewhere.person_geo_unit = 999;

  struct Case {
    std::string what;
    std::vector<Home> homes;
    VaccinationCampaignConfig campaign;
  };
  // A scope on the venue type would leave out anyone with no venue at all
  VaccinationTier open_scope = children();
  shareResidence(open_scope, DeciderCombine::Mean);
  open_scope.decision_rules.scope = {criterion("age", ">=", 0.0)};
  VaccinationTier gap = adults();
  gap.uptake.rows = {row({criterion("age", "<", 35.0)}, 0.5)};
  const std::vector<Case> cases = {
      {"no uptake row", homesWith(Home{{30, 40}}, 3, 2), makeCampaign({gap})},
      {"venue count", homesWith(none, 3, 2), makeCampaign({open_scope})},
      {"no decider", homesWith(Home{{5, 7}}, 3, 2), childrenDecidedByParents()},
      {"split unit", homesWith(elsewhere, 3, 2), childrenDecidedByParents()},
  };
  for (const Case& c : cases) {
    CAPTURE(c.what);
    std::string serial;
    try {
      uptake_detail::throwOnFailure(
          c.campaign, settleUptakeLocal(makeWorld(c.homes), c.campaign,
                                        kCampaignIndex, kSeed));
    } catch (const std::runtime_error& e) {
      serial = e.what();
    }
    REQUIRE_FALSE(serial.empty());

    const WorldState mine = myWorld(c.homes);
    std::string text;
    try {
      settleUptake(mine, c.campaign, kCampaignIndex, kSeed);
    } catch (const std::runtime_error& e) {
      text = e.what();
    }
    CHECK(text == serial);
  }
}

#else
TEST_CASE("MPI vaccination uptake requires MPI") {}
#endif
