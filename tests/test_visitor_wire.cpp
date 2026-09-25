// Unit tests for the visitor wire seam (include/parallel/visitor_wire.h):
// per-record size, pack and unpack of Domain::VisitorData. Pure memcpy, no
// MPI runtime, so this runs as a plain ctest binary, not under mpirun.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#ifdef USE_MPI
#include <stdexcept>
#include <vector>

#include "parallel/visitor_wire.h"

using june::Domain;
namespace visitor_wire = june::visitor_wire;

namespace {

Domain::VisitorData makeVisitor(int num_modes, int fomite_sub_bins) {
  Domain::VisitorData visitor{};
  visitor.person_id = 4242;
  visitor.home_rank = 1;
  visitor.venue_id = 777;
  visitor.subset_idx = 3;
  visitor.is_infected = true;
  visitor.is_infectious = true;
  visitor.immunity_level = 0.25f;
  visitor.encounter_type_id = 5;
  visitor.symptom_id = 9;
  for (int mode = 0; mode < num_modes; ++mode)
    visitor.integrated_infectiousness.push_back(0.1 * (mode + 1) + 1e-17);
  for (int bin = 0; bin < fomite_sub_bins; ++bin)
    visitor.fomite_deposition_sub.push_back(1.0 / 3.0 * (bin + 1));
  return visitor;
}

void checkRoundTrip(const Domain::VisitorData& sent,
                    const visitor_wire::TailCounts& tails) {
  const int size = visitor_wire::recordSize(sent, tails);
  std::vector<char> buffer(size);
  char* packed_end = visitor_wire::pack(buffer.data(), sent, tails);
  CHECK(packed_end - buffer.data() == size);

  Domain::VisitorData received{};
  const char* unpacked_end =
      visitor_wire::unpack(buffer.data(), received, tails);
  CHECK(unpacked_end - buffer.data() == size);

  CHECK(received.person_id == sent.person_id);
  CHECK(received.home_rank == sent.home_rank);
  CHECK(received.venue_id == sent.venue_id);
  CHECK(received.subset_idx == sent.subset_idx);
  CHECK(received.is_infected == sent.is_infected);
  CHECK(received.is_infectious == sent.is_infectious);
  CHECK(received.immunity_level == sent.immunity_level);
  CHECK(received.encounter_type_id == sent.encounter_type_id);
  CHECK(received.symptom_id == sent.symptom_id);
  // Exact equality: tails must arrive bit-identical.
  CHECK(received.integrated_infectiousness == sent.integrated_infectiousness);
  CHECK(received.fomite_deposition_sub == sent.fomite_deposition_sub);
}

}  // namespace

// Uninfected, as the sender builds it: no tails.
Domain::VisitorData makeUninfectedVisitor() {
  Domain::VisitorData visitor = makeVisitor(0, 0);
  visitor.is_infected = false;
  visitor.is_infectious = false;
  visitor.symptom_id = 0;
  return visitor;
}

// Infected, not yet infectious, as the sender builds it: deposits only.
Domain::VisitorData makeIncubatingVisitor(int fomite_sub_bins) {
  Domain::VisitorData visitor = makeVisitor(0, fomite_sub_bins);
  visitor.is_infectious = false;
  return visitor;
}

TEST_CASE("visitor wire: record round-trips header and both tails") {
  const visitor_wire::TailCounts tails{2, 10};
  checkRoundTrip(makeVisitor(2, 10), tails);
}

TEST_CASE("visitor wire: uninfected record sends the header alone") {
  const visitor_wire::TailCounts tails{2, 10};
  const Domain::VisitorData uninfected = makeUninfectedVisitor();
  checkRoundTrip(uninfected, tails);
  CHECK(visitor_wire::recordSize(uninfected, tails) +
            (2 + 10) * static_cast<int>(sizeof(double)) ==
        visitor_wire::recordSize(makeVisitor(2, 10), tails));
}

TEST_CASE("visitor wire: pack throws when a tail's length differs from its "
          "count") {
  const visitor_wire::TailCounts tails{2, 10};
  std::vector<char> buffer(visitor_wire::recordSize(makeVisitor(2, 10), tails) +
                           64);

  SUBCASE("integrated_infectiousness") {
    CHECK_THROWS_AS(visitor_wire::pack(buffer.data(), makeVisitor(3, 10), tails),
                    std::runtime_error);
  }
  SUBCASE("fomite_deposition_sub") {
    CHECK_THROWS_AS(visitor_wire::pack(buffer.data(), makeVisitor(2, 9), tails),
                    std::runtime_error);
  }
}

TEST_CASE("visitor wire: zero fomite sub-bins round-trips") {
  // No fomite mode, as in config_2021.
  const visitor_wire::TailCounts tails{2, 0};
  checkRoundTrip(makeVisitor(2, 0), tails);
}

TEST_CASE("visitor wire: slice of records round-trips in order") {
  const visitor_wire::TailCounts tails{2, 10};
  std::vector<Domain::VisitorData> sent;
  for (int index = 0; index < 3; ++index) {
    sent.push_back(makeVisitor(2, 10));
    sent.back().person_id = 100 + index;
  }

  const int size = visitor_wire::sliceSize(sent, tails);
  std::vector<char> buffer(size);
  char* ptr = buffer.data();
  for (const auto& visitor : sent) ptr = visitor_wire::pack(ptr, visitor, tails);
  CHECK(ptr - buffer.data() == size);

  std::vector<Domain::VisitorData> received;
  visitor_wire::unpackSlice(
      buffer.data(), buffer.data() + buffer.size(), tails,
      [&](Domain::VisitorData&& visitor) { received.push_back(visitor); });
  REQUIRE(received.size() == sent.size());
  for (std::size_t index = 0; index < sent.size(); ++index) {
    CHECK(received[index].person_id == sent[index].person_id);
    CHECK(received[index].fomite_deposition_sub ==
          sent[index].fomite_deposition_sub);
  }
}

TEST_CASE("visitor wire: unpacking a slice throws unless records end exactly "
          "at its end") {
  const visitor_wire::TailCounts tails{2, 10};
  const std::vector<Domain::VisitorData> sent{makeVisitor(2, 10),
                                              makeVisitor(2, 10)};
  std::vector<char> buffer(visitor_wire::sliceSize(sent, tails));
  char* ptr = buffer.data();
  for (const auto& visitor : sent) ptr = visitor_wire::pack(ptr, visitor, tails);
  auto ignore = [](Domain::VisitorData&&) {};

  SUBCASE("truncated: last record's tail cut short") {
    buffer.pop_back();
    CHECK_THROWS_AS(visitor_wire::unpackSlice(buffer.data(),
                                              buffer.data() + buffer.size(),
                                              tails, ignore),
                    std::runtime_error);
  }
  SUBCASE("leftover: bytes after the last record") {
    buffer.insert(buffer.end(), 3, '\0');
    CHECK_THROWS_AS(visitor_wire::unpackSlice(buffer.data(),
                                              buffer.data() + buffer.size(),
                                              tails, ignore),
                    std::runtime_error);
  }
}

TEST_CASE("visitor wire: incubating record sends deposits, not ii") {
  const visitor_wire::TailCounts tails{2, 10};
  const Domain::VisitorData incubating = makeIncubatingVisitor(10);
  checkRoundTrip(incubating, tails);
  CHECK(visitor_wire::recordSize(incubating, tails) +
            2 * static_cast<int>(sizeof(double)) ==
        visitor_wire::recordSize(makeVisitor(2, 10), tails));
}

TEST_CASE("visitor wire: pack throws when a tail it skips is nonzero") {
  const visitor_wire::TailCounts tails{2, 10};
  std::vector<char> buffer(visitor_wire::recordSize(makeVisitor(2, 10), tails));

  SUBCASE("uninfected with a deposit") {
    Domain::VisitorData uninfected = makeUninfectedVisitor();
    uninfected.fomite_deposition_sub.assign(10, 0.0);
    uninfected.fomite_deposition_sub[4] = 1e-9;
    CHECK_THROWS_AS(visitor_wire::pack(buffer.data(), uninfected, tails),
                    std::runtime_error);
  }
  SUBCASE("not infectious with integrated infectiousness") {
    Domain::VisitorData incubating = makeIncubatingVisitor(10);
    incubating.integrated_infectiousness = {0.0, 1e-9};
    CHECK_THROWS_AS(visitor_wire::pack(buffer.data(), incubating, tails),
                    std::runtime_error);
  }
}

TEST_CASE("visitor wire: an all-zero skipped tail packs and arrives empty") {
  const visitor_wire::TailCounts tails{2, 10};
  Domain::VisitorData uninfected = makeUninfectedVisitor();
  uninfected.integrated_infectiousness.assign(2, 0.0);
  uninfected.fomite_deposition_sub.assign(10, 0.0);
  std::vector<char> buffer(visitor_wire::recordSize(uninfected, tails));
  visitor_wire::pack(buffer.data(), uninfected, tails);

  Domain::VisitorData received{};
  visitor_wire::unpack(buffer.data(), received, tails);
  CHECK(received.integrated_infectiousness.empty());
  CHECK(received.fomite_deposition_sub.empty());
}

#endif  // USE_MPI
