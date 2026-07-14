// Unit tests for WireRecord (include/parallel/domain_communicator_detail.h):
// the shared seam that derives wire size, pack, and unpack from a single
// field list per record type. No MPI runtime involved (pure memcpy), so this
// runs as a plain ctest binary, not under mpirun.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#ifdef USE_MPI
#include <cstdint>
#include <vector>

#include "parallel/domain_communicator_detail.h"

using june::domain_comm_detail::makeWireRecord;
using june::domain_comm_detail::packField;
using june::domain_comm_detail::unpackField;

namespace {

struct DummyRecord {
  int32_t id;
  uint8_t small_flag;
  double value;
};

}  // namespace

TEST_CASE("WireRecord: size() sums the listed members' sizeof") {
  constexpr auto kDummyWire = makeWireRecord(
      &DummyRecord::id, &DummyRecord::small_flag, &DummyRecord::value);
  CHECK(kDummyWire.size() ==
        static_cast<int>(sizeof(int32_t) + sizeof(uint8_t) + sizeof(double)));
}

TEST_CASE("WireRecord: pack/unpack round-trips the listed fields") {
  constexpr auto kDummyWire = makeWireRecord(
      &DummyRecord::id, &DummyRecord::small_flag, &DummyRecord::value);

  DummyRecord original{42, 7, 3.14159};
  std::vector<char> buf(kDummyWire.size());

  char* end = kDummyWire.pack(buf.data(), original);
  CHECK(end == buf.data() + buf.size());

  DummyRecord restored{};
  const char* read_end = kDummyWire.unpack(buf.data(), restored);
  CHECK(read_end == buf.data() + buf.size());

  CHECK(restored.id == original.id);
  CHECK(restored.small_flag == original.small_flag);
  CHECK(restored.value == doctest::Approx(original.value));
}

TEST_CASE("WireRecord: field order in the list is the wire order") {
  // Deliberately list fields out of declaration order to confirm the wire
  // layout follows the list, not the struct's own member order.
  constexpr auto kReorderedWire = makeWireRecord(
      &DummyRecord::small_flag, &DummyRecord::id, &DummyRecord::value);

  DummyRecord original{99, 3, 2.5};
  std::vector<char> buf(kReorderedWire.size());
  kReorderedWire.pack(buf.data(), original);

  // small_flag (1 byte) packed first, then id (4 bytes).
  uint8_t first_byte;
  unpackField(buf.data(), first_byte);
  CHECK(first_byte == original.small_flag);

  int32_t second_field;
  unpackField(buf.data() + sizeof(uint8_t), second_field);
  CHECK(second_field == original.id);
}

TEST_CASE(
    "WireRecord: composes with a manually-appended special field, "
    "mirroring the EncounterReply status_byte pattern") {
  // WireRecord only ever covers plain, directly packField-able members.
  // Fields needing custom transcoding (e.g. an enum deliberately narrowed to
  // a fixed-width wire type) stay outside WireRecord and are packed/unpacked
  // as an explicit extra step around it.
  constexpr auto kDummyWire =
      makeWireRecord(&DummyRecord::id, &DummyRecord::value);

  enum class Status : uint32_t { OK = 0, FAILED = 1 };

  DummyRecord original{5, 0, 1.0};
  Status status = Status::FAILED;

  const int total_size = kDummyWire.size() + static_cast<int>(sizeof(uint8_t));
  std::vector<char> buf(total_size);

  char* ptr = kDummyWire.pack(buf.data(), original);
  uint8_t status_byte = static_cast<uint8_t>(status);
  ptr = packField(ptr, status_byte);
  CHECK(ptr == buf.data() + buf.size());

  DummyRecord restored{};
  const char* read_ptr = kDummyWire.unpack(buf.data(), restored);
  uint8_t restored_status_byte;
  read_ptr = unpackField(read_ptr, restored_status_byte);
  Status restored_status = static_cast<Status>(restored_status_byte);

  CHECK(read_ptr == buf.data() + buf.size());
  CHECK(restored.id == original.id);
  CHECK(restored.value == doctest::Approx(original.value));
  CHECK(restored_status == status);
}

#else

TEST_CASE("WireRecord tests skipped: USE_MPI not enabled") {}

#endif  // USE_MPI
