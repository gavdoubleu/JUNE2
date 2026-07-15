#pragma once

#ifdef USE_MPI

#include <cstring>
#include <tuple>
#include <type_traits>

// Field-level pack/unpack used by both translation units of
// DomainCommunicator (the visitor-phase exchanges in
// src/parallel/domain_communicator.cpp and the post-interaction-phase
// exchanges in src/parallel/domain_communicator_encounters.cpp). The
// helpers copy sizeof(T) bytes and advance the cursor; sender and
// receiver must agree on field order and types. The wire size of a
// record is the sum of sizeof(field) across its packField calls.
namespace june::domain_comm_detail {

template <typename T>
inline char* packField(char* ptr, const T& v) {
  static_assert(std::is_trivially_copyable_v<T>,
                "packField requires a trivially copyable type");
  std::memcpy(ptr, &v, sizeof(T));
  return ptr + sizeof(T);
}

template <typename T>
inline const char* unpackField(const char* ptr, T& v) {
  static_assert(std::is_trivially_copyable_v<T>,
                "unpackField requires a trivially copyable type");
  std::memcpy(&v, ptr, sizeof(T));
  return ptr + sizeof(T);
}

// A record type's wire field list written once, deriving size, pack, and
// unpack from it. Only covers plain packField-able members; fields needing
// custom transcoding (e.g. an enum narrowed to a fixed-width wire type) or
// variable-length tails stay outside and are composed around pack()/unpack().
template <typename T, typename... Members>
class WireRecord {
 public:
  constexpr explicit WireRecord(Members T::*... members)
      : members_(members...) {}

  static constexpr int size() { return (static_cast<int>(sizeof(Members)) + ... + 0); }

  char* pack(char* ptr, const T& obj) const {
    return packMembers(ptr, obj, std::index_sequence_for<Members...>{});
  }

  const char* unpack(const char* ptr, T& obj) const {
    return unpackMembers(ptr, obj, std::index_sequence_for<Members...>{});
  }

 private:
  template <std::size_t... I>
  char* packMembers(char* ptr, const T& obj, std::index_sequence<I...>) const {
    ((ptr = packField(ptr, obj.*std::get<I>(members_))), ...);
    return ptr;
  }

  template <std::size_t... I>
  const char* unpackMembers(const char* ptr, T& obj,
                             std::index_sequence<I...>) const {
    ((ptr = unpackField(ptr, obj.*std::get<I>(members_))), ...);
    return ptr;
  }

  std::tuple<Members T::*...> members_;
};

template <typename T, typename... Members>
constexpr auto makeWireRecord(Members T::*... members) {
  return WireRecord<T, Members...>(members...);
}

}  // namespace june::domain_comm_detail

#endif  // USE_MPI
