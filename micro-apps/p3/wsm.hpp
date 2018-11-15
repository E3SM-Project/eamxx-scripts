#ifndef WSM_HPP
#define WSM_HPP

#include "types.hpp"

#include <memory>

namespace unit_test {
struct UnitWrap;
}

namespace util {

template <typename T, typename D=DefaultDevice>
class WorkspaceManager
{
  class Impl;
  class WorkspaceImpl;

  std::shared_ptr<Impl> impl;

  friend struct unit_test::UnitWrap;

 public:
  using TeamPolicy = typename KokkosTypes<D>::TeamPolicy;
  using MemberType = typename KokkosTypes<D>::MemberType;

  template <typename S>
  using view_1d = typename KokkosTypes<D>::template view_1d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KokkosTypes<D>::template view_1d_ptr_array<S, N>;

  //
  // public API
  //

  WorkspaceManager(int size, int max_used, TeamPolicy policy);

  void report() const;

  class Workspace;

  KOKKOS_INLINE_FUNCTION
  Workspace get_workspace(const MemberType& team) const;

  class Workspace {
   public:
    template <typename S=T>
    KOKKOS_INLINE_FUNCTION
    Unmanaged<view_1d<S> > take(const char* name) const;

    template <size_t N, typename S=T>
    KOKKOS_INLINE_FUNCTION
    void take_many_contiguous_unsafe(const Kokkos::Array<const char*, N>& names,
                                     const view_1d_ptr_array<S, N>& ptrs) const;

    template <size_t N, typename S=T>
    KOKKOS_INLINE_FUNCTION
    void take_many(const Kokkos::Array<const char*, N>& names,
                   const view_1d_ptr_array<S, N>& ptrs) const;

    template <size_t N, typename S=T>
    KOKKOS_INLINE_FUNCTION
    void take_many_and_reset(const Kokkos::Array<const char*, N>& names,
                             const view_1d_ptr_array<S, N>& ptrs) const;

    // Wrapper so caller doesn't have to specify scalar type.
    template <typename View>
    KOKKOS_FORCEINLINE_FUNCTION
    void release(const View& space, std::enable_if<View::rank == 1>* = 0) const;

#ifndef NDEBUG
    template <typename View>
    KOKKOS_INLINE_FUNCTION
    const char* get_name(const View& space, std::enable_if<View::rank == 1>* = 0) const;
#endif

    KOKKOS_INLINE_FUNCTION
    void reset() const;

    // Print the linked list. Obviously not a device function.
    void print() const;

   private:

    KOKKOS_INLINE_FUNCTION
    Workspace(const Impl& parent, int ws_idx, const MemberType& team);

    std::shared_ptr<WorkspaceImpl> impl;

    friend struct unit_test::UnitWrap;

    friend Impl;
  }; // class Workspace
}; // class WorkspaceManager

} // namespace util

#endif
