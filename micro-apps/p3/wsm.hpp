#ifndef WSM_HPP
#define WSM_HPP

#include "types.hpp"

namespace unit_test {
struct UnitWrap;
}

namespace util {

template <typename T, typename DeviceT=DefaultDevice>
class WorkspaceManager
{
 public:

  using Device = DeviceT;

  using TeamPolicy = typename KokkosTypes<Device>::TeamPolicy;
  using MemberType = typename KokkosTypes<Device>::MemberType;
  using ExeSpace   = typename KokkosTypes<Device>::ExeSpace;

  template <typename S>
  using view_1d = typename KokkosTypes<Device>::template view_1d<S>;
  template <typename S>
  using view_2d = typename KokkosTypes<Device>::template view_2d<S>;
  template <typename S>
  using view_3d = typename KokkosTypes<Device>::template view_3d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KokkosTypes<Device>::template view_1d_ptr_array<S, N>;

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
    void release(const View& space, std::enable_if<View::rank == 1>* = 0) const
    { release_impl<typename View::value_type>(space); }

#ifndef NDEBUG
    template <typename View>
    KOKKOS_INLINE_FUNCTION
    const char* get_name(const View& space, std::enable_if<View::rank == 1>* = 0) const
    { return get_name_impl<typename View::value_type>(space); }
#endif

    KOKKOS_INLINE_FUNCTION
    void reset() const;

    // Print the linked list. Obviously not a device function.
    void print() const;

    //
    // Private
    //

#ifndef KOKKOS_ENABLE_CUDA
   private: // for CUDA
#endif

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    void release_impl(const Unmanaged<view_1d<S> >& space) const;

#ifndef NDEBUG
    template <typename S>
    KOKKOS_INLINE_FUNCTION
    const char* get_name_impl(const Unmanaged<view_1d<S> >& space) const;

    KOKKOS_INLINE_FUNCTION
    void change_num_used(int change_by) const;

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    void change_indv_meta(const Unmanaged<view_1d<S> >& space, const char* name, bool release=false) const;

    KOKKOS_INLINE_FUNCTION
    int get_name_idx(const char* name, bool add) const;

    KOKKOS_INLINE_FUNCTION
    int get_alloc_count(const char* name) const
    { return m_parent.m_counts(m_ws_idx, get_name_idx(name), 0); }

    KOKKOS_INLINE_FUNCTION
    int get_release_count(const char* name) const
    { return m_parent.m_counts(m_ws_idx, get_name_idx(name), 1); }

    KOKKOS_INLINE_FUNCTION
    int get_num_used() const
    { return m_parent.m_num_used(m_ws_idx); }

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    bool is_active(const Unmanaged<view_1d<S> >& space) const
    { return m_parent.m_active(m_ws_idx, m_parent.template get_index<S>(space));}
#endif

    KOKKOS_INLINE_FUNCTION
    Workspace(const WorkspaceManager& parent, int ws_idx, const MemberType& team);

    friend struct unit_test::UnitWrap;
    friend class WorkspaceManager;

    const WorkspaceManager& m_parent;
    const MemberType& m_team;
    const int m_ws_idx;
    int& m_next_slot;
  }; // class Workspace

#ifndef KOKKOS_ENABLE_CUDA
 private:
#endif

  friend struct unit_test::UnitWrap;

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int get_index(const Unmanaged<view_1d<S> >& space) const
  { return reinterpret_cast<const int*>(reinterpret_cast<const T*>(space.data()) - m_reserve)[0]; }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int get_next(const Unmanaged<view_1d<S> >& space) const
  { return reinterpret_cast<const int*>(reinterpret_cast<const T*>(space.data()) - m_reserve)[1]; }

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  int set_next_and_get_index(const Unmanaged<view_1d<S> >& space, int next) const;

  template <typename S=T>
  KOKKOS_FORCEINLINE_FUNCTION
  Unmanaged<view_1d<S> > get_space_in_slot(const int team_idx, const int slot) const;

  KOKKOS_INLINE_FUNCTION
  void init_metadata(const int ws_idx, const int slot) const;

  static void init(const WorkspaceManager& wm, const view_2d<T>& data,
                   const int concurrent_teams, const int max_used, const int total);

  //
  // data
  //

  enum { m_pad_factor   = OnGpu<ExeSpace>::value ? 1 : 32,
         m_max_name_len = 128,
         m_max_names    = 256
  };

  util::TeamUtils<ExeSpace> m_tu;
  int m_concurrent_teams, m_reserve, m_size, m_total, m_max_used;
#ifndef NDEBUG
  view_1d<int> m_num_used;
  view_1d<int> m_high_water;
  view_2d<bool> m_active;
  view_3d<char> m_curr_names;
  view_3d<char> m_all_names;
  view_3d<int> m_counts;
#endif
  view_1d<int> m_next_slot;
  view_2d<T> m_data;
}; // class WorkspaceManager

} // namespace util

#include "wsm_impl.hpp"

#endif
