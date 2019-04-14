#ifndef LI_KOKKOS_HPP
#define LI_KOKKOS_HPP

#include "types.hpp"
#include "util.hpp"
#include "scream_assert.hpp"
#include "kokkos_util.hpp"

#include <algorithm>

namespace li {

template <typename ScalarT, typename DeviceT=DefaultDevice>
struct LiKokkos
{
  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  using KT = KokkosTypes<Device>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  using ExeSpace    = typename KT::ExeSpace;
  using MemberType  = typename KT::MemberType;
  using TeamPolicy  = typename KT::TeamPolicy;

  using Pack = Scalar;

  //
  // ------ public API -------
  //

  static constexpr const char* NAME = "kokkos";

  LiKokkos(int ncol, int km1, int km2, Scalar minthresh) :
    m_ncol(ncol),
    m_km1(km1),
    m_km2(km2),
    m_minthresh(minthresh),
    m_init(false),
    m_policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, km2)),
    m_indx_map("m_indx_map", ncol, km2)
#ifndef NDEBUG
    , m_indx_map_dbg("m_indx_map_dbg", ncol, km2)
#endif
  {}

  int km1_pack() const { return m_km1; }
  int km2_pack() const { return m_km2; }

  void setup(const view_2d<const Scalar>& x1, const view_2d<const Scalar>& x2)
  {
    if (!m_init) {
      setup_nlogn(*this, x1, x2);
#ifndef NDEBUG
      setup_n2(*this, x1, x2);
#endif
    }

    m_init = true;
  }

  KOKKOS_INLINE_FUNCTION
  void lin_interp(const view_1d<const Scalar>& x1, const view_1d<const Scalar>& x2, const view_1d<const Scalar>& y1,
                  const view_1d<Scalar>& y2, const MemberType& team) const
  {
    lin_interp_impl(*this, x1, x2, y1, y2, team);
  }

  // Linearly interpolate y(x1) onto coordinates x2
  KOKKOS_INLINE_FUNCTION
  static void lin_interp_impl(
    const LiKokkos& lik,
    const view_1d<const Scalar>& x1, const view_1d<const Scalar>& x2, const view_1d<const Scalar>& y1,
    const view_1d<Scalar>& y2, const MemberType& team)
  {
    micro_kassert_msg(lik.m_init, "Not set up");

    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, lik.m_km2), [&] (Int k2) {
      micro_kassert(lik.m_indx_map(i, k2) == lik.m_indx_map_dbg(i, k2));
      const int k1 = lik.m_indx_map(i, k2);
      if (k1+1 == lik.m_km1) {
        y2(k2) = y1(k1) + (y1(k1)-y1(k1-1))*(x2(k2)-x1(k1))/(x1(k1)-x1(k1-1));
      }
      else {
        y2(k2) = y1(k1) + (y1(k1+1)-y1(k1))*(x2(k2)-x1(k1))/(x1(k1+1)-x1(k1));
      }

      if (y2(k2) < lik.m_minthresh) {
        y2(k2) = lik.m_minthresh;
      }
    });
  }

#ifndef NDEBUG
  static void setup_n2(LiKokkos& lik, const view_2d<const Scalar>& x1, const view_2d<const Scalar>& x2)
  {
    Kokkos::parallel_for("setup_n2", lik.m_policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, lik.m_km2), [&] (Int k2) {
        if( x2(i, k2) <= x1(i, 0) ) { // x2[k2] comes before x1[0]
          lik.m_indx_map_dbg(i, k2) = 0;
        }
        else if( x2(i, k2) >= x1(i, lik.m_km1-1) ) { // x2[k2] comes after x1[-1]
          lik.m_indx_map_dbg(i, k2) = lik.m_km1-1;
        }
        else {
          for (int k1 = 1; k1 < lik.m_km1; ++k1) { // scan over x1
            if( (x2(i, k2)>=x1(i, k1-1)) && (x2(i, k2)<x1(i, k1)) ) { // check if x2[k2] lies within x1[k1-1] and x1[k1]
              lik.m_indx_map_dbg(i, k2) = k1-1;
            }
          }
        }
      });
    });
  }
#endif

  static void setup_nlogn(LiKokkos& lik, const view_2d<const Scalar>& x1, const view_2d<const Scalar>& x2)
  {
    Kokkos::parallel_for("setup_nlogn", lik.m_policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, lik.m_km2), [&] (Int k2) {
        const Scalar x1_indv = x2(i, k2);
        auto begin = &x1(i, 0);
        auto upper = begin + lik.m_km1;

        auto ub = util::upper_bound(begin, upper, x1_indv);
        int x1_idx = ub - begin;
        if (x1_idx > 0) {
          --x1_idx;
        }
        lik.m_indx_map(i, k2) = x1_idx;
      });
    });
  }

  int m_ncol;
  int m_km1;
  int m_km2;
  Scalar m_minthresh;
  bool m_init;
  TeamPolicy m_policy;
  view_2d<int> m_indx_map; // [x2-idx] -> x1-idx
#ifndef NDEBUG
  view_2d<int> m_indx_map_dbg; // [x2-idx] -> x1-idx
#endif
};

} //namespace li

#endif
