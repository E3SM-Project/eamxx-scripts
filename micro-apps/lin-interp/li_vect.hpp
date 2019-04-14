#ifndef LI_VECT_HPP
#define LI_VECT_HPP

#include "types.hpp"
#include "util.hpp"
#include "scream_assert.hpp"
#include "kokkos_util.hpp"
#include "scream_pack.hpp"

#include <algorithm>

namespace li {

using scream::pack::IntPack;
using scream::pack::IntSmallPack;
using scream::pack::smallize;
using scream::pack::scalarize;
using scream::pack::BigPack;
using scream::pack::SmallPack;

template <typename ScalarT, typename DeviceT=DefaultDevice>
struct LiVect
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

  using Pack = BigPack<Scalar>;

  //
  // ------ public API -------
  //

  static constexpr const char* NAME = "vect";

  LiVect(int ncol, int km1, int km2, Scalar minthresh) :
    m_ncol(ncol),
    m_km1(km1),
    m_km2(km2),
    m_km1_pack(scream::pack::npack<Pack>(km1)),
    m_km2_pack(scream::pack::npack<Pack>(km2)),
    m_minthresh(minthresh),
    m_init(false),
    m_policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, m_km2_pack)),
    m_indx_map("m_indx_map", ncol, scream::pack::npack<IntPack>(km2))
#ifndef NDEBUG
    , m_indx_map_dbg("m_indx_map_dbg", ncol, scream::pack::npack<IntPack>(km2))
#endif
  {}

  int km1_pack() const { return m_km1_pack; }
  int km2_pack() const { return m_km2_pack; }

  void setup(const view_2d<const Pack>& x1, const view_2d<const Pack>& x2)
  {
    if (!m_init) {
      setup_nlogn(*this, x1, x2);
#ifndef NDEBUG
      setup_n2(*this, x1, x2);
#endif
    }

    m_init = true;
  }

  // Linearly interpolate y(x1) onto coordinates x2
  KOKKOS_INLINE_FUNCTION
  void lin_interp(const view_1d<const Pack>& x1, const view_1d<const Pack>& x2, const view_1d<const Pack>& y1,
                  const view_1d<Pack>& y2, const MemberType& team) const
  {
    lin_interp_impl(*this, x1, x2, y1, y2, team);
  }

  KOKKOS_INLINE_FUNCTION
  static void lin_interp_impl(const LiVect& liv,
                              const view_1d<const Pack>& x1, const view_1d<const Pack>& x2, const view_1d<const Pack>& y1,
                              const view_1d<Pack>& y2, const MemberType& team)
  {
    micro_kassert_msg(liv.m_init, "Not set up");

    auto x1s = scalarize(x1);
    auto y1s = scalarize(y1);

    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, liv.m_km2_pack), [&] (Int k2) {
      const auto indx_pk = liv.m_indx_map(i, k2);
      micro_kassert((indx_pk == liv.m_indx_map_dbg(i, k2)).all());
      const auto end_mask = indx_pk == liv.m_km1 - 1;
      const auto not_end = !end_mask;
      if (end_mask.any()) {
        scream_masked_loop(end_mask, s) {
          int k1 = indx_pk[s];
          y2(k2)[s] = y1s(k1) + (y1s(k1)-y1s(k1-1))*(x2(k2)[s]-x1s(k1))/(x1s(k1)-x1s(k1-1));
        }
        scream_masked_loop(not_end, s) {
          int k1 = indx_pk[s];
          y2(k2)[s] = y1s(k1) + (y1s(k1+1)-y1s(k1))*(x2(k2)[s]-x1s(k1))/(x1s(k1+1)-x1s(k1));
        }
      }
      else {
        for (int s = 0; s < indx_pk.n; ++s) {
          int k1 = indx_pk[s];
          y2(k2)[s] = y1s(k1) + (y1s(k1+1)-y1s(k1))*(x2(k2)[s]-x1s(k1))/(x1s(k1+1)-x1s(k1));
        }
      }

      y2(k2).set(y2(k2) < liv.m_minthresh, liv.m_minthresh);
    });
  }

#ifndef NDEBUG
  static void setup_n2(LiVect& liv, const view_2d<const Pack>& x1, const view_2d<const Pack>& x2)
  {
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(liv.m_ncol, liv.m_km2));

    auto x1s = scalarize(x1);
    auto x2s = scalarize(x2);
    auto idxs = scalarize(liv.m_indx_map_dbg);

    Kokkos::parallel_for("setup_n2", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, liv.m_km2), [&] (Int k2) {

        if( x2s(i, k2) <= x1s(i, 0) ) { // x2[k2] comes before x1[0]
          idxs(i, k2) = 0;
        }
        else if( x2s(i, k2) >= x1s(i, liv.m_km1-1) ) { // x2[k2] comes after x1[-1]
          idxs(i, k2) = liv.m_km1-1;
        }
        else {
          for (int k1 = 1; k1 < liv.m_km1; ++k1) { // scan over x1
            if( (x2s(i, k2)>=x1s(i, k1-1)) && (x2s(i, k2)<x1s(i, k1)) ) { // check if x2[k2] lies within x1[k1-1] and x1[k1]
              idxs(i, k2) = k1-1;
            }
          }
        }
      });
    });
  }
#endif

  static void setup_nlogn(LiVect& liv, const view_2d<const Pack>& x1, const view_2d<const Pack>& x2)
  {
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(liv.m_ncol, liv.m_km2));
    auto x1s = scalarize(x1);
    auto x2s = scalarize(x2);
    auto idxs = scalarize(liv.m_indx_map);

    Kokkos::parallel_for("setup_nlogn", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, liv.m_km2), [&] (Int k2) {

        const Scalar x1_indv = x2s(i, k2);
        auto begin = &x1s(i, 0);
        auto upper = begin + liv.m_km1;

        auto ub = util::upper_bound(begin, upper, x1_indv);
        int x1_idx = ub - begin;
        if (x1_idx > 0) {
          --x1_idx;
        }
        idxs(i, k2) = x1_idx;
      });
    });
  }

  int m_ncol;
  int m_km1;
  int m_km2;
  int m_km1_pack;
  int m_km2_pack;
  Scalar m_minthresh;
  bool m_init;
  TeamPolicy m_policy;
  view_2d<IntPack> m_indx_map; // [x2-idx] -> x1-idx
#ifndef NDEBUG
  view_2d<IntPack> m_indx_map_dbg; // [x2-idx] -> x1-idx
#endif
};

} //namespace li

#endif
