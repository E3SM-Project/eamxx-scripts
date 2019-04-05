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
    m_indx_map("m_indx_map", scream::pack::npack<IntPack>(km2))
#ifndef NDEBUG
    , m_indx_map_dbg("m_indx_map_dbg", scream::pack::npack<IntPack>(km2))
#endif
  {}

  int km1_pack() const { return m_km1_pack; }
  int km2_pack() const { return m_km2_pack; }

  // Linearly interpolate y(x1) onto coordinates x2
  void lin_interp(const view_2d<const Pack>& x1, const view_2d<const Pack>& x2, const view_2d<const Pack>& y1,
                  const view_2d<Pack>& y2)
  {
    if (!m_init) {
      // assumes all cols have the same coords
      setup_nlogn(util::subview(x1, 0), util::subview(x2, 0));
#ifndef NDEBUG
      setup_n2(util::subview(x1, 0), util::subview(x2, 0));
#endif
    }

    Kokkos::parallel_for("lin_interp", m_policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();
      auto x1s = scalarize(util::subview(x1, i));
      auto x2s = util::subview(x2, i);
      auto y1s = scalarize(util::subview(y1, i));
      auto y2s = util::subview(y2, i);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, m_km2_pack), [&] (Int k2) {
        const auto indx_pk = m_indx_map(k2);
        micro_kassert((indx_pk == m_indx_map_dbg(k2)).all());
        const auto end_mask = indx_pk == m_km1 - 1;
        const auto not_end = !end_mask;
        if (end_mask.any()) {
          scream_masked_loop(end_mask, s) {
            int k1 = indx_pk[s];
            y2s(k2)[s] = y1s(k1) + (y1s(k1)-y1s(k1-1))*(x2s(k2)[s]-x1s(k1))/(x1s(k1)-x1s(k1-1));
          }
          scream_masked_loop(not_end, s) {
            int k1 = indx_pk[s];
            y2s(k2)[s] = y1s(k1) + (y1s(k1+1)-y1s(k1))*(x2s(k2)[s]-x1s(k1))/(x1s(k1+1)-x1s(k1));
          }
        }
        else {
          for (int s = 0; s < indx_pk.n; ++s) {
            int k1 = indx_pk[s];
            y2s(k2)[s] = y1s(k1) + (y1s(k1+1)-y1s(k1))*(x2s(k2)[s]-x1s(k1))/(x1s(k1+1)-x1s(k1));
          }
        }

        y2s(k2).set(y2s(k2) < m_minthresh, m_minthresh);
      });
    });
  }

 private:

  void setup_n2(const view_1d<const Pack>& x1, const view_1d<const Pack>& x2)
  {
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_km2, m_km1-1));

    auto x1s = scalarize(x1);
    auto x2s = scalarize(x2);
    auto idxs = scalarize(m_indx_map_dbg);

    Kokkos::parallel_for("setup_n2", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int k2 = team.league_rank();
      if( x2s(k2) <= x1s(0) ) { // x2[k2] comes before x1[0]
        idxs(k2) = 0;
      }
      else if( x2s(k2) >= x1s(m_km1-1) ) { // x2[k2] comes after x1[-1]
        idxs(k2) = m_km1-1;
      }
      else {
        int k1_max = 0;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, m_km1), [&] (int k1_arg, int& k1_max) {
          const int k1 = k1_arg + 1;
          if( (x2s(k2)>=x1s(k1-1)) && (x2s(k2)<x1s(k1)) ) { // check if x2[k2] lies within x1[k1-1] and x1[k1]
            k1_max = k1-1;
          }
        }, Kokkos::Max<int>(k1_max));
        idxs(k2) = k1_max;
      }
    });
  }

  void setup_nlogn(const view_1d<const Pack>& x1, const view_1d<const Pack>& x2)
  {
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_km2, 1));
    auto x1s = scalarize(x1);
    auto x2s = scalarize(x2);
    auto idxs = scalarize(m_indx_map);

    Kokkos::parallel_for("setup_nlogn", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int k2 = team.league_rank();
      const Scalar x1_indv = x2s(k2);
      auto begin = &x1s(0);
      auto upper = begin + m_km1;

      auto ub = std::upper_bound(begin, upper, x1_indv);
      int x1_idx = ub - begin;
      if (x1_idx > 0) {
        --x1_idx;
      }
      idxs(k2) = x1_idx;
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
  view_1d<IntPack> m_indx_map; // [x2-idx] -> x1-idx
#ifndef NDEBUG
  view_1d<IntPack> m_indx_map_dbg; // [x2-idx] -> x1-idx
#endif
};

} //namespace li

#endif
