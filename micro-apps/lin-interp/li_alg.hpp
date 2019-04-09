#ifndef LI_ALG_HPP
#define LI_ALG_HPP

#include "types.hpp"
#include "util.hpp"
#include "scream_assert.hpp"

#include <algorithm>

namespace li {

template <typename Scalar>
struct LiAlg
{
  static constexpr const char* NAME = "alg";

  LiAlg(int ncol, int km1, int km2, Scalar minthresh) :
    m_ncol(ncol),
    m_km1(km1),
    m_km2(km2),
    m_minthresh(minthresh),
    m_init(false),
    m_indx_map(km2)
#ifndef NDEBUG
    , m_indx_map_dbg(km2)
#endif
  {}

  // Linearly interpolate y(x1) onto coordinates x2
  void lin_interp(const vector_2d_t<Scalar>& x1, const vector_2d_t<Scalar>& x2, const vector_2d_t<Scalar>& y1,
                  vector_2d_t<Scalar>& y2)
  {
    if (!m_init) {
      setup_nlogn(x1[0], x2[0]);
#ifndef NDEBUG
      setup_n2(x1[0], x2[0]);
#endif
    }

    for (int i = 0; i < m_ncol; ++i) {
      for (int k2 = 0; k2 < m_km2; ++k2) {
        micro_assert_msg(m_indx_map[k2] == m_indx_map_dbg[k2],
                         "For k2=" << k2 << ", " << m_indx_map[k2] << " != " << m_indx_map_dbg[k2] << std::endl);
        const int k1 = m_indx_map[k2];
        if (k1+1 == m_km1) {
          y2[i][k2] = y1[i][k1] + (y1[i][k1]-y1[i][k1-1])*(x2[i][k2]-x1[i][k1])/(x1[i][k1]-x1[i][k1-1]);
        }
        else {
          y2[i][k2] = y1[i][k1] + (y1[i][k1+1]-y1[i][k1])*(x2[i][k2]-x1[i][k1])/(x1[i][k1+1]-x1[i][k1]);
        }

        if (y2[i][k2] < m_minthresh) {
          y2[i][k2] = m_minthresh;
        }
      }
    }
  }

 private:

#ifndef NDEBUG
  void setup_n2(const std::vector<Scalar>& x1, const std::vector<Scalar>& x2)
  {
    for (int k2 = 0; k2 < m_km2; ++k2) {
      if( x2[k2] <= x1[0] ) { // x2[k2] comes before x1[0]
        m_indx_map_dbg[k2] = 0;
      }
      else if( x2[k2] >= x1[m_km1-1] ) { // x2[k2] comes after x1[-1]
        m_indx_map_dbg[k2] = m_km1-1;
      }
      else {
        for (int k1 = 1; k1 < m_km1; ++k1) { // scan over x1
          if( (x2[k2]>=x1[k1-1]) && (x2[k2]<x1[k1]) ) { // check if x2[k2] lies within x1[k1-1] and x1[k1]
            m_indx_map_dbg[k2] = k1-1;
          }
        }
      }
    }
  }
#endif

  void setup_nlogn(const std::vector<Scalar>& x1, const std::vector<Scalar>& x2)
  {
    auto begin = x1.begin();
    auto current = begin;
    auto upper = x1.end();
    for (int i = 0; i < m_km2; ++i) {
      const Scalar x1_indv = x2[i];
      auto ub = std::upper_bound(current, upper, x1_indv);
      int x1_idx = ub - begin;
      if (x1_idx > 0) {
        --x1_idx;
      }
      m_indx_map[i] = x1_idx;
      current = ub;
    }
  }

  int m_ncol;
  int m_km1;
  int m_km2;
  Scalar m_minthresh;
  bool m_init;
  std::vector<int> m_indx_map; // [x2-idx] -> x1-idx
#ifndef NDEBUG
  std::vector<int> m_indx_map_dbg; // [x2-idx] -> x1-idx
#endif
};

} //namespace li

#endif
