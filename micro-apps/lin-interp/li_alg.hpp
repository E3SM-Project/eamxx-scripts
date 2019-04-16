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
    m_indx_map(ncol, std::vector<int>(km2))
#ifndef NDEBUG
    , m_indx_map_dbg(ncol, std::vector<int>(km2))
#endif
  {}

  void setup(const std::vector<Scalar>& x1, const std::vector<Scalar>& x2, int i)
  {
    setup_nlogn(x1, x2, i);
#ifndef NDEBUG
    setup_n2(x1, x2, i);
#endif
  }

  // Linearly interpolate y(x1) onto coordinates x2
  void lin_interp(const std::vector<Scalar>& x1, const std::vector<Scalar>& x2, const std::vector<Scalar>& y1,
                  std::vector<Scalar>& y2, int i)
  {
    for (int k2 = 0; k2 < m_km2; ++k2) {
      micro_assert_msg(m_indx_map[i][k2] == m_indx_map_dbg[i][k2],
                       "For k2=" << k2 << ", " << m_indx_map[i][k2] << " != " << m_indx_map_dbg[i][k2] << std::endl);
      const int k1 = m_indx_map[i][k2];
      if (k1+1 == m_km1) {
        y2[k2] = y1[k1] + (y1[k1]-y1[k1-1])*(x2[k2]-x1[k1])/(x1[k1]-x1[k1-1]);
      }
      else {
        y2[k2] = y1[k1] + (y1[k1+1]-y1[k1])*(x2[k2]-x1[k1])/(x1[k1+1]-x1[k1]);
      }

      if (y2[k2] < m_minthresh) {
        y2[k2] = m_minthresh;
      }
    }
  }

 private:

#ifndef NDEBUG
  void setup_n2(const std::vector<Scalar>& x1, const std::vector<Scalar>& x2, int i)
  {
    for (int k2 = 0; k2 < m_km2; ++k2) {
      if( x2[k2] <= x1[0] ) { // x2[k2] comes before x1[0]
        m_indx_map_dbg[i][k2] = 0;
      }
      else if( x2[k2] >= x1[m_km1-1] ) { // x2[k2] comes after x1[-1]
        m_indx_map_dbg[i][k2] = m_km1-1;
      }
      else {
        for (int k1 = 1; k1 < m_km1; ++k1) { // scan over x1
          if( (x2[k2]>=x1[k1-1]) && (x2[k2]<x1[k1]) ) { // check if x2[k2] lies within x1[k1-1] and x1[k1]
            m_indx_map_dbg[i][k2] = k1-1;
          }
        }
      }
    }
  }
#endif

  void setup_nlogn(const std::vector<Scalar>& x1, const std::vector<Scalar>& x2, int i)
  {
    auto begin = x1.begin();
    auto current = begin;
    auto upper = x1.end();
    for (int k = 0; k < m_km2; ++k) {
      const Scalar x1_indv = x2[k];
      auto ub = std::upper_bound(current, upper, x1_indv);
      int x1_idx = ub - begin;
      if (x1_idx > 0) {
        --x1_idx;
      }
      m_indx_map[i][k] = x1_idx;
      current = ub;
    }
  }

  int m_ncol;
  int m_km1;
  int m_km2;
  Scalar m_minthresh;
  vector_2d_t<int> m_indx_map; // [i][x2-idx] -> x1-idx
#ifndef NDEBUG
  vector_2d_t<int> m_indx_map_dbg; // [i][x2-idx] -> x1-idx
#endif
};

} //namespace li

#endif
