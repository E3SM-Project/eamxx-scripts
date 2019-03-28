#ifndef LI_VANILLA_HPP
#define LI_VANILLA_HPP

#include "types.hpp"
#include "util.hpp"

namespace li {

template <typename Scalar>
struct LiVanilla
{
  static constexpr const char* NAME = "vanilla";

  LiVanilla(int ncol, int km1, int km2, Scalar minthresh) :
    m_ncol(ncol),
    m_km1(km1),
    m_km2(km2),
    m_minthresh(minthresh)
  {}

  // Linearly interpolate y(x1) onto coordinates x2
  void lin_interp(const vector_2d_t<Scalar>& x1, const vector_2d_t<Scalar>& x2, const vector_2d_t<Scalar>& y1,
                  vector_2d_t<Scalar>& y2) const
  {
    for (int i = 0; i < m_ncol; ++i) {
      for (int k2 = 0; k2 < m_km2; ++k2) {
        if( x2[i][k2] <= x1[i][0] ) {
          y2[i][k2] = y1[i][0] + (y1[i][1]-y1[i][0])*(x2[i][k2]-x1[i][0])/(x1[i][1]-x1[i][0]);
        }
        else if( x2[i][k2] >= x1[i][m_km1-1] ) {
          y2[i][k2] = y1[i][m_km1-1] + (y1[i][m_km1-1]-y1[i][m_km1-2])*(x2[i][k2]-x1[i][m_km1-1])/(x1[i][m_km1-1]-x1[i][m_km1-2]);
        }
        else {
          for (int k1 = 1; k1 < m_km1; ++k1) {
            if( (x2[i][k2]>=x1[i][k1-1]) && (x2[i][k2]<x1[i][k1]) ) {
              y2[i][k2] = y1[i][k1-1] + (y1[i][k1]-y1[i][k1-1])*(x2[i][k2]-x1[i][k1-1])/(x1[i][k1]-x1[i][k1-1]);
            }
          }
        }

        if (y2[i][k2] < m_minthresh) {
          y2[i][k2] = m_minthresh;
        }
      }
    }
  }

 private:
  int m_ncol;
  int m_km1;
  int m_km2;
  Scalar m_minthresh;
};

} //namespace li

#endif
