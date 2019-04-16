#ifndef LI_VANILLA_HPP
#define LI_VANILLA_HPP

#include "types.hpp"
#include "util.hpp"

namespace li {

template <typename Scalar>
struct LiVanilla
{
  static constexpr const char* NAME = "vanilla";

  LiVanilla(int /*ncol*/, int km1, int km2, Scalar minthresh) :
    m_km1(km1),
    m_km2(km2),
    m_minthresh(minthresh)
  {}

  void setup(const std::vector<Scalar>& /*x1*/, const std::vector<Scalar>& /*x2*/, int /*i*/)
  {}

  // Linearly interpolate y(x1) onto coordinates x2
  void lin_interp(const std::vector<Scalar>& x1, const std::vector<Scalar>& x2, const std::vector<Scalar>& y1,
                  std::vector<Scalar>& y2, int /*i*/) const
  {
    for (int k2 = 0; k2 < m_km2; ++k2) {
      if( x2[k2] <= x1[0] ) {
        y2[k2] = y1[0] + (y1[1]-y1[0])*(x2[k2]-x1[0])/(x1[1]-x1[0]);
      }
      else if( x2[k2] >= x1[m_km1-1] ) {
        y2[k2] = y1[m_km1-1] + (y1[m_km1-1]-y1[m_km1-2])*(x2[k2]-x1[m_km1-1])/(x1[m_km1-1]-x1[m_km1-2]);
      }
      else {
        for (int k1 = 1; k1 < m_km1; ++k1) {
          if( (x2[k2]>=x1[k1-1]) && (x2[k2]<x1[k1]) ) {
            y2[k2] = y1[k1-1] + (y1[k1]-y1[k1-1])*(x2[k2]-x1[k1-1])/(x1[k1]-x1[k1-1]);
          }
        }
      }

      if (y2[k2] < m_minthresh) {
        y2[k2] = m_minthresh;
      }
    }
  }

 private:
  int m_km1;
  int m_km2;
  Scalar m_minthresh;
};

} //namespace li

#endif
