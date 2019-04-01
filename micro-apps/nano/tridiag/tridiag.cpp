#define AMB_NO_MPI
#include "../hommexx/dbg.hpp"

#include <cassert>
#include <limits>

#include <Kokkos_Core.hpp>

#ifdef HAVE_CUSPARSE
# include <cusparse.h>
#endif

/*
  A = [0 b0 c0; a1 b1 c1; ...; an bn 0]
  X fast index is j

  todo:
  x get test1 working on gpu
  x get a perf test going
  x try shared mem
  x impl pcr
  x make a perf test that covers thomas, cr, pcr
  x try LayoutStride, esp with pcr. nope: LayoutStride has overhead
    (revealed by order = {2,1,0}), so it's unusable.
  x make a cusparse perf test
  x 2x2 opt suggests i should try thomas in cr, after all. no.
  x cr_a1x1: hand-code dl,d,u layout since LayoutStride won't cut it.
  x switch interface to take dl,d,du
  x homme gpu version that packs A differently. also can put A and X updates
    together. with this version, nrhs = 1 case should call homme version. then
    check perf against cusparse nrhs=1 case.
  x raw cuda
  - perf test for cr_amxm vs cusparse
  - then opt cr_amxm, following a1x?? examples
  - make a blas perf test: dttrfb, dttrsb
  - call thomas with pack
  - make an analysis routine that returns a small struct of POD indicating which
    alg and with what parms to run for a given set of prob parms. then make a
    top-level routine that switches based on the struct.
 */

namespace util {
static inline double gettime () {
  timeval t;
  gettimeofday(&t, 0);
  static const double us = 1e6;
  return (t.tv_sec*us + t.tv_usec)/us;
}

std::string active_avx_string () {
  std::string s;
#if defined __AVX512F__
  s += "-AVX512F";
#endif
#if defined __AVX2__
  s += "-AVX2";
#endif
#if defined __AVX__
  s += "-AVX";
#endif
  return s;
}

bool eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

void expect_another_arg (int i, int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}
}

namespace ko {
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }
}

template <typename TeamMember>
KOKKOS_INLINE_FUNCTION
int get_thread_id_within_team (const TeamMember& team) {
  return team.team_rank();
}

template <typename TeamMember>
KOKKOS_INLINE_FUNCTION
int get_team_nthr (const TeamMember& team) {
  return team.team_size();
}

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
Kokkos::TeamPolicy<ExeSpace>
get_default_team_policy (const int nprob, const int nrow, const int nrhs) {
#ifdef KOKKOS_ENABLE_OPENMP
# ifdef MIMIC_GPU
  return Kokkos::TeamPolicy<ExeSpace>(nprob, omp_get_max_threads(), 1);
# else
  const int per = (omp_get_max_threads() + nprob - 1)/nprob;
  return Kokkos::TeamPolicy<ExeSpace>(nprob, per, 1);
# endif
#else
  return Kokkos::TeamPolicy<ExeSpace>(nprob, 1, 1);
#endif
}

#ifdef KOKKOS_ENABLE_CUDA
KOKKOS_INLINE_FUNCTION
int get_thread_id_within_team (const Kokkos::Impl::CudaTeamMember& team) {
#ifdef __CUDA_ARCH__
  // Can't use team.team_rank() here because vector direction also uses physical
  // threads but TeamMember types don't expose that information.
  return blockDim.x * threadIdx.y + threadIdx.x;
#else
  assert(0);
  return -1;
#endif
}

KOKKOS_INLINE_FUNCTION
int get_team_nthr (const Kokkos::Impl::CudaTeamMember& team) {
#ifdef __CUDA_ARCH__
  return blockDim.x * blockDim.y;
#else
  assert(0);
  return -1;
#endif
}

// This is not intended to be optimal for the tridiag problem; rather,
// it is intended to reflect what a column physics param would likely
// use.
template <>
Kokkos::TeamPolicy<Kokkos::Cuda>
get_default_team_policy<Kokkos::Cuda> (const int nprob, const int nrow, const int nrhs) {
  const int nwarp_per_gpu = 2560;
  const int tpw = 32;
  const int nwarp_per_team =
    std::min(// max warp/team
#ifdef WHITE
      16,
#else
      32,
#endif
      std::max(1,
      std::max((nrow + tpw - 1)/tpw,
      (nwarp_per_gpu + nprob - 1)/nprob)));
      // switch to (nprob, nwarp_per_gpu, tpw) to test vec
      return Kokkos::TeamPolicy<Kokkos::Cuda>(nprob, tpw*nwarp_per_team, 1);
      }
#endif

// NB: The caller must provide the team_barrier after this function
// returns before A is accessed.
template <typename TeamMember, typename TridiagDiag>
KOKKOS_INLINE_FUNCTION
void thomas_factorize (const TeamMember& team,
                       TridiagDiag dl, TridiagDiag d, TridiagDiag du) {
  const int nrow = d.extent_int(0);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  const auto f = [&] () {
    for (int i = 1; i < nrow; ++i) {
      dl(i) /= d(i-1);
      d (i) -= dl(i) * du(i-1);
    }
  };
  Kokkos::single(Kokkos::PerTeam(team), f);
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void thomas_solve (const TeamMember& team,
                   TridiagDiag dl, TridiagDiag d, TridiagDiag du,
                   DataArray X) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  const auto f = [&] (const int& j) {
    const auto g = [&] () {
      for (int i = 1; i < nrow; ++i)
        X(i,j) -= dl(i) * X(i-1,j);
      X(nrow-1,j) /= d(nrow-1);
      for (int i = nrow-1; i > 0; --i)
        X(i-1,j) = (X(i-1,j) - du(i-1) * X(i,j)) / d(i-1);
    };
    Kokkos::single(Kokkos::PerThread(team), g);
  };
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nrhs), f);
}

// Cyclic reduction at the Kokkos team level. Any (thread, vector)
// parameterization is intended to work.
template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void cr_a1xm (const TeamMember& team,
              TridiagDiag dl, TridiagDiag d, TridiagDiag du,
              DataArray X) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  assert(X.extent_int(0) == nrow);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  const int nrhs = X.extent_int(1);
  const int tid = get_thread_id_within_team(team);
  const int nthr = get_team_nthr(team);
  const int team_size = ko::min(nrhs, nthr);
  const int nteam = nthr / team_size;
  const int team_id = tid / team_size;
  const int team_tid = tid % team_size;
  const bool team_lead = tid % team_size == 0;
  int os = 1, stride;
  while ((stride = (os << 1)) < nrow) {
    const int sub_nrow = (nrow + stride - 1) / stride;
    const int nit = (sub_nrow + nteam - 1) / nteam;
    for (int i = stride*team_id, it = 0; it < nit; i += stride*nteam, ++it) {
      int im = i - os;
      int ip = i + os;
      Scalar f1 = 0, f2 = 0;
      const bool run = team_id < nteam && i < nrow;
      assert(team_id != 0 || run);
      if (run) {
        f1 = im >= 0   ? -dl(i)/d(im) : 0;
        f2 = ip < nrow ? -du(i)/d(ip) : 0;
        im = im >= 0   ? im : i;
        ip = ip < nrow ? ip : i;
        for (int j = team_tid; j < nrhs; j += team_size)
          X(i,j) += f1*X(im,j) + f2*X(ip,j);
      }
      // Update A only after all threads are done using current values.
      team.team_barrier();
      if (team_lead && run) {
        dl(i)  = f1*dl(im);
        du(i)  = f2*du(ip);
        d (i) += f1*du(im) + f2*dl(ip);
      }
    }
    os <<= 1;
    // Go down in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
  if (team_id == 0) {
    if (os >= nrow) {
      for (int j = team_tid; j < nrhs; j += team_size)
        X(0,j) /= d(0);
    } else {
      for (int j = team_tid; j < nrhs; j += team_size) {
        const auto
          det = d(0)*d(os) - du(0)*dl(os),
          x0 = X(0,j), x1 = X(os,j);
        X( 0,j) = (d(os)*x0 - du( 0)*x1)/det;
        X(os,j) = (d( 0)*x1 - dl(os)*x0)/det;
      }
    }
  }
  team.team_barrier();
  os >>= 1;
  assert(os < nrow);
  while (os) {
    stride = os << 1;
    if (team_id < nteam) {
      for (int i = stride*team_id + os; i < nrow; i += stride*nteam) {
        const int im = i - os;
        const int ip = i + os;
        assert(im >= 0 || ip < nrow);
        for (int j = team_tid; j < nrhs; j += team_size) {
          Scalar f = 0;
          f += im >=   0 ? dl(i)*X(im,j) : 0;
          f += ip < nrow ? du(i)*X(ip,j) : 0;
          X(i,j) = (X(i,j) - f)/d(i);
        }
      }
    }
    os >>= 1;
    // Go up in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void cr_amxm (const TeamMember& team,
              TridiagDiag dl, TridiagDiag d, TridiagDiag du,
              DataArray X) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(dl.extent_int(1) == nrhs);
  assert(d. extent_int(1) == nrhs);
  assert(du.extent_int(1) == nrhs);
  assert(dl.extent_int(0) == nrow);
  assert(d. extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X. extent_int(0) == nrow);
  const int tid = get_thread_id_within_team(team);
  const int nthr = get_team_nthr(team);
  const int team_size = ko::min(nrhs, nthr);
  const int nteam = nthr / team_size;
  const int team_id = tid / team_size;
  const int team_tid = tid % team_size;
  int os = 1, stride;
  while ((stride = (os << 1)) < nrow) {
    const int sub_nrow = (nrow + stride - 1) / stride;
    const int nit = (sub_nrow + nteam - 1) / nteam;
    for (int i = stride*team_id, it = 0; it < nit; i += stride*nteam, ++it) {
      const int im = i - os;
      const int ip = i + os;
      const bool run = team_id < nteam && i < nrow;
      assert(team_id != 0 || run);
      if (run) {
        if (im < 0) {
          assert(ip < nrow);
          for (int j = team_tid; j < nrhs; j += team_size) {
            const auto f = -du(i,j)/d(ip,j);
            du(i,j)  = f*du(ip,j);
            d (i,j) += f*dl(ip,j);
            X (i,j) += f*X (ip,j);
          }
        } else if (ip >= nrow) {
          for (int j = team_tid; j < nrhs; j += team_size) {
            const auto f = -dl(i,j)/d(im,j);
            dl(i,j)  = f*dl(im,j);
            d (i,j) += f*du(im,j);
            X (i,j) += f*X (im,j);
          }
        } else {
          for (int j = team_tid; j < nrhs; j += team_size) {
            const auto f1 = -dl(i,j)/d(im,j);
            dl(i,j)  = f1*dl(im,j);
            const auto f2 = -du(i,j)/d(ip,j);
            du(i,j)  = f2*du(ip,j);
            d (i,j) += f1*du(im,j) + f2*dl(ip,j);
            X (i,j) += f1*X (im,j) + f2*X(ip,j);
          }
        }
      }
    }
    os <<= 1;
    // Go down in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
  if (team_id == 0) {
    if (os >= nrow) {
      for (int j = team_tid; j < nrhs; j += team_size)
        X(0,j) /= d(0,j);
    } else {
      for (int j = team_tid; j < nrhs; j += team_size) {
        const auto
          det = d(0,j)*d(os,j) - du(0,j)*dl(os,j),
          x0 = X(0,j), x1 = X(os,j);
        X( 0,j) = (d(os,j)*x0 - du( 0,j)*x1)/det;
        X(os,j) = (d( 0,j)*x1 - dl(os,j)*x0)/det;
      }
    }
  }
  team.team_barrier();
  os >>= 1;
  assert(os < nrow);
  while (os) {
    stride = os << 1;
    if (team_id < nteam) {
      for (int i = stride*team_id + os; i < nrow; i += stride*nteam) {
        const int im = i - os;
        const int ip = i + os;
        if (im < 0) {
          assert(ip < nrow);
          for (int j = team_tid; j < nrhs; j += team_size) {
            const auto f = du(i,j)*X(ip,j);
            X(i,j) = (X(i,j) - f)/d(i,j);
          }
        } else if (ip >= nrow) {
          for (int j = team_tid; j < nrhs; j += team_size) {
            const auto f = dl(i,j)*X(im,j);
            X(i,j) = (X(i,j) - f)/d(i,j);
          }        
        } else {
          for (int j = team_tid; j < nrhs; j += team_size) {
            const auto f = dl(i,j)*X(im,j) + du(i,j)*X(ip,j);
            X(i,j) = (X(i,j) - f)/d(i,j);
          }
        }
      }
    }
    os >>= 1;
    // Go up in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void cr_a1x1 (const TeamMember& team,
              TridiagDiag dl, TridiagDiag d, TridiagDiag du,
              DataArray X) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  assert(dl.extent_int(0) == nrow);
  assert(d. extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X. extent_int(0) == nrow);
  const int team_id = get_thread_id_within_team(team);
  const int nteam = get_team_nthr(team);
  int os = 1, stride;
  // Go down reduction.
  while ((stride = (os << 1)) < nrow) {
    const int sub_nrow = (nrow + stride - 1) / stride;
    const int nit = (sub_nrow + nteam - 1) / nteam;
    for (int i = stride*team_id, it = 0; it < nit; i += stride*nteam, ++it) {
      int im = i - os;
      int ip = i + os;
      const bool run = i < nrow;
      assert(team_id != 0 || run);
      if (run) {
        // GPU does well with ternary ?: op. Use it throughout this
        // impl. It requires the trick noted in a few lines.
        const auto f1 = im >= 0   ? -dl(i)/d(im) : 0;
        const auto f2 = ip < nrow ? -du(i)/d(ip) : 0;
        // Trick to keep im, ip in bounds; the index is modified only
        // when the corresponding f is 0, so the resulting invalid
        // value is multipled by 0.
        im = im >= 0   ? im : i;
        ip = ip < nrow ? ip : i;
        dl(i)  = f1*dl(im);
        du(i)  = f2*du(ip);
        d (i) += f1*du(im) + f2*dl(ip);
        X (i) += f1*X (im) + f2*X (ip);
      }
    }
    os <<= 1;
    // Go down in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
  // Bottom 1 or 2 levels of the reduction. This could be folded into
  // the previous loop, but it's a slight opt to handle these cases
  // separately.
  if (team_id == 0) {
    if (os >= nrow) {
      X(0) /= d(0);
    } else {
      const auto
        det = d(0)*d(os) - du(0)*dl(os),
        x0 = X(0), x1 = X(os);
      X( 0) = (d(os)*x0 - du( 0)*x1)/det;
      X(os) = (d( 0)*x1 - dl(os)*x0)/det;
    }
  }
  team.team_barrier();
  os >>= 1;
  assert(os < nrow);
  // Go up reduction.
  while (os) {
    stride = os << 1;
    for (int i = stride*team_id + os; i < nrow; i += stride*nteam) {
      const int im = i - os;
      const int ip = i + os;
      assert(im >= 0 || ip < nrow);
      Scalar f = 0;
      f += im >=   0 ? dl(i)*X(im) : 0;
      f += ip < nrow ? du(i)*X(ip) : 0;
      X(i) = (X(i) - f)/d(i);
    }
    os >>= 1;
    // Go up in cyclic reduction level only when this level is complete.
    team.team_barrier();
  }
}

// Pure Cuda version to measure the Kokkos overhead. For this kind of
// computation, every little bit counts, so we need to know that
// overhead.
#ifdef KOKKOS_ENABLE_CUDA
template <typename Scalar>
__device__
void cr_a1x1p (const int nrow,
               Scalar* const dl, Scalar* const d, Scalar* const du,
               Scalar* const X) {
  const int team_id = threadIdx.x;
  const int nteam = blockDim.x;
  int os = 1, stride;
  while ((stride = (os << 1)) < nrow) {
    const int sub_nrow = (nrow + stride - 1) / stride;
    const int nit = (sub_nrow + nteam - 1) / nteam;
    for (int i = stride*team_id, it = 0; it < nit; i += stride*nteam, ++it) {
      int im = i - os;
      int ip = i + os;
      const bool run = i < nrow;
      assert(team_id != 0 || run);
      if (run) {
        const auto f1 = im >= 0   ? -dl[i]/d[im] : 0;
        const auto f2 = ip < nrow ? -du[i]/d[ip] : 0;
        im = im >= 0   ? im : i;
        ip = ip < nrow ? ip : i;
        dl[i]  = f1*dl[im];
        du[i]  = f2*du[ip];
        d [i] += f1*du[im] + f2*dl[ip];
        X [i] += f1*X [im] + f2*X [ip];
      }
    }
    os <<= 1;
    __syncthreads();
  }
  if (team_id == 0) {
    if (os >= nrow) {
      X[0] /= d[0];
    } else {
      const auto
        det = d[0]*d[os] - du[0]*dl[os],
        x0 = X[0], x1 = X[os];
      X[ 0] = (d[os]*x0 - du[ 0]*x1)/det;
      X[os] = (d[ 0]*x1 - dl[os]*x0)/det;
    }
  }
  __syncthreads();
  os >>= 1;
  assert(os < nrow);
  while (os) {
    stride = os << 1;
    for (int i = stride*team_id + os; i < nrow; i += stride*nteam) {
      const int im = i - os;
      const int ip = i + os;
      assert(im >= 0 || ip < nrow);
      Scalar f = 0;
      f += im >= 0   ? dl[i]*X[im] : 0;
      f += ip < nrow ? du[i]*X[ip] : 0;
      X[i] = (X[i] - f)/d[i];
    }
    os >>= 1;
    __syncthreads();
  }
}
#endif

template <typename TridiagDiag, typename XArray, typename YArray>
KOKKOS_INLINE_FUNCTION
int matvec (TridiagDiag dl, TridiagDiag d, TridiagDiag du, XArray X, YArray Y) {
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X.extent_int(0) == nrow);
  assert(Y.extent_int(0) == nrow);
  assert(Y.extent_int(1) == nrhs);
  if (nrow == 1) {
    for (int j = 0; j < nrhs; ++j)
      Y(0,j) = d(0) * X(0,j);
    return 0;
  }
  for (int j = 0; j < nrhs; ++j)
    Y(0,j) = d(0) * X(0,j) + du(0) * X(1,j);
  for (int i = 1; i < nrow-1; ++i)
    for (int j = 0; j < nrhs; ++j)
      Y(i,j) = (dl(i) * X(i-1,j) +
                d (i) * X(i  ,j) +
                du(i) * X(i+1,j));
  const int i = nrow-1;
  for (int j = 0; j < nrhs; ++j)
    Y(i,j) = dl(i) * X(i-1,j) + d(i) * X(i,j);
  return 0;  
}

KOKKOS_INLINE_FUNCTION
static unsigned int nextpow2 (unsigned int v) {
  // From https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

template <typename TeamMember, typename TridiagDiag, typename DataArray>
KOKKOS_INLINE_FUNCTION
void pcr (const TeamMember& team,
          TridiagDiag dl, TridiagDiag d, TridiagDiag du,
          DataArray X) {
  using Scalar = typename TridiagDiag::non_const_value_type;
  const int nrow = d.extent_int(0);
  const int nrhs = X.extent_int(1);
  const int nthr = get_team_nthr(team);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  assert(X.extent_int(0) == nrow);
  assert(nrow*nrhs <= nthr);
  const int tid = get_thread_id_within_team(team);
  const int team_size = nrhs;
  const int team_id = tid / team_size;
  const int team_tid = tid % team_size;
  const bool team_lead = tid % team_size == 0;
  const int i = team_id;
  const int j = team_tid;
  int n = nextpow2(nrow);
  int os = 1;
  while (n) {
    const int im = i - os;
    const int ip = i + os;
    Scalar a = 0, b, c = 0, f, y;
    if (i < nrow) {
      if (im >= 0 || ip < nrow) {
        b = d(i);
        y = X(i,j);
        if (im >= 0) {
          f  = -dl(i)/d(im);
          a  = f*dl(im);
          b += f*du(im);
          y += f*X (im,j);
        }
        if (ip < nrow) {
          f  = -du(i)/d(ip);
          b += f*dl(ip);
          c  = f*du(ip);
          y += f*X (ip,j);
        }
      } else {      
        b = 1;
        y = X(i,j)/d(i);
      }
    }
    team.team_barrier();
    if (i < nrow) {
      X(i,j) = y;
      if (team_lead) {
        dl(i) = a;
        d (i) = b;
        du(i) = c;
      }
    }
    team.team_barrier();
    os <<= 1;
    n >>= 1;
  }
}

#ifdef HAVE_CUSPARSE
/*
  cusparseStatus_t cusparseCreate(cusparseHandle_t *handle)
  cusparseStatus_t cusparseDestroy(cusparseHandle_t handle)
  cusparseStatus_t cusparseDgtsv2StridedBatch_bufferSizeExt(
    cusparseHandle_t handle,
    int m,
    const double *dl,
    const double  *d,
    const double *du,
    const double *x,
    int batchCount,
    int batchStride,
    size_t *bufferSizeInBytes)
  cusparseStatus_t cusparseDgtsv2StridedBatch(
    cusparseHandle_t handle,
    int m,
    const double *dl,
    const double  *d,
    const double *du,
    double *x,
    int batchCount,
    int batchStride,
    void *pBuffer)
 */
template <typename T>
cusparseStatus_t cusparseTgtsv2StridedBatch_bufferSizeExt(
  cusparseHandle_t, int, const T*, const T*, const T*, const T*, int, int, size_t*);
template <typename T>
cusparseStatus_t cusparseTgtsv2StridedBatch(
  cusparseHandle_t, int, const T*, const T*, const T*, T*, int, int, T*);
template <> cusparseStatus_t cusparseTgtsv2StridedBatch_bufferSizeExt<float> (
  cusparseHandle_t h, int m, const float* dl, const float* d, const float* du, const float* x,
  int bc, int bs, size_t* sz)
{ return cusparseSgtsv2StridedBatch_bufferSizeExt(h, m, dl, d, du, x, bc, bs, sz); }
template <> cusparseStatus_t cusparseTgtsv2StridedBatch<float> (
  cusparseHandle_t h, int m, const float* dl, const float* d, const float* du, float* x,
  int bc, int bs, float* b)
{ return cusparseSgtsv2StridedBatch(h, m, dl, d, du, x, bc, bs, b); }
template <> cusparseStatus_t cusparseTgtsv2StridedBatch_bufferSizeExt<double> (
  cusparseHandle_t h, int m, const double* dl, const double* d, const double* du, const double* x,
  int bc, int bs, size_t* sz)
{ return cusparseDgtsv2StridedBatch_bufferSizeExt(h, m, dl, d, du, x, bc, bs, sz); }
template <> cusparseStatus_t cusparseTgtsv2StridedBatch<double> (
  cusparseHandle_t h, int m, const double* dl, const double* d, const double* du, double* x,
  int bc, int bs, double* b)
{ return cusparseDgtsv2StridedBatch(h, m, dl, d, du, x, bc, bs, b); }

template <typename Scalar>
class CusparseSolver {
  int nrow_, nprob_;
  bool ok_;
  cusparseHandle_t handle_;
  Kokkos::View<Scalar*> buffer_;

  void init (const int& nrow, const int& nprob) {
    nrow_ = nrow;
    nprob_ = nprob;
    ok_ = true;
    if (nrow_ < 3) ok_ = false;
    auto status = cusparseCreate(&handle_);
    if (status != CUSPARSE_STATUS_SUCCESS) {
      pr("cusparseCreate" pu(status));
      ok_ = false;
      return;
    }
    size_t sz;
    status = cusparseTgtsv2StridedBatch_bufferSizeExt<Scalar>(
      handle_, nrow_, nullptr, nullptr, nullptr, nullptr, nprob_, nrow_, &sz);
    if (status != CUSPARSE_STATUS_SUCCESS) {
      pr("bufferSizeExt" pu(status));
      ok_ = false;
      return;
    }
    buffer_ = Kokkos::View<Scalar*>("buffer", sz/sizeof(Scalar));
  }

public:
  CusparseSolver (const int& nrow, const int& nprob) { init(nrow, nprob); }

  ~CusparseSolver () {
    cusparseDestroy(handle_);
  }

  template <typename Acs, typename Dcs>
  bool solve (const Acs& A, Dcs& X) {
    if ( ! ok_) return false;
    auto status = cusparseTgtsv2StridedBatch(
      handle_, nrow_, &A.impl_map().reference(0,0,0), &A.impl_map().reference(1,0,0),
      &A.impl_map().reference(2,0,0), X.data(), nprob_, nrow_, buffer_.data());
    if (status != CUSPARSE_STATUS_SUCCESS) {
      pr("solve" pu(status));
    }
    return status = CUSPARSE_STATUS_SUCCESS;
  }

  template <typename Astd, typename Acs>
  void to_matrix (const Astd& a, const Acs& b, const int nrhs) const {
    const int nprob = nprob_;
    const int nrow = nrow_;
    const int nrp = nrow*nprob;
    assert(a.extent_int(0) == nprob/nrhs);
    assert(a.extent_int(2) == nrow);
    assert(a.extent_int(1) == 3);
    assert(b.extent_int(0) == 3);
    assert(b.extent_int(1) == nprob);
    assert(b.extent_int(2) == nrow);
    auto f = KOKKOS_LAMBDA (const int i) {
      const int c = i/nrp;
      const int prob = (i % nrp) / nrow;
      const int row = i % nrow;
      b(c,prob,row) = a(prob/nrhs,c,row);
    };
    Kokkos::parallel_for(3*nrow*nprob, f);
  }

  template <typename Dstd, typename Dcs>
  void to_data (const Dstd& a, const Dcs& b) const {
    const int nrhs = a.extent_int(2);
    const int nprob = nprob_;
    const int nrow = nrow_;
    assert(a.extent_int(0) == nprob/nrhs);
    assert(a.extent_int(1) == nrow);
    assert(a.extent_int(2) == nrhs);
    assert(b.extent_int(0) == nprob);
    assert(b.extent_int(1) == nrow);
    auto f = KOKKOS_LAMBDA (const int i) {
      const int prob = i / nrow;
      const int row = i % nrow;
      const int rhs = prob % nrhs;
      b(prob,row) = a(prob/nrhs,row,rhs);
    };
    Kokkos::parallel_for(nrow*nprob, f);
  }

  template <typename Dcs, typename Dstd>
  void from_data (const Dcs& a, const Dstd& b) const {
    const int nrhs = b.extent_int(2);
    const int nprob = nprob_;
    const int nrow = nrow_;
    assert(b.extent_int(0) == nprob/nrhs);
    assert(b.extent_int(1) == nrow);
    assert(b.extent_int(2) == nrhs);
    assert(a.extent_int(0) == nprob);
    assert(a.extent_int(1) == nrow);
    auto f = KOKKOS_LAMBDA (const int i) {
      const int prob = i / nrow;
      const int row = i % nrow;
      const int rhs = prob % nrhs;
      b(prob/nrhs,row,rhs) = a(prob,row);
    };
    Kokkos::parallel_for(nrow*nprob, f);
  }
};
#endif // HAVE_CUSPARSE

template <typename TridiagDiag>
KOKKOS_INLINE_FUNCTION
void fill_tridiag_matrix (TridiagDiag dl, TridiagDiag d, TridiagDiag du,
                          const int& seed) {
  const int nrow = d.extent_int(0);
  assert(dl.extent_int(0) == nrow);
  assert(du.extent_int(0) == nrow);
  for (int i = 0; i < nrow; ++i) {
    const int k = seed + i;
    dl(i) = (k % 5 == 0 ? -1 : 1) * 1.3 * (0.1 + ((k*k) % 11));
    du(i) = (k % 7 == 0 ? -1 : 1) * 1.7 * (0.2 + ((k*k) % 13));
    d (i) = ((k % 3 == 0 ? -1 : 1) *
             (0.7 + std::abs(dl(i)) + std::abs(du(i)) + (k % 17)));
  }
}

template <typename DataArray>
KOKKOS_INLINE_FUNCTION
void fill_data_matrix (DataArray X, const int& seed) {
  const int nrow = X.extent_int(0);
  const int nrhs = X.extent_int(1);
  for (int i = 0; i < nrow; ++i)
    for (int j = 0; j < nrhs; ++j)
      X(i,j) = (((7*i + 11*j + 3*i*j) % 3 == 0 ? -1 : 1) *
                1.7 * ((17*(i - 19) + 13*(j - 11) + 5*(i - 5)*(j - 7) + seed) % 47));
}

typedef double Real;

template <typename Array>
Real reldif (const Array& a, const Array& b) {
  assert(a.extent_int(0) == b.extent_int(0));
  assert(a.extent_int(1) == b.extent_int(1));
  assert(a.rank == 2);
  assert(b.rank == 2);
  Real num = 0, den = 0;
  for (int i = 0; i < a.extent_int(0); ++i)
    for (int j = 0; j < a.extent_int(1); ++j) {
      num = std::max(num, std::abs(a(i,j) - b(i,j)));
      den = std::max(den, std::abs(a(i,j)));
    }
  return num/den;
}

struct Solver {
  enum Enum { thomas, cr_a1xm, cr_amxm, cr_a1x1, cr_a1x1p, pcr, cusparse, error };
  static std::string convert (Enum e) {
    switch (e) {
    case thomas: return "thomas";
    case cr_a1xm: return "cr_a1xm";
    case cr_amxm: return "cr_amxm";
    case cr_a1x1: return "cr_a1x1";
    case cr_a1x1p: return "cr_a1x1p";
    case pcr: return "pcr";
    case cusparse: return "cusparse";
    default:
      assert(0);
      return "";
    }
  }
  static Enum convert (const std::string& s) {
    if (s == "thomas") return thomas;
    if (s == "cr_a1xm") return cr_a1xm;
    if (s == "cr_amxm") return cr_amxm;
    if (s == "cr_a1x1") return cr_a1x1;
    if (s == "cr_a1x1p") return cr_a1x1p;
    if (s == "pcr") return pcr;
    if (s == "cusparse") return cusparse;
    return error;
  }
  static Enum all[];
};

Solver::Enum Solver::all[] = {
  thomas, cr_a1xm, cr_amxm, cr_a1x1, pcr
#ifdef HAVE_CUSPARSE
  //, cusparse
#endif
};

// LayoutStride has too much overhead to use it. Instead we use the
// (dl,d,du) triple rather than a matrix A to get the layout
// flexibility we need.
//#define STRIDE
#ifdef STRIDE
using BulkLayout = Kokkos::LayoutStride;
using TeamLayout = Kokkos::LayoutLeft;
#else
using BulkLayout = Kokkos::LayoutRight;
using TeamLayout = Kokkos::LayoutRight;
#endif

template <typename ScalarType>
using TridiagArray = Kokkos::View<ScalarType*[3], TeamLayout>;
template <typename ScalarType>
using DataArray = Kokkos::View<ScalarType**, TeamLayout>;

static void test1 () {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using MT = typename TeamPolicy::member_type;
  const auto policy = get_default_team_policy<>(1, 128, 1);
  const int nthr = policy.team_size();
  for (const auto solver : Solver::all) {
    for (int nrow = 1; nrow <= 129; ++nrow) {
      const int nrhs_max = solver != Solver::pcr ? 60 : 8;
      const int nrhs_inc = solver != Solver::pcr ? 11 : 1;
      for (int nrhs = 1; nrhs <= nrhs_max; nrhs += nrhs_inc) {
        if (solver == Solver::pcr && nrhs*nrow >= nthr) continue;
        if (solver == Solver::cr_amxm && nrhs > 1) continue;
        if (solver == Solver::cr_a1x1 && nrhs > 1) continue;
        TridiagArray<Real> A("A", nrow), Acopy("A", A.extent(0));
        DataArray<Real> B("B", A.extent(0), nrhs), X("X", B.extent(0), B.extent(1)),
          Y("Y", X.extent(0), X.extent(1));
        auto Am = Kokkos::create_mirror_view(A);
        auto Bm = Kokkos::create_mirror_view(B);
        {
          const auto dl = Kokkos::subview(Am, Kokkos::ALL(), 0);
          const auto d  = Kokkos::subview(Am, Kokkos::ALL(), 1);
          const auto du = Kokkos::subview(Am, Kokkos::ALL(), 2);
          fill_tridiag_matrix(dl, d, du, nrhs);
        }
        fill_data_matrix(Bm, nrhs);
        Kokkos::deep_copy(A, Am);
        Kokkos::deep_copy(B, Bm);
        Kokkos::deep_copy(Acopy, A);
        Kokkos::deep_copy(X, B);
        const auto policy = get_default_team_policy<>(1, nrow, nrhs);
        switch (solver) {
        case Solver::thomas: {
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            const auto dl = Kokkos::subview(A, Kokkos::ALL(), 0);
            const auto d  = Kokkos::subview(A, Kokkos::ALL(), 1);
            const auto du = Kokkos::subview(A, Kokkos::ALL(), 2);
            thomas_factorize(team, dl, d, du);
            team.team_barrier();
            thomas_solve(team, dl, d, du, X);
          };
          Kokkos::parallel_for(policy, f);
        } break;
        case Solver::cr_a1xm: {
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            const auto dl = Kokkos::subview(A, Kokkos::ALL(), 0);
            const auto d  = Kokkos::subview(A, Kokkos::ALL(), 1);
            const auto du = Kokkos::subview(A, Kokkos::ALL(), 2);
            cr_a1xm(team, dl, d, du, X);
          };
          Kokkos::parallel_for(policy, f);
        } break;
        case Solver::cr_amxm: {
          Kokkos::View<Real***, TeamLayout> Am(A.data(), nrow, 3, nrhs);
          const auto dl = Kokkos::subview(Am, Kokkos::ALL(), 0, Kokkos::ALL());
          const auto d  = Kokkos::subview(Am, Kokkos::ALL(), 1, Kokkos::ALL());
          const auto du = Kokkos::subview(Am, Kokkos::ALL(), 2, Kokkos::ALL());
          const auto f = KOKKOS_LAMBDA (const MT& team) { cr_amxm(team, dl, d, du, X); };
          Kokkos::parallel_for(policy, f);
        } break;
        case Solver::cr_a1x1: {
          Kokkos::View<Real**, TeamLayout> Am(A.data(), nrow, 3);
          const auto dl = Kokkos::subview(Am, Kokkos::ALL(), 0);
          const auto d  = Kokkos::subview(Am, Kokkos::ALL(), 1);
          const auto du = Kokkos::subview(Am, Kokkos::ALL(), 2);
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            cr_a1x1(team, dl, d, du, Kokkos::subview(X, Kokkos::ALL(), 0));
          };
          Kokkos::parallel_for(policy, f);
        } break;
        case Solver::pcr: {
          Kokkos::View<Real**, TeamLayout> Am(A.data(), nrow, 3);
          const auto dl = Kokkos::subview(Am, Kokkos::ALL(), 0);
          const auto d  = Kokkos::subview(Am, Kokkos::ALL(), 1);
          const auto du = Kokkos::subview(Am, Kokkos::ALL(), 2);
          const auto f = KOKKOS_LAMBDA (const MT& team) {
            pcr(team, dl, d, du, X);
          };
          Kokkos::parallel_for(policy, f);
        } break;
        }
        Real re;
        {
          auto Acopym = Kokkos::create_mirror_view(Acopy);
          const auto dl = Kokkos::subview(Acopym, Kokkos::ALL(), 0);
          const auto d  = Kokkos::subview(Acopym, Kokkos::ALL(), 1);
          const auto du = Kokkos::subview(Acopym, Kokkos::ALL(), 2);
          auto Xm = Kokkos::create_mirror_view(X);
          auto Ym = Kokkos::create_mirror_view(Y);
          Kokkos::deep_copy(Acopym, Acopy);
          Kokkos::deep_copy(Xm, X);
          matvec(dl, d, du, Xm, Ym);
          re = reldif(Bm, Ym);
        }
        if (re > 50*std::numeric_limits<Real>::epsilon())
          std::cout << "test1: solver " << Solver::convert(solver)
                    << " nrow " << nrow << " nrhs " << nrhs
                    << " re " << re << "\n";
      }
    }
  }
}

struct Input {
  Solver::Enum method;
  int nprob, nrow, nrhs, nwarp;
  bool use_scratch;

  Input ()
    : method(Solver::thomas), nprob(2048), nrow(128), nrhs(43), nwarp(-1),
      use_scratch(false)
  {}

  bool parse (int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
      if (util::eq(argv[i], "-m", "--method")) {
        util::expect_another_arg(i, argc);
        method = Solver::convert(argv[++i]);
        if (method == Solver::error) {
          std::cout << "Not a solver: " << argv[i] << "\n";
          return false;
        }
      } else if (util::eq(argv[i], "-np", "--nprob")) {
        util::expect_another_arg(i, argc);
        nprob = std::atoi(argv[++i]);
      } else if (util::eq(argv[i], "-nr", "--nrow")) {
        util::expect_another_arg(i, argc);
        nrow = std::atoi(argv[++i]);
      } else if (util::eq(argv[i], "-nc", "--nrhs")) {
        util::expect_another_arg(i, argc);
        nrhs = std::atoi(argv[++i]);
      } else if (util::eq(argv[i], "-nw", "--nwarp")) {
        util::expect_another_arg(i, argc);
        nwarp = std::atoi(argv[++i]);
#if 0
      } else if (util::eq(argv[i], "-s", "--scratch")) {
        use_scratch = true;
#endif
      } else {
        std::cout << "Unexpected arg: " << argv[i] << "\n";
        return false;
      }
    }
    return true;
  }
};

#ifdef KOKKOS_ENABLE_CUDA
template <typename Real>
__global__
void call_cr_a1x1p (const int nrow, Real* const Adata, Real* const Xdata) {
  const int ip = blockIdx.x;
  Real *dl = Adata + ip*3*nrow, *d = dl+nrow, *du = d+nrow, *x = Xdata + ip*nrow;
  cr_a1x1p(nrow, dl, d, du, x);
}

template <typename AT, typename XT>
double dispatch_cr_a1x1p (const Input& in, AT& A, XT& X) {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  TeamPolicy policy(in.nprob, 32*in.nwarp, 1);
  Kokkos::fence();
  double t0 = util::gettime();
  const int nrow = in.nrow;
  Real* const Adata = A.data();
  Real* const Xdata = X.data();
  call_cr_a1x1p<<<in.nprob,32*in.nwarp>>>(nrow, Adata, Xdata);
  Kokkos::fence();
  double t1 = util::gettime();
  return t1 - t0;
}
#endif

template <typename ScalarType>
using TridiagArrays = Kokkos::View<ScalarType***, BulkLayout>;
template <typename ScalarType>
using DataArrays = Kokkos::View<ScalarType***, BulkLayout>;

static void run (const Input& in) {
  using Kokkos::subview;
  using Kokkos::ALL;
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  std::cout << "run: solver " << Solver::convert(in.method) << " nprob " << in.nprob
            << " nrow " << in.nrow << " nrhs " << in.nrhs << "\n";
#ifdef STRIDE
  const int order[] = {2,0,1};
  const int Adim[] = {in.nprob, in.nrow, 3};
  const int Xdim[] = {in.nprob, in.nrow, in.nrhs};
  const auto Astride = Kokkos::LayoutStride::order_dimensions(3, order, Adim);
  const auto Xstride = Kokkos::LayoutStride::order_dimensions(3, order, Xdim);
  TridiagArrays<Real> A("A", Astride), Acopy("Acopy", Astride);
  DataArrays<Real> B("B", Xstride), X("X", Xstride), Y("Y", Xstride);
#else
  TridiagArrays<Real> A("A", in.nprob, 3, in.nrow),
    Acopy("Acopy", A.extent(0), A.extent(1), A.extent(2));
  DataArrays<Real> B("B", in.nprob, in.nrow, in.nrhs), X("X", in.nprob, in.nrow, in.nrhs),
    Y("Y", in.nprob, in.nrow, in.nrhs);
#endif
  auto Am = Kokkos::create_mirror_view(A);
  auto Bm = Kokkos::create_mirror_view(B);
  for (int i = 0; i < in.nprob; ++i) {
    const auto dl = Kokkos::subview(Am, i, 0, Kokkos::ALL());
    const auto d  = Kokkos::subview(Am, i, 1, Kokkos::ALL());
    const auto du = Kokkos::subview(Am, i, 2, Kokkos::ALL());
    fill_tridiag_matrix(dl, d, du, i);
    fill_data_matrix(subview(Bm, i, ALL(), ALL()), i);
  }
  Kokkos::deep_copy(A, Am);
  Kokkos::deep_copy(B, Bm);
  Kokkos::deep_copy(Acopy, A);
  Kokkos::deep_copy(X, B);
  auto policy = in.nwarp < 0 ? get_default_team_policy<>(in.nprob, in.nrow, in.nrhs) :
    TeamPolicy(in.nprob, 32*in.nwarp, 1);
  pr(puf(policy.league_size()) pu(policy.team_size()));
  if (in.method == Solver::pcr && policy.team_size() < in.nrow*in.nrhs) {
    std::stringstream ss;
    ss << "PCR requires nthr >= nrow*nrhs but nthr " << policy.team_size()
       << " nrow " << in.nrow << " nrhs " << in.nrhs << "\n";
    throw std::runtime_error(ss.str());
  }
  if ((in.method == Solver::cr_amxm || in.method == Solver::cr_a1x1) && in.nrhs != 1)
    throw std::runtime_error("In this performance test, cr_amxm and cr_a1x1 require nrhs 1.");
  Kokkos::fence();
  double t0, t1;
  switch (in.method) {
  case Solver::thomas: {
    t0 = util::gettime();
    const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
      const int ip = team.league_rank();
      const auto dl = Kokkos::subview(A, ip, 0, Kokkos::ALL());
      const auto d  = Kokkos::subview(A, ip, 1, Kokkos::ALL());
      const auto du = Kokkos::subview(A, ip, 2, Kokkos::ALL());
      thomas_factorize(team, dl, d, du);
      team.team_barrier();
      thomas_solve(team, dl, d, du, subview(X, team.league_rank(), ALL(), ALL()));
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    t1 = util::gettime();
  } break;
  case Solver::cr_a1xm: {
    t0 = util::gettime();
    const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
      const int ip = team.league_rank();
      const auto dl = Kokkos::subview(A, ip, 0, Kokkos::ALL());
      const auto d  = Kokkos::subview(A, ip, 1, Kokkos::ALL());
      const auto du = Kokkos::subview(A, ip, 2, Kokkos::ALL());
      cr_a1xm(team,
              dl, d, du,
              subview(X, ip, ALL(), ALL()));
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    t1 = util::gettime();
  } break;
  case Solver::cr_amxm: {
    t0 = util::gettime();
    const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
      const int ip = team.league_rank();
      Kokkos::View<Real***, TeamLayout, Kokkos::MemoryUnmanaged>
        Am(&A(ip,0,0), 3, in.nrow, 1);
      const auto dl = Kokkos::subview(Am, 0, Kokkos::ALL(), Kokkos::ALL());
      const auto d  = Kokkos::subview(Am, 1, Kokkos::ALL(), Kokkos::ALL());
      const auto du = Kokkos::subview(Am, 2, Kokkos::ALL(), Kokkos::ALL());
      cr_amxm(team, dl, d, du, subview(X, ip, ALL(), ALL()));
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    t1 = util::gettime();
  } break;
  case Solver::cr_a1x1: {
    t0 = util::gettime();
    Real* const Adata = A.data();
    Real* const Xdata = X.data();
    const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
      const int ip = team.league_rank();
      const auto dl = Kokkos::subview(A, ip, 0, Kokkos::ALL());
      const auto d  = Kokkos::subview(A, ip, 1, Kokkos::ALL());
      const auto du = Kokkos::subview(A, ip, 2, Kokkos::ALL());
      const auto x = subview(X, ip, ALL(), 0);
      cr_a1x1(team, dl, d, du, x);
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    t1 = util::gettime();
  } break;
#ifdef KOKKOS_ENABLE_CUDA
  case Solver::cr_a1x1p: {
    t0 = 0;
    t1 = dispatch_cr_a1x1p(in, A, X);
  } break;
#endif
  case Solver::pcr: {
    t0 = util::gettime();
    const auto f = KOKKOS_LAMBDA (const typename TeamPolicy::member_type& team) {
      const int ip = team.league_rank();
      const auto dl = Kokkos::subview(A, ip, 0, Kokkos::ALL());
      const auto d  = Kokkos::subview(A, ip, 1, Kokkos::ALL());
      const auto du = Kokkos::subview(A, ip, 2, Kokkos::ALL());
      pcr(team, dl, d, du, subview(X, ip, ALL(), 0));
    };
    Kokkos::parallel_for(policy, f);
    Kokkos::fence();
    t1 = util::gettime();
  } break;
#ifdef HAVE_CUSPARSE
  case Solver::cusparse: {
    Kokkos::View<Real***, Kokkos::LayoutRight> Acs("Acs", 3, in.nrhs*in.nprob, in.nrow);
    Kokkos::View<Real**, Kokkos::LayoutRight> Xcs("Xcs", in.nrhs*in.nprob, in.nrow);
    CusparseSolver<Real> cs(in.nrow, in.nrhs*in.nprob);
    cs.to_matrix(A, Acs, in.nrhs);
    cs.to_data(X, Xcs);
    Kokkos::fence();
    t0 = util::gettime();
    cs.solve(Acs, Xcs);
    Kokkos::fence();
    t1 = util::gettime();
    cs.from_data(Xcs, X);
  } break;
#endif
  }
  const auto et = t1 - t0;
  printf("run: et %1.3e et/datum %1.3e\n", et, (et/(in.nprob*in.nrow*in.nrhs)));
  Real re;
  {
    auto Acopym = Kokkos::create_mirror_view(Acopy);
    auto Xm = Kokkos::create_mirror_view(X);
    auto Ym = Kokkos::create_mirror_view(Y);
    Kokkos::deep_copy(Acopym, Acopy);
    Kokkos::deep_copy(Xm, X);
    const auto ip = in.nprob-1;
    const auto dl = Kokkos::subview(Acopym, ip, 0, Kokkos::ALL());
    const auto d  = Kokkos::subview(Acopym, ip, 1, Kokkos::ALL());
    const auto du = Kokkos::subview(Acopym, ip, 2, Kokkos::ALL());
    matvec(dl, d, du,
           subview(Xm, in.nprob-1, ALL(), ALL()),
           subview(Ym, in.nprob-1, ALL(), ALL()));
    re = reldif(
      subview(Bm, in.nprob-1, ALL(), ALL()),
      subview(Ym, in.nprob-1, ALL(), ALL()));
  }
  if (re > 50*std::numeric_limits<Real>::epsilon())
    std::cout << "run: " << " re " << re << "\n";
}

int main (int argc, char** argv) {
  int stat = 0;
  Kokkos::initialize(argc, argv); {
    if (argc > 1) {
      Input in;
      stat = in.parse(argc, argv);
      if (stat)
        run(in);
    } else
      test1();  
  } Kokkos::finalize();
  return stat;
}
