#include "types.hpp"
#include "util.hpp"
#include "initial_conditions.hpp"
#include "micro_sed_vanilla_kokkos.hpp"
#include "micro_kokkos.hpp"

#include <vector>
#include <iostream>
#include <exception>

#ifdef FPE
# include <xmmintrin.h>
#endif

// Nothing in this file is intended to be indicative of the eventual perf
// portable impl. E.g., Kokkos is not used here. The objective is simply to test
// for input/output regression errors.

extern "C" {
  void p3_init();
  void micro_sed_func_c(
    Int kts, Int kte, Int kdir, Int ni, Int nk, Int its, Int ite, Real dt,
    Real* qr, Real* nr, const Real* th, const Real* dzq, const Real* pres,
    Real* prt_liq);
}

struct TransposeDirection {
  enum Enum { c2f, f2c };
};

template <TransposeDirection::Enum direction, typename Scalar>
void transpose_layout (const ic::MicroSedData<Scalar>& s,
                        ic::MicroSedData<Scalar>& d) {
  const auto transpose = [&] (const Scalar* sv, Scalar* dv) {
    for (Int k = 0; k < s.nk; ++k)
      for (Int i = 0; i < s.ni; ++i)
        if (direction == TransposeDirection::c2f)
          dv[d.ni*k + i] = sv[d.nk*i + k];
        else
          dv[d.nk*i + k] = sv[d.ni*k + i];
  };
  transpose(s.qr, d.qr);
  transpose(s.nr, d.nr);
  transpose(s.th, d.th);
  transpose(s.dzq, d.dzq);
  transpose(s.pres, d.pres);
  for (Int i = 0; i < s.ni; ++i) d.prt_liq[i] = s.prt_liq[i];
}

// MicroSedData <-> micro_sed.f90 format.
template <typename Scalar>
struct FortranBridge : public ic::MicroSedData<Scalar> {
  FortranBridge (const ic::MicroSedData<Scalar>& d)
    : ic::MicroSedData<Scalar>(d)
  { transpose_layout<TransposeDirection::c2f>(d, *this); }

  void sync_to (ic::MicroSedData<Scalar>& d) {
    transpose_layout<TransposeDirection::f2c>(*this, d);
  }
};

template <typename Scalar>
void micro_sed_func (ic::MicroSedData<Scalar>& d, FortranBridge<Scalar>& b) {
  micro_sed_func_c(1, d.nk,
                   d.reverse ? -1 : 1,
                   d.ni, d.nk, 1, d.ni, d.dt, b.qr, b.nr, b.th,
                   b.dzq, b.pres, b.prt_liq);
  b.sync_to(d);
}

// MicroSedData <-> micro_sed_vanilla format.
template <typename Scalar>
class VanillaCppBridge
{
 public:
  VanillaCppBridge(const ic::MicroSedData<Scalar>& d) :
    qr(d.ni, std::vector<Scalar>(d.nk)),
    nr(d.ni, std::vector<Scalar>(d.nk)),
    th(d.ni, std::vector<Scalar>(d.nk)),
    dzq(d.ni, std::vector<Scalar>(d.nk)),
    pres(d.ni, std::vector<Scalar>(d.nk)),
    prt_liq(d.ni)
  {
    sync_from(d);
  }

  void sync_from(const ic::MicroSedData<Scalar>& d)
  {
    for (int i = 0; i < d.ni; ++i) {
      for (int k = 0; k < d.nk; ++k) {
        qr[i][k]   = d.qr[i*d.nk + k];
        nr[i][k]   = d.nr[i*d.nk + k];
        th[i][k]   = d.th[i*d.nk + k];
        dzq[i][k]  = d.dzq[i*d.nk + k];
        pres[i][k] = d.pres[i*d.nk + k];
      }
      prt_liq[i] = d.prt_liq[i];
    }
  }

  void sync_to(ic::MicroSedData<Scalar>& d)
  {
    for (int i = 0; i < d.ni; ++i) {
      for (int k = 0; k < d.nk; ++k) {
        d.qr[i*d.nk + k]   = qr[i][k];
        d.nr[i*d.nk + k]   = nr[i][k];
        d.th[i*d.nk + k]   = th[i][k];
        d.dzq[i*d.nk + k]  = dzq[i][k];
        d.pres[i*d.nk + k] = pres[i][k];
      }
      d.prt_liq[i] = prt_liq[i];
    }
  }

  vector_2d_t<Scalar> qr, nr, th, dzq, pres;
  std::vector<Scalar> prt_liq;
};

template <typename Scalar>
class KokkosCppBridge {
public:
  kokkos_2d_t<Real> qr, nr, th, dzq, pres;
  kokkos_1d_t<Real> prt_liq;

  KokkosCppBridge(const ic::MicroSedData<Scalar>& d)
    : qr("qr", d.ni, d.nk),
      nr("nr", d.ni, d.nk),
      th("th", d.ni, d.nk),
      dzq("dzq", d.ni, d.nk),
      pres("pres", d.ni, d.nk),
      prt_liq("prt_liq", d.ni)
  {
    sync_from(d);
  }

  void sync_from(const ic::MicroSedData<Scalar>& d)
  {
    auto mirror_qr      = Kokkos::create_mirror_view(qr);
    auto mirror_nr      = Kokkos::create_mirror_view(nr);
    auto mirror_th      = Kokkos::create_mirror_view(th);
    auto mirror_dzq     = Kokkos::create_mirror_view(dzq);
    auto mirror_pres    = Kokkos::create_mirror_view(pres);
    auto mirror_prt_liq = Kokkos::create_mirror_view(prt_liq);

    for (int i = 0; i < d.ni; ++i) {
      for (int k = 0; k < d.nk; ++k) {
        mirror_qr  (i, k) = d.qr  [i*d.nk + k];
        mirror_nr  (i, k) = d.nr  [i*d.nk + k];
        mirror_th  (i, k) = d.th  [i*d.nk + k];
        mirror_dzq (i, k) = d.dzq [i*d.nk + k];
        mirror_pres(i, k) = d.pres[i*d.nk + k];
      }
      mirror_prt_liq(i) = d.prt_liq[i];
    }

    Kokkos::deep_copy(qr, mirror_qr);
    Kokkos::deep_copy(nr, mirror_nr);
    Kokkos::deep_copy(th, mirror_th);
    Kokkos::deep_copy(dzq, mirror_dzq);
    Kokkos::deep_copy(pres, mirror_pres);
  }

  void sync_to(ic::MicroSedData<Scalar>& d) const
  {
    auto mirror_qr      = Kokkos::create_mirror_view(qr);
    auto mirror_nr      = Kokkos::create_mirror_view(nr);
    auto mirror_th      = Kokkos::create_mirror_view(th);
    auto mirror_dzq     = Kokkos::create_mirror_view(dzq);
    auto mirror_pres    = Kokkos::create_mirror_view(pres);
    auto mirror_prt_liq = Kokkos::create_mirror_view(prt_liq);

    Kokkos::deep_copy(mirror_qr, qr);
    Kokkos::deep_copy(mirror_nr, nr);
    Kokkos::deep_copy(mirror_th, th);
    Kokkos::deep_copy(mirror_dzq, dzq);
    Kokkos::deep_copy(mirror_pres, pres);

    for (int i = 0; i < d.ni; ++i) {
      for (int k = 0; k < d.nk; ++k) {
        d.qr  [i*d.nk + k] = mirror_qr  (i, k);
        d.nr  [i*d.nk + k] = mirror_nr  (i, k);
        d.th  [i*d.nk + k] = mirror_th  (i, k);
        d.dzq [i*d.nk + k] = mirror_dzq (i, k);
        d.pres[i*d.nk + k] = mirror_pres(i, k);
      }
      d.prt_liq[i] = mirror_prt_liq(i);
    }

  }
};

template <typename Scalar>
void micro_sed_func_cpp (ic::MicroSedData<Scalar>& d, VanillaCppBridge<Scalar>& bridge)
{
  p3::micro_sed_vanilla::micro_sed_func_vanilla<Scalar>( d.reverse ? d.nk : 1, d.reverse ? 1 : d.nk,
                                                        d.ni, d.nk, 1, d.ni, d.dt, bridge.qr, bridge.nr, bridge.th,
                                                        bridge.dzq, bridge.pres, bridge.prt_liq);

  bridge.sync_to(d);
}

template <typename Scalar>
void micro_sed_func_cpp_kokkos (ic::MicroSedData<Scalar>& d, KokkosCppBridge<Scalar>& bridge, p3::micro_sed_vanilla::MicroSedFuncVanillaKokkos<Real>& msvk)
{
  msvk.micro_sed_func_vanilla_kokkos( d.reverse ? d.nk : 1, d.reverse ? 1 : d.nk,
                                      d.ni, d.nk, 1, d.ni, d.dt, bridge.qr, bridge.nr, bridge.th,
                                      bridge.dzq, bridge.pres, bridge.prt_liq);

  bridge.sync_to(d);
}

template <typename Scalar>
struct MicroSedObserver {
  virtual ~MicroSedObserver () {}
  virtual void observe (const ic::MicroSedData<Scalar>& d) = 0;
};

struct BaselineConsts {
  // Take nstep outer steps so we can record output throughout the d.dt run.
  static constexpr Int nstep = 30;
};

template <typename Scalar>
void run_over_parameter_sets (MicroSedObserver<Scalar>& o, const Int ncol) {
  /* Co_max, the Courant number CN, sets the number of substeps. This can be
     estimated as follows:
       Given:
         dz [m]. Let's say clouds go up to ~10km. 100 uniform lvls => dz ~ 100m.
         V_qr [m/s]. Max is ~9 m/s.
     Then
       CN ~ V_qr dt / dz   [CN] = 1
       CN = 1 => dt = dz/V_qr ~ 10s
     So to get ~30 substeps, set dt to ~300s. This is very rough since dz is not
     uniform, etc.
       Here, I'm choosing dt to be really large, to get almost all the rain to
     sediment out. Another way to think of this is we can run the model with
     this dt and have useful computations all the way out to this dt. Obviously
     this dt must be adjusted if any of the other ICs are changed.
       There's also no reason to run this dt all at once. We can check against
     baselines throughout, e.g., to see growth in error.
   */
  const Real dt_tot = 300 * BaselineConsts::nstep; // s

  // will init both fortran and c
  p3::micro_sed_vanilla::p3_init_cpp_kokkos<Scalar>();

  ic::MicroSedData<Scalar> d(ncol, 111);
  d.dt = dt_tot;
  populate(d);

  o.observe(d);

  const auto d_rev(reverse_k(d));
  o.observe(d_rev);

  p3::micro_sed_vanilla::p3_deinit_cpp_kokkos<Scalar>();
}

struct BaselineObserver : public MicroSedObserver<Real> {
  Int nerr;

  BaselineObserver () : nerr(0) {}

  Int run (const Int ncolumns) {
    run_over_parameter_sets(*this, ncolumns);
    return nerr;
  }
};

template <typename Scalar>
static Int generate_baseline (const std::string& bfn) {
  struct Observer : public BaselineObserver {
    util::FILEPtr fid;

    Observer (const std::string& bfn) {
      fid = util::FILEPtr(fopen(bfn.c_str(), "w"));
      micro_throw_if( ! fid, "generate_baseline can't write " << bfn);
    }

    virtual void observe (const ic::MicroSedData<Scalar>& d_ic) override {
      // Record only the final column.
      auto d(d_ic);
      d.dt /= BaselineConsts::nstep;
      FortranBridge<Scalar> bridge(d);
      for (Int step = 0; step < BaselineConsts::nstep; ++step) {
        // The original Fortran serves as the baseline.
        micro_sed_func(d, bridge);
        // Record its output.
        util::write(&d.ni, 1, fid);
        util::write(&d.nk, 1, fid);
        util::write(&d.dt, 1, fid);
        const auto n = d.ni*d.nk;
        util::write(d.qr, n, fid);
        util::write(d.nr, n, fid);
        util::write(d.th, n, fid);
        util::write(d.dzq, n, fid);
        util::write(d.pres, n, fid);
        util::write(d.prt_liq, d.ni, fid);
      }
    }
  };
  return Observer(bfn).run(1);
}

template <typename Scalar>
static Int compare (const std::string& label, const Scalar* a,
                    const Scalar* b, const Int& n, const Real& tol) {
  Int nerr = 0;
  Real den = 0;
  for (Int i = 0; i < n; ++i)
    den = std::max(den, std::abs(a[i]));
  Real worst = 0;
  for (Int i = 0; i < n; ++i) {
    const auto num = std::abs(a[i] - b[i]);
    if (num > tol*den) {
      ++nerr;
#if 0
      std::cout << label << " bad idx: " << i << std::fixed << std::setprecision(12)
                << std::setw(20) << a[i] << " " << b[i] << std::endl;
#endif
      worst = std::max(worst, num);
    }
  }
  if (nerr)
    std::cout << label << " nerr " << nerr << " worst " << (worst/den)
              << " with denominator " << den << "\n";
  return nerr;
}

template <typename Scalar>
static Int compare (const std::string& label, const ic::MicroSedData<Scalar>& d_ref,
                    const ic::MicroSedData<Scalar>& d, const Real& tol) {
  assert(d_ref.ni == 1 && d.nk == d_ref.nk);
  // Compare just the last column.
  const auto os = d.nk*(d.ni-1);
  const auto n = d.nk;
  return (compare(label + " qr", d_ref.qr, d.qr + os, n, tol) +
          compare(label + " nr", d_ref.nr, d.nr + os, n, tol) +
          compare(label + " prt_liq", d_ref.prt_liq, d.prt_liq + (d.ni-1), d_ref.ni, tol) +
          // The rest should not be written, so check that they are BFB.
          compare(label + " th", d_ref.th, d.th + os, n, 0) +
          compare(label + " dzq", d_ref.dzq, d.dzq + os, n, 0) +
          compare(label + " pres", d_ref.pres, d.pres + os, n, 0));
}

template <typename Scalar>
static Int run_and_cmp (const std::string& bfn, const Real& tol) {
  struct Observer : public BaselineObserver {
    util::FILEPtr fid;
    const Real tol;

    Observer (const std::string& bfn, const Real& tol_)
      : tol(tol_)
    {
      fid = util::FILEPtr(fopen(bfn.c_str(), "r"));
      micro_throw_if( ! fid, "run_and_cmp can't read " << bfn);
    }

    virtual void observe (const ic::MicroSedData<Scalar>& d_ic) override {
      auto d_orig_fortran(d_ic);
      d_orig_fortran.dt /= BaselineConsts::nstep;
      FortranBridge<Scalar> f_bridge(d_orig_fortran);

      auto d_vanilla_cpp(d_orig_fortran);
      VanillaCppBridge<Scalar> vcpp_bridge(d_vanilla_cpp);

      auto d_kokkos_cpp(d_orig_fortran);
      KokkosCppBridge<Scalar> kcpp_bridge(d_kokkos_cpp);

      p3::micro_sed_vanilla::MicroSedFuncVanillaKokkos<Scalar> msvk(d_ic.ni, d_ic.nk);

      for (Int step = 0; step < BaselineConsts::nstep; ++step) {
        // Read the baseline.
        Int ni, nk;
        util::read(&ni, 1, fid);
        util::read(&nk, 1, fid);
        micro_throw_if(ni != 1 || nk != d_ic.nk,
                       "Baseline file has (ni,nk) = " << ni << ", " << nk
                       << " but we expected " << 1 << ", " << d_ic.nk);
        const auto n = ni*nk;
        ic::MicroSedData<Scalar> d_ref(ni, nk);
        util::read(&d_ref.dt, 1, fid);
        util::read(d_ref.qr, n, fid);
        util::read(d_ref.nr, n, fid);
        util::read(d_ref.th, n, fid);
        util::read(d_ref.dzq, n, fid);
        util::read(d_ref.pres, n, fid);
        util::read(d_ref.prt_liq, d_ref.ni, fid);

        // Run various models and compare to baseline.

        { // Compare the Fortran code in case we need to change it, as we will
          // for handling the issue of single vs double precision. In this case,
          // we expect the baseline file to be generated from the master-branch
          // version, and we're comparing a modified version in a branch.
          micro_sed_func(d_orig_fortran, f_bridge);
          std::stringstream ss;
          ss << "Original Fortran step " << step;
          nerr += compare(ss.str(), d_ref, d_orig_fortran, tol);
        }

        { // Super-vanilla C++.
          micro_sed_func_cpp(d_vanilla_cpp, vcpp_bridge);
          std::stringstream ss;
          ss << "Super-vanilla C++ step " << step;
          nerr += compare(ss.str(), d_ref, d_vanilla_cpp,
                          util::is_single_precision<Real>::value ? 2e-5 : tol);
        }

        { // Vanilla C++ kokkos.
          micro_sed_func_cpp_kokkos(d_kokkos_cpp, kcpp_bridge, msvk);
          std::stringstream ss;
          ss << "Vanilla Kokkos C++ step " << step;
          nerr += compare(ss.str(), d_ref, d_kokkos_cpp,
                          util::is_single_precision<Real>::value ? 2e-5 : tol);
        }
      }
    }
  };

  Int nerr = 0;
  Int ne = Observer(bfn, tol).run(1);
  if (ne) std::cout << "1-column test failed.\n";
  nerr += ne;
  ne = Observer(bfn, tol).run(7);
  if (ne) std::cout << "Multiple-column test failed.\n";
  nerr += ne;
  return nerr;
}

static void expect_another_arg (Int i, Int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}

static Int unittest () {
  Int nerr = 0;
  if (util::is_single_precision<double>::value) nerr++;
  if ( ! util::is_single_precision<float>::value) nerr++;
  return nerr;
}

int main (int argc, char** argv) {
#ifdef FPE
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
                         ~( _MM_MASK_INVALID |
                            _MM_MASK_DIV_ZERO |
                            _MM_MASK_OVERFLOW |
                            _MM_MASK_UNDERFLOW ));
#endif

  if (argc == 1) {
    std::cout <<
      argv[0] << " [options] baseline-filename\n"
      "Options:\n"
      "  -g        Generate baseline file.\n"
      "  -t <tol>  Tolerance for relative error.\n";
    return -1;
  }

  bool generate = false;
  Real tol = 0;
  for (Int i = 1; i < argc-1; ++i) {
    if (util::eq(argv[i], "-g", "--generate")) generate = true;
    if (util::eq(argv[i], "-t", "--tol")) {
      expect_another_arg(i, argc);
      ++i;
      tol = std::atof(argv[i]);
    }
  }

  // Always decorate baseline name with precision info
  std::string baseline_fn(argv[argc-1]);
  baseline_fn += std::to_string(sizeof(Real));

  Int out = 0;
  Kokkos::initialize(argc, argv); {
    out = unittest();
    if (generate) {
      // Generate a single-column baseline file.
      std::cout << "Generating to " << baseline_fn << "\n";
      out += generate_baseline<Real>(baseline_fn);
    } else {
      // Run with multiple columns, but compare only the last one to the baseline.
      printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
      out += run_and_cmp<Real>(baseline_fn, tol);
    }
  } Kokkos::finalize_all();

  return out;
}
