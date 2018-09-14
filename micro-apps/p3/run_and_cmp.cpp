#include "types.hpp"
#include "util.hpp"
#include "initial_conditions.hpp"
#include "micro_sed_vanilla_kokkos.hpp"
#include "micro_sed_workspace_kokkos.hpp"
#include "micro_sed_pack_kokkos.hpp"
#include "micro_kokkos.hpp"
#include "cmp.hpp"

#include <vector>
#include <iostream>
#include <exception>

// Nothing in this file is intended to be indicative of the eventual perf
// portable impl. E.g., Kokkos is not used here. The objective is simply to test
// for input/output regression errors.

extern "C" {
  void p3_init();
  void micro_sed_func_c(
    Int kts, Int kte, Int kdir, Int its, Int ite, Real dt,
    Real* qr, Real* nr, const Real* th, const Real* dzq, const Real* pres,
    Real* prt_liq);
}


template <cmp::TransposeDirection::Enum direction, typename Scalar>
void transpose_layout (const ic::MicroSedData<Scalar>& s,
                        ic::MicroSedData<Scalar>& d) {
  cmp::transpose<direction>(s.qr, d.qr, d.ni, d.nk);
  cmp::transpose<direction>(s.nr, d.nr, d.ni, d.nk);
  cmp::transpose<direction>(s.th, d.th, d.ni, d.nk);
  cmp::transpose<direction>(s.dzq, d.dzq, d.ni, d.nk);
  cmp::transpose<direction>(s.pres, d.pres, d.ni, d.nk);
  for (Int i = 0; i < s.ni; ++i) d.prt_liq[i] = s.prt_liq[i];
}

// MicroSedData <-> micro_sed.f90 format.
template <typename Scalar>
struct FortranBridge : public ic::MicroSedData<Scalar> {
  FortranBridge (const ic::MicroSedData<Scalar>& d)
    : ic::MicroSedData<Scalar>(d)
  { transpose_layout<cmp::TransposeDirection::c2f>(d, *this); }

  void sync_to (ic::MicroSedData<Scalar>& d) {
    transpose_layout<cmp::TransposeDirection::f2c>(*this, d);
  }
};

template <typename Scalar>
void micro_sed_func (ic::MicroSedData<Scalar>& d, FortranBridge<Scalar>& b) {
  micro_sed_func_c(1, d.nk,
                   d.reverse ? -1 : 1,
                   1, d.ni, d.dt, b.qr, b.nr, b.th,
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
    Kokkos::deep_copy(prt_liq, mirror_prt_liq);
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
    Kokkos::deep_copy(mirror_prt_liq, prt_liq);

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
class KokkosPackBridge {
public:
  using RealPack = p3::micro_sed_pack::RealPack;

  int np;
  kokkos_2d_t<RealPack> qr, nr, th, dzq, pres;
  kokkos_1d_t<Real> prt_liq;

  KokkosPackBridge(const ic::MicroSedData<Scalar>& d)
    : np(scream::pack::npack<RealPack>(d.nk)),
      qr("qr", d.ni, np),
      nr("nr", d.ni, np),
      th("th", d.ni, np),
      dzq("dzq", d.ni, np),
      pres("pres", d.ni, np),
      prt_liq("prt_liq", d.ni)
  {
    sync_from(d);
  }

  void sync_from(const ic::MicroSedData<Scalar>& d)
  {
    using p3::micro_sed_pack::scalarize;

    const auto
      qr_m = Kokkos::create_mirror_view(qr),
      nr_m = Kokkos::create_mirror_view(nr),
      th_m = Kokkos::create_mirror_view(th),
      dzq_m = Kokkos::create_mirror_view(dzq),
      pres_m = Kokkos::create_mirror_view(pres);
    const auto prt_liq_m = Kokkos::create_mirror_view(prt_liq);

    const auto
      sqr = scalarize(qr_m),
      snr = scalarize(nr_m),
      sth = scalarize(th_m),
      sdzq = scalarize(dzq_m),
      spres = scalarize(pres_m);

    for (int i = 0; i < d.ni; ++i) {
      for (int k = 0; k < d.nk; ++k) {
        sqr  (i, k) = d.qr  [i*d.nk + k];
        snr  (i, k) = d.nr  [i*d.nk + k];
        sth  (i, k) = d.th  [i*d.nk + k];
        sdzq (i, k) = d.dzq [i*d.nk + k];
        spres(i, k) = d.pres[i*d.nk + k];
      }
      prt_liq_m(i) = d.prt_liq[i];
    }

    Kokkos::deep_copy(qr, qr_m);
    Kokkos::deep_copy(nr, nr_m);
    Kokkos::deep_copy(th, th_m);
    Kokkos::deep_copy(dzq, dzq_m);
    Kokkos::deep_copy(pres, pres_m);
    Kokkos::deep_copy(prt_liq, prt_liq_m);
  }

  void sync_to(ic::MicroSedData<Scalar>& d) const
  {
    using p3::micro_sed_pack::scalarize;

    const auto
      qr_m = Kokkos::create_mirror_view(qr),
      nr_m = Kokkos::create_mirror_view(nr),
      th_m = Kokkos::create_mirror_view(th),
      dzq_m = Kokkos::create_mirror_view(dzq),
      pres_m = Kokkos::create_mirror_view(pres);
    const auto prt_liq_m = Kokkos::create_mirror_view(prt_liq);

    Kokkos::deep_copy(qr_m, qr);
    Kokkos::deep_copy(nr_m, nr);
    Kokkos::deep_copy(th_m, th);
    Kokkos::deep_copy(dzq_m, dzq);
    Kokkos::deep_copy(pres_m, pres);
    Kokkos::deep_copy(prt_liq_m, prt_liq);

    const auto
      sqr = scalarize(qr_m),
      snr = scalarize(nr_m),
      sth = scalarize(th_m),
      sdzq = scalarize(dzq_m),
      spres = scalarize(pres_m);

    for (int i = 0; i < d.ni; ++i) {
      for (int k = 0; k < d.nk; ++k) {
        d.qr  [i*d.nk + k] = sqr  (i, k);
        d.nr  [i*d.nk + k] = snr  (i, k);
        d.th  [i*d.nk + k] = sth  (i, k);
        d.dzq [i*d.nk + k] = sdzq (i, k);
        d.pres[i*d.nk + k] = spres(i, k);
      }
      d.prt_liq[i] = prt_liq_m(i);
    }
  }
};

template <typename Scalar>
void micro_sed_func_cpp (ic::MicroSedData<Scalar>& d, VanillaCppBridge<Scalar>& bridge)
{
  p3::micro_sed::micro_sed_func<Scalar>( d.reverse ? d.nk : 1, d.reverse ? 1 : d.nk,
                                         1, d.ni, d.dt, bridge.qr, bridge.nr, bridge.th,
                                         bridge.dzq, bridge.pres, bridge.prt_liq);

  bridge.sync_to(d);
}

template <typename Scalar, typename MSK>
void micro_sed_func_cpp_kokkos (ic::MicroSedData<Scalar>& d, KokkosCppBridge<Scalar>& bridge, MSK& msk)
{
  p3::micro_sed::micro_sed_func( msk, d.reverse ? d.nk : 1, d.reverse ? 1 : d.nk,
                                 1, d.ni, d.dt, bridge.qr, bridge.nr, bridge.th,
                                 bridge.dzq, bridge.pres, bridge.prt_liq);

  bridge.sync_to(d);
}

template <typename Scalar>
void micro_sed_func_cpp_pack_kokkos (ic::MicroSedData<Scalar>& d, KokkosPackBridge<Scalar>& bridge, p3::micro_sed_pack::MicroSedFuncPackKokkos& msvk)
{
  micro_sed_func_pack_kokkos( msvk, d.reverse ? d.nk : 1, d.reverse ? 1 : d.nk,
                              1, d.ni, d.dt, bridge.qr, bridge.nr, bridge.th,
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
  p3::micro_sed::p3_init_cpp<Scalar>();

  ic::MicroSedData<Scalar> d(ncol, 111);
  d.dt = dt_tot;
  populate(d, 1);

  o.observe(d);

  const auto d_rev(reverse_k(d));
  o.observe(d_rev);
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
static Int compare (const std::string& label, const ic::MicroSedData<Scalar>& d_ref,
                    const ic::MicroSedData<Scalar>& d, const Real& tol, bool verbose) {
  assert(d_ref.ni == 1 && d.nk == d_ref.nk);
  // Compare just the last column.
  const auto os = d.nk*(d.ni-1);
  const auto n = d.nk;
  return (cmp::compare(label + " qr", d_ref.qr, d.qr + os, n, tol, verbose) +
          cmp::compare(label + " nr", d_ref.nr, d.nr + os, n, tol, verbose) +
          cmp::compare(label + " prt_liq", d_ref.prt_liq, d.prt_liq + (d.ni-1), d_ref.ni, tol, verbose) +
          // The rest should not be written, so check that they are BFB.
          cmp::compare(label + " th", d_ref.th, d.th + os, n, 0, verbose) +
          cmp::compare(label + " dzq", d_ref.dzq, d.dzq + os, n, 0, verbose) +
          cmp::compare(label + " pres", d_ref.pres, d.pres + os, n, 0, verbose));
}

template <typename Bridge, typename Lambda, typename Scalar>
int do_compare(const Lambda& func, ic::MicroSedData<Scalar>& d, const ic::MicroSedData<Scalar>& d_ref, const Real& tol,
               const char* label, int step, bool verbose)
{
  Bridge bridge(d);

  func(d, bridge);

  std::stringstream ss;
  ss << label << " step " << step;

  return compare(ss.str(), d_ref, d, tol, verbose);
}

template <typename Scalar>
static Int run_and_cmp (const std::string& bfn, const Real& tol, bool verbose) {
  struct Observer : public BaselineObserver {
    util::FILEPtr fid;
    const Real tol;
    bool verbose;

    Observer (const std::string& bfn, const Real& tol_, bool verbose_)
      : tol(tol_), verbose(verbose_)
    {
      fid = util::FILEPtr(fopen(bfn.c_str(), "r"));
      micro_throw_if( ! fid, "run_and_cmp can't read " << bfn);

      // Sanity check.
      micro_throw_if( ! util::is_single_precision<Real>::value && tol != 0,
                      "We want BFB in double precision, at least in DEBUG builds.");
    }

    virtual void observe (const ic::MicroSedData<Scalar>& d_ic) override {
      auto d_ic_cp(d_ic);
      d_ic_cp.dt /= BaselineConsts::nstep;

      std::vector<ic::MicroSedData<Scalar> > ds = {ic::MicroSedData<Scalar>(d_ic_cp), ic::MicroSedData<Scalar>(d_ic_cp), ic::MicroSedData<Scalar>(d_ic_cp), ic::MicroSedData<Scalar>(d_ic_cp)};

      p3::micro_sed::MicroSedFuncWorkspaceKokkos<Scalar> mswk(d_ic.ni, d_ic.nk);
      p3::micro_sed::MicroSedFuncVanillaKokkos<Scalar> msvk(d_ic.ni, d_ic.nk);

      p3::micro_sed_pack::MicroSedFuncPackKokkos mspk(d_ic.ni, d_ic.nk);
      auto d_pack_kokkos_cpp(d_ic_cp);
      KokkosPackBridge<Scalar> kcpp_pack_bridge(d_pack_kokkos_cpp);

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

        const Real sptol = 2e-5;
        Real cpp_tol = (util::is_single_precision<Real>::value && tol < sptol) ? sptol : tol;

        // Compare the Fortran code in case we need to change it, as we will
        // for handling the issue of single vs double precision. In this case,
        // we expect the baseline file to be generated from the master-branch
        // version, and we're comparing a modified version in a branch.
        nerr += do_compare<FortranBridge<Scalar> >([] (ic::MicroSedData<Scalar>& d, FortranBridge<Scalar>& b) { micro_sed_func(d, b); },
                                                   ds[0], d_ref, tol, "Original Fortran", step, verbose);

        nerr += do_compare<VanillaCppBridge<Scalar> >([] (ic::MicroSedData<Scalar>& d, VanillaCppBridge<Scalar>& b) { micro_sed_func_cpp(d, b); },
                                                      ds[1], d_ref, cpp_tol, "Super-vanilla C++", step, verbose);

        nerr += do_compare<KokkosCppBridge<Scalar> >([&] (ic::MicroSedData<Scalar>& d, KokkosCppBridge<Scalar>& b) { micro_sed_func_cpp_kokkos(d, b, msvk); },
                                                     ds[2], d_ref, cpp_tol, "Vanilla Kokkos C++", step, verbose);

        nerr += do_compare<KokkosCppBridge<Scalar> >([&] (ic::MicroSedData<Scalar>& d, KokkosCppBridge<Scalar>& b) { micro_sed_func_cpp_kokkos(d, b, mswk); },
                                                     ds[3], d_ref, cpp_tol, "Workspace Kokkos C++", step, verbose);

        { // C++ kokkos with packs.
          micro_sed_func_cpp_pack_kokkos(d_pack_kokkos_cpp, kcpp_pack_bridge, mspk);
          std::stringstream ss;
          ss << "Pack Kokkos C++ step " << step;
          const Real sptol = 2e-5;
          Real kokkos_tol = (util::is_single_precision<Real>::value && tol < sptol) ? sptol : tol;
          // Sanity check.
          micro_throw_if( ! util::is_single_precision<Real>::value && tol != 0,
                          "We want BFB in double precision, at least in DEBUG builds.");
          nerr += compare(ss.str(), d_ref, d_pack_kokkos_cpp, kokkos_tol, verbose);
        }
      }
    }
  };

  Int nerr = 0;
  Int ne = Observer(bfn, tol, verbose).run(1);
  if (ne) std::cout << "1-column test failed.\n";
  nerr += ne;
  ne = Observer(bfn, tol, verbose).run(7);
  if (ne) std::cout << "Multiple-column test failed.\n";
  nerr += ne;
  return nerr;
}

static void expect_another_arg (Int i, Int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}

int main (int argc, char** argv) {
  util::initialize();

  if (argc == 1) {
    std::cout <<
      argv[0] << " [options] baseline-filename\n"
      "Options:\n"
      "  -g        Generate baseline file.\n"
      "  -v        Run with extra verbose output.\n"
      "  -t <tol>  Tolerance for relative error.\n";
    return -1;
  }

  bool generate = false, verbose=false;
  Real tol = 0;
  for (Int i = 1; i < argc-1; ++i) {
    if (util::eq(argv[i], "-g", "--generate")) generate = true;
    if (util::eq(argv[i], "-v", "--verbose")) verbose = true;
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
    if (generate) {
      // Generate a single-column baseline file.
      std::cout << "Generating to " << baseline_fn << "\n";
      out += generate_baseline<Real>(baseline_fn);
    } else {
      // Run with multiple columns, but compare only the last one to the baseline.
      printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
      out += run_and_cmp<Real>(baseline_fn, tol, verbose);
    }
  } Kokkos::finalize();

  return out;
}
