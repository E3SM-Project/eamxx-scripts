#include "types.hpp"
#include "util.hpp"
#include "initial_conditions.hpp"

#include <vector>
#include <iostream>
#include <exception>

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

template <typename Scalar>
struct MicroSedData {
  const Int ni, nk;

private:
  std::vector<Scalar> buf_;

  void init_ptrs () {
    qr = buf_.data();
    nr = qr + ni*nk;
    th = nr + ni*nk;
    dzq = th + ni*nk;
    pres = dzq + ni*nk;
    prt_liq = pres + ni;
  }

public:
  Real dt;
  bool reverse;
  Real* qr;
  Real* nr;
  Real* th;
  Real* dzq;  // Not sure what 'q' means here, but this is dz [m].
  Real* pres;
  Real* prt_liq;

  MicroSedData(Int ni_, Int nk_)
    : ni(ni_), nk(nk_),
      buf_(5*ni*nk + ni),
      dt(0), reverse(false)
  {
    // For now, don't deal with the issue of fast index.
    micro_throw_if(ni > 1, "We're not yet handling chunks.");
    init_ptrs();
  }

  MicroSedData (const MicroSedData& d)
    : ni(d.ni), nk(d.nk), dt(d.dt), reverse(d.reverse), buf_(d.buf_)
  { init_ptrs(); }
};

template <typename Scalar>
static void duplicate_columns (MicroSedData<Scalar>& d) {
  const auto copy = [&] (Scalar* v, const Int& i) {
    std::copy(v, v + d.nk, v + i*d.nk);
  };
  for (Int i = 1; i < d.ni; ++i) {
    copy(d.qr, i);
    copy(d.nr, i);
    copy(d.th, i);
    copy(d.dzq, i);
    copy(d.pres, i);
    d.prt_liq[i] = d.prt_liq[0];
  }
}

template <typename Scalar>
static MicroSedData<Scalar> reverse_k (const MicroSedData<Scalar>& msd) {
  const auto reverse = [&] (const Scalar* s, Scalar* d) {
    for (Int i = 0; i < msd.ni; ++i) {
      const Scalar* si = s + msd.nk*i;
      Scalar* di = d + msd.nk*i;
      for (Int k = 0; k < msd.nk; ++k)
        di[k] = si[msd.nk - k - 1];
    }
  };
  MicroSedData<Scalar> r(msd);
  r.reverse = ! msd.reverse;
  reverse(msd.qr, r.qr);
  reverse(msd.nr, r.nr);
  reverse(msd.th, r.th);
  reverse(msd.dzq, r.dzq);
  reverse(msd.pres, r.pres);
  return r;
}

template <typename Scalar>
void micro_sed_func (MicroSedData<Scalar>& d) {
  micro_sed_func_c(1, d.nk,
                   d.reverse ? -1 : 1,
                   d.ni, d.nk, 1, d.ni, d.dt, d.qr, d.nr, d.th,
                   d.dzq, d.pres, d.prt_liq);
}

template <typename Scalar>
struct MicroSedObserver {
  virtual ~MicroSedObserver () {}
  virtual void observe (const MicroSedData<Scalar>& d) = 0;
};

template <typename Scalar>
void run_over_parameter_sets (MicroSedObserver<Scalar>& o) {
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
  const Real dt_tot = 9000; // s

  p3_init();

  MicroSedData<Scalar> d(1, 111);
  d.dt = dt_tot;
  ic::set_hydrostatic(d.nk, d.dzq, d.pres);
  for (Int k = 0; k < d.nk; ++k)
    d.th[k] = ic::consts::th_ref;
  ic::set_rain(d.nk, d.dzq, d.qr, d.nr);
  for (Int i = 0; i < d.ni; ++i)
    d.prt_liq[i] = 0;
  duplicate_columns(d);

  o.observe(d);

  const auto d_rev(reverse_k(d));
  o.observe(d_rev);
}

struct BaselineObserver : public MicroSedObserver<Real> {
  Int nerr;

  BaselineObserver () : nerr(0) {}

  Int run () {
    run_over_parameter_sets(*this);
    return nerr;
  }
};

struct BaselineConsts {
  // Take nstep outer steps so we can record output throughout the d.dt run.
  static constexpr Int nstep = 30;
};

template <typename Scalar>
static Int generate_baseline (const std::string& bfn) {
  struct Observer : public BaselineObserver {
    util::FILEPtr fid;

    Observer (const std::string& bfn) {
      fid = util::FILEPtr(fopen(bfn.c_str(), "w"));
      micro_throw_if( ! fid, "generate_baseline can't write " << bfn);
    }

    virtual void observe (const MicroSedData<Scalar>& d_ic) override {
      auto d(d_ic);
      d.dt /= BaselineConsts::nstep;
      for (Int step = 0; step < BaselineConsts::nstep; ++step) {
        // The original Fortran serves as the baseline.
        micro_sed_func(d);
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
  return Observer(bfn).run();
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
      worst = std::max(worst, num);
    }
  }
  if (nerr)
    std::cout << label << " nerr " << nerr << " worst " << (worst/den)
              << " with denominator " << den << "\n";
  return nerr;
}

template <typename Scalar>
static Int compare (const std::string& label, const MicroSedData<Scalar>& d_ref,
                    const MicroSedData<Scalar>& d, const Real& tol) {
  assert(d.ni == d_ref.ni && d.nk == d_ref.nk);
  const auto n = d.ni*d.nk;
  return (compare(label + " qr", d_ref.qr, d.qr, n, tol) +
          compare(label + " nr", d_ref.nr, d.nr, n, tol) +
          compare(label + " prt_liq", d_ref.prt_liq, d.prt_liq, d.ni, tol) +
          // The rest should not be written, so check that they are BFB.
          compare(label + " th", d_ref.th, d.th, n, 0) +
          compare(label + " dzq", d_ref.dzq, d.dzq, n, 0) +
          compare(label + " pres", d_ref.pres, d.pres, n, 0));
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

    virtual void observe (const MicroSedData<Scalar>& d_ic) override {
      auto d_orig_fortran(d_ic);
      d_orig_fortran.dt /= BaselineConsts::nstep;
      auto d_vanilla_cpp(d_orig_fortran);
      for (Int step = 0; step < BaselineConsts::nstep; ++step) {
        // Read the baseline.
        Int ni, nk;
        util::read(&ni, 1, fid);
        util::read(&nk, 1, fid);
        micro_throw_if(ni != d_ic.ni || nk != d_ic.nk,
                       "Baseline file has (ni,nk) = " << ni << ", " << nk
                       << " but we expected " << d_ic.ni << ", " << d_ic.nk);
        const auto n = ni*nk;
        MicroSedData<Scalar> d_ref(ni, nk);
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
          micro_sed_func(d_orig_fortran);
          std::stringstream ss;
          ss << "Original Fortran step " << step;
          nerr += compare(ss.str(), d_ref, d_orig_fortran, tol);
        }

        if (false) { // Super-vanilla C++.
          std::stringstream ss;
          ss << "Super-vanilla C++ step " << step;
          nerr += compare(ss.str(), d_ref, d_vanilla_cpp, tol);          
        }
      }
    }
  };
  return Observer(bfn, tol).run();
}

static void expect_another_arg (Int i, Int argc) {
  if (i == argc-1)
    throw std::runtime_error("Expected another cmd-line arg.");
}

int main (int argc, char** argv) {
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

  const std::string baseline_fn(argv[argc-1]);

  Int out = 0;
  if (generate) {
    std::cout << "Generating to " << baseline_fn << "\n";
    out = generate_baseline<Real>(baseline_fn);
  } else {
    printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
    out = run_and_cmp<Real>(baseline_fn, tol);
  }

  return out;
}
