#include "types.hpp"
#include "util.hpp"

#include <vector>
#include <iostream>
#include <exception>

// Nothing in this file is intended to be indicative of the eventual perf
// portable impl. E.g., Kokkos is not used here. The objective is simply to test
// for input/output regression errors.

extern "C" void
micro_sed_func_c(
  Int kts, Int kte, Int ni, Int nk, Int its, Int ite, Real dt,
  Real* qr, Real* nr, const Real* th, const Real* dzq, const Real* pres,
  Real* prt_liq);

template <typename Scalar>
struct MicroSedData {
  const Int ni, nk, nit;

private:
  std::vector<Scalar> buf_;

public:
  Real dt;
  Real* const qr;
  Real* const nr;
  Real* const th;
  Real* const dzq;
  Real* const pres;
  Real* const prt_liq;

  MicroSedData(Int ni_, Int nk_, Int nit_)
    : ni(ni_), nk(nk_), nit(nit_),
      buf_(6*nit*nk),
      dt(0),
      qr(buf_.data()),
      nr(qr + nit*nk),
      th(nr + nit*nk),
      dzq(th + nit*nk),
      pres(dzq + nit*nk),
      prt_liq(pres + nit*nk)
  {}
};

template <typename Scalar>
void micro_sed_func (MicroSedData<Scalar>& d) {
  micro_sed_func_c(1, d.nk, d.ni, d.nk, 1, d.nit, d.dt, d.qr, d.nr, d.th, d.dzq,
                   d.pres, d.prt_liq);
}

template <typename Scalar>
struct MicroSedObserver {
  virtual ~MicroSedObserver () {}
  virtual void observe (const MicroSedData<Scalar>& d) = 0;
};

template <typename Scalar>
void run_over_parameter_sets (MicroSedObserver<Scalar>& o) {
  MicroSedData<Scalar> d(10,10,10);
  o.observe(d);
}

struct BaselineObserver : public MicroSedObserver<Real> {
  Int nerr;

  BaselineObserver () : nerr(0) {}

  Int run () {
    run_over_parameter_sets(*this);
    return nerr;
  }
};

static Int generate_baseline (const std::string& bfn) {
  struct Observer : public BaselineObserver {
    virtual void observe (const MicroSedData<Real>& d) override {
    }
  };
  return Observer().run();
}

static Int run_and_cmp (const std::string& bfn) {
  struct Observer : public BaselineObserver {
    virtual void observe (const MicroSedData<Real>& d) override {
    }
  };
  return Observer().run();
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
    out = generate_baseline(baseline_fn);
  } else {
    printf("Comparing with %s at tol %1.1e\n", baseline_fn.c_str(), tol);
    out = run_and_cmp(baseline_fn);
  }

  return out;
}
