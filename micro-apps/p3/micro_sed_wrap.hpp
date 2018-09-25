#ifndef INCLUDE_MICRO_SED_WRAP
#define INCLUDE_MICRO_SED_WRAP

#include "util.hpp"
#include "micro_kokkos.hpp"

#include <vector>
#include <chrono>

namespace p3 {
namespace micro_sed {

template <typename Real>
void dump_to_file_k(const char* basename,
                    const kokkos_2d_t<Real>& qr, const kokkos_2d_t<Real>& nr, const kokkos_2d_t<Real>& th, const kokkos_2d_t<Real>& dzq,
                    const kokkos_2d_t<Real>& pres, const kokkos_1d_t<Real>& prt_liq, const int ni, const int nk, const Real dt, const int ts)
{
  auto qr_m      = Kokkos::create_mirror_view(qr);
  auto nr_m      = Kokkos::create_mirror_view(nr);
  auto th_m      = Kokkos::create_mirror_view(th);
  auto dzq_m     = Kokkos::create_mirror_view(dzq);
  auto pres_m    = Kokkos::create_mirror_view(pres);
  auto prt_liq_m = Kokkos::create_mirror_view(prt_liq);

  Kokkos::deep_copy(qr_m,      qr);
  Kokkos::deep_copy(nr_m,      nr);
  Kokkos::deep_copy(th_m,      th);
  Kokkos::deep_copy(dzq_m,     dzq);
  Kokkos::deep_copy(pres_m,    pres);
  Kokkos::deep_copy(prt_liq_m, prt_liq);

  util::dump_to_file(basename, qr_m.data(), nr_m.data(), th_m.data(), dzq_m.data(), pres_m.data(), prt_liq_m.data(), ni, nk, dt, ts);
}

template <typename Real>
void dump_to_file_v(const vector_2d_t<Real>& qr, const vector_2d_t<Real>& nr, const vector_2d_t<Real>& th, const vector_2d_t<Real>& dzq,
                    const vector_2d_t<Real>& pres, const std::vector<Real>& prt_liq, const Real dt, const int ts)
{
  const int ni = qr.size();
  const int nk = qr[0].size();

  std::vector<Real> qr_1d(ni*nk), nr_1d(ni*nk), th_1d(ni*nk), dzq_1d(ni*nk), pres_1d(ni*nk);

  for (int i = 0; i < ni; ++i) {
    for (auto item : { std::make_pair(&qr, &qr_1d), std::make_pair(&nr, &nr_1d), std::make_pair(&th, &th_1d),
          std::make_pair(&dzq, &dzq_1d), std::make_pair(&pres, &pres_1d) }) {
      const Real* data = (*item.first)[i].data();
      std::copy( data, data + nk, item.second->data() + i*nk);
    }
  }

  util::dump_to_file("vanilla", qr_1d.data(), nr_1d.data(), th_1d.data(), dzq_1d.data(), pres_1d.data(), prt_liq.data(), ni, nk, dt, ts);
}

template <typename Real>
void populate_kokkos_from_vec(const int num_horz, const int num_vert, vector_2d_t<Real> const& vec, kokkos_2d_t<Real>& device)
{
  typename kokkos_2d_t<Real>::HostMirror mirror = Kokkos::create_mirror_view(device);

  for (int i = 0; i < num_horz; ++i) {
    for (int k = 0; k < num_vert; ++k) {
      mirror(i, k) = vec[i][k];
    }
  }

  Kokkos::deep_copy(device, mirror);
}

template <typename Real>
void micro_sed_func_vanilla_wrap(const int ni, const int nk, const Real dt, const int ts, const int kdir, const int repeat)
{
  vector_2d_t<Real> qr(ni,    std::vector<Real>(nk)),
                    nr(ni,    std::vector<Real>(nk)),
                    th(ni,    std::vector<Real>(nk)),
                    dzq(ni,   std::vector<Real>(nk)),
                    pres(ni,  std::vector<Real>(nk));

  std::vector<Real> prt_liq(ni);

  util::dump_arch();
  std::cout << "Running micro_sed_vanilla with ni=" << ni << ", nk=" << nk
            << ", dt=" << dt << ", ts=" << ts << ", kdir=" << kdir << std::endl;

  populate_input(ni, nk, kdir, qr, nr, th, dzq, pres);

  auto start = std::chrono::steady_clock::now();

  for (int i = 0; i < ts; ++i) {
    micro_sed_func<Real>(kdir == 1 ? 1 : nk, kdir == 1 ? nk : 1,
                         1, ni, dt, qr, nr, th, dzq, pres, prt_liq);
  }

  auto finish = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
  printf("Time = %1.3e seconds\n", 1e-6*duration.count());

  dump_to_file_v(qr, nr, th, dzq, pres, prt_liq, dt, ts);
}

template <typename Real, typename MSK>
void micro_sed_func_kokkos_wrap(const int ni, const int nk, const Real dt, const int ts, const int kdir, const int repeat)
{
  util::dump_arch();

  std::cout << "Running " << MSK::NAME <<  " with ni=" << ni << ", nk=" << nk
            << ", dt=" << dt << ", ts=" << ts << ", kdir=" << kdir << MSK::custom_msg() << std::endl;

  MSK msk(ni, nk);

  typename MSK::msk_2d_kokkos_t qr("qr", ni, msk.get_num_vert()),
    nr("nr", ni, msk.get_num_vert()),
    th("th", ni, msk.get_num_vert()),
    dzq("dzq", ni, msk.get_num_vert()),
    pres("pres", ni, msk.get_num_vert());

  kokkos_1d_t<Real> prt_liq("prt_liq", ni);

  {
    std::vector<Real> v(nk);
    vector_2d_t<Real> qr_v(ni, v), nr_v(ni, v), th_v(ni, v), dzq_v(ni, v), pres_v(ni, v);

    populate_input(ni, nk, kdir, qr_v, nr_v, th_v, dzq_v, pres_v);

    for (auto item : { std::make_pair(&qr_v, &qr), std::make_pair(&nr_v, &nr), std::make_pair(&th_v, &th),
          std::make_pair(&dzq_v, &dzq), std::make_pair(&pres_v, &pres)}) {
      vector_2d_t<Real>& arg1 = *(item.first);
      typename MSK::msk_2d_kokkos_t& arg2 = *(item.second);
      populate_kokkos_from_vec(ni, nk, arg1, arg2);
    }
  }

  auto start = std::chrono::steady_clock::now();

  for (int i = 0; i < ts; ++i) {
    micro_sed_func(msk,
                   kdir == 1 ? 1 : nk, kdir == 1 ? nk : 1,
                   1, ni, dt, qr, nr, th, dzq, pres, prt_liq);
    Kokkos::fence();
  }

  auto finish = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
  printf("Time = %1.3e seconds\n", 1e-6*duration.count());

  dump_to_file_k(MSK::NAME, qr, nr, th, dzq, pres, prt_liq, ni, nk, dt, ts);
}

} // namespace micro_sed
} // namespace p3

#endif
