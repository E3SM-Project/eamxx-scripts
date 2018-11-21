#ifndef P3_COMMON_HPP
#define P3_COMMON_HPP

#include "initial_conditions.hpp"
#include "types.hpp"
#include "scream_pack.hpp"
#include "kokkos_util.hpp"

#include <vector>

extern "C" {

// Things we get from fortran
void p3_init();
Real* c_get_vn_table();
Real* c_get_vm_table();
Real* c_get_mu_r_table();

} // extern C

namespace p3 {
namespace micro_sed {

using scream::pack::scalarize;
using scream::pack::BigPack;

template <typename Scalar>
struct Globals
{
  static constexpr Scalar INV_RHOW = 1.e-3;
  static constexpr Scalar RHOW     = 997.0;
  static constexpr Scalar THRD     = 1.0/3.0;
  static constexpr Scalar SXTH     = 1.0/6.0;
  static constexpr Scalar PI       = 3.14159265;
  static constexpr Scalar PIOV6    = PI*SXTH;
  static constexpr Scalar CONS1    = PIOV6*RHOW;
  static constexpr Scalar QSMALL   = 1.e-14;
  static constexpr Scalar NSMALL   = 1.e-16;
  static constexpr Scalar RD       = 287.15;
  static constexpr Scalar RHOSUR   = 100000.0/(RD*273.15);
  static constexpr Scalar CP       = 1005.0;
  static constexpr Scalar INV_CP   = 1.0/CP;

  static vector_2d_t<Scalar> VN_TABLE, VM_TABLE;
  static std::vector<Scalar> MU_R_TABLE;
};

template <typename Scalar>
vector_2d_t<Scalar> Globals<Scalar>::VN_TABLE;

template <typename Scalar>
vector_2d_t<Scalar> Globals<Scalar>::VM_TABLE;

template <typename Scalar>
std::vector<Scalar> Globals<Scalar>::MU_R_TABLE;

template <typename Scalar>
constexpr Scalar Globals<Scalar>::NSMALL;

template <typename Scalar>
void populate_input(const int nk, const int kdir,
                    std::vector<Scalar> & qr, std::vector<Scalar> & nr, std::vector<Scalar> & th, std::vector<Scalar> & dzq, std::vector<Scalar> & pres)
{
  ic::MicroSedData<Scalar> data(1, nk);
  populate(data, kdir);

  for (int k = 0; k < nk; ++k) {
    qr[k]   = data.qr[k];
    nr[k]   = data.nr[k];
    th[k]   = data.th[k];
    dzq[k]  = data.dzq[k];
    pres[k] = data.pres[k];
  }
}

/**
 * Generate lookup table for rain fallspeed and ventilation parameters
 * the lookup table is two dimensional as a function of number-weighted mean size
 * proportional to qr/Nr and shape parameter mu_r
 */
template <typename Scalar>
void p3_init_cpp()
{
  static bool is_init = false;
  if (is_init) {
    return;
  }
  is_init = true;

  Globals<Scalar>::VN_TABLE.resize(300, std::vector<Scalar>(10));
  Globals<Scalar>::VM_TABLE.resize(300, std::vector<Scalar>(10));
  Globals<Scalar>::MU_R_TABLE.resize(150);

  p3_init();

  Scalar* vn_table   = c_get_vn_table();
  Scalar* vm_table   = c_get_vm_table();
  Scalar* mu_r_table = c_get_mu_r_table();

  for (int i = 0; i < 300; ++i) {
    for (int k = 0; k < 10; ++k) {
      Globals<Scalar>::VN_TABLE[i][k] = vn_table[300*k + i];
      Globals<Scalar>::VM_TABLE[i][k] = vm_table[300*k + i];
    }
  }

  for (int i = 0; i < 150; ++i) {
    Globals<Scalar>::MU_R_TABLE[i] = mu_r_table[i];
  }
}

template <typename Scalar, typename D=DefaultDevice>
void dump_to_file_k(const char* basename,
                    const typename KokkosTypes<D>::template view_2d<Scalar>& qr,
                    const typename KokkosTypes<D>::template view_2d<Scalar>& nr,
                    const typename KokkosTypes<D>::template view_2d<Scalar>& th,
                    const typename KokkosTypes<D>::template view_2d<Scalar>& dzq,
                    const typename KokkosTypes<D>::template view_2d<Scalar>& pres,
                    const typename KokkosTypes<D>::template view_1d<Scalar>& prt_liq,
                    const int ni, const int nk, const Scalar dt, const int ts)
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

  const Scalar* qr_md   = reinterpret_cast<const Scalar*>(qr_m.data());
  const Scalar* nr_md   = reinterpret_cast<const Scalar*>(nr_m.data());
  const Scalar* th_md   = reinterpret_cast<const Scalar*>(th_m.data());
  const Scalar* dzq_md  = reinterpret_cast<const Scalar*>(dzq_m.data());
  const Scalar* pres_md = reinterpret_cast<const Scalar*>(pres_m.data());

  util::dump_to_file(
    basename, qr_md, nr_md, th_md, dzq_md, pres_md, prt_liq_m.data(),
    ni, nk, dt, ts);
}

template <typename Scalar, typename D=DefaultDevice>
void dump_to_file_k (
  const char* basename,
  const typename KokkosTypes<D>::template view_2d<BigPack<Scalar> >& qr,
  const typename KokkosTypes<D>::template view_2d<BigPack<Scalar> >& nr,
  const typename KokkosTypes<D>::template view_2d<BigPack<Scalar> >& th,
  const typename KokkosTypes<D>::template view_2d<BigPack<Scalar> >& dzq,
  const typename KokkosTypes<D>::template view_2d<BigPack<Scalar> >& pres,
  const typename KokkosTypes<D>::template view_1d<Scalar>& prt_liq,
  const int ni, const int nk, const Scalar dt, const int ts)
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

  const auto
    sqr = scalarize(qr_m),
    snr = scalarize(nr_m),
    sth = scalarize(th_m),
    sdzq = scalarize(dzq_m),
    spres = scalarize(pres_m);
  const int ldk = sqr.extent_int(1);

  const Scalar* qr_md   = reinterpret_cast<const Scalar*>(sqr.data());
  const Scalar* nr_md   = reinterpret_cast<const Scalar*>(snr.data());
  const Scalar* th_md   = reinterpret_cast<const Scalar*>(sth.data());
  const Scalar* dzq_md  = reinterpret_cast<const Scalar*>(sdzq.data());
  const Scalar* pres_md = reinterpret_cast<const Scalar*>(spres.data());

  util::dump_to_file(
    basename, qr_md, nr_md, th_md, dzq_md, pres_md, prt_liq_m.data(),
    ni, nk, dt, ts, ldk);
}

template <typename T, typename Scalar, typename D=DefaultDevice>
void populate_kokkos_from_vec(const int num_vert, std::vector<Scalar> const& vec,
                              typename KokkosTypes<D>::template view_1d<T>& device)
{
  const auto mirror = Kokkos::create_mirror_view(device);

  for (int k = 0; k < num_vert; ++k) {
    reinterpret_cast<Scalar*>(mirror.data())[k] = vec[k];
  }

  Kokkos::deep_copy(device, mirror);
}

template <typename Scalar, typename MSK>
void micro_sed_func_kokkos_wrap(const int ni, const int nk, const Scalar dt, const int ts, const int kdir, const int repeat)
{
  util::dump_arch();

  std::cout << "Running " << MSK::NAME <<  " with ni=" << ni << ", nk=" << nk
            << ", dt=" << dt << ", ts=" << ts << ", kdir=" << kdir << MSK::custom_msg() << std::endl;

  MSK msk(ni, nk);

  const int num_vert = msk.get_num_vert();
  typename MSK::template view_2d<typename MSK::Pack> qr("qr", ni, num_vert),
    nr("nr", ni, num_vert),
    th("th", ni, num_vert),
    dzq("dzq", ni, num_vert),
    pres("pres", ni, num_vert);

  typename MSK::template view_1d<typename MSK::Pack> qr_i("qr_i", num_vert),
    nr_i("nr_i", num_vert),
    th_i("th_i", num_vert),
    dzq_i("dzq_i", num_vert),
    pres_i("pres_i", num_vert);

  typename MSK::template view_1d<Scalar> prt_liq("prt_liq", ni), prt_liq_i("prt_liq", ni);

  {
    std::vector<Scalar> qr_v(nk), nr_v(nk), th_v(nk), dzq_v(nk), pres_v(nk);

    populate_input(nk, kdir, qr_v, nr_v, th_v, dzq_v, pres_v);

    for (auto item : { std::make_pair(&qr_v, &qr_i), std::make_pair(&nr_v, &nr_i), std::make_pair(&th_v, &th_i),
          std::make_pair(&dzq_v, &dzq_i), std::make_pair(&pres_v, &pres_i)}) {
      populate_kokkos_from_vec<typename MSK::Pack, Scalar>(nk, *(item.first), *(item.second));
    }
  }

  // This time is thrown out, I just wanted to be able to use auto
  auto start = std::chrono::steady_clock::now();

  for (int r = 0; r < repeat+1; ++r) {
    Kokkos::parallel_for("Re-init",
                         util::ExeSpaceUtils<typename MSK::ExeSpace>::get_default_team_policy(ni, num_vert),
                         KOKKOS_LAMBDA(typename MSK::MemberType team_member) {
      const int i = team_member.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, num_vert), [=] (int k) {
        qr(i, k)   = qr_i(k);
        nr(i, k)   = nr_i(k);
        th(i, k)   = th_i(k);
        dzq(i, k)  = dzq_i(k);
        pres(i, k) = pres_i(k);
      });
      prt_liq(i) = prt_liq_i(i);
    });

    for (int i = 0; i < ts; ++i) {
      msk.micro_sed_func(kdir == 1 ? 1 : nk, kdir == 1 ? nk : 1,
                         1, ni, dt, qr, nr, th, dzq, pres, prt_liq);
      Kokkos::fence();
    }

    if (r == 0) {
      start = std::chrono::steady_clock::now();
    }
  }

  auto finish = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);

  const double report_time = (1e-6*duration.count()) / repeat;

  printf("Time = %1.3e seconds\n", report_time);

  dump_to_file_k(MSK::NAME, qr, nr, th, dzq, pres, prt_liq, ni, nk, dt, ts);
}

} // namespace p3
} // namespace micro_sed

#endif
