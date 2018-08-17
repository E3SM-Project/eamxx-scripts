#include "micro_kokkos.hpp"

using Layout = Kokkos::LayoutRight;

#ifdef KOKKOS_ENABLE_CUDA
//using MemSpace = Kokkos::HostSpace;
using MemSpace = Kokkos::CudaSpace;
#else
using MemSpace = Kokkos::HostSpace;
#endif

using ExecSpace = Kokkos::Cuda;

using kokkos_2d_t = Kokkos::View<double**, Layout, MemSpace>;

using kokkos_1d_t = Kokkos::View<double*, Layout, MemSpace>;

using vector_2d_t = std::vector<std::vector<double> >;

using kokkos_2d_table_t = Kokkos::View<double[300][10], Layout, MemSpace>;

using kokkos_1d_table_t = Kokkos::View<double[150], Layout, MemSpace>;


class MicroSedFuncVanillaKokkos
{
 private:
  kokkos_1d_t V_qr;//, V_nr, flux_qx, flux_nx;
  // kokkos_2d_t mu_r, lamr, rhofacr, inv_dzq, rho, inv_rho, t, tmparr1;
  // kokkos_2d_table_t vn_table, vm_table;
  // kokkos_1d_table_t mu_r_table;
  int _num_horz, _num_vert;

public:

  MicroSedFuncVanillaKokkos(int num_horz, int num_vert) :
    V_qr("V_qr", num_vert),
    // V_nr("V_nr", num_vert),
    // flux_qx("flux_qx", num_vert),
    // flux_nx("flux_nx", num_vert),
    // mu_r("mu_r", num_horz, num_vert),
    // lamr("lamr", num_horz, num_vert),
    // rhofacr("rhofacr", num_horz, num_vert),
    // inv_dzq("inv_dzq", num_horz, num_vert),
    // rho("rho", num_horz, num_vert),
    // inv_rho("inv_rho", num_horz, num_vert),
    // t("t", num_horz, num_vert),
    // tmparr1("tmparr1", num_horz, num_vert),
    // vn_table("VN_TABLE"), vm_table("VM_TABLE"),
    // mu_r_table("MU_R_TABLE"),
    _num_horz(num_horz), _num_vert(num_vert)
  {

    // // initialize on host
    reset();

    // auto mirror_vn_table = Kokkos::create_mirror_view(vn_table);
    // auto mirror_vm_table = Kokkos::create_mirror_view(vm_table);
    // auto mirror_mu_table = Kokkos::create_mirror_view(mu_r_table);

    // for (int i = 0; i < 300; ++i) {
    //   for (int k = 0; k < 10; ++k) {
    //     mirror_vn_table(i, k) = Globals::VN_TABLE[i][k];
    //     mirror_vm_table(i, k) = Globals::VM_TABLE[i][k];
    //   }
    // }

    // for (int i = 0; i < 150; ++i) {
    //   mirror_mu_table(i) = Globals::MU_R_TABLE[i];
    // }

    // // deep copy to device
    // Kokkos::deep_copy(vn_table, mirror_vn_table);
    // Kokkos::deep_copy(vm_table, mirror_vm_table);
    // Kokkos::deep_copy(mu_r_table, mirror_mu_table);
  }

  void reset()
  {
    //kokkos_1d_t V_qr("V_qr", _num_vert);//, V_nr, flux_qx, flux_nx;

    Kokkos::parallel_for("1d reset", Kokkos::RangePolicy<ExecSpace>(0, _num_vert), KOKKOS_LAMBDA(const int k) {
        //Kokkos::parallel_for("1d reset", _num_vert, KOKKOS_LAMBDA(int k) {
      V_qr(k) = 0.0;
      // V_nr(k) = 0.0;
      // flux_qx(k) = 0.0;
      // flux_nx(k) = 0.0;
    });


    //Kokkos::parallel_for("1d reset", Kokkos::RangePolicy<ExecSpace>(0, _num_vert), KOKKOS_LAMBDA(const int k) {
        //Kokkos::parallel_for("1d reset", _num_vert, KOKKOS_LAMBDA(int k) {
        //V_qr(k) = 0.0;
      // V_nr(k) = 0.0;
      // flux_qx(k) = 0.0;
      // flux_nx(k) = 0.0;
    //});

    // Kokkos::parallel_for("2d reset", _num_horz, KOKKOS_LAMBDA(int i) {
    //   for (int k = 0; k < _num_vert; ++k) {
    //     mu_r(i, k)    = 0.0;
    //     lamr(i, k)    = 0.0;
    //     rhofacr(i, k) = 0.0;
    //     inv_dzq(i, k) = 0.0;
    //     rho(i, k)     = 0.0;
    //     inv_rho(i, k) = 0.0;
    //     t(i, k)       = 0.0;
    //     tmparr1(i, k) = 0.0;
    //   }
    // });
  }
};

// using MemSpace = Kokkos::CudaSpace;
// using Layout = Kokkos::LayoutRight;
// using ExecSpace = Kokkos::Cuda;

// using kokkos_1d_t1 = Kokkos::View<double*, Layout, MemSpace>;
// using kokkos_1d_t2 = Kokkos::View<double*, MemSpace>;
// using kokkos_1d_t3 = Kokkos::View<double*, Layout>;
// using kokkos_1d_t4 = Kokkos::View<double*>;

int main(int argc, char** argv)
{

  Kokkos::initialize(argc, argv); {

    MicroSedFuncVanillaKokkos msvk(1, 111);

    const int num_vert = 111;
    kokkos_1d_t V_qr("V_qr", num_vert);//, V_nr, flux_qx, flux_nx;

    Kokkos::parallel_for("1d reset", Kokkos::RangePolicy<ExecSpace>(0, num_vert), KOKKOS_LAMBDA(const int k) {
        //Kokkos::parallel_for("1d reset", _num_vert, KOKKOS_LAMBDA(int k) {
      V_qr(k) = 0.0;
      // V_nr(k) = 0.0;
      // flux_qx(k) = 0.0;
      // flux_nx(k) = 0.0;
    });


//     const int size = 111;
//     kokkos_1d_t4 test_view("test view", size);

// //    Kokkos::parallel_for("test ||for", Kokkos::RangePolicy<ExecSpace>(0, size), KOKKOS_LAMBDA(int i) {
//     Kokkos::parallel_for("test ||for", size, KOKKOS_LAMBDA(int i) {
//       test_view(i) = 0.0;
    //});

  } Kokkos::finalize();

  return 0;
}
