#include <Kokkos_Core.hpp>

#include <mpi.h>

#include <iostream>
#include <string>

template<typename VT>
std::string v2str (const VT& v) {
  if (v.size()==0) {
    return "";
  }

  auto vh = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(vh,v);
  std::string s = std::to_string(vh.data()[0]);

  for (size_t i=1; i<vh.size(); ++i) {
    s += " " + std::to_string(vh.data()[i]);
  }
  return s;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  Kokkos::initialize(argc,argv);

  MPI_Comm comm(MPI_COMM_WORLD);
  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  // Wrap in scope, so view is destroyed before Kokkos::finalize
  {

    Kokkos::View<int*> v("",10);
    Kokkos::deep_copy(v,rank);
    std::cout << "rank=" << rank << ", pre-bcast: v: " << v2str(v) << "\n";

    MPI_Bcast (v.data(),v.size(),MPI_INT,0,comm);

    std::cout << "rank=" << rank << ", post-bcast: v: " << v2str(v) << "\n";
  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
