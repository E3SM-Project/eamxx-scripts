KOKKOS = $(HOME)/lib/kokkos/cpu
#CXXFLAGS = -g -DTRIDIAG_MIMIC_GPU -std=c++11 -fPIC -fopenmp -I$(KOKKOS)/include -lblas -DTRIDIAG_HAVE_LAPACK
CXXFLAGS = -O3 -std=c++11 -fPIC -fopenmp -I$(KOKKOS)/include -fopenmp -DTRIDIAG_HAVE_LAPACK
LDFLAGS = -L$(KOKKOS)/lib -lkokkos -ldl -fopenmp -llapack -lblas 
