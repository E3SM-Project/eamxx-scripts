KOKKOS = $(HOME)/lib/kokkos/cpu
CXXFLAGS = -g -DMIMIC_GPU -std=c++11 -fPIC -fopenmp -I$(KOKKOS)/include -fopenmp
#CXXFLAGS = -O3 -std=c++11 -fPIC -fopenmp -I$(KOKKOS)/include -fopenmp
LDFLAGS = -L$(KOKKOS)/lib -lkokkos -ldl -fopenmp
