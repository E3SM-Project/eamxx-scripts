F90=gfortran
CPP=g++

CFLAGS=-fPIC -O3
#CFLAGS=-fPIC -g -O0 -DWITH_FPE

LDFLAGS=-shared

.SUFFIXES: .cpp .F90 .o

CPP_SRC = fpe.cpp
F90_SRC = shoc.F90 shoc_c_intr.F90
OBJECTS = $(F90_SRC:.F90=.o) $(CPP_SRC:.cpp=.o)

%.o: %.mod

.cpp.o:
	$(CXX) $(CFLAGS) -c $< -o $@

.F90.o:
	$(F90) $(CFLAGS) -c $< -o $@

libshoc: $(OBJECTS)
	$(F90) $(OBJECTS) $(LDFLAGS) -o libshoc.so

clean:
	rm -f *.o *.mod libshoc.so

libshoc: $(OBJECTS)
shoc_c_intr.o: shoc_c_intr.F90
shoc.o: shoc.F90
