F90=gfortran
F90FLAGS=-fPIC -O3
LDFLAGS=-shared

.SUFFIXES: .cpp .F90 .o

F90_SRC = shoc.F90 shoc_c_intr.F90
OBJECTS = $(F90_SRC:.F90=.o)

%.o: %.mod

.F90.o:
	$(F90) $(F90FLAGS) -c $< -o $@

libshoc: $(OBJECTS)
	$(F90) $(OBJECTS) $(LDFLAGS) -o libshoc.so

clean:
	rm -f *.o *.mod libshoc.so

libshoc: $(OBJECTS)
shoc_c_intr.o: shoc_c_intr.F90
shoc.o: shoc.F90
