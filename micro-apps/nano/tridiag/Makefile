include make.inc

SOURCES = tridiag.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

all: exe

libtridiag: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -shared -o libtridiag.so

exe: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o tridiag

clean:
	rm -f *.o *.so
