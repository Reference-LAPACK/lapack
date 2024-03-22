CXX=g++
# CXXFLAGS=-std=c++17 -Wall -Wextra -pedantic -g
CXXFLAGS=-std=c++17 -O3 -march=native -ffast-math -g -DNDEBUG

FC=gfortran
FFLAGS=-Wall -Wextra -pedantic -g

all: ./example/test_lasr3 ./example/profile_lasr3 ./example/profile_steqr3 $(OBJFILES)

HEADERS = $(wildcard include/lapack_c/*.h) \
		  $(wildcard include/lapack_c/*/*.h) \
		  $(wildcard include/lapack_c/*/*/*.h) \
		  $(wildcard include/lapack_cpp/*.hpp) \
		  $(wildcard include/lapack_cpp/*/*.hpp) \
		  $(wildcard include/lapack_cpp/*/*/*.hpp)

SRCFILES = $(wildcard src/lapack_c/*.cpp) \
		   $(wildcard src/lapack_c/*.f90) \
		   $(wildcard src/lapack_cpp/*.cpp) \
		   $(wildcard src/lapack_fortran/*.f90)

OBJFILES = $(patsubst %.f90, %_f90.o ,$(patsubst %.cpp,%_cpp.o,$(SRCFILES)))

# Fortran files
src/lapack_fortran/%_f90.o: src/lapack_fortran/%.f90
	$(FC) $(FFLAGS) -cpp -c -o $@ $<

# C files

src/lapack_c/%_cpp.o: src/lapack_c/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -I./include -c -o $@ $<

src/lapack_c/%_f90.o: src/lapack_c/%.f90
	$(FC) $(FFLAGS) -cpp -c -o $@ $<

# C++ files
src/lapack_cpp/%_cpp.o: src/lapack_cpp/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -I./include -c -o $@ $<

# Test files
example/%: example/%.cpp $(OBJFILES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJFILES) -I./include -lstdc++ -lgfortran -lblas -llapack

clean:
	find . -type f -name '*.o' -delete