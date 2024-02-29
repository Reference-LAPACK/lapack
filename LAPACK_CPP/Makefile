CXX=g++
# CXXFLAGS=-std=c++17 -Wall -Wextra -pedantic -g -Og
CXXFLAGS=-std=c++17 -O3 -march=native -ffast-math -g -DNDEBUG

FC=gfortran
FFLAGS=-Wall -Wextra -pedantic -g

all: ./example/test_lasr3 ./example/profile_lasr3

HEADERS = $(wildcard lapack_c/include/*.h) \
		  $(wildcard lapack_c/include/*/*.h) \
		  $(wildcard lapack_c/include/*/*/*.h) \
		  $(wildcard lapack_cpp/include/*.hpp) \
		  $(wildcard lapack_cpp/include/*/*.hpp) \
		  $(wildcard lapack_cpp/include/*/*/*.hpp)

OBJFILES = src/lapack_cpp/rot_cpp.o \
		   src/lapack_cpp/lartg_cpp.o \
		   src/lapack_c/rot_c.o \
		   src/lapack_c/lartg_c.o \
		   src/lapack_c/lasr3_c.o \
		   src/lapack_fortran/slasr3.o \
		   src/lapack_fortran/dlasr3.o \
		   src/lapack_fortran/clasr3.o \
		   src/lapack_fortran/zlasr3.o \

# Fortran files
src/lapack_fortran/%.o: src/lapack_fortran/%.f90
	$(FC) $(FFLAGS) -cpp -c -o $@ $<

# C files

src/lapack_c/%.o: src/lapack_c/%.cpp
	$(CXX) $(CXXFLAGS) -I./include -c -o $@ $<

src/lapack_c/%.o: src/lapack_c/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

# C++ files
src/lapack_cpp/%.o: src/lapack_cpp/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -I./include -c -o $@ $<

# Test files
example/test_lasr3: example/test_lasr3.cpp $(OBJFILES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJFILES) -I./include -lstdc++ -lgfortran -lblas -llapack

example/profile_lasr3: example/profile_lasr3.cpp $(OBJFILES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJFILES) -I./include -lstdc++ -lgfortran -lblas -llapack

clean:
	find . -type f -name '*.o' -delete