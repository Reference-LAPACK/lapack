using BinaryBuilder, Pkg

haskey(ENV, "BLAS_LAPACK_RELEASE") || error("The environment variable BLAS_LAPACK_RELEASE is not defined.")
haskey(ENV, "BLAS_LAPACK_COMMIT") || error("The environment variable BLAS_LAPACK_COMMIT is not defined.")
haskey(ENV, "BLAS_LAPACK_URL") || error("The environment variable BLAS_LAPACK_URL is not defined.")

name = "blas_lapack"
version = VersionNumber(ENV["BLAS_LAPACK_RELEASE"])

# Collection of sources required to complete build
sources = [
    GitSource(ENV["BLAS_LAPACK_URL"], ENV["BLAS_LAPACK_COMMIT"])
]

# Bash recipe for building across all platforms
script = raw"""
cd ${WORKSPACE}/srcdir/lapack

# FortranCInterface_VERIFY fails on macOS, but it's not actually needed for the current build
sed -i 's/FortranCInterface_VERIFY/# FortranCInterface_VERIFY/g' ./CBLAS/CMakeLists.txt
sed -i 's/FortranCInterface_VERIFY/# FortranCInterface_VERIFY/g' ./LAPACKE/include/CMakeLists.txt

mkdir build && cd build
cmake .. \
   -DCBLAS=ON \
   -DLAPACKE=ON \
   -DCMAKE_INSTALL_PREFIX="$prefix" \
   -DCMAKE_FIND_ROOT_PATH="$prefix" \
   -DCMAKE_TOOLCHAIN_FILE="${CMAKE_TARGET_TOOLCHAIN}" \
   -DCMAKE_BUILD_TYPE=Release \
   -DBUILD_SHARED_LIBS=OFF \
   -DBUILD_INDEX64_EXT_API=OFF \
   -DTEST_FORTRAN_COMPILER=OFF \
   -DLAPACKE_WITH_TMG=OFF

make -j${nproc}
make install

install_license $WORKSPACE/srcdir/lapack/LICENSE
"""

# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line
platforms = supported_platforms()
platforms = expand_gfortran_versions(platforms)

# The products that we will ensure are always built
products = [
    FileProduct("lib/libblas.a", :libblas_a),
    FileProduct("lib/libcblas.a", :libcblas_a),
    FileProduct("lib/liblapack.a", :liblapack_a),
    FileProduct("lib/liblapacke.a", :liblapacke_a),
    # LibraryProduct("libblas", :libblas),
    # LibraryProduct("libcblas", :libcblas),
    # LibraryProduct("liblapack", :liblapack),
    # LibraryProduct("liblapacke", :liblapacke),
]

# Dependencies that must be installed before this package can be built
dependencies = [
    Dependency(PackageSpec(name="CompilerSupportLibraries_jll", uuid="e66e0078-7015-5450-92f7-15fbd957f2ae")),
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies; julia_compat="1.6")
