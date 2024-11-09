#ifdef USE_FORTRAN_BLAS
    #include <cassert>
    #include <complex>
    #include <type_traits>

    #include "lapack_c.h"
    #include "lapack_cpp/base.hpp"
    #include "lapack_cpp/blas/rot.hpp"

namespace lapack_cpp {

template<>
void rot(const Vector<float, lapack_idx_t>& x,
         const Vector<float, lapack_idx_t>& y,
         float c,
         float s)
{
    assert(x.size() == y.size());
    lapack_c_srot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
}

template<>
void rot(const Vector<double, lapack_idx_t>& x,
         const Vector<double, lapack_idx_t>& y,
         double c,
         double s)
{
    assert(x.size() == y.size());
    lapack_c_drot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
}

template<>
void rot(const Vector<std::complex<float>, lapack_idx_t>& x,
         const Vector<std::complex<float>, lapack_idx_t>& y,
         float c,
         std::complex<float> s)
{
    assert(x.size() == y.size());
    lapack_c_crot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
}

template<>
void rot(const Vector<std::complex<double>, lapack_idx_t>& x,
         const Vector<std::complex<double>, lapack_idx_t>& y,
         double c,
         std::complex<double> s)
{
    assert(x.size() == y.size());
    lapack_c_zrot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
}

template<>
void rot(const Vector<float, int>& x,
         const Vector<float, int>& y,
         float c,
         float s)
{
    assert(x.size() == y.size());
    lapack_c_srot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
}

template<>
void rot(const Vector<double, int>& x,
         const Vector<double, int>& y,
         double c,
         double s)
{
    assert(x.size() == y.size());
    lapack_c_drot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
}

template<>
void rot(const Vector<std::complex<float>, int>& x,
         const Vector<std::complex<float>, int>& y,
         float c,
         std::complex<float> s)
{
    assert(x.size() == y.size());
    lapack_c_crot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
}

template<>
void rot(const Vector<std::complex<double>, int>& x,
         const Vector<std::complex<double>, int>& y,
         double c,
         std::complex<double> s)
{
    assert(x.size() == y.size());
    lapack_c_zrot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
}

}  // namespace lapack_cpp
#endif  // USE_FORTRAN_BLAS