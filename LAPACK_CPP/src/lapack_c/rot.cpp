#ifndef USE_FORTRAN_BLAS
#include "lapack_c/rot.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/blas/rot.hpp"

using namespace lapack_cpp;

// Actual definitions
extern "C"
{

    void lapack_c_drot(lapack_idx n,
                       double *x,
                       lapack_idx incx,
                       double *y,
                       lapack_idx incy,
                       double c,
                       double s)
    {
        Vector<double, lapack_idx> x_(n, x, incx);
        Vector<double, lapack_idx> y_(n, y, incy);
        rot(x_, y_, c, s);
    }

    void lapack_c_srot(lapack_idx n,
                       float *x,
                       lapack_idx incx,
                       float *y,
                       lapack_idx incy,
                       float c,
                       float s)
    {
        Vector<float, lapack_idx> x_(n, x, incx);
        Vector<float, lapack_idx> y_(n, y, incy);
        rot(x_, y_, c, s);
    }

    void lapack_c_crot(lapack_idx n,
                       lapack_float_complex *x,
                       lapack_idx incx,
                       lapack_float_complex *y,
                       lapack_idx incy,
                       float c,
                       lapack_float_complex s)
    {
        Vector<lapack_float_complex, lapack_idx> x_(n, x, incx);
        Vector<lapack_float_complex, lapack_idx> y_(n, y, incy);
        rot(x_, y_, c, s);
    }

    void lapack_c_zrot(lapack_idx n,
                       lapack_double_complex *x,
                       lapack_idx incx,
                       lapack_double_complex *y,
                       lapack_idx incy,
                       double c,
                       lapack_double_complex s)
    {
        Vector<lapack_double_complex, lapack_idx> x_(n, x, incx);
        Vector<lapack_double_complex, lapack_idx> y_(n, y, incy);
        rot(x_, y_, c, s);
    }

} // extern "C"
#endif