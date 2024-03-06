#include "lapack_c/steqr3.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/lapack/steqr3.hpp"

using namespace lapack_cpp;

// Actual definitions
extern "C"
{

    void lapack_c_ssteqr3(char compz,
                          lapack_idx n,
                          float *d,
                          float *e,
                          float *Z,
                          lapack_idx ldz,
                          float *work,
                          lapack_idx lwork,
                          float *rwork,
                          lapack_idx lrwork,
                          lapack_idx *info)
    {
        CompQ compz_ = char2compq(compz);

        Vector<float, lapack_idx> d_(n, d);
        Vector<float, lapack_idx> e_(n-1, e);
        Matrix<float, Layout::ColMajor, lapack_idx> Z_(n, n, Z, ldz);

        if (lwork == -1)
        {
            lapack_idx lwork_ = steqr3_workquery(compz_, d_, e_, Z_);
            *work = lwork_;
        }

        if (lrwork == -1)
        {
            lapack_idx lrwork_ = steqr3_rworkquery(compz_, d_, e_, Z_);
            *rwork = lrwork_;
        }

        if( lwork == -1 || lrwork == -1)
            return;

        MemoryBlock<float, lapack_idx> work_(lwork, work);
        MemoryBlock<float, lapack_idx> rwork_(lrwork, rwork);

        *info = steqr3(compz_, d_, e_, Z_, work_, rwork_);
    }


    void lapack_c_dsteqr3(char compz,
                          lapack_idx n,
                          double *d,
                          double *e,
                          double *Z,
                          lapack_idx ldz,
                          double *work,
                          lapack_idx lwork,
                          double *rwork,
                          lapack_idx lrwork,
                          lapack_idx *info)
    {
        CompQ compz_ = char2compq(compz);

        Vector<double, lapack_idx> d_(n, d);
        Vector<double, lapack_idx> e_(n-1, e);
        Matrix<double, Layout::ColMajor, lapack_idx> Z_(n, n, Z, ldz);

        if (lwork == -1)
        {
            lapack_idx lwork_ = steqr3_workquery(compz_, d_, e_, Z_);
            *work = lwork_;
        }

        if (lrwork == -1)
        {
            lapack_idx lrwork_ = steqr3_rworkquery(compz_, d_, e_, Z_);
            *rwork = lrwork_;
        }

        if( lwork == -1 || lrwork == -1)
            return;

        MemoryBlock<double, lapack_idx> work_(lwork, work);
        MemoryBlock<double, lapack_idx> rwork_(lrwork, rwork);

        *info = steqr3(compz_, d_, e_, Z_, work_, rwork_);
    }

    void lapack_c_csteqr3(char compz,
                          lapack_idx n,
                          float *d,
                          float *e,
                          lapack_float_complex *Z,
                          lapack_idx ldz,
                          lapack_float_complex *work,
                          lapack_idx lwork,
                          float *rwork,
                          lapack_idx lrwork,
                          lapack_idx *info)
    {
        CompQ compz_ = char2compq(compz);

        Vector<float, lapack_idx> d_(n, d);
        Vector<float, lapack_idx> e_(n-1, e);
        Matrix<lapack_float_complex, Layout::ColMajor, lapack_idx> Z_(n, n, Z, ldz);

        if (lwork == -1)
        {
            lapack_idx lwork_ = steqr3_workquery(compz_, d_, e_, Z_);
            *work = lwork_;
        }

        if (lrwork == -1)
        {
            lapack_idx lrwork_ = steqr3_rworkquery(compz_, d_, e_, Z_);
            *rwork = lrwork_;
        }

        if( lwork == -1 || lrwork == -1)
            return;

        MemoryBlock<lapack_float_complex, lapack_idx> work_(lwork, work);
        MemoryBlock<float, lapack_idx> rwork_(lrwork, rwork);

        *info = steqr3(compz_, d_, e_, Z_, work_, rwork_);
    }


    void lapack_c_zsteqr3(char compz,
                          lapack_idx n,
                          double *d,
                          double *e,
                          lapack_double_complex *Z,
                          lapack_idx ldz,
                          lapack_double_complex *work,
                          lapack_idx lwork,
                          double *rwork,
                          lapack_idx lrwork,
                          lapack_idx *info)
    {
        CompQ compz_ = char2compq(compz);

        Vector<double, lapack_idx> d_(n, d);
        Vector<double, lapack_idx> e_(n-1, e);
        Matrix<lapack_double_complex, Layout::ColMajor, lapack_idx> Z_(n, n, Z, ldz);

        if (lwork == -1)
        {
            lapack_idx lwork_ = steqr3_workquery(compz_, d_, e_, Z_);
            *work = lwork_;
        }

        if (lrwork == -1)
        {
            lapack_idx lrwork_ = steqr3_rworkquery(compz_, d_, e_, Z_);
            *rwork = lrwork_;
        }

        if( lwork == -1 || lrwork == -1)
            return;

        MemoryBlock<lapack_double_complex, lapack_idx> work_(lwork, work);
        MemoryBlock<double, lapack_idx> rwork_(lrwork, rwork);

        *info = steqr3(compz_, d_, e_, Z_, work_, rwork_);
    }

} // extern "C"