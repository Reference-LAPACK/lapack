#include "lapack_c/lasr3.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/lapack/lasr3.hpp"

using namespace lapack_cpp;

// Actual definitions
extern "C"
{

    void lapack_c_slasr3(char side,
                         char direct,
                         lapack_idx m,
                         lapack_idx n,
                         lapack_idx k,
                         const float *C,
                         lapack_idx ldc,
                         const float *S,
                         lapack_idx lds,
                         float *A,
                         lapack_idx lda,
                         float *work,
                         lapack_idx lwork)
    {
        Side side_ = char2side(side);
        Direction direct_ = char2direct(direct);
        const lapack_idx n_rot = side_ == Side::Left ? m - 1 : n - 1;

        ConstMatrix<float, Layout::ColMajor, lapack_idx> C_(n_rot, k, C, ldc);
        ConstMatrix<float, Layout::ColMajor, lapack_idx> S_(n_rot, k, S, lds);
        Matrix<float, Layout::ColMajor, lapack_idx> A_(m, n, A, lda);

        if (lwork == -1)
        {
            lapack_idx lwork_ = lasr3_workquery(side_, direct_, C_, S_, A_);
            *work = lwork_;
            return;
        }

        MemoryBlock<float, lapack_idx> work_(lwork, work);
        lasr3(side_, direct_, C_, S_, A_, work_);
    }

    void lapack_c_dlasr3(char side,
                         char direct,
                         lapack_idx m,
                         lapack_idx n,
                         lapack_idx k,
                         const double *C,
                         lapack_idx ldc,
                         const double *S,
                         lapack_idx lds,
                         double *A,
                         lapack_idx lda,
                         double *work,
                         lapack_idx lwork)
    {
        Side side_ = char2side(side);
        Direction direct_ = char2direct(direct);
        const lapack_idx n_rot = side_ == Side::Left ? m - 1 : n - 1;

        ConstMatrix<double, Layout::ColMajor, lapack_idx> C_(n_rot, k, C, ldc);
        ConstMatrix<double, Layout::ColMajor, lapack_idx> S_(n_rot, k, S, lds);
        Matrix<double, Layout::ColMajor, lapack_idx> A_(m, n, A, lda);

        if (lwork == -1)
        {
            lapack_idx lwork_ = lasr3_workquery(side_, direct_, C_, S_, A_);
            *work = lwork_;
            return;
        }

        MemoryBlock<double, lapack_idx> work_(lwork, work);
        lasr3(side_, direct_, C_, S_, A_, work_);
    }

    void lapack_c_clasr3(char side,
                         char direct,
                         lapack_idx m,
                         lapack_idx n,
                         lapack_idx k,
                         const float *C,
                         lapack_idx ldc,
                         const lapack_float_complex *S,
                         lapack_idx lds,
                         lapack_float_complex *A,
                         lapack_idx lda,
                         lapack_float_complex *work,
                         lapack_idx lwork)
    {
        Side side_ = char2side(side);
        Direction direct_ = char2direct(direct);
        const lapack_idx n_rot = side_ == Side::Left ? m - 1 : n - 1;

        ConstMatrix<float, Layout::ColMajor, lapack_idx> C_(n_rot, k, C, ldc);
        ConstMatrix<lapack_float_complex, Layout::ColMajor, lapack_idx> S_(n_rot, k,
                                                                           S, lds);
        Matrix<lapack_float_complex, Layout::ColMajor, lapack_idx> A_(m, n, A, lda);

        if (lwork == -1)
        {
            lapack_idx lwork_ = lasr3_workquery(side_, direct_, C_, S_, A_);
            *work = lwork_;
            return;
        }

        MemoryBlock<lapack_float_complex, lapack_idx> work_(lwork, work);
        lasr3(side_, direct_, C_, S_, A_, work_);
    }

    void lapack_c_zlasr3(char side,
                         char direct,
                         lapack_idx m,
                         lapack_idx n,
                         lapack_idx k,
                         const double *C,
                         lapack_idx ldc,
                         const lapack_double_complex *S,
                         lapack_idx lds,
                         lapack_double_complex *A,
                         lapack_idx lda,
                         lapack_double_complex *work,
                         lapack_idx lwork)
    {
        Side side_ = char2side(side);
        Direction direct_ = char2direct(direct);
        const lapack_idx n_rot = side_ == Side::Left ? m - 1 : n - 1;

        ConstMatrix<double, Layout::ColMajor, lapack_idx> C_(n_rot, k, C, ldc);
        ConstMatrix<lapack_double_complex, Layout::ColMajor, lapack_idx> S_(
            n_rot, k, S, lds);
        Matrix<lapack_double_complex, Layout::ColMajor, lapack_idx> A_(m, n, A,
                                                                       lda);

        if (lwork == -1)
        {
            lapack_idx lwork_ = lasr3_workquery(side_, direct_, C_, S_, A_);
            *work = lwork_;
            return;
        }

        MemoryBlock<lapack_double_complex, lapack_idx> work_(lwork, work);
        lasr3(side_, direct_, C_, S_, A_, work_);
    }

    void lapack_c_sclasr3(char side,
                          char direct,
                          lapack_idx m,
                          lapack_idx n,
                          lapack_idx k,
                          const float *C,
                          lapack_idx ldc,
                          const float *S,
                          lapack_idx lds,
                          lapack_float_complex *A,
                          lapack_idx lda,
                          lapack_float_complex *work,
                          lapack_idx lwork)
    {
        Side side_ = char2side(side);
        Direction direct_ = char2direct(direct);
        const lapack_idx n_rot = side_ == Side::Left ? m - 1 : n - 1;

        ConstMatrix<float, Layout::ColMajor, lapack_idx> C_(n_rot, k, C, ldc);
        ConstMatrix<float, Layout::ColMajor, lapack_idx> S_(n_rot, k, S, lds);
        Matrix<lapack_float_complex, Layout::ColMajor, lapack_idx> A_(m, n, A, lda);

        if (lwork == -1)
        {
            lapack_idx lwork_ = lasr3_workquery(side_, direct_, C_, S_, A_);
            *work = lwork_;
            return;
        }

        MemoryBlock<lapack_float_complex, lapack_idx> work_(lwork, work);
        lasr3(side_, direct_, C_, S_, A_, work_);
    }

    void lapack_c_dzlasr3(char side,
                          char direct,
                          lapack_idx m,
                          lapack_idx n,
                          lapack_idx k,
                          const double *C,
                          lapack_idx ldc,
                          const double *S,
                          lapack_idx lds,
                          lapack_double_complex *A,
                          lapack_idx lda,
                          lapack_double_complex *work,
                          lapack_idx lwork)
    {
        Side side_ = char2side(side);
        Direction direct_ = char2direct(direct);
        const lapack_idx n_rot = side_ == Side::Left ? m - 1 : n - 1;

        ConstMatrix<double, Layout::ColMajor, lapack_idx> C_(n_rot, k, C, ldc);
        ConstMatrix<double, Layout::ColMajor, lapack_idx> S_(n_rot, k, S, lds);
        Matrix<lapack_double_complex, Layout::ColMajor, lapack_idx> A_(m, n, A,
                                                                       lda);

        if (lwork == -1)
        {
            lapack_idx lwork_ = lasr3_workquery(side_, direct_, C_, S_, A_);
            *work = lwork_;
            return;
        }

        MemoryBlock<lapack_double_complex, lapack_idx> work_(lwork, work);
        lasr3(side_, direct_, C_, S_, A_, work_);
    }

} // extern "C"