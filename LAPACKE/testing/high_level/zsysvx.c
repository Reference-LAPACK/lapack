#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zsysvx. */
#define LAPACKE_ZSYSVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_sym(layout, N, a, LD);                              \
        lapacke_test_zfill_sym(layout, N, af, LD);                             \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_zfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zsysvx)(                         \
                               layout, 'N', 'U', N, NRHS, a, LD, af, LD, ipiv, \
                               b, LD, x, LD, &rcond, ferr, berr),              \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zsysvx)
{
    lapack_complex_double a[LD * LD];
    lapack_complex_double af[LD * LD];
    lapack_int ipiv[LD * LD];
    lapack_complex_double b[LD * LD];
    lapack_complex_double x[LD * LD];
    double rcond;
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx a fact=N uplo=U", l, N, N, a, LD, lapacke_test_region_upper,
            -6,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'N', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx b fact=N uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'N', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx a fact=N uplo=L", l, N, N, a, LD, lapacke_test_region_lower,
            -6,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'N', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx b fact=N uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'N', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx a fact=F uplo=U", l, N, N, a, LD, lapacke_test_region_upper,
            -6,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zsytrf)(layout, 'U', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'F', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx af fact=F uplo=U", l, N, N, af, LD,
            lapacke_test_region_upper, -8,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zsytrf)(layout, 'U', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'F', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx b fact=F uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zsytrf)(layout, 'U', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'F', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx a fact=F uplo=L", l, N, N, a, LD, lapacke_test_region_lower,
            -6,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zsytrf)(layout, 'L', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'F', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx af fact=F uplo=L", l, N, N, af, LD,
            lapacke_test_region_lower, -8,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zsytrf)(layout, 'L', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'F', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zsysvx b fact=F uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zsytrf)(layout, 'L', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'F', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_zfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_zfill_nan(layout, N, N, a, LD);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_zfill_nan(layout, N, N, af, LD);
        lapacke_test_check(
            "zsysvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zsysvx)(layout, 'N', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr,
                                       berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZSYSVX_ALLOC_TEST(0, 0, "zsysvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZSYSVX_ALLOC_TEST(0, 1, "zsysvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZSYSVX_ALLOC_TEST(0, 2, "zsysvx allocation count", 0);
    lapacke_test_check_alloc_count("zsysvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZSYSVX_ALLOC_TEST(1, 0, "zsysvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZSYSVX_ALLOC_TEST(1, 1, "zsysvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZSYSVX_ALLOC_TEST(1, 2, "zsysvx transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZSYSVX_ALLOC_TEST(1, 3, "zsysvx transpose alloc failure (af_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZSYSVX_ALLOC_TEST(1, 4, "zsysvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZSYSVX_ALLOC_TEST(1, 5, "zsysvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZSYSVX_ALLOC_TEST(1, 6, "zsysvx allocation count", 0);
    lapacke_test_check_alloc_count("zsysvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZSYSVX_ALLOC_TEST(2, 0, "zsysvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zsysvx invalid layout allocation count");
}
