#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call chesvx. */
#define LAPACKE_CHESVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_sym(layout, N, a, LD);                              \
        lapacke_test_cfill_sym(layout, N, af, LD);                             \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_cfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_chesvx)(                         \
                               layout, 'N', 'U', N, NRHS, a, LD, af, LD, ipiv, \
                               b, LD, x, LD, &rcond, ferr, berr),              \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(chesvx)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float af[LD * LD];
    lapack_int ipiv[LD * LD];
    lapack_complex_float b[LD * LD];
    lapack_complex_float x[LD * LD];
    float rcond;
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx a fact=N uplo=U", l, N, N, a, LD, lapacke_test_region_upper,
            -6,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'N', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx b fact=N uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'N', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx a fact=N uplo=L", l, N, N, a, LD, lapacke_test_region_lower,
            -6,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'N', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx b fact=N uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'N', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx a fact=F uplo=U", l, N, N, a, LD, lapacke_test_region_upper,
            -6,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_chetrf)(layout, 'U', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'F', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx af fact=F uplo=U", l, N, N, af, LD,
            lapacke_test_region_upper, -8,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_chetrf)(layout, 'U', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'F', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx b fact=F uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_chetrf)(layout, 'U', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'F', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx a fact=F uplo=L", l, N, N, a, LD, lapacke_test_region_lower,
            -6,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_chetrf)(layout, 'L', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'F', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx af fact=F uplo=L", l, N, N, af, LD,
            lapacke_test_region_lower, -8,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_chetrf)(layout, 'L', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'F', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "chesvx b fact=F uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_cfill_sym(layout, N, a, LD),
             lapacke_test_cfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_chetrf)(layout, 'L', N, af, LD, ipiv)),
            API_SUFFIX(LAPACKE_chesvx)(layout, 'F', 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_cfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_cfill_nan(layout, N, N, af, LD);
        lapacke_test_check(
            "chesvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_chesvx)(layout, 'N', 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, &rcond, ferr,
                                       berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CHESVX_ALLOC_TEST(0, 0, "chesvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHESVX_ALLOC_TEST(0, 1, "chesvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHESVX_ALLOC_TEST(0, 2, "chesvx allocation count", 0);
    lapacke_test_check_alloc_count("chesvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CHESVX_ALLOC_TEST(1, 0, "chesvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHESVX_ALLOC_TEST(1, 1, "chesvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHESVX_ALLOC_TEST(1, 2, "chesvx transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CHESVX_ALLOC_TEST(1, 3, "chesvx transpose alloc failure (af_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CHESVX_ALLOC_TEST(1, 4, "chesvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CHESVX_ALLOC_TEST(1, 5, "chesvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CHESVX_ALLOC_TEST(1, 6, "chesvx allocation count", 0);
    lapacke_test_check_alloc_count("chesvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CHESVX_ALLOC_TEST(2, 0, "chesvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("chesvx invalid layout allocation count");
}
