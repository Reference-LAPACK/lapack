#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cspsvx. */
#define LAPACKE_CSPSVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_sp(layout, 'U', N, ap);                             \
        lapacke_test_cfill_sp(layout, 'U', N, afp);                            \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_cfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cspsvx)(                         \
                               layout, 'N', 'U', N, NRHS, ap, afp, ipiv, b,    \
                               LD, x, LD, &rcond, ferr, berr),                 \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cspsvx)
{
    lapack_complex_float ap[LD * LD];
    lapack_complex_float afp[LD * LD];
    lapack_int ipiv[LD * LD];
    lapack_complex_float b[LD * LD];
    lapack_complex_float x[LD * LD];
    float rcond;
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx ap fact=N uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_cfill_sp(layout, 'U', N, ap),
             lapacke_test_cfill_sp(layout, 'U', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'N', 'U', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx b fact=N uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_sp(layout, 'U', N, ap),
             lapacke_test_cfill_sp(layout, 'U', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'N', 'U', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx ap fact=N uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_cfill_sp(layout, 'L', N, ap),
             lapacke_test_cfill_sp(layout, 'L', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'N', 'L', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx b fact=N uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_sp(layout, 'L', N, ap),
             lapacke_test_cfill_sp(layout, 'L', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'N', 'L', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx afp fact=F uplo=U", l, N * (N + 1) / 2, 1, afp,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -7,
            (lapacke_test_cfill_sp(layout, 'U', N, ap),
             lapacke_test_cfill_sp(layout, 'U', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_csptrf)(layout, 'U', N, afp, ipiv)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'F', 'U', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx ap fact=F uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_cfill_sp(layout, 'U', N, ap),
             lapacke_test_cfill_sp(layout, 'U', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_csptrf)(layout, 'U', N, afp, ipiv)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'F', 'U', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx b fact=F uplo=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_sp(layout, 'U', N, ap),
             lapacke_test_cfill_sp(layout, 'U', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_csptrf)(layout, 'U', N, afp, ipiv)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'F', 'U', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx afp fact=F uplo=L", l, N * (N + 1) / 2, 1, afp,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -7,
            (lapacke_test_cfill_sp(layout, 'L', N, ap),
             lapacke_test_cfill_sp(layout, 'L', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_csptrf)(layout, 'L', N, afp, ipiv)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'F', 'L', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx ap fact=F uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -6,
            (lapacke_test_cfill_sp(layout, 'L', N, ap),
             lapacke_test_cfill_sp(layout, 'L', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_csptrf)(layout, 'L', N, afp, ipiv)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'F', 'L', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cspsvx b fact=F uplo=L", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_sp(layout, 'L', N, ap),
             lapacke_test_cfill_sp(layout, 'L', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_csptrf)(layout, 'L', N, afp, ipiv)),
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'F', 'L', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_cfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_cfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_cfill_nan(layout, N * (N + 1) / 2, 1, afp,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_check(
            "cspsvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cspsvx)(layout, 'N', 'U', N, NRHS, ap, afp, ipiv,
                                       b, LD, x, LD, &rcond, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CSPSVX_ALLOC_TEST(0, 0, "cspsvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CSPSVX_ALLOC_TEST(0, 1, "cspsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CSPSVX_ALLOC_TEST(0, 2, "cspsvx allocation count", 0);
    lapacke_test_check_alloc_count("cspsvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CSPSVX_ALLOC_TEST(1, 0, "cspsvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CSPSVX_ALLOC_TEST(1, 1, "cspsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CSPSVX_ALLOC_TEST(1, 2, "cspsvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CSPSVX_ALLOC_TEST(1, 3, "cspsvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CSPSVX_ALLOC_TEST(1, 4, "cspsvx transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CSPSVX_ALLOC_TEST(1, 5, "cspsvx transpose alloc failure (afp_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CSPSVX_ALLOC_TEST(1, 6, "cspsvx allocation count", 0);
    lapacke_test_check_alloc_count("cspsvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CSPSVX_ALLOC_TEST(2, 0, "cspsvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cspsvx invalid layout allocation count");
}
