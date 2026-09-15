#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dspcon. */
#define LAPACKE_DSPCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_sp(layout, 'U', N, ap);                             \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        anorm[0] = 1.0;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dspcon)(layout, 'U', N, ap,      \
                                                      ipiv, anorm[0], &rcond), \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dspcon)
{
    double ap[LD * LD];
    lapack_int ipiv[LD * LD];
    double anorm[1];
    double rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dspcon anorm uplo=U", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_dfill_sp(layout, 'U', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_dspcon)(layout, 'U', N, ap, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "dspcon ap uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_dfill_sp(layout, 'U', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_dspcon)(layout, 'U', N, ap, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "dspcon anorm uplo=L", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_dfill_sp(layout, 'L', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_dspcon)(layout, 'L', N, ap, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "dspcon ap uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_dfill_sp(layout, 'L', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_dspcon)(layout, 'L', N, ap, ipiv, anorm[0],
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        anorm[0] = 1.0;
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_dfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_check("dspcon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_dspcon)(layout, 'U', N, ap, ipiv,
                                                      anorm[0], &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DSPCON_ALLOC_TEST(0, 0, "dspcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DSPCON_ALLOC_TEST(0, 1, "dspcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DSPCON_ALLOC_TEST(0, 2, "dspcon allocation count", 0);
    lapacke_test_check_alloc_count("dspcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DSPCON_ALLOC_TEST(1, 0, "dspcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DSPCON_ALLOC_TEST(1, 1, "dspcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DSPCON_ALLOC_TEST(1, 2, "dspcon transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DSPCON_ALLOC_TEST(1, 3, "dspcon allocation count", 0);
    lapacke_test_check_alloc_count("dspcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DSPCON_ALLOC_TEST(2, 0, "dspcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dspcon invalid layout allocation count");
}
