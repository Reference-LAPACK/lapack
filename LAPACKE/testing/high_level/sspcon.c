#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sspcon. */
#define LAPACKE_SSPCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_sp(layout, 'U', N, ap);                             \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        anorm[0] = 1.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sspcon)(layout, 'U', N, ap,      \
                                                      ipiv, anorm[0], &rcond), \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sspcon)
{
    float ap[LD * LD];
    lapack_int ipiv[LD * LD];
    float anorm[1];
    float rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sspcon anorm uplo=U", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_sfill_sp(layout, 'U', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_sspcon)(layout, 'U', N, ap, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "sspcon ap uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_sfill_sp(layout, 'U', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_sspcon)(layout, 'U', N, ap, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "sspcon anorm uplo=L", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_sfill_sp(layout, 'L', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_sspcon)(layout, 'L', N, ap, ipiv, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "sspcon ap uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_sfill_sp(layout, 'L', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_sspcon)(layout, 'L', N, ap, ipiv, anorm[0],
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        anorm[0] = 1.0f;
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_sfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_check("sspcon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_sspcon)(layout, 'U', N, ap, ipiv,
                                                      anorm[0], &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SSPCON_ALLOC_TEST(0, 0, "sspcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSPCON_ALLOC_TEST(0, 1, "sspcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSPCON_ALLOC_TEST(0, 2, "sspcon allocation count", 0);
    lapacke_test_check_alloc_count("sspcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SSPCON_ALLOC_TEST(1, 0, "sspcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSPCON_ALLOC_TEST(1, 1, "sspcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSPCON_ALLOC_TEST(1, 2, "sspcon transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SSPCON_ALLOC_TEST(1, 3, "sspcon allocation count", 0);
    lapacke_test_check_alloc_count("sspcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SSPCON_ALLOC_TEST(2, 0, "sspcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sspcon invalid layout allocation count");
}
