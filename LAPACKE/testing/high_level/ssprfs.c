#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call ssprfs. */
#define LAPACKE_SSPRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_sp(layout, 'U', N, ap);                             \
        lapacke_test_sfill_sp(layout, 'U', N, afp);                            \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_sfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_ssprfs)(layout, 'U', N, NRHS,    \
                                                      ap, afp, ipiv, b, LD, x, \
                                                      LD, ferr, berr),         \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(ssprfs)
{
    float ap[LD * LD];
    float afp[LD * LD];
    lapack_int ipiv[LD * LD];
    float b[LD * LD];
    float x[LD * LD];
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP("ssprfs afp uplo=U", l, N * (N + 1) / 2, 1, afp,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -6,
                                (lapacke_test_sfill_sp(layout, 'U', N, ap),
                                 lapacke_test_sfill_sp(layout, 'U', N, afp),
                                 lapacke_test_fill_ipiv(LD * LD, ipiv),
                                 lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
                                 lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
                                 lapacke_test_sfill_vec(LD * LD, ferr),
                                 lapacke_test_sfill_vec(LD * LD, berr)),
                                API_SUFFIX(LAPACKE_ssprfs)(layout, 'U', N, NRHS,
                                                           ap, afp, ipiv, b, LD,
                                                           x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP("ssprfs ap uplo=U", l, N * (N + 1) / 2, 1, ap,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -5,
                                (lapacke_test_sfill_sp(layout, 'U', N, ap),
                                 lapacke_test_sfill_sp(layout, 'U', N, afp),
                                 lapacke_test_fill_ipiv(LD * LD, ipiv),
                                 lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
                                 lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
                                 lapacke_test_sfill_vec(LD * LD, ferr),
                                 lapacke_test_sfill_vec(LD * LD, berr)),
                                API_SUFFIX(LAPACKE_ssprfs)(layout, 'U', N, NRHS,
                                                           ap, afp, ipiv, b, LD,
                                                           x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "ssprfs b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_sfill_sp(layout, 'U', N, ap),
             lapacke_test_sfill_sp(layout, 'U', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ssprfs)(layout, 'U', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "ssprfs x uplo=U", l, N, NRHS, x, LD, lapacke_test_region_full, -10,
            (lapacke_test_sfill_sp(layout, 'U', N, ap),
             lapacke_test_sfill_sp(layout, 'U', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ssprfs)(layout, 'U', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP("ssprfs afp uplo=L", l, N * (N + 1) / 2, 1, afp,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -6,
                                (lapacke_test_sfill_sp(layout, 'L', N, ap),
                                 lapacke_test_sfill_sp(layout, 'L', N, afp),
                                 lapacke_test_fill_ipiv(LD * LD, ipiv),
                                 lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
                                 lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
                                 lapacke_test_sfill_vec(LD * LD, ferr),
                                 lapacke_test_sfill_vec(LD * LD, berr)),
                                API_SUFFIX(LAPACKE_ssprfs)(layout, 'L', N, NRHS,
                                                           ap, afp, ipiv, b, LD,
                                                           x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP("ssprfs ap uplo=L", l, N * (N + 1) / 2, 1, ap,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -5,
                                (lapacke_test_sfill_sp(layout, 'L', N, ap),
                                 lapacke_test_sfill_sp(layout, 'L', N, afp),
                                 lapacke_test_fill_ipiv(LD * LD, ipiv),
                                 lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
                                 lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
                                 lapacke_test_sfill_vec(LD * LD, ferr),
                                 lapacke_test_sfill_vec(LD * LD, berr)),
                                API_SUFFIX(LAPACKE_ssprfs)(layout, 'L', N, NRHS,
                                                           ap, afp, ipiv, b, LD,
                                                           x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "ssprfs b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_sfill_sp(layout, 'L', N, ap),
             lapacke_test_sfill_sp(layout, 'L', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ssprfs)(layout, 'L', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "ssprfs x uplo=L", l, N, NRHS, x, LD, lapacke_test_region_full, -10,
            (lapacke_test_sfill_sp(layout, 'L', N, ap),
             lapacke_test_sfill_sp(layout, 'L', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ssprfs)(layout, 'L', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_sfill_nan(layout, N * (N + 1) / 2, 1, afp,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_sfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_sfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_sfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "ssprfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_ssprfs)(layout, 'U', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SSPRFS_ALLOC_TEST(0, 0, "ssprfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSPRFS_ALLOC_TEST(0, 1, "ssprfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSPRFS_ALLOC_TEST(0, 2, "ssprfs allocation count", 0);
    lapacke_test_check_alloc_count("ssprfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SSPRFS_ALLOC_TEST(1, 0, "ssprfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSPRFS_ALLOC_TEST(1, 1, "ssprfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSPRFS_ALLOC_TEST(1, 2, "ssprfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SSPRFS_ALLOC_TEST(1, 3, "ssprfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SSPRFS_ALLOC_TEST(1, 4, "ssprfs transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SSPRFS_ALLOC_TEST(1, 5, "ssprfs transpose alloc failure (afp_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SSPRFS_ALLOC_TEST(1, 6, "ssprfs allocation count", 0);
    lapacke_test_check_alloc_count("ssprfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SSPRFS_ALLOC_TEST(2, 0, "ssprfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("ssprfs invalid layout allocation count");
}
