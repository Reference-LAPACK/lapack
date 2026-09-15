#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dsprfs. */
#define LAPACKE_DSPRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_sp(layout, 'U', N, ap);                             \
        lapacke_test_dfill_sp(layout, 'U', N, afp);                            \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_dfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dsprfs)(layout, 'U', N, NRHS,    \
                                                      ap, afp, ipiv, b, LD, x, \
                                                      LD, ferr, berr),         \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dsprfs)
{
    double ap[LD * LD];
    double afp[LD * LD];
    lapack_int ipiv[LD * LD];
    double b[LD * LD];
    double x[LD * LD];
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP("dsprfs afp uplo=U", l, N * (N + 1) / 2, 1, afp,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -6,
                                (lapacke_test_dfill_sp(layout, 'U', N, ap),
                                 lapacke_test_dfill_sp(layout, 'U', N, afp),
                                 lapacke_test_fill_ipiv(LD * LD, ipiv),
                                 lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
                                 lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
                                 lapacke_test_dfill_vec(LD * LD, ferr),
                                 lapacke_test_dfill_vec(LD * LD, berr)),
                                API_SUFFIX(LAPACKE_dsprfs)(layout, 'U', N, NRHS,
                                                           ap, afp, ipiv, b, LD,
                                                           x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP("dsprfs ap uplo=U", l, N * (N + 1) / 2, 1, ap,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -5,
                                (lapacke_test_dfill_sp(layout, 'U', N, ap),
                                 lapacke_test_dfill_sp(layout, 'U', N, afp),
                                 lapacke_test_fill_ipiv(LD * LD, ipiv),
                                 lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
                                 lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
                                 lapacke_test_dfill_vec(LD * LD, ferr),
                                 lapacke_test_dfill_vec(LD * LD, berr)),
                                API_SUFFIX(LAPACKE_dsprfs)(layout, 'U', N, NRHS,
                                                           ap, afp, ipiv, b, LD,
                                                           x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dsprfs b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_dfill_sp(layout, 'U', N, ap),
             lapacke_test_dfill_sp(layout, 'U', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dsprfs)(layout, 'U', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dsprfs x uplo=U", l, N, NRHS, x, LD, lapacke_test_region_full, -10,
            (lapacke_test_dfill_sp(layout, 'U', N, ap),
             lapacke_test_dfill_sp(layout, 'U', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dsprfs)(layout, 'U', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP("dsprfs afp uplo=L", l, N * (N + 1) / 2, 1, afp,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -6,
                                (lapacke_test_dfill_sp(layout, 'L', N, ap),
                                 lapacke_test_dfill_sp(layout, 'L', N, afp),
                                 lapacke_test_fill_ipiv(LD * LD, ipiv),
                                 lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
                                 lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
                                 lapacke_test_dfill_vec(LD * LD, ferr),
                                 lapacke_test_dfill_vec(LD * LD, berr)),
                                API_SUFFIX(LAPACKE_dsprfs)(layout, 'L', N, NRHS,
                                                           ap, afp, ipiv, b, LD,
                                                           x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP("dsprfs ap uplo=L", l, N * (N + 1) / 2, 1, ap,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -5,
                                (lapacke_test_dfill_sp(layout, 'L', N, ap),
                                 lapacke_test_dfill_sp(layout, 'L', N, afp),
                                 lapacke_test_fill_ipiv(LD * LD, ipiv),
                                 lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
                                 lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
                                 lapacke_test_dfill_vec(LD * LD, ferr),
                                 lapacke_test_dfill_vec(LD * LD, berr)),
                                API_SUFFIX(LAPACKE_dsprfs)(layout, 'L', N, NRHS,
                                                           ap, afp, ipiv, b, LD,
                                                           x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dsprfs b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_dfill_sp(layout, 'L', N, ap),
             lapacke_test_dfill_sp(layout, 'L', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dsprfs)(layout, 'L', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dsprfs x uplo=L", l, N, NRHS, x, LD, lapacke_test_region_full, -10,
            (lapacke_test_dfill_sp(layout, 'L', N, ap),
             lapacke_test_dfill_sp(layout, 'L', N, afp),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dsprfs)(layout, 'L', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_dfill_nan(layout, N * (N + 1) / 2, 1, afp,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_dfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "dsprfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dsprfs)(layout, 'U', N, NRHS, ap, afp, ipiv, b,
                                       LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DSPRFS_ALLOC_TEST(0, 0, "dsprfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DSPRFS_ALLOC_TEST(0, 1, "dsprfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DSPRFS_ALLOC_TEST(0, 2, "dsprfs allocation count", 0);
    lapacke_test_check_alloc_count("dsprfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DSPRFS_ALLOC_TEST(1, 0, "dsprfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DSPRFS_ALLOC_TEST(1, 1, "dsprfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DSPRFS_ALLOC_TEST(1, 2, "dsprfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DSPRFS_ALLOC_TEST(1, 3, "dsprfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DSPRFS_ALLOC_TEST(1, 4, "dsprfs transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DSPRFS_ALLOC_TEST(1, 5, "dsprfs transpose alloc failure (afp_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DSPRFS_ALLOC_TEST(1, 6, "dsprfs allocation count", 0);
    lapacke_test_check_alloc_count("dsprfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DSPRFS_ALLOC_TEST(2, 0, "dsprfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dsprfs invalid layout allocation count");
}
