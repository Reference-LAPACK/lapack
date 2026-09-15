#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

#define REGION_AB(i, j) lapacke_test_region_band(i, j, N, KU)
#define REGION_AFB(i, j) lapacke_test_region_band(i, j, N, KL + KU)

/* Refill the inputs, schedule the malloc failure, call dgbrfs. */
#define LAPACKE_DGBRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD);      \
        lapacke_test_dfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb, LD); \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_dfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dgbrfs)(                         \
                               layout, 'N', N, KL, KU, NRHS, ab, LD, afb, LD,  \
                               ipiv, b, LD, x, LD, ferr, berr),                \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dgbrfs)
{
    double ab[LD * LD];
    double afb[LD * LD];
    lapack_int ipiv[LD * LD];
    double b[LD * LD];
    double x[LD * LD];
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgbrfs ab", l, KL + KU + 1, N, ab, LD, REGION_AB, -7,
            (lapacke_test_dfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_dfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgbrfs)(layout, 'N', N, KL, KU, NRHS, ab, LD,
                                       afb, LD, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgbrfs afb", l, 2 * KL + KU + 1, N, afb, LD, REGION_AFB, -9,
            (lapacke_test_dfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_dfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgbrfs)(layout, 'N', N, KL, KU, NRHS, ab, LD,
                                       afb, LD, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgbrfs b", l, N, NRHS, b, LD, lapacke_test_region_full, -12,
            (lapacke_test_dfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_dfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgbrfs)(layout, 'N', N, KL, KU, NRHS, ab, LD,
                                       afb, LD, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgbrfs x", l, N, NRHS, x, LD, lapacke_test_region_full, -14,
            (lapacke_test_dfill_gb(layout, N, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_dfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, afb,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgbrfs)(layout, 'N', N, KL, KU, NRHS, ab, LD,
                                       afb, LD, ipiv, b, LD, x, LD, ferr,
                                       berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_dfill_nan(layout, KL + KU + 1, N, ab, LD);
        lapacke_test_dfill_nan(layout, 2 * KL + KU + 1, N, afb, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "dgbrfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgbrfs)(layout, 'N', N, KL, KU, NRHS, ab, LD,
                                       afb, LD, ipiv, b, LD, x, LD, ferr,
                                       berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DGBRFS_ALLOC_TEST(0, 0, "dgbrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGBRFS_ALLOC_TEST(0, 1, "dgbrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGBRFS_ALLOC_TEST(0, 2, "dgbrfs allocation count", 0);
    lapacke_test_check_alloc_count("dgbrfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGBRFS_ALLOC_TEST(1, 0, "dgbrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGBRFS_ALLOC_TEST(1, 1, "dgbrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGBRFS_ALLOC_TEST(1, 2, "dgbrfs transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGBRFS_ALLOC_TEST(1, 3, "dgbrfs transpose alloc failure (afb_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGBRFS_ALLOC_TEST(1, 4, "dgbrfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGBRFS_ALLOC_TEST(1, 5, "dgbrfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGBRFS_ALLOC_TEST(1, 6, "dgbrfs allocation count", 0);
    lapacke_test_check_alloc_count("dgbrfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGBRFS_ALLOC_TEST(2, 0, "dgbrfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgbrfs invalid layout allocation count");
}
