#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AFB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)
#define REGION_AFB_L(i, j) lapacke_test_region_band(i, j, N, 0)

/* Refill the inputs, schedule the malloc failure, call dpbrfs. */
#define LAPACKE_DPBRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_pb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_dfill_pb(layout, 'U', N, KD, afb, LD);                    \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_dfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,  \
                                       LD, b, LD, x, LD, ferr, berr),          \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(dpbrfs)
{
    double ab[LD * LD];
    double afb[LD * LD];
    double b[LD * LD];
    double x[LD * LD];
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dpbrfs ab uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -6,
            (lapacke_test_dfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_dfill_pb(layout, 'U', N, KD, afb, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpbrfs afb uplo=U", l, KD + 1, N, afb, LD, REGION_AFB_U, -8,
            (lapacke_test_dfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_dfill_pb(layout, 'U', N, KD, afb, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpbrfs b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_dfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_dfill_pb(layout, 'U', N, KD, afb, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpbrfs x uplo=U", l, N, NRHS, x, LD, lapacke_test_region_full, -12,
            (lapacke_test_dfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_dfill_pb(layout, 'U', N, KD, afb, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpbrfs ab uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -6,
            (lapacke_test_dfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_dfill_pb(layout, 'L', N, KD, afb, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'L', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpbrfs afb uplo=L", l, KD + 1, N, afb, LD, REGION_AFB_L, -8,
            (lapacke_test_dfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_dfill_pb(layout, 'L', N, KD, afb, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'L', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpbrfs b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_dfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_dfill_pb(layout, 'L', N, KD, afb, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'L', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpbrfs x uplo=L", l, N, NRHS, x, LD, lapacke_test_region_full, -12,
            (lapacke_test_dfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_dfill_pb(layout, 'L', N, KD, afb, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'L', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_dfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_dfill_nan(layout, KD + 1, N, afb, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "dpbrfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DPBRFS_ALLOC_TEST(0, 0, "dpbrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPBRFS_ALLOC_TEST(0, 1, "dpbrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPBRFS_ALLOC_TEST(0, 2, "dpbrfs allocation count", 0);
    lapacke_test_check_alloc_count("dpbrfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DPBRFS_ALLOC_TEST(1, 0, "dpbrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPBRFS_ALLOC_TEST(1, 1, "dpbrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPBRFS_ALLOC_TEST(1, 2, "dpbrfs transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPBRFS_ALLOC_TEST(1, 3, "dpbrfs transpose alloc failure (afb_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPBRFS_ALLOC_TEST(1, 4, "dpbrfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPBRFS_ALLOC_TEST(1, 5, "dpbrfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPBRFS_ALLOC_TEST(1, 6, "dpbrfs allocation count", 0);
    lapacke_test_check_alloc_count("dpbrfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DPBRFS_ALLOC_TEST(2, 0, "dpbrfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dpbrfs invalid layout allocation count");
}
