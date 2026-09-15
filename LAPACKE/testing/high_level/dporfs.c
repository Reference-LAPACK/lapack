#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dporfs. */
#define LAPACKE_DPORFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_spd(layout, N, a, LD);                              \
        lapacke_test_dfill_spd(layout, N, af, LD);                             \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_dfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dporfs)(layout, 'U', N, NRHS, a, \
                                                      LD, af, LD, b, LD, x,    \
                                                      LD, ferr, berr),         \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dporfs)
{
    double a[LD * LD];
    double af[LD * LD];
    double b[LD * LD];
    double x[LD * LD];
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dporfs a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -5,
            (lapacke_test_dfill_spd(layout, N, a, LD),
             lapacke_test_dfill_spd(layout, N, af, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dporfs)(layout, 'U', N, NRHS, a, LD, af, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dporfs af uplo=U", l, N, N, af, LD, lapacke_test_region_upper, -7,
            (lapacke_test_dfill_spd(layout, N, a, LD),
             lapacke_test_dfill_spd(layout, N, af, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dporfs)(layout, 'U', N, NRHS, a, LD, af, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dporfs b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -9,
            (lapacke_test_dfill_spd(layout, N, a, LD),
             lapacke_test_dfill_spd(layout, N, af, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dporfs)(layout, 'U', N, NRHS, a, LD, af, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dporfs x uplo=U", l, N, NRHS, x, LD, lapacke_test_region_full, -11,
            (lapacke_test_dfill_spd(layout, N, a, LD),
             lapacke_test_dfill_spd(layout, N, af, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dporfs)(layout, 'U', N, NRHS, a, LD, af, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dporfs a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -5,
            (lapacke_test_dfill_spd(layout, N, a, LD),
             lapacke_test_dfill_spd(layout, N, af, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dporfs)(layout, 'L', N, NRHS, a, LD, af, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dporfs af uplo=L", l, N, N, af, LD, lapacke_test_region_lower, -7,
            (lapacke_test_dfill_spd(layout, N, a, LD),
             lapacke_test_dfill_spd(layout, N, af, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dporfs)(layout, 'L', N, NRHS, a, LD, af, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dporfs b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -9,
            (lapacke_test_dfill_spd(layout, N, a, LD),
             lapacke_test_dfill_spd(layout, N, af, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dporfs)(layout, 'L', N, NRHS, a, LD, af, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dporfs x uplo=L", l, N, NRHS, x, LD, lapacke_test_region_full, -11,
            (lapacke_test_dfill_spd(layout, N, a, LD),
             lapacke_test_dfill_spd(layout, N, af, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dporfs)(layout, 'L', N, NRHS, a, LD, af, LD, b,
                                       LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_dfill_nan(layout, N, N, a, LD);
        lapacke_test_dfill_nan(layout, N, N, af, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "dporfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dporfs)(layout, 'U', N, NRHS, a, LD, af, LD, b,
                                       LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DPORFS_ALLOC_TEST(0, 0, "dporfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPORFS_ALLOC_TEST(0, 1, "dporfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPORFS_ALLOC_TEST(0, 2, "dporfs allocation count", 0);
    lapacke_test_check_alloc_count("dporfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DPORFS_ALLOC_TEST(1, 0, "dporfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPORFS_ALLOC_TEST(1, 1, "dporfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DPORFS_ALLOC_TEST(1, 2, "dporfs transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPORFS_ALLOC_TEST(1, 3, "dporfs transpose alloc failure (af_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPORFS_ALLOC_TEST(1, 4, "dporfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPORFS_ALLOC_TEST(1, 5, "dporfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPORFS_ALLOC_TEST(1, 6, "dporfs allocation count", 0);
    lapacke_test_check_alloc_count("dporfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DPORFS_ALLOC_TEST(2, 0, "dporfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dporfs invalid layout allocation count");
}
