#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AFB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)
#define REGION_AFB_L(i, j) lapacke_test_region_band(i, j, N, 0)

/* Refill the inputs, schedule the malloc failure, call cpbrfs. */
#define LAPACKE_CPBRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_cfill_pb(layout, 'U', N, KD, afb, LD);                    \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_cfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,  \
                                       LD, b, LD, x, LD, ferr, berr),          \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(cpbrfs)
{
    lapack_complex_float ab[LD * LD];
    lapack_complex_float afb[LD * LD];
    lapack_complex_float b[LD * LD];
    lapack_complex_float x[LD * LD];
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbrfs ab uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -6,
            (lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_cfill_pb(layout, 'U', N, KD, afb, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbrfs afb uplo=U", l, KD + 1, N, afb, LD, REGION_AFB_U, -8,
            (lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_cfill_pb(layout, 'U', N, KD, afb, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbrfs b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_cfill_pb(layout, 'U', N, KD, afb, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbrfs x uplo=U", l, N, NRHS, x, LD, lapacke_test_region_full, -12,
            (lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_cfill_pb(layout, 'U', N, KD, afb, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbrfs ab uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -6,
            (lapacke_test_cfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_cfill_pb(layout, 'L', N, KD, afb, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'L', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbrfs afb uplo=L", l, KD + 1, N, afb, LD, REGION_AFB_L, -8,
            (lapacke_test_cfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_cfill_pb(layout, 'L', N, KD, afb, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'L', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbrfs b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_cfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_cfill_pb(layout, 'L', N, KD, afb, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'L', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbrfs x uplo=L", l, N, NRHS, x, LD, lapacke_test_region_full, -12,
            (lapacke_test_cfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_cfill_pb(layout, 'L', N, KD, afb, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'L', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_cfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_cfill_nan(layout, KD + 1, N, afb, LD);
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_cfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "cpbrfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cpbrfs)(layout, 'U', N, KD, NRHS, ab, LD, afb,
                                       LD, b, LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CPBRFS_ALLOC_TEST(0, 0, "cpbrfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPBRFS_ALLOC_TEST(0, 1, "cpbrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPBRFS_ALLOC_TEST(0, 2, "cpbrfs allocation count", 0);
    lapacke_test_check_alloc_count("cpbrfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CPBRFS_ALLOC_TEST(1, 0, "cpbrfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPBRFS_ALLOC_TEST(1, 1, "cpbrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CPBRFS_ALLOC_TEST(1, 2, "cpbrfs transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPBRFS_ALLOC_TEST(1, 3, "cpbrfs transpose alloc failure (afb_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPBRFS_ALLOC_TEST(1, 4, "cpbrfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPBRFS_ALLOC_TEST(1, 5, "cpbrfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPBRFS_ALLOC_TEST(1, 6, "cpbrfs allocation count", 0);
    lapacke_test_check_alloc_count("cpbrfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CPBRFS_ALLOC_TEST(2, 0, "cpbrfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cpbrfs invalid layout allocation count");
}
