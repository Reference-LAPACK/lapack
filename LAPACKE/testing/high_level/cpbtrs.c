#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)

/* Refill the inputs, schedule the malloc failure, call cpbtrs. */
#define LAPACKE_CPBTRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cpbtrs)(layout, 'U', N, KD,      \
                                                      NRHS, ab, LD, b, LD),    \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cpbtrs)
{
    lapack_complex_float ab[LD * LD];
    lapack_complex_float b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbtrs ab uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -6,
            (lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cpbtrs)(layout, 'U', N, KD, NRHS, ab, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbtrs b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cpbtrs)(layout, 'U', N, KD, NRHS, ab, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbtrs ab uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -6,
            (lapacke_test_cfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cpbtrs)(layout, 'L', N, KD, NRHS, ab, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbtrs b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_cfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cpbtrs)(layout, 'L', N, KD, NRHS, ab, LD, b,
                                       LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check("cpbtrs NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_cpbtrs)(layout, 'U', N, KD, NRHS,
                                                      ab, LD, b, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CPBTRS_ALLOC_TEST(0, 0, "cpbtrs allocation count", 0);
    lapacke_test_check_alloc_count("cpbtrs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CPBTRS_ALLOC_TEST(1, 0, "cpbtrs transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPBTRS_ALLOC_TEST(1, 1, "cpbtrs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPBTRS_ALLOC_TEST(1, 2, "cpbtrs allocation count", 0);
    lapacke_test_check_alloc_count("cpbtrs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CPBTRS_ALLOC_TEST(2, 0, "cpbtrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cpbtrs invalid layout allocation count");
}
