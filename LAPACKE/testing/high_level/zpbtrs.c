#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)

/* Refill the inputs, schedule the malloc failure, call zpbtrs. */
#define LAPACKE_ZPBTRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zpbtrs)(layout, 'U', N, KD,      \
                                                      NRHS, ab, LD, b, LD),    \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zpbtrs)
{
    lapack_complex_double ab[LD * LD];
    lapack_complex_double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbtrs ab uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -6,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zpbtrs)(layout, 'U', N, KD, NRHS, ab, LD, b,
                                       LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbtrs b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zpbtrs)(layout, 'U', N, KD, NRHS, ab, LD, b,
                                       LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbtrs ab uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -6,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zpbtrs)(layout, 'L', N, KD, NRHS, ab, LD, b,
                                       LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbtrs b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zpbtrs)(layout, 'L', N, KD, NRHS, ab, LD, b,
                                       LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check("zpbtrs NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_zpbtrs)(layout, 'U', N, KD, NRHS,
                                                      ab, LD, b, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZPBTRS_ALLOC_TEST(0, 0, "zpbtrs allocation count", 0);
    lapacke_test_check_alloc_count("zpbtrs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPBTRS_ALLOC_TEST(1, 0, "zpbtrs transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPBTRS_ALLOC_TEST(1, 1, "zpbtrs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPBTRS_ALLOC_TEST(1, 2, "zpbtrs allocation count", 0);
    lapacke_test_check_alloc_count("zpbtrs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPBTRS_ALLOC_TEST(2, 0, "zpbtrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zpbtrs invalid layout allocation count");
}
