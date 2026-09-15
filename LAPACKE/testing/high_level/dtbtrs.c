#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_U_UNIT(i, j)                                                 \
    (lapacke_test_region_band(i, j, N, KD) && (i) != KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)
#define REGION_AB_L_UNIT(i, j)                                                 \
    (lapacke_test_region_band(i, j, N, 0) && (i) != 0)

/* Refill the inputs, schedule the malloc failure, call dtbtrs. */
#define LAPACKE_DTBTRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_tb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dtbtrs)(layout, 'U', 'N', 'N',   \
                                                      N, KD, NRHS, ab, LD, b,  \
                                                      LD),                     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dtbtrs)
{
    double ab[LD * LD];
    double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dtbtrs ab uplo=U diag=N", l, KD + 1, N, ab, LD, REGION_AB_U, -8,
            (lapacke_test_dfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dtbtrs)(layout, 'U', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtbtrs b uplo=U diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_dfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dtbtrs)(layout, 'U', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtbtrs ab uplo=U diag=U", l, KD + 1, N, ab, LD, REGION_AB_U_UNIT,
            -8,
            (lapacke_test_dfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dtbtrs)(layout, 'U', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtbtrs b uplo=U diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_dfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dtbtrs)(layout, 'U', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtbtrs ab uplo=L diag=N", l, KD + 1, N, ab, LD, REGION_AB_L, -8,
            (lapacke_test_dfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dtbtrs)(layout, 'L', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtbtrs b uplo=L diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_dfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dtbtrs)(layout, 'L', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtbtrs ab uplo=L diag=U", l, KD + 1, N, ab, LD, REGION_AB_L_UNIT,
            -8,
            (lapacke_test_dfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dtbtrs)(layout, 'L', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtbtrs b uplo=L diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_dfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dtbtrs)(layout, 'L', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check(
            "dtbtrs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dtbtrs)(layout, 'U', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_DTBTRS_ALLOC_TEST(0, 0, "dtbtrs allocation count", 0);
    lapacke_test_check_alloc_count("dtbtrs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DTBTRS_ALLOC_TEST(1, 0, "dtbtrs transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DTBTRS_ALLOC_TEST(1, 1, "dtbtrs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DTBTRS_ALLOC_TEST(1, 2, "dtbtrs allocation count", 0);
    lapacke_test_check_alloc_count("dtbtrs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DTBTRS_ALLOC_TEST(2, 0, "dtbtrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dtbtrs invalid layout allocation count");
}
