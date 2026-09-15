#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zposv. */
#define LAPACKE_ZPOSV_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_spd(layout, N, a, LD);                              \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zposv)(layout, 'U', N, NRHS, a, LD, b, LD),     \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(zposv)
{
    lapack_complex_double a[LD * LD];
    lapack_complex_double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zposv a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -5,
            (lapacke_test_zfill_spd(layout, N, a, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zposv)(layout, 'U', N, NRHS, a, LD, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zposv b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -7,
            (lapacke_test_zfill_spd(layout, N, a, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zposv)(layout, 'U', N, NRHS, a, LD, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zposv a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -5,
            (lapacke_test_zfill_spd(layout, N, a, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zposv)(layout, 'L', N, NRHS, a, LD, b, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zposv b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -7,
            (lapacke_test_zfill_spd(layout, N, a, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_zposv)(layout, 'L', N, NRHS, a, LD, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, N, N, a, LD);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check(
            "zposv NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zposv)(layout, 'U', N, NRHS, a, LD, b, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZPOSV_ALLOC_TEST(0, 0, "zposv allocation count", 0);
    lapacke_test_check_alloc_count("zposv col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPOSV_ALLOC_TEST(1, 0, "zposv transpose alloc failure (a_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPOSV_ALLOC_TEST(1, 1, "zposv transpose alloc failure (b_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPOSV_ALLOC_TEST(1, 2, "zposv allocation count", 0);
    lapacke_test_check_alloc_count("zposv row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPOSV_ALLOC_TEST(2, 0, "zposv invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zposv invalid layout allocation count");
}
