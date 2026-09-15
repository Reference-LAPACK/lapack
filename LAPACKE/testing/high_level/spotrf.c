#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call spotrf. */
#define LAPACKE_SPOTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_spd(layout, N, a, LD);                              \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_spotrf)(layout, 'U', N, a, LD),  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(spotrf)
{
    float a[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "spotrf a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -4,
            (lapacke_test_sfill_spd(layout, N, a, LD)),
            API_SUFFIX(LAPACKE_spotrf)(layout, 'U', N, a, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "spotrf a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -4,
            (lapacke_test_sfill_spd(layout, N, a, LD)),
            API_SUFFIX(LAPACKE_spotrf)(layout, 'L', N, a, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "spotrf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_spotrf)(layout, 'U', N, a, LD) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_SPOTRF_ALLOC_TEST(0, 0, "spotrf allocation count", 0);
    lapacke_test_check_alloc_count("spotrf col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPOTRF_ALLOC_TEST(1, 0, "spotrf transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPOTRF_ALLOC_TEST(1, 1, "spotrf allocation count", 0);
    lapacke_test_check_alloc_count("spotrf row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPOTRF_ALLOC_TEST(2, 0, "spotrf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("spotrf invalid layout allocation count");
}
