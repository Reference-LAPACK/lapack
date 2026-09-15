#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zpotrf. */
#define LAPACKE_ZPOTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_spd(layout, N, a, LD);                              \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zpotrf)(layout, 'U', N, a, LD),  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zpotrf)
{
    lapack_complex_double a[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpotrf a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -4,
            (lapacke_test_zfill_spd(layout, N, a, LD)),
            API_SUFFIX(LAPACKE_zpotrf)(layout, 'U', N, a, LD));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpotrf a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -4,
            (lapacke_test_zfill_spd(layout, N, a, LD)),
            API_SUFFIX(LAPACKE_zpotrf)(layout, 'L', N, a, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "zpotrf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zpotrf)(layout, 'U', N, a, LD) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZPOTRF_ALLOC_TEST(0, 0, "zpotrf allocation count", 0);
    lapacke_test_check_alloc_count("zpotrf col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPOTRF_ALLOC_TEST(1, 0, "zpotrf transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPOTRF_ALLOC_TEST(1, 1, "zpotrf allocation count", 0);
    lapacke_test_check_alloc_count("zpotrf row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPOTRF_ALLOC_TEST(2, 0, "zpotrf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zpotrf invalid layout allocation count");
}
