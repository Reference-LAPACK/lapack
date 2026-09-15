#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zlacgv. */
#define LAPACKE_ZLACGV_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_zfill_vec(LD * LD, x);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zlacgv)(N, x, 1), expected);     \
    } while (0)

LAPACKE_TEST(zlacgv)
{
    lapack_complex_double x[LD * LD];

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zlacgv x", l, N, 1, x, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -2, (lapacke_test_zfill_vec(LD * LD, x)),
            API_SUFFIX(LAPACKE_zlacgv)(N, x, 1));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, N, 1, x, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_check("zlacgv NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_zlacgv)(N, x, 1) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZLACGV_ALLOC_TEST(0, 0, "zlacgv allocation count", 0);
    lapacke_test_check_alloc_count("zlacgv col-major allocation count");
}
