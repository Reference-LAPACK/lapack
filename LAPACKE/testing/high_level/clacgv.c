#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call clacgv. */
#define LAPACKE_CLACGV_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_cfill_vec(LD * LD, x);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_clacgv)(N, x, 1), expected);     \
    } while (0)

LAPACKE_TEST(clacgv)
{
    lapack_complex_float x[LD * LD];

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "clacgv x", l, N, 1, x, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -2, (lapacke_test_cfill_vec(LD * LD, x)),
            API_SUFFIX(LAPACKE_clacgv)(N, x, 1));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, N, 1, x, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_check("clacgv NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_clacgv)(N, x, 1) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CLACGV_ALLOC_TEST(0, 0, "clacgv allocation count", 0);
    lapacke_test_check_alloc_count("clacgv col-major allocation count");
}
