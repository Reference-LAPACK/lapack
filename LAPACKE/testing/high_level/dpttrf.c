#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dpttrf. */
#define LAPACKE_DPTTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_dfill_pos(LD * LD, d);                                    \
        lapacke_test_dfill_vec(LD * LD, e);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dpttrf)(N, d, e), expected);     \
    } while (0)

LAPACKE_TEST(dpttrf)
{
    double d[LD * LD];
    double e[LD * LD];

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP("dpttrf d", l, N, 1, d,
                                LAPACKE_TEST_VLD(layout, N),
                                lapacke_test_region_full, -2,
                                (lapacke_test_dfill_pos(LD * LD, d),
                                 lapacke_test_dfill_vec(LD * LD, e)),
                                API_SUFFIX(LAPACKE_dpttrf)(N, d, e));

        LAPACKE_TEST_DNAN_SWEEP("dpttrf e", l, N - 1, 1, e,
                                LAPACKE_TEST_VLD(layout, N - 1),
                                lapacke_test_region_full, -3,
                                (lapacke_test_dfill_pos(LD * LD, d),
                                 lapacke_test_dfill_vec(LD * LD, e)),
                                API_SUFFIX(LAPACKE_dpttrf)(N, d, e));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_dfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check("dpttrf NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_dpttrf)(N, d, e) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_DPTTRF_ALLOC_TEST(0, 0, "dpttrf allocation count", 0);
    lapacke_test_check_alloc_count("dpttrf col-major allocation count");
}
