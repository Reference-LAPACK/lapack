#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cpttrf. */
#define LAPACKE_CPTTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_sfill_pos(LD * LD, d);                                    \
        lapacke_test_cfill_vec(LD * LD, e);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cpttrf)(N, d, e), expected);     \
    } while (0)

LAPACKE_TEST(cpttrf)
{
    float d[LD * LD];
    lapack_complex_float e[LD * LD];

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP("cpttrf d", l, N, 1, d,
                                LAPACKE_TEST_VLD(layout, N),
                                lapacke_test_region_full, -2,
                                (lapacke_test_sfill_pos(LD * LD, d),
                                 lapacke_test_cfill_vec(LD * LD, e)),
                                API_SUFFIX(LAPACKE_cpttrf)(N, d, e));

        LAPACKE_TEST_CNAN_SWEEP("cpttrf e", l, N - 1, 1, e,
                                LAPACKE_TEST_VLD(layout, N - 1),
                                lapacke_test_region_full, -3,
                                (lapacke_test_sfill_pos(LD * LD, d),
                                 lapacke_test_cfill_vec(LD * LD, e)),
                                API_SUFFIX(LAPACKE_cpttrf)(N, d, e));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_cfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check("cpttrf NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_cpttrf)(N, d, e) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CPTTRF_ALLOC_TEST(0, 0, "cpttrf allocation count", 0);
    lapacke_test_check_alloc_count("cpttrf col-major allocation count");
}
