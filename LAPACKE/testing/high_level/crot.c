#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call crot. */
#define LAPACKE_CROT_ALLOC_TEST(layout_index, countdown, name, expected)       \
    do {                                                                       \
        lapacke_test_cfill_vec(LD * LD, cx);                                   \
        lapacke_test_cfill_vec(LD * LD, cy);                                   \
        c[0] = 0.5f;                                                           \
        s[0] = lapack_make_complex_float(0.5f, 0.0f);                          \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_crot)(N, cx, 1, cy, 1, c[0], s[0]), expected);  \
    } while (0)

LAPACKE_TEST(crot)
{
    lapack_complex_float cx[LD * LD];
    lapack_complex_float cy[LD * LD];
    float c[1];
    lapack_complex_float s[1];

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "crot cx", l, N, 1, cx, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -2,
            (lapacke_test_cfill_vec(LD * LD, cx),
             lapacke_test_cfill_vec(LD * LD, cy), (c[0] = 0.5f),
             (s[0] = lapack_make_complex_float(0.5f, 0.0f))),
            API_SUFFIX(LAPACKE_crot)(N, cx, 1, cy, 1, c[0], s[0]));

        LAPACKE_TEST_CNAN_SWEEP(
            "crot cy", l, N, 1, cy, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -4,
            (lapacke_test_cfill_vec(LD * LD, cx),
             lapacke_test_cfill_vec(LD * LD, cy), (c[0] = 0.5f),
             (s[0] = lapack_make_complex_float(0.5f, 0.0f))),
            API_SUFFIX(LAPACKE_crot)(N, cx, 1, cy, 1, c[0], s[0]));

        LAPACKE_TEST_SNAN_SWEEP(
            "crot c", l, 1, 1, c, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_cfill_vec(LD * LD, cx),
             lapacke_test_cfill_vec(LD * LD, cy), (c[0] = 0.5f),
             (s[0] = lapack_make_complex_float(0.5f, 0.0f))),
            API_SUFFIX(LAPACKE_crot)(N, cx, 1, cy, 1, c[0], s[0]));

        LAPACKE_TEST_CNAN_SWEEP(
            "crot s", l, 1, 1, s, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -7,
            (lapacke_test_cfill_vec(LD * LD, cx),
             lapacke_test_cfill_vec(LD * LD, cy), (c[0] = 0.5f),
             (s[0] = lapack_make_complex_float(0.5f, 0.0f))),
            API_SUFFIX(LAPACKE_crot)(N, cx, 1, cy, 1, c[0], s[0]));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, N, 1, cx, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_cfill_nan(layout, N, 1, cy, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_sfill_nan(layout, 1, 1, c, LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_cfill_nan(layout, 1, 1, s, LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_check(
            "crot NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_crot)(N, cx, 1, cy, 1, c[0], s[0]) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CROT_ALLOC_TEST(0, 0, "crot allocation count", 0);
    lapacke_test_check_alloc_count("crot col-major allocation count");
}
