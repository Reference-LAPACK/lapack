#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zrot. */
#define LAPACKE_ZROT_ALLOC_TEST(layout_index, countdown, name, expected)       \
    do {                                                                       \
        lapacke_test_zfill_vec(LD * LD, cx);                                   \
        lapacke_test_zfill_vec(LD * LD, cy);                                   \
        c[0] = 0.5;                                                            \
        s[0] = lapack_make_complex_double(0.5, 0.0);                           \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zrot)(N, cx, 1, cy, 1, c[0], s[0]), expected);  \
    } while (0)

LAPACKE_TEST(zrot)
{
    lapack_complex_double cx[LD * LD];
    lapack_complex_double cy[LD * LD];
    double c[1];
    lapack_complex_double s[1];

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zrot cx", l, N, 1, cx, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -2,
            (lapacke_test_zfill_vec(LD * LD, cx),
             lapacke_test_zfill_vec(LD * LD, cy), (c[0] = 0.5),
             (s[0] = lapack_make_complex_double(0.5, 0.0))),
            API_SUFFIX(LAPACKE_zrot)(N, cx, 1, cy, 1, c[0], s[0]));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zrot cy", l, N, 1, cy, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -4,
            (lapacke_test_zfill_vec(LD * LD, cx),
             lapacke_test_zfill_vec(LD * LD, cy), (c[0] = 0.5),
             (s[0] = lapack_make_complex_double(0.5, 0.0))),
            API_SUFFIX(LAPACKE_zrot)(N, cx, 1, cy, 1, c[0], s[0]));

        LAPACKE_TEST_DNAN_SWEEP(
            "zrot c", l, 1, 1, c, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            (lapacke_test_zfill_vec(LD * LD, cx),
             lapacke_test_zfill_vec(LD * LD, cy), (c[0] = 0.5),
             (s[0] = lapack_make_complex_double(0.5, 0.0))),
            API_SUFFIX(LAPACKE_zrot)(N, cx, 1, cy, 1, c[0], s[0]));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zrot s", l, 1, 1, s, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -7,
            (lapacke_test_zfill_vec(LD * LD, cx),
             lapacke_test_zfill_vec(LD * LD, cy), (c[0] = 0.5),
             (s[0] = lapack_make_complex_double(0.5, 0.0))),
            API_SUFFIX(LAPACKE_zrot)(N, cx, 1, cy, 1, c[0], s[0]));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, N, 1, cx, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_zfill_nan(layout, N, 1, cy, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_dfill_nan(layout, 1, 1, c, LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_zfill_nan(layout, 1, 1, s, LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_check(
            "zrot NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zrot)(N, cx, 1, cy, 1, c[0], s[0]) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZROT_ALLOC_TEST(0, 0, "zrot allocation count", 0);
    lapacke_test_check_alloc_count("zrot col-major allocation count");
}
