#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cgttrf. */
#define LAPACKE_CGTTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_cfill_vec(LD * LD, dl);                                   \
        lapacke_test_cfill_pos(LD * LD, d);                                    \
        lapacke_test_cfill_vec(LD * LD, du);                                   \
        lapacke_test_cfill_vec(LD * LD, du2);                                  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_cgttrf)(N, dl, d, du, du2, ipiv), expected);    \
    } while (0)

LAPACKE_TEST(cgttrf)
{
    lapack_complex_float dl[LD * LD];
    lapack_complex_float d[LD * LD];
    lapack_complex_float du[LD * LD];
    lapack_complex_float du2[LD * LD];
    lapack_int ipiv[LD * LD];

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cgttrf d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -3,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_cgttrf)(N, dl, d, du, du2, ipiv));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgttrf dl", l, N - 1, 1, dl, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -2,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_cgttrf)(N, dl, d, du, du2, ipiv));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgttrf du", l, N - 1, 1, du, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -4,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_cgttrf)(N, dl, d, du, du2, ipiv));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_vec(LD * LD, du2);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_cfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_cfill_nan(layout, N - 1, 1, dl,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_cfill_nan(layout, N - 1, 1, du,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check(
            "cgttrf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgttrf)(N, dl, d, du, du2, ipiv) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CGTTRF_ALLOC_TEST(0, 0, "cgttrf allocation count", 0);
    lapacke_test_check_alloc_count("cgttrf col-major allocation count");
}
