#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sgttrf. */
#define LAPACKE_SGTTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_sfill_vec(LD * LD, dl);                                   \
        lapacke_test_sfill_pos(LD * LD, d);                                    \
        lapacke_test_sfill_vec(LD * LD, du);                                   \
        lapacke_test_sfill_vec(LD * LD, du2);                                  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_sgttrf)(N, dl, d, du, du2, ipiv), expected);    \
    } while (0)

LAPACKE_TEST(sgttrf)
{
    float dl[LD * LD];
    float d[LD * LD];
    float du[LD * LD];
    float du2[LD * LD];
    lapack_int ipiv[LD * LD];

    for (size_t l = 0; l < 1; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgttrf d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -3,
            (lapacke_test_sfill_vec(LD * LD, dl),
             lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, du),
             lapacke_test_sfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_sgttrf)(N, dl, d, du, du2, ipiv));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgttrf dl", l, N - 1, 1, dl, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -2,
            (lapacke_test_sfill_vec(LD * LD, dl),
             lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, du),
             lapacke_test_sfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_sgttrf)(N, dl, d, du, du2, ipiv));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgttrf du", l, N - 1, 1, du, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -4,
            (lapacke_test_sfill_vec(LD * LD, dl),
             lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, du),
             lapacke_test_sfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_sgttrf)(N, dl, d, du, du2, ipiv));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_vec(LD * LD, du2);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_sfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_sfill_nan(layout, N - 1, 1, dl,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_sfill_nan(layout, N - 1, 1, du,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check(
            "sgttrf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sgttrf)(N, dl, d, du, du2, ipiv) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_SGTTRF_ALLOC_TEST(0, 0, "sgttrf allocation count", 0);
    lapacke_test_check_alloc_count("sgttrf col-major allocation count");
}
