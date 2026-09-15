#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sptsv. */
#define LAPACKE_SPTSV_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_pos(LD * LD, d);                                    \
        lapacke_test_sfill_vec(LD * LD, e);                                    \
        lapacke_test_sfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_sptsv)(layout, N, NRHS, d, e, b, LD),           \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(sptsv)
{
    float d[LD * LD];
    float e[LD * LD];
    float b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsv b", l, N, NRHS, b, LD, lapacke_test_region_full, -6,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_sptsv)(layout, N, NRHS, d, e, b, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsv d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -4,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_sptsv)(layout, N, NRHS, d, e, b, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsv e", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -5,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_sptsv)(layout, N, NRHS, d, e, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_sfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_sfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check(
            "sptsv NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sptsv)(layout, N, NRHS, d, e, b, LD) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_SPTSV_ALLOC_TEST(0, 0, "sptsv allocation count", 0);
    lapacke_test_check_alloc_count("sptsv col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPTSV_ALLOC_TEST(1, 0, "sptsv transpose alloc failure (b_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPTSV_ALLOC_TEST(1, 1, "sptsv allocation count", 0);
    lapacke_test_check_alloc_count("sptsv row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPTSV_ALLOC_TEST(2, 0, "sptsv invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sptsv invalid layout allocation count");
}
