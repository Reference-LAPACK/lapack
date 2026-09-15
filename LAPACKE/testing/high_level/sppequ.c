#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sppequ. */
#define LAPACKE_SPPEQU_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_pp(layout, 'U', N, ap);                             \
        lapacke_test_sfill_pos(LD * LD, s);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_sppequ)(layout, 'U', N, ap, s, &scond, &amax),  \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(sppequ)
{
    float ap[LD * LD];
    float s[LD * LD];
    float scond;
    float amax;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sppequ ap uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_sfill_pp(layout, 'U', N, ap),
             lapacke_test_sfill_pos(LD * LD, s)),
            API_SUFFIX(LAPACKE_sppequ)(layout, 'U', N, ap, s, &scond, &amax));

        LAPACKE_TEST_SNAN_SWEEP(
            "sppequ ap uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_sfill_pp(layout, 'L', N, ap),
             lapacke_test_sfill_pos(LD * LD, s)),
            API_SUFFIX(LAPACKE_sppequ)(layout, 'L', N, ap, s, &scond, &amax));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_pos(LD * LD, s);
        lapacke_test_sfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_check("sppequ NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_sppequ)(layout, 'U', N, ap, s,
                                                      &scond, &amax) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_SPPEQU_ALLOC_TEST(0, 0, "sppequ allocation count", 0);
    lapacke_test_check_alloc_count("sppequ col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPPEQU_ALLOC_TEST(1, 0, "sppequ transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPPEQU_ALLOC_TEST(1, 1, "sppequ allocation count", 0);
    lapacke_test_check_alloc_count("sppequ row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPPEQU_ALLOC_TEST(2, 0, "sppequ invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sppequ invalid layout allocation count");
}
