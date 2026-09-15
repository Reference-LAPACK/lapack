#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call spoequ. */
#define LAPACKE_SPOEQU_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, N, N, a, LD);                               \
        lapacke_test_sfill_pos(LD * LD, s);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_spoequ)(layout, N, a, LD, s, &scond, &amax),    \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(spoequ)
{
    float a[LD * LD];
    float s[LD * LD];
    float scond;
    float amax;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "spoequ a", l, N, N, a, LD, lapacke_test_region_full, -3,
            (lapacke_test_sfill(layout, N, N, a, LD),
             lapacke_test_sfill_pos(LD * LD, s)),
            API_SUFFIX(LAPACKE_spoequ)(layout, N, a, LD, s, &scond, &amax));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_pos(LD * LD, s);
        lapacke_test_sfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "spoequ NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_spoequ)(layout, N, a, LD, s, &scond, &amax) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_SPOEQU_ALLOC_TEST(0, 0, "spoequ allocation count", 0);
    lapacke_test_check_alloc_count("spoequ col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPOEQU_ALLOC_TEST(1, 0, "spoequ transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPOEQU_ALLOC_TEST(1, 1, "spoequ allocation count", 0);
    lapacke_test_check_alloc_count("spoequ row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPOEQU_ALLOC_TEST(2, 0, "spoequ invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("spoequ invalid layout allocation count");
}
