#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zpoequ. */
#define LAPACKE_ZPOEQU_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, N, N, a, LD);                               \
        lapacke_test_dfill_pos(LD * LD, s);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zpoequ)(layout, N, a, LD, s, &scond, &amax),    \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(zpoequ)
{
    lapack_complex_double a[LD * LD];
    double s[LD * LD];
    double scond;
    double amax;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpoequ a", l, N, N, a, LD, lapacke_test_region_full, -3,
            (lapacke_test_zfill(layout, N, N, a, LD),
             lapacke_test_dfill_pos(LD * LD, s)),
            API_SUFFIX(LAPACKE_zpoequ)(layout, N, a, LD, s, &scond, &amax));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_pos(LD * LD, s);
        lapacke_test_zfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "zpoequ NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zpoequ)(layout, N, a, LD, s, &scond, &amax) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZPOEQU_ALLOC_TEST(0, 0, "zpoequ allocation count", 0);
    lapacke_test_check_alloc_count("zpoequ col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPOEQU_ALLOC_TEST(1, 0, "zpoequ transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPOEQU_ALLOC_TEST(1, 1, "zpoequ allocation count", 0);
    lapacke_test_check_alloc_count("zpoequ row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPOEQU_ALLOC_TEST(2, 0, "zpoequ invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zpoequ invalid layout allocation count");
}
