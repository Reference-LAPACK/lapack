#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zgeequ. */
#define LAPACKE_ZGEEQU_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, M, N, a, LD);                               \
        lapacke_test_dfill_pos(LD * LD, r);                                    \
        lapacke_test_dfill_pos(LD * LD, c);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zgeequ)(layout, M, N, a, LD, r,  \
                                                      c, &rowcnd, &colcnd,     \
                                                      &amax),                  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zgeequ)
{
    lapack_complex_double a[LD * LD];
    double r[LD * LD];
    double c[LD * LD];
    double rowcnd;
    double colcnd;
    double amax;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgeequ a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_zfill(layout, M, N, a, LD),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c)),
            API_SUFFIX(LAPACKE_zgeequ)(layout, M, N, a, LD, r, c, &rowcnd,
                                       &colcnd, &amax));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_pos(LD * LD, r);
        lapacke_test_dfill_pos(LD * LD, c);
        lapacke_test_zfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "zgeequ NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgeequ)(layout, M, N, a, LD, r, c, &rowcnd,
                                       &colcnd, &amax) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZGEEQU_ALLOC_TEST(0, 0, "zgeequ allocation count", 0);
    lapacke_test_check_alloc_count("zgeequ col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZGEEQU_ALLOC_TEST(1, 0, "zgeequ transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGEEQU_ALLOC_TEST(1, 1, "zgeequ allocation count", 0);
    lapacke_test_check_alloc_count("zgeequ row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZGEEQU_ALLOC_TEST(2, 0, "zgeequ invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgeequ invalid layout allocation count");
}
