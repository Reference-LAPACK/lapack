#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dgeequ. */
#define LAPACKE_DGEEQU_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, M, N, a, LD);                               \
        lapacke_test_dfill_pos(LD * LD, r);                                    \
        lapacke_test_dfill_pos(LD * LD, c);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dgeequ)(layout, M, N, a, LD, r,  \
                                                      c, &rowcnd, &colcnd,     \
                                                      &amax),                  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dgeequ)
{
    double a[LD * LD];
    double r[LD * LD];
    double c[LD * LD];
    double rowcnd;
    double colcnd;
    double amax;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgeequ a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_dfill(layout, M, N, a, LD),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c)),
            API_SUFFIX(LAPACKE_dgeequ)(layout, M, N, a, LD, r, c, &rowcnd,
                                       &colcnd, &amax));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_pos(LD * LD, r);
        lapacke_test_dfill_pos(LD * LD, c);
        lapacke_test_dfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "dgeequ NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgeequ)(layout, M, N, a, LD, r, c, &rowcnd,
                                       &colcnd, &amax) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_DGEEQU_ALLOC_TEST(0, 0, "dgeequ allocation count", 0);
    lapacke_test_check_alloc_count("dgeequ col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGEEQU_ALLOC_TEST(1, 0, "dgeequ transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGEEQU_ALLOC_TEST(1, 1, "dgeequ allocation count", 0);
    lapacke_test_check_alloc_count("dgeequ row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGEEQU_ALLOC_TEST(2, 0, "dgeequ invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgeequ invalid layout allocation count");
}
