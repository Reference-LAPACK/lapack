#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dgeqlf. */
#define LAPACKE_DGEQLF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, M, N, a, LD);                               \
        lapacke_test_dfill_vec(LD * LD, tau);                                  \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dgeqlf)(layout, M, N, a, LD, tau), expected);   \
    } while (0)

LAPACKE_TEST(dgeqlf)
{
    double a[LD * LD];
    double tau[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgeqlf a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_dfill(layout, M, N, a, LD),
             lapacke_test_dfill_vec(LD * LD, tau)),
            API_SUFFIX(LAPACKE_dgeqlf)(layout, M, N, a, LD, tau));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_vec(LD * LD, tau);
        lapacke_test_dfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "dgeqlf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgeqlf)(layout, M, N, a, LD, tau) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DGEQLF_ALLOC_TEST(0, 0, "dgeqlf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGEQLF_ALLOC_TEST(0, 1, "dgeqlf allocation count", 0);
    lapacke_test_check_alloc_count("dgeqlf col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGEQLF_ALLOC_TEST(1, 0, "dgeqlf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGEQLF_ALLOC_TEST(1, 1, "dgeqlf transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGEQLF_ALLOC_TEST(1, 2, "dgeqlf allocation count", 0);
    lapacke_test_check_alloc_count("dgeqlf row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGEQLF_ALLOC_TEST(2, 0, "dgeqlf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgeqlf invalid layout allocation count");
}
