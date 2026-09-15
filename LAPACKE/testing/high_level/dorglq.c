#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define K N

/* Refill the inputs, schedule the malloc failure, call dorglq. */
#define LAPACKE_DORGLQ_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, N, M, a, LD);                               \
        lapacke_test_dfill_vec(LD * LD, tau);                                  \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dorglq)(layout, N, M, K, a, LD, tau),           \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(dorglq)
{
    double a[LD * LD];
    double tau[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dorglq a", l, N, M, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_dfill(layout, N, M, a, LD),
             lapacke_test_dfill_vec(LD * LD, tau)),
            API_SUFFIX(LAPACKE_dorglq)(layout, N, M, K, a, LD, tau));

        LAPACKE_TEST_DNAN_SWEEP(
            "dorglq tau", l, K, 1, tau, LAPACKE_TEST_VLD(layout, K),
            lapacke_test_region_full, -7,
            (lapacke_test_dfill(layout, N, M, a, LD),
             lapacke_test_dfill_vec(LD * LD, tau)),
            API_SUFFIX(LAPACKE_dorglq)(layout, N, M, K, a, LD, tau));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_nan(layout, N, M, a, LD);
        lapacke_test_dfill_nan(layout, K, 1, tau, LAPACKE_TEST_VLD(layout, K));
        lapacke_test_check(
            "dorglq NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dorglq)(layout, N, M, K, a, LD, tau) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DORGLQ_ALLOC_TEST(0, 0, "dorglq work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DORGLQ_ALLOC_TEST(0, 1, "dorglq allocation count", 0);
    lapacke_test_check_alloc_count("dorglq col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DORGLQ_ALLOC_TEST(1, 0, "dorglq work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DORGLQ_ALLOC_TEST(1, 1, "dorglq transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DORGLQ_ALLOC_TEST(1, 2, "dorglq allocation count", 0);
    lapacke_test_check_alloc_count("dorglq row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DORGLQ_ALLOC_TEST(2, 0, "dorglq invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dorglq invalid layout allocation count");
}
