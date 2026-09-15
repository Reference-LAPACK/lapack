#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define NB 2

/* Refill the inputs, schedule the malloc failure, call dgeqrt. */
#define LAPACKE_DGEQRT_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, M, N, a, LD);                               \
        lapacke_test_dfill_rhs(layout, NB, N, t, LD);                          \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dgeqrt)(layout, M, N, NB, a, LD, t, LD),        \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(dgeqrt)
{
    double a[LD * LD];
    double t[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgeqrt a", l, M, N, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_dfill(layout, M, N, a, LD),
             lapacke_test_dfill_rhs(layout, NB, N, t, LD)),
            API_SUFFIX(LAPACKE_dgeqrt)(layout, M, N, NB, a, LD, t, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_rhs(layout, NB, N, t, LD);
        lapacke_test_dfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "dgeqrt NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgeqrt)(layout, M, N, NB, a, LD, t, LD) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DGEQRT_ALLOC_TEST(0, 0, "dgeqrt work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGEQRT_ALLOC_TEST(0, 1, "dgeqrt allocation count", 0);
    lapacke_test_check_alloc_count("dgeqrt col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGEQRT_ALLOC_TEST(1, 0, "dgeqrt work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGEQRT_ALLOC_TEST(1, 1, "dgeqrt transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGEQRT_ALLOC_TEST(1, 2, "dgeqrt transpose alloc failure (t_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGEQRT_ALLOC_TEST(1, 3, "dgeqrt allocation count", 0);
    lapacke_test_check_alloc_count("dgeqrt row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGEQRT_ALLOC_TEST(2, 0, "dgeqrt invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgeqrt invalid layout allocation count");
}
