#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call slaswp. */
#define LAPACKE_SLASWP_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, N, N, a, LD);                               \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_slaswp)(layout, N, a, LD, 1, M, ipiv, 1),       \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(slaswp)
{
    float a[LD * LD];
    lapack_int ipiv[LD * LD];

    /* column-major: no allocation at all */
    LAPACKE_SLASWP_ALLOC_TEST(0, 0, "slaswp allocation count", 0);
    lapacke_test_check_alloc_count("slaswp col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SLASWP_ALLOC_TEST(1, 0, "slaswp transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SLASWP_ALLOC_TEST(1, 1, "slaswp allocation count", 0);
    lapacke_test_check_alloc_count("slaswp row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SLASWP_ALLOC_TEST(2, 0, "slaswp invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("slaswp invalid layout allocation count");
}
