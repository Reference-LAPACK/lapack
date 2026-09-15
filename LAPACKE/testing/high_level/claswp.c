#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call claswp. */
#define LAPACKE_CLASWP_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, N, N, a, LD);                               \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_claswp)(layout, N, a, LD, 1, M, ipiv, 1),       \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(claswp)
{
    lapack_complex_float a[LD * LD];
    lapack_int ipiv[LD * LD];

    /* column-major: no allocation at all */
    LAPACKE_CLASWP_ALLOC_TEST(0, 0, "claswp allocation count", 0);
    lapacke_test_check_alloc_count("claswp col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CLASWP_ALLOC_TEST(1, 0, "claswp transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CLASWP_ALLOC_TEST(1, 1, "claswp allocation count", 0);
    lapacke_test_check_alloc_count("claswp row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CLASWP_ALLOC_TEST(2, 0, "claswp invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("claswp invalid layout allocation count");
}
