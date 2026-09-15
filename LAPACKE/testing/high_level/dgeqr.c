#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define TSIZE (LD * LD)

/* Refill the inputs, schedule the malloc failure, call dgeqr. */
#define LAPACKE_DGEQR_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, M, N, a, LD);                               \
        lapacke_test_dfill_vec(LD * LD, t);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dgeqr)(layout, M, N, a, LD, t, TSIZE),          \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(dgeqr)
{
    double a[LD * LD];
    double t[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgeqr a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_dfill(layout, M, N, a, LD),
             lapacke_test_dfill_vec(LD * LD, t)),
            API_SUFFIX(LAPACKE_dgeqr)(layout, M, N, a, LD, t, TSIZE));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_vec(LD * LD, t);
        lapacke_test_dfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "dgeqr NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgeqr)(layout, M, N, a, LD, t, TSIZE) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DGEQR_ALLOC_TEST(0, 0, "dgeqr work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGEQR_ALLOC_TEST(0, 1, "dgeqr allocation count", 0);
    lapacke_test_check_alloc_count("dgeqr col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGEQR_ALLOC_TEST(1, 0, "dgeqr work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGEQR_ALLOC_TEST(1, 1, "dgeqr transpose alloc failure (a_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGEQR_ALLOC_TEST(1, 2, "dgeqr allocation count", 0);
    lapacke_test_check_alloc_count("dgeqr row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGEQR_ALLOC_TEST(2, 0, "dgeqr invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgeqr invalid layout allocation count");
}
