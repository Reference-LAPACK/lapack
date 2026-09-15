#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define TSIZE (LD * LD)

/* Refill the inputs, schedule the malloc failure, call cgeqr. */
#define LAPACKE_CGEQR_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, M, N, a, LD);                               \
        lapacke_test_cfill_vec(LD * LD, t);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_cgeqr)(layout, M, N, a, LD, t, TSIZE),          \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(cgeqr)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float t[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cgeqr a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_cfill(layout, M, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, t)),
            API_SUFFIX(LAPACKE_cgeqr)(layout, M, N, a, LD, t, TSIZE));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_vec(LD * LD, t);
        lapacke_test_cfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "cgeqr NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgeqr)(layout, M, N, a, LD, t, TSIZE) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CGEQR_ALLOC_TEST(0, 0, "cgeqr work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGEQR_ALLOC_TEST(0, 1, "cgeqr allocation count", 0);
    lapacke_test_check_alloc_count("cgeqr col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CGEQR_ALLOC_TEST(1, 0, "cgeqr work alloc failure (work)",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGEQR_ALLOC_TEST(1, 1, "cgeqr transpose alloc failure (a_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGEQR_ALLOC_TEST(1, 2, "cgeqr allocation count", 0);
    lapacke_test_check_alloc_count("cgeqr row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CGEQR_ALLOC_TEST(2, 0, "cgeqr invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cgeqr invalid layout allocation count");
}
