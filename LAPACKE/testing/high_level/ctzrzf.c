#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call ctzrzf. */
#define LAPACKE_CTZRZF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, N, M, a, LD);                               \
        lapacke_test_cfill_vec(LD * LD, tau);                                  \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_ctzrzf)(layout, N, M, a, LD, tau), expected);   \
    } while (0)

LAPACKE_TEST(ctzrzf)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float tau[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "ctzrzf a", l, N, M, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_cfill(layout, N, M, a, LD),
             lapacke_test_cfill_vec(LD * LD, tau)),
            API_SUFFIX(LAPACKE_ctzrzf)(layout, N, M, a, LD, tau));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_vec(LD * LD, tau);
        lapacke_test_cfill_nan(layout, N, M, a, LD);
        lapacke_test_check(
            "ctzrzf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_ctzrzf)(layout, N, M, a, LD, tau) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CTZRZF_ALLOC_TEST(0, 0, "ctzrzf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTZRZF_ALLOC_TEST(0, 1, "ctzrzf allocation count", 0);
    lapacke_test_check_alloc_count("ctzrzf col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CTZRZF_ALLOC_TEST(1, 0, "ctzrzf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTZRZF_ALLOC_TEST(1, 1, "ctzrzf transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CTZRZF_ALLOC_TEST(1, 2, "ctzrzf allocation count", 0);
    lapacke_test_check_alloc_count("ctzrzf row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CTZRZF_ALLOC_TEST(2, 0, "ctzrzf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("ctzrzf invalid layout allocation count");
}
