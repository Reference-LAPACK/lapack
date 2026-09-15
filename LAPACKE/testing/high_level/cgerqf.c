#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cgerqf. */
#define LAPACKE_CGERQF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, M, N, a, LD);                               \
        lapacke_test_cfill_vec(LD * LD, tau);                                  \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_cgerqf)(layout, M, N, a, LD, tau), expected);   \
    } while (0)

LAPACKE_TEST(cgerqf)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float tau[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cgerqf a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_cfill(layout, M, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, tau)),
            API_SUFFIX(LAPACKE_cgerqf)(layout, M, N, a, LD, tau));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_vec(LD * LD, tau);
        lapacke_test_cfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "cgerqf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgerqf)(layout, M, N, a, LD, tau) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CGERQF_ALLOC_TEST(0, 0, "cgerqf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGERQF_ALLOC_TEST(0, 1, "cgerqf allocation count", 0);
    lapacke_test_check_alloc_count("cgerqf col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CGERQF_ALLOC_TEST(1, 0, "cgerqf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGERQF_ALLOC_TEST(1, 1, "cgerqf transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGERQF_ALLOC_TEST(1, 2, "cgerqf allocation count", 0);
    lapacke_test_check_alloc_count("cgerqf row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CGERQF_ALLOC_TEST(2, 0, "cgerqf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cgerqf invalid layout allocation count");
}
