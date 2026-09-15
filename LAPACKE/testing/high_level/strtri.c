#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call strtri. */
#define LAPACKE_STRTRI_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_tri(layout, 'U', N, a, LD);                         \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_strtri)(layout, 'U', 'N', N, a, LD), expected); \
    } while (0)

LAPACKE_TEST(strtri)
{
    float a[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "strtri a uplo=U diag=N", l, N, N, a, LD, lapacke_test_region_upper,
            -5, (lapacke_test_sfill_tri(layout, 'U', N, a, LD)),
            API_SUFFIX(LAPACKE_strtri)(layout, 'U', 'N', N, a, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "strtri a uplo=U diag=U", l, N, N, a, LD,
            lapacke_test_region_strict_upper, -5,
            (lapacke_test_sfill_tri(layout, 'U', N, a, LD)),
            API_SUFFIX(LAPACKE_strtri)(layout, 'U', 'U', N, a, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "strtri a uplo=L diag=N", l, N, N, a, LD, lapacke_test_region_lower,
            -5, (lapacke_test_sfill_tri(layout, 'L', N, a, LD)),
            API_SUFFIX(LAPACKE_strtri)(layout, 'L', 'N', N, a, LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "strtri a uplo=L diag=U", l, N, N, a, LD,
            lapacke_test_region_strict_lower, -5,
            (lapacke_test_sfill_tri(layout, 'L', N, a, LD)),
            API_SUFFIX(LAPACKE_strtri)(layout, 'L', 'U', N, a, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "strtri NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_strtri)(layout, 'U', 'N', N, a, LD) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_STRTRI_ALLOC_TEST(0, 0, "strtri allocation count", 0);
    lapacke_test_check_alloc_count("strtri col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_STRTRI_ALLOC_TEST(1, 0, "strtri transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_STRTRI_ALLOC_TEST(1, 1, "strtri allocation count", 0);
    lapacke_test_check_alloc_count("strtri row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_STRTRI_ALLOC_TEST(2, 0, "strtri invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("strtri invalid layout allocation count");
}
