#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call ssytrf. */
#define LAPACKE_SSYTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_sym(layout, N, a, LD);                              \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_ssytrf)(layout, 'U', N, a, LD, ipiv),           \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(ssytrf)
{
    float a[LD * LD];
    lapack_int ipiv[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "ssytrf a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -4,
            (lapacke_test_sfill_sym(layout, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_ssytrf)(layout, 'U', N, a, LD, ipiv));

        LAPACKE_TEST_SNAN_SWEEP(
            "ssytrf a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -4,
            (lapacke_test_sfill_sym(layout, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_ssytrf)(layout, 'L', N, a, LD, ipiv));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_sfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "ssytrf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_ssytrf)(layout, 'U', N, a, LD, ipiv) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SSYTRF_ALLOC_TEST(0, 0, "ssytrf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSYTRF_ALLOC_TEST(0, 1, "ssytrf allocation count", 0);
    lapacke_test_check_alloc_count("ssytrf col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SSYTRF_ALLOC_TEST(1, 0, "ssytrf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SSYTRF_ALLOC_TEST(1, 1, "ssytrf transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SSYTRF_ALLOC_TEST(1, 2, "ssytrf allocation count", 0);
    lapacke_test_check_alloc_count("ssytrf row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SSYTRF_ALLOC_TEST(2, 0, "ssytrf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("ssytrf invalid layout allocation count");
}
