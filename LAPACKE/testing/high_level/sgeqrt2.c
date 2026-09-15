#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define NB 2

/* Refill the inputs, schedule the malloc failure, call sgeqrt2. */
#define LAPACKE_SGEQRT2_ALLOC_TEST(layout_index, countdown, name, expected)    \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, M, N, a, LD);                               \
        lapacke_test_sfill_rhs(layout, NB, N, t, LD);                          \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_sgeqrt2)(layout, M, N, a, LD, t, LD),           \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(sgeqrt2)
{
    float a[LD * LD];
    float t[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgeqrt2 a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_sfill_rhs(layout, NB, N, t, LD)),
            API_SUFFIX(LAPACKE_sgeqrt2)(layout, M, N, a, LD, t, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_rhs(layout, NB, N, t, LD);
        lapacke_test_sfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "sgeqrt2 NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sgeqrt2)(layout, M, N, a, LD, t, LD) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_SGEQRT2_ALLOC_TEST(0, 0, "sgeqrt2 allocation count", 0);
    lapacke_test_check_alloc_count("sgeqrt2 col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SGEQRT2_ALLOC_TEST(1, 0, "sgeqrt2 transpose alloc failure (a_t)",
                               LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGEQRT2_ALLOC_TEST(1, 1, "sgeqrt2 transpose alloc failure (t_t)",
                               LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGEQRT2_ALLOC_TEST(1, 2, "sgeqrt2 allocation count", 0);
    lapacke_test_check_alloc_count("sgeqrt2 row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SGEQRT2_ALLOC_TEST(2, 0, "sgeqrt2 invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgeqrt2 invalid layout allocation count");
}
