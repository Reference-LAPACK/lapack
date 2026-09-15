#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define K N

/* Refill the inputs, schedule the malloc failure, call sorgql. */
#define LAPACKE_SORGQL_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, M, N, a, LD);                               \
        lapacke_test_sfill_vec(LD * LD, tau);                                  \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_sorgql)(layout, M, N, K, a, LD, tau),           \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(sorgql)
{
    float a[LD * LD];
    float tau[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sorgql a", l, M, N, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_sfill_vec(LD * LD, tau)),
            API_SUFFIX(LAPACKE_sorgql)(layout, M, N, K, a, LD, tau));

        LAPACKE_TEST_SNAN_SWEEP(
            "sorgql tau", l, K, 1, tau, LAPACKE_TEST_VLD(layout, K),
            lapacke_test_region_full, -7,
            (lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_sfill_vec(LD * LD, tau)),
            API_SUFFIX(LAPACKE_sorgql)(layout, M, N, K, a, LD, tau));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, M, N, a, LD);
        lapacke_test_sfill_nan(layout, K, 1, tau, LAPACKE_TEST_VLD(layout, K));
        lapacke_test_check(
            "sorgql NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sorgql)(layout, M, N, K, a, LD, tau) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SORGQL_ALLOC_TEST(0, 0, "sorgql work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SORGQL_ALLOC_TEST(0, 1, "sorgql allocation count", 0);
    lapacke_test_check_alloc_count("sorgql col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SORGQL_ALLOC_TEST(1, 0, "sorgql work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SORGQL_ALLOC_TEST(1, 1, "sorgql transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SORGQL_ALLOC_TEST(1, 2, "sorgql allocation count", 0);
    lapacke_test_check_alloc_count("sorgql row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SORGQL_ALLOC_TEST(2, 0, "sorgql invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sorgql invalid layout allocation count");
}
