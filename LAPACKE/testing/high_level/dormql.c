#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define K N

/* Refill the inputs, schedule the malloc failure, call dormql. */
#define LAPACKE_DORMQL_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, M, K, a, LD);                               \
        lapacke_test_dfill_vec(LD * LD, tau);                                  \
        lapacke_test_dfill_rhs(layout, M, N, c, LD);                           \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dormql)(layout, 'L', 'N', M, N,  \
                                                      K, a, LD, tau, c, LD),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dormql)
{
    double a[LD * LD];
    double tau[LD * LD];
    double c[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dormql a side=L", l, M, K, a, LD, lapacke_test_region_full, -7,
            (lapacke_test_dfill(layout, M, K, a, LD),
             lapacke_test_dfill_vec(LD * LD, tau),
             lapacke_test_dfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_dormql)(layout, 'L', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dormql c side=L", l, M, N, c, LD, lapacke_test_region_full, -10,
            (lapacke_test_dfill(layout, M, K, a, LD),
             lapacke_test_dfill_vec(LD * LD, tau),
             lapacke_test_dfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_dormql)(layout, 'L', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dormql tau side=L", l, K, 1, tau, LAPACKE_TEST_VLD(layout, K),
            lapacke_test_region_full, -9,
            (lapacke_test_dfill(layout, M, K, a, LD),
             lapacke_test_dfill_vec(LD * LD, tau),
             lapacke_test_dfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_dormql)(layout, 'L', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dormql a side=R", l, N, K, a, LD, lapacke_test_region_full, -7,
            (lapacke_test_dfill(layout, N, K, a, LD),
             lapacke_test_dfill_vec(LD * LD, tau),
             lapacke_test_dfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_dormql)(layout, 'R', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dormql c side=R", l, M, N, c, LD, lapacke_test_region_full, -10,
            (lapacke_test_dfill(layout, N, K, a, LD),
             lapacke_test_dfill_vec(LD * LD, tau),
             lapacke_test_dfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_dormql)(layout, 'R', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dormql tau side=R", l, K, 1, tau, LAPACKE_TEST_VLD(layout, K),
            lapacke_test_region_full, -9,
            (lapacke_test_dfill(layout, N, K, a, LD),
             lapacke_test_dfill_vec(LD * LD, tau),
             lapacke_test_dfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_dormql)(layout, 'R', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_nan(layout, M, K, a, LD);
        lapacke_test_dfill_nan(layout, M, N, c, LD);
        lapacke_test_dfill_nan(layout, K, 1, tau, LAPACKE_TEST_VLD(layout, K));
        lapacke_test_check("dormql NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_dormql)(layout, 'L', 'N', M, N, K,
                                                      a, LD, tau, c, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DORMQL_ALLOC_TEST(0, 0, "dormql work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DORMQL_ALLOC_TEST(0, 1, "dormql allocation count", 0);
    lapacke_test_check_alloc_count("dormql col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DORMQL_ALLOC_TEST(1, 0, "dormql work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DORMQL_ALLOC_TEST(1, 1, "dormql transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DORMQL_ALLOC_TEST(1, 2, "dormql transpose alloc failure (c_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DORMQL_ALLOC_TEST(1, 3, "dormql allocation count", 0);
    lapacke_test_check_alloc_count("dormql row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DORMQL_ALLOC_TEST(2, 0, "dormql invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dormql invalid layout allocation count");
}
