#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define K N

/* Refill the inputs, schedule the malloc failure, call cunmlq. */
#define LAPACKE_CUNMLQ_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, K, M, a, LD);                               \
        lapacke_test_cfill_vec(LD * LD, tau);                                  \
        lapacke_test_cfill_rhs(layout, M, N, c, LD);                           \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cunmlq)(layout, 'L', 'N', M, N,  \
                                                      K, a, LD, tau, c, LD),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cunmlq)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float tau[LD * LD];
    lapack_complex_float c[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cunmlq a side=L", l, K, M, a, LD, lapacke_test_region_full, -7,
            (lapacke_test_cfill(layout, K, M, a, LD),
             lapacke_test_cfill_vec(LD * LD, tau),
             lapacke_test_cfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_cunmlq)(layout, 'L', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cunmlq c side=L", l, M, N, c, LD, lapacke_test_region_full, -10,
            (lapacke_test_cfill(layout, K, M, a, LD),
             lapacke_test_cfill_vec(LD * LD, tau),
             lapacke_test_cfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_cunmlq)(layout, 'L', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cunmlq tau side=L", l, K, 1, tau, LAPACKE_TEST_VLD(layout, K),
            lapacke_test_region_full, -9,
            (lapacke_test_cfill(layout, K, M, a, LD),
             lapacke_test_cfill_vec(LD * LD, tau),
             lapacke_test_cfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_cunmlq)(layout, 'L', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cunmlq a side=R", l, K, N, a, LD, lapacke_test_region_full, -7,
            (lapacke_test_cfill(layout, K, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, tau),
             lapacke_test_cfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_cunmlq)(layout, 'R', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cunmlq c side=R", l, M, N, c, LD, lapacke_test_region_full, -10,
            (lapacke_test_cfill(layout, K, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, tau),
             lapacke_test_cfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_cunmlq)(layout, 'R', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cunmlq tau side=R", l, K, 1, tau, LAPACKE_TEST_VLD(layout, K),
            lapacke_test_region_full, -9,
            (lapacke_test_cfill(layout, K, N, a, LD),
             lapacke_test_cfill_vec(LD * LD, tau),
             lapacke_test_cfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_cunmlq)(layout, 'R', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, K, M, a, LD);
        lapacke_test_cfill_nan(layout, M, N, c, LD);
        lapacke_test_cfill_nan(layout, K, 1, tau, LAPACKE_TEST_VLD(layout, K));
        lapacke_test_check("cunmlq NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_cunmlq)(layout, 'L', 'N', M, N, K,
                                                      a, LD, tau, c, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CUNMLQ_ALLOC_TEST(0, 0, "cunmlq work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CUNMLQ_ALLOC_TEST(0, 1, "cunmlq allocation count", 0);
    lapacke_test_check_alloc_count("cunmlq col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CUNMLQ_ALLOC_TEST(1, 0, "cunmlq work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CUNMLQ_ALLOC_TEST(1, 1, "cunmlq transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CUNMLQ_ALLOC_TEST(1, 2, "cunmlq transpose alloc failure (c_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CUNMLQ_ALLOC_TEST(1, 3, "cunmlq allocation count", 0);
    lapacke_test_check_alloc_count("cunmlq row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CUNMLQ_ALLOC_TEST(2, 0, "cunmlq invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cunmlq invalid layout allocation count");
}
