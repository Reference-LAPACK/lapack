#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dlaset. */
#define LAPACKE_DLASET_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        alpha[0] = 1.0;                                                        \
        beta[0] = 2.0;                                                         \
        lapacke_test_dfill(layout, M, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dlaset)(                         \
                               layout, 'U', M, N, alpha[0], beta[0], a, LD),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dlaset)
{
    double alpha[1];
    double beta[1];
    double a[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dlaset alpha", l, 1, 1, alpha, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -5,
            ((alpha[0] = 1.0), (beta[0] = 2.0),
             lapacke_test_dfill(layout, M, N, a, LD)),
            API_SUFFIX(LAPACKE_dlaset)(layout, 'U', M, N, alpha[0], beta[0], a,
                                       LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dlaset beta", l, 1, 1, beta, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            ((alpha[0] = 1.0), (beta[0] = 2.0),
             lapacke_test_dfill(layout, M, N, a, LD)),
            API_SUFFIX(LAPACKE_dlaset)(layout, 'U', M, N, alpha[0], beta[0], a,
                                       LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill(layout, M, N, a, LD);
        lapacke_test_dfill_nan(layout, 1, 1, alpha,
                               LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_dfill_nan(layout, 1, 1, beta, LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_check("dlaset NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_dlaset)(
                               layout, 'U', M, N, alpha[0], beta[0], a, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_DLASET_ALLOC_TEST(0, 0, "dlaset allocation count", 0);
    lapacke_test_check_alloc_count("dlaset col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DLASET_ALLOC_TEST(1, 0, "dlaset transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DLASET_ALLOC_TEST(1, 1, "dlaset allocation count", 0);
    lapacke_test_check_alloc_count("dlaset row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DLASET_ALLOC_TEST(2, 0, "dlaset invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dlaset invalid layout allocation count");
}
