#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call claset. */
#define LAPACKE_CLASET_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        alpha[0] = lapack_make_complex_float(1.0f, 0.0f);                      \
        beta[0] = lapack_make_complex_float(2.0f, 0.0f);                       \
        lapacke_test_cfill(layout, M, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_claset)(                         \
                               layout, 'U', M, N, alpha[0], beta[0], a, LD),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(claset)
{
    lapack_complex_float alpha[1];
    lapack_complex_float beta[1];
    lapack_complex_float a[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "claset alpha", l, 1, 1, alpha, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -5,
            ((alpha[0] = lapack_make_complex_float(1.0f, 0.0f)),
             (beta[0] = lapack_make_complex_float(2.0f, 0.0f)),
             lapacke_test_cfill(layout, M, N, a, LD)),
            API_SUFFIX(LAPACKE_claset)(layout, 'U', M, N, alpha[0], beta[0], a,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "claset beta", l, 1, 1, beta, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -6,
            ((alpha[0] = lapack_make_complex_float(1.0f, 0.0f)),
             (beta[0] = lapack_make_complex_float(2.0f, 0.0f)),
             lapacke_test_cfill(layout, M, N, a, LD)),
            API_SUFFIX(LAPACKE_claset)(layout, 'U', M, N, alpha[0], beta[0], a,
                                       LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill(layout, M, N, a, LD);
        lapacke_test_cfill_nan(layout, 1, 1, alpha,
                               LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_cfill_nan(layout, 1, 1, beta, LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_check("claset NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_claset)(
                               layout, 'U', M, N, alpha[0], beta[0], a, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CLASET_ALLOC_TEST(0, 0, "claset allocation count", 0);
    lapacke_test_check_alloc_count("claset col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CLASET_ALLOC_TEST(1, 0, "claset transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CLASET_ALLOC_TEST(1, 1, "claset allocation count", 0);
    lapacke_test_check_alloc_count("claset row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CLASET_ALLOC_TEST(2, 0, "claset invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("claset invalid layout allocation count");
}
