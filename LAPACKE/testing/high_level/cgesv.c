#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cgesv. */
#define LAPACKE_CGESV_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, N, N, a, LD);                               \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_cgesv)(layout, N, NRHS, a, LD, ipiv, b, LD),    \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(cgesv)
{
    lapack_complex_float a[LD * LD];
    lapack_int ipiv[LD * LD];
    lapack_complex_float b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cgesv a", l, N, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_cfill(layout, N, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cgesv)(layout, N, NRHS, a, LD, ipiv, b, LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgesv b", l, N, NRHS, b, LD, lapacke_test_region_full, -7,
            (lapacke_test_cfill(layout, N, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cgesv)(layout, N, NRHS, a, LD, ipiv, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check(
            "cgesv NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgesv)(layout, N, NRHS, a, LD, ipiv, b, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CGESV_ALLOC_TEST(0, 0, "cgesv allocation count", 0);
    lapacke_test_check_alloc_count("cgesv col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CGESV_ALLOC_TEST(1, 0, "cgesv transpose alloc failure (a_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGESV_ALLOC_TEST(1, 1, "cgesv transpose alloc failure (b_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGESV_ALLOC_TEST(1, 2, "cgesv allocation count", 0);
    lapacke_test_check_alloc_count("cgesv row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CGESV_ALLOC_TEST(2, 0, "cgesv invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cgesv invalid layout allocation count");
}
