#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cppsv. */
#define LAPACKE_CPPSV_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_pp(layout, 'U', N, ap);                             \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_cppsv)(layout, 'U', N, NRHS, ap, b, LD),        \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(cppsv)
{
    lapack_complex_float ap[LD * LD];
    lapack_complex_float b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cppsv ap uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -5,
            (lapacke_test_cfill_pp(layout, 'U', N, ap),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cppsv)(layout, 'U', N, NRHS, ap, b, LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cppsv b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -6,
            (lapacke_test_cfill_pp(layout, 'U', N, ap),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cppsv)(layout, 'U', N, NRHS, ap, b, LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cppsv ap uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -5,
            (lapacke_test_cfill_pp(layout, 'L', N, ap),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cppsv)(layout, 'L', N, NRHS, ap, b, LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cppsv b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -6,
            (lapacke_test_cfill_pp(layout, 'L', N, ap),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_cppsv)(layout, 'L', N, NRHS, ap, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check(
            "cppsv NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cppsv)(layout, 'U', N, NRHS, ap, b, LD) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CPPSV_ALLOC_TEST(0, 0, "cppsv allocation count", 0);
    lapacke_test_check_alloc_count("cppsv col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CPPSV_ALLOC_TEST(1, 0, "cppsv transpose alloc failure (b_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPPSV_ALLOC_TEST(1, 1, "cppsv transpose alloc failure (ap_t)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPPSV_ALLOC_TEST(1, 2, "cppsv allocation count", 0);
    lapacke_test_check_alloc_count("cppsv row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CPPSV_ALLOC_TEST(2, 0, "cppsv invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cppsv invalid layout allocation count");
}
