#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dpptrs. */
#define LAPACKE_DPPTRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_pp(layout, 'U', N, ap);                             \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dpptrs)(layout, 'U', N, NRHS, ap, b, LD),       \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(dpptrs)
{
    double ap[LD * LD];
    double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dpptrs ap uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -5,
            (lapacke_test_dfill_pp(layout, 'U', N, ap),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dpptrs)(layout, 'U', N, NRHS, ap, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpptrs b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -6,
            (lapacke_test_dfill_pp(layout, 'U', N, ap),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dpptrs)(layout, 'U', N, NRHS, ap, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpptrs ap uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -5,
            (lapacke_test_dfill_pp(layout, 'L', N, ap),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dpptrs)(layout, 'L', N, NRHS, ap, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dpptrs b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -6,
            (lapacke_test_dfill_pp(layout, 'L', N, ap),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dpptrs)(layout, 'L', N, NRHS, ap, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check(
            "dpptrs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dpptrs)(layout, 'U', N, NRHS, ap, b, LD) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_DPPTRS_ALLOC_TEST(0, 0, "dpptrs allocation count", 0);
    lapacke_test_check_alloc_count("dpptrs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DPPTRS_ALLOC_TEST(1, 0, "dpptrs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPPTRS_ALLOC_TEST(1, 1, "dpptrs transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPPTRS_ALLOC_TEST(1, 2, "dpptrs allocation count", 0);
    lapacke_test_check_alloc_count("dpptrs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DPPTRS_ALLOC_TEST(2, 0, "dpptrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dpptrs invalid layout allocation count");
}
