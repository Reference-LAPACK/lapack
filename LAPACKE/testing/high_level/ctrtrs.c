#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call ctrtrs. */
#define LAPACKE_CTRTRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_tri(layout, 'U', N, a, LD);                         \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_ctrtrs)(layout, 'U', 'N', 'N',   \
                                                      N, NRHS, a, LD, b, LD),  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(ctrtrs)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrtrs a uplo=U diag=N", l, N, N, a, LD, lapacke_test_region_upper,
            -7,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_ctrtrs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrtrs b uplo=U diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_ctrtrs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrtrs a uplo=U diag=U", l, N, N, a, LD,
            lapacke_test_region_strict_upper, -7,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_ctrtrs)(layout, 'U', 'N', 'U', N, NRHS, a, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrtrs b uplo=U diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_ctrtrs)(layout, 'U', 'N', 'U', N, NRHS, a, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrtrs a uplo=L diag=N", l, N, N, a, LD, lapacke_test_region_lower,
            -7,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_ctrtrs)(layout, 'L', 'N', 'N', N, NRHS, a, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrtrs b uplo=L diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_ctrtrs)(layout, 'L', 'N', 'N', N, NRHS, a, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrtrs a uplo=L diag=U", l, N, N, a, LD,
            lapacke_test_region_strict_lower, -7,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_ctrtrs)(layout, 'L', 'N', 'U', N, NRHS, a, LD, b,
                                       LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrtrs b uplo=L diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_ctrtrs)(layout, 'L', 'N', 'U', N, NRHS, a, LD, b,
                                       LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check("ctrtrs NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_ctrtrs)(layout, 'U', 'N', 'N', N,
                                                      NRHS, a, LD, b, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CTRTRS_ALLOC_TEST(0, 0, "ctrtrs allocation count", 0);
    lapacke_test_check_alloc_count("ctrtrs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CTRTRS_ALLOC_TEST(1, 0, "ctrtrs transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CTRTRS_ALLOC_TEST(1, 1, "ctrtrs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CTRTRS_ALLOC_TEST(1, 2, "ctrtrs allocation count", 0);
    lapacke_test_check_alloc_count("ctrtrs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CTRTRS_ALLOC_TEST(2, 0, "ctrtrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("ctrtrs invalid layout allocation count");
}
