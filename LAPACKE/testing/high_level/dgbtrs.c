#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

#define REGION_AB(i, j) lapacke_test_region_band(i, j, N, KL + KU)

/* Refill the inputs, schedule the malloc failure, call dgbtrs. */
#define LAPACKE_DGBTRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab, LD);  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dgbtrs)(layout, 'N', N, KL, KU,  \
                                                      NRHS, ab, LD, ipiv, b,   \
                                                      LD),                     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dgbtrs)
{
    double ab[LD * LD];
    lapack_int ipiv[LD * LD];
    double b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgbtrs ab", l, 2 * KL + KU + 1, N, ab, LD, REGION_AB, -7,
            (lapacke_test_dfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgbtrs)(layout, 'N', N, KL, KU, NRHS, ab, LD,
                                       ipiv, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgbtrs b", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_dfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgbtrs)(layout, 'N', N, KL, KU, NRHS, ab, LD,
                                       ipiv, b, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_dfill_nan(layout, 2 * KL + KU + 1, N, ab, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check(
            "dgbtrs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgbtrs)(layout, 'N', N, KL, KU, NRHS, ab, LD,
                                       ipiv, b, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_DGBTRS_ALLOC_TEST(0, 0, "dgbtrs allocation count", 0);
    lapacke_test_check_alloc_count("dgbtrs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGBTRS_ALLOC_TEST(1, 0, "dgbtrs transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGBTRS_ALLOC_TEST(1, 1, "dgbtrs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGBTRS_ALLOC_TEST(1, 2, "dgbtrs allocation count", 0);
    lapacke_test_check_alloc_count("dgbtrs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGBTRS_ALLOC_TEST(2, 0, "dgbtrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgbtrs invalid layout allocation count");
}
