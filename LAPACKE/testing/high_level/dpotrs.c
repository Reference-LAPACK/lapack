#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* One allocation test: refill the inputs, schedule the malloc failure
 * countdown and check the info result of a dpotrs call in the indexed
 * layout. */
#define LAPACKE_DPOTRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        lapacke_test_fill_spd(lapacke_test_layouts[layout_index], N, a, LD);   \
        lapacke_test_fill_rhs(lapacke_test_layouts[layout_index], N, NRHS, b,  \
                              LD);                                             \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           LAPACKE_dpotrs(lapacke_test_layouts[layout_index],  \
                                          'U', N, NRHS, a, LD, b, LD),         \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dpotrs)
{
    double a[LD * LD], b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        /* Only the uplo triangle of the factor is a documented input. */
        LAPACKE_TEST_NAN_SWEEP(
            "dpotrs a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -5,
            (lapacke_test_fill_spd(layout, N, a, LD),
             lapacke_test_fill_rhs(layout, N, NRHS, b, LD)),
            LAPACKE_dpotrs(layout, 'U', N, NRHS, a, LD, b, LD));

        LAPACKE_TEST_NAN_SWEEP(
            "dpotrs a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -5,
            (lapacke_test_fill_spd(layout, N, a, LD),
             lapacke_test_fill_rhs(layout, N, NRHS, b, LD)),
            LAPACKE_dpotrs(layout, 'L', N, NRHS, a, LD, b, LD));

        LAPACKE_TEST_NAN_SWEEP(
            "dpotrs b", l, N, NRHS, b, LD, lapacke_test_region_full, -7,
            (lapacke_test_fill_spd(layout, N, a, LD),
             lapacke_test_fill_rhs(layout, N, NRHS, b, LD)),
            LAPACKE_dpotrs(layout, 'U', N, NRHS, a, LD, b, LD));

        /* With NaN checking disabled even all-NaN input must go through to
         * the Fortran routine (valid arguments, so info must not be
         * negative). */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_nan(layout, N, N, a, LD);
        lapacke_test_fill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check(
            "dpotrs NaN with nancheck off", lapacke_test_layout_names[l],
            LAPACKE_dpotrs(layout, 'U', N, NRHS, a, LD, b, LD) >= 0 ? 0 : -999,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* Row-major allocates the transposed copies of A, then B. */
    LAPACKE_DPOTRS_ALLOC_TEST(1, 0, "dpotrs transpose alloc failure (a)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DPOTRS_ALLOC_TEST(1, 1, "dpotrs transpose alloc failure (b)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Recovery, scheduled one past the last allocation. */
    LAPACKE_DPOTRS_ALLOC_TEST(1, 2, "dpotrs recovers after alloc failure", 0);
    lapacke_test_check_alloc_count("dpotrs allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_DPOTRS_ALLOC_TEST(2, 0, "dpotrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dpotrs invalid layout allocation count");
}
