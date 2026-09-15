#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sgerfs. */
#define LAPACKE_SGERFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, N, N, a, LD);                               \
        lapacke_test_sfill(layout, N, N, af, LD);                              \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_sfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sgerfs)(layout, 'N', N, NRHS, a, \
                                                      LD, af, LD, ipiv, b, LD, \
                                                      x, LD, ferr, berr),      \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sgerfs)
{
    float a[LD * LD];
    float af[LD * LD];
    lapack_int ipiv[LD * LD];
    float b[LD * LD];
    float x[LD * LD];
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgerfs a", l, N, N, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_sfill(layout, N, N, a, LD),
             lapacke_test_sfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sgerfs)(layout, 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgerfs af", l, N, N, af, LD, lapacke_test_region_full, -7,
            (lapacke_test_sfill(layout, N, N, a, LD),
             lapacke_test_sfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sgerfs)(layout, 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgerfs b", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_sfill(layout, N, N, a, LD),
             lapacke_test_sfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sgerfs)(layout, 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgerfs x", l, N, NRHS, x, LD, lapacke_test_region_full, -12,
            (lapacke_test_sfill(layout, N, N, a, LD),
             lapacke_test_sfill(layout, N, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sgerfs)(layout, 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_sfill_nan(layout, N, N, a, LD);
        lapacke_test_sfill_nan(layout, N, N, af, LD);
        lapacke_test_sfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_sfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "sgerfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sgerfs)(layout, 'N', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SGERFS_ALLOC_TEST(0, 0, "sgerfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGERFS_ALLOC_TEST(0, 1, "sgerfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGERFS_ALLOC_TEST(0, 2, "sgerfs allocation count", 0);
    lapacke_test_check_alloc_count("sgerfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SGERFS_ALLOC_TEST(1, 0, "sgerfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGERFS_ALLOC_TEST(1, 1, "sgerfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGERFS_ALLOC_TEST(1, 2, "sgerfs transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGERFS_ALLOC_TEST(1, 3, "sgerfs transpose alloc failure (af_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGERFS_ALLOC_TEST(1, 4, "sgerfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGERFS_ALLOC_TEST(1, 5, "sgerfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGERFS_ALLOC_TEST(1, 6, "sgerfs allocation count", 0);
    lapacke_test_check_alloc_count("sgerfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SGERFS_ALLOC_TEST(2, 0, "sgerfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgerfs invalid layout allocation count");
}
