#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call ctrrfs. */
#define LAPACKE_CTRRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_tri(layout, 'U', N, a, LD);                         \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_cfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_ctrrfs)(layout, 'U', 'N', 'N',   \
                                                      N, NRHS, a, LD, b, LD,   \
                                                      x, LD, ferr, berr),      \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(ctrrfs)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float b[LD * LD];
    lapack_complex_float x[LD * LD];
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs a uplo=U diag=N", l, N, N, a, LD, lapacke_test_region_upper,
            -7,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs b uplo=U diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs x uplo=U diag=N", l, N, NRHS, x, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs a uplo=U diag=U", l, N, N, a, LD,
            lapacke_test_region_strict_upper, -7,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'U', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs b uplo=U diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'U', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs x uplo=U diag=U", l, N, NRHS, x, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_cfill_tri(layout, 'U', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'U', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs a uplo=L diag=N", l, N, N, a, LD, lapacke_test_region_lower,
            -7,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'L', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs b uplo=L diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'L', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs x uplo=L diag=N", l, N, NRHS, x, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'L', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs a uplo=L diag=U", l, N, N, a, LD,
            lapacke_test_region_strict_lower, -7,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'L', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs b uplo=L diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'L', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "ctrrfs x uplo=L diag=U", l, N, NRHS, x, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_cfill_tri(layout, 'L', N, a, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'L', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_cfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "ctrrfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_ctrrfs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CTRRFS_ALLOC_TEST(0, 0, "ctrrfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTRRFS_ALLOC_TEST(0, 1, "ctrrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTRRFS_ALLOC_TEST(0, 2, "ctrrfs allocation count", 0);
    lapacke_test_check_alloc_count("ctrrfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CTRRFS_ALLOC_TEST(1, 0, "ctrrfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTRRFS_ALLOC_TEST(1, 1, "ctrrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CTRRFS_ALLOC_TEST(1, 2, "ctrrfs transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CTRRFS_ALLOC_TEST(1, 3, "ctrrfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CTRRFS_ALLOC_TEST(1, 4, "ctrrfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CTRRFS_ALLOC_TEST(1, 5, "ctrrfs allocation count", 0);
    lapacke_test_check_alloc_count("ctrrfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CTRRFS_ALLOC_TEST(2, 0, "ctrrfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("ctrrfs invalid layout allocation count");
}
