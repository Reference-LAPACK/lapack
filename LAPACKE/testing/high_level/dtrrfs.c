#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dtrrfs. */
#define LAPACKE_DTRRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_tri(layout, 'U', N, a, LD);                         \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_dfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dtrrfs)(layout, 'U', 'N', 'N',   \
                                                      N, NRHS, a, LD, b, LD,   \
                                                      x, LD, ferr, berr),      \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dtrrfs)
{
    double a[LD * LD];
    double b[LD * LD];
    double x[LD * LD];
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs a uplo=U diag=N", l, N, N, a, LD, lapacke_test_region_upper,
            -7,
            (lapacke_test_dfill_tri(layout, 'U', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs b uplo=U diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_dfill_tri(layout, 'U', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs x uplo=U diag=N", l, N, NRHS, x, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_dfill_tri(layout, 'U', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs a uplo=U diag=U", l, N, N, a, LD,
            lapacke_test_region_strict_upper, -7,
            (lapacke_test_dfill_tri(layout, 'U', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'U', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs b uplo=U diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_dfill_tri(layout, 'U', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'U', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs x uplo=U diag=U", l, N, NRHS, x, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_dfill_tri(layout, 'U', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'U', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs a uplo=L diag=N", l, N, N, a, LD, lapacke_test_region_lower,
            -7,
            (lapacke_test_dfill_tri(layout, 'L', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'L', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs b uplo=L diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_dfill_tri(layout, 'L', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'L', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs x uplo=L diag=N", l, N, NRHS, x, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_dfill_tri(layout, 'L', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'L', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs a uplo=L diag=U", l, N, N, a, LD,
            lapacke_test_region_strict_lower, -7,
            (lapacke_test_dfill_tri(layout, 'L', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'L', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs b uplo=L diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -9,
            (lapacke_test_dfill_tri(layout, 'L', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'L', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dtrrfs x uplo=L diag=U", l, N, NRHS, x, LD,
            lapacke_test_region_full, -11,
            (lapacke_test_dfill_tri(layout, 'L', N, a, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'L', 'N', 'U', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_dfill_nan(layout, N, N, a, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "dtrrfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dtrrfs)(layout, 'U', 'N', 'N', N, NRHS, a, LD, b,
                                       LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DTRRFS_ALLOC_TEST(0, 0, "dtrrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DTRRFS_ALLOC_TEST(0, 1, "dtrrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DTRRFS_ALLOC_TEST(0, 2, "dtrrfs allocation count", 0);
    lapacke_test_check_alloc_count("dtrrfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DTRRFS_ALLOC_TEST(1, 0, "dtrrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DTRRFS_ALLOC_TEST(1, 1, "dtrrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DTRRFS_ALLOC_TEST(1, 2, "dtrrfs transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DTRRFS_ALLOC_TEST(1, 3, "dtrrfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DTRRFS_ALLOC_TEST(1, 4, "dtrrfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DTRRFS_ALLOC_TEST(1, 5, "dtrrfs allocation count", 0);
    lapacke_test_check_alloc_count("dtrrfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DTRRFS_ALLOC_TEST(2, 0, "dtrrfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dtrrfs invalid layout allocation count");
}
