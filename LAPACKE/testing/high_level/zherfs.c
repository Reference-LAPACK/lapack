#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zherfs. */
#define LAPACKE_ZHERFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_sym(layout, N, a, LD);                              \
        lapacke_test_zfill_sym(layout, N, af, LD);                             \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_zfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zherfs)(layout, 'U', N, NRHS, a, \
                                                      LD, af, LD, ipiv, b, LD, \
                                                      x, LD, ferr, berr),      \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zherfs)
{
    lapack_complex_double a[LD * LD];
    lapack_complex_double af[LD * LD];
    lapack_int ipiv[LD * LD];
    lapack_complex_double b[LD * LD];
    lapack_complex_double x[LD * LD];
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zherfs a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -5,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zherfs)(layout, 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zherfs af uplo=U", l, N, N, af, LD, lapacke_test_region_upper, -7,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zherfs)(layout, 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zherfs b uplo=U", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zherfs)(layout, 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zherfs x uplo=U", l, N, NRHS, x, LD, lapacke_test_region_full, -12,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zherfs)(layout, 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zherfs a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -5,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zherfs)(layout, 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zherfs af uplo=L", l, N, N, af, LD, lapacke_test_region_lower, -7,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zherfs)(layout, 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zherfs b uplo=L", l, N, NRHS, b, LD, lapacke_test_region_full, -10,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zherfs)(layout, 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zherfs x uplo=L", l, N, NRHS, x, LD, lapacke_test_region_full, -12,
            (lapacke_test_zfill_sym(layout, N, a, LD),
             lapacke_test_zfill_sym(layout, N, af, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zherfs)(layout, 'L', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_zfill_nan(layout, N, N, a, LD);
        lapacke_test_zfill_nan(layout, N, N, af, LD);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_zfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "zherfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zherfs)(layout, 'U', N, NRHS, a, LD, af, LD,
                                       ipiv, b, LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZHERFS_ALLOC_TEST(0, 0, "zherfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZHERFS_ALLOC_TEST(0, 1, "zherfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZHERFS_ALLOC_TEST(0, 2, "zherfs allocation count", 0);
    lapacke_test_check_alloc_count("zherfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZHERFS_ALLOC_TEST(1, 0, "zherfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZHERFS_ALLOC_TEST(1, 1, "zherfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZHERFS_ALLOC_TEST(1, 2, "zherfs transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZHERFS_ALLOC_TEST(1, 3, "zherfs transpose alloc failure (af_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZHERFS_ALLOC_TEST(1, 4, "zherfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZHERFS_ALLOC_TEST(1, 5, "zherfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZHERFS_ALLOC_TEST(1, 6, "zherfs allocation count", 0);
    lapacke_test_check_alloc_count("zherfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZHERFS_ALLOC_TEST(2, 0, "zherfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zherfs invalid layout allocation count");
}
