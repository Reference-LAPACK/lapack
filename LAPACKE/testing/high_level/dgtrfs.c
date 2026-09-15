#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dgtrfs. */
#define LAPACKE_DGTRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_vec(LD * LD, dl);                                   \
        lapacke_test_dfill_pos(LD * LD, d);                                    \
        lapacke_test_dfill_vec(LD * LD, du);                                   \
        lapacke_test_dfill_vec(LD * LD, dlf);                                  \
        lapacke_test_dfill_pos(LD * LD, df);                                   \
        lapacke_test_dfill_vec(LD * LD, duf);                                  \
        lapacke_test_dfill_vec(LD * LD, du2);                                  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_dfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_dfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dgtrfs)(                         \
                               layout, 'N', N, NRHS, dl, d, du, dlf, df, duf,  \
                               du2, ipiv, b, LD, x, LD, ferr, berr),           \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dgtrfs)
{
    double dl[LD * LD];
    double d[LD * LD];
    double du[LD * LD];
    double dlf[LD * LD];
    double df[LD * LD];
    double duf[LD * LD];
    double du2[LD * LD];
    lapack_int ipiv[LD * LD];
    double b[LD * LD];
    double x[LD * LD];
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgtrfs b", l, N, NRHS, b, LD, lapacke_test_region_full, -13,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, dlf),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_dfill_vec(LD * LD, duf),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgtrfs d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -6,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, dlf),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_dfill_vec(LD * LD, duf),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgtrfs df", l, N, 1, df, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -9,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, dlf),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_dfill_vec(LD * LD, duf),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgtrfs dl", l, N - 1, 1, dl, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -5,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, dlf),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_dfill_vec(LD * LD, duf),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgtrfs dlf", l, N - 1, 1, dlf, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -8,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, dlf),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_dfill_vec(LD * LD, duf),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgtrfs du", l, N - 1, 1, du, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -7,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, dlf),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_dfill_vec(LD * LD, duf),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgtrfs du2", l, N - 2, 1, du2, LAPACKE_TEST_VLD(layout, N - 2),
            lapacke_test_region_full, -11,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, dlf),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_dfill_vec(LD * LD, duf),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgtrfs duf", l, N - 1, 1, duf, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -10,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, dlf),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_dfill_vec(LD * LD, duf),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgtrfs x", l, N, NRHS, x, LD, lapacke_test_region_full, -15,
            (lapacke_test_dfill_vec(LD * LD, dl),
             lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_dfill_vec(LD * LD, du),
             lapacke_test_dfill_vec(LD * LD, dlf),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_dfill_vec(LD * LD, duf),
             lapacke_test_dfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_dfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_dfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_dfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_dfill_nan(layout, N, 1, df, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_dfill_nan(layout, N - 1, 1, dl,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_dfill_nan(layout, N - 1, 1, dlf,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_dfill_nan(layout, N - 1, 1, du,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_dfill_nan(layout, N - 2, 1, du2,
                               LAPACKE_TEST_VLD(layout, N - 2));
        lapacke_test_dfill_nan(layout, N - 1, 1, duf,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_dfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "dgtrfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DGTRFS_ALLOC_TEST(0, 0, "dgtrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGTRFS_ALLOC_TEST(0, 1, "dgtrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGTRFS_ALLOC_TEST(0, 2, "dgtrfs allocation count", 0);
    lapacke_test_check_alloc_count("dgtrfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGTRFS_ALLOC_TEST(1, 0, "dgtrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGTRFS_ALLOC_TEST(1, 1, "dgtrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGTRFS_ALLOC_TEST(1, 2, "dgtrfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGTRFS_ALLOC_TEST(1, 3, "dgtrfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGTRFS_ALLOC_TEST(1, 4, "dgtrfs allocation count", 0);
    lapacke_test_check_alloc_count("dgtrfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGTRFS_ALLOC_TEST(2, 0, "dgtrfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgtrfs invalid layout allocation count");
}
