#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cgtrfs. */
#define LAPACKE_CGTRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_vec(LD * LD, dl);                                   \
        lapacke_test_cfill_pos(LD * LD, d);                                    \
        lapacke_test_cfill_vec(LD * LD, du);                                   \
        lapacke_test_cfill_vec(LD * LD, dlf);                                  \
        lapacke_test_cfill_pos(LD * LD, df);                                   \
        lapacke_test_cfill_vec(LD * LD, duf);                                  \
        lapacke_test_cfill_vec(LD * LD, du2);                                  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_cfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_cfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cgtrfs)(                         \
                               layout, 'N', N, NRHS, dl, d, du, dlf, df, duf,  \
                               du2, ipiv, b, LD, x, LD, ferr, berr),           \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cgtrfs)
{
    lapack_complex_float dl[LD * LD];
    lapack_complex_float d[LD * LD];
    lapack_complex_float du[LD * LD];
    lapack_complex_float dlf[LD * LD];
    lapack_complex_float df[LD * LD];
    lapack_complex_float duf[LD * LD];
    lapack_complex_float du2[LD * LD];
    lapack_int ipiv[LD * LD];
    lapack_complex_float b[LD * LD];
    lapack_complex_float x[LD * LD];
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtrfs b", l, N, NRHS, b, LD, lapacke_test_region_full, -13,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, dlf),
             lapacke_test_cfill_pos(LD * LD, df),
             lapacke_test_cfill_vec(LD * LD, duf),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtrfs d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -6,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, dlf),
             lapacke_test_cfill_pos(LD * LD, df),
             lapacke_test_cfill_vec(LD * LD, duf),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtrfs df", l, N, 1, df, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, dlf),
             lapacke_test_cfill_pos(LD * LD, df),
             lapacke_test_cfill_vec(LD * LD, duf),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtrfs dl", l, N - 1, 1, dl, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -5,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, dlf),
             lapacke_test_cfill_pos(LD * LD, df),
             lapacke_test_cfill_vec(LD * LD, duf),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtrfs dlf", l, N - 1, 1, dlf, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -8,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, dlf),
             lapacke_test_cfill_pos(LD * LD, df),
             lapacke_test_cfill_vec(LD * LD, duf),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtrfs du", l, N - 1, 1, du, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -7,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, dlf),
             lapacke_test_cfill_pos(LD * LD, df),
             lapacke_test_cfill_vec(LD * LD, duf),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtrfs du2", l, N - 2, 1, du2, LAPACKE_TEST_VLD(layout, N - 2),
            lapacke_test_region_full, -11,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, dlf),
             lapacke_test_cfill_pos(LD * LD, df),
             lapacke_test_cfill_vec(LD * LD, duf),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtrfs duf", l, N - 1, 1, duf, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -10,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, dlf),
             lapacke_test_cfill_pos(LD * LD, df),
             lapacke_test_cfill_vec(LD * LD, duf),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        LAPACKE_TEST_CNAN_SWEEP(
            "cgtrfs x", l, N, NRHS, x, LD, lapacke_test_region_full, -15,
            (lapacke_test_cfill_vec(LD * LD, dl),
             lapacke_test_cfill_pos(LD * LD, d),
             lapacke_test_cfill_vec(LD * LD, du),
             lapacke_test_cfill_vec(LD * LD, dlf),
             lapacke_test_cfill_pos(LD * LD, df),
             lapacke_test_cfill_vec(LD * LD, duf),
             lapacke_test_cfill_vec(LD * LD, du2),
             lapacke_test_fill_ipiv(LD * LD, ipiv),
             lapacke_test_cfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_cfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_cfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_cfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_cfill_nan(layout, N, 1, df, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_cfill_nan(layout, N - 1, 1, dl,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_cfill_nan(layout, N - 1, 1, dlf,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_cfill_nan(layout, N - 1, 1, du,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_cfill_nan(layout, N - 2, 1, du2,
                               LAPACKE_TEST_VLD(layout, N - 2));
        lapacke_test_cfill_nan(layout, N - 1, 1, duf,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_cfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "cgtrfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgtrfs)(layout, 'N', N, NRHS, dl, d, du, dlf, df,
                                       duf, du2, ipiv, b, LD, x, LD, ferr,
                                       berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CGTRFS_ALLOC_TEST(0, 0, "cgtrfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGTRFS_ALLOC_TEST(0, 1, "cgtrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGTRFS_ALLOC_TEST(0, 2, "cgtrfs allocation count", 0);
    lapacke_test_check_alloc_count("cgtrfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CGTRFS_ALLOC_TEST(1, 0, "cgtrfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGTRFS_ALLOC_TEST(1, 1, "cgtrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGTRFS_ALLOC_TEST(1, 2, "cgtrfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGTRFS_ALLOC_TEST(1, 3, "cgtrfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGTRFS_ALLOC_TEST(1, 4, "cgtrfs allocation count", 0);
    lapacke_test_check_alloc_count("cgtrfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CGTRFS_ALLOC_TEST(2, 0, "cgtrfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cgtrfs invalid layout allocation count");
}
