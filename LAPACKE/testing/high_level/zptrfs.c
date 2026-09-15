#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zptrfs. */
#define LAPACKE_ZPTRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_pos(LD * LD, d);                                    \
        lapacke_test_zfill_vec(LD * LD, e);                                    \
        lapacke_test_dfill_pos(LD * LD, df);                                   \
        lapacke_test_zfill_vec(LD * LD, ef);                                   \
        lapacke_test_zfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_zfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_dfill_vec(LD * LD, ferr);                                 \
        lapacke_test_dfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zptrfs)(layout, 'U', N, NRHS, d, \
                                                      e, df, ef, b, LD, x, LD, \
                                                      ferr, berr),             \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zptrfs)
{
    double d[LD * LD];
    lapack_complex_double e[LD * LD];
    double df[LD * LD];
    lapack_complex_double ef[LD * LD];
    lapack_complex_double b[LD * LD];
    lapack_complex_double x[LD * LD];
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptrfs b", l, N, NRHS, b, LD, lapacke_test_region_full, -9,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zptrfs)(layout, 'U', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "zptrfs d", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -5,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zptrfs)(layout, 'U', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "zptrfs df", l, N, 1, df, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -7,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zptrfs)(layout, 'U', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptrfs e", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -6,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zptrfs)(layout, 'U', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptrfs ef", l, N - 1, 1, ef, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -8,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zptrfs)(layout, 'U', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptrfs x", l, N, NRHS, x, LD, lapacke_test_region_full, -11,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zptrfs)(layout, 'U', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_dfill_nan(layout, N, 1, df, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_zfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_zfill_nan(layout, N - 1, 1, ef,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_zfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "zptrfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zptrfs)(layout, 'U', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZPTRFS_ALLOC_TEST(0, 0, "zptrfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPTRFS_ALLOC_TEST(0, 1, "zptrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPTRFS_ALLOC_TEST(0, 2, "zptrfs allocation count", 0);
    lapacke_test_check_alloc_count("zptrfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPTRFS_ALLOC_TEST(1, 0, "zptrfs work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPTRFS_ALLOC_TEST(1, 1, "zptrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPTRFS_ALLOC_TEST(1, 2, "zptrfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPTRFS_ALLOC_TEST(1, 3, "zptrfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPTRFS_ALLOC_TEST(1, 4, "zptrfs allocation count", 0);
    lapacke_test_check_alloc_count("zptrfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPTRFS_ALLOC_TEST(2, 0, "zptrfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zptrfs invalid layout allocation count");
}
