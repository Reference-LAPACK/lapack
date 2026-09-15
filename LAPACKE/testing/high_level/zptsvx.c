#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zptsvx. */
#define LAPACKE_ZPTSVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
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
                           API_SUFFIX(LAPACKE_zptsvx)(layout, 'N', N, NRHS, d, \
                                                      e, df, ef, b, LD, x, LD, \
                                                      &rcond, ferr, berr),     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zptsvx)
{
    double d[LD * LD];
    lapack_complex_double e[LD * LD];
    double df[LD * LD];
    lapack_complex_double ef[LD * LD];
    lapack_complex_double b[LD * LD];
    lapack_complex_double x[LD * LD];
    double rcond;
    double ferr[LD * LD];
    double berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptsvx b fact=N", l, N, NRHS, b, LD, lapacke_test_region_full, -9,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zptsvx)(layout, 'N', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "zptsvx d fact=N", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -5,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zptsvx)(layout, 'N', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptsvx e fact=N", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -6,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_zptsvx)(layout, 'N', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptsvx b fact=F", l, N, NRHS, b, LD, lapacke_test_region_full, -9,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_zptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "zptsvx d fact=F", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -5,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_zptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_DNAN_SWEEP(
            "zptsvx df fact=F", l, N, 1, df, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -7,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_zptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptsvx e fact=F", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -6,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_zptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zptsvx ef fact=F", l, N - 1, 1, ef,
            LAPACKE_TEST_VLD(layout, N - 1), lapacke_test_region_full, -8,
            (lapacke_test_dfill_pos(LD * LD, d),
             lapacke_test_zfill_vec(LD * LD, e),
             lapacke_test_dfill_pos(LD * LD, df),
             lapacke_test_zfill_vec(LD * LD, ef),
             lapacke_test_zfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_zfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_dfill_vec(LD * LD, ferr),
             lapacke_test_dfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_zpttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_zptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_dfill_vec(LD * LD, ferr);
        lapacke_test_dfill_vec(LD * LD, berr);
        lapacke_test_zfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_dfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_zfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_dfill_nan(layout, N, 1, df, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_zfill_nan(layout, N - 1, 1, ef,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check(
            "zptsvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zptsvx)(layout, 'N', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZPTSVX_ALLOC_TEST(0, 0, "zptsvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPTSVX_ALLOC_TEST(0, 1, "zptsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPTSVX_ALLOC_TEST(0, 2, "zptsvx allocation count", 0);
    lapacke_test_check_alloc_count("zptsvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPTSVX_ALLOC_TEST(1, 0, "zptsvx work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPTSVX_ALLOC_TEST(1, 1, "zptsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPTSVX_ALLOC_TEST(1, 2, "zptsvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPTSVX_ALLOC_TEST(1, 3, "zptsvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPTSVX_ALLOC_TEST(1, 4, "zptsvx allocation count", 0);
    lapacke_test_check_alloc_count("zptsvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPTSVX_ALLOC_TEST(2, 0, "zptsvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zptsvx invalid layout allocation count");
}
