#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sptsvx. */
#define LAPACKE_SPTSVX_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_pos(LD * LD, d);                                    \
        lapacke_test_sfill_vec(LD * LD, e);                                    \
        lapacke_test_sfill_pos(LD * LD, df);                                   \
        lapacke_test_sfill_vec(LD * LD, ef);                                   \
        lapacke_test_sfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sptsvx)(layout, 'N', N, NRHS, d, \
                                                      e, df, ef, b, LD, x, LD, \
                                                      &rcond, ferr, berr),     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sptsvx)
{
    float d[LD * LD];
    float e[LD * LD];
    float df[LD * LD];
    float ef[LD * LD];
    float b[LD * LD];
    float x[LD * LD];
    float rcond;
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsvx b fact=N", l, N, NRHS, b, LD, lapacke_test_region_full, -9,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_pos(LD * LD, df),
             lapacke_test_sfill_vec(LD * LD, ef),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sptsvx)(layout, 'N', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsvx d fact=N", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -5,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_pos(LD * LD, df),
             lapacke_test_sfill_vec(LD * LD, ef),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sptsvx)(layout, 'N', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsvx e fact=N", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -6,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_pos(LD * LD, df),
             lapacke_test_sfill_vec(LD * LD, ef),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_sptsvx)(layout, 'N', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsvx b fact=F", l, N, NRHS, b, LD, lapacke_test_region_full, -9,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_pos(LD * LD, df),
             lapacke_test_sfill_vec(LD * LD, ef),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_sptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsvx d fact=F", l, N, 1, d, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -5,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_pos(LD * LD, df),
             lapacke_test_sfill_vec(LD * LD, ef),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_sptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsvx df fact=F", l, N, 1, df, LAPACKE_TEST_VLD(layout, N),
            lapacke_test_region_full, -7,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_pos(LD * LD, df),
             lapacke_test_sfill_vec(LD * LD, ef),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_sptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsvx e fact=F", l, N - 1, 1, e, LAPACKE_TEST_VLD(layout, N - 1),
            lapacke_test_region_full, -6,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_pos(LD * LD, df),
             lapacke_test_sfill_vec(LD * LD, ef),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_sptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "sptsvx ef fact=F", l, N - 1, 1, ef,
            LAPACKE_TEST_VLD(layout, N - 1), lapacke_test_region_full, -8,
            (lapacke_test_sfill_pos(LD * LD, d),
             lapacke_test_sfill_vec(LD * LD, e),
             lapacke_test_sfill_pos(LD * LD, df),
             lapacke_test_sfill_vec(LD * LD, ef),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr),
             API_SUFFIX(LAPACKE_spttrf)(N, df, ef)),
            API_SUFFIX(LAPACKE_sptsvx)(layout, 'F', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_sfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_sfill_nan(layout, N, 1, d, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_sfill_nan(layout, N - 1, 1, e,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_sfill_nan(layout, N, 1, df, LAPACKE_TEST_VLD(layout, N));
        lapacke_test_sfill_nan(layout, N - 1, 1, ef,
                               LAPACKE_TEST_VLD(layout, N - 1));
        lapacke_test_check(
            "sptsvx NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sptsvx)(layout, 'N', N, NRHS, d, e, df, ef, b,
                                       LD, x, LD, &rcond, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SPTSVX_ALLOC_TEST(0, 0, "sptsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPTSVX_ALLOC_TEST(0, 1, "sptsvx allocation count", 0);
    lapacke_test_check_alloc_count("sptsvx col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPTSVX_ALLOC_TEST(1, 0, "sptsvx work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPTSVX_ALLOC_TEST(1, 1, "sptsvx transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPTSVX_ALLOC_TEST(1, 2, "sptsvx transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPTSVX_ALLOC_TEST(1, 3, "sptsvx allocation count", 0);
    lapacke_test_check_alloc_count("sptsvx row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPTSVX_ALLOC_TEST(2, 0, "sptsvx invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sptsvx invalid layout allocation count");
}
