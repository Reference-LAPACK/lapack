#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_U_UNIT(i, j)                                                 \
    (lapacke_test_region_band(i, j, N, KD) && (i) != KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)
#define REGION_AB_L_UNIT(i, j)                                                 \
    (lapacke_test_region_band(i, j, N, 0) && (i) != 0)

/* Refill the inputs, schedule the malloc failure, call stbrfs. */
#define LAPACKE_STBRFS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_sfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_sfill_rhs(layout, N, NRHS, x, LD);                        \
        lapacke_test_sfill_vec(LD * LD, ferr);                                 \
        lapacke_test_sfill_vec(LD * LD, berr);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_stbrfs)(layout, 'U', 'N', 'N',   \
                                                      N, KD, NRHS, ab, LD, b,  \
                                                      LD, x, LD, ferr, berr),  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(stbrfs)
{
    float ab[LD * LD];
    float b[LD * LD];
    float x[LD * LD];
    float ferr[LD * LD];
    float berr[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs ab uplo=U diag=N", l, KD + 1, N, ab, LD, REGION_AB_U, -8,
            (lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'U', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs b uplo=U diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'U', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs x uplo=U diag=N", l, N, NRHS, x, LD,
            lapacke_test_region_full, -12,
            (lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'U', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs ab uplo=U diag=U", l, KD + 1, N, ab, LD, REGION_AB_U_UNIT,
            -8,
            (lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'U', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs b uplo=U diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'U', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs x uplo=U diag=U", l, N, NRHS, x, LD,
            lapacke_test_region_full, -12,
            (lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'U', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs ab uplo=L diag=N", l, KD + 1, N, ab, LD, REGION_AB_L, -8,
            (lapacke_test_sfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'L', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs b uplo=L diag=N", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_sfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'L', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs x uplo=L diag=N", l, N, NRHS, x, LD,
            lapacke_test_region_full, -12,
            (lapacke_test_sfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'L', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs ab uplo=L diag=U", l, KD + 1, N, ab, LD, REGION_AB_L_UNIT,
            -8,
            (lapacke_test_sfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'L', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs b uplo=L diag=U", l, N, NRHS, b, LD,
            lapacke_test_region_full, -10,
            (lapacke_test_sfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'L', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbrfs x uplo=L diag=U", l, N, NRHS, x, LD,
            lapacke_test_region_full, -12,
            (lapacke_test_sfill_tb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, x, LD),
             lapacke_test_sfill_vec(LD * LD, ferr),
             lapacke_test_sfill_vec(LD * LD, berr)),
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'L', 'N', 'U', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_vec(LD * LD, ferr);
        lapacke_test_sfill_vec(LD * LD, berr);
        lapacke_test_sfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_sfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_sfill_nan(layout, N, NRHS, x, LD);
        lapacke_test_check(
            "stbrfs NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_stbrfs)(layout, 'U', 'N', 'N', N, KD, NRHS, ab,
                                       LD, b, LD, x, LD, ferr, berr) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_STBRFS_ALLOC_TEST(0, 0, "stbrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_STBRFS_ALLOC_TEST(0, 1, "stbrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_STBRFS_ALLOC_TEST(0, 2, "stbrfs allocation count", 0);
    lapacke_test_check_alloc_count("stbrfs col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_STBRFS_ALLOC_TEST(1, 0, "stbrfs work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_STBRFS_ALLOC_TEST(1, 1, "stbrfs work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_STBRFS_ALLOC_TEST(1, 2, "stbrfs transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_STBRFS_ALLOC_TEST(1, 3, "stbrfs transpose alloc failure (b_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_STBRFS_ALLOC_TEST(1, 4, "stbrfs transpose alloc failure (x_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_STBRFS_ALLOC_TEST(1, 5, "stbrfs allocation count", 0);
    lapacke_test_check_alloc_count("stbrfs row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_STBRFS_ALLOC_TEST(2, 0, "stbrfs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("stbrfs invalid layout allocation count");
}
