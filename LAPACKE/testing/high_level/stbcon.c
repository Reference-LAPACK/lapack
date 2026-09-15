#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_U_UNIT(i, j)                                                 \
    (lapacke_test_region_band(i, j, N, KD) && (i) != KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)
#define REGION_AB_L_UNIT(i, j)                                                 \
    (lapacke_test_region_band(i, j, N, 0) && (i) != 0)

/* Refill the inputs, schedule the malloc failure, call stbcon. */
#define LAPACKE_STBCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_stbcon)(layout, '1', 'U', 'N',   \
                                                      N, KD, ab, LD, &rcond),  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(stbcon)
{
    float ab[LD * LD];
    float rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "stbcon ab uplo=U diag=N", l, KD + 1, N, ab, LD, REGION_AB_U, -7,
            (lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD)),
            API_SUFFIX(LAPACKE_stbcon)(layout, '1', 'U', 'N', N, KD, ab, LD,
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbcon ab uplo=U diag=U", l, KD + 1, N, ab, LD, REGION_AB_U_UNIT,
            -7, (lapacke_test_sfill_tb(layout, 'U', N, KD, ab, LD)),
            API_SUFFIX(LAPACKE_stbcon)(layout, '1', 'U', 'U', N, KD, ab, LD,
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbcon ab uplo=L diag=N", l, KD + 1, N, ab, LD, REGION_AB_L, -7,
            (lapacke_test_sfill_tb(layout, 'L', N, KD, ab, LD)),
            API_SUFFIX(LAPACKE_stbcon)(layout, '1', 'L', 'N', N, KD, ab, LD,
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "stbcon ab uplo=L diag=U", l, KD + 1, N, ab, LD, REGION_AB_L_UNIT,
            -7, (lapacke_test_sfill_tb(layout, 'L', N, KD, ab, LD)),
            API_SUFFIX(LAPACKE_stbcon)(layout, '1', 'L', 'U', N, KD, ab, LD,
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_check("stbcon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_stbcon)(layout, '1', 'U', 'N', N,
                                                      KD, ab, LD, &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_STBCON_ALLOC_TEST(0, 0, "stbcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_STBCON_ALLOC_TEST(0, 1, "stbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_STBCON_ALLOC_TEST(0, 2, "stbcon allocation count", 0);
    lapacke_test_check_alloc_count("stbcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_STBCON_ALLOC_TEST(1, 0, "stbcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_STBCON_ALLOC_TEST(1, 1, "stbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_STBCON_ALLOC_TEST(1, 2, "stbcon transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_STBCON_ALLOC_TEST(1, 3, "stbcon allocation count", 0);
    lapacke_test_check_alloc_count("stbcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_STBCON_ALLOC_TEST(2, 0, "stbcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("stbcon invalid layout allocation count");
}
