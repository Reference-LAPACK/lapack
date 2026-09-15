#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)

/* Refill the inputs, schedule the malloc failure, call spbcon. */
#define LAPACKE_SPBCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD);                     \
        anorm[0] = 1.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_spbcon)(layout, 'U', N, KD, ab,  \
                                                      LD, anorm[0], &rcond),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(spbcon)
{
    float ab[LD * LD];
    float anorm[1];
    float rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "spbcon ab uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -5,
            (lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD),
             (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_spbcon)(layout, 'U', N, KD, ab, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbcon anorm uplo=U", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -7,
            (lapacke_test_sfill_pb(layout, 'U', N, KD, ab, LD),
             (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_spbcon)(layout, 'U', N, KD, ab, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbcon ab uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -5,
            (lapacke_test_sfill_pb(layout, 'L', N, KD, ab, LD),
             (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_spbcon)(layout, 'L', N, KD, ab, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "spbcon anorm uplo=L", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -7,
            (lapacke_test_sfill_pb(layout, 'L', N, KD, ab, LD),
             (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_spbcon)(layout, 'L', N, KD, ab, LD, anorm[0],
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        anorm[0] = 1.0f;
        lapacke_test_sfill_nan(layout, KD + 1, N, ab, LD);
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_check("spbcon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_spbcon)(layout, 'U', N, KD, ab,
                                                      LD, anorm[0], &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SPBCON_ALLOC_TEST(0, 0, "spbcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPBCON_ALLOC_TEST(0, 1, "spbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPBCON_ALLOC_TEST(0, 2, "spbcon allocation count", 0);
    lapacke_test_check_alloc_count("spbcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPBCON_ALLOC_TEST(1, 0, "spbcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPBCON_ALLOC_TEST(1, 1, "spbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPBCON_ALLOC_TEST(1, 2, "spbcon transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPBCON_ALLOC_TEST(1, 3, "spbcon allocation count", 0);
    lapacke_test_check_alloc_count("spbcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPBCON_ALLOC_TEST(2, 0, "spbcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("spbcon invalid layout allocation count");
}
