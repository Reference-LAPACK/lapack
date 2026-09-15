#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)

/* Refill the inputs, schedule the malloc failure, call zpbcon. */
#define LAPACKE_ZPBCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD);                     \
        anorm[0] = 1.0;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zpbcon)(layout, 'U', N, KD, ab,  \
                                                      LD, anorm[0], &rcond),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zpbcon)
{
    lapack_complex_double ab[LD * LD];
    double anorm[1];
    double rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbcon ab uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -5,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zpbcon)(layout, 'U', N, KD, ab, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "zpbcon anorm uplo=U", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -7,
            (lapacke_test_zfill_pb(layout, 'U', N, KD, ab, LD),
             (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zpbcon)(layout, 'U', N, KD, ab, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zpbcon ab uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -5,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zpbcon)(layout, 'L', N, KD, ab, LD, anorm[0],
                                       &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "zpbcon anorm uplo=L", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -7,
            (lapacke_test_zfill_pb(layout, 'L', N, KD, ab, LD),
             (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zpbcon)(layout, 'L', N, KD, ab, LD, anorm[0],
                                       &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        anorm[0] = 1.0;
        lapacke_test_zfill_nan(layout, KD + 1, N, ab, LD);
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_check("zpbcon NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_zpbcon)(layout, 'U', N, KD, ab,
                                                      LD, anorm[0], &rcond) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZPBCON_ALLOC_TEST(0, 0, "zpbcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPBCON_ALLOC_TEST(0, 1, "zpbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPBCON_ALLOC_TEST(0, 2, "zpbcon allocation count", 0);
    lapacke_test_check_alloc_count("zpbcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPBCON_ALLOC_TEST(1, 0, "zpbcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPBCON_ALLOC_TEST(1, 1, "zpbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZPBCON_ALLOC_TEST(1, 2, "zpbcon transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPBCON_ALLOC_TEST(1, 3, "zpbcon allocation count", 0);
    lapacke_test_check_alloc_count("zpbcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPBCON_ALLOC_TEST(2, 0, "zpbcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zpbcon invalid layout allocation count");
}
