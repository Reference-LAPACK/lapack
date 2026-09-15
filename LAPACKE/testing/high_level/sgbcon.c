#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

#define REGION_AB(i, j) lapacke_test_region_band(i, j, N, KL + KU)

/* Refill the inputs, schedule the malloc failure, call sgbcon. */
#define LAPACKE_SGBCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab, LD);  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        anorm[0] = 1.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sgbcon)(layout, '1', N, KL, KU,  \
                                                      ab, LD, ipiv, anorm[0],  \
                                                      &rcond),                 \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sgbcon)
{
    float ab[LD * LD];
    lapack_int ipiv[LD * LD];
    float anorm[1];
    float rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgbcon ab", l, 2 * KL + KU + 1, N, ab, LD, REGION_AB, -6,
            (lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_sgbcon)(layout, '1', N, KL, KU, ab, LD, ipiv,
                                       anorm[0], &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgbcon anorm", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -9,
            (lapacke_test_sfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_sgbcon)(layout, '1', N, KL, KU, ab, LD, ipiv,
                                       anorm[0], &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        anorm[0] = 1.0f;
        lapacke_test_sfill_nan(layout, 2 * KL + KU + 1, N, ab, LD);
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_check(
            "sgbcon NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sgbcon)(layout, '1', N, KL, KU, ab, LD, ipiv,
                                       anorm[0], &rcond) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SGBCON_ALLOC_TEST(0, 0, "sgbcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGBCON_ALLOC_TEST(0, 1, "sgbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGBCON_ALLOC_TEST(0, 2, "sgbcon allocation count", 0);
    lapacke_test_check_alloc_count("sgbcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SGBCON_ALLOC_TEST(1, 0, "sgbcon work alloc failure (iwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGBCON_ALLOC_TEST(1, 1, "sgbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGBCON_ALLOC_TEST(1, 2, "sgbcon transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGBCON_ALLOC_TEST(1, 3, "sgbcon allocation count", 0);
    lapacke_test_check_alloc_count("sgbcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SGBCON_ALLOC_TEST(2, 0, "sgbcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgbcon invalid layout allocation count");
}
