#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

#define REGION_AB(i, j) lapacke_test_region_band(i, j, N, KL + KU)

/* Refill the inputs, schedule the malloc failure, call zgbcon. */
#define LAPACKE_ZGBCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab, LD);  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        anorm[0] = 1.0;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zgbcon)(layout, '1', N, KL, KU,  \
                                                      ab, LD, ipiv, anorm[0],  \
                                                      &rcond),                 \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zgbcon)
{
    lapack_complex_double ab[LD * LD];
    lapack_int ipiv[LD * LD];
    double anorm[1];
    double rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgbcon ab", l, 2 * KL + KU + 1, N, ab, LD, REGION_AB, -6,
            (lapacke_test_zfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zgbcon)(layout, '1', N, KL, KU, ab, LD, ipiv,
                                       anorm[0], &rcond));

        LAPACKE_TEST_DNAN_SWEEP(
            "zgbcon anorm", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -9,
            (lapacke_test_zfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0)),
            API_SUFFIX(LAPACKE_zgbcon)(layout, '1', N, KL, KU, ab, LD, ipiv,
                                       anorm[0], &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        anorm[0] = 1.0;
        lapacke_test_zfill_nan(layout, 2 * KL + KU + 1, N, ab, LD);
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_check(
            "zgbcon NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgbcon)(layout, '1', N, KL, KU, ab, LD, ipiv,
                                       anorm[0], &rcond) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZGBCON_ALLOC_TEST(0, 0, "zgbcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGBCON_ALLOC_TEST(0, 1, "zgbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGBCON_ALLOC_TEST(0, 2, "zgbcon allocation count", 0);
    lapacke_test_check_alloc_count("zgbcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZGBCON_ALLOC_TEST(1, 0, "zgbcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGBCON_ALLOC_TEST(1, 1, "zgbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGBCON_ALLOC_TEST(1, 2, "zgbcon transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGBCON_ALLOC_TEST(1, 3, "zgbcon allocation count", 0);
    lapacke_test_check_alloc_count("zgbcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZGBCON_ALLOC_TEST(2, 0, "zgbcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgbcon invalid layout allocation count");
}
