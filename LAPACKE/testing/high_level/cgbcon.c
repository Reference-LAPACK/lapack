#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

#define REGION_AB(i, j) lapacke_test_region_band(i, j, N, KL + KU)

/* Refill the inputs, schedule the malloc failure, call cgbcon. */
#define LAPACKE_CGBCON_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab, LD);  \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        anorm[0] = 1.0f;                                                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cgbcon)(layout, '1', N, KL, KU,  \
                                                      ab, LD, ipiv, anorm[0],  \
                                                      &rcond),                 \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cgbcon)
{
    lapack_complex_float ab[LD * LD];
    lapack_int ipiv[LD * LD];
    float anorm[1];
    float rcond;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cgbcon ab", l, 2 * KL + KU + 1, N, ab, LD, REGION_AB, -6,
            (lapacke_test_cfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cgbcon)(layout, '1', N, KL, KU, ab, LD, ipiv,
                                       anorm[0], &rcond));

        LAPACKE_TEST_SNAN_SWEEP(
            "cgbcon anorm", l, 1, 1, anorm, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -9,
            (lapacke_test_cfill_gb(layout, N, N, KL, KU, 2 * KL + KU + 1, ab,
                                   LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv), (anorm[0] = 1.0f)),
            API_SUFFIX(LAPACKE_cgbcon)(layout, '1', N, KL, KU, ab, LD, ipiv,
                                       anorm[0], &rcond));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        anorm[0] = 1.0f;
        lapacke_test_cfill_nan(layout, 2 * KL + KU + 1, N, ab, LD);
        /* anorm stays finite: LAPACK rejects it as an illegal argument */
        lapacke_test_check(
            "cgbcon NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgbcon)(layout, '1', N, KL, KU, ab, LD, ipiv,
                                       anorm[0], &rcond) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CGBCON_ALLOC_TEST(0, 0, "cgbcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGBCON_ALLOC_TEST(0, 1, "cgbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGBCON_ALLOC_TEST(0, 2, "cgbcon allocation count", 0);
    lapacke_test_check_alloc_count("cgbcon col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CGBCON_ALLOC_TEST(1, 0, "cgbcon work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGBCON_ALLOC_TEST(1, 1, "cgbcon work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGBCON_ALLOC_TEST(1, 2, "cgbcon transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGBCON_ALLOC_TEST(1, 3, "cgbcon allocation count", 0);
    lapacke_test_check_alloc_count("cgbcon row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CGBCON_ALLOC_TEST(2, 0, "cgbcon invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cgbcon invalid layout allocation count");
}
