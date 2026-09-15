#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

#define REGION_AB(i, j) lapacke_test_region_band(i, j, M, KU)

/* Refill the inputs, schedule the malloc failure, call dgbequ. */
#define LAPACKE_DGBEQU_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_gb(layout, M, N, KL, KU, KL + KU + 1, ab, LD);      \
        lapacke_test_dfill_pos(LD * LD, r);                                    \
        lapacke_test_dfill_pos(LD * LD, c);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dgbequ)(layout, M, N, KL, KU,    \
                                                      ab, LD, r, c, &rowcnd,   \
                                                      &colcnd, &amax),         \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dgbequ)
{
    double ab[LD * LD];
    double r[LD * LD];
    double c[LD * LD];
    double rowcnd;
    double colcnd;
    double amax;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgbequ ab", l, KL + KU + 1, N, ab, LD, REGION_AB, -6,
            (lapacke_test_dfill_gb(layout, M, N, KL, KU, KL + KU + 1, ab, LD),
             lapacke_test_dfill_pos(LD * LD, r),
             lapacke_test_dfill_pos(LD * LD, c)),
            API_SUFFIX(LAPACKE_dgbequ)(layout, M, N, KL, KU, ab, LD, r, c,
                                       &rowcnd, &colcnd, &amax));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_pos(LD * LD, r);
        lapacke_test_dfill_pos(LD * LD, c);
        lapacke_test_dfill_nan(layout, KL + KU + 1, N, ab, LD);
        lapacke_test_check(
            "dgbequ NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgbequ)(layout, M, N, KL, KU, ab, LD, r, c,
                                       &rowcnd, &colcnd, &amax) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_DGBEQU_ALLOC_TEST(0, 0, "dgbequ allocation count", 0);
    lapacke_test_check_alloc_count("dgbequ col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGBEQU_ALLOC_TEST(1, 0, "dgbequ transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGBEQU_ALLOC_TEST(1, 1, "dgbequ allocation count", 0);
    lapacke_test_check_alloc_count("dgbequ row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGBEQU_ALLOC_TEST(2, 0, "dgbequ invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgbequ invalid layout allocation count");
}
