#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KD 1

#define REGION_AB_U(i, j) lapacke_test_region_band(i, j, N, KD)
#define REGION_AB_L(i, j) lapacke_test_region_band(i, j, N, 0)

/* Refill the inputs, schedule the malloc failure, call cpbequ. */
#define LAPACKE_CPBEQU_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD);                     \
        lapacke_test_sfill_pos(LD * LD, s);                                    \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cpbequ)(layout, 'U', N, KD, ab,  \
                                                      LD, s, &scond, &amax),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cpbequ)
{
    lapack_complex_float ab[LD * LD];
    float s[LD * LD];
    float scond;
    float amax;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbequ ab uplo=U", l, KD + 1, N, ab, LD, REGION_AB_U, -5,
            (lapacke_test_cfill_pb(layout, 'U', N, KD, ab, LD),
             lapacke_test_sfill_pos(LD * LD, s)),
            API_SUFFIX(LAPACKE_cpbequ)(layout, 'U', N, KD, ab, LD, s, &scond,
                                       &amax));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpbequ ab uplo=L", l, KD + 1, N, ab, LD, REGION_AB_L, -5,
            (lapacke_test_cfill_pb(layout, 'L', N, KD, ab, LD),
             lapacke_test_sfill_pos(LD * LD, s)),
            API_SUFFIX(LAPACKE_cpbequ)(layout, 'L', N, KD, ab, LD, s, &scond,
                                       &amax));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_pos(LD * LD, s);
        lapacke_test_cfill_nan(layout, KD + 1, N, ab, LD);
        lapacke_test_check("cpbequ NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_cpbequ)(layout, 'U', N, KD, ab,
                                                      LD, s, &scond, &amax) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CPBEQU_ALLOC_TEST(0, 0, "cpbequ allocation count", 0);
    lapacke_test_check_alloc_count("cpbequ col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CPBEQU_ALLOC_TEST(1, 0, "cpbequ transpose alloc failure (ab_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CPBEQU_ALLOC_TEST(1, 1, "cpbequ allocation count", 0);
    lapacke_test_check_alloc_count("cpbequ row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CPBEQU_ALLOC_TEST(2, 0, "cpbequ invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cpbequ invalid layout allocation count");
}
