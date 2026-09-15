#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

/* Refill the inputs, schedule the malloc failure, call slascl. */
#define LAPACKE_SLASCL_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        cfrom[0] = 1.0f;                                                       \
        cto[0] = 2.0f;                                                         \
        lapacke_test_sfill(layout, M, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_slascl)(layout, 'G', KL, KU,     \
                                                      cfrom[0], cto[0], M, N,  \
                                                      a, LD),                  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(slascl)
{
    float cfrom[1];
    float cto[1];
    float a[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "slascl a", l, M, N, a, LD, lapacke_test_region_full, -9,
            ((cfrom[0] = 1.0f), (cto[0] = 2.0f),
             lapacke_test_sfill(layout, M, N, a, LD)),
            API_SUFFIX(LAPACKE_slascl)(layout, 'G', KL, KU, cfrom[0], cto[0], M,
                                       N, a, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        cfrom[0] = 1.0f;
        cto[0] = 2.0f;
        lapacke_test_sfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "slascl NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_slascl)(layout, 'G', KL, KU, cfrom[0], cto[0], M,
                                       N, a, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_SLASCL_ALLOC_TEST(0, 0, "slascl allocation count", 0);
    lapacke_test_check_alloc_count("slascl col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SLASCL_ALLOC_TEST(1, 0, "slascl transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SLASCL_ALLOC_TEST(1, 1, "slascl allocation count", 0);
    lapacke_test_check_alloc_count("slascl row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SLASCL_ALLOC_TEST(2, 0, "slascl invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("slascl invalid layout allocation count");
}
