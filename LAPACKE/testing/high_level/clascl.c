#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

/* Refill the inputs, schedule the malloc failure, call clascl. */
#define LAPACKE_CLASCL_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        cfrom[0] = 1.0f;                                                       \
        cto[0] = 2.0f;                                                         \
        lapacke_test_cfill(layout, M, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_clascl)(layout, 'G', KL, KU,     \
                                                      cfrom[0], cto[0], M, N,  \
                                                      a, LD),                  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(clascl)
{
    float cfrom[1];
    float cto[1];
    lapack_complex_float a[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "clascl a", l, M, N, a, LD, lapacke_test_region_full, -9,
            ((cfrom[0] = 1.0f), (cto[0] = 2.0f),
             lapacke_test_cfill(layout, M, N, a, LD)),
            API_SUFFIX(LAPACKE_clascl)(layout, 'G', KL, KU, cfrom[0], cto[0], M,
                                       N, a, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        cfrom[0] = 1.0f;
        cto[0] = 2.0f;
        lapacke_test_cfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "clascl NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_clascl)(layout, 'G', KL, KU, cfrom[0], cto[0], M,
                                       N, a, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_CLASCL_ALLOC_TEST(0, 0, "clascl allocation count", 0);
    lapacke_test_check_alloc_count("clascl col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CLASCL_ALLOC_TEST(1, 0, "clascl transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CLASCL_ALLOC_TEST(1, 1, "clascl allocation count", 0);
    lapacke_test_check_alloc_count("clascl row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CLASCL_ALLOC_TEST(2, 0, "clascl invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("clascl invalid layout allocation count");
}
