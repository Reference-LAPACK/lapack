#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define KL 1
#define KU 1

/* Refill the inputs, schedule the malloc failure, call zlascl. */
#define LAPACKE_ZLASCL_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        cfrom[0] = 1.0;                                                        \
        cto[0] = 2.0;                                                          \
        lapacke_test_zfill(layout, M, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zlascl)(layout, 'G', KL, KU,     \
                                                      cfrom[0], cto[0], M, N,  \
                                                      a, LD),                  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zlascl)
{
    double cfrom[1];
    double cto[1];
    lapack_complex_double a[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zlascl a", l, M, N, a, LD, lapacke_test_region_full, -9,
            ((cfrom[0] = 1.0), (cto[0] = 2.0),
             lapacke_test_zfill(layout, M, N, a, LD)),
            API_SUFFIX(LAPACKE_zlascl)(layout, 'G', KL, KU, cfrom[0], cto[0], M,
                                       N, a, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        cfrom[0] = 1.0;
        cto[0] = 2.0;
        lapacke_test_zfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "zlascl NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zlascl)(layout, 'G', KL, KU, cfrom[0], cto[0], M,
                                       N, a, LD) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZLASCL_ALLOC_TEST(0, 0, "zlascl allocation count", 0);
    lapacke_test_check_alloc_count("zlascl col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZLASCL_ALLOC_TEST(1, 0, "zlascl transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZLASCL_ALLOC_TEST(1, 1, "zlascl allocation count", 0);
    lapacke_test_check_alloc_count("zlascl row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZLASCL_ALLOC_TEST(2, 0, "zlascl invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zlascl invalid layout allocation count");
}
