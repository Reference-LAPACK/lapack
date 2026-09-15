#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zgetf2. */
#define LAPACKE_ZGETF2_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, M, N, a, LD);                               \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zgetf2)(layout, M, N, a, LD, ipiv), expected);  \
    } while (0)

LAPACKE_TEST(zgetf2)
{
    lapack_complex_double a[LD * LD];
    lapack_int ipiv[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgetf2 a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_zfill(layout, M, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_zgetf2)(layout, M, N, a, LD, ipiv));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_zfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "zgetf2 NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgetf2)(layout, M, N, a, LD, ipiv) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZGETF2_ALLOC_TEST(0, 0, "zgetf2 allocation count", 0);
    lapacke_test_check_alloc_count("zgetf2 col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZGETF2_ALLOC_TEST(1, 0, "zgetf2 transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGETF2_ALLOC_TEST(1, 1, "zgetf2 allocation count", 0);
    lapacke_test_check_alloc_count("zgetf2 row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZGETF2_ALLOC_TEST(2, 0, "zgetf2 invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgetf2 invalid layout allocation count");
}
