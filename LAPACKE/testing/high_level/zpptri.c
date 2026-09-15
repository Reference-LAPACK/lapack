#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zpptri. */
#define LAPACKE_ZPPTRI_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_pp(layout, 'U', N, ap);                             \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zpptri)(layout, 'U', N, ap),     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zpptri)
{
    lapack_complex_double ap[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP("zpptri ap uplo=U", l, N * (N + 1) / 2, 1, ap,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -4,
                                (lapacke_test_zfill_pp(layout, 'U', N, ap)),
                                API_SUFFIX(LAPACKE_zpptri)(layout, 'U', N, ap));

        LAPACKE_TEST_ZNAN_SWEEP("zpptri ap uplo=L", l, N * (N + 1) / 2, 1, ap,
                                LAPACKE_TEST_VLD(layout, N * (N + 1) / 2),
                                lapacke_test_region_full, -4,
                                (lapacke_test_zfill_pp(layout, 'L', N, ap)),
                                API_SUFFIX(LAPACKE_zpptri)(layout, 'L', N, ap));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_check(
            "zpptri NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zpptri)(layout, 'U', N, ap) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZPPTRI_ALLOC_TEST(0, 0, "zpptri allocation count", 0);
    lapacke_test_check_alloc_count("zpptri col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZPPTRI_ALLOC_TEST(1, 0, "zpptri transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZPPTRI_ALLOC_TEST(1, 1, "zpptri allocation count", 0);
    lapacke_test_check_alloc_count("zpptri row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZPPTRI_ALLOC_TEST(2, 0, "zpptri invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zpptri invalid layout allocation count");
}
