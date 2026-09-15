#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zhptrf. */
#define LAPACKE_ZHPTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_sp(layout, 'U', N, ap);                             \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zhptrf)(layout, 'U', N, ap, ipiv), expected);   \
    } while (0)

LAPACKE_TEST(zhptrf)
{
    lapack_complex_double ap[LD * LD];
    lapack_int ipiv[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zhptrf ap uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_zfill_sp(layout, 'U', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_zhptrf)(layout, 'U', N, ap, ipiv));

        LAPACKE_TEST_ZNAN_SWEEP(
            "zhptrf ap uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_zfill_sp(layout, 'L', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_zhptrf)(layout, 'L', N, ap, ipiv));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_zfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_check(
            "zhptrf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zhptrf)(layout, 'U', N, ap, ipiv) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: no allocation at all */
    LAPACKE_ZHPTRF_ALLOC_TEST(0, 0, "zhptrf allocation count", 0);
    lapacke_test_check_alloc_count("zhptrf col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZHPTRF_ALLOC_TEST(1, 0, "zhptrf transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZHPTRF_ALLOC_TEST(1, 1, "zhptrf allocation count", 0);
    lapacke_test_check_alloc_count("zhptrf row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZHPTRF_ALLOC_TEST(2, 0, "zhptrf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zhptrf invalid layout allocation count");
}
