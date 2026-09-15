#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call chptri. */
#define LAPACKE_CHPTRI_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_sp(layout, 'U', N, ap);                             \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_chptri)(layout, 'U', N, ap, ipiv), expected);   \
    } while (0)

LAPACKE_TEST(chptri)
{
    lapack_complex_float ap[LD * LD];
    lapack_int ipiv[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "chptri ap uplo=U", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_cfill_sp(layout, 'U', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_chptri)(layout, 'U', N, ap, ipiv));

        LAPACKE_TEST_CNAN_SWEEP(
            "chptri ap uplo=L", l, N * (N + 1) / 2, 1, ap,
            LAPACKE_TEST_VLD(layout, N * (N + 1) / 2), lapacke_test_region_full,
            -4,
            (lapacke_test_cfill_sp(layout, 'L', N, ap),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_chptri)(layout, 'L', N, ap, ipiv));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_cfill_nan(layout, N * (N + 1) / 2, 1, ap,
                               LAPACKE_TEST_VLD(layout, N * (N + 1) / 2));
        lapacke_test_check(
            "chptri NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_chptri)(layout, 'U', N, ap, ipiv) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CHPTRI_ALLOC_TEST(0, 0, "chptri work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHPTRI_ALLOC_TEST(0, 1, "chptri allocation count", 0);
    lapacke_test_check_alloc_count("chptri col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CHPTRI_ALLOC_TEST(1, 0, "chptri work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHPTRI_ALLOC_TEST(1, 1, "chptri transpose alloc failure (ap_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CHPTRI_ALLOC_TEST(1, 2, "chptri allocation count", 0);
    lapacke_test_check_alloc_count("chptri row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CHPTRI_ALLOC_TEST(2, 0, "chptri invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("chptri invalid layout allocation count");
}
