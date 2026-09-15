#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sgetri. */
#define LAPACKE_SGETRI_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, N, N, a, LD);                               \
        lapacke_test_fill_ipiv(LD * LD, ipiv);                                 \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sgetri)(layout, N, a, LD, ipiv), \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sgetri)
{
    float a[LD * LD];
    lapack_int ipiv[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgetri a", l, N, N, a, LD, lapacke_test_region_full, -3,
            (lapacke_test_sfill(layout, N, N, a, LD),
             lapacke_test_fill_ipiv(LD * LD, ipiv)),
            API_SUFFIX(LAPACKE_sgetri)(layout, N, a, LD, ipiv));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_ipiv(LD * LD, ipiv);
        lapacke_test_sfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "sgetri NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sgetri)(layout, N, a, LD, ipiv) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SGETRI_ALLOC_TEST(0, 0, "sgetri work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGETRI_ALLOC_TEST(0, 1, "sgetri allocation count", 0);
    lapacke_test_check_alloc_count("sgetri col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SGETRI_ALLOC_TEST(1, 0, "sgetri work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGETRI_ALLOC_TEST(1, 1, "sgetri transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGETRI_ALLOC_TEST(1, 2, "sgetri allocation count", 0);
    lapacke_test_check_alloc_count("sgetri row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SGETRI_ALLOC_TEST(2, 0, "sgetri invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgetri invalid layout allocation count");
}
