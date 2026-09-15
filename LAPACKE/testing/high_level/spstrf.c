#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call spstrf. */
#define LAPACKE_SPSTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill_spd(layout, N, a, LD);                              \
        lapacke_test_fill_int(LD * LD, piv, 0);                                \
        tol[0] = -1.0f;                                                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_spstrf)(layout, 'U', N, a, LD,   \
                                                      piv, &rank, tol[0]),     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(spstrf)
{
    float a[LD * LD];
    lapack_int piv[LD * LD];
    lapack_int rank;
    float tol[1];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "spstrf a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -4,
            (lapacke_test_sfill_spd(layout, N, a, LD),
             lapacke_test_fill_int(LD * LD, piv, 0), (tol[0] = -1.0f)),
            API_SUFFIX(LAPACKE_spstrf)(layout, 'U', N, a, LD, piv, &rank,
                                       tol[0]));

        LAPACKE_TEST_SNAN_SWEEP(
            "spstrf tol uplo=U", l, 1, 1, tol, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -8,
            (lapacke_test_sfill_spd(layout, N, a, LD),
             lapacke_test_fill_int(LD * LD, piv, 0), (tol[0] = -1.0f)),
            API_SUFFIX(LAPACKE_spstrf)(layout, 'U', N, a, LD, piv, &rank,
                                       tol[0]));

        LAPACKE_TEST_SNAN_SWEEP(
            "spstrf a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -4,
            (lapacke_test_sfill_spd(layout, N, a, LD),
             lapacke_test_fill_int(LD * LD, piv, 0), (tol[0] = -1.0f)),
            API_SUFFIX(LAPACKE_spstrf)(layout, 'L', N, a, LD, piv, &rank,
                                       tol[0]));

        LAPACKE_TEST_SNAN_SWEEP(
            "spstrf tol uplo=L", l, 1, 1, tol, LAPACKE_TEST_VLD(layout, 1),
            lapacke_test_region_full, -8,
            (lapacke_test_sfill_spd(layout, N, a, LD),
             lapacke_test_fill_int(LD * LD, piv, 0), (tol[0] = -1.0f)),
            API_SUFFIX(LAPACKE_spstrf)(layout, 'L', N, a, LD, piv, &rank,
                                       tol[0]));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_int(LD * LD, piv, 0);
        lapacke_test_sfill_nan(layout, N, N, a, LD);
        lapacke_test_sfill_nan(layout, 1, 1, tol, LAPACKE_TEST_VLD(layout, 1));
        lapacke_test_check("spstrf NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_spstrf)(layout, 'U', N, a, LD,
                                                      piv, &rank, tol[0]) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SPSTRF_ALLOC_TEST(0, 0, "spstrf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPSTRF_ALLOC_TEST(0, 1, "spstrf allocation count", 0);
    lapacke_test_check_alloc_count("spstrf col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SPSTRF_ALLOC_TEST(1, 0, "spstrf work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SPSTRF_ALLOC_TEST(1, 1, "spstrf transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SPSTRF_ALLOC_TEST(1, 2, "spstrf allocation count", 0);
    lapacke_test_check_alloc_count("spstrf row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SPSTRF_ALLOC_TEST(2, 0, "spstrf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("spstrf invalid layout allocation count");
}
