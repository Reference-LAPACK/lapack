#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define K N

/* Refill the inputs, schedule the malloc failure, call sormqr. */
#define LAPACKE_SORMQR_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, M, K, a, LD);                               \
        lapacke_test_sfill_vec(LD * LD, tau);                                  \
        lapacke_test_sfill_rhs(layout, M, N, c, LD);                           \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sormqr)(layout, 'L', 'N', M, N,  \
                                                      K, a, LD, tau, c, LD),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sormqr)
{
    float a[LD * LD];
    float tau[LD * LD];
    float c[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sormqr a side=L", l, M, K, a, LD, lapacke_test_region_full, -7,
            (lapacke_test_sfill(layout, M, K, a, LD),
             lapacke_test_sfill_vec(LD * LD, tau),
             lapacke_test_sfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_sormqr)(layout, 'L', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sormqr c side=L", l, M, N, c, LD, lapacke_test_region_full, -10,
            (lapacke_test_sfill(layout, M, K, a, LD),
             lapacke_test_sfill_vec(LD * LD, tau),
             lapacke_test_sfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_sormqr)(layout, 'L', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sormqr tau side=L", l, K, 1, tau, LAPACKE_TEST_VLD(layout, K),
            lapacke_test_region_full, -9,
            (lapacke_test_sfill(layout, M, K, a, LD),
             lapacke_test_sfill_vec(LD * LD, tau),
             lapacke_test_sfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_sormqr)(layout, 'L', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sormqr a side=R", l, N, K, a, LD, lapacke_test_region_full, -7,
            (lapacke_test_sfill(layout, N, K, a, LD),
             lapacke_test_sfill_vec(LD * LD, tau),
             lapacke_test_sfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_sormqr)(layout, 'R', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sormqr c side=R", l, M, N, c, LD, lapacke_test_region_full, -10,
            (lapacke_test_sfill(layout, N, K, a, LD),
             lapacke_test_sfill_vec(LD * LD, tau),
             lapacke_test_sfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_sormqr)(layout, 'R', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sormqr tau side=R", l, K, 1, tau, LAPACKE_TEST_VLD(layout, K),
            lapacke_test_region_full, -9,
            (lapacke_test_sfill(layout, N, K, a, LD),
             lapacke_test_sfill_vec(LD * LD, tau),
             lapacke_test_sfill_rhs(layout, M, N, c, LD)),
            API_SUFFIX(LAPACKE_sormqr)(layout, 'R', 'N', M, N, K, a, LD, tau, c,
                                       LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, M, K, a, LD);
        lapacke_test_sfill_nan(layout, M, N, c, LD);
        lapacke_test_sfill_nan(layout, K, 1, tau, LAPACKE_TEST_VLD(layout, K));
        lapacke_test_check("sormqr NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_sormqr)(layout, 'L', 'N', M, N, K,
                                                      a, LD, tau, c, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SORMQR_ALLOC_TEST(0, 0, "sormqr work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SORMQR_ALLOC_TEST(0, 1, "sormqr allocation count", 0);
    lapacke_test_check_alloc_count("sormqr col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SORMQR_ALLOC_TEST(1, 0, "sormqr work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SORMQR_ALLOC_TEST(1, 1, "sormqr transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SORMQR_ALLOC_TEST(1, 2, "sormqr transpose alloc failure (c_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SORMQR_ALLOC_TEST(1, 3, "sormqr allocation count", 0);
    lapacke_test_check_alloc_count("sormqr row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SORMQR_ALLOC_TEST(2, 0, "sormqr invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sormqr invalid layout allocation count");
}
