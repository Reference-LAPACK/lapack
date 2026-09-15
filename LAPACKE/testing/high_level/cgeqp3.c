#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call cgeqp3. */
#define LAPACKE_CGEQP3_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, M, N, a, LD);                               \
        lapacke_test_fill_int(LD * LD, jpvt, 0);                               \
        lapacke_test_cfill_vec(LD * LD, tau);                                  \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_cgeqp3)(layout, M, N, a, LD, jpvt, tau),        \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(cgeqp3)
{
    lapack_complex_float a[LD * LD];
    lapack_int jpvt[LD * LD];
    lapack_complex_float tau[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cgeqp3 a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_cfill(layout, M, N, a, LD),
             lapacke_test_fill_int(LD * LD, jpvt, 0),
             lapacke_test_cfill_vec(LD * LD, tau)),
            API_SUFFIX(LAPACKE_cgeqp3)(layout, M, N, a, LD, jpvt, tau));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_fill_int(LD * LD, jpvt, 0);
        lapacke_test_cfill_vec(LD * LD, tau);
        lapacke_test_cfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "cgeqp3 NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cgeqp3)(layout, M, N, a, LD, jpvt, tau) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CGEQP3_ALLOC_TEST(0, 0, "cgeqp3 work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGEQP3_ALLOC_TEST(0, 1, "cgeqp3 work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGEQP3_ALLOC_TEST(0, 2, "cgeqp3 allocation count", 0);
    lapacke_test_check_alloc_count("cgeqp3 col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CGEQP3_ALLOC_TEST(1, 0, "cgeqp3 work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGEQP3_ALLOC_TEST(1, 1, "cgeqp3 work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGEQP3_ALLOC_TEST(1, 2, "cgeqp3 transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGEQP3_ALLOC_TEST(1, 3, "cgeqp3 allocation count", 0);
    lapacke_test_check_alloc_count("cgeqp3 row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CGEQP3_ALLOC_TEST(2, 0, "cgeqp3 invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cgeqp3 invalid layout allocation count");
}
