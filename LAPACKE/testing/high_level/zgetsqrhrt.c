#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define NB 2

/* Refill the inputs, schedule the malloc failure, call zgetsqrhrt. */
#define LAPACKE_ZGETSQRHRT_ALLOC_TEST(layout_index, countdown, name, expected) \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, M, N, a, LD);                               \
        lapacke_test_zfill_rhs(layout, NB, N, t, LD);                          \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zgetsqrhrt)(layout, M, N, M, 2,  \
                                                          2, a, LD, t, LD),    \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zgetsqrhrt)
{
    lapack_complex_double a[LD * LD];
    lapack_complex_double t[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP("zgetsqrhrt a", l, M, N, a, LD,
                                lapacke_test_region_full, -7,
                                (lapacke_test_zfill(layout, M, N, a, LD),
                                 lapacke_test_zfill_rhs(layout, NB, N, t, LD)),
                                API_SUFFIX(LAPACKE_zgetsqrhrt)(
                                    layout, M, N, M, 2, 2, a, LD, t, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_rhs(layout, NB, N, t, LD);
        lapacke_test_zfill_nan(layout, M, N, a, LD);
        lapacke_test_check("zgetsqrhrt NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_zgetsqrhrt)(layout, M, N, M, 2, 2,
                                                          a, LD, t, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZGETSQRHRT_ALLOC_TEST(0, 0, "zgetsqrhrt work alloc failure (work)",
                                  LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGETSQRHRT_ALLOC_TEST(0, 1, "zgetsqrhrt allocation count", 0);
    lapacke_test_check_alloc_count("zgetsqrhrt col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZGETSQRHRT_ALLOC_TEST(1, 0, "zgetsqrhrt work alloc failure (work)",
                                  LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGETSQRHRT_ALLOC_TEST(1, 1,
                                  "zgetsqrhrt transpose alloc failure (a_t)",
                                  LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGETSQRHRT_ALLOC_TEST(1, 2,
                                  "zgetsqrhrt transpose alloc failure (t_t)",
                                  LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGETSQRHRT_ALLOC_TEST(1, 3, "zgetsqrhrt allocation count", 0);
    lapacke_test_check_alloc_count("zgetsqrhrt row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZGETSQRHRT_ALLOC_TEST(2, 0, "zgetsqrhrt invalid matrix_layout", -1);
    lapacke_test_check_alloc_count(
        "zgetsqrhrt invalid layout allocation count");
}
