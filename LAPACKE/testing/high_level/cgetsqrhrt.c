#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define NB 2

/* Refill the inputs, schedule the malloc failure, call cgetsqrhrt. */
#define LAPACKE_CGETSQRHRT_ALLOC_TEST(layout_index, countdown, name, expected) \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, M, N, a, LD);                               \
        lapacke_test_cfill_rhs(layout, NB, N, t, LD);                          \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cgetsqrhrt)(layout, M, N, M, 2,  \
                                                          2, a, LD, t, LD),    \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cgetsqrhrt)
{
    lapack_complex_float a[LD * LD];
    lapack_complex_float t[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP("cgetsqrhrt a", l, M, N, a, LD,
                                lapacke_test_region_full, -7,
                                (lapacke_test_cfill(layout, M, N, a, LD),
                                 lapacke_test_cfill_rhs(layout, NB, N, t, LD)),
                                API_SUFFIX(LAPACKE_cgetsqrhrt)(
                                    layout, M, N, M, 2, 2, a, LD, t, LD));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_rhs(layout, NB, N, t, LD);
        lapacke_test_cfill_nan(layout, M, N, a, LD);
        lapacke_test_check("cgetsqrhrt NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_cgetsqrhrt)(layout, M, N, M, 2, 2,
                                                          a, LD, t, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_CGETSQRHRT_ALLOC_TEST(0, 0, "cgetsqrhrt work alloc failure (work)",
                                  LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGETSQRHRT_ALLOC_TEST(0, 1, "cgetsqrhrt allocation count", 0);
    lapacke_test_check_alloc_count("cgetsqrhrt col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_CGETSQRHRT_ALLOC_TEST(1, 0, "cgetsqrhrt work alloc failure (work)",
                                  LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CGETSQRHRT_ALLOC_TEST(1, 1,
                                  "cgetsqrhrt transpose alloc failure (a_t)",
                                  LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGETSQRHRT_ALLOC_TEST(1, 2,
                                  "cgetsqrhrt transpose alloc failure (t_t)",
                                  LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_CGETSQRHRT_ALLOC_TEST(1, 3, "cgetsqrhrt allocation count", 0);
    lapacke_test_check_alloc_count("cgetsqrhrt row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_CGETSQRHRT_ALLOC_TEST(2, 0, "cgetsqrhrt invalid matrix_layout", -1);
    lapacke_test_check_alloc_count(
        "cgetsqrhrt invalid layout allocation count");
}
