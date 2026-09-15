#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call dgesvd. */
#define LAPACKE_DGESVD_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, M, N, a, LD);                               \
        lapacke_test_dfill_pos(LD * LD, s);                                    \
        lapacke_test_dfill_rhs(layout, M, M, u, LD);                           \
        lapacke_test_dfill_rhs(layout, N, N, vt, LD);                          \
        lapacke_test_dfill_vec(LD * LD, superb);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dgesvd)(layout, 'A', 'A', M, N,  \
                                                      a, LD, s, u, LD, vt, LD, \
                                                      superb),                 \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dgesvd)
{
    double a[LD * LD];
    double s[LD * LD];
    double u[LD * LD];
    double vt[LD * LD];
    double superb[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgesvd a", l, M, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_dfill(layout, M, N, a, LD),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_dfill_rhs(layout, M, M, u, LD),
             lapacke_test_dfill_rhs(layout, N, N, vt, LD),
             lapacke_test_dfill_vec(LD * LD, superb)),
            API_SUFFIX(LAPACKE_dgesvd)(layout, 'A', 'A', M, N, a, LD, s, u, LD,
                                       vt, LD, superb));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_pos(LD * LD, s);
        lapacke_test_dfill_rhs(layout, M, M, u, LD);
        lapacke_test_dfill_rhs(layout, N, N, vt, LD);
        lapacke_test_dfill_vec(LD * LD, superb);
        lapacke_test_dfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "dgesvd NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_dgesvd)(layout, 'A', 'A', M, N, a, LD, s, u, LD,
                                       vt, LD, superb) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_DGESVD_ALLOC_TEST(0, 0, "dgesvd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGESVD_ALLOC_TEST(0, 1, "dgesvd allocation count", 0);
    lapacke_test_check_alloc_count("dgesvd col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_DGESVD_ALLOC_TEST(1, 0, "dgesvd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGESVD_ALLOC_TEST(1, 1, "dgesvd transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGESVD_ALLOC_TEST(1, 2, "dgesvd transpose alloc failure (u_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGESVD_ALLOC_TEST(1, 3, "dgesvd transpose alloc failure (vt_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGESVD_ALLOC_TEST(1, 4, "dgesvd allocation count", 0);
    lapacke_test_check_alloc_count("dgesvd row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_DGESVD_ALLOC_TEST(2, 0, "dgesvd invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgesvd invalid layout allocation count");
}
