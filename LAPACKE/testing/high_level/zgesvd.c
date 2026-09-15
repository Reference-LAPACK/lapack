#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call zgesvd. */
#define LAPACKE_ZGESVD_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, M, N, a, LD);                               \
        lapacke_test_dfill_pos(LD * LD, s);                                    \
        lapacke_test_zfill_rhs(layout, M, M, u, LD);                           \
        lapacke_test_zfill_rhs(layout, N, N, vt, LD);                          \
        lapacke_test_dfill_vec(LD * LD, superb);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zgesvd)(layout, 'A', 'A', M, N,  \
                                                      a, LD, s, u, LD, vt, LD, \
                                                      superb),                 \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zgesvd)
{
    lapack_complex_double a[LD * LD];
    double s[LD * LD];
    lapack_complex_double u[LD * LD];
    lapack_complex_double vt[LD * LD];
    double superb[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_ZNAN_SWEEP(
            "zgesvd a", l, M, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_zfill(layout, M, N, a, LD),
             lapacke_test_dfill_pos(LD * LD, s),
             lapacke_test_zfill_rhs(layout, M, M, u, LD),
             lapacke_test_zfill_rhs(layout, N, N, vt, LD),
             lapacke_test_dfill_vec(LD * LD, superb)),
            API_SUFFIX(LAPACKE_zgesvd)(layout, 'A', 'A', M, N, a, LD, s, u, LD,
                                       vt, LD, superb));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_pos(LD * LD, s);
        lapacke_test_zfill_rhs(layout, M, M, u, LD);
        lapacke_test_zfill_rhs(layout, N, N, vt, LD);
        lapacke_test_dfill_vec(LD * LD, superb);
        lapacke_test_zfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "zgesvd NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgesvd)(layout, 'A', 'A', M, N, a, LD, s, u, LD,
                                       vt, LD, superb) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_ZGESVD_ALLOC_TEST(0, 0, "zgesvd work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGESVD_ALLOC_TEST(0, 1, "zgesvd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGESVD_ALLOC_TEST(0, 2, "zgesvd allocation count", 0);
    lapacke_test_check_alloc_count("zgesvd col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_ZGESVD_ALLOC_TEST(1, 0, "zgesvd work alloc failure (rwork)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGESVD_ALLOC_TEST(1, 1, "zgesvd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGESVD_ALLOC_TEST(1, 2, "zgesvd transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGESVD_ALLOC_TEST(1, 3, "zgesvd transpose alloc failure (u_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGESVD_ALLOC_TEST(1, 4, "zgesvd transpose alloc failure (vt_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZGESVD_ALLOC_TEST(1, 5, "zgesvd allocation count", 0);
    lapacke_test_check_alloc_count("zgesvd row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_ZGESVD_ALLOC_TEST(2, 0, "zgesvd invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgesvd invalid layout allocation count");
}
