#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* Refill the inputs, schedule the malloc failure, call sgesvd. */
#define LAPACKE_SGESVD_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, M, N, a, LD);                               \
        lapacke_test_sfill_pos(LD * LD, s);                                    \
        lapacke_test_sfill_rhs(layout, M, M, u, LD);                           \
        lapacke_test_sfill_rhs(layout, N, N, vt, LD);                          \
        lapacke_test_sfill_vec(LD * LD, superb);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sgesvd)(layout, 'A', 'A', M, N,  \
                                                      a, LD, s, u, LD, vt, LD, \
                                                      superb),                 \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sgesvd)
{
    float a[LD * LD];
    float s[LD * LD];
    float u[LD * LD];
    float vt[LD * LD];
    float superb[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_SNAN_SWEEP(
            "sgesvd a", l, M, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_sfill(layout, M, N, a, LD),
             lapacke_test_sfill_pos(LD * LD, s),
             lapacke_test_sfill_rhs(layout, M, M, u, LD),
             lapacke_test_sfill_rhs(layout, N, N, vt, LD),
             lapacke_test_sfill_vec(LD * LD, superb)),
            API_SUFFIX(LAPACKE_sgesvd)(layout, 'A', 'A', M, N, a, LD, s, u, LD,
                                       vt, LD, superb));

        /* NaN checks off: all-NaN input must reach the Fortran routine. */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_pos(LD * LD, s);
        lapacke_test_sfill_rhs(layout, M, M, u, LD);
        lapacke_test_sfill_rhs(layout, N, N, vt, LD);
        lapacke_test_sfill_vec(LD * LD, superb);
        lapacke_test_sfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "sgesvd NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_sgesvd)(layout, 'A', 'A', M, N, a, LD, s, u, LD,
                                       vt, LD, superb) < 0,
            0);
        LAPACKE_set_nancheck(1);
    }

    /* column-major: the high-level workspaces */
    LAPACKE_SGESVD_ALLOC_TEST(0, 0, "sgesvd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGESVD_ALLOC_TEST(0, 1, "sgesvd allocation count", 0);
    lapacke_test_check_alloc_count("sgesvd col-major allocation count");

    /* row-major: the workspaces, then the transposed copies */
    LAPACKE_SGESVD_ALLOC_TEST(1, 0, "sgesvd work alloc failure (work)",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_SGESVD_ALLOC_TEST(1, 1, "sgesvd transpose alloc failure (a_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGESVD_ALLOC_TEST(1, 2, "sgesvd transpose alloc failure (u_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGESVD_ALLOC_TEST(1, 3, "sgesvd transpose alloc failure (vt_t)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGESVD_ALLOC_TEST(1, 4, "sgesvd allocation count", 0);
    lapacke_test_check_alloc_count("sgesvd row-major allocation count");

    /* invalid matrix_layout: rejected before any allocation */
    LAPACKE_SGESVD_ALLOC_TEST(2, 0, "sgesvd invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgesvd invalid layout allocation count");
}
