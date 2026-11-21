#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One eigenvalue-only call: refill the band matrix with kd superdiagonals,
 * schedule the malloc failure countdown and check the info result of a
 * chbevx_2stage call with the given range in the indexed layout. */
#define LAPACKE_CHBEVX_2STAGE_TEST(layout_index, kd, range, countdown, name,   \
                                   expected)                                   \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, (kd) + 1, N, ab, LD);                       \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_chbevx_2stage)(                  \
                               layout, 'N', range, 'U', N, kd, ab, LD, q, LD,  \
                               -10.0, 10.0, 1, N - 1, 0.0, &m, w, z, LD,       \
                               ifail),                                         \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(csbevx_2stage)
{
    lapack_complex_float ab[LD * LD], q[LD * LD], z[LD * LD];
    float w[N];
    lapack_int m, ifail[N];

    /* CHBEVX_2STAGE is the complex counterpart of SSBEVX_2STAGE (a real
     * symmetric band matrix becomes a Hermitian one). LAPACKE allocates the
     * WORK size of its workspace query and the fixed RWORK(7*N) and
     * IWORK(5*N). The test allocator fails the test when a call writes past
     * the end of a LAPACKE allocation. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_CHBEVX_2STAGE_TEST(l, 2, 'A', -1, "chbevx_2stage range = A", 0);
        LAPACKE_CHBEVX_2STAGE_TEST(l, 2, 'I', -1, "chbevx_2stage range = I", 0);
        LAPACKE_CHBEVX_2STAGE_TEST(l, 2, 'V', -1, "chbevx_2stage range = V", 0);
    }

    /* Column-major: IWORK, RWORK and the workspace are the only
     * allocations. */
    LAPACKE_CHBEVX_2STAGE_TEST(0, 2, 'A', 0,
                               "chbevx_2stage iwork alloc failure",
                               LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHBEVX_2STAGE_TEST(0, 2, 'A', 1,
                               "chbevx_2stage rwork alloc failure",
                               LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHBEVX_2STAGE_TEST(0, 2, 'A', 2, "chbevx_2stage work alloc failure",
                               LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_CHBEVX_2STAGE_TEST(0, 2, 'A', 3, "chbevx_2stage allocation count",
                               0);
    lapacke_test_check_alloc_count("chbevx_2stage col-major allocation count");

    /* Row-major allocates the same workspace, then the transposed copy of AB
     * (Q and Z are not referenced without eigenvectors). The workspace query
     * must not allocate. */
    LAPACKE_CHBEVX_2STAGE_TEST(1, 2, 'A', 0,
                               "chbevx_2stage iwork alloc failure",
                               LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHBEVX_2STAGE_TEST(1, 2, 'A', 1,
                               "chbevx_2stage rwork alloc failure",
                               LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHBEVX_2STAGE_TEST(1, 2, 'A', 2, "chbevx_2stage work alloc failure",
                               LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_CHBEVX_2STAGE_TEST(1, 2, 'A', 3,
                               "chbevx_2stage transpose alloc failure (ab)",
                               LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_CHBEVX_2STAGE_TEST(1, 2, 'A', 4, "chbevx_2stage allocation count",
                               0);
    lapacke_test_check_alloc_count("chbevx_2stage row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_CHBEVX_2STAGE_TEST(2, 2, 'A', 0,
                               "chbevx_2stage invalid matrix_layout", -1);
    lapacke_test_check_alloc_count(
        "chbevx_2stage invalid layout allocation count");
}
