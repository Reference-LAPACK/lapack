#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One call computing only the singular values: refill the rows-by-cols A,
 * schedule the malloc failure countdown and check the info result of a
 * dgejsv call in the indexed layout. */
#define LAPACKE_DGEJSV_TEST(layout_index, rows, cols, countdown, name,         \
                            expected)                                          \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, rows, cols, a, LD);                         \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dgejsv)(                         \
                               layout, 'C', 'N', 'N', 'N', 'N', 'N', rows,     \
                               cols, a, LD, sva, u, LD, v, LD, stat, istat),   \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dgejsv)
{
    double a[LD * LD], u[LD * LD], v[LD * LD], sva[LD], stat[7];
    lapack_int istat[3];

    /* DGEJSV sets IWORK(1:3) and needs MAX(3,M+3*N) entries, which LAPACKE
     * allocates. The test allocator fails the test when a call writes past
     * the end of a LAPACKE allocation. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_DGEJSV_TEST(l, 1, 1, -1, "dgejsv 1-by-1", 0);
        LAPACKE_DGEJSV_TEST(l, M, N, -1, "dgejsv M-by-N", 0);
    }

    /* Column-major: IWORK and the workspace are the only allocations. */
    LAPACKE_DGEJSV_TEST(0, M, N, 0, "dgejsv iwork alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGEJSV_TEST(0, M, N, 1, "dgejsv work alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_DGEJSV_TEST(0, M, N, 2, "dgejsv allocation count", 0);
    lapacke_test_check_alloc_count("dgejsv col-major allocation count");

    /* Row-major allocates the same workspace, then the transposed copy of A
     * (U and V are not referenced). */
    LAPACKE_DGEJSV_TEST(1, M, N, 2, "dgejsv transpose alloc failure (a)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_DGEJSV_TEST(1, M, N, 3, "dgejsv allocation count", 0);
    lapacke_test_check_alloc_count("dgejsv row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_DGEJSV_TEST(2, M, N, 0, "dgejsv invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgejsv invalid layout allocation count");
}
