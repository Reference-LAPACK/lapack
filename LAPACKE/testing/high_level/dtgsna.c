#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One call for all eigenpairs: refill the upper triangular pair (A,B), which
 * is in generalized Schur form, and the eigenvector matrices, schedule the
 * malloc failure countdown and check the info result of a dtgsna call with
 * the given job in the indexed layout. */
#define LAPACKE_DTGSNA_TEST(layout_index, job, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill_tri(layout, 'U', N, a, LD);                         \
        lapacke_test_dfill_tri(layout, 'U', N, b, LD);                         \
        lapacke_test_dfill(layout, N, N, vl, LD);                              \
        lapacke_test_dfill(layout, N, N, vr, LD);                              \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dtgsna)(                         \
                               layout, job, 'A', select, N, a, LD, b, LD, vl,  \
                               LD, vr, LD, s, dif, N, &m),                     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dtgsna)
{
    double a[LD * LD], b[LD * LD], vl[LD * LD], vr[LD * LD], s[N], dif[N];
    const lapack_logical select[N] = {1, 1, 1};
    lapack_int m;

    /* DTGSNA uses WORK for every job: with job = 'E' it holds the
     * products with the eigenvectors that form the eigenvalue condition
     * numbers. LAPACKE must allocate it for job = 'E' too. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_DTGSNA_TEST(l, 'E', -1, "dtgsna job = E", 0);
        LAPACKE_DTGSNA_TEST(l, 'V', -1, "dtgsna job = V", 0);
        LAPACKE_DTGSNA_TEST(l, 'B', -1, "dtgsna job = B", 0);
    }

    /* Column-major, job = 'E': the workspace is the only allocation (IWORK
     * is not referenced). */
    LAPACKE_DTGSNA_TEST(0, 'E', 0, "dtgsna work alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_DTGSNA_TEST(0, 'E', 1, "dtgsna allocation count", 0);
    lapacke_test_check_alloc_count("dtgsna col-major allocation count");

    /* Row-major, job = 'E': the workspace, then the transposed copies of A,
     * B, VL and VR. */
    LAPACKE_DTGSNA_TEST(1, 'E', 0, "dtgsna work alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DTGSNA_TEST(1, 'E', 1, "dtgsna transpose alloc failure (a)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DTGSNA_TEST(1, 'E', 2, "dtgsna transpose alloc failure (b)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DTGSNA_TEST(1, 'E', 3, "dtgsna transpose alloc failure (vl)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DTGSNA_TEST(1, 'E', 4, "dtgsna transpose alloc failure (vr)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_DTGSNA_TEST(1, 'E', 5, "dtgsna allocation count", 0);
    lapacke_test_check_alloc_count("dtgsna row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_DTGSNA_TEST(2, 'E', 0, "dtgsna invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dtgsna invalid layout allocation count");
}
