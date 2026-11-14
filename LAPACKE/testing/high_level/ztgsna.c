#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One call for all eigenpairs: refill the upper triangular pair (A,B), which
 * is in generalized Schur form, and the eigenvector matrices, schedule the
 * malloc failure countdown and check the info result of a ztgsna call with
 * the given job in the indexed layout. */
#define LAPACKE_ZTGSNA_TEST(layout_index, job, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill_tri(layout, 'U', N, a, LD);                         \
        lapacke_test_zfill_tri(layout, 'U', N, b, LD);                         \
        lapacke_test_zfill(layout, N, N, vl, LD);                              \
        lapacke_test_zfill(layout, N, N, vr, LD);                              \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_ztgsna)(                         \
                               layout, job, 'A', select, N, a, LD, b, LD, vl,  \
                               LD, vr, LD, s, dif, N, &m),                     \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(ztgsna)
{
    lapack_complex_double a[LD * LD], b[LD * LD], vl[LD * LD], vr[LD * LD];
    double s[N], dif[N];
    const lapack_logical select[N] = {1, 1, 1};
    lapack_int m;

    /* ZTGSNA uses WORK for every job: with job = 'E' it holds the
     * products with the eigenvectors that form the eigenvalue condition
     * numbers. LAPACKE must allocate it for job = 'E' too. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_ZTGSNA_TEST(l, 'E', -1, "ztgsna job = E", 0);
        LAPACKE_ZTGSNA_TEST(l, 'V', -1, "ztgsna job = V", 0);
        LAPACKE_ZTGSNA_TEST(l, 'B', -1, "ztgsna job = B", 0);
    }

    /* Column-major, job = 'E': the workspace is the only allocation (IWORK
     * is not referenced). */
    LAPACKE_ZTGSNA_TEST(0, 'E', 0, "ztgsna work alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_ZTGSNA_TEST(0, 'E', 1, "ztgsna allocation count", 0);
    lapacke_test_check_alloc_count("ztgsna col-major allocation count");

    /* Row-major, job = 'E': the workspace, then the transposed copies of A,
     * B, VL and VR. */
    LAPACKE_ZTGSNA_TEST(1, 'E', 0, "ztgsna work alloc failure",
                        LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZTGSNA_TEST(1, 'E', 1, "ztgsna transpose alloc failure (a)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZTGSNA_TEST(1, 'E', 2, "ztgsna transpose alloc failure (b)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZTGSNA_TEST(1, 'E', 3, "ztgsna transpose alloc failure (vl)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_ZTGSNA_TEST(1, 'E', 4, "ztgsna transpose alloc failure (vr)",
                        LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_ZTGSNA_TEST(1, 'E', 5, "ztgsna allocation count", 0);
    lapacke_test_check_alloc_count("ztgsna row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_ZTGSNA_TEST(2, 'E', 0, "ztgsna invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("ztgsna invalid layout allocation count");
}
