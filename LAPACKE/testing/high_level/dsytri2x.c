#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define NB 2

/* A diagonal matrix with ipiv(k) = k is its own factorization by DSYTRF
 * (all pivots 1-by-1, U = L = I), so it is valid input for DSYTRI2X. */
static void lapacke_test_dsytri2x_fill(double *a, lapack_int *ipiv)
{
    for (lapack_int i = 0; i < LD * LD; i++)
        a[i] = 0.0;
    for (lapack_int i = 0; i < N; i++) {
        a[i * LD + i] = 2.0 + i;
        ipiv[i] = i + 1;
    }
}

/* One call: refill the factored matrix, schedule the malloc failure
 * countdown and check the info result of a dsytri2x call in the indexed
 * layout. */
#define LAPACKE_DSYTRI2X_TEST(layout_index, uplo, countdown, name, expected)   \
    do {                                                                       \
        lapacke_test_dsytri2x_fill(a, ipiv);                                   \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_dsytri2x)(                       \
                               lapacke_test_layouts[layout_index], uplo, N, a, \
                               LD, ipiv, NB),                                  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(dsytri2x)
{
    double a[LD * LD];
    lapack_int ipiv[N];

    /* DSYTRI2X takes WORK(N+NB+1,NB+3), which LAPACKE allocates; it keeps
     * the inverted diagonal blocks in columns NB+2 and NB+3. The test
     * allocator fails the test when a call writes past the end of a LAPACKE
     * allocation. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_DSYTRI2X_TEST(l, 'U', -1, "dsytri2x uplo = U", 0);
        LAPACKE_DSYTRI2X_TEST(l, 'L', -1, "dsytri2x uplo = L", 0);
    }

    /* Column-major: the workspace is the only allocation. */
    LAPACKE_DSYTRI2X_TEST(0, 'U', 0, "dsytri2x work alloc failure",
                          LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_DSYTRI2X_TEST(0, 'U', 1, "dsytri2x allocation count", 0);
    lapacke_test_check_alloc_count("dsytri2x col-major allocation count");

    /* Row-major allocates the workspace, then the transposed copy of A. */
    LAPACKE_DSYTRI2X_TEST(1, 'U', 0, "dsytri2x work alloc failure",
                          LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DSYTRI2X_TEST(1, 'U', 1, "dsytri2x transpose alloc failure (a)",
                          LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_DSYTRI2X_TEST(1, 'U', 2, "dsytri2x allocation count", 0);
    lapacke_test_check_alloc_count("dsytri2x row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_DSYTRI2X_TEST(2, 'U', 0, "dsytri2x invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dsytri2x invalid layout allocation count");
}
