#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD
#define NB 2

/* A real diagonal matrix with ipiv(k) = k is its own factorization by
 * CSYTRF and by CHETRF (all pivots 1-by-1, U = L = I), so it is valid input
 * for both CSYTRI2X and CHETRI2X. */
static void lapacke_test_csytri2x_fill(lapack_complex_float *a,
                                       lapack_int *ipiv)
{
    for (lapack_int i = 0; i < LD * LD; i++)
        a[i] = lapack_make_complex_float(0.0, 0.0);
    for (lapack_int i = 0; i < N; i++) {
        a[i * LD + i] = lapack_make_complex_float(2.0 + i, 0.0);
        ipiv[i] = i + 1;
    }
}

/* One call: refill the factored matrix, schedule the malloc failure
 * countdown and check the info result of a call of routine (csytri2x or
 * chetri2x) in the indexed layout. */
#define LAPACKE_CSYTRI2X_TEST(routine, layout_index, uplo, countdown, name,    \
                              expected)                                        \
    do {                                                                       \
        lapacke_test_csytri2x_fill(a, ipiv);                                   \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_##routine)(                      \
                               lapacke_test_layouts[layout_index], uplo, N, a, \
                               LD, ipiv, NB),                                  \
                           expected);                                          \
    } while (0)

/* The allocation checks for one routine. Column-major allocates only the
 * workspace; row-major allocates the workspace, then the transposed copy of
 * A. A countdown one past the last allocation must not fire, so the call
 * must succeed and consume the countdown exactly. An invalid matrix_layout
 * must be rejected as an error in argument 1, before any allocation. */
#define LAPACKE_CSYTRI2X_ALLOC_TESTS(routine)                                  \
    do {                                                                       \
        LAPACKE_CSYTRI2X_TEST(routine, 0, 'U', 0,                              \
                              #routine " work alloc failure",                  \
                              LAPACK_WORK_MEMORY_ERROR);                       \
        LAPACKE_CSYTRI2X_TEST(routine, 0, 'U', 1,                              \
                              #routine " allocation count", 0);                \
        lapacke_test_check_alloc_count(#routine                                \
                                       " col-major allocation count");         \
        LAPACKE_CSYTRI2X_TEST(routine, 1, 'U', 0,                              \
                              #routine " work alloc failure",                  \
                              LAPACK_WORK_MEMORY_ERROR);                       \
        LAPACKE_CSYTRI2X_TEST(routine, 1, 'U', 1,                              \
                              #routine " transpose alloc failure (a)",         \
                              LAPACK_TRANSPOSE_MEMORY_ERROR);                  \
        LAPACKE_CSYTRI2X_TEST(routine, 1, 'U', 2,                              \
                              #routine " allocation count", 0);                \
        lapacke_test_check_alloc_count(#routine                                \
                                       " row-major allocation count");         \
        LAPACKE_CSYTRI2X_TEST(routine, 2, 'U', 0,                              \
                              #routine " invalid matrix_layout", -1);          \
        lapacke_test_check_alloc_count(#routine                                \
                                       " invalid layout allocation count");    \
    } while (0)

LAPACKE_TEST(csytri2x)
{
    lapack_complex_float a[LD * LD];
    lapack_int ipiv[N];

    /* CSYTRI2X and CHETRI2X take WORK(N+NB+1,NB+3), which LAPACKE
     * allocates; they keep the inverted diagonal blocks in columns NB+2 and
     * NB+3. The test allocator fails the test when a call writes past the
     * end of a LAPACKE allocation. */
    for (size_t l = 0; l < 2; l++) {
        LAPACKE_CSYTRI2X_TEST(csytri2x, l, 'U', -1, "csytri2x uplo = U", 0);
        LAPACKE_CSYTRI2X_TEST(csytri2x, l, 'L', -1, "csytri2x uplo = L", 0);
        LAPACKE_CSYTRI2X_TEST(chetri2x, l, 'U', -1, "chetri2x uplo = U", 0);
        LAPACKE_CSYTRI2X_TEST(chetri2x, l, 'L', -1, "chetri2x uplo = L", 0);
    }

    LAPACKE_CSYTRI2X_ALLOC_TESTS(csytri2x);
    LAPACKE_CSYTRI2X_ALLOC_TESTS(chetri2x);
}
