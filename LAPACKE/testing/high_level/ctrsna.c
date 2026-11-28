#include <math.h>
#include "lapacke_test.h"

LAPACKE_TEST(ctrsna)
{
    lapack_complex_float t[4] = {0}, unused[1] = {0};
    const lapack_logical select[2] = {1, 1};

    t[0] = lapack_make_complex_float(1, 0);
    t[3] = lapack_make_complex_float(3, 0);

    /* With JOB = 'V', neither eigenvector matrix is referenced. */
    for (size_t l = 0; l < 2; l++) {
        float s[2] = {0}, sep[2] = {0};
        lapack_int m = 0;
        lapack_int info = API_SUFFIX(LAPACKE_ctrsna)(
            lapacke_test_layouts[l], 'V', 'A', select, 2, t, 2,
            unused, 1, unused, 1, s, sep, 2, &m);

        lapacke_test_check("ctrsna", lapacke_test_layout_names[l], info, 0);
        if (info == 0) {
            lapacke_test_check("eigenpair count", lapacke_test_layout_names[l], m, 2);
            for (size_t i = 0; i < 2; i++) {
                lapacke_test_check("eigenvalue separation", lapacke_test_layout_names[l],
                                   fabs(sep[i] - 2) < 1e-5, 1);
            }
        }
    }

    /* Check the row-major stride limits with MM = N, including N = 0. */
    {
        static const struct {
            const char *label;
            char job;
            lapack_int n, ldvl, ldvr, expected;
        } cases[] = {
            {"both read", 'E', 2, 2, 2, 0},
            {"short VL", 'E', 2, 1, 2, -9},
            {"short VR", 'B', 2, 2, 1, -11},
            {"zero unused LDVL", 'V', 2, 0, 1, -9},
            {"zero unused LDVR", 'V', 2, 1, 0, -11},
            {"empty", 'V', 0, 0, 0, 0},
        };
        lapack_complex_float vec[4] = {0};
        float s[2], sep[2];
        lapack_int m;

        vec[0] = vec[3] = t[0];  /* Eigenvectors of the diagonal T. */
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
            lapack_int info = API_SUFFIX(LAPACKE_ctrsna)(
                LAPACK_ROW_MAJOR, cases[i].job, 'A', select, cases[i].n,
                t, 2, vec, cases[i].ldvl, vec, cases[i].ldvr,
                s, sep, cases[i].n, &m);
            lapacke_test_check("ctrsna", cases[i].label, info, cases[i].expected);
        }
    }
}
