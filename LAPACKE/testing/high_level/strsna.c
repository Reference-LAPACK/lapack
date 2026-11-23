#include <math.h>
#include "lapacke_test.h"

LAPACKE_TEST(strsna)
{
    float t[4] = {0}, unused[1] = {0};
    const lapack_logical select[2] = {1, 1};

    t[0] = 1;
    t[3] = 3;

    /* With JOB = 'V', neither eigenvector matrix is referenced. */
    for (size_t l = 0; l < 2; l++) {
        float s[2] = {0}, sep[2] = {0};
        lapack_int m = 0;
        lapack_int info = API_SUFFIX(LAPACKE_strsna)(
            lapacke_test_layouts[l], 'V', 'A', select, 2, t, 2,
            unused, 1, unused, 1, s, sep, 2, &m);

        lapacke_test_check("strsna", lapacke_test_layout_names[l], info, 0);
        if (info == 0) {
            lapacke_test_check("eigenpair count", lapacke_test_layout_names[l], m, 2);
            for (size_t i = 0; i < 2; i++) {
                lapacke_test_check("eigenvalue separation", lapacke_test_layout_names[l],
                                   fabs(sep[i] - 2) < 1e-5, 1);
            }
        }
    }
}
