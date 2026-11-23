#include "lapacke_test.h"

LAPACKE_TEST(dtrevc)
{
    double t[4] = {0};
    lapack_logical select[2] = {1, 1};

    t[0] = 1;
    t[3] = 3;

    /* A diagonal matrix with distinct eigenvalues has coordinate eigenvectors. */
    for (size_t l = 0; l < 2; l++) {
        for (int left = 0; left <= 1; left++) {
            const char *side = left ? "left" : "right";
            double v[4] = {0}, unused[1] = {0};
            lapack_int m = 0;
            lapack_int info = API_SUFFIX(LAPACKE_dtrevc)(
                lapacke_test_layouts[l], left ? 'L' : 'R', 'A', select, 2,
                t, 2, left ? v : unused, left ? 2 : 1,
                left ? unused : v, left ? 1 : 2, 2, &m);

            lapacke_test_check(side, lapacke_test_layout_names[l], info, 0);
            if (info == 0) {
                lapacke_test_check("eigenvector count", side, m, 2);
                for (size_t i = 0; i < 4; i++) {
                    const int nonzero = v[i] != 0;
                    lapacke_test_check("eigenvector entry", side, nonzero,
                                       i == 0 || i == 3);
                }
            }
        }
    }
}
