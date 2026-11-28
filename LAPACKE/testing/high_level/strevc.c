#include "lapacke_test.h"

LAPACKE_TEST(strevc)
{
    float t[4] = {0};
    lapack_logical select[2] = {1, 1};

    t[0] = 1;
    t[3] = 3;

    /* A diagonal matrix with distinct eigenvalues has coordinate eigenvectors. */
    for (size_t l = 0; l < 2; l++) {
        for (int left = 0; left <= 1; left++) {
            const char *side = left ? "left" : "right";
            float v[4] = {0}, unused[1] = {0};
            lapack_int m = 0;
            lapack_int info = API_SUFFIX(LAPACKE_strevc)(
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

    /* Check the row-major stride limits with MM = N, including N = 0. */
    {
        static const struct {
            const char *label;
            char side;
            lapack_int n, ldvl, ldvr, expected;
        } cases[] = {
            {"both sides", 'B', 2, 2, 2, 0},
            {"short VL", 'L', 2, 1, 1, -9},
            {"short VR", 'R', 2, 1, 1, -11},
            {"zero unused LDVL", 'R', 2, 0, 2, -9},
            {"zero unused LDVR", 'L', 2, 2, 0, -11},
            {"empty left", 'L', 0, 0, 0, 0},
            {"empty right", 'R', 0, 0, 0, 0},
        };
        float vl[4] = {0}, vr[4] = {0};
        lapack_int m;

        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
            lapack_int info = API_SUFFIX(LAPACKE_strevc)(
                LAPACK_ROW_MAJOR, cases[i].side, 'A', select, cases[i].n,
                t, 2, vl, cases[i].ldvl, vr, cases[i].ldvr, cases[i].n, &m);
            lapacke_test_check("strevc", cases[i].label, info, cases[i].expected);
        }
    }
}
