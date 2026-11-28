#include "lapacke_test.h"

LAPACKE_TEST(ztgevc)
{
    lapack_complex_double s[4] = {0};
    lapack_complex_double p[4] = {0};
    lapack_logical select[2] = {1, 1};

    s[0] = lapack_make_complex_double(1, 0);
    s[3] = lapack_make_complex_double(6, 0);
    p[0] = lapack_make_complex_double(1, 0);
    p[3] = lapack_make_complex_double(2, 0);

    /* A diagonal pair with distinct eigenvalues has coordinate eigenvectors. */
    for (size_t l = 0; l < 2; l++) {
        for (int left = 0; left <= 1; left++) {
            const char *side = left ? "left" : "right";
            lapack_complex_double v[4] = {0}, unused[1] = {0};
            lapack_int m = 0;
            lapack_int info = API_SUFFIX(LAPACKE_ztgevc)(
                lapacke_test_layouts[l], left ? 'L' : 'R', 'A', select, 2,
                s, 2, p, 2, left ? v : unused, left ? 2 : 1,
                left ? unused : v, left ? 1 : 2, 2, &m);

            lapacke_test_check(side, lapacke_test_layout_names[l], info, 0);
            if (info == 0) {
                lapacke_test_check("eigenvector count", side, m, 2);
                for (size_t i = 0; i < 4; i++) {
                    const int nonzero = lapack_complex_double_real(v[i]) != 0 ||
                                    lapack_complex_double_imag(v[i]) != 0;
                    lapacke_test_check("eigenvector entry", side, nonzero,
                                       i == 0 || i == 3);
                }
            }
        }
    }

    /* With MM = 2, the requested side needs two columns and the unused side needs one. */
    {
        static const struct {
            const char *label;
            char side;
            lapack_int ldvl, ldvr, expected;
        } cases[] = {
            {"short requested VL", 'L', 1, 1, -11},
            {"short requested VR", 'R', 1, 1, -13},
            {"zero unused LDVL", 'R', 0, 2, -11},
            {"zero unused LDVR", 'L', 2, 0, -13},
        };
        lapack_complex_double vl[4] = {0}, vr[4] = {0};
        lapack_int m = 0;
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
            lapack_int info = API_SUFFIX(LAPACKE_ztgevc)(
                LAPACK_ROW_MAJOR, cases[i].side, 'A', select, 2, s, 2, p, 2,
                vl, cases[i].ldvl, vr, cases[i].ldvr, 2, &m);
            lapacke_test_check("ztgevc", cases[i].label, info, cases[i].expected);
        }
    }
}
