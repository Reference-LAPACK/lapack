#include <lapack.h>
#ifdef CHECK_INTERFACES
#include <cblas.h>
#include <lapacke.h>
#ifdef CHECK_EXTENDED
#include <cblas_64.h>
#include <lapacke_64.h>
#endif
#endif

typedef char lapack_integer_abi[sizeof(lapack_int) == CHECK_INTEGER_BYTES ? 1 : -1];
#ifdef CHECK_INTERFACES
typedef char cblas_integer_abi[sizeof(CBLAS_INT) == CHECK_INTEGER_BYTES ? 1 : -1];
#endif

int main(void)
{
    double a = 2.0, b = 6.0;
#ifndef CHECK_ONLY_EXTENDED
    lapack_int n = 1, nrhs = 1, ipiv = 0;
#ifndef CHECK_INTERFACES
    lapack_int info = -1;
    LAPACK_dgesv(&n, &nrhs, &a, &n, &ipiv, &b, &n, &info);
    if (info != 0 || ipiv != 1 || b != 3.0) return 1;
#else
    if (LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, &a, n, &ipiv, &b, n)
        || ipiv != 1 || b != 3.0) return 2;
    if (cblas_ddot(1, &a, 1, &b, 1) != 6.0) return 3;
#endif
#endif
#ifdef CHECK_EXTENDED
    int64_t ipiv64 = 0;
    a = 2.0;
    b = 6.0;
    if (LAPACKE_dgesv_64(LAPACK_COL_MAJOR, 1, 1, &a, 1, &ipiv64, &b, 1)
        || ipiv64 != 1 || b != 3.0) return 4;
    if (cblas_ddot_64(1, &a, 1, &b, 1) != 6.0) return 5;
#endif
    return 0;
}
