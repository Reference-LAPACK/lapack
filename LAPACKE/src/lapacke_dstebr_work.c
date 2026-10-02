#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_dstebr_work)( lapack_int n, double* d, double* e,
                                double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork )
{
    lapack_int info = 0;
    /* Call LAPACK function and adjust info */
    LAPACK_dstebr( &n, d, e, work, &lwork, iwork, &liwork, &info );
    return info;
}
