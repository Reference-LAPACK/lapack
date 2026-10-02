#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_sstebr_work)( lapack_int n, float* d, float* e,
                                float* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork )
{
    lapack_int info = 0;
    /* Call LAPACK function and adjust info */
    LAPACK_sstebr( &n, d, e, work, &lwork, iwork, &liwork, &info );
    return info;
}
