#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_sstebr)( lapack_int n, float* d, float* e )
{
    lapack_int info = 0;
    lapack_int liwork = -1;
    lapack_int lwork = -1;
    lapack_int* iwork = NULL;
    float* work = NULL;
    lapack_int iwork_query;
    float work_query;
#ifndef LAPACK_DISABLE_NAN_CHECK
    if( LAPACKE_get_nancheck() ) {
        /* Optionally check input matrices for NaNs */
        if( API_SUFFIX(LAPACKE_s_nancheck)( n, d, 1 ) ) {
            return -2;
        }
        if( API_SUFFIX(LAPACKE_s_nancheck)( n-1, e, 1 ) ) {
            return -3;
        }
    }
#endif
    /* Query optimal working array(s) size */
    info = API_SUFFIX(LAPACKE_sstebr_work)( n, d, e, &work_query, lwork,
                                &iwork_query, liwork );
    if( info != 0 ) {
        goto exit_level_0;
    }
    liwork = iwork_query;
    lwork = (lapack_int)work_query;
    /* Allocate memory for work arrays */
    iwork = (lapack_int*)LAPACKE_malloc( sizeof(lapack_int) * liwork );
    if( iwork == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_0;
    }
    work = (float*)LAPACKE_malloc( sizeof(float) * lwork );
    if( work == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_1;
    }
    /* Call middle-level interface */
    info = API_SUFFIX(LAPACKE_sstebr_work)( n, d, e, work, lwork, iwork,
                                liwork );
    /* Release memory and exit */
    LAPACKE_free( work );
exit_level_1:
    LAPACKE_free( iwork );
exit_level_0:
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_sstebr", info );
    }
    return info;
}
