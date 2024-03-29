#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_strsyl3_work)( int matrix_layout, char trana, char tranb,
                                 lapack_int isgn, lapack_int m, lapack_int n,
                                 const float* a, lapack_int lda,
                                 const float* b, lapack_int ldb, float* c,
                                 lapack_int ldc, float* scale,
                                 lapack_int* iwork, lapack_int liwork,
                                 float* swork, lapack_int ldswork )
{
    lapack_int info = 0;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        LAPACK_strsyl3( &trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc,
                        scale, iwork, &liwork, swork, &ldswork, &info );
        if( info < 0 ) {
            info = info - 1;
        }
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        lapack_int lda_t = MAX(1,m);
        lapack_int ldb_t = MAX(1,n);
        lapack_int ldc_t = MAX(1,m);
        float* a_t = NULL;
        float* b_t = NULL;
        float* c_t = NULL;
        /* Check leading dimension(s) */
        if( lda < m ) {
            info = -8;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_strsyl3_work", info );
            return info;
        }
        if( ldb < n ) {
            info = -10;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_strsyl3_work", info );
            return info;
        }
        if( ldc < n ) {
            info = -12;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_strsyl3_work", info );
            return info;
        }
        /* Allocate memory for temporary array(s) */
        a_t = (float*)LAPACKE_malloc( sizeof(float) * lda_t * MAX(1,m) );
        if( a_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_0;
        }
        b_t = (float*)LAPACKE_malloc( sizeof(float) * ldb_t * MAX(1,n) );
        if( b_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_1;
        }
        c_t = (float*)LAPACKE_malloc( sizeof(float) * ldc_t * MAX(1,n) );
        if( c_t == NULL ) {
            info = LAPACK_TRANSPOSE_MEMORY_ERROR;
            goto exit_level_2;
        }
        /* Transpose input matrices */
        API_SUFFIX(LAPACKE_sge_trans)( matrix_layout, m, m, a, lda, a_t, lda_t );
        API_SUFFIX(LAPACKE_sge_trans)( matrix_layout, n, n, b, ldb, b_t, ldb_t );
        API_SUFFIX(LAPACKE_sge_trans)( matrix_layout, m, n, c, ldc, c_t, ldc_t );
        /* Call LAPACK function and adjust info */
        LAPACK_strsyl3( &trana, &tranb, &isgn, &m, &n, a_t, &lda_t, b_t, &ldb_t,
                        c_t, &ldc_t, scale, iwork, &liwork, swork, &ldswork,
                        &info );
        if( info < 0 ) {
            info = info - 1;
        }
        /* Transpose output matrices */
        API_SUFFIX(LAPACKE_sge_trans)( LAPACK_COL_MAJOR, m, n, c_t, ldc_t, c, ldc );
        /* Release memory and exit */
        LAPACKE_free( c_t );
exit_level_2:
        LAPACKE_free( b_t );
exit_level_1:
        LAPACKE_free( a_t );
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_strsyl3_work", info );
        }
    } else {
        info = -1;
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_strsyl3_work", info );
    }
    return info;
}
