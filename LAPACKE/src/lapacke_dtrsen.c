/*****************************************************************************
  Copyright (c) 2014, Intel Corp.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************
* Contents: Native high-level C interface to LAPACK function dtrsen
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_dtrsen)( int matrix_layout, char job, char compq,
                           const lapack_logical* select, lapack_int n,
                           double* t, lapack_int ldt, double* q, lapack_int ldq,
                           double* wr, double* wi, lapack_int* m, double* s,
                           double* sep )
{
    lapack_int info = 0;
    lapack_int liwork = -1;
    lapack_int lwork = -1;
    lapack_int* iwork = NULL;
    double* work = NULL;
    lapack_int iwork_query;
    double work_query;
    if( matrix_layout != LAPACK_COL_MAJOR && matrix_layout != LAPACK_ROW_MAJOR ) {
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_dtrsen", -1 );
        return -1;
    }
#ifndef LAPACK_DISABLE_NAN_CHECK
    if( LAPACKE_get_nancheck() ) {
        /* Optionally check input matrices for NaNs */
        if( API_SUFFIX(LAPACKE_lsame)( compq, 'v' ) ) {
            if( API_SUFFIX(LAPACKE_dge_nancheck)( matrix_layout, n, n, q, ldq ) ) {
                return -8;
            }
        }
        if( API_SUFFIX(LAPACKE_dge_nancheck)( matrix_layout, n, n, t, ldt ) ) {
            return -6;
        }
    }
#endif
    /* Query optimal working array(s) size */
    info = API_SUFFIX(LAPACKE_dtrsen_work)( matrix_layout, job, compq, select, n, t, ldt, q,
                                ldq, wr, wi, m, s, sep, &work_query, lwork,
                                &iwork_query, liwork );
    if( info != 0 ) {
        goto exit_level_0;
    }
    liwork = iwork_query;
    lwork = (lapack_int)work_query;
    /* Allocate memory for work arrays */
    if( API_SUFFIX(LAPACKE_lsame)( job, 'b' ) || API_SUFFIX(LAPACKE_lsame)( job, 'v' ) ) {
        iwork = (lapack_int*)LAPACKE_malloc( sizeof(lapack_int) * liwork );
        if( iwork == NULL ) {
            info = LAPACK_WORK_MEMORY_ERROR;
            goto exit_level_0;
        }
    } else {
        iwork = (lapack_int*)LAPACKE_malloc( sizeof(lapack_int) );
        if( iwork == NULL ) {
            info = LAPACK_WORK_MEMORY_ERROR;
            goto exit_level_0;
        }
    }
    work = (double*)LAPACKE_malloc( sizeof(double) * lwork );
    if( work == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_1;
    }
    /* Call middle-level interface */
    info = API_SUFFIX(LAPACKE_dtrsen_work)( matrix_layout, job, compq, select, n, t, ldt, q,
                                ldq, wr, wi, m, s, sep, work, lwork, iwork,
                                liwork );
    /* Release memory and exit */
    LAPACKE_free( work );
exit_level_1:
    if( API_SUFFIX(LAPACKE_lsame)( job, 'b' ) || API_SUFFIX(LAPACKE_lsame)( job, 'v' ) ) {
        LAPACKE_free( iwork );
    }
exit_level_0:
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_dtrsen", info );
    }
    return info;
}
