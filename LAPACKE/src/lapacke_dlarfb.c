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
* Contents: Native high-level C interface to LAPACK function dlarfb
* Author: Intel Corporation
*****************************************************************************/

#include "lapacke_utils.h"

lapack_int API_SUFFIX(LAPACKE_dlarfb)( int matrix_layout, char side, char trans, char direct,
                           char storev, lapack_int m, lapack_int n,
                           lapack_int k, const double* v, lapack_int ldv,
                           const double* t, lapack_int ldt, double* c,
                           lapack_int ldc )
{
    lapack_int info = 0;
    lapack_int ldwork;
    double* work = NULL;
    lapack_int nrows_v, ncols_v;
    lapack_logical left, col, forward;
    char uplo;
    if( matrix_layout != LAPACK_COL_MAJOR && matrix_layout != LAPACK_ROW_MAJOR ) {
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_dlarfb", -1 );
        return -1;
    }
#ifndef LAPACK_DISABLE_NAN_CHECK
    if( LAPACKE_get_nancheck() ) {
        /* Optionally check input matrices for NaNs */
        left = API_SUFFIX(LAPACKE_lsame)( side, 'l' );
        col = API_SUFFIX(LAPACKE_lsame)( storev, 'c' );
        forward = API_SUFFIX(LAPACKE_lsame)( direct, 'f' );

        nrows_v = ( col && left ) ? m : ( ( col && !left ) ? n : ( !col ? k : 1) );
        ncols_v = ( !col && left ) ? m : ( ( !col && !left ) ? n : ( col ? k : 1 ) );
        uplo = ( ( forward && col ) || !( forward || col ) ) ? 'l' : 'u';

        if( ( col && k > nrows_v ) || ( !col && k > ncols_v ) ) {
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_dlarfb", -8 );
            return -8;
        }
        if( API_SUFFIX(LAPACKE_dtz_nancheck)( matrix_layout, direct, uplo, 'u',
                                  nrows_v, ncols_v, v, ldv ) ) {
            return -9;
        }
        if( API_SUFFIX(LAPACKE_dge_nancheck)( matrix_layout, k, k, t, ldt ) ) {
            return -11;
        }
        if( API_SUFFIX(LAPACKE_dge_nancheck)( matrix_layout, m, n, c, ldc ) ) {
            return -13;
        }
    }
#endif
    if( API_SUFFIX(LAPACKE_lsame)( side, 'l' ) ) {
        ldwork = n;
    } else if( API_SUFFIX(LAPACKE_lsame)( side, 'r' ) ) {
        ldwork = m;
    } else {
        ldwork = 1;
    }
    /* Allocate memory for working array(s) */
    work = (double*)LAPACKE_malloc( sizeof(double) * ldwork * MAX(1,k) );
    if( work == NULL ) {
        info = LAPACK_WORK_MEMORY_ERROR;
        goto exit_level_0;
    }
    /* Call middle-level interface */
    info = API_SUFFIX(LAPACKE_dlarfb_work)( matrix_layout, side, trans, direct, storev, m, n,
                                k, v, ldv, t, ldt, c, ldc, work, ldwork );
    /* Release memory and exit */
    LAPACKE_free( work );
exit_level_0:
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_dlarfb", info );
    }
    return info;
}
