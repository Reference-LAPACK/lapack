/*****************************************************************************
  Copyright (c) 2022, Intel Corp.
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
* Contents: Native middle-level C interface to LAPACK function zlangb
* Author: Simon Märtens
*****************************************************************************/

#include "lapacke_utils.h"

double API_SUFFIX(LAPACKE_zlangb_work)( int matrix_layout, char norm, lapack_int n,
                            lapack_int kl, lapack_int ku,
                            const lapack_complex_double* ab, lapack_int ldab,
                            double* work )
{
    lapack_int info = 0;
    double res = 0.;
    if( matrix_layout == LAPACK_COL_MAJOR ) {
        /* Call LAPACK function and adjust info */
        res = LAPACK_zlangb( &norm, &n, &kl, &ku, ab, &ldab, work );
    } else if( matrix_layout == LAPACK_ROW_MAJOR ) {
        char norm_lapack;
        double* work_lapack = NULL;
        /* Check leading dimension(s) */
        if( ldab < kl+ku+1 ) {
            info = -7;
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_zlangb_work", info );
            return info;
        }
        if( API_SUFFIX(LAPACKE_lsame)( norm, '1' ) || API_SUFFIX(LAPACKE_lsame)( norm, 'o' ) ) {
            norm_lapack = 'i';
        } else if( API_SUFFIX(LAPACKE_lsame)( norm, 'i' ) ) {
            norm_lapack = '1';
        } else {
            norm_lapack = norm;
        }
        /* Allocate memory for work array(s) */
        if( API_SUFFIX(LAPACKE_lsame)( norm_lapack, 'i' ) ) {
            work_lapack = (double*)LAPACKE_malloc( sizeof(double) * MAX(1,n) );
            if( work_lapack == NULL ) {
                info = LAPACK_WORK_MEMORY_ERROR;
                goto exit_level_0;
            }
        }
        /* Call LAPACK function */
        res = LAPACK_zlangb( &norm, &n, &ku, &kl, ab, &ldab, work );
        /* Release memory and exit */
        if( work_lapack ) {
            LAPACKE_free( work_lapack );
        }
exit_level_0:
        if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
            API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_zlangb_work", info );
        }
    } else {
        info = -1;
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_zlangb_work", info );
    }
    return res;
}
