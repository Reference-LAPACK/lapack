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
* Contents: Native high-level C interface to LAPACK function clangb
* Author: Simon Märtens
*****************************************************************************/

#include "lapacke_utils.h"

float API_SUFFIX(LAPACKE_clangb)( int matrix_layout, char norm, lapack_int n,
                      lapack_int kl, lapack_int ku,
                      const lapack_complex_float* ab, lapack_int ldab )
{
    lapack_int info = 0;
    float res = 0.;
    float* work = NULL;
    if( matrix_layout != LAPACK_COL_MAJOR && matrix_layout != LAPACK_ROW_MAJOR ) {
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_clangb", -1 );
        return -1;
    }
#ifndef LAPACK_DISABLE_NAN_CHECK
    if( LAPACKE_get_nancheck() ) {
        /* Optionally check input matrices for NaNs */
        if( API_SUFFIX(LAPACKE_cgb_nancheck)( matrix_layout, n, n, kl, ku, ab, ldab ) ) {
            return -6;
        }
    }
#endif
    /* Allocate memory for working array(s) */
    if( API_SUFFIX(LAPACKE_lsame)( norm, 'i' ) ) {
        work = (float*)LAPACKE_malloc( sizeof(float) * MAX(1,n) );
        if( work == NULL ) {
            info = LAPACK_WORK_MEMORY_ERROR;
            goto exit_level_0;
        }
    }
    /* Call middle-level interface */
    res = API_SUFFIX(LAPACKE_clangb_work)( matrix_layout, norm, n, kl, ku, ab, ldab, work );
    /* Release memory and exit */
    if( API_SUFFIX(LAPACKE_lsame)( norm, 'i' ) ) {
        LAPACKE_free( work );
    }
exit_level_0:
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        API_SUFFIX(LAPACKE_xerbla)( "LAPACKE_clangb", info );
    }
    return res;
}
