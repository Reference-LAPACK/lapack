*> \brief \b ZTRMMOOP computes an out of place triangular times general matrix multiplication
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     RECURSIVE SUBROUTINE ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB,
*    $ DIAG, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*        .. Scalar Arguments ..
*        COMPLEX*16 ALPHA, BETA
*        INTEGER M, N, LDA, LDB, LDC
*        CHARACTER SIDE, UPLO, TRANSA, TRANSB, DIAG
*        ..
*        .. Array Arguments ..
*        COMPLEX*16 A(LDA,*), B(LDB,*), C(LDC,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRMMOOP performs one of the matrix-matrix operations
*>
*>       C = \alpha op(A) * op(B) + \beta C
*>                      or
*>       C = \alpha op(B) * op(A) + \beta C
*>
*> where \alpha and \beta are scalars, C is an m-by-n matrix, A is
*> a unit, or non-unit, upper or lower triangular matrix, and op(A) is
*> is one of
*>
*>       op(A) = A or op(A) = A**T op(A) = A**H
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether op(A) multiplies op(B) from
*>           the left or right as follows:
*>
*>             SIDE = 'L' or 'l' C = \alpha op(A) * op(B) + \beta C
*>
*>             SIDE = 'R' or 'r' C = \alpha op(B) * op(A) + \beta C
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>             UPLO = 'U' or 'u' A is upper triangular
*>
*>             UPLO = 'L' or 'l' A is lower triangular
*> \Endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op(A) to be used in
*>           the matrix multiplication as follows:
*>             TRANSA = 'N' or 'n' op(A) = A
*>
*>             TRANSA = 'T' or 't' op(A) = A**T
*>
*>             TRANSA = 'C' or 'c' op(A) = A**H
*> \endverbatim
*>
*> \param[in] TRANSB
*> \verbatim
*>          TRANSB is CHARACTER*1
*>           On entry, TRANSB specifies the form of op(B) to be used in
*>           the matrix multiplication as follows:
*>             TRANSB = 'N' or 'n' op(B) = B
*>
*>             TRANSB = 'T' or 't' op(B) = B**T
*>
*>             TRANSB = 'C' or 'c' op(B) = B**H
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u' A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n' A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of C. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of C. N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16.
*>           On entry, ALPHA specifies the scalar alpha. When alpha is
*>           zero then A and B are not referenced, and A and B need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, K ) where
*>           K is M when SIDE = 'L' and K is N when SIDE='R'
*>           Before entry with UPLO = 'U' or 'u', the leading k-by-k
*>           upper triangular part of the array A must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading k-by-k
*>           lower triangular part of the array A must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when DIAG = 'U' or 'u', the diagonal elements of
*>           A are not referenced either, but are assumed to be unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. When SIDE = 'L' or 'l' then
*>           LDA must be at least max( 1, m ), when SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>           B is COMPLEX*16 array, dimension ( LDB, K ), where K is M
*>           If SIDE='R' and TRANSA='N', or SIDE='L' and TRANSA='T' and N
*>           otherwise. On entry, the leading k-by-k sub-matrix must contain
*>           B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in the calling (sub) program. When SIDE = 'R' and TRANSB='N'
*>           then LDB must be at least max( 1, m ), when SIDE = 'R'
*>           and TRANSB = 'T' then LDB must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16.
*>           On entry, BETA specifies the scalar beta. When beta is
*>           zero then C is not referenced on entry, and C need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension ( LDC, N )
*>           Before entry, the leading m-by-n part of the array C must
*>           contain the matrix C, and on exit is overwritten by the
*>           transformed matrix.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in the calling (sub) program. LDC must be at least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*  =====================================================================
      RECURSIVE SUBROUTINE ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB,
     $         DIAG, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, BETA
      INTEGER           M, N, LDA, LDB, LDC
      CHARACTER         SIDE, UPLO, TRANSA, TRANSB, DIAG
*     ..
*     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL           LSAME
      COMPLEX*16        ZDOTC, ZDOTU
      EXTERNAL          LSAME, ZDOTC, ZDOTU
*     ..
*     .. External Subroutines ..
      EXTERNAL          ZTRMMOOP_LVL2, ZGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC         CONJG, MIN
*     ..
*     .. Local Scalars ..
      INTEGER           I, J, L, INCB
      LOGICAL           LSIDE, UPPER, UNIT, TRANST, TRANSG,
     $                  CONJA, CONJB
*     ..
*     .. Local Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER(ONE=(1.0D+0,0.0D+0), ZERO=(0.0D+0,0.0D+0))
*     ..
*
*     Beginning of Executable Statements
*
      LSIDE = LSAME(SIDE, 'L')
      UPPER = LSAME(UPLO, 'U')
*
*     If we are transposing the triangular matrix (A)
*
      CONJA = LSAME(TRANSA, 'C')
      TRANST= LSAME(TRANSA, 'T').OR.CONJA
*
*     If we are transposing the general matrix (B)
*
      CONJB = LSAME(TRANSB, 'C')
      TRANSG= LSAME(TRANSB, 'T').OR.CONJB
*
*     Terminating Case
*
      UNIT  = LSAME(DIAG, 'U')
      IF (M.EQ.1.OR.N.EQ.1) THEN
            ! Let the level 2 implementation handle this termination case
            CALL ZTRMMOOP_LVL2(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $         M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
         RETURN
      END IF
*
*     Recursive Case
*
      L = MIN(M,N)/2
      IF (LSIDE) THEN
*
*        We are multiplying A from the left IE we are computing
*        C = \alpha op(A)*op(B) + \beta C
*
         IF (UPPER) THEN
*
*           A is upper triangular
*
            IF (TRANST) THEN
*
*              We are transposing A
*
               IF (TRANSG) THEN
*
*                 We are transposing
*
*                 So we are computing
*                 C = \alpha A**T * B**T + \beta C. We break this down as follows
*
*                       |-------------|         |-------------------|
*                 C =   |C_{11} C_{12}| A**T =  |A_{11}**T 0        |
*                       |C_{21} C_{22}|         |A_{12}**T A_{22}**T|
*                       |-------------|         |-------------------|
*
*                       |-------------------|
*                 B**T =|B_{11}**T B_{21}**T|
*                       |B_{12}**T B_{22}**T|
*                       |-------------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times m-\ell}
*                 A_{21}\in\R^{m-\ell\times\ell} A_{22}\in\R^{m-\ell\times m-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times m-\ell}
*                 B_{21}\in\R^{n-\ell\times\ell} B_{22}\in\R^{n-\ell\times m-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha A_{11}**T * B_{11}**T + \beta C_{11}
*                 C_{12} = \alpha A_{11}**T * B_{21}**T + \beta C_{12}
*                 C_{21} = \alpha A_{12}**T * B_{11}**T + \alpha A_{22}**T * B_{12}**T + \beta C_{21}
*                 C_{22} = \alpha A_{12}**T * B_{21}**T + \alpha A_{22}**T * B_{22}**T + \beta C_{22}
*
*                 Computing C_{12} and C_{12} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{21} and C_{22} as follows
*
*                 C_{21} = \alpha A_{12}**T * B_{11}**T + \beta C_{21} (GEMM call)
*                 C_{21} = \alpha A_{22}**T * B_{12}**T + C_{21} (This routine)
*
*                 C_{22} = \alpha A_{12}**T * B_{21}**T + \beta C_{22} (GEMM call)
*                 C_{22} = \alpha A_{22}**T * B_{22}**T + C_{22} (This routine)
*
*                 C_{11}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                     L, L, ALPHA, A, LDA, B, LDB, BETA, C,
     $                     LDC)
*
*                 C_{12}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                     L, N-L, ALPHA, A, LDA, B(L+1, 1), LDB,
     $                     BETA, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZGEMM(TRANSA, TRANSB, M-L, L, L, ALPHA,
     $                     A(1, L+1), LDA, B, LDB, BETA, C(L+1,1),
     $                     LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                     M-L, L, ALPHA, A(L+1,L+1), LDA, B(1,L+1),
     $                     LDB, ONE, C(L+1,1), LDC)
*
*                 C_{22}
*
                  CALL ZGEMM(TRANSA, TRANSB, M-L, N-L, L, ALPHA,
     $                     A(1, L+1), LDA, B(L+1,1), LDB, BETA,
     $                     C(L+1,L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                     M-L, N-L, ALPHA, A(L+1,L+1), LDA,
     $                     B(L+1,L+1), LDB, ONE, C(L+1,L+1), LDC)
               ELSE
*
*                 We are not transposing B.
*
*                 So we are computing
*                 C = \alpha A**T * B + \beta C. We break this down as follows
*
*                       |-------------|         |-------------------|
*                 C =   |C_{11} C_{12}| A**T =  |A_{11}**T 0        |
*                       |C_{21} C_{22}|         |A_{12}**T A_{22}**T|
*                       |-------------|         |-------------------|
*
*                       |-------------|
*                 B =   |B_{11} B_{12}|
*                       |B_{21} B_{22}|
*                       |-------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times m-\ell}
*                 A_{21}\in\R^{m-\ell\times\ell} A_{22}\in\R^{m-\ell\times m-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times n-\ell}
*                 B_{21}\in\R^{m-\ell\times\ell} B_{22}\in\R^{m-\ell\times n-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha A_{11}**T * B_{11} + \beta C_{11}
*                 C_{12} = \alpha A_{11}**T * B_{12} + \beta C_{12}
*                 C_{21} = \alpha A_{12}**T * B_{11} + \alpha A_{22}**T * B_{21} + \beta C_{21}
*                 C_{22} = \alpha A_{12}**T * B_{12} + \alpha A_{22}**T * B_{22} + \beta C_{22}
*
*                 Computing C_{11} and C_{12} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{21} and C_{22} as follows
*
*                 C_{21} = \alpha A_{12}**T * B_{11} + \beta C_{21} (GEMM call)
*                 C_{21} = \alpha A_{22}**T * B_{21} + C_{21} (This routine)
*
*                 C_{22} = \alpha A_{12}**T * B_{12} + \beta C_{22} (GEMM call)
*                 C_{22} = \alpha A_{22}**T * B_{22} + C_{22} (This routine)
*
*                 C_{11}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*                 C_{12}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A, LDA, B(1, L+1), LDB, BETA,
     $                  C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZGEMM(TRANSA, TRANSB, M-L, L, L, ALPHA,
     $                  A(1, L+1), LDA, B, LDB, BETA, C(L+1, 1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A(L+1, L+1), LDA, B(L+1, 1),
     $                  LDB, ONE, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZGEMM(TRANSA, TRANSB, M-L, N-L, L,
     $                  ALPHA, A(1, L+1), LDA, B(1, L+1), LDB, BETA,
     $                  C(L+1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1,L+1), LDB, ONE, C(L+1,L+1), LDC)
               ENDIF
            ELSE
*
*              We are not transposing A
*
               IF (TRANSG) THEN
*
*                 We are transposing B.
*
*                 So we are computing
*                 C = \alpha A * B**T + \beta C. We break this down as follows
*
*                       |-------------|      |-------------|
*                 C =   |C_{11} C_{12}| A =  |A_{11} A_{12}|
*                       |C_{21} C_{22}|      |0      A_{22}|
*                       |-------------|      |-------------|
*
*                       |-------------------|
*                 B**T =|B_{11}**T B_{21}**T|
*                       |B_{12}**T B_{22}**T|
*                       |-------------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times m-\ell}
*                 A_{21}\in\R^{m-\ell\times\ell} A_{22}\in\R^{m-\ell\times m-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times m-\ell}
*                 B_{21}\in\R^{n-\ell\times\ell} B_{22}\in\R^{n-\ell\times m-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha A_{11} * B_{11}**T + \alpha A_{12} * B_{12}**T + \beta C_{11}
*                 C_{12} = \alpha A_{11} * B_{21}**T + \alpha A_{12} * B_{22}**T + \beta C_{12}
*                 C_{21} = \alpha A_{22} * B_{12}**T + \beta C_{21}
*                 C_{22} = \alpha A_{22} * B_{22}**T + \beta C_{22}
*
*                 Computing C_{21} and C_{22} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{11} and C_{12} as follows
*
*                 C_{11} = \alpha A_{12} * B_{12}**T + \beta C_{11} (GEMM call)
*                 C_{11} = \alpha A_{11} * B_{11}**T + C_{11} (This routine)
*
*                 C_{12} = \alpha A_{12} * B_{22}**T + \beta C_{12} (GEMM call)
*                 C_{12} = \alpha A_{11} * B_{21}**T + C_{12} (This routine)
*
*                 C_{11}
*
                  CALL ZGEMM(TRANSA, TRANSB, L, L, M-L, ALPHA,
     $                  A(1, L+1), LDA, B(1, L+1), LDB, BETA, C, LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, ONE, C, LDC)
*
*                 C_{12}
*
                  CALL ZGEMM(TRANSA, TRANSB, L, N-L, M-L, ALPHA,
     $                  A(1, L+1), LDA, B(L+1, L+1), LDB, BETA,
     $                  C(1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A, LDA, B(L+1,1), LDB, ONE,
     $                  C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A(L+1, L+1), LDA, B(1, L+1),
     $                  LDB, BETA, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, BETA, C(L+1, L+1), LDC)
               ELSE
*
*                 We are not transposing B.
*
*                 So we are computing
*                 C = \alpha A * B + \beta C. We break this down as follows
*
*                       |-------------|      |-------------|
*                 C =   |C_{11} C_{12}| A =  |A_{11} A_{12}|
*                       |C_{21} C_{22}|      |0      A_{22}|
*                       |-------------|      |-------------|
*
*                       |-------------|
*                 B =   |B_{11} B_{12}|
*                       |B_{21} B_{22}|
*                       |-------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times m-\ell}
*                 A_{21}\in\R^{m-\ell\times\ell} A_{22}\in\R^{m-\ell\times m-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times n-\ell}
*                 B_{21}\in\R^{m-\ell\times\ell} B_{22}\in\R^{m-\ell\times n-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha A_{11} * B_{11} + \alpha A_{12} * B_{21} + \beta C_{11}
*                 C_{12} = \alpha A_{11} * B_{12} + \alpha A_{12} * B_{22} + \beta C_{12}
*                 C_{21} = \alpha A_{22} * B_{21} + \beta C_{21}
*                 C_{22} = \alpha A_{22} * B_{22} + \beta C_{22}
*
*                 Computing C_{21} and C_{22} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{11} and C_{12} as follows
*
*                 C_{11} = \alpha A_{12} * B_{21} + \beta C_{11} (GEMM call)
*                 C_{11} = \alpha A_{11} * B_{11} + C_{11} (This routine)
*
*                 C_{12} = \alpha A_{12} * B_{22} + \beta C_{12} (GEMM call)
*                 C_{12} = \alpha A_{11} * B_{12} + C_{12} (This routine)
*
*                 C_{11}
*
                  CALL ZGEMM(TRANSA, TRANSB, L, L, M-L, ALPHA,
     $                  A(1, L+1), LDA, B(L+1, 1), LDB, BETA, C, LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, ONE, C, LDC)
*
*                 C_{12}
*
                  CALL ZGEMM(TRANSB, TRANSA, L, N-L, M-L, ALPHA,
     $                  A(1, L+1), LDA, B(L+1, L+1), LDB, BETA,
     $                  C(1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A, LDA, B(1, L+1), LDB,
     $                  ONE, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A(L+1, L+1), LDA, B(L+1, 1),
     $                  LDB, BETA, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, BETA, C(L+1, L+1), LDC)
               ENDIF
            END IF
         ELSE
*
*           A is lower triangular
*
            IF (TRANST) THEN
*
*              We are transposing A
*
               IF (TRANSG) THEN
*
*                 We are transposing B.
*
*                 So we are computing
*                 C = \alpha A**T * B**T + \beta C. We break this down as follows
*
*                       |-------------|         |-------------------|
*                 C =   |C_{11} C_{12}| A**T =  |A_{11}**T A_{21}**T|
*                       |C_{21} C_{22}|         |0         A_{22}**T|
*                       |-------------|         |-------------------|
*
*                       |-------------------|
*                 B**T =|B_{11}**T B_{21}**T|
*                       |B_{12}**T B_{22}**T|
*                       |-------------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times m-\ell}
*                 A_{21}\in\R^{m-\ell\times\ell} A_{22}\in\R^{m-\ell\times m-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times m-\ell}
*                 B_{21}\in\R^{n-\ell\times\ell} B_{22}\in\R^{n-\ell\times m-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha A_{11}**T * B_{11}**T + \alpha A_{21}**T * B_{12}**T + \beta C_{11}
*                 C_{12} = \alpha A_{11}**T * B_{21}**T + \alpha A_{21}**T * B_{22}**T + \beta C_{12}
*                 C_{21} = \alpha A_{22}**T * B_{12}**T + \beta C_{21}
*                 C_{22} = \alpha A_{22}**T * B_{22}**T + \beta C_{22}
*
*                 Computing C_{21} and C_{22} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{11} and C_{12} as follows
*
*                 C_{11} = \alpha A_{21}**T * B_{12}**T + \beta C_{11} (GEMM call)
*                 C_{11} = \alpha A_{11}**T * B_{11}**T + C_{11} (This routine)
*
*                 C_{12} = \alpha A_{21}**T * B_{22}**T + \beta C_{12} (GEMM call)
*                 C_{12} = \alpha A_{11}**T * B_{21}**T + C_{12} (This routine)
*
*                 C_{11}
*
                  CALL ZGEMM(TRANSA, TRANSB, L, L, M-L, ALPHA,
     $                  A(L+1, 1), LDA, B(1, L+1), LDB, BETA, C, LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, ONE, C, LDC)
*
*                 C_{12}
*
                  CALL ZGEMM(TRANSA, TRANSB, L, N-L, M-L, ALPHA,
     $                  A(L+1, 1), LDA, B(L+1, L+1), LDB, BETA,
     $                  C(1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A, LDA, B(L+1, 1), LDB, ONE,
     $                  C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A(L+1, L+1), LDA, B(1, L+1),
     $                  LDB, BETA, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, BETA, C(L+1, L+1), LDC)
               ELSE
*
*                 We are not transposing B.
*
*                 So we are computing
*                 C = \alpha A**T * B + \beta C. We break this down as follows
*
*                       |-------------|         |-------------------|
*                 C =   |C_{11} C_{12}| A**T =  |A_{11}**T A_{21}**T|
*                       |C_{21} C_{22}|         |0         A_{22}**T|
*                       |-------------|         |-------------------|
*
*                       |-------------|
*                 B =   |B_{11} B_{12}|
*                       |B_{21} B_{22}|
*                       |-------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times m-\ell}
*                 A_{21}\in\R^{m-\ell\times\ell} A_{22}\in\R^{m-\ell\times m-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times n-\ell}
*                 B_{21}\in\R^{m-\ell\times\ell} B_{22}\in\R^{m-\ell\times n-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha A_{11}**T * B_{11} + \alpha A_{21}**T * B_{21} + \beta C_{11}
*                 C_{12} = \alpha A_{11}**T * B_{12} + \alpha A_{21}**T * B_{22} + \beta C_{12}
*                 C_{21} = \alpha A_{22}**T * B_{21} + \beta C_{21}
*                 C_{22} = \alpha A_{22}**T * B_{22} + \beta C_{22}
*
*                 Computing C_{21} and C_{22} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{11} and C_{12} as follows
*
*                 C_{11} = \alpha A_{21}**T * B_{21} + \beta C_{11} (GEMM call)
*                 C_{11} = \alpha A_{11}**T * B_{11} + C_{11} (This routine)
*
*                 C_{12} = \alpha A_{21}**T * B_{22} + \beta C_{12} (GEMM call)
*                 C_{12} = \alpha A_{11}**T * B_{12} + C_{12} (This routine)
*
*                 C_{11}
*
                  CALL ZGEMM(TRANSA, TRANSB, L, L, M-L, ALPHA,
     $                  A(L+1, 1), LDA, B(L+1, 1), LDB, BETA, C, LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, ONE, C, LDC)
*
*                 C_{12}
*
                  CALL ZGEMM(TRANSA, TRANSB, L, N-L, M-L, ALPHA,
     $                  A(L+1, 1), LDA, B(L+1, L+1), LDB, BETA,
     $                  C(1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A, LDA, B(1, L+1), LDB, ONE,
     $                  C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A(L+1, L+1), LDA, B(L+1, 1),
     $                  LDB, BETA, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, BETA, C(L+1, L+1), LDC)
               ENDIF
            ELSE
*
*              We are not transposing A
*
               IF (TRANSG) THEN
*
*                 We are transposing B.
*
*                 So we are computing
*                 C = \alpha A * B**T + \beta C. We break this down as follows
*
*                       |-------------|      |-------------|
*                 C =   |C_{11} C_{12}| A =  |A_{11} 0     |
*                       |C_{21} C_{22}|      |A_{21} A_{22}|
*                       |-------------|      |-------------|
*
*                       |-------------------|
*                 B**T =|B_{11}**T B_{21}**T|
*                       |B_{12}**T B_{22}**T|
*                       |-------------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times m-\ell}
*                 A_{21}\in\R^{m-\ell\times\ell} A_{22}\in\R^{m-\ell\times m-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times m-\ell}
*                 B_{21}\in\R^{n-\ell\times\ell} B_{22}\in\R^{n-\ell\times m-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha A_{11} * B_{11}**T + \beta C_{11}
*                 C_{12} = \alpha A_{11} * B_{21}**T + \beta C_{12}
*                 C_{21} = \alpha A_{21} * B_{11}**T + \alpha A_{22} * B_{12}**T + \beta * C_{21}
*                 C_{22} = \alpha A_{21} * B_{21}**T + \alpha A_{22} * B_{22}**T + \beta * C_{22}
*
*                 Computing C_{11} and C_{12} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{21} and C_{22} as follows
*
*                 C_{21} = \alpha A_{21} * B_{11}**T + \beta C_{21} (GEMM call)
*                 C_{21} = \alpha A_{22} * B_{12}**T + C_{21} (This routine)
*
*                 C_{22} = \alpha A_{21} * B_{21}**T + \beta C_{22} (GEMM call)
*                 C_{22} = \alpha A_{22} * B_{22}**T + C_{22} (This routine)
*
*                 C_{11}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*                 C_{12}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A, LDA, B(L+1, 1), LDB,
     $                  BETA, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZGEMM(TRANSA, TRANSB, M-L, L, L, ALPHA,
     $                  A(L+1, 1), LDA, B, LDB, BETA, C(L+1, 1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A(L+1, L+1), LDA, B(1, L+1),
     $                  LDB, ONE, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZGEMM(TRANSA, TRANSB, M-L, N-L, L,
     $                  ALPHA, A(L+1, 1), LDA, B(L+1, 1), LDB, BETA,
     $                  C(L+1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, ONE, C(L+1, L+1), LDC)
               ELSE
*
*                 We are not transposing B.
*
*                 So we are computing
*                 C = \alpha A * B + \beta C. We break this down as follows
*
*                       |-------------|      |-------------|
*                 C =   |C_{11} C_{12}| A =  |A_{11} 0     |
*                       |C_{21} C_{22}|      |A_{21} A_{22}|
*                       |-------------|      |-------------|
*
*                       |-------------|
*                 B =   |B_{11} B_{12}|
*                       |B_{21} B_{22}|
*                       |-------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times m-\ell}
*                 A_{21}\in\R^{m-\ell\times\ell} A_{22}\in\R^{m-\ell\times m-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times n-\ell}
*                 B_{21}\in\R^{m-\ell\times\ell} B_{22}\in\R^{m-\ell\times n-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha A_{11} * B_{11} + \beta C_{11}
*                 C_{12} = \alpha A_{11} * B_{12} + \beta C_{12}
*                 C_{21} = \alpha A_{21} * B_{11} + \alpha A_{22} * B_{21} + \beta * C_{21}
*                 C_{22} = \alpha A_{21} * B_{12} + \alpha A_{22} * B_{22} + \beta * C_{22}
*
*                 Computing C_{11} and C_{12} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{21} and C_{22} as follows
*
*                 C_{21} = \alpha A_{21} * B_{11} + \beta C_{21} (GEMM call)
*                 C_{21} = \alpha A_{22} * B_{21} + C_{21} (This routine)
*
*                 C_{22} = \alpha A_{21} * B_{12} + \beta C_{22} (GEMM call)
*                 C_{22} = \alpha A_{22} * B_{22} + C_{22} (This routine)
*
*                 C_{11}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*                 C_{12}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A, LDA, B(1, L+1), LDB,
     $                  BETA, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZGEMM(TRANSA, TRANSB, M-L, L, L, ALPHA,
     $                  A(L+1, 1), LDA, B, LDB, BETA, C(L+1, 1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A(L+1, L+1), LDA, B(L+1, 1),
     $                  LDB, ONE, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZGEMM(TRANSB, TRANSA, M-L, N-L, L,
     $                  ALPHA, A(L+1, 1), LDA, B(1, L+1), LDB, BETA,
     $                  C(L+1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, ONE, C(L+1, L+1), LDC)
               ENDIF
            END IF
         END IF
      ELSE
*
*        We are multiplying A from the right IE we are computing
*        C = \alpha op(B)*op(A) + \beta C
*
         IF (UPPER) THEN
*
*           A is upper triangular
*
            IF (TRANST) THEN
*
*              We are transposing A
*
               IF (TRANSG) THEN
*
*                 We are transposing B.
*
*                 So we are computing
*                 C = \alpha  B**T * A**T + \beta C. We break this down as follows
*
*                       |-------------|         |-------------------|
*                 C =   |C_{11} C_{12}| A**T =  |A_{11}**T 0        |
*                       |C_{21} C_{22}|         |A_{12}**T A_{22}**T|
*                       |-------------|         |-------------------|
*
*                       |-------------------|
*                 B**T =|B_{11}**T B_{21}**T|
*                       |B_{12}**T B_{22}**T|
*                       |-------------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times n-\ell}
*                 A_{21}\in\R^{n-\ell\times\ell} A_{22}\in\R^{n-\ell\times n-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times m-\ell}
*                 B_{21}\in\R^{n-\ell\times\ell} B_{22}\in\R^{n-\ell\times m-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha B_{11}**T * A_{11}**T + \alpha B_{21}**T * A_{12}**T + \beta C_{11}
*                 C_{12} = \alpha B_{21}**T * A_{22}**T + \beta C_{12}
*                 C_{21} = \alpha B_{12}**T * A_{11}**T + \alpha B_{22}**T * A_{12}**T + \beta C_{21}
*                 C_{22} = \alpha B_{22}**T * A_{22}**T + \beta C_{22}
*
*                 Computing C_{12} and C_{22} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{11} and C_{21} as follows
*
*                 C_{11} = \alpha B_{21}**T * A_{12}**T + \beta C_{11} (GEMM call)
*                 C_{11} = \alpha B_{11}**T * A_{11}**T + C_{11} (This routine)
*
*                 C_{21} = \alpha B_{22}**T * A_{12}**T + \beta C_{21} (GEMM call)
*                 C_{21} = \alpha B_{12}**T * A_{11}**T + C_{21} (This routine)
*
*                 C_{11}
*
                  CALL ZGEMM(TRANSB, TRANSA, L, L, N-L, ALPHA,
     $                  B(L+1, 1), LDB, A(1, L+1), LDA, BETA, C, LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, ONE, C, LDC)
*
*                 C_{12}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A(L+1, L+1), LDA, B(L+1, 1),
     $                  LDB, BETA, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZGEMM(TRANSB, TRANSA, M-L, L, N-L, ALPHA,
     $                  B(L+1, L+1), LDB, A(1, L+1), LDA, BETA,
     $                  C(L+1, 1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A, LDA, B(1, L+1), LDB,
     $                  ONE, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, BETA, C(L+1, L+1), LDC)
               ELSE
*
*                 We are not transposing B.
*
*                 So we are computing
*                 C = \alpha B * A**T + \beta C. We break this down as follows
*
*                       |-------------|         |-------------------|
*                 C =   |C_{11} C_{12}| A**T =  |A_{11}**T 0        |
*                       |C_{21} C_{22}|         |A_{12}**T A_{22}**T|
*                       |-------------|         |-------------------|
*
*                       |-------------|
*                 B =   |B_{11} B_{12}|
*                       |B_{21} B_{22}|
*                       |-------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times n-\ell}
*                 A_{21}\in\R^{n-\ell\times\ell} A_{22}\in\R^{n-\ell\times n-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times n-\ell}
*                 B_{21}\in\R^{m-\ell\times\ell} B_{22}\in\R^{m-\ell\times n-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha B_{11} * A_{11}**T + \alpha B_{12} * A_{12}**T + \beta C_{11}
*                 C_{12} = \alpha B_{12} * A_{22}**T + \beta C_{12}
*                 C_{21} = \alpha B_{21} * A_{11}**T + \alpha B_{22} * A_{12}**T + \beta C_{21}
*                 C_{22} = \alpha B_{22} * A_{22}**T + \beta C_{22}
*
*                 Computing C_{12} and C_{22} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{11} and C_{21} as follows
*
*                 C_{11} = \alpha B_{12} * A_{12}**T + \beta C_{11} (GEMM call)
*                 C_{11} = \alpha B_{11} * A_{11}**T + C_{11} (This routine)
*
*                 C_{21} = \alpha B_{22} * A_{12}**T + \beta C_{21} (GEMM call)
*                 C_{21} = \alpha B_{21} * A_{11}**T + C_{21} (This routine)
*
*                 C_{11}
*
                  CALL ZGEMM(TRANSB, TRANSA, L, L, N-L, ALPHA,
     $                  B(1,L+1), LDB, A(1,L+1), LDA, BETA, C, LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, ONE, C, LDC)
*
*                 C_{12}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A(L+1, L+1), LDA, B(1, L+1),
     $                  LDB, BETA, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZGEMM(TRANSB, TRANSA, M-L, L, N-L, ALPHA,
     $                  B(L+1, L+1), LDB, A(1, L+1), LDA, BETA,
     $                  C(L+1, 1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A, LDA, B(L+1, 1), LDB,
     $                  ONE, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, BETA, C(L+1, L+1), LDC)
               ENDIF
            ELSE
*
*              We are not transposing A
*
               IF (TRANSG) THEN
*
*                 We are transposing B.
*
*                 So we are computing
*                 C = \alpha B**T * A + \beta C. We break this down as follows
*
*                       |-------------|      |-------------|
*                 C =   |C_{11} C_{12}| A =  |A_{11} A_{12}|
*                       |C_{21} C_{22}|      |0      A_{22}|
*                       |-------------|      |-------------|
*
*                       |-------------------|
*                 B**T =|B_{11}**T B_{21}**T|
*                       |B_{12}**T B_{22}**T|
*                       |-------------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times n-\ell}
*                 A_{21}\in\R^{n-\ell\times\ell} A_{22}\in\R^{n-\ell\times n-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times m-\ell}
*                 B_{21}\in\R^{n-\ell\times\ell} B_{22}\in\R^{n-\ell\times m-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha B_{11}**T * A_{11} + \beta C_{11}
*                 C_{12} = \alpha B_{11}**T * A_{12} + \alpha B_{21}**T * A_{22} + \beta C_{12}
*                 C_{21} = \alpha B_{12}**T * A_{11} + \beta C_{21}
*                 C_{22} = \alpha B_{12}**T * A_{12} + \alpha B_{22}**T * A_{22} + \beta C_{22}
*
*                 Computing C_{11} and C_{21} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{12} and C_{22} as follows
*
*                 C_{12} = \alpha B_{11}**T * A_{12} + \beta C_{12} (GEMM call)
*                 C_{12} = \alpha B_{21}**T * A_{22} + C_{12} (This routine)
*
*                 C_{22} = \alpha B_{12}**T * A_{12} + \beta C_{22} (GEMM call)
*                 C_{22} = \alpha B_{22}**T * A_{22} + C_{22} (This routine)
*
*                 C_{11}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*                 C_{12}
*
                  CALL ZGEMM(TRANSB, TRANSA, L, N-L, L, ALPHA,
     $                  B, LDB, A(1, L+1), LDA, BETA, C(1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A(L+1, L+1), LDA, B(L+1, 1),
     $                  LDB, ONE, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A, LDA, B(1, L+1), LDB,
     $                  BETA, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZGEMM(TRANSB, TRANSA, M-L, N-L, L,
     $                  ALPHA, B(1, L+1), LDB, A(1, L+1), LDA, BETA,
     $                  C(L+1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, ONE, C(L+1, L+1), LDC)
               ELSE
*
*                 We are not transposing B.
*
*                 So we are computing
*                 C = \alpha B * A + \beta C. We break this down as follows
*
*                       |-------------|      |-------------|
*                 C =   |C_{11} C_{12}| A =  |A_{11} A_{12}|
*                       |C_{21} C_{22}|      |0      A_{22}|
*                       |-------------|      |-------------|
*
*                       |-------------|
*                 B =   |B_{11} B_{12}|
*                       |B_{21} B_{22}|
*                       |-------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times n-\ell}
*                 A_{21}\in\R^{n-\ell\times\ell} A_{22}\in\R^{n-\ell\times n-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times n-\ell}
*                 B_{21}\in\R^{m-\ell\times\ell} B_{22}\in\R^{m-\ell\times n-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha B_{11} * A_{11} + \beta C_{11}
*                 C_{12} = \alpha B_{11} * A_{12} + \alpha B_{12} * A_{22} + \beta C_{12}
*                 C_{21} = \alpha B_{21} * A_{11} + \beta C_{21}
*                 C_{22} = \alpha B_{21} * A_{12} + \alpha B_{22} * A_{22} + \beta C_{22}
*
*                 Computing C_{11} and C_{21} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{12} and C_{22} as follows
*
*                 C_{12} = \alpha B_{11} * A_{12} + \beta C_{12} (GEMM call)
*                 C_{12} = \alpha B_{12} * A_{22} + C_{12} (This routine)
*
*                 C_{22} = \alpha B_{21} * A_{12} + \beta C_{22} (GEMM call)
*                 C_{22} = \alpha B_{22} * A_{22} + C_{22} (This routine)
*
*                 C_{11}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*                 C_{12}
*
                  CALL ZGEMM(TRANSB, TRANSA, L, N-L, L, ALPHA,
     $                  B, LDB, A(1, L+1), LDA, BETA, C(1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A(L+1, L+1), LDA, B(1, L+1),
     $                  LDB, ONE, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A, LDA, B(L+1, 1), LDB, BETA,
     $                  C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZGEMM(TRANSB, TRANSA, M-L, N-L, L,
     $                  ALPHA, B(L+1, 1), LDB, A(1, L+1), LDA,
     $                  BETA, C(L+1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, ONE, C(L+1, L+1), LDC)
               ENDIF
            END IF
         ELSE
*
*           A is lower triangular
*
            IF (TRANST) THEN
*
*              We are transposing A
*
               IF (TRANSG) THEN
*
*                 We are transposing B.
*
*                 So we are computing
*                 C = \alpha B**T * A**T + \beta C. We break this down as follows
*
*                       |-------------|         |-------------------|
*                 C =   |C_{11} C_{12}| A**T =  |A_{11}**T A_{21}**T|
*                       |C_{21} C_{22}|         |0         A_{22}**T|
*                       |-------------|         |-------------------|
*
*                       |-------------------|
*                 B**T =|B_{11}**T B_{21}**T|
*                       |B_{12}**T B_{22}**T|
*                       |-------------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times n-\ell}
*                 A_{21}\in\R^{n-\ell\times\ell} A_{22}\in\R^{n-\ell\times n-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times m-\ell}
*                 B_{21}\in\R^{n-\ell\times\ell} B_{22}\in\R^{n-\ell\times m-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha B_{11}**T * A_{11}**T + \beta C_{11}
*                 C_{12} = \alpha B_{11}**T * A_{21}**T + \alpha B_{21}**T * A_{22}**T + \beta C_{12}
*                 C_{21} = \alpha B_{12}**T * A_{11}**T + \beta C_{21}
*                 C_{22} = \alpha B_{12}**T * A_{21}**T + \alpha B_{22}**T * A_{22}**T + \beta C_{22}
*
*                 Computing C_{11} and C_{21} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{12} and C_{22} as follows
*
*                 C_{12} = \alpha B_{11}**T * A_{21}**T + \beta C_{12} (GEMM call)
*                 C_{12} = \alpha B_{21}**T * A_{22}**T + C_{12} (This routine)
*
*                 C_{22} = \alpha B_{12}**T * A_{21}**T + \beta C_{22} (GEMM call)
*                 C_{22} = \alpha B_{22}**T * A_{22}**T + C_{22} (This routine)
*
*                 C_{11}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*                 C_{12}
*
                  CALL ZGEMM(TRANSB, TRANSA, L, N-L, L, ALPHA,
     $                  B, LDB, A(L+1, 1), LDA, BETA, C(1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A(L+1, L+1), LDA, B(L+1, 1),
     $                  LDB, ONE, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A, LDA, B(1, L+1), LDB,
     $                  BETA, C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZGEMM(TRANSB, TRANSA, M-L, N-L, L, ALPHA,
     $                  B(1, L+1), LDB, A(L+1, 1), LDA, BETA,
     $                  C(L+1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, ONE, C(L+1, L+1), LDC)
               ELSE
*
*                 We are not transposing B.
*
*                 So we are computing
*                 C = \alpha B * A**T + \beta C. We break this down as follows
*
*                       |-------------|         |-------------------|
*                 C =   |C_{11} C_{12}| A**T =  |A_{11}**T A_{21}**T|
*                       |C_{21} C_{22}|         |0         A_{22}**T|
*                       |-------------|         |-------------------|
*
*                       |-------------|
*                 B =   |B_{11} B_{12}|
*                       |B_{21} B_{22}|
*                       |-------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times n-\ell}
*                 A_{21}\in\R^{n-\ell\times\ell} A_{22}\in\R^{n-\ell\times n-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times n-\ell}
*                 B_{21}\in\R^{m-\ell\times\ell} B_{22}\in\R^{m-\ell\times n-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha B_{11} * A_{11} + \beta C_{11}
*                 C_{12} = \alpha B_{11} * A_{21}**T + \alpha A_{12} * B_{22}**T + \beta C_{12}
*                 C_{21} = \alpha B_{21} * A_{11}**T + \beta C_{21}
*                 C_{22} = \alpha B_{21} * A_{21}**T + \alpha A_{22} * B_{22}**T + \beta C_{22}
*
*                 Computing C_{11} and C_{21} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{12} and C_{22} as follows
*
*                 C_{12} = \alpha B_{11} * A_{21}**T + \beta C_{12} (GEMM call)
*                 C_{12} = \alpha B_{12} * A_{22}**T + C_{12} (This routine)
*
*                 C_{22} = \alpha B_{21} * A_{21}**T + \beta C_{22} (GEMM call)
*                 C_{22} = \alpha B_{22} * A_{22}**T + C_{22} (This routine)
*
*                 C_{11}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*                 C_{12}
*
                  CALL ZGEMM(TRANSB, TRANSA, L, N-L, L, ALPHA,
     $                  B, LDB, A(L+1, 1), LDA, BETA, C(1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A(L+1, L+1), LDA, B(1, L+1),
     $                  LDB, ONE, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A, LDA, B(L+1, 1), LDB, BETA,
     $                  C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZGEMM(TRANSB, TRANSA, M-L, N-L, L, ALPHA,
     $                  B(L+1, 1), LDB, A(L+1, 1), LDA, BETA,
     $                  C(L+1, L+1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, ONE, C(L+1, L+1), LDC)
               ENDIF
            ELSE
*
*              We are not transposing A
*
               IF (TRANSG) THEN
*
*                 We are transposing B.
*
*                 So we are computing
*                 C = \alpha B**T * A + \beta C. We break this down as follows
*
*                       |-------------|      |-------------|
*                 C =   |C_{11} C_{12}| A =  |A_{11} 0     |
*                       |C_{21} C_{22}|      |A_{21} A_{22}|
*                       |-------------|      |-------------|
*
*                       |-------------------|
*                 B**T =|B_{11}**T B_{21}**T|
*                       |B_{12}**T B_{22}**T|
*                       |-------------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times n-\ell}
*                 A_{21}\in\R^{n-\ell\times\ell} A_{22}\in\R^{n-\ell\times n-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times m-\ell}
*                 B_{21}\in\R^{n-\ell\times\ell} B_{22}\in\R^{n-\ell\times m-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha B_{11}**T * A_{11} + \alpha B_{21}**T * A_{21} + \beta C_{11}
*                 C_{12} = \alpha B_{21}**T * A_{22} + \beta C_{12}
*                 C_{21} = \alpha B_{12}**T * A_{11} + \alpha B_{22}**T * A_{21} + \beta C_{21}
*                 C_{22} = \alpha B_{22}**T * A_{22} + \beta C_{22}
*
*                 Computing C_{12} and C_{22} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{11} and C_{21} as follows
*
*                 C_{11} = \alpha B_{21}**T * A_{21} + \beta C_{11} (GEMM call)
*                 C_{11} = \alpha B_{11}**T * A_{11} + C_{11}(This routine)
*
*                 C_{21} = \alpha B_{22}**T * A_{21} + \beta C_{21} (GEMM call)
*                 C_{21} = \alpha B_{12}**T * A_{11} + C_{21} (This routine)
*
*                 C_{11}
*
                  CALL ZGEMM(TRANSB, TRANSA, L, L, N-L, ALPHA,
     $                  B(L+1, 1), LDB, A(L+1, 1), LDA, BETA, C, LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, ONE, C, LDC)
*
*                 C_{12}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A(L+1, L+1), LDA, B(L+1, 1),
     $                  LDB, BETA, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZGEMM(TRANSB, TRANSA, M-L, L, N-L, ALPHA,
     $                  B(L+1, L+1), LDB, A(L+1, 1), LDA, BETA,
     $                  C(L+1, 1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A, LDA, B(1, L+1), LDB, ONE,
     $                  C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, BETA, C(L+1, L+1), LDC)
               ELSE
*
*                 We are not transposing B.
*
*                 So we are computing
*                 C = \alpha B * A + \beta C. We break this down as follows
*
*                       |-------------|      |-------------|
*                 C =   |C_{11} C_{12}| A =  |A_{11} 0     |
*                       |C_{21} C_{22}|      |A_{21} A_{22}|
*                       |-------------|      |-------------|
*
*                       |-------------|
*                 B =   |B_{11} B_{12}|
*                       |B_{21} B_{22}|
*                       |-------------|
*
*                 Where
*                 C_{11}\in\R^{\ell\times\ell}   C_{12}\in\R^{\ell\times n-\ell}
*                 C_{21}\in\R^{m-\ell\times\ell} C_{22}\in\R^{m-\ell\times n-\ell}
*
*                 A_{11}\in\R^{\ell\times\ell}   A_{12}\in\R^{\ell\times n-\ell}
*                 A_{21}\in\R^{n-\ell\times\ell} A_{22}\in\R^{n-\ell\times n-\ell}
*
*                 B_{11}\in\R^{\ell\times\ell}   B_{12}\in\R^{\ell\times n-\ell}
*                 B_{21}\in\R^{m-\ell\times\ell} B_{22}\in\R^{m-\ell\times n-\ell}
*
*                 Which means that we get
*                 C_{11} = \alpha B_{11} * A_{11} + \alpha B_{12} * A_{21} + \beta C_{11}
*                 C_{12} = \alpha B_{12} * A_{22} + \beta C_{12}
*                 C_{21} = \alpha B_{21} * A_{11} + \alpha B_{22} * A_{21} + \beta C_{21}
*                 C_{22} = \alpha B_{22} * A_{22} + \beta C_{22}
*
*                 Computing C_{12} and C_{22} is just a recursive call to
*                 this routine but we can break down computing
*                 C_{11} and C_{21} as follows
*
*                 C_{11} = \alpha B_{12} * A_{21} + \beta C_{11} (GEMM call)
*                 C_{11} = \alpha B_{11} * A_{11} + C_{11}(This routine)
*
*                 C_{21} = \alpha B_{22} * A_{21} + \beta C_{21} (GEMM call)
*                 C_{21} = \alpha B_{21} * A_{11} + C_{21} (This routine)
*
*                 C_{11}
*
                  CALL ZGEMM(TRANSB, TRANSA, L, L, N-L, ALPHA,
     $                  B(1, L+1), LDB, A(L+1, 1), LDA, BETA, C, LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, L, ALPHA, A, LDA, B, LDB, ONE, C, LDC)
*
*                 C_{12}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  L, N-L, ALPHA, A(L+1, L+1), LDA, B(1, L+1),
     $                  LDB, BETA, C(1, L+1), LDC)
*
*                 C_{21}
*
                  CALL ZGEMM(TRANSB, TRANSA, M-L, L, N-L, ALPHA,
     $                  B(L+1, L+1), LDB, A(L+1, 1), LDA, BETA,
     $                  C(L+1, 1), LDC)
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, L, ALPHA, A, LDA, B(L+1, 1), LDB, ONE,
     $                  C(L+1, 1), LDC)
*
*                 C_{22}
*
                  CALL ZTRMMOOP(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $                  M-L, N-L, ALPHA, A(L+1, L+1), LDA,
     $                  B(L+1, L+1), LDB, BETA, C(L+1, L+1), LDC)
               ENDIF
            END IF
         END IF
      END IF
      END SUBROUTINE
