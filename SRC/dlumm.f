*> \brief \b DLUMM computes an in place triangular times triangluar matrix multiplication
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     RECURSIVE SUBROUTINE DLUMM(SIDEL, DIAGL, DIAGU, N, ALPHA,
*    $                        A, LDA)
*
*     .. Scalar Arguments ..
*     INTEGER           N, LDA
*     CHARACTER         SIDEL, DIAGL, DIAGU
*     DOUBLE PRECISION  ALPHA
*
*     .. Array Arguments ..
*     DOUBLE PRECISION  A(LDA,*)
*     ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLUMM performs one of the matrix-matrix operations
*>
*>                C = \alpha L * U
*>                      or
*>                C = \alpha U * L
*>
*> where \alpha is a scalar, L is a unit, or non-unit, lower triangular matrix, and U is a unit, or
*> non-unit, upper triangular matrix, and at most one of L and U are non-unit
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDEL
*> \verbatim
*>          SIDEL is CHARACTER*1
*>           On entry, SIDE specifies whether L multiplies U from
*>           the left or right as follows:
*>
*>             SIDE = 'L' or 'l'    A = \alpha L * U
*>
*>             SIDE = 'R' or 'r'    A = \alpha U * L
*> \endverbatim
*>
*> \param[in] DIAGL
*> \verbatim
*>          DIAGL is CHARACTER*1
*>           On entry, DIAGL specifies whether or not L is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      L is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      L is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] DIAGU
*> \verbatim
*>          DIAGU is CHARACTER*1
*>           On entry, DIAGU specifies whether or not U is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      U is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      U is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          M is INTEGER
*>           On entry, N specifies the number of rows and columns of L and U. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha. When alpha is
*>           zero then A is not referenced, and A need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, N ) where
*>           Before entry the leading n-by-n strictly upper triangular part of the array
*>           A must contain the upper triangular matrix U and the strictly lower triangular part of
*>           the leading n-by-n submatrix must contain the lower triangular matrix L.
*>           If DIAGL != 'U', then the diagonal is assumed to be part of L, and if
*>           DIAGU != 'U', then the diagonal is assumed to be part of U.
*>           Note: At most one of DIAGL and DIAGU can be not equal to 'U'.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least max( 1, n ).
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
c     Cost: 2/3 * (n^3 - n)
      RECURSIVE SUBROUTINE DLUMM(SIDEL, DIAGL, DIAGU, N, ALPHA,
     $                        A, LDA)
*
*        .. Scalar Arguments ..
         INTEGER           N, LDA
         CHARACTER         SIDEL, DIAGL, DIAGU
         DOUBLE PRECISION  ALPHA
*
*        .. Array Arguments ..
         DOUBLE PRECISION  A(LDA,*)
*        ..
*
*  =====================================================================
*
*        .. External Functions ..
         LOGICAL           LSAME
         EXTERNAL          LSAME
*        ..
*        .. External Subroutines ..
         EXTERNAL          DGEMM, DTRMM, DLASET
*        ..
*        .. Local Scalars ..
         INTEGER           K
         LOGICAL           LLEFT, LUNIT, UUNIT
*        ..
*        .. Local Parameters ..
         DOUBLE PRECISION  ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)
*        ..
*
*        Determine if our flags are valid or not. We can have at
*        most one of DIAGU, DIAGL not equal to 'U'
*
         LUNIT = LSAME(DIAGL, 'U')
         UUNIT = LSAME(DIAGU, 'U')
*
*        If both of the above are false, then it is impossible to have the
*        structure that we are exploiting in this routine
*        Note: It is possible to allow the matrices to share a non-unit
*        diagonal as long as the values are the exact same, but there is
*        currently no use case for this that I am aware of.
*
         IF ((.NOT.LUNIT).AND.(.NOT.UUNIT)) THEN
*
*        We say the error is in the last set DIAG value as we cannot know
*        what the user actually meant.
*
            CALL XERBLA( 'DLUMM', 3 )
            RETURN
         END IF
*
*        Determine which side L is on
*
         LLEFT = LSAME(SIDEL, 'L')
*
*        Early exit if possible
*
         IF (N.EQ.0) THEN
            RETURN
         END IF
         IF (ALPHA.EQ.ZERO) THEN
            CALL DLASET('All', N, N, ZERO, ZERO, A, LDA)
            RETURN
         END IF
*
*     Terminating Case
*
         IF (N.EQ.1) THEN
*
*           Since at most one of L and U are non-unit triangular, whatever side L is on, we are still
*           always computing one of
*
*           1) A(1,1) = ALPHA*A(1,1)
*           2) A(1,1) = ALPHA
*
*           Where the first case happens when exactly one of L and U are unit triangular, while the
*           second case happens when both L and U are unit triangular
*
            IF (LUNIT.AND.UUNIT) THEN
               A(1,1) = ALPHA
            ELSE
               A(1,1) = ALPHA*A(1,1)
            END IF
            RETURN
         END IF
*
*     Recursive Case
*
         K = N/2
*
*     Regardless of us computing A = L*U or A = U*L, break break A apart as follows:
*
*           |---|
*     A =   | U |
*           | L |
*           |---|
*
*     Further break down L as
*           |---------------|
*     L =   | L_{11} 0      |
*           | L_{21} L_{22} |
*           |---------------|
*
*     Where:
*
*     L_{11}\in\R^{k\times k} is lower triangular (assumed unit iff DIAGL == 'U')
*     L_{21}\in\R^{n-k\times n} is rectangular
*     L_{22}\in\R^{n-k\times n-k} is lower triangular (assumed unit iff DIAGL == 'U')
*
*     Further break down U as
*           |---------------|
*     U =   | U_{11} U_{21} |
*           | 0      U_{22} |
*           |---------------|
*
*     Where:
*
*     U_{11}\in\R^{k\times k} is upper triangular (assumed unit iff DIAGU == 'U')
*     U_{12}\in\R^{n\times n-k} is rectangular
*     U_{22}\in\R^{n-k\times n-k} is upper triangular (assumed unit iff DIAGU == 'U')
         IF (LLEFT) THEN
*
*        This means we are computing
*                         |---------------|   |---------------|
*        A = L*U = \alpha | L_{11} 0      | * | U_{11} U_{12} |
*                         | L_{21} L_{22} |   | 0      U_{22} |
*                         |---------------|   |---------------|
*
*                 |---------------------------------------------|
*        = \alpha | L_{11}*U_{11} L_{11}*U_{12}                 |
*                 | L_{21}*U_{11} L_{21}*U_{12} + L_{22}*U_{22} |
*                 |---------------------------------------------|
*
*        We compute these in the following order
*
*        A_{22} = \alpha*L_{22}*U_{22}          (This routine)
*        A_{22} = \alpha*L_{21}*U_{12} + A_{22} (GEMM)
*
*        A_{12} = \alpha*L_{11}*U_{12}          (TRMM)
*        A_{21} = \alpha*L_{21}*U_{11}          (TRMM)
*
*        A_{11} = \alpha*L_{11}*U_{11}          (This routine)
*
*        Compute A_{22}
*
*        A_{22} = \alpha*L_{22}*U_{22}
*
            CALL DLUMM(SIDEL, DIAGL, DIAGU, N-K, ALPHA,
     $               A(K+1, K+1), LDA)
*
*        A_{22} = \alpha L_{21}*U_{12} + A_{22}
*
            CALL DGEMM('No Transpose', 'No Transpose', N-K, N-K, K,
     $               ALPHA, A(K+1,1), LDA, A(1,K+1), LDA, ONE,
     $               A(K+1,K+1), LDA)
*
*        Compute A_{12}
*
*        A_{12} = \alpha*L_{11}*U_{12}
*
            CALL DTRMM('Left', 'Lower', 'No Transpose', DIAGL, K,
     $               N-K, ALPHA, A, LDA, A(1,K+1), LDA)
*
*        Compute A_{21}
*
*        A_{21} = \alpha*L_{21}*U_{11}
*
            CALL DTRMM('Right', 'Upper', 'No Transpose', DIAGU, N-K,
     $               K, ALPHA, A, LDA, A(K+1,1), LDA)
*
*        Compute A_{11}
*
*        A_{11} = \alpha*L_{11}*U_{11}
*
            CALL DLUMM(SIDEL, DIAGL, DIAGU, K, ALPHA, A, LDA)
         ELSE
*
*        This means we are computing
*                         |---------------|   |---------------|
*        A = U*L = \alpha | U_{11} U_{12} | * | L_{11} 0      |
*                         | 0      U_{22} |   | L_{21} L_{22} |
*                         |---------------|   |---------------|
*
*                 |---------------------------------------------|
*        = \alpha | U_{11}*L_{11} + U_{12}*L_{21} U_{12}*L_{22} |
*                 | U_{22}*L_{21}                 U_{22}*L_{22} |
*                 |---------------------------------------------|
*
*        We compute these in the following order
*
*        A_{11} = \alpha*U_{11}*L_{11}          (This routine)
*        A_{11} = \alpha*U_{12}*L_{21} + A_{11} (GEMM)
*
*        A_{12} = \alpha*U_{12}*L_{22}          (TRMM)
*        A_{21} = \alpha*U_{22}*L_{21}          (TRMM)
*
*        A_{22} = \alpha*U_{22}*L_{22}          (This routine)
*
*        Compute A_{11}
*
*        A_{11} = \alpha*U_{11}*L_{11}
*
            CALL DLUMM(SIDEL, DIAGL, DIAGU, K, ALPHA, A, LDA)
*
*        A_{11} = \alpha*U_{12}*L_{21} + A_{11}
*
            CALL DGEMM('No Transpose', 'No Transpose', K, K, N-K,
     $               ALPHA, A(1,K+1), LDA, A(K+1,1), LDA, ONE, A, LDA)
*
*        Compute A_{12}
*
*        A_{12} = \alpha*U_{12}*L_{22}
*
            CALL DTRMM('Right', 'Lower', 'No Transpose', DIAGL, K,
     $               N-K, ALPHA, A(K+1,K+1), LDA, A(1,K+1), LDA)
*
*        Compute A_{21}
*
*        A_{21} = \alpha*U_{22}*L_{21}
*
            CALL DTRMM('Left', 'Upper', 'No Transpose', DIAGU, N-K,
     $               K, ALPHA, A(K+1, K+1), LDA, A(K+1,1), LDA)
*
*        Compute A_{22}
*
*        A_{22} = \alpha*U_{22}*L_{22}
*
            CALL DLUMM(SIDEL, DIAGL, DIAGU, N-K, ALPHA,
     $               A(K+1, K+1), LDA)
      END IF
      END SUBROUTINE
