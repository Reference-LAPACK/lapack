*> \brief \b DTRMMOOP_LVL2 computes an out of place triangular times general matrix multiplication
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE DTRMMOOP_LVL2(SIDE, UPLO, TRANSA, TRANSB, DIAG
*    $         M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*        .. Scalar Arguments ..
*        DOUBLE PRECISION  ALPHA, BETA
*        INTEGER           M, N, LDA, LDB, LDC
*        CHARACTER         SIDE, UPLO, TRANSA, TRANSB, DIAG
*        ..
*        .. Array Arguments ..
*        DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRMMOOP performs one of the matrix-matrix operations
*>
*>       C = \alpha op(A) * op(B) + \beta C
*>                      or
*>       C = \alpha op(B) * op(A) + \beta C
*>
*> where \alpha and \beta are scalars, C is an m-by-n matrix, A is
*> a unit, or non-unit, upper or lower triangular matrix, and op(A) is
*> is one of
*>
*>       op(A) = A      or       op(A) = A**T
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
*>             SIDE = 'L' or 'l'    C = \alpha op(A) * op(B) + \beta C
*>
*>             SIDE = 'R' or 'r'    C = \alpha op(B) * op(A) + \beta C
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>             UPLO = 'U' or 'u'    A is upper triangular
*>
*>             UPLO = 'L' or 'l'    A is lower triangular
*> \Endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op(A) to be used in
*>           the matrix multiplication as follows:
*>             TRANSA = 'N' or 'n'    op(A) = A
*>
*>             TRANSA = 'T' or 't'    op(A) = A**T
*>
*>             TRANSA = 'C' or 'c'    op(A) = A**T
*> \endverbatim
*>
*> \param[in] TRANSB
*> \verbatim
*>          TRANSB is CHARACTER*1
*>           On entry, TRANSB specifies the form of op(B) to be used in
*>           the matrix multiplication as follows:
*>             TRANSB = 'N' or 'n'     op(B) = B
*>
*>             TRANSB = 'T' or 't'     op(B) = B**T
*>
*>             TRANSB = 'C' or 'c'     op(B) = B**T
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      A is not assumed to be unit
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
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha. When alpha is
*>           zero then A and B are not referenced, and A and B need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, K ) where
*>           K is M when SIDE = 'L' and K is N when SIDE='R'
*>           Before entry with UPLO = 'U' or 'u', the leading k-by-k
*>           upper triangular part of the array A must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l', the leading k-by-k
*>           lower triangular part of the array A must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           A  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. When SIDE = 'L' or 'l' then
*>           LDA must be at least max( 1, m ), when  SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>           B is DOUBLE PRECISION array, dimension ( LDB, K ), where K is M
*>           If SIDE='R' and TRANSA='N', or SIDE='L' and TRANSA='T' and N
*>           otherwise. On entry, the leading k-by-k submatrix must contain
*>           B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in the calling (sub) program.  When  SIDE = 'R' and TRANSB='N'
*>           then LDB  must be at least  max( 1, m ), when SIDE = 'R'
*>           and TRANSB = 'T' then LDB must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta. When beta is
*>           zero then C is not referenced on entry, and C need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension ( LDC, N )
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
      SUBROUTINE DTRMMOOP_LVL2(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $         M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*        .. Scalar Arguments ..
         DOUBLE PRECISION  ALPHA, BETA
         INTEGER           M, N, LDA, LDB, LDC
         CHARACTER         SIDE, UPLO, TRANSA, TRANSB, DIAG
*        ..
*        .. Array Arguments ..
         DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*)
*        ..
*
*  =====================================================================
*
*        .. External Functions ..
         LOGICAL           LSAME
         EXTERNAL          LSAME
*        ..
*        .. External Subroutines ..
         EXTERNAL          DTRMVOOP
*        ..
*        .. Intrinsic Functions ..
         INTRINSIC         MIN
*        ..
*        .. Local Scalars ..
         INTEGER           I
         LOGICAL           LSIDE, TRANSG
         CHARACTER         TRANSF
*        ..
*
*        Beginning of Executable Statements
*
         LSIDE = LSAME(SIDE, 'L')
*
*        Determine if we are transposing the general matrix B
*
         TRANSG = LSAME(TRANSB, 'T').OR.LSAME(TRANSB, 'C')
         IF (LSIDE) THEN
*
*           This means we are computing C = \alpha op(A)*op(B) + \beta * C
*
            IF (TRANSG) THEN
*
*              op(B) = B**T
*
*              So, we are computing
*
*              C = \alpha op(A)*B**T + \beta C
*
*              Which, we compute columnwise
*
               DO I = 1, N
                  CALL DTRMVOOP(UPLO, TRANSA, DIAG, M, ALPHA,
     $                  A, LDA, B(I,1), LDB, BETA, C(1,I), 1)
               END DO
            ELSE
*
*              op(B) = B
*
*              So, we are computing
*
*              C = \alpha op(A)*B + \beta C
*
*              Which, we compute columnwise
*
               DO I = 1, N
                  CALL DTRMVOOP(UPLO, TRANSA, DIAG, M, ALPHA,
     $                  A, LDA, B(1,I), 1, BETA, C(1,I), 1)
               END DO
            END IF
         ELSE
*
*           This means we are computing C = \alpha op(B)*op(A) + \beta * C
*
            TRANSF = 'T'
            IF (LSAME(TRANSA, 'T').OR.LSAME(TRANSA,'C')) THEN
               TRANSF = 'N'
            END IF
            IF (TRANSG) THEN
*
*              op(B) = B**T
*
*              So, we are computing
*
*              C = \alpha B**T*op(A) + \beta C
*
*              Which, we compute rowwise. This means we compute the
*              I-th row of C as follows
*
*              C(I,1:N) = \alpha (B(I,1:N)**T)*op(A))**T + \beta C(I,1:N)
*                       = \alpha op(A)**T*B(I,1:N) + \beta C(I,1:N)
*
*              We transpose the main operation here as our output is
*              meant to be a vector and our kernel TRMVOOP requires A
*              to be on the left
*
               DO I = 1, M
                  CALL DTRMVOOP(UPLO, TRANSF, DIAG, N, ALPHA,
     $                  A, LDA, B(1,I), 1, BETA, C(I,1), LDC)
               END DO
            ELSE
*
*              op(B) = B
*
*              So, we are computing
*
*              C = \alpha B*op(A) + \beta C
*
*              Which, we compute rowwise. This means we compute the
*              I-th row of C as follows
*
*              C(I,1:N) = \alpha (B(I,1:N))*op(A))**T + \beta C(I,1:N)
*                       = \alpha op(A)**T*B(I,1:N)**T + \beta C(I,1:N)
*
*              We transpose the main operation here as our output is
*              meant to be a vector and our kernel TRMVOOP requires A
*              to be on the left
*
               DO I = 1, M
                  CALL DTRMVOOP(UPLO, TRANSF, DIAG, N, ALPHA,
     $                  A, LDA, B(I,1), LDB, BETA, C(I,1), LDC)
               END DO
            END IF
         END IF
      END SUBROUTINE
