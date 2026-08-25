*> \brief \b CTRMMOOP_LVL2 computes an out of place triangular times general matrix multiplication
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE CTRMMOOP_LVL2(SIDE, UPLO, TRANSA, TRANSB, DIAG
*    $         M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*        .. Scalar Arguments ..
*        COMPLEX           ALPHA, BETA
*        INTEGER           M, N, LDA, LDB, LDC
*        CHARACTER         SIDE, UPLO, TRANSA, TRANSB, DIAG
*        ..
*        .. Array Arguments ..
*        COMPLEX           A(LDA,*), B(LDB,*), C(LDC,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CTRMMOOP performs one of the matrix-matrix operations
*>
*>       C = \alpha op(A) * op(B) + \beta C
*>                      or
*>       C = \alpha op(B) * op(A) + \beta C
*>
*> where \alpha and \beta are scalars, C is an m-by-n matrix, A is
*> a unit, or non-unit, upper or lower triangular matrix, and op(A) is
*> is one of
*>
*>          op(A) = A  or  op(A) = A**T  or  op(A) = A**H
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
*>             TRANSA = 'C' or 'c'    op(A) = A**H
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
*>             TRANSB = 'C' or 'c'     op(B) = B**H
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
*>          ALPHA is COMPLEX.
*>           On entry, ALPHA specifies the scalar alpha. When alpha is
*>           zero then A and B are not referenced, and A and B need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension ( LDA, K ) where
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
*>           B is COMPLEX array, dimension ( LDB, K ), where K is M
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
*>          BETA is COMPLEX.
*>           On entry, BETA specifies the scalar beta. When beta is
*>           zero then C is not referenced on entry, and C need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX array, dimension ( LDC, N )
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
      SUBROUTINE CTRMMOOP_LVL2(SIDE, UPLO, TRANSA, TRANSB, DIAG,
     $         M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
*
*        .. Scalar Arguments ..
         COMPLEX           ALPHA, BETA
         INTEGER           M, N, LDA, LDB, LDC
         CHARACTER         SIDE, UPLO, TRANSA, TRANSB, DIAG
*        ..
*        .. Array Arguments ..
         COMPLEX           A(LDA,*), B(LDB,*), C(LDC,*)
*        ..
*
*  =====================================================================
*
*        .. External Functions ..
         LOGICAL           LSAME
         COMPLEX           CDOTC,CDOTU
         EXTERNAL          LSAME,CDOTC,CDOTU
*        ..
*        .. External Subroutines ..
         EXTERNAL          CTRMVOOP
*        ..
*        .. Intrinsic Functions ..
         INTRINSIC         MIN, CONJG
*        ..
*        .. Local Scalars ..
         INTEGER           I, J
         LOGICAL           LSIDE, TRANSG, CONJB, CONJA, TRANST, UNIT,
     $                     UPPER
         CHARACTER         TRANSF
*        ..
*
*        Beginning of Executable Statements
*
         LSIDE = LSAME(SIDE, 'L')
*
*        Determine if we are transposing the general matrix B
*
         CONJB = LSAME(TRANSB,'C')
         CONJA = LSAME(TRANSA,'C')

         TRANSG = CONJB.OR.LSAME(TRANSB, 'T')
         TRANST = CONJA.OR.LSAME(TRANSA, 'T')

         UNIT = LSAME(DIAG, 'U')
         UPPER = LSAME(UPLO,'U')
         IF (LSIDE) THEN
*
*           This means we are computing C = \alpha op(A)*op(B) + \beta * C
*
            IF (CONJB) THEN
*
*              op(B) = B**H
*
*              So, we are computing
*
*              C = \alpha op(A)*B**H + \beta C
*
*              Which, we compute columnwise
*
               DO I = 1, N
                  CALL CTRMCVOOP(UPLO, TRANSA, DIAG, M, ALPHA,
     $                  A, LDA, B(I,1), LDB, BETA, C(1,I), 1)
               END DO
            ELSE IF (TRANSG) THEN
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
                  CALL CTRMVOOP(UPLO, TRANSA, DIAG, M, ALPHA,
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
                  CALL CTRMVOOP(UPLO, TRANSA, DIAG, M, ALPHA,
     $                  A, LDA, B(1,I), 1, BETA, C(1,I), 1)
               END DO
            END IF
         ELSE
*
*           This means we are computing C = \alpha op(B)*op(A) + \beta * C
*
            IF (CONJB) THEN
*
*              op(B) = B**H
*
*              So, we are computing
*
*              C = \alpha B**H*op(A) + \beta C
*
*              Which we compute componentwise.
*              TODO: Consider extracting this into another subroutine
               IF (CONJA) THEN
                  IF (UPPER) THEN
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CONJG(B(J,I)) +
     $                           ALPHA*CONJG(CDOTU(N-J, B(J+1,I), 1,
     $                           A(J,J+1), LDA)) + BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CONJG(CDOTU(N-J+1,
     $                           B(J,I), 1, A(J,J), LDA)) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  ELSE
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CONJG(B(J,I)) +
     $                           ALPHA*CONJG(CDOTU(J-1, B(1,I), 1,
     $                            A(J,1), LDA)) + BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CONJG(CDOTU(J,
     $                           B(1,I), 1, A(J,1), LDA)) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  END IF
               ELSE IF (TRANST) THEN
                  ! op(A) = A**T
                  IF (UPPER) THEN
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CONJG(B(J,I)) +
     $                           ALPHA*CDOTC(N-J, B(J+1,I), 1,
     $                           A(J,J+1), LDA) + BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CDOTC(N-J+1,
     $                           B(J,I), 1, A(J,J), LDA) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  ELSE
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CONJG(B(J,I)) +
     $                           ALPHA*CDOTC(J-1, B(1,I), 1,
     $                            A(J,1), LDA) + BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CDOTC(J, B(1,I), 1,
     $                           A(J,1), LDA) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  END IF
               ELSE! op(A) = A
                  IF (UPPER) THEN
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CONJG(B(J,I)) +
     $                           ALPHA*CDOTC(J-1, B(1,I), 1,
     $                           A(1,J), 1) + BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CDOTC(J,
     $                           B(1,I), 1, A(1,J), 1) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  ELSE ! A is lower
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CONJG(B(J,I)) +
     $                           ALPHA*CDOTC(N-J, B(J+1,I), 1,
     $                            A(J+1,J), 1) + BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CDOTC(N-J+1,
     $                           B(J,I), 1, A(J,J), 1) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  END IF
               END IF
            ELSE IF (TRANSG) THEN
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
               IF (CONJA) THEN
                  IF (UPPER) THEN
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*B(J,I) +
     $                           ALPHA*CDOTC(N-J, A(J,J+1), LDA,
     $                           B(J+1,I), 1) + BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CDOTC(N-J+1,
     $                           A(J,J), LDA, B(J,I), 1) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  ELSE
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*B(J,I) +
     $                           ALPHA*CDOTC(J-1, A(J,1), LDA,
     $                            B(1,I), 1) + BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CDOTC(J,
     $                           A(J,1), LDA, B(1,I), 1) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  END IF
               ELSE
                  TRANSF = 'T'
                  IF (TRANST) THEN
                     TRANSF = 'N'
                  END IF
                  DO I = 1, M
                     CALL CTRMVOOP(UPLO, TRANSF, DIAG, N, ALPHA,
     $                     A, LDA, B(1,I), 1, BETA, C(I,1), LDC)
                  END DO
               END IF
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
               IF (CONJA) THEN
                  IF (UPPER) THEN
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*B(I,J) +
     $                           ALPHA*CDOTC(N-J, A(J,J+1), LDA,
     $                           B(I,J+1), LDB) + BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CDOTC(N-J+1,
     $                           A(J,J), LDA, B(I,J), LDB) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  ELSE
                     IF (UNIT) THEN
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*B(I,J) +
     $                           ALPHA*CDOTC(J-1, A(J,1), LDA,
     $                           B(I,1), LDB) +  BETA*C(I,J)
                           END DO
                        END DO
                     ELSE
                        DO J = 1, N
                           DO I = 1, M
                              C(I,J) = ALPHA*CDOTC(J, A(J,1), LDA,
     $                           B(I,1), LDB) + BETA*C(I,J)
                           END DO
                        END DO
                     END IF
                  END IF
               ELSE
                  TRANSF = 'T'
                  IF (TRANST) THEN
                     TRANSF = 'N'
                  END IF
                  DO I = 1, M
                     CALL CTRMVOOP(UPLO, TRANSF, DIAG, N, ALPHA,
     $                     A, LDA, B(I,1), LDB, BETA, C(I,1), LDC)
                  END DO
               END IF
            END IF
         END IF
      END SUBROUTINE
