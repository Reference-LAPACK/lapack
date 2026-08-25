*> \brief \b CTRSM_LVL2_MOD solves a singular triangular system
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CTRSM_LVL2_MOD(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*       .. Scalar Arguments ..
*       COMPLEX ALPHA
*       INTEGER LDA,LDB,M,N
*       CHARACTER DIAG,SIDE,TRANSA,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX A(LDA,*),B(LDB,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CTRSM_LVL2_MOD  solves one of the matrix equations
*>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T.
*>
*> The matrix X is overwritten on B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether op( A ) appears on the left
*>           or right of X as follows:
*>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n'   op( A ) = A.
*>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c'   op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX.
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*>           zero then  A is not referenced and  B need not be set before
*>           entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension ( LDA, k ),
*>           where k is m when SIDE = 'L' or 'l'
*>             and k is n when SIDE = 'R' or 'r'.
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*>           upper triangular part of the array  A must contain the upper
*>           triangular matrix  and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*>           lower triangular part of the array  A must contain the lower
*>           triangular matrix  and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           A  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension ( LDB, N )
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain  the  right-hand  side  matrix  B,  and  on exit  is
*>           overwritten by the solution matrix  X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
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
*> \ingroup trsm
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*> Modified from ctrsm originally written in 1989
*> This routine solves singular triangular systems with explicit
*> 0s stored only on the diagonal of A on input. We assume relevant
*> parts of B are 0 without checking the entry values
*>
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CTRSM_LVL2_MOD(SIDE, UPLO, TRANSA, DIAG, M, N,
     $      ALPHA, A, LDA, B, LDB)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      COMPLEX ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX A(LDA,*),B(LDB,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      COMPLEX TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOUNIT,UPPER,NOTRAN,CONJA
*     ..
*     .. Parameters ..
      COMPLEX ONE,ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
         INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
         INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
         INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND.
     +         (.NOT.LSAME(DIAG,'N'))) THEN
         INFO = 4
      ELSE IF (M.LT.0) THEN
         INFO = 5
      ELSE IF (N.LT.0) THEN
         INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
         INFO = 11
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('CTRSM ',INFO)
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN

         CALL CLASET('All', M, N, ZERO, ZERO, B, LDB)
         RETURN
      END IF
*
*     Start the operations.
*
      LSIDE = LSAME(SIDE,'L')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
      NOTRAN = LSAME(TRANSA,'N')
      CONJA = LSAME(TRANSA,'C')
      ! Pre-scale by alpha if necessary
      IF (ALPHA.NE.ONE) THEN
         DO J = 1,N
            CALL CSCAL(M, ALPHA, B(1,J), 1)
         END DO
      END IF

      IF (LSIDE) THEN
         IF (NOTRAN) THEN
            IF (UPPER) THEN
               DO J = 1,N
                  DO K = M,1,-1
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(K,J) = ZERO
                        ELSE
                           B(K,J) = B(K,J)/A(K,K)
                           TEMP = B(K,J)
                           DO I = 1,K-1
                              B(I,J) = B(I,J)-TEMP*A(I,K)
                           END DO
                        END IF
                     ELSE
                        TEMP = B(K,J)
                        DO I = 1,K-1
                           B(I,J) = B(I,J)-TEMP*A(I,K)
                        END DO
                     END IF
                  END DO
               END DO
            ELSE ! A is lower triangular
               DO J = 1,N
                  DO K = 1,M
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(K,J) = ZERO
                        ELSE
                           B(K,J) = B(K,J)/A(K,K)
                           TEMP = B(K,J)

                           DO I = K+1,M
                              B(I,J) = B(I,J)-TEMP*A(I,K)
                           END DO
                        END IF
                     ELSE
                        TEMP = B(K,J)
                        DO I = K+1,M
                           B(I,J) = B(I,J)-TEMP*A(I,K)
                        END DO
                     END IF
                  END DO
               END DO
            END IF
         ELSE IF( CONJA ) THEN
            IF (UPPER) THEN
               DO J = 1,N
                  DO K = 1,M
                     TEMP = B(K,J)
                     DO I = K+1,M
                        B(I,J) = B(I,J)-TEMP*CONJG(A(K,I))
                     END DO
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(K,J) = ZERO
                        ELSE
                           B(K,J) = B(K,J)/CONJG(A(K,K))
                        END IF
                     END IF
                  END DO
               END DO
            ELSE ! A is lower triangular
               DO J = 1,N
                  DO K = M,1,-1
                     TEMP = B(K,J)
                     DO I = 1,K-1
                        B(I,J) = B(I,J)-TEMP*CONJG(A(K,I))
                     END DO
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(K,J) = ZERO
                        ELSE
                           B(K,J) = B(K,J)/CONJG(A(K,K))
                        END IF
                     END IF
                  END DO
               END DO
            END IF
         ELSE ! A is transposed
            IF (UPPER) THEN
               DO J = 1,N
                  DO K = 1,M
                     TEMP = B(K,J)
                     DO I = K+1,M
                        B(I,J) = B(I,J)-TEMP*A(K,I)
                     END DO
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(K,J) = ZERO
                        ELSE
                           B(K,J) = B(K,J)/A(K,K)
                        END IF
                     END IF
                  END DO
               END DO
            ELSE ! A is lower triangular
               DO J = 1,N
                  DO K = M,1,-1
                     TEMP = B(K,J)
                     DO I = 1,K-1
                        B(I,J) = B(I,J)-TEMP*A(K,I)
                     END DO
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(K,J) = ZERO
                        ELSE
                           B(K,J) = B(K,J)/A(K,K)
                        END IF
                     END IF
                  END DO
               END DO
            END IF
         END IF
      ELSE ! A is on the right
         IF (NOTRAN) THEN
            IF (UPPER) THEN
               DO I = 1,M
                  DO K = N,1,-1
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(I,K) = ZERO
                        ELSE
                           B(I,K) = B(I,K)/A(K,K)
                           TEMP = B(I,K)
                           DO J = 1,K-1
                              B(I,J) = B(I,J)-TEMP*A(J,K)
                           END DO
                        END IF
                     ELSE
                        TEMP = B(I,K)
                        DO J = 1,K-1
                           B(I,J) = B(I,J)-TEMP*A(J,K)
                        END DO
                     END IF
                  END DO
               END DO
            ELSE
               DO I = 1,M
                  DO K = 1,N
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(I,K) = ZERO
                        ELSE
                           B(I,K) = B(I,K)/A(K,K)
                           TEMP = B(I,K)
                           DO J = K+1,N
                              B(I,J) = B(I,J)-TEMP*A(J,K)
                           END DO
                        END IF
                     ELSE
                        TEMP = B(I,K)
                        DO J = K+1,N
                           B(I,J) = B(I,J)-TEMP*A(J,K)
                        END DO
                     END IF
                  END DO
               END DO
            END IF
         ELSE IF( CONJA ) THEN
            IF (UPPER) THEN
               DO I = 1,M
                  DO K = 1,N
                     TEMP = B(I,K)
                     DO J = K+1,N
                        B(I,J) = B(I,J)-TEMP*CONJG(A(K,J))
                     END DO
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(I,K) = ZERO
                        ELSE
                           B(I,K) = B(I,K)/CONJG(A(K,K))
                        END IF
                     END IF
                  END DO
               END DO
            ELSE
               DO I = 1,M
                  DO K = N,1,-1
                     TEMP = B(I,K)
                     DO J = 1,K-1
                        B(I,J) = B(I,J)-TEMP*CONJG(A(K,J))
                     END DO
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(I,K) = ZERO
                        ELSE
                           B(I,K) = B(I,K)/CONJG(A(K,K))
                        END IF
                     END IF
                  END DO
               END DO
            END IF
         ELSE ! A is transposed
            IF (UPPER) THEN
               DO I = 1,M
                  DO K = 1,N
                     TEMP = B(I,K)
                     DO J = K+1,N
                        B(I,J) = B(I,J)-TEMP*A(K,J)
                     END DO
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(I,K) = ZERO
                        ELSE
                           B(I,K) = B(I,K)/A(K,K)
                        END IF
                     END IF
                  END DO
               END DO
            ELSE
               DO I = 1,M
                  DO K = N,1,-1
                     TEMP = B(I,K)
                     DO J = 1,K-1
                        B(I,J) = B(I,J)-TEMP*A(K,J)
                     END DO
                     IF (NOUNIT) THEN
                        IF (A(K,K).EQ.ZERO) THEN
                           B(I,K) = ZERO
                        ELSE
                           B(I,K) = B(I,K)/A(K,K)
                        END IF
                     END IF
                  END DO
               END DO
            END IF
         END IF
      END IF
      RETURN
*
*     End of CTRSM
*
      END
