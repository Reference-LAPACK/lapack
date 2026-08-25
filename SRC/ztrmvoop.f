*> \brief \b ZTRMVOOP
*
*  =========== DOCUMENTATION ===========
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRMVOOP( UPLO, TRANS, DIAG, N, ALPHA, A, LDA,
*    $               X, INCX, BETA, Y, INCY )
*
*     .. Scalar Arguments ..
*     INTEGER           N, LDA, INCX, INCY
*     CHARACTER         UPLO, TRANS, DIAG
*     COMPLEX*16        ALPHA, BETA
*     ..
*     .. Array Arguments ..
*     COMPLEX*16        A(LDA,*),X(*),Y(*)
*     ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRMVOOP  performs one of the matrix-vector operations
*>
*>                y := alpha*op(A)*x + beta*y,
*>
*> where op(A) = A, A**T, or A**H, alpha and beta are scalars, x and y
*> are n element vectors, and A is an n by n unit, or non-unit, upper
*> or lower triangular matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   y := A*x + y.
*>
*>              TRANS = 'T' or 't'   y := A**T*x + y.
*>
*>              TRANS = 'C' or 'c'   y := A**H*x + y.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit
*>           triangular as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, N )
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*>           A are not referenced either, but are assumed to be unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must be positive.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16.
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is COMPLEX*16 array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ).
*>           Before entry with BETA non-zero, the incremented array Y
*>           must contain the vector y. On exit, Y is overwritten by the
*>           updated vector y.
*>           If n is zero, then Y is not referenced and the function
*>           performs a quick return.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must be positive.
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
*> \ingroup trmv
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0
*> \endverbatim
*>
*  =====================================================================
        SUBROUTINE ZTRMVOOP( UPLO, TRANS, DIAG, N, ALPHA, A, LDA,
     $               X, INCX, BETA, Y, INCY )
*
*  -- Reference LAPACK level2 routine --
*  -- Reference LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           N, LDA, INCX, INCY
      CHARACTER         UPLO, TRANS, DIAG
      COMPLEX*16        ALPHA, BETA
*     ..
*     .. Array Arguments ..
      COMPLEX*16        A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16        ZERO
      PARAMETER (ZERO=(0.0D+0,0.0D+0))
*     ..
*     .. Local Scalars ..
      INTEGER           INFO,KX,KY,IY,IX,I
      LOGICAL           NOUNIT,UPPER
*     ..
*     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL          XERBLA, ZACXPY, ZAXPY, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC         MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     $    .NOT.LSAME(TRANS,'C')) THEN
         INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND.
     $         .NOT.LSAME(DIAG,'N')) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = 7
      ELSE IF (INCX.EQ.0) THEN
         INFO = 9
      ELSE IF (INCY.EQ.0) THEN
         INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZTRMVOOP',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
*     Determine the starting point for X and Y if their increments are not 1
*
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
         ! Did not find a good way to test, but not needed by me. Can be
         ! implemented later if needed
         !KX = 1 - (N-1)*INCX
         CALL XERBLA('ZTRMVOOP',9)
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
         ! Did not find a good way to test, but not needed by me. Can be
         ! implemented later if needed
         !KY = 1 - (N-1)*INCY
         CALL XERBLA('ZTRMVOOP',12)
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (BETA.EQ.ZERO) THEN
         IF (INCY.EQ.1) THEN
            DO I = 1, N
               Y(I) = ZERO
            END DO
         ELSE
            IY = KY
            DO I = 1, N
               Y(IY) = ZERO
               IY = IY + INCY
            END DO
         END IF
      ELSE
*
*        A sane implementation will do nothing on BETA.EQ.DCMPLX(1.0D+0)
*
         CALL ZSCAL(N, BETA, Y, INCY)
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
         IF (UPPER) THEN
            IF (INCX.EQ.1) THEN
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(1:I,I) + y
*
                     CALL ZAXPY(I, ALPHA*X(I), A(1,I), 1, Y, INCY)
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(1:I-1,I) + y
*
                     CALL ZAXPY(I-1, ALPHA*X(I), A(1,I), 1, Y, INCY)
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            ELSE
               IX = KX
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(1:I,I) + y
*
                     CALL ZAXPY(I, ALPHA*X(IX), A(1,I), 1, Y, INCY)
                     IX = IX + INCX
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(1:I-1,I) + y
*
                     CALL ZAXPY(I-1, ALPHA*X(IX), A(1,I), 1, Y, INCY)
                     IX = IX + INCX
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            END IF
         ELSE
            IF (INCX.EQ.1) THEN
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I:N,I) + y
*
                     CALL ZAXPY(N-I+1, ALPHA*X(I), A(I,I), 1,
     $                     Y(KY + (I-1)*INCY), INCY)
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I+1:N,I) + y
*
                     CALL ZAXPY(N-I, ALPHA*X(I), A(I+1,I), 1,
     $                     Y(KY + I*INCY), INCY)
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            ELSE
               IX = KX
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I:N,I) + y
*
                     CALL ZAXPY(N-I+1, ALPHA*X(IX), A(I,I), 1,
     $                     Y(KY + (I-1)*INCY), INCY)
                     IX = IX + INCX
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I+1:N,I) + y
*
                     CALL ZAXPY(N-I, ALPHA*X(IX), A(I+1,I), 1,
     $                     Y(KY + I*INCY), INCY)
                     IX = IX + INCX
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            END IF
         END IF
      ELSE IF (LSAME(TRANS,'T')) THEN
*
*        Form  y := alpha*A**T*x + y.
*
         IF (UPPER) THEN
            IF (INCX.EQ.1) THEN
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I,I:N) + y
*
                     CALL ZAXPY(N-I+1, ALPHA*X(I), A(I,I), LDA,
     $                     Y(KY + (I-1)*INCY), INCY)
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I,I+1:N) + y
*
                     CALL ZAXPY(N-I, ALPHA*X(I), A(I,I+1), LDA,
     $                     Y(KY + I*INCY), INCY)
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            ELSE
               IX = KX
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I,I:N) + y
*
                     CALL ZAXPY(N-I+1, ALPHA*X(IX), A(I,I), LDA,
     $                     Y(KY + (I-1)*INCY), INCY)
                     IX = IX + INCX
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I,I+1:N) + y
*
                     CALL ZAXPY(N-I, ALPHA*X(IX), A(I,I+1), LDA,
     $                     Y(KY + I*INCY), INCY)
                     IX = IX + INCX
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            END IF
         ELSE
            IF (INCX.EQ.1) THEN
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I,1:I) + y
*
                     CALL ZAXPY(I, ALPHA*X(I), A(I,1), LDA,
     $                     Y, INCY)
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I,1:I-1) + y
*
                     CALL ZAXPY(I-1, ALPHA*X(I), A(I,1), LDA,
     $                     Y, INCY)
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            ELSE
               IX = KX
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I,1:I) + y
*
                     CALL ZAXPY(I, ALPHA*X(IX), A(I,1), LDA,
     $                     Y, INCY)
                     IX = IX + INCX
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*A(I,1:I-1) + y
*
                     CALL ZAXPY(I-1, ALPHA*X(IX), A(I,1), LDA,
     $                     Y, INCY)
                     IX = IX + INCX
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            END IF
         END IF
      ELSE !IF (LSAME(TRANS,'C')) THEN
*
*        Form  y := alpha*A**H*x + y.
*
         IF (UPPER) THEN
            IF (INCX.EQ.1) THEN
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*DCONJG(A(I,I:N)) + y
*
                     CALL ZACXPY(N-I+1, ALPHA*X(I), A(I,I), LDA,
     $                     Y(KY + (I-1)*INCY), INCY)
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*DCONJG(A(I,I+1:N)) + y
*
                     CALL ZACXPY(N-I, ALPHA*X(I), A(I,I+1), LDA,
     $                     Y(KY + I*INCY), INCY)
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            ELSE
               IX = KX
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*DCONJG(A(I,I:N)) + y
*
                     CALL ZACXPY(N-I+1, ALPHA*X(IX), A(I,I), LDA,
     $                     Y(KY + (I-1)*INCY), INCY)
                     IX = IX + INCX
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*DCONJG(A(I,I+1:N)) + y
*
                     CALL ZACXPY(N-I, ALPHA*X(IX), A(I,I+1), LDA,
     $                     Y(KY + I*INCY), INCY)
                     IX = IX + INCX
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            END IF
         ELSE
            IF (INCX.EQ.1) THEN
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*DCONJG(A(I,1:I)) + y
*
                     CALL ZACXPY(I, ALPHA*X(I), A(I,1), LDA,
     $                     Y, INCY)
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*DCONJG(A(I,1:I-1)) + y
*
                     CALL ZACXPY(I-1, ALPHA*X(I), A(I,1), LDA,
     $                     Y, INCY)
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            ELSE
               IX = KX
               IF (NOUNIT) THEN
                  DO I = 1, N
*
*                    y = (alpha*x(i))*DCONJG(A(I,1:I)) + y
*
                     CALL ZACXPY(I, ALPHA*X(IX), A(I,1), LDA,
     $                     Y, INCY)
                     IX = IX + INCX
                  END DO
               ELSE
                  DO I = 1, N
*
*                    y = (alpha*x(i))*DCONJG(A(I,1:I-1)) + y
*
                     CALL ZACXPY(I-1, ALPHA*X(IX), A(I,1), LDA,
     $                     Y, INCY)
                     IX = IX + INCX
                  END DO
*
*                 y = alpha*x + y
*
                  CALL ZAXPY(N, ALPHA, X, INCX, Y, INCY)
               END IF
            END IF
         END IF
      END IF
      END SUBROUTINE
