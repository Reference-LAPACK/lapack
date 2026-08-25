*> \brief \b DTRMVOOP
*
*  =========== DOCUMENTATION ===========
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTRMVOOP( UPLO, TRANS, DIAG, N, ALPHA, A, LDA,
*    $               X, INCX, BETA, Y, INCY )
*
*     .. Scalar Arguments ..
*     INTEGER           N, LDA, INCX, INCY
*     CHARACTER         UPLO, TRANS, DIAG
*     DOUBLE PRECISION  ALPHA, BETA
*     ..
*     .. Array Arguments ..
*     DOUBLE PRECISION  A(LDA,*),X(*),Y(*)
*     ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTRMVOOP  performs one of the matrix-vector operations
*>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are n element vectors, and
*> A is an n by n unit, or non-unit, upper or lower triangular matrix.
*>
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
*>              TRANS = 'C' or 'c'   y := A**T*x + y.
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
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension ( LDA, N )
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
*>          X is DOUBLE PRECISION array, dimension at least
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
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array, dimension at least
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
        SUBROUTINE DTRMVOOP( UPLO, TRANS, DIAG, N, ALPHA, A, LDA,
     $               X, INCX, BETA, Y, INCY )
*
*  -- Reference LAPACK level2 routine --
*  -- Reference LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           N, LDA, INCX, INCY
      CHARACTER         UPLO, TRANS, DIAG
      DOUBLE PRECISION  ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER (ZERO=0.0D+0, ONE=1.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER           INFO,KX,KY,IY,IX,I
      DOUBLE PRECISION  TEMP
      LOGICAL           NOUNIT,UPPER
*     ..
*     .. External Functions ..
      DOUBLE PRECISION  DDOT
      LOGICAL           LSAME
      EXTERNAL          LSAME, DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL          XERBLA, DAXPY, DSCAL
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
          CALL XERBLA('DTRMVOOP',INFO)
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
         ! Did not find a good way to test this, but not needed by me. Can be
         ! implemented later if needed
         !KX = 1 - (N-1)*INCX
         CALL XERBLA('DTRMVOOP',9)
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
         ! Did not find a good way to test this, but not needed by me. Can be
         ! implemented later if needed
         !KY = 1 - (N-1)*INCY
         CALL XERBLA('DTRMVOOP',12)
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ZERO.AND.BETA.NE.ONE) THEN
         CALL DSCAL(N, BETA, Y, INCY)
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      ! This implementation computes y one component at a time
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
         IF (UPPER) THEN
            IF (NOUNIT) THEN
               ! A is NOT assumed unit triangular, so we use A(I,I)
               IF (BETA.EQ.ZERO) THEN
                  ! We overwrite y with alpha*A*x
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * DDOT(N-I+1, A(I,I), LDA,
     $                           X(I), INCX)
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * DDOT(N-I+1, A(I,I), LDA,
     $                           X(I), INCX)
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * DDOT(N-I+1, A(I,I), LDA,
     $                           X(IX), INCX)
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * DDOT(N-I+1, A(I,I), LDA,
     $                           X(IX), INCX)
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               ELSE
                  ! We add y to alpha*A*x
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * DDOT(N-I+1, A(I,I), LDA,
     $                           X(I), INCX) + Y(I)
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * DDOT(N-I+1, A(I,I), LDA,
     $                           X(I), INCX) + Y(IY)
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * DDOT(N-I+1, A(I,I), LDA,
     $                           X(IX), INCX) + Y(I)
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * DDOT(N-I+1, A(I,I), LDA,
     $                           X(IX), INCX) + Y(IY)
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               END IF
            ELSE
               ! A is assumed unit triangular, so we do NOT use A(I,I)
               IF (BETA.EQ.ZERO) THEN
                  ! We overwrite y with alpha*A*x
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(N-I, A(I,I+1), LDA,
     $                           X(I+1), INCX) + X(I))
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(N-I, A(I,I+1), LDA,
     $                           X(I+1), INCX) + X(I))
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(N-I, A(I,I+1), LDA,
     $                           X(IX+INCX), INCX) + X(IX))
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(N-I, A(I,I+1), LDA,
     $                           X(IX+INCX), INCX) + X(IX))
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               ELSE
                  ! We add y to alpha*A*x
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(N-I, A(I,I+1), LDA,
     $                           X(I+1), INCX) + X(I)) + Y(I)
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(N-I, A(I,I+1), LDA,
     $                           X(I+1), INCX) + X(I)) + Y(IY)
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(N-I, A(I,I+1), LDA,
     $                           X(IX+INCX), INCX) + X(IX)) + Y(I)
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(N-I, A(I,I+1), LDA,
     $                           X(IX+INCX), INCX) + X(IX)) + Y(IY)
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               END IF
            END IF
         ELSE ! A is lower triangular
            IF (NOUNIT) THEN
               IF (BETA.EQ.ZERO) THEN
                  ! Replace y with alpha*A*x
                  IF (INCY.EQ.1) THEN
                     DO I = 1, N
                        Y(I) = ALPHA * DDOT(I, A(I,1), LDA,
     $                        X, INCX)
                     END DO
                  ELSE
                     IY = KY
                     DO I = 1, N
                        Y(IY) = ALPHA * DDOT(I, A(I,1), LDA,
     $                        X, INCX)
                        IY = IY + INCY
                     END DO
                  END IF
               ELSE
                  ! add alpha*A*x to y
                  IF (INCY.EQ.1) THEN
                     DO I = 1, N
                        Y(I) = ALPHA * DDOT(I, A(I,1), LDA,
     $                        X, INCX) + Y(I)
                     END DO
                  ELSE
                     IY = KY
                     DO I = 1, N
                        Y(IY) = ALPHA * DDOT(I, A(I,1), LDA,
     $                        X, INCX) + Y(IY)
                        IY = IY + INCY
                     END DO
                  END IF
               END IF
            ELSE ! A is assumed unit
               IF (BETA.EQ.ZERO) THEN
                  ! Replace y with alpha*A*x
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(I-1, A(I,1), LDA,
     $                           X, INCX) + X(I))
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(I-1, A(I,1), LDA,
     $                           X, INCX) + X(I))
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(I-1, A(I,1), LDA,
     $                           X, INCX) + X(IX))
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(I-1, A(I,1), LDA,
     $                           X, INCX) + X(IX))
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               ELSE
                  ! add alpha*A*x to y
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(I-1, A(I,1), LDA,
     $                           X, INCX) + X(I)) + Y(I)
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(I-1, A(I,1), LDA,
     $                           X, INCX) + X(I)) + Y(IY)
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(I-1, A(I,1), LDA,
     $                           X, INCX) + X(IX)) + Y(I)
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(I-1, A(I,1), LDA,
     $                           X, INCX) + X(IX)) + Y(IY)
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               END IF
            END IF
         END IF
      ELSE
*
*        Form  y := alpha*A**T*x + y.
*
         IF (UPPER) THEN
            IF (NOUNIT) THEN
               IF (BETA.EQ.ZERO) THEN
                  ! replace y with alpha A**T x
                  IF (INCY.EQ.1) THEN
                     DO I = 1, N
                        Y(I) = ALPHA * DDOT(I, A(1,I), 1, X, INCX)
                     END DO
                  ELSE
                     IY = KY
                     DO I = 1, N
                        Y(IY) = ALPHA * DDOT(I, A(1,I), 1, X, INCX)
                        IY = IY + INCY
                     END DO
                  END IF
               ELSE
                   ! add alpha A**T x to y
                  IF (INCY.EQ.1) THEN
                     DO I = 1, N
                        Y(I) = ALPHA * DDOT(I, A(1,I), 1, X, INCX)
     $                        + Y(I)
                     END DO
                  ELSE
                     IY = KY
                     DO I = 1, N
                        Y(IY) = ALPHA * DDOT(I, A(1,I), 1, X, INCX)
     $                        + Y(IY)
                        IY = IY + INCY
                     END DO
                  END IF
               END IF
            ELSE ! A is assumed unit triangular
               IF (BETA.EQ.ZERO) THEN
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(I-1, A(1,I), 1,
     $                           X, INCX) + X(I))
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(I-1, A(1,I), 1,
     $                           X, INCX) + X(I))
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(I-1, A(1,I), 1,
     $                           X, INCX) + X(IX))
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(I-1, A(1,I), 1,
     $                           X, INCX) + X(IX))
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               ELSE
                  ! add alpha*A**T*x to y
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(I-1, A(1,I), 1,
     $                           X, INCX) + X(I)) + Y(I)
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(I-1, A(1,I), 1,
     $                           X, INCX) + X(I)) + Y(IY)
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(I-1, A(1,I), 1,
     $                           X, INCX) + X(IX)) + Y(I)
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(I-1, A(1,I), 1,
     $                           X, INCX) + X(IX)) + Y(IY)
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               END IF
            END IF
         ELSE ! A is lower triangular
            IF (NOUNIT) THEN
               IF (BETA.EQ.ZERO) THEN
                  ! replace y with alpha*A**T *x
                  ! assume incx,incy.ne.1
                  IX = KX
                  IY = KY
                  DO I = 1, N
                     Y(IY) = ALPHA * DDOT(N-I+1, A(I,I), 1,
     $                     X(IX), INCX)
                     IX = IX + INCX
                     IY = IY + INCY
                  END DO
               ELSE
                  ! add alpha*A**T *x to y
                  ! assume incx,incy.ne.1
                  IX = KX
                  IY = KY
                  DO I = 1, N
                     Y(IY) = ALPHA * DDOT(N-I+1, A(I,I), 1,
     $                     X(IX), INCX) + Y(IY)
                     IX = IX + INCX
                     IY = IY + INCY
                  END DO
               END IF
            ELSE
               IF (BETA.EQ.ZERO) THEN
                  ! replace y with alpha*A**T *x
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(N-I, A(I+1,I), 1,
     $                           X(I+1), 1) + X(I))
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(N-I, A(I+1,I), 1,
     $                           X(I+1), 1) + X(I))
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(N-I, A(I+1,I), 1,
     $                           X(IX+INCX), INCX) + X(IX))
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(N-I, A(I+1,I), 1,
     $                           X(IX+INCX), INCX) + X(IX))
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               ELSE
                  ! add alpha*A**T *x to y
                  IF (INCX.EQ.1) THEN
                     IF (INCY.EQ.1) THEN
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(N-I, A(I+1,I), 1,
     $                           X(I+1), 1) + X(I)) + Y(I)
                        END DO
                     ELSE
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(N-I, A(I+1,I), 1,
     $                           X(I+1), 1) + X(I)) + Y(IY)
                           IY = IY + INCY
                        END DO
                     END IF
                  ELSE
                     IF (INCY.EQ.1) THEN
                        IX = KX
                        DO I = 1, N
                           Y(I) = ALPHA * (DDOT(N-I, A(I+1,I), 1,
     $                           X(IX+INCX), INCX) + X(IX)) + Y(I)
                           IX = IX + INCX
                        END DO
                     ELSE
                        IX = KX
                        IY = KY
                        DO I = 1, N
                           Y(IY) = ALPHA * (DDOT(N-I, A(I+1,I), 1,
     $                           X(IX+INCX), INCX) + X(IX)) + Y(IY)
                           IX = IX + INCX
                           IY = IY + INCY
                        END DO
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
      END SUBROUTINE
