      SUBROUTINE CBDT03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK,
     $                   RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            KD, LDU, LDVT, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), S( * )
      COMPLEX            U( LDU, * ), VT( LDVT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CBDT03 reconstructs a bidiagonal matrix B from its SVD:
*     S = U' * B * V
*  where U and V are orthogonal matrices and S is diagonal.
*
*  The test ratio to test the singular value decomposition is
*     RESID = norm( B - U * S * VT ) / ( n * norm(B) * EPS )
*  where VT = V' and EPS is the machine precision.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix B is upper or lower bidiagonal.
*          = 'U':  Upper bidiagonal
*          = 'L':  Lower bidiagonal
*
*  N       (input) INTEGER
*          The order of the matrix B.
*
*  KD      (input) INTEGER
*          The bandwidth of the bidiagonal matrix B.  If KD = 1, the
*          matrix B is bidiagonal, and if KD = 0, B is diagonal and E is
*          not referenced.  If KD is greater than 1, it is assumed to be
*          1, and if KD is less than 0, it is assumed to be 0.
*
*  D       (input) REAL array, dimension (N)
*          The n diagonal elements of the bidiagonal matrix B.
*
*  E       (input) REAL array, dimension (N-1)
*          The (n-1) superdiagonal elements of the bidiagonal matrix B
*          if UPLO = 'U', or the (n-1) subdiagonal elements of B if
*          UPLO = 'L'.
*
*  U       (input) COMPLEX array, dimension (LDU,N)
*          The n by n orthogonal matrix U in the reduction B = U'*A*P.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= max(1,N)
*
*  S       (input) REAL array, dimension (N)
*          The singular values from the SVD of B, sorted in decreasing
*          order.
*
*  VT      (input) COMPLEX array, dimension (LDVT,N)
*          The n by n orthogonal matrix V' in the reduction
*          B = U * S * V'.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT.
*
*  WORK    (workspace) COMPLEX array, dimension (2*N)
*
*  RESID   (output) REAL
*          The test ratio:  norm(B - U * S * V') / ( n * norm(A) * EPS )
*
* ======================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               BNORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ISAMAX
      REAL               SCASUM, SLAMCH
      EXTERNAL           LSAME, ISAMAX, SCASUM, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      RESID = ZERO
      IF( N.LE.0 )
     $   RETURN
*
*     Compute B - U * S * V' one column at a time.
*
      BNORM = ZERO
      IF( KD.GE.1 ) THEN
*
*        B is bidiagonal.
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           B is upper bidiagonal.
*
            DO 20 J = 1, N
               DO 10 I = 1, N
                  WORK( N+I ) = S( I )*VT( I, J )
   10          CONTINUE
               CALL CGEMV( 'No transpose', N, N, -CMPLX( ONE ), U, LDU,
     $                     WORK( N+1 ), 1, CMPLX( ZERO ), WORK, 1 )
               WORK( J ) = WORK( J ) + D( J )
               IF( J.GT.1 ) THEN
                  WORK( J-1 ) = WORK( J-1 ) + E( J-1 )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J-1 ) ) )
               ELSE
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               END IF
               RESID = MAX( RESID, SCASUM( N, WORK, 1 ) )
   20       CONTINUE
         ELSE
*
*           B is lower bidiagonal.
*
            DO 40 J = 1, N
               DO 30 I = 1, N
                  WORK( N+I ) = S( I )*VT( I, J )
   30          CONTINUE
               CALL CGEMV( 'No transpose', N, N, -CMPLX( ONE ), U, LDU,
     $                     WORK( N+1 ), 1, CMPLX( ZERO ), WORK, 1 )
               WORK( J ) = WORK( J ) + D( J )
               IF( J.LT.N ) THEN
                  WORK( J+1 ) = WORK( J+1 ) + E( J )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J ) ) )
               ELSE
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               END IF
               RESID = MAX( RESID, SCASUM( N, WORK, 1 ) )
   40       CONTINUE
         END IF
      ELSE
*
*        B is diagonal.
*
         DO 60 J = 1, N
            DO 50 I = 1, N
               WORK( N+I ) = S( I )*VT( I, J )
   50       CONTINUE
            CALL CGEMV( 'No transpose', N, N, -CMPLX( ONE ), U, LDU,
     $                  WORK( N+1 ), 1, CMPLX( ZERO ), WORK, 1 )
            WORK( J ) = WORK( J ) + D( J )
            RESID = MAX( RESID, SCASUM( N, WORK, 1 ) )
   60    CONTINUE
         J = ISAMAX( N, D, 1 )
         BNORM = ABS( D( J ) )
      END IF
*
*     Compute norm(B - U * S * V') / ( n * norm(B) * EPS )
*
      EPS = SLAMCH( 'Precision' )
*
      IF( BNORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         IF( BNORM.GE.RESID ) THEN
            RESID = ( RESID / BNORM ) / ( REAL( N )*EPS )
         ELSE
            IF( BNORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REAL( N )*BNORM ) / BNORM ) /
     $                 ( REAL( N )*EPS )
            ELSE
               RESID = MIN( RESID / BNORM, REAL( N ) ) /
     $                 ( REAL( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of CBDT03
*
      END
