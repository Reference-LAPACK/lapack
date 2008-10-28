      SUBROUTINE ZHBT21( UPLO, N, KA, KS, A, LDA, D, E, U, LDU, WORK,
     $                   RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            KA, KS, LDA, LDU, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), RESULT( 2 ), RWORK( * )
      COMPLEX*16         A( LDA, * ), U( LDU, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHBT21  generally checks a decomposition of the form
*
*          A = U S U*
*
*  where * means conjugate transpose, A is hermitian banded, U is
*  unitary, and S is diagonal (if KS=0) or symmetric
*  tridiagonal (if KS=1).
*
*  Specifically:
*
*          RESULT(1) = | A - U S U* | / ( |A| n ulp ) *and*
*          RESULT(2) = | I - UU* | / ( n ulp )
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER
*          If UPLO='U', the upper triangle of A and V will be used and
*          the (strictly) lower triangle will not be referenced.
*          If UPLO='L', the lower triangle of A and V will be used and
*          the (strictly) upper triangle will not be referenced.
*
*  N       (input) INTEGER
*          The size of the matrix.  If it is zero, ZHBT21 does nothing.
*          It must be at least zero.
*
*  KA      (input) INTEGER
*          The bandwidth of the matrix A.  It must be at least zero.  If
*          it is larger than N-1, then max( 0, N-1 ) will be used.
*
*  KS      (input) INTEGER
*          The bandwidth of the matrix S.  It may only be zero or one.
*          If zero, then S is diagonal, and E is not referenced.  If
*          one, then S is symmetric tri-diagonal.
*
*  A       (input) COMPLEX*16 array, dimension (LDA, N)
*          The original (unfactored) matrix.  It is assumed to be
*          hermitian, and only the upper (UPLO='U') or only the lower
*          (UPLO='L') will be referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  It must be at least 1
*          and at least min( KA, N-1 ).
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The diagonal of the (symmetric tri-) diagonal matrix S.
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal of the (symmetric tri-) diagonal matrix S.
*          E(1) is the (1,2) and (2,1) element, E(2) is the (2,3) and
*          (3,2) element, etc.
*          Not referenced if KS=0.
*
*  U       (input) COMPLEX*16 array, dimension (LDU, N)
*          The unitary matrix in the decomposition, expressed as a
*          dense matrix (i.e., not as a product of Householder
*          transformations, Givens transformations, etc.)
*
*  LDU     (input) INTEGER
*          The leading dimension of U.  LDU must be at least N and
*          at least 1.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N**2)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (2)
*          The values computed by the two tests described above.  The
*          values are currently limited to 1/ulp, to avoid overflow.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER
      CHARACTER          CUPLO
      INTEGER            IKA, J, JC, JR
      DOUBLE PRECISION   ANORM, ULP, UNFL, WNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANHB, ZLANHP
      EXTERNAL           LSAME, DLAMCH, ZLANGE, ZLANHB, ZLANHP
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMM, ZHPR, ZHPR2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Constants
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      IKA = MAX( 0, MIN( N-1, KA ) )
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         LOWER = .FALSE.
         CUPLO = 'U'
      ELSE
         LOWER = .TRUE.
         CUPLO = 'L'
      END IF
*
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
*
*     Some Error Checks
*
*     Do Test 1
*
*     Norm of A:
*
      ANORM = MAX( ZLANHB( '1', CUPLO, N, IKA, A, LDA, RWORK ), UNFL )
*
*     Compute error matrix:    Error = A - U S U*
*
*     Copy A from SB to SP storage format.
*
      J = 0
      DO 50 JC = 1, N
         IF( LOWER ) THEN
            DO 10 JR = 1, MIN( IKA+1, N+1-JC )
               J = J + 1
               WORK( J ) = A( JR, JC )
   10       CONTINUE
            DO 20 JR = IKA + 2, N + 1 - JC
               J = J + 1
               WORK( J ) = ZERO
   20       CONTINUE
         ELSE
            DO 30 JR = IKA + 2, JC
               J = J + 1
               WORK( J ) = ZERO
   30       CONTINUE
            DO 40 JR = MIN( IKA, JC-1 ), 0, -1
               J = J + 1
               WORK( J ) = A( IKA+1-JR, JC )
   40       CONTINUE
         END IF
   50 CONTINUE
*
      DO 60 J = 1, N
         CALL ZHPR( CUPLO, N, -D( J ), U( 1, J ), 1, WORK )
   60 CONTINUE
*
      IF( N.GT.1 .AND. KS.EQ.1 ) THEN
         DO 70 J = 1, N - 1
            CALL ZHPR2( CUPLO, N, -DCMPLX( E( J ) ), U( 1, J ), 1,
     $                  U( 1, J+1 ), 1, WORK )
   70    CONTINUE
      END IF
      WNORM = ZLANHP( '1', CUPLO, N, WORK, RWORK )
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         ELSE
            RESULT( 1 ) = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
         END IF
      END IF
*
*     Do Test 2
*
*     Compute  UU* - I
*
      CALL ZGEMM( 'N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK,
     $            N )
*
      DO 80 J = 1, N
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
   80 CONTINUE
*
      RESULT( 2 ) = MIN( ZLANGE( '1', N, N, WORK, N, RWORK ),
     $              DBLE( N ) ) / ( N*ULP )
*
      RETURN
*
*     End of ZHBT21
*
      END
