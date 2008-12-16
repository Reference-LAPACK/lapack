      SUBROUTINE DSTT21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK,
     $                   RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            KBAND, LDU, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AD( * ), AE( * ), RESULT( 2 ), SD( * ),
     $                   SE( * ), U( LDU, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSTT21 checks a decomposition of the form
*
*     A = U S U'
*
*  where ' means transpose, A is symmetric tridiagonal, U is orthogonal,
*  and S is diagonal (if KBAND=0) or symmetric tridiagonal (if KBAND=1).
*  Two tests are performed:
*
*     RESULT(1) = | A - U S U' | / ( |A| n ulp )
*
*     RESULT(2) = | I - UU' | / ( n ulp )
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The size of the matrix.  If it is zero, DSTT21 does nothing.
*          It must be at least zero.
*
*  KBAND   (input) INTEGER
*          The bandwidth of the matrix S.  It may only be zero or one.
*          If zero, then S is diagonal, and SE is not referenced.  If
*          one, then S is symmetric tri-diagonal.
*
*  AD      (input) DOUBLE PRECISION array, dimension (N)
*          The diagonal of the original (unfactored) matrix A.  A is
*          assumed to be symmetric tridiagonal.
*
*  AE      (input) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal of the original (unfactored) matrix A.  A
*          is assumed to be symmetric tridiagonal.  AE(1) is the (1,2)
*          and (2,1) element, AE(2) is the (2,3) and (3,2) element, etc.
*
*  SD      (input) DOUBLE PRECISION array, dimension (N)
*          The diagonal of the (symmetric tri-) diagonal matrix S.
*
*  SE      (input) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal of the (symmetric tri-) diagonal matrix S.
*          Not referenced if KBSND=0.  If KBAND=1, then AE(1) is the
*          (1,2) and (2,1) element, SE(2) is the (2,3) and (3,2)
*          element, etc.
*
*  U       (input) DOUBLE PRECISION array, dimension (LDU, N)
*          The orthogonal matrix in the decomposition.
*
*  LDU     (input) INTEGER
*          The leading dimension of U.  LDU must be at least N.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N*(N+1))
*
*  RESULT  (output) DOUBLE PRECISION array, dimension (2)
*          The values computed by the two tests described above.  The
*          values are currently limited to 1/ulp, to avoid overflow.
*          RESULT(1) is always modified.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   ANORM, TEMP1, TEMP2, ULP, UNFL, WNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANSY
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLASET, DSYR, DSYR2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     1)      Constants
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 )
     $   RETURN
*
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Precision' )
*
*     Do Test 1
*
*     Copy A & Compute its 1-Norm:
*
      CALL DLASET( 'Full', N, N, ZERO, ZERO, WORK, N )
*
      ANORM = ZERO
      TEMP1 = ZERO
*
      DO 10 J = 1, N - 1
         WORK( ( N+1 )*( J-1 )+1 ) = AD( J )
         WORK( ( N+1 )*( J-1 )+2 ) = AE( J )
         TEMP2 = ABS( AE( J ) )
         ANORM = MAX( ANORM, ABS( AD( J ) )+TEMP1+TEMP2 )
         TEMP1 = TEMP2
   10 CONTINUE
*
      WORK( N**2 ) = AD( N )
      ANORM = MAX( ANORM, ABS( AD( N ) )+TEMP1, UNFL )
*
*     Norm of A - USU'
*
      DO 20 J = 1, N
         CALL DSYR( 'L', N, -SD( J ), U( 1, J ), 1, WORK, N )
   20 CONTINUE
*
      IF( N.GT.1 .AND. KBAND.EQ.1 ) THEN
         DO 30 J = 1, N - 1
            CALL DSYR2( 'L', N, -SE( J ), U( 1, J ), 1, U( 1, J+1 ), 1,
     $                  WORK, N )
   30    CONTINUE
      END IF
*
      WNORM = DLANSY( '1', 'L', N, WORK, N, WORK( N**2+1 ) )
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
*     Compute  UU' - I
*
      CALL DGEMM( 'N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK,
     $            N )
*
      DO 40 J = 1, N
         WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
   40 CONTINUE
*
      RESULT( 2 ) = MIN( DBLE( N ), DLANGE( '1', N, N, WORK, N,
     $              WORK( N**2+1 ) ) ) / ( N*ULP )
*
      RETURN
*
*     End of DSTT21
*
      END
