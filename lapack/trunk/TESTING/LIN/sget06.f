      REAL             FUNCTION SGET06( RCOND, RCONDC )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      REAL               RCOND, RCONDC
*     ..
*
*  Purpose
*  =======
*
*  SGET06 computes a test ratio to compare two values for RCOND.
*
*  Arguments
*  ==========
*
*  RCOND   (input) REAL
*          The estimate of the reciprocal of the condition number of A,
*          as computed by SGECON.
*
*  RCONDC  (input) REAL
*          The reciprocal of the condition number of A, computed as
*          ( 1/norm(A) ) / norm(inv(A)).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      REAL               EPS, RAT
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
      IF( RCOND.GT.ZERO ) THEN
         IF( RCONDC.GT.ZERO ) THEN
            RAT = MAX( RCOND, RCONDC ) / MIN( RCOND, RCONDC ) -
     $            ( ONE-EPS )
         ELSE
            RAT = RCOND / EPS
         END IF
      ELSE
         IF( RCONDC.GT.ZERO ) THEN
            RAT = RCONDC / EPS
         ELSE
            RAT = ZERO
         END IF
      END IF
      SGET06 = RAT
      RETURN
*
*     End of SGET06
*
      END
