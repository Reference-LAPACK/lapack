      REAL             FUNCTION SMFLOP( OPS, TIME, INFO )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO
      REAL               OPS, TIME
*     ..
*
*  Purpose
*  =======
*
*  SMFLOP computes the megaflop rate given the number of operations
*  and time in seconds.  This is basically just a divide operation,
*  but care is taken not to divide by zero.
*
*  Arguments
*  =========
*
*  OPS     (input) REAL
*          The number of floating point operations.
*          performed by the timed routine.
*
*  TIME    (input) REAL
*          The total time in seconds.
*
*  INFO    (input) INTEGER
*          The return code from the timed routine.  If INFO is not 0,
*          then SMFLOP returns a negative value, indicating an error.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL
*     ..
*     .. Executable Statements ..
*
      IF( TIME.LE.ZERO ) THEN
         SMFLOP = ZERO
      ELSE
         SMFLOP = OPS / ( 1.0E6*TIME )
      END IF
      IF( INFO.NE.0 )
     $   SMFLOP = -ABS( REAL( INFO ) )
      RETURN
*
*     End of SMFLOP
*
      END
