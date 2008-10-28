      DOUBLE PRECISION FUNCTION DMFLOP( OPS, TIME, INFO )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO
      DOUBLE PRECISION   OPS, TIME
*     ..
*
*  Purpose
*  =======
*
*     DMFLOP computes the megaflop rate given the number of operations
*     and time in seconds.  This is basically just a divide operation,
*     but care is taken not to divide by zero.
*
*  Arguments
*  =========
*
*  OPS    - DOUBLE PRECISION
*           On entry, OPS is the number of floating point operations
*           performed by the timed routine.
*
*  TIME   - DOUBLE PRECISION
*           On entry, TIME is the total time in seconds used by the
*           timed routine.
*
*  INFO   - INTEGER
*           On entry, INFO specifies the return code from the timed
*           routine.  If INFO is not 0, then DMFLOP returns a negative
*           value, indicating an error.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE
*     ..
*     .. Executable Statements ..
*
      IF( TIME.LE.ZERO ) THEN
         DMFLOP = ZERO
      ELSE
         DMFLOP = OPS / ( 1.0D6*TIME )
      END IF
      IF( INFO.NE.0 )
     $   DMFLOP = -ABS( DBLE( INFO ) )
      RETURN
*
*     End of DMFLOP
*
      END
