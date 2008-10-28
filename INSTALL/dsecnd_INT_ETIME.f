      DOUBLE PRECISION FUNCTION DSECND( )
*
*  -- LAPACK auxiliary routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     February 2007
*
*  Purpose
*  =======
*
*  DSECND returns the user time for a process in seconds.
*  This version gets the time from the INTERNAL function ETIME.
*
* =====================================================================
*
*     .. Local Scalars ..
      REAL               T1
*     ..
*     .. Local Arrays ..
      REAL               TARRAY( 2 )
*     ..
*     .. Intrinsic Functions ..
      REAL               ETIME
      INTRINSIC          ETIME
*     ..
*     .. Executable Statements ..
*
      T1 = ETIME( TARRAY )
      DSECND = TARRAY( 1 )
      RETURN
*
*     End of DSECND
*
      END
