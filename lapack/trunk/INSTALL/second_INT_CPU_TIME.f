      REAL             FUNCTION SECOND( )
*
*  -- LAPACK auxiliary routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     February 2007
*
*  Purpose
*  =======
*
*  SECOND returns the user time for a process in seconds.
*  This version gets the time from the INTERNAL function CPU_TIME.
*
* =====================================================================
*
*     .. Local Scalars ..
* 
      REAL T
* 
* .. Intrinsic Functions ..
* 
      INTRINSIC CPU_TIME
* 
* .. Executable Statements .. *
* 
      CALL CPU_TIME( T )
      SECOND = T
      RETURN
*
*     End of SECOND
*
      END
