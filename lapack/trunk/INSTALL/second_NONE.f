      REAL             FUNCTION SECOND( )
*
*  -- LAPACK auxiliary routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     February 2007
*
*  Purpose
*  =======
*
*  SECOND returns nothing instead of returning the user time for a process in seconds.
*  If you are using that routine, it means that neither EXTERNAL ETIME,
*  EXTERNAL ETIME_, INTERNAL ETIME, INTERNAL CPU_TIME is available  on
*  your machine.
*
* =====================================================================
*
      SECOND = 0.0E+0
      RETURN
*
*     End of SECOND
*
      END
