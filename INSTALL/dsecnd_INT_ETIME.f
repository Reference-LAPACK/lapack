*> \brief \b DSECND Using the INTERNAL function ETIME.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*      DOUBLE PRECISION FUNCTION DSECND( )
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>  DSECND returns the user time for a process in seconds.
*>  This version gets the time from the INTERNAL function ETIME.
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
*> \ingroup second
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DSECND( )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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
