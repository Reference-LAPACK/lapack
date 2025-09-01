*> \brief \b DSECND  Using ETIME_
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
*>  This version gets the time from the system function ETIME_.
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
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
* =====================================================================
*
*     .. Local Scalars ..
      REAL               T1
*     ..
*     .. Local Arrays ..
      REAL               TARRAY( 2 )
*     ..
*     .. External Functions ..
      REAL               ETIME_
      EXTERNAL           ETIME_
*     ..
*     .. Executable Statements ..
*
      T1 = ETIME_( TARRAY )
      DSECND = TARRAY( 1 )
      RETURN
*
*     End of DSECND
*
      END
