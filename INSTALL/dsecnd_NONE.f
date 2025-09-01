*> \brief \b DSECND returns nothing
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
*>  DSECND returns nothing instead of returning the user time for a process in seconds.
*>  If you are using that routine, it means that neither EXTERNAL ETIME,
*>  EXTERNAL ETIME_, INTERNAL ETIME, INTERNAL CPU_TIME is available  on
*>  your machine.
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
      DSECND = 0.0D+0
      RETURN
*
*     End of DSECND
*
      END
