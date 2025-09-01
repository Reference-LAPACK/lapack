*> \brief \b DLASUM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASUM( TYPE, IOUNIT, IE, NRUN )
*
*       .. Scalar Arguments ..
*       CHARACTER*3        TYPE
*       INTEGER            IE, IOUNIT, NRUN
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLASUM prints a summary of the results from one of the test routines.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE DLASUM( TYPE, IOUNIT, IE, NRUN )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*3        TYPE
      INTEGER            IE, IOUNIT, NRUN
*     ..
*
*  =====================================================================
*
*     .. Executable Statements ..
*
      IF( IE.GT.0 ) THEN
         WRITE( IOUNIT, FMT = 9999 )TYPE, ': ', IE, ' out of ', NRUN,
     $      ' tests failed to pass the threshold'
      ELSE
         WRITE( IOUNIT, FMT = 9998 )'All tests for ', TYPE,
     $      ' passed the threshold ( ', NRUN, ' tests run)'
      END IF
 9999 FORMAT( 1X, A3, A2, I4, A8, I5, A35 )
 9998 FORMAT( / 1X, A14, A3, A24, I5, A11 )
      RETURN
*
*     End of DLASUM
*
      END
