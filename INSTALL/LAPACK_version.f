*> \brief \b LAPACK_VERSION
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*      PROGRAM LAPACK_VERSION
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      PROGRAM LAPACK_VERSION
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      CHARACTER*30 LIBNAME
      INTEGER MAJOR, MINOR, PATCH, OUT_LEN
*     ..
*     .. External Subroutines ..
      EXTERNAL ILAVER, GET_LAPACK_LIBNAME
*
      CALL GET_LAPACK_LIBNAME ( LIBNAME, OUT_LEN )
      CALL ILAVER ( MAJOR, MINOR, PATCH )
      WRITE ( *, FMT = 9972 ) LIBNAME(1:OUT_LEN-1),
     $    MAJOR, MINOR, PATCH
*
 9972 FORMAT( A, ' version ', I0, '.', I0, '.', I0 )
*
      END
