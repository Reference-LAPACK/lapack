*> \brief \b IBLAVER returns the BLAS version.
**
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE IBLAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
*
*     INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>  This subroutine returns the BLAS version.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*>  \param[out] VERS_MAJOR
*>      VERS_MAJOR is INTEGER
*>      return the BLAS major version
*>
*>  \param[out] VERS_MINOR
*>      VERS_MINOR is INTEGER
*>      return the BLAS minor version from the major version
*>
*>  \param[out] VERS_PATCH
*>      VERS_PATCH is INTEGER
*>      return the BLAS patch version from the minor version
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup iblaver
*
*  =====================================================================
      SUBROUTINE IBLAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
      IMPLICIT NONE
*
*  -- Reference BLAS auxiliary routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,--
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*  =====================================================================
*
      INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH
*  =====================================================================
      VERS_MAJOR = 3
      VERS_MINOR = 12
      VERS_PATCH = 1
*  =====================================================================
*
      RETURN
      END
