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
*>  This subroutine returns the BLAS API version, that is the version
*>  of LAPACK whose reference BLAS it matches.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*>  \param[out] VERS_MAJOR
*>      VERS_MAJOR is INTEGER
*>      return the BLAS API major version
*>
*>  \param[out] VERS_MINOR
*>      VERS_MINOR is INTEGER
*>      return the BLAS API minor version from the major version
*>
*>  \param[out] VERS_PATCH
*>      VERS_PATCH is INTEGER
*>      return the BLAS API patch version from the minor version
*
*  Authors:
*  ========
*
*> \author Ed J
*
*> \ingroup iblaver
*
*  =====================================================================
      SUBROUTINE IBLAVER(VERS_MAJOR, VERS_MINOR, VERS_PATCH)
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
