      SUBROUTINE ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
*     
*  -- LAPACK routine (version 3.2.1)                                  --
*
*  -- April 2009                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     ..
*
*  Purpose
*  =======
*
*  This subroutine return the Lapack version.
*
*  Arguments
*  =========
*  VERS_MAJOR   (output) INTEGER
*      return the lapack major version
*  VERS_MINOR   (output) INTEGER
*      return the lapack minor version from the major version
*  VERS_PATCH   (output) INTEGER
*      return the lapack patch version from the minor version
*  =====================================================================
*
      INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH
*  =====================================================================
      VERS_MAJOR = 3
      VERS_MINOR = 2
      VERS_PATCH = 1
*  =====================================================================
*
      RETURN
      END
