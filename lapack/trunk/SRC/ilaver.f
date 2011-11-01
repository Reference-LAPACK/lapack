*> \brief \b ILAVER returns the LAPACK version.
*>\details
*> \b Purpose:
*>\verbatim
*>
*>  This subroutine returns the LAPACK version.
*>
*>\endverbatim
*>
*> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*> \date November 2011
*> \ingroup auxOTHERcomputational
*>
*>  \param[out] VERS_MAJOR
*>      return the lapack major version
*>
*>  \param[out] VERS_MINOR
*>      return the lapack minor version from the major version
*>
*>  \param[out] VERS_PATCH
*>      return the lapack patch version from the minor version
*>
      SUBROUTINE ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
C
C  -- LAPACK computational routine (version 3.3.1) --
C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     November 2011
C
C  =====================================================================
C
      INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH
C  =====================================================================
      VERS_MAJOR = 3
      VERS_MINOR = 3
      VERS_PATCH = 1
C  =====================================================================
C
      RETURN
      END
*> \htmlonly
*> Download ILAVER + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaver.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaver.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaver.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
