*> \brief \b DLABAD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DLABAD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlabad.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlabad.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlabad.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLABAD( SMALL, LARGE )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   LARGE, SMALL
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLABAD is a no-op and kept for compatibility reasons. It used
*> to correct the overflow/underflow behavior of machines that
*> are not IEEE-754 compliant.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in,out] SMALL
*> \verbatim
*>          SMALL is DOUBLE PRECISION
*>          On entry, the underflow threshold as computed by DLAMCH.
*>          On exit, the unchanged value SMALL.
*> \endverbatim
*>
*> \param[in,out] LARGE
*> \verbatim
*>          LARGE is DOUBLE PRECISION
*>          On entry, the overflow threshold as computed by DLAMCH.
*>          On exit, the unchanged value LARGE.
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
*> \ingroup labad
*
*  =====================================================================
      SUBROUTINE DLABAD( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
*     ..
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
*      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
*         SMALL = SQRT( SMALL )
*         LARGE = SQRT( LARGE )
*      END IF
*
      RETURN
*
*     End of DLABAD
*
      END
