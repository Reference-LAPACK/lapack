*> \brief \b SROUNDUP_LWORK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*      REAL             FUNCTION SROUNDUP_LWORK( LWORK )
*
*     .. Scalar Arguments ..
*      INTEGER          LWORK
*     ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SROUNDUP_LWORK deals with a subtle bug with returning LWORK as a Float.
*> If LWORK > 2**24, then it will get rounded as a Float.
*> This routine guarantees it is rounded up instead of down by
*> multiplying LWORK by 1+eps, where eps is the relative machine precision.
*>
*>        float( 16777217            ) == 16777216
*>        float( 16777217 * (1.+eps) ) == 16777218
*>
*> \return SROUNDUP_LWORK
*> \verbatim
*>         SROUNDUP_LWORK >= LWORK.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] LWORK Workspace size.
*
*  Authors:
*  ========
*
*> \author Weslley Pereira, University of Colorado Denver, USA
*
*> \ingroup auxOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>  This routine was inspired in the method `magma_zmake_lwork` from MAGMA.
*>  \see https://bitbucket.org/icl/magma/src/master/control/magma_zauxiliary.cpp
*> \endverbatim
*
*  =====================================================================
      REAL             FUNCTION SROUNDUP_LWORK( LWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER          LWORK
*     ..
*
* =====================================================================
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Executable Statements ..
*     ..
      SROUNDUP_LWORK = LWORK * ( 1 + SLAMCH('EPS') )
      RETURN
*
*     End of SROUNDUP_LWORK
*
      END
