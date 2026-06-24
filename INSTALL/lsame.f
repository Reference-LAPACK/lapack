*> \brief \b LSAME
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*      LOGICAL FUNCTION LSAME( CA, CB )
*
*     .. Scalar Arguments ..
*      CHARACTER          CA, CB
*     ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> LSAME returns .TRUE. if CA is the same letter as CB regardless of
*> case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CA
*> \param[in] CB
*> \verbatim
*>          CA and CB specify the single characters to be compared.
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
*> \ingroup lsame
*
*  =====================================================================
      LOGICAL FUNCTION LSAME( CA, CB )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          IACHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB
*     ..
*     .. Executable Statements ..
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Use IACHAR to get ASCII codes for case-insensitive comparison,
*     regardless of the platform's native character encoding.
*
      INTA = IACHAR( CA )
      INTB = IACHAR( CB )
      IF( INTA.GE.65 .AND. INTA.LE.90 ) INTA = INTA + 32
      IF( INTB.GE.65 .AND. INTB.LE.90 ) INTB = INTB + 32
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
