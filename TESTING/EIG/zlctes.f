*> \brief \b ZLCTES
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       LOGICAL          FUNCTION ZLCTES( Z, D )
*
*       .. Scalar Arguments ..
*       COMPLEX*16         D, Z
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLCTES returns .TRUE. if the eigenvalue Z/D is to be selected
*> (specifically, in this subroutine, if the real part of the
*> eigenvalue is negative), and otherwise it returns .FALSE..
*>
*> It is used by the test routine ZDRGES to test whether the driver
*> routine ZGGES successfully sorts eigenvalues.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] Z
*> \verbatim
*>          Z is COMPLEX*16
*>          The numerator part of a complex eigenvalue Z/D.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is COMPLEX*16
*>          The denominator part of a complex eigenvalue Z/D.
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
*> \ingroup complex16_eig
*
*  =====================================================================
      LOGICAL          FUNCTION ZLCTES( Z, D )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16         D, Z
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ZMAX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, DBLE, MAX, SIGN
*     ..
*     .. Executable Statements ..
*
      IF( D.EQ.CZERO ) THEN
         ZLCTES = ( DBLE( Z ).LT.ZERO )
      ELSE
         IF( DBLE( Z ).EQ.ZERO .OR. DBLE( D ).EQ.ZERO ) THEN
            ZLCTES = ( SIGN( ONE, AIMAG( Z ) ).NE.
     $               SIGN( ONE, AIMAG( D ) ) )
         ELSE IF( AIMAG( Z ).EQ.ZERO .OR. AIMAG( D ).EQ.ZERO ) THEN
            ZLCTES = ( SIGN( ONE, DBLE( Z ) ).NE.
     $               SIGN( ONE, DBLE( D ) ) )
         ELSE
            ZMAX = MAX( ABS( DBLE( Z ) ), ABS( AIMAG( Z ) ) )
            ZLCTES = ( ( DBLE( Z ) / ZMAX )*DBLE( D )+
     $               ( AIMAG( Z ) / ZMAX )*AIMAG( D ).LT.ZERO )
         END IF
      END IF
*
      RETURN
*
*     End of ZLCTES
*
      END
*
*  =====================================================================
*> \brief \b ZLCTSM
*>
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLCTSM returns .TRUE. if the eigenvalue Z/D is to be selected
*> (specifically, in this subroutine, if the real part of Z is
*> negative, D is nonzero (the eigenvalue is finite) and Z or D is far
*> outside the range [sqrt(safe minimum), 1/sqrt(safe minimum)] that a
*> well scaled problem lies in), and otherwise it returns .FALSE..
*>
*> Unlike ZLCTES, whose choice depends only on signs, this choice
*> depends on the magnitudes of Z and D. It is used by the test
*> routines ZDRGES and ZDRGES3 to test that ZGGES and ZGGES3 pass
*> the selection function the eigenvalues in the units of the input
*> matrices when they have scaled the matrices.
*> \endverbatim
*>
*> \param[in] Z
*> \verbatim
*>          Z is COMPLEX*16
*>          The numerator part of a complex eigenvalue Z/D.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is COMPLEX*16
*>          The denominator part of a complex eigenvalue Z/D.
*> \endverbatim
*
*> \ingroup complex16_eig
*
*  =====================================================================
      LOGICAL          FUNCTION ZLCTSM( Z, D )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16         D, Z
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   BIG, DABS1, SMALL, ZABS1
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, DBLE, SQRT
*     ..
*     .. Executable Statements ..
*
      SMALL = SQRT( DLAMCH( 'Safe minimum' ) )
      BIG = ONE / SMALL
      ZABS1 = ABS( DBLE( Z ) ) + ABS( AIMAG( Z ) )
      DABS1 = ABS( DBLE( D ) ) + ABS( AIMAG( D ) )
*
      ZLCTSM = ( DBLE( Z ).LT.ZERO ) .AND. ( DABS1.GT.ZERO ) .AND.
     $         ( ZABS1.GT.BIG .OR. DABS1.GT.BIG .OR.
     $         ( ZABS1.GT.ZERO .AND. ZABS1.LT.SMALL ) .OR.
     $         ( DABS1.GT.ZERO .AND. DABS1.LT.SMALL ) )
*
      RETURN
*
*     End of ZLCTSM
*
      END
