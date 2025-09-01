*> \brief \b CSROT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CSROT( N, CX, INCX, CY, INCY, C, S )
*
*       .. Scalar Arguments ..
*       INTEGER           INCX, INCY, N
*       REAL              C, S
*       ..
*       .. Array Arguments ..
*       COMPLEX           CX( * ), CY( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CSROT applies a plane rotation, where the cos and sin (c and s) are real
*> and the vectors cx and cy are complex.
*> jack dongarra, linpack, 3/11/78.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the vectors cx and cy.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in,out] CX
*> \verbatim
*>          CX is COMPLEX array, dimension at least
*>           ( 1 + ( N - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array CX must contain the n
*>           element vector cx. On exit, CX is overwritten by the updated
*>           vector cx.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           CX. INCX must not be zero.
*> \endverbatim
*>
*> \param[in,out] CY
*> \verbatim
*>          CY is COMPLEX array, dimension at least
*>           ( 1 + ( N - 1 )*abs( INCY ) ).
*>           Before entry, the incremented array CY must contain the n
*>           element vector cy. On exit, CY is overwritten by the updated
*>           vector cy.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           CY. INCY must not be zero.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is REAL
*>           On entry, C specifies the cosine, cos.
*> \endverbatim
*>
*> \param[in] S
*> \verbatim
*>          S is REAL
*>           On entry, S specifies the sine, sin.
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
*> \ingroup rot
*
*  =====================================================================
      SUBROUTINE CSROT( N, CX, INCX, CY, INCY, C, S )
      IMPLICIT NONE
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           INCX, INCY, N
      REAL              C, S
*     ..
*     .. Array Arguments ..
      COMPLEX           CX( * ), CY( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER           I, IX, IY
      COMPLEX           CTEMP
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
*
*        code for both increments equal to 1
*
         DO I = 1, N
            CTEMP = C*CX( I ) + S*CY( I )
            CY( I ) = C*CY( I ) - S*CX( I )
            CX( I ) = CTEMP
         END DO
      ELSE
*
*        code for unequal increments or equal increments not equal
*          to 1
*
         IX = 1
         IY = 1
         IF( INCX.LT.0 )
     $      IX = ( -N+1 )*INCX + 1
         IF( INCY.LT.0 )
     $      IY = ( -N+1 )*INCY + 1
         DO I = 1, N
            CTEMP = C*CX( IX ) + S*CY( IY )
            CY( IY ) = C*CY( IY ) - S*CX( IX )
            CX( IX ) = CTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
*
*     End of CSROT
*
      END
