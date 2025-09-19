*> \brief \b CAXPBY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CAXPBY(N,CA,CX,INCX,CB,CY,INCY)
*
*       .. Scalar Arguments ..
*       COMPLEX CA,CB
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX CX(*),CY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    CAXPBY constant times a vector plus constant times a vector.
*>
*>    Y = ALPHA * X + BETA * Y
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] CA
*> \verbatim
*>          CA is COMPLEX
*>          On entry, CA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] CX
*> \verbatim
*>          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          storage spacing between elements of CX
*> \endverbatim
*>
*> \param[in] CB
*> \verbatim
*>          CB is COMPLEX
*>          On entry, CB specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] CY
*> \verbatim
*>          CY is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>          storage spacing between elements of CY
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*> \author Martin Koehler, MPI Magdeburg
*
*> \ingroup axpby
*
*  =====================================================================
      SUBROUTINE CAXPBY(N,CA,CX,INCX,CB,CY,INCY)
      IMPLICIT NONE
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX CA, CB
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX CX(*),CY(*)
*     ..
*     .. External Subroutines ..
      EXTERNAL CSCAL
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN

      IF (CA .EQ. (0.0,0.0) .AND. CB.NE.(0.0,0.0)) THEN
         CALL CSCAL(N,CB, CY, INCY)
         RETURN
      END IF

      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO I = 1,N
            CY(I) = CB*CY(I) + CA*CX(I)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CY(IY) = CB*CY(IY) + CA*CX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
*
      RETURN
*
*     End of CAXBPY
*
      END
