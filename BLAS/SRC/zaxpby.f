*> \brief \b ZAXPBY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZAXPBY(N,ZA,ZX,INCX,ZB,ZY,INCY)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ZA,ZB
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*),ZY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZAXPBY constant times a vector plus constant times a vector.
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
*> \param[in] ZA
*> \verbatim
*>          ZA is COMPLEX*16
*>          On entry, ZA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          storage spacing between elements of ZX
*> \endverbatim
*>
*> \param[in] ZB
*> \verbatim
*>          ZB is COMPLEX*16
*>          On entry, ZB specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] ZY
*> \verbatim
*>          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>          storage spacing between elements of ZY
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
      SUBROUTINE ZAXPBY(N,ZA,ZX,INCX,ZB,ZY,INCY)
      IMPLICIT NONE
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*16 ZA,ZB
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*     .. External Subroutines ..
      EXTERNAL ZSCAL
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN

*     Scale if ZA .EQ. 0
      IF ( ZA.EQ.(0.0D0,0.0D0) .AND. ZB.NE.(0.0D0,0.0D0)) THEN
        CALL ZSCAL(N, ZB, ZY, INCY)
        RETURN
      END IF

      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO I = 1,N
            ZY(I) = ZB*ZY(I) + ZA*ZX(I)
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
            ZY(IY) = ZB*ZY(IY) + ZA*ZX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
*
      RETURN
*
*     End of ZAXBPY
*
      END
