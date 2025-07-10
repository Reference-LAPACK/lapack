*> \brief \b SAXPBY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SAXPBY(N,SA,SX,INCX,SB,SY,INCY)
*
*       .. Scalar Arguments ..
*       REAL SA,SB
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       REAL SX(*),SY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    SAXPBY constant times a vector plus constant times a vector.
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
*> \param[in] SA
*> \verbatim
*>           SA is REAL
*>           On entry, SA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] SX
*> \verbatim
*>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          storage spacing between elements of SX
*> \endverbatim
*>
*> \param[in] SB
*> \verbatim
*>           SB is REAL
*>           On entry, SB specifies the scalar beta.
*> \endverbatim
*>
*> \param[in,out] SY
*> \verbatim
*>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of SY
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
      SUBROUTINE SAXPBY(N,SA,SX,INCX,SB,SY,INCY)
      IMPLICIT NONE
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL SA,SB
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL SX(*),SY(*)
*     ..
*     .. External Subroutines ..
      EXTERNAL SSCAL
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN

*     Scale if SA.EQ.0
      IF (SA.EQ.0.0E0 .AND. SB.NE.0.0E0) THEN
          CALL SSCAL(N, SB, SY, INCY)
          RETURN
      END IF


      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO I = 1,N
            SY(I) = SB*SY(I) + SA*SX(I)
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
          SY(IY) = SB*SY(IY) + SA*SX(IX)
          IX = IX + INCX
          IY = IY + INCY
         END DO
      END IF
      RETURN
*
*     End of SAXPBY
*
      END
