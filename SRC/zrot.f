      SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   C
      COMPLEX*16         S
*     ..
*     .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
*     ..
*
*  Purpose
*  =======
*
*  ZROT   applies a plane rotation, where the cos (C) is real and the
*  sin (S) is complex, and the vectors CX and CY are complex.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements in the vectors CX and CY.
*
*  CX      (input/output) COMPLEX*16 array, dimension (N)
*          On input, the vector X.
*          On output, CX is overwritten with C*X + S*Y.
*
*  INCX    (input) INTEGER
*          The increment between successive values of CY.  INCX <> 0.
*
*  CY      (input/output) COMPLEX*16 array, dimension (N)
*          On input, the vector Y.
*          On output, CY is overwritten with -CONJG(S)*X + C*Y.
*
*  INCY    (input) INTEGER
*          The increment between successive values of CY.  INCX <> 0.
*
*  C       (input) DOUBLE PRECISION
*  S       (input) COMPLEX*16
*          C and S define a rotation
*             [  C          S  ]
*             [ -conjg(S)   C  ]
*          where C*C + S*CONJG(S) = 1.0.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX*16         STEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 )
     $   GO TO 20
*
*     Code for unequal increments or equal increments not equal to 1
*
      IX = 1
      IY = 1
      IF( INCX.LT.0 )
     $   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 )
     $   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         STEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - DCONJG( S )*CX( IX )
         CX( IX ) = STEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
*
*     Code for both increments equal to 1
*
   20 CONTINUE
      DO 30 I = 1, N
         STEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - DCONJG( S )*CX( I )
         CX( I ) = STEMP
   30 CONTINUE
      RETURN
      END
