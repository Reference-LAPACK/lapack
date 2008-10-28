      SUBROUTINE ZLACRT( N, CX, INCX, CY, INCY, C, S )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      COMPLEX*16         C, S
*     ..
*     .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLACRT performs the operation
*
*     (  c  s )( x )  ==> ( x )
*     ( -s  c )( y )      ( y )
*
*  where c and s are complex and the vectors x and y are complex.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements in the vectors CX and CY.
*
*  CX      (input/output) COMPLEX*16 array, dimension (N)
*          On input, the vector x.
*          On output, CX is overwritten with c*x + s*y.
*
*  INCX    (input) INTEGER
*          The increment between successive values of CX.  INCX <> 0.
*
*  CY      (input/output) COMPLEX*16 array, dimension (N)
*          On input, the vector y.
*          On output, CY is overwritten with -s*x + c*y.
*
*  INCY    (input) INTEGER
*          The increment between successive values of CY.  INCY <> 0.
*
*  C       (input) COMPLEX*16
*  S       (input) COMPLEX*16
*          C and S define the matrix
*             [  C   S  ].
*             [ -S   C  ]
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX*16         CTEMP
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
         CTEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - S*CX( IX )
         CX( IX ) = CTEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
*
*     Code for both increments equal to 1
*
   20 CONTINUE
      DO 30 I = 1, N
         CTEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - S*CX( I )
         CX( I ) = CTEMP
   30 CONTINUE
      RETURN
      END
