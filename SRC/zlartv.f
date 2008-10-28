      SUBROUTINE ZLARTV( N, X, INCX, Y, INCY, C, S, INCC )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCC, INCX, INCY, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( * )
      COMPLEX*16         S( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARTV applies a vector of complex plane rotations with real cosines
*  to elements of the complex vectors x and y. For i = 1,2,...,n
*
*     ( x(i) ) := (        c(i)   s(i) ) ( x(i) )
*     ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of plane rotations to be applied.
*
*  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
*          The vector x.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCY)
*          The vector y.
*
*  INCY    (input) INTEGER
*          The increment between elements of Y. INCY > 0.
*
*  C       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
*          The cosines of the plane rotations.
*
*  S       (input) COMPLEX*16 array, dimension (1+(N-1)*INCC)
*          The sines of the plane rotations.
*
*  INCC    (input) INTEGER
*          The increment between elements of C and S. INCC > 0.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IX, IY
      COMPLEX*16         XI, YI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. Executable Statements ..
*
      IX = 1
      IY = 1
      IC = 1
      DO 10 I = 1, N
         XI = X( IX )
         YI = Y( IY )
         X( IX ) = C( IC )*XI + S( IC )*YI
         Y( IY ) = C( IC )*YI - DCONJG( S( IC ) )*XI
         IX = IX + INCX
         IY = IY + INCY
         IC = IC + INCC
   10 CONTINUE
      RETURN
*
*     End of ZLARTV
*
      END
