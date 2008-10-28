      SUBROUTINE ICOPY( N, SX, INCX, SY, INCY )
*
*  -- LAPACK auxiliary test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
*     ..
*     .. Array Arguments ..
      INTEGER            SX( * ), SY( * )
*     ..
*
*  Purpose
*  =======
*
*  ICOPY copies an integer vector x to an integer vector y.
*  Uses unrolled loops for increments equal to 1.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vectors SX and SY.
*
*  SX      (input) INTEGER array, dimension (1+(N-1)*abs(INCX))
*          The vector X.
*
*  INCX    (input) INTEGER
*          The spacing between consecutive elements of SX.
*
*  SY      (output) INTEGER array, dimension (1+(N-1)*abs(INCY))
*          The vector Y.
*
*  INCY    (input) INTEGER
*          The spacing between consecutive elements of SY.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX, IY, M, MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
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
         SY( IY ) = SX( IX )
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
*
*     Code for both increments equal to 1
*
*     Clean-up loop
*
   20 CONTINUE
      M = MOD( N, 7 )
      IF( M.EQ.0 )
     $   GO TO 40
      DO 30 I = 1, M
         SY( I ) = SX( I )
   30 CONTINUE
      IF( N.LT.7 )
     $   RETURN
   40 CONTINUE
      MP1 = M + 1
      DO 50 I = MP1, N, 7
         SY( I ) = SX( I )
         SY( I+1 ) = SX( I+1 )
         SY( I+2 ) = SX( I+2 )
         SY( I+3 ) = SX( I+3 )
         SY( I+4 ) = SX( I+4 )
         SY( I+5 ) = SX( I+5 )
         SY( I+6 ) = SX( I+6 )
   50 CONTINUE
      RETURN
*
*     End of ICOPY
*
      END
