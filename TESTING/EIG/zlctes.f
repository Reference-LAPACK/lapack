      LOGICAL          FUNCTION ZLCTES( Z, D )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      COMPLEX*16         D, Z
*     ..
*
*  Purpose
*  =======
*
*  ZLCTES returns .TRUE. if the eigenvalue Z/D is to be selected
*  (specifically, in this subroutine, if the real part of the
*  eigenvalue is negative), and otherwise it returns .FALSE..
*
*  It is used by the test routine ZDRGES to test whether the driver
*  routine ZGGES succesfully sorts eigenvalues.
*
*  Arguments
*  =========
*
*  Z       (input) COMPLEX*16
*          The numerator part of a complex eigenvalue Z/D.
*
*  D       (input) COMPLEX*16
*          The denominator part of a complex eigenvalue Z/D.
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
      INTRINSIC          ABS, DBLE, DIMAG, MAX, SIGN
*     ..
*     .. Executable Statements ..
*
      IF( D.EQ.CZERO ) THEN
         ZLCTES = ( DBLE( Z ).LT.ZERO )
      ELSE
         IF( DBLE( Z ).EQ.ZERO .OR. DBLE( D ).EQ.ZERO ) THEN
            ZLCTES = ( SIGN( ONE, DIMAG( Z ) ).NE.
     $               SIGN( ONE, DIMAG( D ) ) )
         ELSE IF( DIMAG( Z ).EQ.ZERO .OR. DIMAG( D ).EQ.ZERO ) THEN
            ZLCTES = ( SIGN( ONE, DBLE( Z ) ).NE.
     $               SIGN( ONE, DBLE( D ) ) )
         ELSE
            ZMAX = MAX( ABS( DBLE( Z ) ), ABS( DIMAG( Z ) ) )
            ZLCTES = ( ( DBLE( Z ) / ZMAX )*DBLE( D )+
     $               ( DIMAG( Z ) / ZMAX )*DIMAG( D ).LT.ZERO )
         END IF
      END IF
*
      RETURN
*
*     End of ZLCTES
*
      END
