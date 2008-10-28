      LOGICAL          FUNCTION DLCTES( ZR, ZI, D )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   D, ZI, ZR
*     ..
*
*  Purpose
*  =======
*
*  DLCTES returns .TRUE. if the eigenvalue (ZR/D) + sqrt(-1)*(ZI/D)
*  is to be selected (specifically, in this subroutine, if the real
*  part of the eigenvalue is negative), and otherwise it returns
*  .FALSE..
*
*  It is used by the test routine DDRGES to test whether the driver
*  routine DGGES succesfully sorts eigenvalues.
*
*  Arguments
*  =========
*
*  ZR      (input) DOUBLE PRECISION
*          The numerator of the real part of a complex eigenvalue
*          (ZR/D) + i*(ZI/D).
*
*  ZI      (input) DOUBLE PRECISION
*          The numerator of the imaginary part of a complex eigenvalue
*          (ZR/D) + i*(ZI).
*
*  D       (input) DOUBLE PRECISION
*          The denominator part of a complex eigenvalue
*          (ZR/D) + i*(ZI/D).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SIGN
*     ..
*     .. Executable Statements ..
*
      IF( D.EQ.ZERO ) THEN
         DLCTES = ( ZR.LT.ZERO )
      ELSE
         DLCTES = ( SIGN( ONE, ZR ).NE.SIGN( ONE, D ) )
      END IF
*
      RETURN
*
*     End of DLCTES
*
      END
