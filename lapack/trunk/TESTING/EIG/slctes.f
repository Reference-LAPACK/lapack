      LOGICAL          FUNCTION SLCTES( ZR, ZI, D )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      REAL               D, ZI, ZR
*     ..
*
*  Purpose
*  =======
*
*  SLCTES returns .TRUE. if the eigenvalue (ZR/D) + sqrt(-1)*(ZI/D)
*  is to be selected (specifically, in this subroutine, if the real
*  part of the eigenvalue is negative), and otherwise it returns
*  .FALSE..
*
*  It is used by the test routine SDRGES to test whether the driver
*  routine SGGES succesfully sorts eigenvalues.
*
*  Arguments
*  =========
*
*  ZR      (input) REAL
*          The numerator of the real part of a complex eigenvalue
*          (ZR/D) + i*(ZI/D).
*
*  ZI      (input) REAL
*          The numerator of the imaginary part of a complex eigenvalue
*          (ZR/D) + i*(ZI).
*
*  D       (input) REAL
*          The denominator part of a complex eigenvalue
*          (ZR/D) + i*(ZI/D).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SIGN
*     ..
*     .. Executable Statements ..
*
      IF( D.EQ.ZERO ) THEN
         SLCTES = ( ZR.LT.ZERO )
      ELSE
         SLCTES = ( SIGN( ONE, ZR ).NE.SIGN( ONE, D ) )
      END IF
*
      RETURN
*
*     End of SLCTES
*
      END
