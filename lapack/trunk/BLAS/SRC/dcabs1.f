      DOUBLE PRECISION FUNCTION DCABS1(Z)
*     .. Scalar Arguments ..
      COMPLEX*16 Z
*     ..
*     ..
*  Purpose
*  =======
*
*  DCABS1 computes absolute value of a double complex number 
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG
*
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN
      END
