      SUBROUTINE SROTG(SA,SB,C,S)
*     .. Scalar Arguments ..
      REAL C,S,SA,SB
*     ..
*
*  Purpose
*  =======
*
*     SROTG construct givens plane rotation.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL R,ROE,SCALE,Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,SQRT
*     ..
      ROE = SB
      IF (ABS(SA).GT.ABS(SB)) ROE = SA
      SCALE = ABS(SA) + ABS(SB)
      IF (SCALE.EQ.0.0) THEN
         C = 1.0
         S = 0.0
         R = 0.0
         Z = 0.0
      ELSE
         R = SCALE*SQRT((SA/SCALE)**2+ (SB/SCALE)**2)
         R = SIGN(1.0,ROE)*R
         C = SA/R
         S = SB/R
         Z = 1.0
         IF (ABS(SA).GT.ABS(SB)) Z = S
         IF (ABS(SB).GE.ABS(SA) .AND. C.NE.0.0) Z = 1.0/C
      END IF
      SA = R
      SB = Z
      RETURN
      END
