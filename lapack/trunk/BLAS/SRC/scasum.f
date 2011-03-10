      REAL FUNCTION SCASUM(N,CX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX CX(*)
*     ..
*
*  Purpose
*  =======
*
*     SCASUM takes the sum of the absolute values of a complex vector and
*     returns a single precision result.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL STEMP
      INTEGER I,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,REAL
*     ..
      SCASUM = 0.0e0
      STEMP = 0.0e0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO I = 1,N
            STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
         END DO
      END IF
      SCASUM = STEMP
      RETURN
      END
