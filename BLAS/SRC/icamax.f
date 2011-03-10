      INTEGER FUNCTION ICAMAX(N,CX,INCX)
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
*     ICAMAX finds the index of element having max. absolute value.
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
      REAL SMAX
      INTEGER I,IX
*     ..
*     .. External Functions ..
      REAL SCABS1
      EXTERNAL SCABS1
*     ..
      ICAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ICAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         SMAX = SCABS1(CX(1))
         DO I = 2,N
            IF (SCABS1(CX(I)).GT.SMAX) THEN
               ICAMAX = I
               SMAX = SCABS1(CX(I))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         IX = 1
         SMAX = SCABS1(CX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (SCABS1(CX(IX)).GT.SMAX) THEN
               ICAMAX = I
               SMAX = SCABS1(CX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
