      REAL FUNCTION SCEIL( A )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     June 2008
*
*     .. Scalar Arguments ..*
      REAL A
*     ..
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
	      INTRINSIC          INT
*     ..
*     .. Executable Statements ..*
*      
      IF (A-INT(A).EQ.0) THEN
          SCEIL = A
      ELSE IF (A.GT.0) THEN
          SCEIL = INT(A)+1;
      ELSE
          SCEIL = INT(A)
      END IF

      RETURN
*
      END
