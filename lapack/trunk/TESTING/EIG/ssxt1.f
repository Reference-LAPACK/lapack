      REAL             FUNCTION SSXT1( IJOB, D1, N1, D2, N2, ABSTOL,
     $                 ULP, UNFL )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IJOB, N1, N2
      REAL               ABSTOL, ULP, UNFL
*     ..
*     .. Array Arguments ..
      REAL               D1( * ), D2( * )
*     ..
*
*  Purpose
*  =======
*
*  SSXT1  computes the difference between a set of eigenvalues.
*  The result is returned as the function value.
*
*  IJOB = 1:   Computes   max { min | D1(i)-D2(j) | }
*                          i     j
*
*  IJOB = 2:   Computes   max { min | D1(i)-D2(j) | /
*                          i     j
*                               ( ABSTOL + |D1(i)|*ULP ) }
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          Specifies the type of tests to be performed.  (See above.)
*
*  D1      (input) REAL array, dimension (N1)
*          The first array.  D1 should be in increasing order, i.e.,
*          D1(j) <= D1(j+1).
*
*  N1      (input) INTEGER
*          The length of D1.
*
*  D2      (input) REAL array, dimension (N2)
*          The second array.  D2 should be in increasing order, i.e.,
*          D2(j) <= D2(j+1).
*
*  N2      (input) INTEGER
*          The length of D2.
*
*  ABSTOL  (input) REAL
*          The absolute tolerance, used as a measure of the error.
*
*  ULP     (input) REAL
*          Machine precision.
*
*  UNFL    (input) REAL
*          The smallest positive number whose reciprocal does not
*          overflow.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               TEMP1, TEMP2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      TEMP1 = ZERO
*
      J = 1
      DO 20 I = 1, N1
   10    CONTINUE
         IF( D2( J ).LT.D1( I ) .AND. J.LT.N2 ) THEN
            J = J + 1
            GO TO 10
         END IF
         IF( J.EQ.1 ) THEN
            TEMP2 = ABS( D2( J )-D1( I ) )
            IF( IJOB.EQ.2 )
     $         TEMP2 = TEMP2 / MAX( UNFL, ABSTOL+ULP*ABS( D1( I ) ) )
         ELSE
            TEMP2 = MIN( ABS( D2( J )-D1( I ) ),
     $              ABS( D1( I )-D2( J-1 ) ) )
            IF( IJOB.EQ.2 )
     $         TEMP2 = TEMP2 / MAX( UNFL, ABSTOL+ULP*ABS( D1( I ) ) )
         END IF
         TEMP1 = MAX( TEMP1, TEMP2 )
   20 CONTINUE
*
      SSXT1 = TEMP1
      RETURN
*
*     End of SSXT1
*
      END
