      SUBROUTINE ZLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS1, RT1, RT2
      COMPLEX*16         A, B, C, SN1
*     ..
*
*  Purpose
*  =======
*
*  ZLAEV2 computes the eigendecomposition of a 2-by-2 Hermitian matrix
*     [  A         B  ]
*     [  CONJG(B)  C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*  eigenvector for RT1, giving the decomposition
*
*  [ CS1  CONJG(SN1) ] [    A     B ] [ CS1 -CONJG(SN1) ] = [ RT1  0  ]
*  [-SN1     CS1     ] [ CONJG(B) C ] [ SN1     CS1     ]   [  0  RT2 ].
*
*  Arguments
*  =========
*
*  A      (input) COMPLEX*16
*         The (1,1) element of the 2-by-2 matrix.
*
*  B      (input) COMPLEX*16
*         The (1,2) element and the conjugate of the (2,1) element of
*         the 2-by-2 matrix.
*
*  C      (input) COMPLEX*16
*         The (2,2) element of the 2-by-2 matrix.
*
*  RT1    (output) DOUBLE PRECISION
*         The eigenvalue of larger absolute value.
*
*  RT2    (output) DOUBLE PRECISION
*         The eigenvalue of smaller absolute value.
*
*  CS1    (output) DOUBLE PRECISION
*  SN1    (output) COMPLEX*16
*         The vector (CS1, SN1) is a unit right eigenvector for RT1.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  CS1 and SN1 are accurate to a few ulps barring over/underflow.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   T
      COMPLEX*16         W
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAEV2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG
*     ..
*     .. Executable Statements ..
*
      IF( ABS( B ).EQ.ZERO ) THEN
         W = ONE
      ELSE
         W = DCONJG( B ) / ABS( B )
      END IF
      CALL DLAEV2( DBLE( A ), ABS( B ), DBLE( C ), RT1, RT2, CS1, T )
      SN1 = W*T
      RETURN
*
*     End of ZLAEV2
*
      END
