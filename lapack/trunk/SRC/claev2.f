      SUBROUTINE CLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      REAL               CS1, RT1, RT2
      COMPLEX            A, B, C, SN1
*     ..
*
*  Purpose
*  =======
*
*  CLAEV2 computes the eigendecomposition of a 2-by-2 Hermitian matrix
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
*  A      (input) COMPLEX
*         The (1,1) element of the 2-by-2 matrix.
*
*  B      (input) COMPLEX
*         The (1,2) element and the conjugate of the (2,1) element of
*         the 2-by-2 matrix.
*
*  C      (input) COMPLEX
*         The (2,2) element of the 2-by-2 matrix.
*
*  RT1    (output) REAL
*         The eigenvalue of larger absolute value.
*
*  RT2    (output) REAL
*         The eigenvalue of smaller absolute value.
*
*  CS1    (output) REAL
*  SN1    (output) COMPLEX
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
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               T
      COMPLEX            W
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAEV2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, REAL
*     ..
*     .. Executable Statements ..
*
      IF( ABS( B ).EQ.ZERO ) THEN
         W = ONE
      ELSE
         W = CONJG( B ) / ABS( B )
      END IF
      CALL SLAEV2( REAL( A ), ABS( B ), REAL( C ), RT1, RT2, CS1, T )
      SN1 = W*T
      RETURN
*
*     End of CLAEV2
*
      END
