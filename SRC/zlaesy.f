      SUBROUTINE ZLAESY( A, B, C, RT1, RT2, EVSCAL, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      COMPLEX*16         A, B, C, CS1, EVSCAL, RT1, RT2, SN1
*     ..
*
*  Purpose
*  =======
*
*  ZLAESY computes the eigendecomposition of a 2-by-2 symmetric matrix
*     ( ( A, B );( B, C ) )
*  provided the norm of the matrix of eigenvectors is larger than
*  some threshold value.
*
*  RT1 is the eigenvalue of larger absolute value, and RT2 of
*  smaller absolute value.  If the eigenvectors are computed, then
*  on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence
*
*  [  CS1     SN1   ] . [ A  B ] . [ CS1    -SN1   ] = [ RT1  0  ]
*  [ -SN1     CS1   ]   [ B  C ]   [ SN1     CS1   ]   [  0  RT2 ]
*
*  Arguments
*  =========
*
*  A       (input) COMPLEX*16
*          The ( 1, 1 ) element of input matrix.
*
*  B       (input) COMPLEX*16
*          The ( 1, 2 ) element of input matrix.  The ( 2, 1 ) element
*          is also given by B, since the 2-by-2 matrix is symmetric.
*
*  C       (input) COMPLEX*16
*          The ( 2, 2 ) element of input matrix.
*
*  RT1     (output) COMPLEX*16
*          The eigenvalue of larger modulus.
*
*  RT2     (output) COMPLEX*16
*          The eigenvalue of smaller modulus.
*
*  EVSCAL  (output) COMPLEX*16
*          The complex value by which the eigenvector matrix was scaled
*          to make it orthonormal.  If EVSCAL is zero, the eigenvectors
*          were not computed.  This means one of two things:  the 2-by-2
*          matrix could not be diagonalized, or the norm of the matrix
*          of eigenvectors before scaling was larger than the threshold
*          value THRESH (set below).
*
*  CS1     (output) COMPLEX*16
*  SN1     (output) COMPLEX*16
*          If EVSCAL .NE. 0,  ( CS1, SN1 ) is the unit right eigenvector
*          for RT1.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
      DOUBLE PRECISION   THRESH
      PARAMETER          ( THRESH = 0.1D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   BABS, EVNORM, TABS, Z
      COMPLEX*16         S, T, TMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*
*     Special case:  The matrix is actually diagonal.
*     To avoid divide by zero later, we treat this case separately.
*
      IF( ABS( B ).EQ.ZERO ) THEN
         RT1 = A
         RT2 = C
         IF( ABS( RT1 ).LT.ABS( RT2 ) ) THEN
            TMP = RT1
            RT1 = RT2
            RT2 = TMP
            CS1 = ZERO
            SN1 = ONE
         ELSE
            CS1 = ONE
            SN1 = ZERO
         END IF
      ELSE
*
*        Compute the eigenvalues and eigenvectors.
*        The characteristic equation is
*           lambda **2 - (A+C) lambda + (A*C - B*B)
*        and we solve it using the quadratic formula.
*
         S = ( A+C )*HALF
         T = ( A-C )*HALF
*
*        Take the square root carefully to avoid over/under flow.
*
         BABS = ABS( B )
         TABS = ABS( T )
         Z = MAX( BABS, TABS )
         IF( Z.GT.ZERO )
     $      T = Z*SQRT( ( T / Z )**2+( B / Z )**2 )
*
*        Compute the two eigenvalues.  RT1 and RT2 are exchanged
*        if necessary so that RT1 will have the greater magnitude.
*
         RT1 = S + T
         RT2 = S - T
         IF( ABS( RT1 ).LT.ABS( RT2 ) ) THEN
            TMP = RT1
            RT1 = RT2
            RT2 = TMP
         END IF
*
*        Choose CS1 = 1 and SN1 to satisfy the first equation, then
*        scale the components of this eigenvector so that the matrix
*        of eigenvectors X satisfies  X * X' = I .  (No scaling is
*        done if the norm of the eigenvalue matrix is less than THRESH.)
*
         SN1 = ( RT1-A ) / B
         TABS = ABS( SN1 )
         IF( TABS.GT.ONE ) THEN
            T = TABS*SQRT( ( ONE / TABS )**2+( SN1 / TABS )**2 )
         ELSE
            T = SQRT( CONE+SN1*SN1 )
         END IF
         EVNORM = ABS( T )
         IF( EVNORM.GE.THRESH ) THEN
            EVSCAL = CONE / T
            CS1 = EVSCAL
            SN1 = SN1*EVSCAL
         ELSE
            EVSCAL = ZERO
         END IF
      END IF
      RETURN
*
*     End of ZLAESY
*
      END
