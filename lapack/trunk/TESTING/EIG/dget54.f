      SUBROUTINE DGET54( N, A, LDA, B, LDB, S, LDS, T, LDT, U, LDU, V,
     $                   LDV, WORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDS, LDT, LDU, LDV, N
      DOUBLE PRECISION   RESULT
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( LDS, * ),
     $                   T( LDT, * ), U( LDU, * ), V( LDV, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGET54 checks a generalized decomposition of the form
*
*           A = U*S*V'  and B = U*T* V'
*
*  where ' means transpose and U and V are orthogonal.
*
*  Specifically,
*
*   RESULT = ||( A - U*S*V', B - U*T*V' )|| / (||( A, B )||*n*ulp )
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The size of the matrix.  If it is zero, DGET54 does nothing.
*          It must be at least zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA, N)
*          The original (unfactored) matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  It must be at least 1
*          and at least N.
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB, N)
*          The original (unfactored) matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  It must be at least 1
*          and at least N.
*
*  S       (input) DOUBLE PRECISION array, dimension (LDS, N)
*          The factored matrix S.
*
*  LDS     (input) INTEGER
*          The leading dimension of S.  It must be at least 1
*          and at least N.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT, N)
*          The factored matrix T.
*
*  LDT     (input) INTEGER
*          The leading dimension of T.  It must be at least 1
*          and at least N.
*
*  U       (input) DOUBLE PRECISION array, dimension (LDU, N)
*          The orthogonal matrix on the left-hand side in the
*          decomposition.
*
*  LDU     (input) INTEGER
*          The leading dimension of U.  LDU must be at least N and
*          at least 1.
*
*  V       (input) DOUBLE PRECISION array, dimension (LDV, N)
*          The orthogonal matrix on the left-hand side in the
*          decomposition.
*
*  LDV     (input) INTEGER
*          The leading dimension of V.  LDV must be at least N and
*          at least 1.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N**2)
*
*  RESULT  (output) DOUBLE PRECISION
*          The value RESULT, It is currently limited to 1/ulp, to
*          avoid overflow. Errors are flagged by RESULT=10/ulp.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ABNORM, ULP, UNFL, WNORM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      RESULT = ZERO
      IF( N.LE.0 )
     $   RETURN
*
*     Constants
*
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
*
*     compute the norm of (A,B)
*
      CALL DLACPY( 'Full', N, N, A, LDA, WORK, N )
      CALL DLACPY( 'Full', N, N, B, LDB, WORK( N*N+1 ), N )
      ABNORM = MAX( DLANGE( '1', N, 2*N, WORK, N, DUM ), UNFL )
*
*     Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)
*
      CALL DLACPY( ' ', N, N, A, LDA, WORK, N )
      CALL DGEMM( 'N', 'N', N, N, N, ONE, U, LDU, S, LDS, ZERO,
     $            WORK( N*N+1 ), N )
*
      CALL DGEMM( 'N', 'C', N, N, N, -ONE, WORK( N*N+1 ), N, V, LDV,
     $            ONE, WORK, N )
*
*     Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)
*
      CALL DLACPY( ' ', N, N, B, LDB, WORK( N*N+1 ), N )
      CALL DGEMM( 'N', 'N', N, N, N, ONE, U, LDU, T, LDT, ZERO,
     $            WORK( 2*N*N+1 ), N )
*
      CALL DGEMM( 'N', 'C', N, N, N, -ONE, WORK( 2*N*N+1 ), N, V, LDV,
     $            ONE, WORK( N*N+1 ), N )
*
*     Compute norm(W)/ ( ulp*norm((A,B)) )
*
      WNORM = DLANGE( '1', N, 2*N, WORK, N, DUM )
*
      IF( ABNORM.GT.WNORM ) THEN
         RESULT = ( WNORM / ABNORM ) / ( 2*N*ULP )
      ELSE
         IF( ABNORM.LT.ONE ) THEN
            RESULT = ( MIN( WNORM, 2*N*ABNORM ) / ABNORM ) / ( 2*N*ULP )
         ELSE
            RESULT = MIN( WNORM / ABNORM, DBLE( 2*N ) ) / ( 2*N*ULP )
         END IF
      END IF
*
      RETURN
*
*     End of DGET54
*
      END
