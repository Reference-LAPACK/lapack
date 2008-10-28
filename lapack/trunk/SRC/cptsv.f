      SUBROUTINE CPTSV( N, NRHS, D, E, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL               D( * )
      COMPLEX            B( LDB, * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  CPTSV computes the solution to a complex system of linear equations
*  A*X = B, where A is an N-by-N Hermitian positive definite tridiagonal
*  matrix, and X and B are N-by-NRHS matrices.
*
*  A is factored as A = L*D*L**H, and the factored form of A is then
*  used to solve the system of equations.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix
*          A.  On exit, the n diagonal elements of the diagonal matrix
*          D from the factorization A = L*D*L**H.
*
*  E       (input/output) COMPLEX array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix A.  On exit, the (n-1) subdiagonal elements of the
*          unit bidiagonal factor L from the L*D*L**H factorization of
*          A.  E can also be regarded as the superdiagonal of the unit
*          bidiagonal factor U from the U**H*D*U factorization of A.
*
*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the solution has not been
*                computed.  The factorization has not been completed
*                unless i = N.
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           CPTTRF, CPTTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPTSV ', -INFO )
         RETURN
      END IF
*
*     Compute the L*D*L' (or U'*D*U) factorization of A.
*
      CALL CPTTRF( N, D, E, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL CPTTRS( 'Lower', N, NRHS, D, E, B, LDB, INFO )
      END IF
      RETURN
*
*     End of CPTSV
*
      END
