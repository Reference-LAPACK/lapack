      SUBROUTINE ZPPT02( UPLO, N, NRHS, A, X, LDX, B, LDB, RWORK,
     $                   RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  ZPPT02 computes the residual in the solution of a Hermitian system
*  of linear equations  A*x = b  when packed storage is used for the
*  coefficient matrix.  The ratio computed is
*
*     RESID = norm(B - A*X) / ( norm(A) * norm(X) * EPS),
*
*  where EPS is the machine precision.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of columns of B, the matrix of right hand sides.
*          NRHS >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (N*(N+1)/2)
*          The original Hermitian matrix A, stored as a packed
*          triangular matrix.
*
*  X       (input) COMPLEX*16 array, dimension (LDX,NRHS)
*          The computed solution vectors for the system of linear
*          equations.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.   LDX >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the right hand side vectors for the system of
*          linear equations.
*          On exit, B is overwritten with the difference B - A*X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RESID   (output) DOUBLE PRECISION
*          The maximum over the number of right hand sides of
*          norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANHP
      EXTERNAL           DLAMCH, DZASUM, ZLANHP
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZHPMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHP( '1', UPLO, N, A, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute  B - A*X  for the matrix of right hand sides B.
*
      DO 10 J = 1, NRHS
         CALL ZHPMV( UPLO, N, -CONE, A, X( 1, J ), 1, CONE, B( 1, J ),
     $               1 )
   10 CONTINUE
*
*     Compute the maximum over the number of right hand sides of
*        norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 20 J = 1, NRHS
         BNORM = DZASUM( N, B( 1, J ), 1 )
         XNORM = DZASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   20 CONTINUE
*
      RETURN
*
*     End of ZPPT02
*
      END
