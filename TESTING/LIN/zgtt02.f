      SUBROUTINE ZGTT02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB,
     $                   RWORK, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ),
     $                   X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGTT02 computes the residual for the solution to a tridiagonal
*  system of equations:
*     RESID = norm(B - op(A)*X) / (norm(A) * norm(X) * EPS),
*  where EPS is the machine epsilon.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER
*          Specifies the form of the residual.
*          = 'N':  B - A * X     (No transpose)
*          = 'T':  B - A**T * X  (Transpose)
*          = 'C':  B - A**H * X  (Conjugate transpose)
*
*  N       (input) INTEGTER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices B and X.  NRHS >= 0.
*
*  DL      (input) COMPLEX*16 array, dimension (N-1)
*          The (n-1) sub-diagonal elements of A.
*
*  D       (input) COMPLEX*16 array, dimension (N)
*          The diagonal elements of A.
*
*  DU      (input) COMPLEX*16 array, dimension (N-1)
*          The (n-1) super-diagonal elements of A.
*
*  X       (input) COMPLEX*16 array, dimension (LDX,NRHS)
*          The computed solution vectors X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the right hand side vectors for the system of
*          linear equations.
*          On exit, B is overwritten with the difference B - op(A)*X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RESID   (output) DOUBLE PRECISION
*          norm(B - op(A)*X) / (norm(A) * norm(X) * EPS)
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANGT
      EXTERNAL           LSAME, DLAMCH, DZASUM, ZLANGT
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLAGTM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0
*
      RESID = ZERO
      IF( N.LE.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     Compute the maximum over the number of right hand sides of
*        norm(B - op(A)*X) / ( norm(A) * norm(X) * EPS ).
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         ANORM = ZLANGT( '1', N, DL, D, DU )
      ELSE
         ANORM = ZLANGT( 'I', N, DL, D, DU )
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute B - op(A)*X.
*
      CALL ZLAGTM( TRANS, N, NRHS, -ONE, DL, D, DU, X, LDX, ONE, B,
     $             LDB )
*
      DO 10 J = 1, NRHS
         BNORM = DZASUM( N, B( 1, J ), 1 )
         XNORM = DZASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of ZGTT02
*
      END
