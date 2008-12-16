      SUBROUTINE CGBT02( TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B,
     $                   LDB, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            KL, KU, LDA, LDB, LDX, M, N, NRHS
      REAL               RESID
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  CGBT02 computes the residual for a solution of a banded system of
*  equations  A*x = b  or  A'*x = b:
*     RESID = norm( B - A*X ) / ( norm(A) * norm(X) * EPS).
*  where EPS is the machine precision.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A *x = b
*          = 'T':  A'*x = b, where A' is the transpose of A
*          = 'C':  A'*x = b, where A' is the transpose of A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  NRHS    (input) INTEGER
*          The number of columns of B.  NRHS >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The original matrix A in band storage, stored in rows 1 to
*          KL+KU+1.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,KL+KU+1).
*
*  X       (input) COMPLEX array, dimension (LDX,NRHS)
*          The computed solution vectors for the system of linear
*          equations.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  If TRANS = 'N',
*          LDX >= max(1,N); if TRANS = 'T' or 'C', LDX >= max(1,M).
*
*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
*          On entry, the right hand side vectors for the system of
*          linear equations.
*          On exit, B is overwritten with the difference B - A*X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  IF TRANS = 'N',
*          LDB >= max(1,M); if TRANS = 'T' or 'C', LDB >= max(1,N).
*
*  RESID   (output) REAL
*          The maximum over the number of right hand sides of
*          norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I1, I2, J, KD, N1
      REAL               ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SCASUM, SLAMCH
      EXTERNAL           LSAME, SCASUM, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGBMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if N = 0 pr NRHS = 0
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = SLAMCH( 'Epsilon' )
      KD = KU + 1
      ANORM = ZERO
      DO 10 J = 1, N
         I1 = MAX( KD+1-J, 1 )
         I2 = MIN( KD+M-J, KL+KD )
         ANORM = MAX( ANORM, SCASUM( I2-I1+1, A( I1, J ), 1 ) )
   10 CONTINUE
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
      IF( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) THEN
         N1 = N
      ELSE
         N1 = M
      END IF
*
*     Compute  B - A*X (or  B - A'*X )
*
      DO 20 J = 1, NRHS
         CALL CGBMV( TRANS, M, N, KL, KU, -CONE, A, LDA, X( 1, J ), 1,
     $               CONE, B( 1, J ), 1 )
   20 CONTINUE
*
*     Compute the maximum over the number of right hand sides of
*        norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
*
      RESID = ZERO
      DO 30 J = 1, NRHS
         BNORM = SCASUM( N1, B( 1, J ), 1 )
         XNORM = SCASUM( N1, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM/ANORM )/XNORM )/EPS )
         END IF
   30 CONTINUE
*
      RETURN
*
*     End of CGBT02
*
      END
