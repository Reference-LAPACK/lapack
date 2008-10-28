      SUBROUTINE ZGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U,
     $                   WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, N, P
      DOUBLE PRECISION   RESULT
*     ..
*     .. Array Arguments ..
*
*  Purpose
*  =======
*
*  ZGLMTS tests ZGGGLM - a subroutine for solving the generalized
*  linear model problem.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of rows of the matrices A and B.  N >= 0.
*
*  M       (input) INTEGER
*          The number of columns of the matrix A.  M >= 0.
*
*  P       (input) INTEGER
*          The number of columns of the matrix B.  P >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,M)
*          The N-by-M matrix A.
*
*  AF      (workspace) COMPLEX*16 array, dimension (LDA,M)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF. LDA >= max(M,N).
*
*  B       (input) COMPLEX*16 array, dimension (LDB,P)
*          The N-by-P matrix A.
*
*  BF      (workspace) COMPLEX*16 array, dimension (LDB,P)
*
*  LDB     (input) INTEGER
*          The leading dimension of the arrays B, BF. LDB >= max(P,N).
*
*  D       (input) COMPLEX*16 array, dimension( N )
*          On input, the left hand side of the GLM.
*
*  DF      (workspace) COMPLEX*16 array, dimension( N )
*
*  X       (output) COMPLEX*16 array, dimension( M )
*          solution vector X in the GLM problem.
*
*  U       (output) COMPLEX*16 array, dimension( P )
*          solution vector U in the GLM problem.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (M)
*
*  RESULT   (output) DOUBLE PRECISION
*          The test ratio:
*                           norm( d - A*x - B*u )
*            RESULT = -----------------------------------------
*                     (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
*
*  ====================================================================
*
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ),
     $                   BF( LDB, * ), D( * ), DF( * ), U( * ),
     $                   WORK( LWORK ), X( * )
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
      DOUBLE PRECISION   ANORM, BNORM, DNORM, EPS, UNFL, XNORM, YNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANGE
      EXTERNAL           DLAMCH, DZASUM, ZLANGE
*     ..
*     .. External Subroutines ..
*
      EXTERNAL           ZCOPY, ZGEMV, ZGGGLM, ZLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
      ANORM = MAX( ZLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( ZLANGE( '1', N, P, B, LDB, RWORK ), UNFL )
*
*     Copy the matrices A and B to the arrays AF and BF,
*     and the vector D the array DF.
*
      CALL ZLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL ZLACPY( 'Full', N, P, B, LDB, BF, LDB )
      CALL ZCOPY( N, D, 1, DF, 1 )
*
*     Solve GLM problem
*
      CALL ZGGGLM( N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK,
     $             INFO )
*
*     Test the residual for the solution of LSE
*
*                       norm( d - A*x - B*u )
*       RESULT = -----------------------------------------
*                (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
*
      CALL ZCOPY( N, D, 1, DF, 1 )
      CALL ZGEMV( 'No transpose', N, M, -CONE, A, LDA, X, 1, CONE, DF,
     $            1 )
*
      CALL ZGEMV( 'No transpose', N, P, -CONE, B, LDB, U, 1, CONE, DF,
     $            1 )
*
      DNORM = DZASUM( N, DF, 1 )
      XNORM = DZASUM( M, X, 1 ) + DZASUM( P, U, 1 )
      YNORM = ANORM + BNORM
*
      IF( XNORM.LE.ZERO ) THEN
         RESULT = ZERO
      ELSE
         RESULT = ( ( DNORM / YNORM ) / XNORM ) / EPS
      END IF
*
      RETURN
*
*     End of ZGLMTS
*
      END
