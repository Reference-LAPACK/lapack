      SUBROUTINE CGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF,
     $                   X, U, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, P, N
      REAL               RESULT
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), B( LDB, * ),
     $                   BF( LDB, * ), D( * ), DF( * ), U( * ),
     $                   WORK( LWORK ), X( * )
*
*  Purpose
*  =======
*
*  CGLMTS tests CGGGLM - a subroutine for solving the generalized
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
*  A       (input) COMPLEX array, dimension (LDA,M)
*          The N-by-M matrix A.
*
*  AF      (workspace) COMPLEX array, dimension (LDA,M)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A, AF. LDA >= max(M,N).
*
*  B       (input) COMPLEX array, dimension (LDB,P)
*          The N-by-P matrix A.
*
*  BF      (workspace) COMPLEX array, dimension (LDB,P)
*
*  LDB     (input) INTEGER
*          The leading dimension of the arrays B, BF. LDB >= max(P,N).
*
*  D       (input) COMPLEX array, dimension( N )
*          On input, the left hand side of the GLM.
*
*  DF      (workspace) COMPLEX array, dimension( N )
*
*  X       (output) COMPLEX array, dimension( M )
*          solution vector X in the GLM problem.
*
*  U       (output) COMPLEX array, dimension( P )
*          solution vector U in the GLM problem.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*  RWORK   (workspace) REAL array, dimension (M)
*
*  RESULT   (output) REAL
*          The test ratio:
*                           norm( d - A*x - B*u )
*            RESULT = -----------------------------------------
*                     (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
*
*  ====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
      REAL               ANORM, BNORM, EPS, XNORM, YNORM, DNORM, UNFL
*     ..
*     .. External Functions ..
      REAL               SCASUM, SLAMCH, CLANGE
      EXTERNAL           SCASUM, SLAMCH, CLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACPY
*
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
      UNFL = SLAMCH( 'Safe minimum' )
      ANORM = MAX( CLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( CLANGE( '1', N, P, B, LDB, RWORK ), UNFL )
*
*     Copy the matrices A and B to the arrays AF and BF,
*     and the vector D the array DF.
*
      CALL CLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL CLACPY( 'Full', N, P, B, LDB, BF, LDB )
      CALL CCOPY( N, D, 1, DF, 1 )
*
*     Solve GLM problem
*
      CALL CGGGLM( N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK,
     $             INFO )
*
*     Test the residual for the solution of LSE
*
*                       norm( d - A*x - B*u )
*       RESULT = -----------------------------------------
*                (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
*
      CALL CCOPY( N, D, 1, DF, 1 )
      CALL CGEMV( 'No transpose', N, M, -CONE, A, LDA, X, 1, CONE,
     $            DF, 1 )
*
      CALL CGEMV( 'No transpose', N, P, -CONE, B, LDB, U, 1, CONE,
     $            DF, 1 )
*
      DNORM = SCASUM( N, DF, 1 )
      XNORM = SCASUM( M, X, 1 ) + SCASUM( P, U, 1 )
      YNORM = ANORM + BNORM
*
      IF( XNORM.LE.ZERO ) THEN
         RESULT = ZERO
      ELSE
         RESULT =  ( ( DNORM / YNORM ) / XNORM ) /EPS
      END IF
*
      RETURN
*
*     End of CGLMTS
*
      END
