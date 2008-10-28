      SUBROUTINE CTPT01( UPLO, DIAG, N, AP, AINVP, RCOND, RWORK, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            N
      REAL               RCOND, RESID
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            AINVP( * ), AP( * )
*     ..
*
*  Purpose
*  =======
*
*  CTPT01 computes the residual for a triangular matrix A times its
*  inverse when A is stored in packed format:
*     RESID = norm(A*AINV - I) / ( N * norm(A) * norm(AINV) * EPS ),
*  where EPS is the machine epsilon.
*
*  Arguments
*  ==========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  AP      (input) COMPLEX array, dimension (N*(N+1)/2)
*          The original upper or lower triangular matrix A, packed
*          columnwise in a linear array.  The j-th column of A is stored
*          in the array AP as follows:
*          if UPLO = 'U', AP((j-1)*j/2 + i) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L',
*             AP((j-1)*(n-j) + j*(j+1)/2 + i-j) = A(i,j) for j<=i<=n.
*
*  AINVP   (input) COMPLEX array, dimension (N*(N+1)/2)
*          On entry, the (triangular) inverse of the matrix A, packed
*          columnwise in a linear array as in AP.
*          On exit, the contents of AINVP are destroyed.
*
*  RCOND   (output) REAL
*          The reciprocal condition number of A, computed as
*          1/(norm(A) * norm(AINV)).
*
*  RWORK   (workspace) REAL array, dimension (N)
*
*  RESID   (output) REAL
*          norm(A*AINV - I) / ( N * norm(A) * norm(AINV) * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UNITD
      INTEGER            J, JC
      REAL               AINVNM, ANORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANTP, SLAMCH
      EXTERNAL           LSAME, CLANTP, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CTPMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RCOND = ONE
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
*
      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANTP( '1', UPLO, DIAG, N, AP, RWORK )
      AINVNM = CLANTP( '1', UPLO, DIAG, N, AINVP, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
*     Compute A * AINV, overwriting AINV.
*
      UNITD = LSAME( DIAG, 'U' )
      IF( LSAME( UPLO, 'U' ) ) THEN
         JC = 1
         DO 10 J = 1, N
            IF( UNITD )
     $         AINVP( JC+J-1 ) = ONE
*
*           Form the j-th column of A*AINV.
*
            CALL CTPMV( 'Upper', 'No transpose', DIAG, J, AP,
     $                  AINVP( JC ), 1 )
*
*           Subtract 1 from the diagonal to form A*AINV - I.
*
            AINVP( JC+J-1 ) = AINVP( JC+J-1 ) - ONE
            JC = JC + J
   10    CONTINUE
      ELSE
         JC = 1
         DO 20 J = 1, N
            IF( UNITD )
     $         AINVP( JC ) = ONE
*
*           Form the j-th column of A*AINV.
*
            CALL CTPMV( 'Lower', 'No transpose', DIAG, N-J+1, AP( JC ),
     $                  AINVP( JC ), 1 )
*
*           Subtract 1 from the diagonal to form A*AINV - I.
*
            AINVP( JC ) = AINVP( JC ) - ONE
            JC = JC + N - J + 1
   20    CONTINUE
      END IF
*
*     Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = CLANTP( '1', UPLO, 'Non-unit', N, AINVP, RWORK )
*
      RESID = ( ( RESID*RCOND ) / REAL( N ) ) / EPS
*
      RETURN
*
*     End of CTPT01
*
      END
