      SUBROUTINE SPPT03( UPLO, N, A, AINV, WORK, LDWORK, RWORK, RCOND,
     $                   RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDWORK, N
      REAL               RCOND, RESID
*     ..
*     .. Array Arguments ..
      REAL               A( * ), AINV( * ), RWORK( * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  SPPT03 computes the residual for a symmetric packed matrix times its
*  inverse:
*     norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS ),
*  where EPS is the machine epsilon.
*
*  Arguments
*  ==========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The number of rows and columns of the matrix A.  N >= 0.
*
*  A       (input) REAL array, dimension (N*(N+1)/2)
*          The original symmetric matrix A, stored as a packed
*          triangular matrix.
*
*  AINV    (input) REAL array, dimension (N*(N+1)/2)
*          The (symmetric) inverse of the matrix A, stored as a packed
*          triangular matrix.
*
*  WORK    (workspace) REAL array, dimension (LDWORK,N)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.  LDWORK >= max(1,N).
*
*  RWORK   (workspace) REAL array, dimension (N)
*
*  RCOND   (output) REAL
*          The reciprocal of the condition number of A, computed as
*          ( 1/norm(A) ) / norm(AINV).
*
*  RESID   (output) REAL
*          norm(I - A*AINV) / ( N * norm(A) * norm(AINV) * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JJ
      REAL               AINVNM, ANORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANGE, SLANSP
      EXTERNAL           LSAME, SLAMCH, SLANGE, SLANSP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SSPMV
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
      ANORM = SLANSP( '1', UPLO, N, A, RWORK )
      AINVNM = SLANSP( '1', UPLO, N, AINV, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.EQ.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE / ANORM ) / AINVNM
*
*     UPLO = 'U':
*     Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
*     expand it to a full matrix, then multiply by A one column at a
*     time, moving the result one column to the left.
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Copy AINV
*
         JJ = 1
         DO 10 J = 1, N - 1
            CALL SCOPY( J, AINV( JJ ), 1, WORK( 1, J+1 ), 1 )
            CALL SCOPY( J-1, AINV( JJ ), 1, WORK( J, 2 ), LDWORK )
            JJ = JJ + J
   10    CONTINUE
         JJ = ( ( N-1 )*N ) / 2 + 1
         CALL SCOPY( N-1, AINV( JJ ), 1, WORK( N, 2 ), LDWORK )
*
*        Multiply by A
*
         DO 20 J = 1, N - 1
            CALL SSPMV( 'Upper', N, -ONE, A, WORK( 1, J+1 ), 1, ZERO,
     $                  WORK( 1, J ), 1 )
   20    CONTINUE
         CALL SSPMV( 'Upper', N, -ONE, A, AINV( JJ ), 1, ZERO,
     $               WORK( 1, N ), 1 )
*
*     UPLO = 'L':
*     Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
*     and multiply by A, moving each column to the right.
*
      ELSE
*
*        Copy AINV
*
         CALL SCOPY( N-1, AINV( 2 ), 1, WORK( 1, 1 ), LDWORK )
         JJ = N + 1
         DO 30 J = 2, N
            CALL SCOPY( N-J+1, AINV( JJ ), 1, WORK( J, J-1 ), 1 )
            CALL SCOPY( N-J, AINV( JJ+1 ), 1, WORK( J, J ), LDWORK )
            JJ = JJ + N - J + 1
   30    CONTINUE
*
*        Multiply by A
*
         DO 40 J = N, 2, -1
            CALL SSPMV( 'Lower', N, -ONE, A, WORK( 1, J-1 ), 1, ZERO,
     $                  WORK( 1, J ), 1 )
   40    CONTINUE
         CALL SSPMV( 'Lower', N, -ONE, A, AINV( 1 ), 1, ZERO,
     $               WORK( 1, 1 ), 1 )
*
      END IF
*
*     Add the identity matrix to WORK .
*
      DO 50 I = 1, N
         WORK( I, I ) = WORK( I, I ) + ONE
   50 CONTINUE
*
*     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
*
      RESID = SLANGE( '1', N, N, WORK, LDWORK, RWORK )
*
      RESID = ( ( RESID*RCOND ) / EPS ) / REAL( N )
*
      RETURN
*
*     End of SPPT03
*
      END
