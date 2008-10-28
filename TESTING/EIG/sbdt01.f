      SUBROUTINE SBDT01( M, N, KD, A, LDA, Q, LDQ, D, E, PT, LDPT, WORK,
     $                   RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            KD, LDA, LDPT, LDQ, M, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), D( * ), E( * ), PT( LDPT, * ),
     $                   Q( LDQ, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SBDT01 reconstructs a general matrix A from its bidiagonal form
*     A = Q * B * P'
*  where Q (m by min(m,n)) and P' (min(m,n) by n) are orthogonal
*  matrices and B is bidiagonal.
*
*  The test ratio to test the reduction is
*     RESID = norm( A - Q * B * PT ) / ( n * norm(A) * EPS )
*  where PT = P' and EPS is the machine precision.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrices A and Q.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and P'.
*
*  KD      (input) INTEGER
*          If KD = 0, B is diagonal and the array E is not referenced.
*          If KD = 1, the reduction was performed by xGEBRD; B is upper
*          bidiagonal if M >= N, and lower bidiagonal if M < N.
*          If KD = -1, the reduction was performed by xGBBRD; B is
*          always upper bidiagonal.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  Q       (input) REAL array, dimension (LDQ,N)
*          The m by min(m,n) orthogonal matrix Q in the reduction
*          A = Q * B * P'.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= max(1,M).
*
*  D       (input) REAL array, dimension (min(M,N))
*          The diagonal elements of the bidiagonal matrix B.
*
*  E       (input) REAL array, dimension (min(M,N)-1)
*          The superdiagonal elements of the bidiagonal matrix B if
*          m >= n, or the subdiagonal elements of B if m < n.
*
*  PT      (input) REAL array, dimension (LDPT,N)
*          The min(m,n) by n orthogonal matrix P' in the reduction
*          A = Q * B * P'.
*
*  LDPT    (input) INTEGER
*          The leading dimension of the array PT.
*          LDPT >= max(1,min(M,N)).
*
*  WORK    (workspace) REAL array, dimension (M+N)
*
*  RESID   (output) REAL
*          The test ratio:  norm(A - Q * B * P') / ( n * norm(A) * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               ANORM, EPS
*     ..
*     .. External Functions ..
      REAL               SASUM, SLAMCH, SLANGE
      EXTERNAL           SASUM, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Compute A - Q * B * P' one column at a time.
*
      RESID = ZERO
      IF( KD.NE.0 ) THEN
*
*        B is bidiagonal.
*
         IF( KD.NE.0 .AND. M.GE.N ) THEN
*
*           B is upper bidiagonal and M >= N.
*
            DO 20 J = 1, N
               CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 10 I = 1, N - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   10          CONTINUE
               WORK( M+N ) = D( N )*PT( N, J )
               CALL SGEMV( 'No transpose', M, N, -ONE, Q, LDQ,
     $                     WORK( M+1 ), 1, ONE, WORK, 1 )
               RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
   20       CONTINUE
         ELSE IF( KD.LT.0 ) THEN
*
*           B is upper bidiagonal and M < N.
*
            DO 40 J = 1, N
               CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 30 I = 1, M - 1
                  WORK( M+I ) = D( I )*PT( I, J ) + E( I )*PT( I+1, J )
   30          CONTINUE
               WORK( M+M ) = D( M )*PT( M, J )
               CALL SGEMV( 'No transpose', M, M, -ONE, Q, LDQ,
     $                     WORK( M+1 ), 1, ONE, WORK, 1 )
               RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
   40       CONTINUE
         ELSE
*
*           B is lower bidiagonal.
*
            DO 60 J = 1, N
               CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
               WORK( M+1 ) = D( 1 )*PT( 1, J )
               DO 50 I = 2, M
                  WORK( M+I ) = E( I-1 )*PT( I-1, J ) +
     $                          D( I )*PT( I, J )
   50          CONTINUE
               CALL SGEMV( 'No transpose', M, M, -ONE, Q, LDQ,
     $                     WORK( M+1 ), 1, ONE, WORK, 1 )
               RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
   60       CONTINUE
         END IF
      ELSE
*
*        B is diagonal.
*
         IF( M.GE.N ) THEN
            DO 80 J = 1, N
               CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 70 I = 1, N
                  WORK( M+I ) = D( I )*PT( I, J )
   70          CONTINUE
               CALL SGEMV( 'No transpose', M, N, -ONE, Q, LDQ,
     $                     WORK( M+1 ), 1, ONE, WORK, 1 )
               RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
   80       CONTINUE
         ELSE
            DO 100 J = 1, N
               CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
               DO 90 I = 1, M
                  WORK( M+I ) = D( I )*PT( I, J )
   90          CONTINUE
               CALL SGEMV( 'No transpose', M, M, -ONE, Q, LDQ,
     $                     WORK( M+1 ), 1, ONE, WORK, 1 )
               RESID = MAX( RESID, SASUM( M, WORK, 1 ) )
  100       CONTINUE
         END IF
      END IF
*
*     Compute norm(A - Q * B * P') / ( n * norm(A) * EPS )
*
      ANORM = SLANGE( '1', M, N, A, LDA, WORK )
      EPS = SLAMCH( 'Precision' )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         IF( ANORM.GE.RESID ) THEN
            RESID = ( RESID / ANORM ) / ( REAL( N )*EPS )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REAL( N )*ANORM ) / ANORM ) /
     $                 ( REAL( N )*EPS )
            ELSE
               RESID = MIN( RESID / ANORM, REAL( N ) ) /
     $                 ( REAL( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of SBDT01
*
      END
