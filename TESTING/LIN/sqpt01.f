      REAL             FUNCTION SQPT01( M, N, K, A, AF, LDA, TAU, JPVT,
     $                 WORK, LWORK )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               A( LDA, * ), AF( LDA, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SQPT01 tests the QR-factorization with pivoting of a matrix A.  The
*  array AF contains the (possibly partial) QR-factorization of A, where
*  the upper triangle of AF(1:k,1:k) is a partial triangular factor,
*  the entries below the diagonal in the first k columns are the
*  Householder vectors, and the rest of AF contains a partially updated
*  matrix.
*
*  This function returns ||A*P - Q*R||/(||norm(A)||*eps*M)
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrices A and AF.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and AF.
*
*  K       (input) INTEGER
*          The number of columns of AF that have been reduced
*          to upper triangular form.
*
*  A       (input) REAL array, dimension (LDA, N)
*          The original matrix A.
*
*  AF      (input) REAL array, dimension (LDA,N)
*          The (possibly partial) output of SGEQPF.  The upper triangle
*          of AF(1:k,1:k) is a partial triangular factor, the entries
*          below the diagonal in the first k columns are the Householder
*          vectors, and the rest of AF contains a partially updated
*          matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A and AF.
*
*  TAU     (input) REAL array, dimension (K)
*          Details of the Householder transformations as returned by
*          SGEQPF.
*
*  JPVT    (input) INTEGER array, dimension (N)
*          Pivot information as returned by SGEQPF.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= M*N+N.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      REAL               NORMA
*     ..
*     .. Local Arrays ..
      REAL               RWORK( 1 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SORMQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      SQPT01 = ZERO
*
*     Test if there is enough workspace
*
      IF( LWORK.LT.M*N+N ) THEN
         CALL XERBLA( 'SQPT01', 10 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      NORMA = SLANGE( 'One-norm', M, N, A, LDA, RWORK )
*
      DO 30 J = 1, K
         DO 10 I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = AF( I, J )
   10    CONTINUE
         DO 20 I = J + 1, M
            WORK( ( J-1 )*M+I ) = ZERO
   20    CONTINUE
   30 CONTINUE
      DO 40 J = K + 1, N
         CALL SCOPY( M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 )
   40 CONTINUE
*
      CALL SORMQR( 'Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK,
     $             M, WORK( M*N+1 ), LWORK-M*N, INFO )
*
      DO 50 J = 1, N
*
*        Compare i-th column of QR and jpvt(i)-th column of A
*
         CALL SAXPY( M, -ONE, A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ),
     $               1 )
   50 CONTINUE
*
      SQPT01 = SLANGE( 'One-norm', M, N, WORK, M, RWORK ) /
     $         ( REAL( MAX( M, N ) )*SLAMCH( 'Epsilon' ) )
      IF( NORMA.NE.ZERO )
     $   SQPT01 = SQPT01 / NORMA
*
      RETURN
*
*     End of SQPT01
*
      END
