      DOUBLE PRECISION FUNCTION ZQPT01( M, N, K, A, AF, LDA, TAU, JPVT,
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
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  ZQPT01 tests the QR-factorization with pivoting of a matrix A.  The
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
*  A       (input) COMPLEX*16 array, dimension (LDA, N)
*          The original matrix A.
*
*  AF      (input) COMPLEX*16 array, dimension (LDA,N)
*          The (possibly partial) output of ZGEQPF.  The upper triangle
*          of AF(1:k,1:k) is a partial triangular factor, the entries
*          below the diagonal in the first k columns are the Householder
*          vectors, and the rest of AF contains a partially updated
*          matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays A and AF.
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          Details of the Householder transformations as returned by
*          ZGEQPF.
*
*  JPVT    (input) INTEGER array, dimension (N)
*          Pivot information as returned by ZGEQPF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= M*N+N.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   NORMA
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZCOPY, ZUNMQR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      ZQPT01 = ZERO
*
*     Test if there is enough workspace
*
      IF( LWORK.LT.M*N+N ) THEN
         CALL XERBLA( 'ZQPT01', 10 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, RWORK )
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
         CALL ZCOPY( M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 )
   40 CONTINUE
*
      CALL ZUNMQR( 'Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK,
     $             M, WORK( M*N+1 ), LWORK-M*N, INFO )
*
      DO 50 J = 1, N
*
*        Compare i-th column of QR and jpvt(i)-th column of A
*
         CALL ZAXPY( M, DCMPLX( -ONE ), A( 1, JPVT( J ) ), 1,
     $               WORK( ( J-1 )*M+1 ), 1 )
   50 CONTINUE
*
      ZQPT01 = ZLANGE( 'One-norm', M, N, WORK, M, RWORK ) /
     $         ( DBLE( MAX( M, N ) )*DLAMCH( 'Epsilon' ) )
      IF( NORMA.NE.ZERO )
     $   ZQPT01 = ZQPT01 / NORMA
*
      RETURN
*
*     End of ZQPT01
*
      END
