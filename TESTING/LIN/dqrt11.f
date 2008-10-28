      DOUBLE PRECISION FUNCTION DQRT11( M, K, A, LDA, TAU, WORK, LWORK )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LWORK, M
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DQRT11 computes the test ratio
*
*        || Q'*Q - I || / (eps * m)
*
*  where the orthogonal matrix Q is represented as a product of
*  elementary transformations.  Each transformation has the form
*
*     H(k) = I - tau(k) v(k) v(k)'
*
*  where tau(k) is stored in TAU(k) and v(k) is an m-vector of the form
*  [ 0 ... 0 1 x(k) ]', where x(k) is a vector of length m-k stored
*  in A(k+1:m,k).
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  K       (input) INTEGER
*          The number of columns of A whose subdiagonal entries
*          contain information about orthogonal transformations.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
*          The (possibly partial) output of a QR reduction routine.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          The scaling factors tau for the elementary transformations as
*          computed by the QR factorization routine.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= M*M + M.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, J
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASET, DORM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RDUMMY( 1 )
*     ..
*     .. Executable Statements ..
*
      DQRT11 = ZERO
*
*     Test for sufficient workspace
*
      IF( LWORK.LT.M*M+M ) THEN
         CALL XERBLA( 'DQRT11', 7 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 )
     $   RETURN
*
      CALL DLASET( 'Full', M, M, ZERO, ONE, WORK, M )
*
*     Form Q
*
      CALL DORM2R( 'Left', 'No transpose', M, M, K, A, LDA, TAU, WORK,
     $             M, WORK( M*M+1 ), INFO )
*
*     Form Q'*Q
*
      CALL DORM2R( 'Left', 'Transpose', M, M, K, A, LDA, TAU, WORK, M,
     $             WORK( M*M+1 ), INFO )
*
      DO 10 J = 1, M
         WORK( ( J-1 )*M+J ) = WORK( ( J-1 )*M+J ) - ONE
   10 CONTINUE
*
      DQRT11 = DLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) /
     $         ( DBLE( M )*DLAMCH( 'Epsilon' ) )
*
      RETURN
*
*     End of DQRT11
*
      END
