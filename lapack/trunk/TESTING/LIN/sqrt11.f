      REAL             FUNCTION SQRT11( M, K, A, LDA, TAU, WORK, LWORK )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LWORK, M
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SQRT11 computes the test ratio
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
*  A       (input) REAL array, dimension (LDA,K)
*          The (possibly partial) output of a QR reduction routine.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  TAU     (input) REAL array, dimension (K)
*          The scaling factors tau for the elementary transformations as
*          computed by the QR factorization routine.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= M*M + M.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, J
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASET, SORM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Local Arrays ..
      REAL               RDUMMY( 1 )
*     ..
*     .. Executable Statements ..
*
      SQRT11 = ZERO
*
*     Test for sufficient workspace
*
      IF( LWORK.LT.M*M+M ) THEN
         CALL XERBLA( 'SQRT11', 7 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 )
     $   RETURN
*
      CALL SLASET( 'Full', M, M, ZERO, ONE, WORK, M )
*
*     Form Q
*
      CALL SORM2R( 'Left', 'No transpose', M, M, K, A, LDA, TAU, WORK,
     $             M, WORK( M*M+1 ), INFO )
*
*     Form Q'*Q
*
      CALL SORM2R( 'Left', 'Transpose', M, M, K, A, LDA, TAU, WORK, M,
     $             WORK( M*M+1 ), INFO )
*
      DO 10 J = 1, M
         WORK( ( J-1 )*M+J ) = WORK( ( J-1 )*M+J ) - ONE
   10 CONTINUE
*
      SQRT11 = SLANGE( 'One-norm', M, M, WORK, M, RDUMMY ) /
     $         ( REAL( M )*SLAMCH( 'Epsilon' ) )
*
      RETURN
*
*     End of SQRT11
*
      END
