      SUBROUTINE DTRTTP( UPLO, N, A, LDA, AP, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- Contributed by Fred Gustavson of the IBM Watson Research Center --
*  --            and Julien Langou of the Univ. of Colorado Denver    --
*  -- November 2008 --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N, LDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), AP( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTTP copies a triangular matrix A from full format (TR) to standard
*  packed format (TP).
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER
*          = 'U':  A is upper triangular.
*          = 'L':  A is lower triangular.
*
*  N       (input) INTEGER
*          The order of the matrices AP and A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          On exit, the triangular matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  AP      (output) DOUBLE PRECISION array, dimension (N*(N+1)/2
*          On exit, the upper or lower triangular matrix A, packed
*          columnwise in a linear array. The j-th column of A is stored
*          in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER
      INTEGER            I, J, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTTP', -INFO )
         RETURN
      END IF
*
      IF( LOWER ) THEN
         K = 0
         DO J = 1, N
            DO I = J, N
               K = K + 1
               AP( K ) = A( I, J )
            END DO
         END DO
      ELSE
         K = 0
         DO J = 1, N
            DO I = 1, J
               K = K + 1
               AP( K ) = A( I, J )
            END DO
         END DO
      END IF
*
*
      RETURN
*
*     End of DTRTTP
*
      END
