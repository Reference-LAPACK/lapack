      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
      IMPLICIT NONE
*
*  -- LAPACK routine (version 3.X) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     May 2008
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This code implements an iterative version of Sivan Toledo's recursive
*  LU algorithm[1].  For square matrices, this iterative versions should
*  be within a factor of two of the optimum number of memory transfers.
*
*  The pattern is as follows, with the large blocks of U being updated
*  in one call to DTRSM, and the dotted lines denoting sections that
*  have had all pending permutations applied:
*
*   1 2 3 4 5 6 7 8
*  +-+-+---+-------+------
*  | |1|   |       |
*  |.+-+ 2 |       |
*  | | |   |       |
*  |.|.+-+-+   4   |
*  | | | |1|       |
*  | | |.+-+       |
*  | | | | |       |
*  |.|.|.|.+-+-+---+  8
*  | | | | | |1|   |
*  | | | | |.+-+ 2 |
*  | | | | | | |   |
*  | | | | |.|.+-+-+
*  | | | | | | | |1|
*  | | | | | | |.+-+
*  | | | | | | | | |
*  |.|.|.|.|.|.|.|.+-----
*  | | | | | | | | |
*
*  The 1-2-1-4-1-2-1-8-... pattern is the position of the last 1 bit in
*  the binary expansion of the current column.  Each Schur update is
*  applied as soon as the necessary portion of U is available.
*
*  [1] Toledo, S. 1997. Locality of Reference in LU Decomposition with
*  Partial Pivoting. SIAM J. Matrix Anal. Appl. 18, 4 (Oct. 1997),
*  1065-1081. http://dx.doi.org/10.1137/S0895479896297744
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, NEGONE
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ONE = (1.0D+0, 0.0D+0) )
      PARAMETER          ( NEGONE = (-1.0D+0, 0.0D+0) )
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN, PIVMAG
      COMPLEX*16         TMP
      INTEGER            I, J, JP, NSTEP, NTOPIV, NPIVED, KAHEAD
      INTEGER            KSTART, IPIVSTART, JPIVSTART, KCOLS
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IZAMAX
      LOGICAL            DISNAN
      EXTERNAL           DLAMCH, IZAMAX, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZTRSM, ZSCAL, XERBLA, ZLASWP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, IAND, ABS
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum
*
      SFMIN = DLAMCH( 'S' )
*
      NSTEP = MIN( M, N )
      DO J = 1, NSTEP
         KAHEAD = IAND( J, -J )
         KSTART = J + 1 - KAHEAD
         KCOLS = MIN( KAHEAD, M-J )
*
*        Find pivot.
*
         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP

!        Permute just this column.
         IF (JP .NE. J) THEN
            TMP = A( J, J )
            A( J, J ) = A( JP, J )
            A( JP, J ) = TMP
         END IF

!        Apply pending permutations to L
         NTOPIV = 1
         IPIVSTART = J
         JPIVSTART = J - NTOPIV
         DO WHILE ( NTOPIV .LT. KAHEAD )
            CALL ZLASWP( NTOPIV, A( 1, JPIVSTART ), LDA, IPIVSTART, J,
     $           IPIV, 1 )
            IPIVSTART = IPIVSTART - NTOPIV;
            NTOPIV = NTOPIV * 2;
            JPIVSTART = JPIVSTART - NTOPIV;
         END DO

!        Permute U block to match L
         CALL ZLASWP( KCOLS, A( 1,J+1 ), LDA, KSTART, J, IPIV, 1 )

!        Factor the current column
         PIVMAG = ABS( A( J, J ) )
         IF( PIVMAG.NE.ZERO .AND. .NOT.DISNAN( PIVMAG ) ) THEN
               IF( PIVMAG .GE. SFMIN ) THEN
                  CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
               ELSE
                 DO I = 1, M-J
                    A( J+I, J ) = A( J+I, J ) / A( J, J )
                 END DO
               END IF
         ELSE IF( PIVMAG .EQ. ZERO .AND. INFO .EQ. 0 ) THEN
            INFO = J
         END IF

!        Solve for U block.
         CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', KAHEAD,
     $        KCOLS, ONE, A( KSTART, KSTART ), LDA,
     $        A( KSTART, J+1 ), LDA )
!        Schur complement.
         CALL ZGEMM( 'No transpose', 'No transpose', M-J,
     $        KCOLS, KAHEAD, NEGONE, A( J+1, KSTART ), LDA,
     $        A( KSTART, J+1 ), LDA, ONE, A( J+1, J+1 ), LDA )
      END DO

!     Handle pivot permutations on the way out of the recursion
      NPIVED = IAND( NSTEP, -NSTEP )
      J = NSTEP - NPIVED
      DO WHILE ( J .GT. 0 )
         NTOPIV = IAND( J, -J )
         CALL ZLASWP( NTOPIV, A( 1, J-NTOPIV+1 ), LDA, J+1, NSTEP,
     $        IPIV, 1 )
         J = J - NTOPIV
      END DO

!     If short and wide, handle the rest of the columns.
      IF ( M .LT. N ) THEN
         CALL ZLASWP( N-M, A( 1, M+KCOLS+1 ), LDA, 1, M, IPIV, 1 )
         CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', M,
     $        N-M, ONE, A, LDA, A( 1,M+KCOLS+1 ), LDA )
      END IF

      RETURN
*
*     End of ZGETRF
*
      END
