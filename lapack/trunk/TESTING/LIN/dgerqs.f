      SUBROUTINE DGERQS( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  Compute a minimum-norm solution
*      min || A*X - B ||
*  using the RQ factorization
*      A = R*Q
*  computed by DGERQF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= M >= 0.
*
*  NRHS    (input) INTEGER
*          The number of columns of B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          Details of the RQ factorization of the original matrix A as
*          returned by DGERQF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= M.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (M)
*          Details of the orthogonal matrix Q.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side vectors for the linear system.
*          On exit, the solution vectors X.  Each solution vector
*          is contained in rows 1:N of a column of B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK must be at least NRHS,
*          and should be at least NRHS*NB, where NB is the block size
*          for this environment.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASET, DORMRQ, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. M.GT.N ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.1 .OR. LWORK.LT.NRHS .AND. M.GT.0 .AND. N.GT.0 )
     $          THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGERQS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Solve R*X = B(n-m+1:n,:)
*
      CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', M, NRHS,
     $            ONE, A( 1, N-M+1 ), LDA, B( N-M+1, 1 ), LDB )
*
*     Set B(1:n-m,:) to zero
*
      CALL DLASET( 'Full', N-M, NRHS, ZERO, ZERO, B, LDB )
*
*     B := Q' * B
*
      CALL DORMRQ( 'Left', 'Transpose', N, NRHS, M, A, LDA, TAU, B, LDB,
     $             WORK, LWORK, INFO )
*
      RETURN
*
*     End of DGERQS
*
      END
