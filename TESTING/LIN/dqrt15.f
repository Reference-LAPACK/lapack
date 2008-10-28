      SUBROUTINE DQRT15( SCALE, RKSEL, M, N, NRHS, A, LDA, B, LDB, S,
     $                   RANK, NORMA, NORMB, ISEED, WORK, LWORK )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, N, NRHS, RANK, RKSEL, SCALE
      DOUBLE PRECISION   NORMA, NORMB
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DQRT15 generates a matrix with full or deficient rank and of various
*  norms.
*
*  Arguments
*  =========
*
*  SCALE   (input) INTEGER
*          SCALE = 1: normally scaled matrix
*          SCALE = 2: matrix scaled up
*          SCALE = 3: matrix scaled down
*
*  RKSEL   (input) INTEGER
*          RKSEL = 1: full rank matrix
*          RKSEL = 2: rank-deficient matrix
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  N       (input) INTEGER
*          The number of columns of A.
*
*  NRHS    (input) INTEGER
*          The number of columns of B.
*
*  A       (output) DOUBLE PRECISION array, dimension (LDA,N)
*          The M-by-N matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  B       (output) DOUBLE PRECISION array, dimension (LDB, NRHS)
*          A matrix that is in the range space of matrix A.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.
*
*  S       (output) DOUBLE PRECISION array, dimension MIN(M,N)
*          Singular values of A.
*
*  RANK    (output) INTEGER
*          number of nonzero singular values of A.
*
*  NORMA   (output) DOUBLE PRECISION
*          one-norm of A.
*
*  NORMB   (output) DOUBLE PRECISION
*          one-norm of B.
*
*  ISEED   (input/output) integer array, dimension (4)
*          seed for random number generator.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          length of work space required.
*          LWORK >= MAX(M+MIN(M,N),NRHS*MIN(M,N),2*N+M)
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, SVMIN
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   SVMIN = 0.1D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, J, MN
      DOUBLE PRECISION   BIGNUM, EPS, SMLNUM, TEMP
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DUMMY( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DASUM, DLAMCH, DLANGE, DLARND, DNRM2
      EXTERNAL           DASUM, DLAMCH, DLANGE, DLARND, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLAORD, DLARF, DLARNV, DLAROR, DLASCL,
     $                   DLASET, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      MN = MIN( M, N )
      IF( LWORK.LT.MAX( M+MN, MN*NRHS, 2*N+M ) ) THEN
         CALL XERBLA( 'DQRT15', 16 )
         RETURN
      END IF
*
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      EPS = DLAMCH( 'Epsilon' )
      SMLNUM = ( SMLNUM / EPS ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Determine rank and (unscaled) singular values
*
      IF( RKSEL.EQ.1 ) THEN
         RANK = MN
      ELSE IF( RKSEL.EQ.2 ) THEN
         RANK = ( 3*MN ) / 4
         DO 10 J = RANK + 1, MN
            S( J ) = ZERO
   10    CONTINUE
      ELSE
         CALL XERBLA( 'DQRT15', 2 )
      END IF
*
      IF( RANK.GT.0 ) THEN
*
*        Nontrivial case
*
         S( 1 ) = ONE
         DO 30 J = 2, RANK
   20       CONTINUE
            TEMP = DLARND( 1, ISEED )
            IF( TEMP.GT.SVMIN ) THEN
               S( J ) = ABS( TEMP )
            ELSE
               GO TO 20
            END IF
   30    CONTINUE
         CALL DLAORD( 'Decreasing', RANK, S, 1 )
*
*        Generate 'rank' columns of a random orthogonal matrix in A
*
         CALL DLARNV( 2, ISEED, M, WORK )
         CALL DSCAL( M, ONE / DNRM2( M, WORK, 1 ), WORK, 1 )
         CALL DLASET( 'Full', M, RANK, ZERO, ONE, A, LDA )
         CALL DLARF( 'Left', M, RANK, WORK, 1, TWO, A, LDA,
     $               WORK( M+1 ) )
*
*        workspace used: m+mn
*
*        Generate consistent rhs in the range space of A
*
         CALL DLARNV( 2, ISEED, RANK*NRHS, WORK )
         CALL DGEMM( 'No transpose', 'No transpose', M, NRHS, RANK, ONE,
     $               A, LDA, WORK, RANK, ZERO, B, LDB )
*
*        work space used: <= mn *nrhs
*
*        generate (unscaled) matrix A
*
         DO 40 J = 1, RANK
            CALL DSCAL( M, S( J ), A( 1, J ), 1 )
   40    CONTINUE
         IF( RANK.LT.N )
     $      CALL DLASET( 'Full', M, N-RANK, ZERO, ZERO, A( 1, RANK+1 ),
     $                   LDA )
         CALL DLAROR( 'Right', 'No initialization', M, N, A, LDA, ISEED,
     $                WORK, INFO )
*
      ELSE
*
*        work space used 2*n+m
*
*        Generate null matrix and rhs
*
         DO 50 J = 1, MN
            S( J ) = ZERO
   50    CONTINUE
         CALL DLASET( 'Full', M, N, ZERO, ZERO, A, LDA )
         CALL DLASET( 'Full', M, NRHS, ZERO, ZERO, B, LDB )
*
      END IF
*
*     Scale the matrix
*
      IF( SCALE.NE.1 ) THEN
         NORMA = DLANGE( 'Max', M, N, A, LDA, DUMMY )
         IF( NORMA.NE.ZERO ) THEN
            IF( SCALE.EQ.2 ) THEN
*
*              matrix scaled up
*
               CALL DLASCL( 'General', 0, 0, NORMA, BIGNUM, M, N, A,
     $                      LDA, INFO )
               CALL DLASCL( 'General', 0, 0, NORMA, BIGNUM, MN, 1, S,
     $                      MN, INFO )
               CALL DLASCL( 'General', 0, 0, NORMA, BIGNUM, M, NRHS, B,
     $                      LDB, INFO )
            ELSE IF( SCALE.EQ.3 ) THEN
*
*              matrix scaled down
*
               CALL DLASCL( 'General', 0, 0, NORMA, SMLNUM, M, N, A,
     $                      LDA, INFO )
               CALL DLASCL( 'General', 0, 0, NORMA, SMLNUM, MN, 1, S,
     $                      MN, INFO )
               CALL DLASCL( 'General', 0, 0, NORMA, SMLNUM, M, NRHS, B,
     $                      LDB, INFO )
            ELSE
               CALL XERBLA( 'DQRT15', 1 )
               RETURN
            END IF
         END IF
      END IF
*
      NORMA = DASUM( MN, S, 1 )
      NORMB = DLANGE( 'One-norm', M, NRHS, B, LDB, DUMMY )
*
      RETURN
*
*     End of DQRT15
*
      END
