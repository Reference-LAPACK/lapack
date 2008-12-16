      DOUBLE PRECISION FUNCTION ZQRT17( TRANS, IRESID, M, N, NRHS, A,
     $                 LDA, X, LDX, B, LDB, C, WORK, LWORK )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDB, * ),
     $                   WORK( LWORK ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  ZQRT17 computes the ratio
*
*     || R'*op(A) ||/(||A||*alpha*max(M,N,NRHS)*eps)
*
*  where R = op(A)*X - B, op(A) is A or A', and
*
*     alpha = ||B|| if IRESID = 1 (zero-residual problem)
*     alpha = ||R|| if IRESID = 2 (otherwise).
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies whether or not the transpose of A is used.
*          = 'N':  No transpose, op(A) = A.
*          = 'C':  Conjugate transpose, op(A) = A'.
*
*  IRESID  (input) INTEGER
*          IRESID = 1 indicates zero-residual problem.
*          IRESID = 2 indicates non-zero residual.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*          If TRANS = 'N', the number of rows of the matrix B.
*          If TRANS = 'C', the number of rows of the matrix X.
*
*  N       (input) INTEGER
*          The number of columns of the matrix  A.
*          If TRANS = 'N', the number of rows of the matrix X.
*          If TRANS = 'C', the number of rows of the matrix B.
*
*  NRHS    (input) INTEGER
*          The number of columns of the matrices X and B.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The m-by-n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= M.
*
*  X       (input) COMPLEX*16 array, dimension (LDX,NRHS)
*          If TRANS = 'N', the n-by-nrhs matrix X.
*          If TRANS = 'C', the m-by-nrhs matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.
*          If TRANS = 'N', LDX >= N.
*          If TRANS = 'C', LDX >= M.
*
*  B       (input) COMPLEX*16 array, dimension (LDB,NRHS)
*          If TRANS = 'N', the m-by-nrhs matrix B.
*          If TRANS = 'C', the n-by-nrhs matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.
*          If TRANS = 'N', LDB >= M.
*          If TRANS = 'C', LDB >= N.
*
*  C       (workspace) COMPLEX*16 array, dimension (LDB,NRHS)
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= NRHS*(M+N).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, ISCL, NCOLS, NROWS
      DOUBLE PRECISION   BIGNUM, ERR, NORMA, NORMB, NORMRS, NORMX,
     $                   SMLNUM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZLACPY, ZLASCL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX
*     ..
*     .. Executable Statements ..
*
      ZQRT17 = ZERO
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         NROWS = M
         NCOLS = N
      ELSE IF( LSAME( TRANS, 'C' ) ) THEN
         NROWS = N
         NCOLS = M
      ELSE
         CALL XERBLA( 'ZQRT17', 1 )
         RETURN
      END IF
*
      IF( LWORK.LT.NCOLS*NRHS ) THEN
         CALL XERBLA( 'ZQRT17', 13 )
         RETURN
      END IF
*
      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 )
     $   RETURN
*
      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, RWORK )
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      ISCL = 0
*
*     compute residual and scale it
*
      CALL ZLACPY( 'All', NROWS, NRHS, B, LDB, C, LDB )
      CALL ZGEMM( TRANS, 'No transpose', NROWS, NRHS, NCOLS,
     $            DCMPLX( -ONE ), A, LDA, X, LDX, DCMPLX( ONE ), C,
     $            LDB )
      NORMRS = ZLANGE( 'Max', NROWS, NRHS, C, LDB, RWORK )
      IF( NORMRS.GT.SMLNUM ) THEN
         ISCL = 1
         CALL ZLASCL( 'General', 0, 0, NORMRS, ONE, NROWS, NRHS, C, LDB,
     $                INFO )
      END IF
*
*     compute R'*A
*
      CALL ZGEMM( 'Conjugate transpose', TRANS, NRHS, NCOLS, NROWS,
     $            DCMPLX( ONE ), C, LDB, A, LDA, DCMPLX( ZERO ), WORK,
     $            NRHS )
*
*     compute and properly scale error
*
      ERR = ZLANGE( 'One-norm', NRHS, NCOLS, WORK, NRHS, RWORK )
      IF( NORMA.NE.ZERO )
     $   ERR = ERR / NORMA
*
      IF( ISCL.EQ.1 )
     $   ERR = ERR*NORMRS
*
      IF( IRESID.EQ.1 ) THEN
         NORMB = ZLANGE( 'One-norm', NROWS, NRHS, B, LDB, RWORK )
         IF( NORMB.NE.ZERO )
     $      ERR = ERR / NORMB
      ELSE
         NORMX = ZLANGE( 'One-norm', NCOLS, NRHS, X, LDX, RWORK )
         IF( NORMX.NE.ZERO )
     $      ERR = ERR / NORMX
      END IF
*
      ZQRT17 = ERR / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N, NRHS ) ) )
      RETURN
*
*     End of ZQRT17
*
      END
