      SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
     $                  INFO )
*
*  -- LAPACK driver routine (instrumented to count ops, version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*     Common block to return operation counts and timings.
*     .. Common blocks ..
      COMMON             / LSTIME / OPCNT, TIMNG
*     ..
*     .. Arrays in Common ..
      DOUBLE PRECISION   OPCNT( 6 ), TIMNG( 6 )
*     ..
*
*  Purpose
*  =======
*
*  ZGELS solves overdetermined or underdetermined complex linear systems
*  involving an M-by-N matrix A, or its conjugate-transpose, using a QR
*  or LQ factorization of A.  It is assumed that A has full rank.
*
*  The following options are provided:
*
*  1. If TRANS = 'N' and m >= n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || B - A*X ||.
*
*  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
*     an underdetermined system A * X = B.
*
*  3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
*     an undetermined system A**H * X = B.
*
*  4. If TRANS = 'C' and m < n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || B - A**H * X ||.
*
*  Several right hand side vectors b and solution vectors x can be
*  handled in a single call; they are stored as the columns of the
*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
*  matrix X.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER
*          = 'N': the linear system involves A;
*          = 'C': the linear system involves A**H.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of
*          columns of the matrices B and X. NRHS >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit,
*            if M >= N, A is overwritten by details of its QR
*                       factorization as returned by ZGEQRF;
*            if M <  N, A is overwritten by details of its LQ
*                       factorization as returned by ZGELQF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the matrix B of right hand side vectors, stored
*          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
*          if TRANS = 'C'.
*          On exit, B is overwritten by the solution vectors, stored
*          columnwise:
*          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
*          squares solution vectors; the residual sum of squares for the
*          solution in each column is given by the sum of squares of
*          elements N+1 to M in that column;
*          if TRANS = 'N' and m < n, rows 1 to N of B contain the
*          minimum norm solution vectors;
*          if TRANS = 'C' and m >= n, rows 1 to M of B contain the
*          minimum norm solution vectors;
*          if TRANS = 'C' and m < n, rows 1 to M of B contain the
*          least squares solution vectors; the residual sum of squares
*          for the solution in each column is given by the sum of
*          squares of elements M+1 to N in that column.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= MAX(1,M,N).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          LWORK >= max( 1, MN + max( MN, NRHS ) ).
*          For optimal performance,
*          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
*          where MN = min(M,N) and NB is the optimum block size.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, TPSD
      INTEGER            BROW, GELQF, GELS, GEQRF, I, IASCL, IBSCL, J,
     $                   MN, NB, SCLLEN, TRSM, UNMLQ, UNMQR, WSIZE
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, SMLNUM, T1, T2
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DOPBL3, DOPLA, DSECND, ZLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, DOPBL3, DOPLA, DSECND,
     $                   ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZGELQF, ZGEQRF, ZLASCL, ZLASET,
     $                   ZTRSM, ZUNMLQ, ZUNMQR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Data statements ..
      DATA               GELQF / 2 / , GELS / 1 / , GEQRF / 2 / ,
     $                   TRSM / 4 / , UNMLQ / 3 / , UNMQR / 3 /
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, MN+MAX( MN, NRHS ) ) .AND. .NOT.LQUERY )
     $          THEN
         INFO = -10
      END IF
*
*     Figure out optimal block size
*
      IF( INFO.EQ.0 .OR. INFO.EQ.-10 ) THEN
*
         TPSD = .TRUE.
         IF( LSAME( TRANS, 'N' ) )
     $      TPSD = .FALSE.
*
         IF( M.GE.N ) THEN
            NB = ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'ZUNMQR', 'LN', M, NRHS, N,
     $              -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'ZUNMQR', 'LC', M, NRHS, N,
     $              -1 ) )
            END IF
         ELSE
            NB = ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'ZUNMLQ', 'LC', N, NRHS, M,
     $              -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'ZUNMLQ', 'LN', N, NRHS, M,
     $              -1 ) )
            END IF
         END IF
*
         WSIZE = MAX( 1, MN+MAX( MN, NRHS )*NB )
         WORK( 1 ) = DBLE( WSIZE )
*
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGELS ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         CALL ZLASET( 'Full', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         RETURN
      END IF
*
*     Get machine parameters
*
      OPCNT( GELS ) = OPCNT( GELS ) + DBLE( 2 )
      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*
*     Scale A, B if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         OPCNT( GELS ) = OPCNT( GELS ) + DBLE( 6*M*N )
         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         OPCNT( GELS ) = OPCNT( GELS ) + DBLE( 6*M*N )
         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL ZLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         GO TO 50
      END IF
*
      BROW = M
      IF( TPSD )
     $   BROW = N
      BNRM = ZLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         OPCNT( GELS ) = OPCNT( GELS ) + DBLE( 6*BROW*NRHS )
         CALL ZLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB,
     $                INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         OPCNT( GELS ) = OPCNT( GELS ) + DBLE( 6*BROW*NRHS )
         CALL ZLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB,
     $                INFO )
         IBSCL = 2
      END IF
*
      IF( M.GE.N ) THEN
*
*        compute QR factorization of A
*
         NB = ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
         OPCNT( GEQRF ) = OPCNT( GEQRF ) +
     $                    DOPLA( 'ZGEQRF', M, N, 0, 0, NB )
         T1 = DSECND( )
         CALL ZGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN,
     $                INFO )
         T2 = DSECND( )
         TIMNG( GEQRF ) = TIMNG( GEQRF ) + ( T2-T1 )
*
*        workspace at least N, optimally N*NB
*
         IF( .NOT.TPSD ) THEN
*
*           Least-Squares Problem min || A * X - B ||
*
*           B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS)
*
            NB = ILAENV( 1, 'ZUNMQR', 'LC', M, NRHS, N, -1 )
            OPCNT( UNMQR ) = OPCNT( UNMQR ) +
     $                       DOPLA( 'ZUNMQR', M, NRHS, N, 0, NB )
            T1 = DSECND( )
            CALL ZUNMQR( 'Left', 'Conjugate transpose', M, NRHS, N, A,
     $                   LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $                   INFO )
            T2 = DSECND( )
            TIMNG( UNMQR ) = TIMNG( UNMQR ) + ( T2-T1 )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
*
            OPCNT( TRSM ) = OPCNT( TRSM ) +
     $                      DOPBL3( 'ZTRSM ', N, NRHS, 0 )
            T1 = DSECND( )
            CALL ZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $                  NRHS, CONE, A, LDA, B, LDB )
            T2 = DSECND( )
            TIMNG( TRSM ) = TIMNG( TRSM ) + ( T2-T1 )
*
            SCLLEN = N
*
         ELSE
*
*           Overdetermined system of equations A' * X = B
*
*           B(1:N,1:NRHS) := inv(R') * B(1:N,1:NRHS)
*
            OPCNT( TRSM ) = OPCNT( TRSM ) +
     $                      DOPBL3( 'ZTRSM ', N, NRHS, 0 )
            T1 = DSECND( )
            CALL ZTRSM( 'Left', 'Upper', 'Conjugate transpose',
     $                  'Non-unit', N, NRHS, CONE, A, LDA, B, LDB )
            T2 = DSECND( )
            TIMNG( TRSM ) = TIMNG( TRSM ) + ( T2-T1 )
*
*           B(N+1:M,1:NRHS) = ZERO
*
            DO 20 J = 1, NRHS
               DO 10 I = N + 1, M
                  B( I, J ) = CZERO
   10          CONTINUE
   20       CONTINUE
*
*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
*
            NB = ILAENV( 1, 'ZUNMQR', 'LN', M, NRHS, N, -1 )
            OPCNT( UNMQR ) = OPCNT( UNMQR ) +
     $                       DOPLA( 'ZUNMQR', M, NRHS, N, 0, NB )
            T1 = DSECND( )
            CALL ZUNMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA,
     $                   WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $                   INFO )
            T2 = DSECND( )
            TIMNG( UNMQR ) = TIMNG( UNMQR ) + ( T2-T1 )
*
*           workspace at least NRHS, optimally NRHS*NB
*
            SCLLEN = M
*
         END IF
*
      ELSE
*
*        Compute LQ factorization of A
*
         NB = ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
         OPCNT( GELQF ) = OPCNT( GELQF ) +
     $                    DOPLA( 'ZGELQF', M, N, 0, 0, NB )
         T1 = DSECND( )
         CALL ZGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN,
     $                INFO )
         T2 = DSECND( )
         TIMNG( GELQF ) = TIMNG( GELQF ) + ( T2-T1 )
*
*        workspace at least M, optimally M*NB.
*
         IF( .NOT.TPSD ) THEN
*
*           underdetermined system of equations A * X = B
*
*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
*
            OPCNT( TRSM ) = OPCNT( TRSM ) +
     $                      DOPBL3( 'ZTRSM ', M, NRHS, 0 )
            T1 = DSECND( )
            CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', M,
     $                  NRHS, CONE, A, LDA, B, LDB )
            T2 = DSECND( )
            TIMNG( TRSM ) = TIMNG( TRSM ) + ( T2-T1 )
*
*           B(M+1:N,1:NRHS) = 0
*
            DO 40 J = 1, NRHS
               DO 30 I = M + 1, N
                  B( I, J ) = CZERO
   30          CONTINUE
   40       CONTINUE
*
*           B(1:N,1:NRHS) := Q(1:N,:)' * B(1:M,1:NRHS)
*
            NB = ILAENV( 1, 'ZUNMLQ', 'LC', N, NRHS, M, -1 )
            OPCNT( UNMLQ ) = OPCNT( UNMLQ ) +
     $                       DOPLA( 'ZUNMLQ', N, NRHS, M, 0, NB )
            T1 = DSECND( )
            CALL ZUNMLQ( 'Left', 'Conjugate transpose', N, NRHS, M, A,
     $                   LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $                   INFO )
            T2 = DSECND( )
            TIMNG( UNMLQ ) = TIMNG( UNMLQ ) + ( T2-T1 )
*
*           workspace at least NRHS, optimally NRHS*NB
*
            SCLLEN = N
*
         ELSE
*
*           overdetermined system min || A' * X - B ||
*
*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
*
            NB = ILAENV( 1, 'ZUNMLQ', 'LN', N, NRHS, M, -1 )
            OPCNT( UNMLQ ) = OPCNT( UNMLQ ) +
     $                       DOPLA( 'ZUNMLQ', N, NRHS, M, 0, NB )
            T1 = DSECND( )
            CALL ZUNMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA,
     $                   WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $                   INFO )
            T2 = DSECND( )
            TIMNG( UNMLQ ) = TIMNG( UNMLQ ) + ( T2-T1 )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(1:M,1:NRHS) := inv(L') * B(1:M,1:NRHS)
*
            OPCNT( TRSM ) = OPCNT( TRSM ) +
     $                      DOPBL3( 'ZTRSM ', M, NRHS, 0 )
            T1 = DSECND( )
            CALL ZTRSM( 'Left', 'Lower', 'Conjugate transpose',
     $                  'Non-unit', M, NRHS, CONE, A, LDA, B, LDB )
            T2 = DSECND( )
            TIMNG( TRSM ) = TIMNG( TRSM ) + ( T2-T1 )
*
            SCLLEN = M
*
         END IF
*
      END IF
*
*     Undo scaling
*
      IF( IASCL.EQ.1 ) THEN
         OPCNT( GELS ) = OPCNT( GELS ) + DBLE( 6*SCLLEN*NRHS )
         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         OPCNT( GELS ) = OPCNT( GELS ) + DBLE( 6*SCLLEN*NRHS )
         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         OPCNT( GELS ) = OPCNT( GELS ) + DBLE( 6*SCLLEN*NRHS )
         CALL ZLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         OPCNT( GELS ) = OPCNT( GELS ) + DBLE( 6*SCLLEN*NRHS )
         CALL ZLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      END IF
*
   50 CONTINUE
      WORK( 1 ) = DBLE( WSIZE )
*
      RETURN
*
*     End of ZGELS
*
      END
