*> \brief <b> CGELS solves overdetermined or underdetermined systems for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CGELS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgels.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgels.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgels.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
*                         INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*       ..
*       .. Array Arguments ..
*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CGELS solves overdetermined or underdetermined complex linear systems
*> involving an M-by-N matrix A, or its conjugate-transpose, using a QR
*> or LQ factorization of A.
*>
*> It is assumed that A has full rank, and only a rudimentary protection
*> against rank-deficient matrices is provided. This subroutine only detects
*> exact rank-deficiency, where a diagonal element of the triangular factor
*> of A is exactly zero.
*>
*> It is conceivable for one (or more) of the diagonal elements of the triangular
*> factor of A to be subnormally tiny numbers without this subroutine signalling
*> an error. The solutions computed for such almost-rank-deficient matrices may
*> be less accurate due to a loss of numerical precision.
*>
*> The following options are provided:
*>
*> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
*>    an overdetermined system, i.e., solve the least squares problem
*>                 minimize || B - A*X ||.
*>
*> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
*>    an underdetermined system A * X = B.
*>
*> 3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
*>    an underdetermined system A**H * X = B.
*>
*> 4. If TRANS = 'C' and m < n:  find the least squares solution of
*>    an overdetermined system, i.e., solve the least squares problem
*>                 minimize || B - A**H * X ||.
*>
*> Several right hand side vectors b and solution vectors x can be
*> handled in a single call; they are stored as the columns of the
*> M-by-NRHS right hand side matrix B and the N-by-NRHS solution
*> matrix X.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': the linear system involves A;
*>          = 'C': the linear system involves A**H.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of
*>          columns of the matrices B and X. NRHS >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>            if M >= N, A is overwritten by details of its QR
*>                       factorization as returned by CGEQRF;
*>            if M <  N, A is overwritten by details of its LQ
*>                       factorization as returned by CGELQF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,NRHS)
*>          On entry, the matrix B of right hand side vectors, stored
*>          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
*>          if TRANS = 'C'.
*>          On exit, if INFO = 0, B is overwritten by the solution
*>          vectors, stored columnwise:
*>          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
*>          squares solution vectors; the residual sum of squares for the
*>          solution in each column is given by the sum of squares of the
*>          modulus of elements N+1 to M in that column;
*>          if TRANS = 'N' and m < n, rows 1 to N of B contain the
*>          minimum norm solution vectors;
*>          if TRANS = 'C' and m >= n, rows 1 to M of B contain the
*>          minimum norm solution vectors;
*>          if TRANS = 'C' and m < n, rows 1 to M of B contain the
*>          least squares solution vectors; the residual sum of squares
*>          for the solution in each column is given by the sum of
*>          squares of the modulus of elements M+1 to N in that column.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= MAX(1,M,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          LWORK >= max( 1, MN + max( MN, NRHS ) ).
*>          For optimal performance,
*>          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
*>          where MN = min(M,N) and NB is the optimum block size.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO =  i, the i-th diagonal element of the
*>                triangular factor of A is exactly zero, so that A does not have
*>                full rank; the least squares solution could not be
*>                computed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup gels
*
*  =====================================================================
      SUBROUTINE CGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK,
     $                  LWORK,
     $                  INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, TPSD
      INTEGER            BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE
      REAL               ANRM, BIGNUM, BNRM, SMLNUM
*     ..
*     .. Local Arrays ..
      REAL               RWORK( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               CLANGE, SLAMCH
      EXTERNAL           LSAME, ILAENV, CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGELQF, CGEQRF, CLASCL, CLASET, CTRTRS,
     $                   CUNMLQ,
     $                   CUNMQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( LSAME( TRANS, 'N' ) .OR.
     $    LSAME( TRANS, 'C' ) ) ) THEN
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
      ELSE IF( LWORK.LT.MAX( 1, MN+MAX( MN, NRHS ) ) .AND.
     $   .NOT.LQUERY ) THEN
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
            NB = ILAENV( 1, 'CGEQRF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'CUNMQR', 'LN', M, NRHS, N,
     $              -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'CUNMQR', 'LC', M, NRHS, N,
     $              -1 ) )
            END IF
         ELSE
            NB = ILAENV( 1, 'CGELQF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'CUNMLQ', 'LC', N, NRHS, M,
     $              -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'CUNMLQ', 'LN', N, NRHS, M,
     $              -1 ) )
            END IF
         END IF
*
         WSIZE = MAX( 1, MN + MAX( MN, NRHS )*NB )
         WORK( 1 ) = REAL( WSIZE )
*
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGELS ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         CALL CLASET( 'Full', MAX( M, N ), NRHS, CZERO, CZERO, B,
     $                LDB )
         RETURN
      END IF
*
*     Get machine parameters
*
      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
*
*     Scale A, B if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = CLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL CLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         GO TO 50
      END IF
*
      BROW = M
      IF( TPSD )
     $   BROW = N
      BNRM = CLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL CLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB,
     $                INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL CLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB,
     $                INFO )
         IBSCL = 2
      END IF
*
      IF( M.GE.N ) THEN
*
*        compute QR factorization of A
*
         CALL CGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ),
     $                LWORK-MN,
     $                INFO )
*
*        workspace at least N, optimally N*NB
*
         IF( .NOT.TPSD ) THEN
*
*           Least-Squares Problem min || A * X - B ||
*
*           B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)
*
            CALL CUNMQR( 'Left', 'Conjugate transpose', M, NRHS, N,
     $                   A,
     $                   LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $                   INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
*
            CALL CTRTRS( 'Upper', 'No transpose', 'Non-unit', N,
     $                   NRHS,
     $                   A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*
            SCLLEN = N
*
         ELSE
*
*           Underdetermined system of equations A**T * X = B
*
*           B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS)
*
            CALL CTRTRS( 'Upper', 'Conjugate transpose','Non-unit',
     $                   N, NRHS, A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
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
            CALL CUNMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA,
     $                   WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $                   INFO )
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
         CALL CGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ),
     $                LWORK-MN,
     $                INFO )
*
*        workspace at least M, optimally M*NB.
*
         IF( .NOT.TPSD ) THEN
*
*           underdetermined system of equations A * X = B
*
*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
*
            CALL CTRTRS( 'Lower', 'No transpose', 'Non-unit', M,
     $                   NRHS,
     $                   A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*
*           B(M+1:N,1:NRHS) = 0
*
            DO 40 J = 1, NRHS
               DO 30 I = M + 1, N
                  B( I, J ) = CZERO
   30          CONTINUE
   40       CONTINUE
*
*           B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS)
*
            CALL CUNMLQ( 'Left', 'Conjugate transpose', N, NRHS, M,
     $                   A,
     $                   LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $                   INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
            SCLLEN = N
*
         ELSE
*
*           overdetermined system min || A**H * X - B ||
*
*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
*
            CALL CUNMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA,
     $                   WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,
     $                   INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS)
*
            CALL CTRTRS( 'Lower', 'Conjugate transpose', 'Non-unit',
     $                   M, NRHS, A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
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
         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL CLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      END IF
*
   50 CONTINUE
      WORK( 1 ) = REAL( WSIZE )
*
      RETURN
*
*     End of CGELS
*
      END
