*  Definition:
*  ===========
*
*       SUBROUTINE SGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB
*     $                   , WORK, LWORK, INFO )

*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*       ..
*       .. Array Arguments ..
*       REAL   A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGETSLS solves overdetermined or underdetermined real linear systems
*> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
*> factorization of A.  It is assumed that A has full rank.
*>
*>
*>
*> The following options are provided:
*>
*> 1. If TRANS = 'N' and m >= n:  find the least squares solution of
*>    an overdetermined system, i.e., solve the least squares problem
*>                 minimize || B - A*X ||.

*> 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
*>    an underdetermined system A * X = B.

*> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
*>    an undetermined system A**T * X = B.

*> 4. If TRANS = 'T' and m < n:  find the least squares solution of
*>    an overdetermined system, i.e., solve the least squares problem
*>                 minimize || B - A**T * X ||.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': the linear system involves A;
*>          = 'T': the linear system involves A**T.
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
*>          columns of the matrices B and X. NRHS >=0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit,
*>          A is overwritten by details of its QR or LQ
*>          factorization as returned by DGETSQR.
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
*>          B is REAL array, dimension (LDB,NRHS)
*>          On entry, the matrix B of right hand side vectors, stored
*>          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
*>          if TRANS = 'T'.
*>          On exit, if INFO = 0, B is overwritten by the solution
*>          vectors, stored columnwise:
*>          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
*>          squares solution vectors.
*>          if TRANS = 'N' and m < n, rows 1 to N of B contain the
*>          minimum norm solution vectors;
*>          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
*>          minimum norm solution vectors;
*>          if TRANS = 'T' and m < n, rows 1 to M of B contain the
*>          least squares solution vectors.
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK,
*>          and WORK(2) returns the minimum LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          IF LWORK=-1, workspace query is assumed, and
*>          WORK(1) returns the optimal LWORK,
*>          and WORK(2) returns the minimum LWORK.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO =  i, the i-th diagonal element of the
*>                triangular factor of A is zero, so that A does not have
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
*> \date November 2011
*
*> \ingroup doubleGEsolve
*
*  =====================================================================
      SUBROUTINE SGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB
     $                   , WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, MB
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
*
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, TRAN
      INTEGER            I, IASCL, IBSCL, J, MINMN, MAXMN, BROW, LW,
     $                   SCLLEN, MNK, WSIZEO, WSIZEM, LW1, LW2, INFO2,
     $                   NB
      REAL               ANRM, BIGNUM, BNRM, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANGE
      EXTERNAL           LSAME, ILAENV, SLABAD, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEQR, SGEMQR, SLASCL, SLASET,
     $                   STRTRS, XERBLA, SGELQ, SGEMLQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO=0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      MNK   = MAX(MINMN,NRHS)
      TRAN  = LSAME( TRANS, 'T' )
*
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( LSAME( TRANS, 'N' ) .OR.
     $    LSAME( TRANS, 'T' ) ) ) THEN
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
      END IF
*
      IF( INFO.EQ.0)  THEN
*
*     Determine the block size and minimum LWORK
*
       IF ( M.GE.N ) THEN
        CALL SGEQR( M, N, A, LDA, WORK(1), -1, WORK(6), -1,
     $   INFO2)
        MB = INT(WORK(4))
        NB = INT(WORK(5))
        LW = INT(WORK(6))
        CALL SGEMQR( 'L', TRANS, M, NRHS, N, A, LDA, WORK(1),
     $        INT(WORK(2)), B, LDB, WORK(6), -1 , INFO2 )
        WSIZEO = INT(WORK(2))+MAX(LW,INT(WORK(6)))
        WSIZEM = INT(WORK(3))+MAX(LW,INT(WORK(6)))
       ELSE
        CALL SGELQ( M, N, A, LDA, WORK(1), -1, WORK(6), -1,
     $   INFO2)
        MB = INT(WORK(4))
        NB = INT(WORK(5))
        LW = INT(WORK(6))
        CALL SGEMLQ( 'L', TRANS, N, NRHS, M, A, LDA, WORK(1),
     $        INT(WORK(2)), B, LDB, WORK(6), -1 , INFO2 )
        WSIZEO = INT(WORK(2))+MAX(LW,INT(WORK(6)))
        WSIZEM = INT(WORK(3))+MAX(LW,INT(WORK(6)))
       END IF
*
       IF((LWORK.LT.WSIZEO).AND.(.NOT.LQUERY)) THEN
          INFO=-10
       END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'SGETSLS', -INFO )
        WORK( 1 ) = REAL( WSIZEO )
        WORK( 2 ) = REAL( WSIZEM )
        RETURN
      ELSE IF (LQUERY) THEN
        WORK( 1 ) = REAL( WSIZEO )
        WORK( 2 ) = REAL( WSIZEM )
        RETURN
      END IF
      IF(LWORK.LT.WSIZEO) THEN
        LW1=INT(WORK(3))
        LW2=MAX(LW,INT(WORK(6)))
      ELSE
        LW1=INT(WORK(2))
        LW2=MAX(LW,INT(WORK(6)))
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
           CALL SLASET( 'FULL', MAX( M, N ), NRHS, ZERO, ZERO,
     $       B, LDB )
           RETURN
      END IF
*
*     Get machine parameters
*
       SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
       BIGNUM = ONE / SMLNUM
       CALL SLABAD( SMLNUM, BIGNUM )
*
*     Scale A, B if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = SLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL SLASET( 'F', MAXMN, NRHS, ZERO, ZERO, B, LDB )
         GO TO 50
      END IF
*
      BROW = M
      IF ( TRAN ) THEN
        BROW = N
      END IF
      BNRM = SLANGE( 'M', BROW, NRHS, B, LDB, WORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL SLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB,
     $                INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL SLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB,
     $                INFO )
         IBSCL = 2
      END IF
*
      IF ( M.GE.N) THEN
*
*        compute QR factorization of A
*
        CALL SGEQR( M, N, A, LDA, WORK(LW2+1), LW1
     $    , WORK(1), LW2, INFO )
        IF (.NOT.TRAN) THEN
*
*           Least-Squares Problem min || A * X - B ||
*
*           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
*
          CALL SGEMQR( 'L' , 'T', M, NRHS, N, A, LDA,
     $         WORK(LW2+1), LW1, B, LDB, WORK(1), LW2, INFO )
*
*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
*
          CALL STRTRS( 'U', 'N', 'N', N, NRHS,
     $                   A, LDA, B, LDB, INFO )
          IF(INFO.GT.0) THEN
            RETURN
          END IF
          SCLLEN = N
        ELSE
*
*           Overdetermined system of equations A**T * X = B
*
*           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
*
            CALL STRTRS( 'U', 'T', 'N', N, NRHS,
     $                   A, LDA, B, LDB, INFO )
*
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
*
*           B(N+1:M,1:NRHS) = ZERO
*
            DO 20 J = 1, NRHS
               DO 10 I = N + 1, M
                  B( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
*
*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
*
            CALL SGEMQR( 'L', 'N', M, NRHS, N, A, LDA,
     $                   WORK( LW2+1), LW1, B, LDB, WORK( 1 ), LW2,
     $                   INFO )
*
            SCLLEN = M
*
         END IF
*
      ELSE
*
*        Compute LQ factorization of A
*
        CALL SGELQ( M, N, A, LDA, WORK(LW2+1), LW1
     $    , WORK(1), LW2, INFO )
*
*        workspace at least M, optimally M*NB.
*
         IF( .NOT.TRAN ) THEN
*
*           underdetermined system of equations A * X = B
*
*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
*
            CALL STRTRS( 'L', 'N', 'N', M, NRHS,
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
                  B( I, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
*
            CALL SGEMLQ( 'L', 'T', N, NRHS, M, A, LDA,
     $                   WORK( LW2+1), LW1, B, LDB, WORK( 1 ), LW2,
     $                   INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
            SCLLEN = N
*
         ELSE
*
*           overdetermined system min || A**T * X - B ||
*
*           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
*
            CALL SGEMLQ( 'L', 'N', N, NRHS, M, A, LDA,
     $                   WORK( LW2+1), LW1, B, LDB, WORK( 1 ), LW2,
     $                   INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
*
            CALL STRTRS( 'Lower', 'Transpose', 'Non-unit', M, NRHS,
     $                   A, LDA, B, LDB, INFO )
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
        CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
        CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $                INFO )
      END IF
*
   50 CONTINUE
       WORK( 1 ) = REAL( WSIZEO )
       WORK( 2 ) = REAL( WSIZEM )
      RETURN
*
*     End of SGETSLS
*
      END
