*> \brief \b DGETSLS
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB,
*     $                     WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGETSLS solves overdetermined or underdetermined real linear systems
*> involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
*> factorization of A.
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
*> 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
*>    an undetermined system A**T * X = B.
*>
*> 4. If TRANS = 'T' and m < n:  find the least squares solution of
*>    an overdetermined system, i.e., solve the least squares problem
*>                 minimize || B - A**T * X ||.
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit,
*>          A is overwritten by details of its QR or LQ
*>          factorization as returned by DGEQR or DGELQ.
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
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
*>          (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) contains optimal (or either minimal
*>          or optimal, if query was assumed) LWORK.
*>          See LWORK for details.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= 1.
*>          If LWORK = -1 or -2, then a workspace query is assumed.
*>          If LWORK = -1, the routine calculates optimal size of WORK for the
*>          optimal performance and returns this value in WORK(1).
*>          If LWORK = -2, the routine calculates minimal size of WORK and 
*>          returns this value in WORK(1).
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
*> \ingroup getsls
*
*  =====================================================================
      SUBROUTINE DGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB,
     $                    WORK, LWORK, INFO )
      IMPLICIT NONE
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
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, TRAN
      INTEGER            I, IASCL, IBSCL, J, MAXMN, BROW,
     $                   SCLLEN, TSZO, TSZM, LWO, LWM, LW1, LW2,
     $                   WSIZEO, WSIZEM, INFO2
      DOUBLE PRECISION   ANRM, BIGNUM, BNRM, SMLNUM, TQ( 5 ), WORKQ( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQR, DGEMQR, DLASCL, DLASET,
     $                   DTRTRS, XERBLA, DGELQ, DGEMLQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN, INT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      MAXMN = MAX( M, N )
      TRAN  = LSAME( TRANS, 'T' )
*
      LQUERY = ( LWORK.EQ.-1 .OR. LWORK.EQ.-2 )
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
      IF( INFO.EQ.0 ) THEN
*
*     Determine the optimum and minimum LWORK
*
       IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         WSIZEM = 1
         WSIZEO = 1
       ELSE IF( M.GE.N ) THEN
         CALL DGEQR( M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 )
         TSZO = INT( TQ( 1 ) )
         LWO  = INT( WORKQ( 1 ) )
         CALL DGEMQR( 'L', TRANS, M, NRHS, N, A, LDA, TQ,
     $                TSZO, B, LDB, WORKQ, -1, INFO2 )
         LWO  = MAX( LWO, INT( WORKQ( 1 ) ) )
         CALL DGEQR( M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 )
         TSZM = INT( TQ( 1 ) )
         LWM  = INT( WORKQ( 1 ) )
         CALL DGEMQR( 'L', TRANS, M, NRHS, N, A, LDA, TQ,
     $                TSZM, B, LDB, WORKQ, -1, INFO2 )
         LWM = MAX( LWM, INT( WORKQ( 1 ) ) )
         WSIZEO = TSZO + LWO
         WSIZEM = TSZM + LWM
       ELSE
         CALL DGELQ( M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 )
         TSZO = INT( TQ( 1 ) )
         LWO  = INT( WORKQ( 1 ) )
         CALL DGEMLQ( 'L', TRANS, N, NRHS, M, A, LDA, TQ,
     $                TSZO, B, LDB, WORKQ, -1, INFO2 )
         LWO  = MAX( LWO, INT( WORKQ( 1 ) ) )
         CALL DGELQ( M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 )
         TSZM = INT( TQ( 1 ) )
         LWM  = INT( WORKQ( 1 ) )
         CALL DGEMLQ( 'L', TRANS, N, NRHS, M, A, LDA, TQ,
     $                TSZM, B, LDB, WORKQ, -1, INFO2 )
         LWM  = MAX( LWM, INT( WORKQ( 1 ) ) )
         WSIZEO = TSZO + LWO
         WSIZEM = TSZM + LWM
       END IF
*
       IF( ( LWORK.LT.WSIZEM ).AND.( .NOT.LQUERY ) ) THEN
          INFO = -10
       END IF
*
       WORK( 1 ) = DBLE( WSIZEO )
*
      END IF
*
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGETSLS', -INFO )
        RETURN
      END IF
      IF( LQUERY ) THEN
        IF( LWORK.EQ.-2 ) WORK( 1 ) = DBLE( WSIZEM )
        RETURN
      END IF
      IF( LWORK.LT.WSIZEO ) THEN
        LW1 = TSZM
        LW2 = LWM
      ELSE
        LW1 = TSZO
        LW2 = LWO
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
           CALL DLASET( 'FULL', MAX( M, N ), NRHS, ZERO, ZERO,
     $                  B, LDB )
           RETURN
      END IF
*
*     Get machine parameters
*
       SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
       BIGNUM = ONE / SMLNUM
*
*     Scale A, B if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = DLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL DLASET( 'F', MAXMN, NRHS, ZERO, ZERO, B, LDB )
         GO TO 50
      END IF
*
      BROW = M
      IF ( TRAN ) THEN
        BROW = N
      END IF
      BNRM = DLANGE( 'M', BROW, NRHS, B, LDB, WORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB,
     $                INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB,
     $                INFO )
         IBSCL = 2
      END IF
*
      IF ( M.GE.N ) THEN
*
*        compute QR factorization of A
*
        CALL DGEQR( M, N, A, LDA, WORK( LW2+1 ), LW1,
     $              WORK( 1 ), LW2, INFO )
        IF ( .NOT.TRAN ) THEN
*
*           Least-Squares Problem min || A * X - B ||
*
*           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
*
          CALL DGEMQR( 'L' , 'T', M, NRHS, N, A, LDA,
     $                 WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2,
     $                 INFO )
*
*           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
*
          CALL DTRTRS( 'U', 'N', 'N', N, NRHS,
     $                  A, LDA, B, LDB, INFO )
          IF( INFO.GT.0 ) THEN
            RETURN
          END IF
          SCLLEN = N
        ELSE
*
*           Overdetermined system of equations A**T * X = B
*
*           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
*
            CALL DTRTRS( 'U', 'T', 'N', N, NRHS,
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
            CALL DGEMQR( 'L', 'N', M, NRHS, N, A, LDA,
     $                   WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2,
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
         CALL DGELQ( M, N, A, LDA, WORK( LW2+1 ), LW1,
     $               WORK( 1 ), LW2, INFO )
*
*        workspace at least M, optimally M*NB.
*
         IF( .NOT.TRAN ) THEN
*
*           underdetermined system of equations A * X = B
*
*           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
*
            CALL DTRTRS( 'L', 'N', 'N', M, NRHS,
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
            CALL DGEMLQ( 'L', 'T', N, NRHS, M, A, LDA,
     $                   WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2,
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
            CALL DGEMLQ( 'L', 'N', N, NRHS, M, A, LDA,
     $                   WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2,
     $                   INFO )
*
*           workspace at least NRHS, optimally NRHS*NB
*
*           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
*
            CALL DTRTRS( 'Lower', 'Transpose', 'Non-unit', M, NRHS,
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
        CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB,
     $               INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
        CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB,
     $               INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
        CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $               INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
        CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB,
     $               INFO )
      END IF
*
   50 CONTINUE
      WORK( 1 ) = DBLE( TSZO + LWO )
      RETURN
*
*     End of DGETSLS
*
      END
