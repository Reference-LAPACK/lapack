*
*  Definition:
*  ===========
*
*      SUBROUTINE DGEMLQ( SIDE, TRANS, M, N, K, A, LDA, WORK1,
*     $                LWORK1, C, LDC, WORK2, LWORK2, INFO )
*
*
*     .. Scalar Arguments ..
*      CHARACTER         SIDE, TRANS
*      INTEGER           INFO, LDA, M, N, K, MB, NB, LWORK1, LWORK2, LDC
*     ..
*     .. Array Arguments ..
*      DOUBLE        A( LDA, * ), WORK1( * ), C(LDC, * ),
*     $                  WORK2( * )
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>     DGEMLQ overwrites the general real M-by-N matrix C with
*>
*>
*>                    SIDE = 'L'     SIDE = 'R'
*>    TRANS = 'N':      Q * C          C * Q
*>    TRANS = 'T':      Q**T * C       C * Q**T
*>    where Q is a real orthogonal matrix defined as the product
*>    of blocked elementary reflectors computed by short wide LQ
*>    factorization (DGELQ)
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**T from the Left;
*>          = 'R': apply Q or Q**T from the Right.
*>
*> \param[in] TRANS
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'T':  Transpose, apply Q**T.
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >=0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= M.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          M >= K >= 0;
*>
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,K)
*>          The i-th row must contain the vector which defines the blocked
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          DLASWLQ in the first k rows of its array argument A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If SIDE = 'L', LDA >= max(1,M);
*>          if SIDE = 'R', LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] WORK1
*> \verbatim
*>          WORK1 is DOUBLE PRECISION array, dimension (MAX(1,LWORK1)) is
*>          returned by GEQR.
*> \endverbatim
*>
*> \param[in] LWORK1
*> \verbatim
*>          LWORK1 is INTEGER
*>          The dimension of the array WORK1.
*> \endverbatim
*>
*> \param[in,out] C
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*> \param[in] LDC
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*>
*> \param[out] WORK2
*> \verbatim
*>         (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK2))
*>
*> \endverbatim
*> \param[in] LWORK2
*> \verbatim
*>          LWORK2 is INTEGER
*>          The dimension of the array WORK2.
*>          If LWORK2 = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK2 array, returns
*>          this value as the third entry of the WORK2 array (WORK2(1)),
*>          and no error message related to LWORK2 is issued by XERBLA.
*>
*> \endverbatim
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>  Depending on the matrix dimensions M and N, and row and column
*>  block sizes MB and NB returned by ILAENV, GELQ will use either
*>  LASWLQ(if the matrix is short-and-wide) or GELQT to compute
*>  the LQ decomposition.
*>  The output of LASWLQ or GELQT representing Q is stored in A and in
*>  array WORK1(6:LWORK1) for later use.
*>  WORK1(2:5) contains the matrix dimensions M,N and block sizes MB, NB
*>  which are needed to interpret A and WORK1(6:LWORK1) for later use.
*>  WORK1(1)=1 indicates that the code needed to take WORK1(2:5) and
*>  decide whether LASWLQ or GELQT was used is the same as used below in
*>  GELQ. For a detailed description of A and WORK1(6:LWORK1), see
*>  Further Details in LASWLQ or GELQT.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEMLQ( SIDE, TRANS, M, N, K, A, LDA, WORK1, LWORK1,
     $      C, LDC, WORK2, LWORK2, INFO )
*
*  -- LAPACK computational routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2013
*
*     .. Scalar Arguments ..
      CHARACTER         SIDE, TRANS
      INTEGER           INFO, LDA, M, N, K, LWORK1, LWORK2, LDC
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A( LDA, * ), C( LDC, * ), WORK1( * ), WORK2( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL    LEFT, RIGHT, TRAN, NOTRAN, LQUERY
      INTEGER    I, II, KK, MB, NB, LW, NBLCKS, MN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           DTPMLQT, DGEMLQT, XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      LQUERY  = (LWORK2.LT.0)
      NOTRAN  = LSAME( TRANS, 'N' )
      TRAN    = LSAME( TRANS, 'T' )
      LEFT    = LSAME( SIDE, 'L' )
      RIGHT   = LSAME( SIDE, 'R' )
*
      MB = INT(WORK1(4))
      NB = INT(WORK1(5))
      IF (LEFT) THEN
        LW = N * MB
        MN = M
      ELSE
        LW = M * MB
        MN = N
      END IF
      IF ((NB.GT.K).AND.(MN.GT.K)) THEN
        IF(MOD(MN-K, NB-K).EQ.0) THEN
          NBLCKS = (MN-K)/(NB-K)
        ELSE
          NBLCKS = (MN-K)/(NB-K) + 1
        END IF
      ELSE
        NBLCKS = 1
      END IF
*
      INFO = 0
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
         INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
        INFO = -3
      ELSE IF( N.LT.0) THEN
        INFO = -4
      ELSE IF( K.LT.0 ) THEN
        INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
        INFO = -7
      ELSE IF( LWORK1.LT.MAX( 1, MB*K*NBLCKS+5 )) THEN
        INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF(( LWORK2.LT.MAX(1,LW)).AND.(.NOT.LQUERY)) THEN
        INFO = -13
      END IF
*
      IF( INFO.EQ.0)  THEN
        WORK2(1) = LW
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGEMLQ', -INFO )
        RETURN
      ELSE IF (LQUERY) THEN
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN(M,N,K).EQ.0 ) THEN
        RETURN
      END IF
*
      IF((LEFT.AND.M.LE.K).OR.(RIGHT.AND.N.LE.K).OR.(NB.LE.K).OR.
     $   (NB.GE.MAX(M,N,K))) THEN
        CALL DGEMLQT( SIDE, TRANS, M, N, K, MB, A, LDA,
     $        WORK1(6), MB, C, LDC, WORK2, INFO)
      ELSE
        CALL DLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, WORK1(6),
     $    MB, C, LDC, WORK2, LWORK2, INFO )
      END IF
*
      WORK2(1) = LW
*
      RETURN
*
*     End of DGEMLQ
*
      END
