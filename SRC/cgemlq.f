*> \brief \b CGEMLQ
*
*  Definition:
*  ===========
*
*      SUBROUTINE CGEMLQ( SIDE, TRANS, M, N, K, A, LDA, T,
*     $                   TSIZE, C, LDC, WORK, LWORK, INFO )
*
*
*     .. Scalar Arguments ..
*     CHARACTER         SIDE, TRANS
*     INTEGER           INFO, LDA, M, N, K, LDT, TSIZE, LWORK, LDC
*     ..
*     .. Array Arguments ..
*     COMPLEX           A( LDA, * ), T( * ), C(LDC, * ), WORK( * )
*     ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>     CGEMLQ overwrites the general real M-by-N matrix C with
*>
*>                      SIDE = 'L'     SIDE = 'R'
*>      TRANS = 'N':      Q * C          C * Q
*>      TRANS = 'C':      Q**H * C       C * Q**H
*>      where Q is a complex unitary matrix defined as the product
*>      of blocked elementary reflectors computed by short wide
*>      LQ factorization (CGELQ)
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**H from the Left;
*>          = 'R': apply Q or Q**H from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'C':  Conjugate transpose, apply Q**H.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >=0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension
*>                               (LDA,M) if SIDE = 'L',
*>                               (LDA,N) if SIDE = 'R'
*>          Part of the data structure to represent Q as returned by CGELQ.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,K).
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is COMPLEX array, dimension (MAX(5,TSIZE)).
*>          Part of the data structure to represent Q as returned by CGELQ.
*> \endverbatim
*>
*> \param[in] TSIZE
*> \verbatim
*>          TSIZE is INTEGER
*>          The dimension of the array T. TSIZE >= 5.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          (workspace) COMPLEX array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the minimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= 1.
*>          If LWORK = -1, then a workspace query is assumed. The routine
*>          only calculates the size of the WORK array, returns this
*>          value as WORK(1), and no error message related to WORK
*>          is issued by XERBLA.
*> \endverbatim
*>
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
*> \par Further Details
*  ====================
*>
*> \verbatim
*>
*> These details are particular for this LAPACK implementation. Users should not
*> take them for granted. These details may change in the future, and are not likely
*> true for another LAPACK implementation. These details are relevant if one wants
*> to try to understand the code. They are not part of the interface.
*>
*> In this version,
*>
*>          T(2): row block size (MB)
*>          T(3): column block size (NB)
*>          T(6:TSIZE): data structure needed for Q, computed by
*>                           CLASWQR or CGELQT
*>
*>  Depending on the matrix dimensions M and N, and row and column
*>  block sizes MB and NB returned by ILAENV, CGELQ will use either
*>  CLASWLQ (if the matrix is wide-and-short) or CGELQT to compute
*>  the LQ factorization.
*>  This version of CGEMLQ will use either CLAMSWLQ or CGEMLQT to
*>  multiply matrix Q by another matrix.
*>  Further Details in CLAMSWLQ or CGEMLQT.
*> \endverbatim
*>
*> \ingroup gemlq
*>
*  =====================================================================
      SUBROUTINE CGEMLQ( SIDE, TRANS, M, N, K, A, LDA, T, TSIZE,
     $                   C, LDC, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, LDA, M, N, K, TSIZE, LWORK, LDC
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), T( * ), C( LDC, * ), WORK( * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN, LQUERY
      INTEGER            MB, NB, LW, NBLCKS, MN, MINMNK, LWMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SROUNDUP_LWORK
      EXTERNAL           LSAME, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLAMSWLQ, CGEMLQT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      LQUERY  = ( LWORK.EQ.-1 )
      NOTRAN  = LSAME( TRANS, 'N' )
      TRAN    = LSAME( TRANS, 'C' )
      LEFT    = LSAME( SIDE, 'L' )
      RIGHT   = LSAME( SIDE, 'R' )
*
      MB = INT( T( 2 ) )
      NB = INT( T( 3 ) )
      IF( LEFT ) THEN
        LW = N * MB
        MN = M
      ELSE
        LW = M * MB
        MN = N
      END IF
*
      MINMNK = MIN( M, N, K )
      IF( MINMNK.EQ.0 ) THEN
         LWMIN = 1
      ELSE
         LWMIN = MAX( 1, LW )
      END IF
*
      IF( ( NB.GT.K ) .AND. ( MN.GT.K ) ) THEN
        IF( MOD( MN - K, NB - K ) .EQ. 0 ) THEN
          NBLCKS = ( MN - K ) / ( NB - K )
        ELSE
          NBLCKS = ( MN - K ) / ( NB - K ) + 1
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
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.MN ) THEN
        INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
        INFO = -7
      ELSE IF( TSIZE.LT.5 ) THEN
        INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -11
      ELSE IF( ( LWORK.LT.LWMIN ) .AND. ( .NOT.LQUERY ) ) THEN
        INFO = -13
      END IF
*
      IF( INFO.EQ.0 ) THEN
        WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      END IF
*
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'CGEMLQ', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( MINMNK.EQ.0 ) THEN
        RETURN
      END IF
*
      IF( ( LEFT .AND. M.LE.K ) .OR. ( RIGHT .AND. N.LE.K )
     $     .OR. ( NB.LE.K ) .OR. ( NB.GE.MAX( M, N, K ) ) ) THEN
        CALL CGEMLQT( SIDE, TRANS, M, N, K, MB, A, LDA,
     $                T( 6 ), MB, C, LDC, WORK, INFO )
      ELSE
        CALL CLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T( 6 ),
     $                 MB, C, LDC, WORK, LWORK, INFO )
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
*
      RETURN
*
*     End of CGEMLQ
*
      END
