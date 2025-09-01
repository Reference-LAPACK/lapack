*> \brief \b DLAMTSQR
*
*  Definition:
*  ===========
*
*      SUBROUTINE DLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T,
*     $                     LDT, C, LDC, WORK, LWORK, INFO )
*
*
*     .. Scalar Arguments ..
*      CHARACTER         SIDE, TRANS
*      INTEGER           INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC
*     ..
*     .. Array Arguments ..
*      DOUBLE        A( LDA, * ), WORK( * ), C(LDC, * ),
*     $                  T( LDT, * )
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>      DLAMTSQR overwrites the general real M-by-N matrix C with
*>
*>
*>                 SIDE = 'L'     SIDE = 'R'
*> TRANS = 'N':      Q * C          C * Q
*> TRANS = 'T':      Q**T * C       C * Q**T
*>      where Q is a real orthogonal matrix defined as the product
*>      of blocked elementary reflectors computed by tall skinny
*>      QR factorization (DLATSQR)
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**T from the Left;
*>          = 'R': apply Q or Q**T from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'T':  Transpose, apply Q**T.
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
*>          the matrix Q. M >= K >= 0;
*>
*> \endverbatim
*>
*> \param[in] MB
*> \verbatim
*>          MB is INTEGER
*>          The block size to be used in the blocked QR.
*>          MB > N. (must be the same as DLATSQR)
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The column block size to be used in the blocked QR.
*>          N >= NB >= 1.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,K)
*>          The i-th column must contain the vector which defines the
*>          blockedelementary reflector H(i), for i = 1,2,...,k, as
*>          returned by DLATSQR in the first k columns of
*>          its array argument A.
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
*> \param[in] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension
*>          ( N * Number of blocks(CEIL(M-K/MB-K)),
*>          The blocked upper triangular block reflectors stored in compact form
*>          as a sequence of upper triangular blocks.  See below
*>          for further details.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
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
*>          (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the minimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If MIN(M,N,K) = 0, LWORK >= 1.
*>          If SIDE = 'L', LWORK >= max(1,N*NB).
*>          If SIDE = 'R', LWORK >= max(1,MB*NB).
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the minimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
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
*> \par Further Details:
*  =====================
*>
*> \verbatim
*> Tall-Skinny QR (TSQR) performs QR by a sequence of orthogonal transformations,
*> representing Q as a product of other orthogonal matrices
*>   Q = Q(1) * Q(2) * . . . * Q(k)
*> where each Q(i) zeros out subdiagonal entries of a block of MB rows of A:
*>   Q(1) zeros out the subdiagonal entries of rows 1:MB of A
*>   Q(2) zeros out the bottom MB-N rows of rows [1:N,MB+1:2*MB-N] of A
*>   Q(3) zeros out the bottom MB-N rows of rows [1:N,2*MB-N+1:3*MB-2*N] of A
*>   . . .
*>
*> Q(1) is computed by GEQRT, which represents Q(1) by Householder vectors
*> stored under the diagonal of rows 1:MB of A, and by upper triangular
*> block reflectors, stored in array T(1:LDT,1:N).
*> For more information see Further Details in GEQRT.
*>
*> Q(i) for i>1 is computed by TPQRT, which represents Q(i) by Householder vectors
*> stored in rows [(i-1)*(MB-N)+N+1:i*(MB-N)+N] of A, and by upper triangular
*> block reflectors, stored in array T(1:LDT,(i-1)*N+1:i*N).
*> The last Q(k) may use fewer rows.
*> For more information see Further Details in TPQRT.
*>
*> For more details of the overall algorithm, see the description of
*> Sequential TSQR in Section 2.2 of [1].
*>
*> [1] “Communication-Optimal Parallel and Sequential QR and LU Factorizations,”
*>     J. Demmel, L. Grigori, M. Hoemmen, J. Langou,
*>     SIAM J. Sci. Comput, vol. 34, no. 1, 2012
*> \endverbatim
*>
*> \ingroup lamtsqr
*>
*  =====================================================================
      SUBROUTINE DLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T,
     $                     LDT, C, LDC, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * ), C( LDC, * ),
     $                   T( LDT, * )
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, RIGHT, TRAN, NOTRAN, LQUERY
      INTEGER            I, II, KK, LW, CTR, Q, MINMNK, LWMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           DGEMQRT, DTPMQRT, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY  = ( LWORK.EQ.-1 )
      NOTRAN  = LSAME( TRANS, 'N' )
      TRAN    = LSAME( TRANS, 'T' )
      LEFT    = LSAME( SIDE, 'L' )
      RIGHT   = LSAME( SIDE, 'R' )
      IF( LEFT ) THEN
        LW = N * NB
        Q = M
      ELSE
        LW = MB * NB
        Q = N
      END IF
*
      MINMNK = MIN( M, N, K )
      IF( MINMNK.EQ.0 ) THEN
        LWMIN = 1
      ELSE
        LWMIN = MAX( 1, LW )
      END IF
*
      IF( .NOT.LEFT .AND. .NOT.RIGHT ) THEN
        INFO = -1
      ELSE IF( .NOT.TRAN .AND. .NOT.NOTRAN ) THEN
        INFO = -2
      ELSE IF( M.LT.K ) THEN
        INFO = -3
      ELSE IF( N.LT.0 ) THEN
        INFO = -4
      ELSE IF( K.LT.0 ) THEN
        INFO = -5
      ELSE IF( K.LT.NB .OR. NB.LT.1 ) THEN
        INFO = -7
      ELSE IF( LDA.LT.MAX( 1, Q ) ) THEN
        INFO = -9
      ELSE IF( LDT.LT.MAX( 1, NB ) ) THEN
        INFO = -11
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
        INFO = -13
      ELSE IF( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) THEN
        INFO = -15
      END IF
*
      IF( INFO.EQ.0 ) THEN
        WORK( 1 ) = LWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DLAMTSQR', -INFO )
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
*     Determine the block size if it is tall skinny or short and wide
*
      IF((MB.LE.K).OR.(MB.GE.MAX(M,N,K))) THEN
        CALL DGEMQRT( SIDE, TRANS, M, N, K, NB, A, LDA,
     $        T, LDT, C, LDC, WORK, INFO )
        RETURN
      END IF
*
      IF(LEFT.AND.NOTRAN) THEN
*
*         Multiply Q to the last block of C
*
         KK = MOD((M-K),(MB-K))
         CTR = (M-K)/(MB-K)
         IF (KK.GT.0) THEN
           II=M-KK+1
           CALL DTPMQRT('L','N',KK , N, K, 0, NB, A(II,1), LDA,
     $       T(1,CTR*K+1),LDT , C(1,1), LDC,
     $       C(II,1), LDC, WORK, INFO )
         ELSE
           II=M+1
         END IF
*
         DO I=II-(MB-K),MB+1,-(MB-K)
*
*         Multiply Q to the current block of C (I:I+MB,1:N)
*
           CTR = CTR - 1
           CALL DTPMQRT('L','N',MB-K , N, K, 0,NB, A(I,1), LDA,
     $         T(1,CTR*K+1),LDT, C(1,1), LDC,
     $         C(I,1), LDC, WORK, INFO )
*
         END DO
*
*         Multiply Q to the first block of C (1:MB,1:N)
*
         CALL DGEMQRT('L','N',MB , N, K, NB, A(1,1), LDA, T
     $            ,LDT ,C(1,1), LDC, WORK, INFO )
*
      ELSE IF (LEFT.AND.TRAN) THEN
*
*         Multiply Q to the first block of C
*
         KK = MOD((M-K),(MB-K))
         II=M-KK+1
         CTR = 1
         CALL DGEMQRT('L','T',MB , N, K, NB, A(1,1), LDA, T
     $            ,LDT ,C(1,1), LDC, WORK, INFO )
*
         DO I=MB+1,II-MB+K,(MB-K)
*
*         Multiply Q to the current block of C (I:I+MB,1:N)
*
          CALL DTPMQRT('L','T',MB-K , N, K, 0,NB, A(I,1), LDA,
     $       T(1,CTR * K + 1),LDT, C(1,1), LDC,
     $       C(I,1), LDC, WORK, INFO )
          CTR = CTR + 1
*
         END DO
         IF(II.LE.M) THEN
*
*         Multiply Q to the last block of C
*
          CALL DTPMQRT('L','T',KK , N, K, 0,NB, A(II,1), LDA,
     $      T(1,CTR * K + 1), LDT, C(1,1), LDC,
     $      C(II,1), LDC, WORK, INFO )
*
         END IF
*
      ELSE IF(RIGHT.AND.TRAN) THEN
*
*         Multiply Q to the last block of C
*
          KK = MOD((N-K),(MB-K))
          CTR = (N-K)/(MB-K)
          IF (KK.GT.0) THEN
            II=N-KK+1
            CALL DTPMQRT('R','T',M , KK, K, 0, NB, A(II,1), LDA,
     $        T(1,CTR*K+1), LDT, C(1,1), LDC,
     $        C(1,II), LDC, WORK, INFO )
          ELSE
            II=N+1
          END IF
*
          DO I=II-(MB-K),MB+1,-(MB-K)
*
*         Multiply Q to the current block of C (1:M,I:I+MB)
*
            CTR = CTR - 1
            CALL DTPMQRT('R','T',M , MB-K, K, 0,NB, A(I,1), LDA,
     $          T(1,CTR*K+1), LDT, C(1,1), LDC,
     $          C(1,I), LDC, WORK, INFO )
*
          END DO
*
*         Multiply Q to the first block of C (1:M,1:MB)
*
          CALL DGEMQRT('R','T',M , MB, K, NB, A(1,1), LDA, T
     $              ,LDT ,C(1,1), LDC, WORK, INFO )
*
      ELSE IF (RIGHT.AND.NOTRAN) THEN
*
*         Multiply Q to the first block of C
*
         KK = MOD((N-K),(MB-K))
         II=N-KK+1
         CTR = 1
         CALL DGEMQRT('R','N', M, MB , K, NB, A(1,1), LDA, T
     $              ,LDT ,C(1,1), LDC, WORK, INFO )
*
         DO I=MB+1,II-MB+K,(MB-K)
*
*         Multiply Q to the current block of C (1:M,I:I+MB)
*
          CALL DTPMQRT('R','N', M, MB-K, K, 0,NB, A(I,1), LDA,
     $         T(1, CTR * K + 1),LDT, C(1,1), LDC,
     $         C(1,I), LDC, WORK, INFO )
          CTR = CTR + 1
*
         END DO
         IF(II.LE.N) THEN
*
*         Multiply Q to the last block of C
*
          CALL DTPMQRT('R','N', M, KK , K, 0,NB, A(II,1), LDA,
     $        T(1, CTR * K + 1),LDT, C(1,1), LDC,
     $        C(1,II), LDC, WORK, INFO )
*
         END IF
*
      END IF
*
      WORK( 1 ) = LWMIN
*
      RETURN
*
*     End of DLAMTSQR
*
      END
