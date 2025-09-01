*> \brief \b ZGELQ
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGELQ( M, N, A, LDA, T, TSIZE, WORK, LWORK,
*                         INFO )
*
*       .. Scalar Arguments ..
*       INTEGER           INFO, LDA, M, N, TSIZE, LWORK
*       ..
*       .. Array Arguments ..
*       COMPLEX*16       A( LDA, * ), T( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGELQ computes an LQ factorization of a complex M-by-N matrix A:
*>
*>    A = ( L 0 ) *  Q
*>
*> where:
*>
*>    Q is a N-by-N orthogonal matrix;
*>    L is a lower-triangular M-by-M matrix;
*>    0 is a M-by-(N-M) zero matrix, if M < N.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
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
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and below the diagonal of the array
*>          contain the M-by-min(M,N) lower trapezoidal matrix L
*>          (L is lower triangular if M <= N);
*>          the elements above the diagonal are used to store part of the 
*>          data structure to represent Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (MAX(5,TSIZE))
*>          On exit, if INFO = 0, T(1) returns optimal (or either minimal 
*>          or optimal, if query is assumed) TSIZE. See TSIZE for details.
*>          Remaining T contains part of the data structure used to represent Q.
*>          If one wants to apply or construct Q, then one needs to keep T 
*>          (in addition to A) and pass it to further subroutines.
*> \endverbatim
*>
*> \param[in] TSIZE
*> \verbatim
*>          TSIZE is INTEGER
*>          If TSIZE >= 5, the dimension of the array T.
*>          If TSIZE = -1 or -2, then a workspace query is assumed. The routine
*>          only calculates the sizes of the T and WORK arrays, returns these
*>          values as the first entries of the T and WORK arrays, and no error
*>          message related to T or WORK is issued by XERBLA.
*>          If TSIZE = -1, the routine calculates optimal size of T for the 
*>          optimum performance and returns this value in T(1).
*>          If TSIZE = -2, the routine calculates minimal size of T and 
*>          returns this value in T(1).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) contains optimal (or either minimal
*>          or optimal, if query was assumed) LWORK.
*>          See LWORK for details.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= 1.
*>          If LWORK = -1 or -2, then a workspace query is assumed. The routine
*>          only calculates the sizes of the T and WORK arrays, returns these
*>          values as the first entries of the T and WORK arrays, and no error
*>          message related to T or WORK is issued by XERBLA.
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
*> The goal of the interface is to give maximum freedom to the developers for
*> creating any LQ factorization algorithm they wish. The triangular 
*> (trapezoidal) L has to be stored in the lower part of A. The lower part of A
*> and the array T can be used to store any relevant information for applying or
*> constructing the Q factor. The WORK array can safely be discarded after exit.
*>
*> Caution: One should not expect the sizes of T and WORK to be the same from one 
*> LAPACK implementation to the other, or even from one execution to the other.
*> A workspace query (for T and WORK) is needed at each execution. However, 
*> for a given execution, the size of T and WORK are fixed and will not change 
*> from one query to the next.
*>
*> \endverbatim
*>
*> \par Further Details particular to this LAPACK implementation:
*  ==============================================================
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
*>                           ZLASWLQ or ZGELQT
*>
*>  Depending on the matrix dimensions M and N, and row and column
*>  block sizes MB and NB returned by ILAENV, ZGELQ will use either
*>  ZLASWLQ (if the matrix is short-and-wide) or ZGELQT to compute
*>  the LQ factorization.
*> \endverbatim
*>
*> \ingroup gelq
*>
*  =====================================================================
      SUBROUTINE ZGELQ( M, N, A, LDA, T, TSIZE, WORK, LWORK,
     $                  INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N, TSIZE, LWORK
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), T( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, LMINWS, MINT, MINW
      INTEGER            MB, NB, MINTSZ, NBLCKS, LWMIN, LWOPT, LWREQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGELQT, ZLASWLQ, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      LQUERY = ( TSIZE.EQ.-1 .OR. TSIZE.EQ.-2 .OR.
     $           LWORK.EQ.-1 .OR. LWORK.EQ.-2 )
*
      MINT = .FALSE.
      MINW = .FALSE.
      IF( TSIZE.EQ.-2 .OR. LWORK.EQ.-2 ) THEN
        IF( TSIZE.NE.-1 ) MINT = .TRUE.
        IF( LWORK.NE.-1 ) MINW = .TRUE.
      END IF
*
*     Determine the block size
*
      IF( MIN( M, N ).GT.0 ) THEN
        MB = ILAENV( 1, 'ZGELQ ', ' ', M, N, 1, -1 )
        NB = ILAENV( 1, 'ZGELQ ', ' ', M, N, 2, -1 )
      ELSE
        MB = 1
        NB = N
      END IF
      IF( MB.GT.MIN( M, N ) .OR. MB.LT.1 ) MB = 1
      IF( NB.GT.N .OR. NB.LE.M ) NB = N
      MINTSZ = M + 5
      IF ( NB.GT.M .AND. N.GT.M ) THEN
        IF( MOD( N - M, NB - M ).EQ.0 ) THEN
          NBLCKS = ( N - M ) / ( NB - M )
        ELSE
          NBLCKS = ( N - M ) / ( NB - M ) + 1
        END IF
      ELSE
        NBLCKS = 1
      END IF
*
*     Determine if the workspace size satisfies minimal size
*
      IF( ( N.LE.M ) .OR. ( NB.LE.M ) .OR. ( NB.GE.N ) ) THEN
         LWMIN = MAX( 1, N )
         LWOPT = MAX( 1, MB*N )
      ELSE
         LWMIN = MAX( 1, M )
         LWOPT = MAX( 1, MB*M )
      END IF
      LMINWS = .FALSE.
      IF( ( TSIZE.LT.MAX( 1, MB*M*NBLCKS + 5 ) .OR. LWORK.LT.LWOPT )
     $    .AND. ( LWORK.GE.LWMIN ) .AND. ( TSIZE.GE.MINTSZ )
     $    .AND. ( .NOT.LQUERY ) ) THEN
        IF( TSIZE.LT.MAX( 1, MB*M*NBLCKS + 5 ) ) THEN
            LMINWS = .TRUE.
            MB = 1
            NB = N
        END IF
        IF( LWORK.LT.LWOPT ) THEN
            LMINWS = .TRUE.
            MB = 1
        END IF
      END IF
      IF( ( N.LE.M ) .OR. ( NB.LE.M ) .OR. ( NB.GE.N ) ) THEN
         LWREQ = MAX( 1, MB*N )
      ELSE
         LWREQ = MAX( 1, MB*M )
      END IF
*
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      ELSE IF( TSIZE.LT.MAX( 1, MB*M*NBLCKS + 5 )
     $   .AND. ( .NOT.LQUERY ) .AND. ( .NOT.LMINWS ) ) THEN
        INFO = -6
      ELSE IF( ( LWORK.LT.LWREQ ) .AND .( .NOT.LQUERY )
     $   .AND. ( .NOT.LMINWS ) ) THEN
        INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
        IF( MINT ) THEN
          T( 1 ) = MINTSZ
        ELSE
          T( 1 ) = MB*M*NBLCKS + 5
        END IF
        T( 2 ) = MB
        T( 3 ) = NB
        IF( MINW ) THEN
          WORK( 1 ) = LWMIN
        ELSE
          WORK( 1 ) = LWREQ
        END IF
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'ZGELQ', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N ).EQ.0 ) THEN
        RETURN
      END IF
*
*     The LQ Decomposition
*
      IF( ( N.LE.M ) .OR. ( NB.LE.M ) .OR. ( NB.GE.N ) ) THEN
        CALL ZGELQT( M, N, MB, A, LDA, T( 6 ), MB, WORK, INFO )
      ELSE
        CALL ZLASWLQ( M, N, MB, NB, A, LDA, T( 6 ), MB, WORK,
     $                LWORK, INFO )
      END IF
*
      WORK( 1 ) = LWREQ
*
      RETURN
*
*     End of ZGELQ
*
      END
