*> \brief <b> DGGQRCS computes the singular value decomposition (SVD) for OTHER matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DGGQRCS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggqrcs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggqrcs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggqrcs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGGQRCS( JOBU1, JOBU2, JOBX, M, N, P, L, SWAPPED,
*                           A, LDA, B, LDB,
*                           ALPHA, BETA,
*                           U1, LDU1, U2, LDU2
*                           WORK, LWORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBU1, JOB2, JOBX
*       INTEGER            INFO, LDA, LDB, LDU1, LDU2, M, N, P, L, LWORK
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ),
*      $                   ALPHA( N ), BETA( N ),
*      $                   U1( LDU1, * ), U2( LDU2, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGGQRCS computes the generalized singular value decomposition (GSVD)
*> of an M-by-N real matrix A and P-by-N real matrix B:
*>
*>       A = U1 * D1 * X,           B = U2 * D2 * X
*>
*> where U1 and U2 are orthogonal matrices. DGGQRCS uses the QR
*> factorization with column pivoting and the 2-by-1 CS decomposition to
*> compute the GSVD.
*>
*> Let L be the effective numerical rank of the matrix (A**T,B**T)**T,
*> then X is a L-by-N nonsingular matrix, D1 and D2 are M-by-L and
*> P-by-L "diagonal" matrices. If SWAPPED is false, then D1 and D2 are
*> of the of the following structures, respectively:
*>
*>                 K1  K
*>            K1 [ I   0   0 ]
*>       D1 = K  [ 0   C   0 ]
*>               [ 0   0   0 ]
*>
*>                     K   K2
*>               [ 0   0   0 ]
*>       D2 = K  [ 0   S   0 ]
*>            K2 [ 0   0   I ]
*>
*> where
*>
*>   K  = MIN(M, P, L, M + P - L),
*>   K1 = MAX(L - P, 0),
*>   K2 = MAX(L - M, 0),
*>   C  = diag( ALPHA(1), ..., ALPHA(K) ),
*>   S  = diag( BETA(1), ..., BETA(K) ), and
*>   C^2 + S^2 = I.
*>
*> If SWAPPED is true, then D1 and D2 are of the of the following
*> structures, respectively:
*>
*>                     K   K1
*>               [ 0   0   0 ]
*>       D1 = K  [ 0   S   0 ]
*>            K1 [ 0   0   I ]
*>
*>                 K2  K
*>            K2 [ I   0   0 ]
*>       D2 = K  [ 0   C   0 ]
*>               [ 0   0   0 ]
*>
*> where
*>
*>   S  = diag( ALPHA(1), ..., ALPHA(K) ),
*>   C  = diag( BETA(1), ..., BETA(K) ), and
*>   C^2 + S^2 = I.
*>
*> The routine computes C, S and optionally the matrices U1, U2, and X.
*> On exit, X is stored in WORK( 2:L*N+1 ).
*>
*> If B is an N-by-N nonsingular matrix, then the GSVD of the matrix
*> pair (A, B) implicitly gives the SVD of A*inv(B):
*>
*>       A*inv(B) = U1*(D1*inv(D2))*U2**T.
*>
*> If (A**T,B**T)**T  has orthonormal columns, then the GSVD of A and B
*> is also equal to the CS decomposition of A and B. Furthermore, the
*> GSVD can be used to derive the solution of the eigenvalue problem:
*>
*>       A**T*A x = lambda * B**T*B x.
*>
*> In some literature, the GSVD of A and B is presented in the form
*>
*>       A = U1*D1*( 0 R )*Q**T,    B = U2*D2*( 0 R )*Q**T
*>
*> where U1, U2, and Q are orthogonal matrices. This latter GSVD form is
*> computed directly by DGGSVD3. It is possible to convert between the
*> two representations by calculating the RQ decomposition of X but this
*> is not recommended for reasons of numerical stability.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBU1
*> \verbatim
*>          JOBU1 is CHARACTER*1
*>          = 'Y':  Orthogonal matrix U1 is computed;
*>          = 'N':  U1 is not computed.
*> \endverbatim
*>
*> \param[in] JOBU2
*> \verbatim
*>          JOBU2 is CHARACTER*1
*>          = 'Y':  Orthogonal matrix U2 is computed;
*>          = 'N':  U2 is not computed.
*> \endverbatim
*>
*> \param[in] JOBX
*> \verbatim
*>          JOBX is CHARACTER*1
*>          = 'Y':  Matrix X is computed;
*>          = 'N':  X is not computed.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 1.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B.  N >= 1.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 1.
*> \endverbatim
*>
*> \param[out] L
*> \verbatim
*>          L is INTEGER
*>          On exit, the effective numerical rank of the matrix
*>          (A**T, B**T)**T.
*> \endverbatim
*>
*> \param[out] SWAPPED
*> \verbatim
*>          L is LOGICAL
*>          On exit, SWAPPED is true if DGGQRCS swapped the input
*>          matrices A, B and computed the GSVD of (B, A); false
*>          otherwise.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,N)
*>          On entry, the P-by-N matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,P).
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>
*>          On exit, ALPHA and BETA contain the K generalized singular
*>          value pairs of A and B.
*> \endverbatim
*>
*> \param[out] U1
*> \verbatim
*>          U1 is DOUBLE PRECISION array, dimension (LDU1,M)
*>          If JOBU1 = 'Y', U1 contains the M-by-M orthogonal matrix U1.
*>          If JOBU1 = 'N', U1 is not referenced.
*> \endverbatim
*>
*> \param[in] LDU1
*> \verbatim
*>          LDU1 is INTEGER
*>          The leading dimension of the array U1. LDU1 >= max(1,M) if
*>          JOBU1 = 'Y'; LDU1 >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] U2
*> \verbatim
*>          U2 is DOUBLE PRECISION array, dimension (LDU2,P)
*>          If JOBU2 = 'Y', U2 contains the P-by-P orthogonal matrix U2.
*>          If JOBU2 = 'N', U2 is not referenced.
*> \endverbatim
*>
*> \param[in] LDU2
*> \verbatim
*>          LDU2 is INTEGER
*>          The leading dimension of the array U2. LDU2 >= max(1,P) if
*>          JOBU2 = 'Y'; LDU2 >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>
*>          If LWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal size of the WORK array,
*>          returns this value as the first entry of the WORK array, and
*>          no error message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (M + N + P)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  DBBCSD did not converge. For further details, see
*>                subroutine DORCSDBY1.
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  W       DOUBLE PRECISION
*>          W is a radix power chosen such that the Frobenius norm of A
*>          and W*B are with SQRT(RADIX) and 1/SQRT(RADIX) of each
*>          other.
*>
*>  TOL     DOUBLE PRECISION
*>          Let G = (A**T,B**T)**T. TOL is the threshold to determine
*>          the effective rank of G. Generally, it is set to
*>                   TOL = MAX( M + P, N ) * norm(G) * MACHEPS,
*>          where norm(G) is the Frobenius norm of G.
*>          The size of TOL may affect the size of backward error of the
*>          decomposition.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Christoph Conrads (https://christoph-conrads.name)
*
*> \ingroup doubleGEsing
*
*> \par Contributors:
*  ==================
*>
*>     Christoph Conrads (https://christoph-conrads.name)
*>
*
*> \par Further Details:
*  =====================
*>
*>  DGGQRCS should be significantly faster than DGGSVD3 for large
*>  matrices because the matrices A and B are reduced to a pair of
*>  well-conditioned bidiagonal matrices instead of pairs of upper
*>  triangular matrices. On the downside, DGGQRCS requires a much larger
*>  workspace whose dimension must be queried at run-time. DGGQRCS also
*>  offers no guarantees which of the two possible diagonal matrices
*>  is used for the matrix factorization.
*>
*  =====================================================================
      RECURSIVE SUBROUTINE DGGQRCS( JOBU1, JOBU2, JOBX, M, N, P, L,
     $                              SWAPPED,
     $                              A, LDA, B, LDB,
     $                              ALPHA, BETA,
     $                              U1, LDU1, U2, LDU2,
     $                              WORK, LWORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      LOGICAL            SWAPPED
      CHARACTER          JOBU1, JOBU2, JOBX
      INTEGER            INFO, LDA, LDB, LDU1, LDU2, L, M, N, P, LWORK
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ),
     $                   ALPHA( N ), BETA( N ),
     $                   U1( LDU1, * ), U2( LDU2, * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            WANTU1, WANTU2, WANTX, LQUERY
      INTEGER            I, J, K, K1, LMAX, IG, IG11, IG21, IG22,
     $                   IVT, IVT12, LDG, LDX, LDVT, LWKMIN, LWKOPT
      DOUBLE PRECISION   BASE, NORMA, NORMB, NORMG, TOL, ULP, UNFL,
     $                   THETA, IOTA, W
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEQP3, DLACPY, DLAPMT, DLASCL,
     $                   DLASET, DORGQR, DORCSD2BY1, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          COS, MAX, MIN, SIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      WANTU1 = LSAME( JOBU1, 'Y' )
      WANTU2 = LSAME( JOBU2, 'Y' )
      WANTX = LSAME( JOBX, 'Y' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     Test the input arguments
*
      INFO = 0
      IF( .NOT.( WANTU1 .OR. LSAME( JOBU1, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTU2 .OR. LSAME( JOBU2, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTX .OR. LSAME( JOBX, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU1.LT.1 .OR. ( WANTU1 .AND. LDU1.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDU2.LT.1 .OR. ( WANTU2 .AND. LDU2.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -20
      END IF
*
*     Make sure A is the matrix smaller in norm
*
      IF( INFO.EQ.0 ) THEN
         NORMA = DLANGE( 'F', M, N, A, LDA, WORK )
         NORMB = DLANGE( 'F', P, N, B, LDB, WORK )
*
         IF( NORMA.GT.SQRT( 2.0D0 ) * NORMB ) THEN
            CALL DGGQRCS( JOBU2, JOBU1, JOBX, P, N, M, L,
     $                    SWAPPED,
     $                    B, LDB, A, LDA,
     $                    BETA, ALPHA,
     $                    U2, LDU2, U1, LDU1,
     $                    WORK, LWORK, IWORK, INFO )
            SWAPPED = .TRUE.
            RETURN
         ENDIF
*
*     Past this point, we know that
*     * NORMA <= NORMB (almost)
*     * W >= 1
*     * ALPHA will contain cosine values at the end
*     * BETA will contain sine values at the end
*
      END IF
*
*     Initialize variables
*
      SWAPPED = .FALSE.
      L = 0
*     The leading dimension must never be zero
      LDG = MAX( M + P, 1 )
      LDVT = N
      LMAX = MIN( M + P, N )
      IG = 1
      IG11 = 1
      IG21 = M + 1
      IG22 = LDG * M + M + 1
      IVT = LDG * N + 2
      IVT12 = IVT + LDVT * M
      THETA = -1
      IOTA = -1
      W = -1
*
*     Compute workspace
*
      IF( INFO.EQ.0 ) THEN
         LWKMIN = 0
         LWKOPT = 0
*
         CALL DGEQP3( M + P, N, WORK( IG ), LDG, IWORK, ALPHA, WORK, -1,
     $                INFO )
         LWKMIN = MAX( LWKMIN, 3 * N + 1 )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
*
         CALL DORGQR( M + P, LMAX, LMAX, WORK( IG ), LDG, ALPHA, WORK,
     $                -1, INFO )
         LWKMIN = MAX( LWKMIN, LMAX )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
*
         CALL DORCSD2BY1( JOBU1, JOBU2, JOBX, M + P, M, LMAX,
     $                    WORK( IG ), LDG, WORK( IG ), LDG,
     $                    ALPHA,
     $                    U1, LDU1, U2, LDU2, WORK( IVT ), LDVT,
     $                    WORK, -1, IWORK, INFO )
         LWKMIN = MAX( LWKMIN, INT( WORK( 1 ) ) )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
*        The matrix (A, B) must be stored sequentially for DORGQR
         LWKMIN = LWKMIN + IVT
         LWKOPT = LWKOPT + IVT
*        2-by-1 CSD matrix V1 must be stored
         IF( WANTX ) THEN
            LWKMIN = LWKMIN + LDVT*N
            LWKOPT = LWKOPT + LDVT*N
         END IF
*        Check for minimum workspace size
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -20
         END IF
*
         WORK( 1 ) = DBLE( LWKOPT )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGQRCS', -INFO )
         RETURN
      END IF
      IF( LQUERY ) THEN
         RETURN
      ENDIF
*     Finish initialization
      IF( .NOT.WANTX ) THEN
         LDVT = 0
      END IF
*
*     Scale matrix A such that norm(A) \approx norm(B)
*
      IF( NORMA.EQ.0.0D0 ) THEN
         W = 1.0D0
      ELSE
         BASE = DLAMCH( 'B' )
         W = BASE ** INT( LOG( NORMB / NORMA ) / LOG( BASE ) )
*
         CALL DLASCL( 'G', -1, -1, 1.0D0, W, M, N, A, LDA, INFO )
         IF ( INFO.NE.0 ) THEN
            RETURN
         END IF
      END IF
*
*     Copy matrices A, B into the (M+P) x N matrix G
*
      CALL DLACPY( 'A', M, N, A, LDA, WORK( IG11 ), LDG )
      CALL DLACPY( 'A', P, N, B, LDB, WORK( IG21 ), LDG )
*
*     Compute the Frobenius norm of matrix G
*
      NORMG = NORMB * SQRT( 1.0D0 + ( ( W * NORMA ) / NORMB )**2 )
*
*     Get machine precision and set up threshold for determining
*     the effective numerical rank of the matrix G.
*
      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe Minimum' )
      TOL = MAX( M + P, N ) * MAX( NORMG, UNFL ) * ULP
*
*     IWORK stores the column permutations computed by DGEQP3.
*     Columns J where IWORK( J ) is non-zero are permuted to the front
*     so we set the all entries to zero here.
*
      IWORK( 1:N ) = 0
*
*     Compute the QR factorization with column pivoting GΠ = Q1 R1
*
      CALL DGEQP3( M + P, N, WORK( IG ), LDG, IWORK, ALPHA,
     $             WORK( IVT ), LWORK - IVT + 1, INFO )
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Determine the rank of G
*
      DO I = 1, MIN( M + P, N )
         IF( ABS( WORK( (I-1) * LDG + I ) ).LE.TOL ) THEN
            EXIT
         END IF
         L = L + 1
      END DO
*
*     Handle rank=0 case
*
      IF( L.EQ.0 ) THEN
         IF( WANTU1 ) THEN
            CALL DLASET( 'A', M, M, 0.0D0, 1.0D0, U1, LDU1 )
         END IF
         IF( WANTU2 ) THEN
            CALL DLASET( 'A', P, P, 0.0D0, 1.0D0, U2, LDU2 )
         END IF
*
         WORK( 1 ) = DBLE( LWKOPT )
         RETURN
      END IF
*
*     Copy R1( 1:L, : ) into A, B and set lower triangular part to zero
*
      IF( WANTX ) THEN
         IF( L.LE.M ) THEN
             CALL DLACPY( 'U', L, N, WORK( IG ), LDG, A, LDA )
             CALL DLASET( 'L', L - 1, N, 0.0D0, 0.0D0, A( 2, 1 ), LDA )
         ELSE
             CALL DLACPY( 'U', M, N, WORK( IG ), LDG, A, LDA )
             CALL DLACPY( 'U', L - M, N - M, WORK( IG22 ), LDG, B, LDB )
*
             CALL DLASET( 'L', M - 1, N, 0.0D0, 0.0D0, A( 2, 1 ), LDA )
             CALL DLASET( 'L', L-M-1, N, 0.0D0, 0.0D0, B( 2, 1 ), LDB )
         END IF
      END IF
*
*     Explicitly form Q1 so that we can compute the CS decomposition
*
      CALL DORGQR( M + P, L, L, WORK( IG ), LDG, ALPHA,
     $             WORK( IVT ), LWORK - IVT + 1, INFO )
      IF ( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Compute the CS decomposition of Q1( :, 1:L )
*
      K = MIN( M, P, L, M + P - L )
      K1 = MAX( L - P, 0 )
      CALL DORCSD2BY1( JOBU1, JOBU2, JOBX, M + P, M, L,
     $                 WORK( IG11 ), LDG, WORK( IG21 ), LDG,
     $                 ALPHA,
     $                 U1, LDU1, U2, LDU2, WORK( IVT ), LDVT,
     $                 WORK( IVT + LDVT*N ), LWORK - IVT - LDVT*N + 1,
     $                 IWORK( N + 1 ), INFO )
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Compute X = V^T R1( 1:L, : ) and adjust for matrix scaling
*
      IF( WANTX ) THEN
         LDX = L
         IF ( L.LE.M ) THEN
            CALL DGEMM( 'N', 'N', L, N, L,
     $                  1.0D0, WORK( IVT ), LDVT, A, LDA,
     $                  0.0D0, WORK( 2 ), LDX )
         ELSE
            CALL DGEMM( 'N', 'N', L, N, M,
     $                  1.0D0, WORK( IVT ), LDVT, A, LDA,
     $                  0.0D0, WORK( 2 ), LDX )
            CALL DGEMM( 'N', 'N', L, N - M, L - M,
     $                  1.0D0, WORK( IVT12 ), LDVT, B, LDB,
     $                  1.0D0, WORK( L*M + 2 ), LDX )
         END IF
*        Revert column permutation Π by permuting the columns of X
         CALL DLAPMT( .FALSE., L, N, WORK( 2 ), LDX, IWORK )
      END IF
*
*     Adjust generalized singular values for matrix scaling
*     Compute sine, cosine values
*     Prepare row scaling of X
*
      DO I = 1, K
         THETA = ALPHA( I )
*        Do not adjust singular value if THETA is greater
*        than pi/2 (infinite singular values won't change)
         IF( COS( THETA ).LE.0.0D0 ) THEN
            ALPHA( I ) = 0.0D0
            BETA( I ) = 1.0D0
            IF( WANTX ) THEN
               WORK( IVT + I ) = 1.0D0
            END IF
         ELSE
*           iota comes in the greek alphabet after theta
            IOTA = ATAN( W * TAN( THETA ) )
*           ensure sine, cosine divisor is far away from zero
*           w is a power of two and will cause no trouble
            IF( SIN( IOTA ) .GE. COS( IOTA ) ) THEN
               ALPHA( I ) =  ( SIN( IOTA ) / TAN( THETA ) ) / W
               BETA( I ) = SIN( IOTA )
               IF( WANTX ) THEN
                  WORK( IVT + I ) = SIN( THETA ) / SIN( IOTA )
               END IF
            ELSE
               ALPHA( I ) = COS( IOTA )
               BETA( I ) = SIN( IOTA )
               IF( WANTX ) THEN
                  WORK( IVT + I ) = COS( THETA ) / COS( IOTA ) / W
               END IF
            END IF
         END IF
      END DO
*     Adjust rows of X for matrix scaling
      IF( WANTX ) THEN
         DO J = 0, N-1
            DO I = 1, K1
               WORK( LDX*J + I + 1 ) = WORK( LDX*J + I + 1 ) / W
            END DO
            DO I = 1, K
               WORK( LDX*J + I + K1 + 1 ) =
     $         WORK( LDX*J + I + K1 + 1 ) * WORK( IVT + I )
            END DO
         END DO
      END IF
*
      WORK( 1 ) = DBLE( LWKOPT )
      RETURN
*
*     End of DGGQRCS
*
      END
