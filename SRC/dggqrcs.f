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
*       SUBROUTINE DGGQRCS( JOBU1, JOBU2, JOBQT, M, N, P, W, L,
*                           A, LDA, B, LDB,
*                           THETA, U1, LDU1, U2, LDU2, QT, LDQT,
*                           WORK, LWORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBU1, JOB2, JOBQT
*       INTEGER            INFO, LDA, LDB, LDU1, LDU2, LDQT,
*      $                   M, N, P, L, LWORK
*       DOUBLE PRECISION   W
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), THETA( * ),
*      $                   U1( LDU1, * ), U2( LDU2, * ), QT( LDQ, * ),
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
*>       U1**T*A*Q = D1*( 0 R ),    U2**T*B*Q = D2*( 0 R )
*>
*> where U1, U2, and Q are orthogonal matrices. DGGQRCS uses the QR
*> factorization with column pivoting and the 2-by-1 CS decomposition to
*> compute the GSVD.
*>
*> Let L be the effective numerical rank of the matrix (A**T,B**T)**T,
*> then R is a L-by-L nonsingular upper triangular matrix, D1 and
*> D2 are M-by-L and P-by-L "diagonal" matrices and of the
*> following structures, respectively:
*>
*>                        K   K1
*>        D1 =     (  0   0   0 )
*>              K  (  0   S   0 )
*>              K1 (  0   0   I )
*>
*>                    K2  K
*>        D2 =  K2 (  I   0   0 )
*>              K  (  0   C   0 )
*>                 (  0   0   0 )
*>
*>                 N-L  L
*>   ( 0 R ) = L (  0   R )
*>
*> where
*>
*>   K  = MIN(M, P, L, M + P - L),
*>   K1 = MAX(L - P, 0),
*>   K2 = MAX(L - M, 0),
*>   C  = diag( COS(THETA(1)), ..., COS(THETA(K)) ),
*>   S  = diag( SIN(THETA(1)), ..., SIN(THETA(K)) ), and
*>   C^2 + S^2 = I.
*>
*> The routine computes C, S, R, and optionally the orthogonal
*> transformation matrices U, V and Q. If L <= M, then R is stored in
*> A(1:L, :) on exit. Otherwise, the first M rows of R are stored in A
*> and the last L-M rows are stored in B(1:L-M, :).
*>
*> In particular, if B is an N-by-N nonsingular matrix, then the GSVD of
*> A and B implicitly gives the SVD of A*inv(B):
*>                      A*inv(B) = U1*(D1*inv(D2))*U2**T.
*> If (A**T,B**T)**T  has orthonormal columns, then the GSVD of A and B
*> is also equal to the CS decomposition of A and B. Furthermore, the
*> GSVD can be used to derive the solution of the eigenvalue problem:
*>                      A**T*A x = lambda * B**T*B x.
*> In some literature, the GSVD of A and B is presented in the form
*>                  U1**T*A*X = ( 0 D1 ),   U2**T*B*X = ( 0 D2 )
*> where U1 and U2 are orthogonal and X is nonsingular, D1 and D2 are
*> ``diagonal''.  The former GSVD form can be converted to the latter
*> form by taking the nonsingular matrix X as
*>
*>                      X = Q*( I   0    )
*>                            ( 0 inv(R) ).
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
*> \param[in] JOBQT
*> \verbatim
*>          JOBQT is CHARACTER*1
*>          = 'Y':  Orthogonal matrix Q is computed;
*>          = 'N':  Q is not computed.
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
*>          The number of columns of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W in DOUBLE PRECISION
*>
*>          On exit, W is a radix power chosen such that the Frobenius
*>          norm of A and W*B are within sqrt(radix) and 1/sqrt(radix)
*>          of each other.
*> \endverbatim
*>
*> \param[out] L
*> \verbatim
*>          L is INTEGER
*>          On exit, the effective numerical rank of the matrix
*>          (A**T, B**T)**T.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, A contains the triangular matrix R or the first M
*>          rows of R, respectively. See Purpose for details.
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
*>          On exit, if L > M, then B contains the last L - M rows of
*>          the triangular matrix R. See Purpose for details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,P).
*> \endverbatim
*>
*> \param[out] THETA
*> \verbatim
*>          THETA is DOUBLE PRECISION array, dimension (N)
*>
*>          On exit, THETA contains K = MIN(M, P, L, M + P - L) values
*>          in radians in ascending order.
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
*> \param[out] QT
*> \verbatim
*>          QT is DOUBLE PRECISION array, dimension (LDQT,N)
*>          If JOBQT = 'Y', QT contains the N-by-N orthogonal matrix
*>          Q**T.
*> \endverbatim
*>
*> \param[in] LDQT
*> \verbatim
*>          LDQT is INTEGER
*>          The leading dimension of the array QT. LDQT >= max(1,N) if
*>          JOBQT = 'Y'; LDQT >= 1 otherwise.
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
*>  TOL     DOUBLE PRECISION
*>          Let G = (A**T,B**T)**T. TOL is the threshold to determine
*>          the effective rank of G. Generally, it is set to
*>                   TOL = MAX(M,P,N) * norm(G) * MACHEPS,
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
*> \date September 2016
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
*>  DGGQRCS should be significantly faster than DGGSVD and DGGSVD3 for
*>  large matrices because the matrices A and B are reduced to a pair of
*>  well-conditioned bidiagonal matrices instead of pairs of upper
*>  triangular matrices. On the downside, DGGQRCS requires a much larger
*>  workspace whose dimension must be queried at run-time.
*>
*  =====================================================================
      SUBROUTINE DGGQRCS( JOBU1, JOBU2, JOBQT, M, N, P, W, L,
     $                    A, LDA, B, LDB,
     $                    THETA, U1, LDU1, U2, LDU2, QT, LDQT,
     $                    WORK, LWORK, IWORK, INFO )
*
*  -- LAPACK driver routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
*     September 2016
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      CHARACTER          JOBU1, JOBU2, JOBQT
      INTEGER            INFO, LDA, LDB, LDU1, LDU2, LDQT,
     $                   L, M, N, P, LWORK
      DOUBLE PRECISION   W
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), THETA( * ),
     $                   U1( LDU1, * ), U2( LDU2, * ), QT( LDQT, * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            WANTU1, WANTU2, WANTQT, LQUERY
      INTEGER            I, J, Z, R, LDG, LWKOPT
      DOUBLE PRECISION   GNORM, TOL, ULP, UNFL, NORMA, NORMB, BASE
*     .. Local Arrays ..
      DOUBLE PRECISION   G( M + P, N )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLASCL, DGEQP3, DORGQR, DGERQF, QORGRQ,
     $                   DORCSD2BY1, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      WANTU1 = LSAME( JOBU1, 'Y' )
      WANTU2 = LSAME( JOBU2, 'Y' )
      WANTQT = LSAME( JOBQT, 'Y' )
      LQUERY = ( LWORK.EQ.-1 )
      LWKOPT = 1
*
*     Initialize variables
*
      L = MIN( M + P, N )
      Z = ( M + P ) * N
      IF ( LQUERY ) THEN
         G = 0
      ELSE
         G = RESHAPE( WORK(1:Z), (/ M + P, N /) )
      END IF
      LDG = M + P
*
*     Test the input arguments
*
      INFO = 0
      IF( .NOT.( WANTU1 .OR. LSAME( JOBU1, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTU2 .OR. LSAME( JOBU2, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTQT .OR. LSAME( JOBQT, 'N' ) ) ) THEN
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
         INFO = -15
      ELSE IF( LDU2.LT.1 .OR. ( WANTU2 .AND. LDU2.LT.P ) ) THEN
         INFO = -17
      ELSE IF( LDQT.LT.1 .OR. ( WANTQT .AND. LDQT.LT.N ) ) THEN
         INFO = -19
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -23
      END IF
*
*     Compute workspace
*
      IF( INFO.EQ.0 ) THEN
         CALL DGEQP3( M+P, N, G, LDG, IWORK, THETA, WORK, -1, INFO )
         LWKOPT = INT( WORK( 1 ) )

         CALL DORGQR( M + P, L, L, G, LDG, THETA, WORK, -1, INFO )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )

         CALL DORCSD2BY1( JOBU2, JOBU1, 'Y', M + P, P, L,
     $                    G, LDG, G, LDG,
     $                    THETA, U2, LDU2, U1, LDU1, QT, LDQT,
     $                    WORK, -1, IWORK, INFO )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         LWKOPT = Z + LWKOPT

*        DGERQF stores L scalar factors for the elementary reflectors
         CALL DGERQF( L, N, QT, LDQT, WORK, WORK, -1, INFO )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + L )

         CALL DORGRQ( N, N, L, QT, LDQT, WORK, WORK, -1, INFO )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) + L )

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
*
*     Scale matrix B such that norm(A) \approx norm(B)
*
      NORMA = DLANGE( 'F', M, N, A, LDA, WORK )
      NORMB = DLANGE( 'F', P, N, B, LDB, WORK )
*
      IF ( NORMB.EQ.0 ) THEN
         W = 1.0D0
      ELSE
         BASE = DLAMCH( 'B' )
         W = BASE ** INT( LOG( NORMA / NORMB ) / LOG( BASE ) )
*
         CALL DLASCL( 'G', -1, -1, 1.0D0, W, P, N, B, LDB, INFO )
         IF ( INFO.NE.0 ) THEN
            RETURN
         END IF
      END IF
*
*     Copy matrices A, B into the (M+P) x n matrix G
*
      CALL DLACPY( 'A', M, N, A, LDA, G( P + 1, 1 ), LDG )
      CALL DLACPY( 'A', P, N, B, LDB, G( 1, 1 ), LDG )
*
*     Compute the Frobenius norm of matrix G
*
      GNORM = DLANGE( 'F', M + P, N, G, LDG, WORK( Z + 1 ) )
*
*     Get machine precision and set up threshold for determining
*     the effective numerical rank of the matrix G.
*
      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe Minimum' )
      TOL = MAX( M + P, N ) * MAX( GNORM, UNFL ) * ULP
*
*     IWORK stores the column permutations computed by DGEQP3.
*     Columns J where IWORK( J ) is non-zero are permuted to the front
*     so we set the all entries to zero here.
*
      DO 10 J = 1, N
         IWORK( J ) = 0
   10 CONTINUE
*
*     Compute the QR factorization with column pivoting GΠ = Q1 R1
*
      CALL DGEQP3( M + P, N, G, LDG, IWORK, THETA,
     $             WORK( Z + 1 ), LWORK - Z, INFO )
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Determine the rank of G
*
      R = 0
      DO 20 I = 1, MIN( M + P, N )
         IF( ABS( G( I, I ) ).LE.TOL ) THEN
            EXIT
         END IF
         R = R + 1
   20 CONTINUE
*
      L = R
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
         IF( WANTQT ) THEN
            CALL DLASET( 'A', N, N, 0.0D0, 1.0D0, QT, LDQT )
         END IF
*
         WORK( 1 ) = DBLE( LWKOPT )
         RETURN
      END IF
*
*     Copy R1 into A
*
      IF( R.LE.M ) THEN
          CALL DLACPY( 'U', R, N, G, LDG, A, LDA )
      ELSE
          CALL DLACPY( 'U', M, N, G, LDG, A, LDA )
          CALL DLACPY( 'A', R - M, N, G( M + 1, 1 ), LDG, B, LDB )
      END IF
*
*     Explicitly form Q1 so that we can compute the CS decomposition
*
      CALL DORGQR( M + P, R, R, G, LDG, THETA,
     $             WORK( Z + 1 ), LWORK - Z, INFO )
      IF ( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Compute the CS decomposition of Q1( :, 1:R )
*
      CALL DORCSD2BY1( JOBU2, JOBU1, 'Y', M + P, P, R,
     $                 G( 1, 1 ), LDG, G( P + 1, 1 ), LDG, THETA,
     $                 U2, LDU2, U1, LDU1, QT, LDQT,
     $                 WORK( Z + 1 ), LWORK - Z, IWORK( N + 1 ), INFO )
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Copy V^T from QT to G
*
      CALL DLACPY( 'A', R, R, QT, LDQT, G, R )
*
*     Compute V^T R1( 1:R, : ) in the last R rows of QT
*
      IF ( R.LE.M ) THEN
         CALL DGEMM( 'N', 'N', R, N, R, 1.0D0, G, R,
     $               A, LDA, 0.0D0, QT( N-R+1, 1 ), LDQT )
      ELSE
         CALL DGEMM( 'N', 'N', R, N, M, 1.0D0, G( 1, 1 ), R,
     $               A, LDA, 0.0D0, QT( N-R+1, 1 ), LDQT )
         CALL DGEMM( 'N', 'N', R, N, R - M, 1.0D0, G( 1, M + 1 ), R,
     $               B, LDB, 1.0D0, QT( N-R+1, 1 ), LDQT )
      END IF
*
*     Compute the RQ decomposition of V^T R1( 1:R, : )
*
      CALL DGERQF( R, N, QT( N-R+1, 1 ), LDQT, WORK,
     $             WORK( L + 1 ), LWORK - L, INFO )
      IF ( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Copy matrix R from QT( 1:R, N-R+1: ) to A, B
*
      IF ( R.LE.M ) THEN
         CALL DLACPY( 'U', R, R, QT( N-R+1, N-R+1 ), LDQT, A, LDA )
      ELSE
         CALL DLACPY( 'U', M,     R, QT( N-R+1, N-R+1 ), LDQT, A, LDA )
         CALL DLACPY( 'U', R - M, R, QT( N-R+M+1, N-R+1 ), LDQT, B, LDB)
      END IF
*
*     Explicitly form Q^T
*
      IF( WANTQT ) THEN
         CALL DORGRQ( N, N, R, QT, LDQT, WORK,
     $                WORK( L + 1 ), LWORK - L, INFO )
         IF ( INFO.NE.0 ) THEN
            RETURN
         END IF
*
*     Revert column permutation Π by permuting the rows of Q^T
*
         CALL DLAPMT( .FALSE., N, N, QT, LDQT, IWORK )
      END IF
*
      WORK( 1 ) = DBLE( LWKOPT )
      RETURN
*
*     End of DGGQRCS
*
      END
