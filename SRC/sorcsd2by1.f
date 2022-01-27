*> \brief \b SORCSD2BY1
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SORCSD2BY1 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorcsd2by1.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorcsd2by1.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorcsd2by1.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SORCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11,
*                              X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T,
*                              LDV1T, WORK, LWORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBU1, JOBU2, JOBV1T
*       INTEGER            INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21,
*      $                   M, P, Q
*       ..
*       .. Array Arguments ..
*       REAL               THETA(*)
*       REAL               U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*),
*      $                   X11(LDX11,*), X21(LDX21,*)
*       INTEGER            IWORK(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*>\verbatim
*>
*> SORCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with
*> orthonormal columns that has been partitioned into a 2-by-1 block
*> structure:
*>
*>                                [  I1 0  0 ]
*>                                [  0  C  0 ]
*>          [ X11 ]   [ U1 |    ] [  0  0  0 ]
*>      X = [-----] = [---------] [----------] V1**T .
*>          [ X21 ]   [    | U2 ] [  0  0  0 ]
*>                                [  0  S  0 ]
*>                                [  0  0  I2]
*>
*> X11 is P-by-Q. The orthogonal matrices U1, U2, and V1 are P-by-P,
*> (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R
*> nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which
*> R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a
*> K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBU1
*> \verbatim
*>          JOBU1 is CHARACTER
*>          = 'Y':      U1 is computed;
*>          otherwise:  U1 is not computed.
*> \endverbatim
*>
*> \param[in] JOBU2
*> \verbatim
*>          JOBU2 is CHARACTER
*>          = 'Y':      U2 is computed;
*>          otherwise:  U2 is not computed.
*> \endverbatim
*>
*> \param[in] JOBV1T
*> \verbatim
*>          JOBV1T is CHARACTER
*>          = 'Y':      V1T is computed;
*>          otherwise:  V1T is not computed.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows in X.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows in X11. 0 <= P <= M.
*> \endverbatim
*>
*> \param[in] Q
*> \verbatim
*>          Q is INTEGER
*>          The number of columns in X11 and X21. 0 <= Q <= M.
*> \endverbatim
*>
*> \param[in,out] X11
*> \verbatim
*>          X11 is REAL array, dimension (LDX11,Q)
*>          On entry, part of the orthogonal matrix whose CSD is desired.
*> \endverbatim
*>
*> \param[in] LDX11
*> \verbatim
*>          LDX11 is INTEGER
*>          The leading dimension of X11. LDX11 >= MAX(1,P).
*> \endverbatim
*>
*> \param[in,out] X21
*> \verbatim
*>          X21 is REAL array, dimension (LDX21,Q)
*>          On entry, part of the orthogonal matrix whose CSD is desired.
*> \endverbatim
*>
*> \param[in] LDX21
*> \verbatim
*>          LDX21 is INTEGER
*>           The leading dimension of X21. LDX21 >= MAX(1,M-P).
*> \endverbatim
*>
*> \param[out] THETA
*> \verbatim
*>          THETA is REAL array, dimension (R), in which R =
*>          MIN(P,M-P,Q,M-Q).
*>          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and
*>          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ).
*> \endverbatim
*>
*> \param[out] U1
*> \verbatim
*>          U1 is REAL array, dimension (P)
*>          If JOBU1 = 'Y', U1 contains the P-by-P orthogonal matrix U1.
*> \endverbatim
*>
*> \param[in] LDU1
*> \verbatim
*>          LDU1 is INTEGER
*>          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >=
*>          MAX(1,P).
*> \endverbatim
*>
*> \param[out] U2
*> \verbatim
*>          U2 is REAL array, dimension (M-P)
*>          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) orthogonal
*>          matrix U2.
*> \endverbatim
*>
*> \param[in] LDU2
*> \verbatim
*>          LDU2 is INTEGER
*>          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >=
*>          MAX(1,M-P).
*> \endverbatim
*>
*> \param[out] V1T
*> \verbatim
*>          V1T is REAL array, dimension (Q)
*>          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix orthogonal
*>          matrix V1**T.
*> \endverbatim
*>
*> \param[in] LDV1T
*> \verbatim
*>          LDV1T is INTEGER
*>          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=
*>          MAX(1,Q).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*>          If INFO > 0 on exit, WORK(2:R) contains the values PHI(1),
*>          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R),
*>          define the matrix in intermediate bidiagonal-block form
*>          remaining after nonconvergence. INFO specifies the number
*>          of nonzero PHI's.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the work array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (M-MIN(P,M-P,Q,M-Q))
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  SBBCSD did not converge. See the description of WORK
*>                above for details.
*> \endverbatim
*
*> \par References:
*  ================
*>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
*>      Algorithms, 50(1):33-65, 2009.
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realOTHERcomputational
*
*  =====================================================================
      SUBROUTINE SORCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11,
     $                       X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T,
     $                       LDV1T, WORK, LWORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBU1, JOBU2, JOBV1T
      INTEGER            INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21,
     $                   M, P, Q
*     ..
*     .. Array Arguments ..
      REAL               THETA(*)
      REAL               U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*),
     $                   X11(LDX11,*), X21(LDX21,*)
      INTEGER            IWORK(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            CHILDINFO, I, IB11D, IB11E, IB12D, IB12E,
     $                   IB21D, IB21E, IB22D, IB22E, IBBCSD, IORBDB,
     $                   IORGLQ, IORGQR, IPHI, ITAUP1, ITAUP2, ITAUQ1,
     $                   J, LBBCSD, LORBDB, LORGLQ, LORGLQMIN,
     $                   LORGLQOPT, LORGQR, LORGQRMIN, LORGQROPT,
     $                   LWORKMIN, LWORKOPT, R
      LOGICAL            LQUERY, WANTU1, WANTU2, WANTV1T
*     ..
*     .. Local Arrays ..
      REAL               DUM1(1), DUM2(1,1)
*     ..
*     .. External Subroutines ..
      EXTERNAL           SBBCSD, SCOPY, SLACPY, SLAPMR, SLAPMT, SORBDB1,
     $                   SORBDB2, SORBDB3, SORBDB4, SORGLQ, SORGQR,
     $                   XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Function ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test input arguments
*
      INFO = 0
      WANTU1 = LSAME( JOBU1, 'Y' )
      WANTU2 = LSAME( JOBU2, 'Y' )
      WANTV1T = LSAME( JOBV1T, 'Y' )
      LQUERY = LWORK .EQ. -1
*
      IF( M .LT. 0 ) THEN
         INFO = -4
      ELSE IF( P .LT. 0 .OR. P .GT. M ) THEN
         INFO = -5
      ELSE IF( Q .LT. 0 .OR. Q .GT. M ) THEN
         INFO = -6
      ELSE IF( LDX11 .LT. MAX( 1, P ) ) THEN
         INFO = -8
      ELSE IF( LDX21 .LT. MAX( 1, M-P ) ) THEN
         INFO = -10
      ELSE IF( WANTU1 .AND. LDU1 .LT. MAX( 1, P ) ) THEN
         INFO = -13
      ELSE IF( WANTU2 .AND. LDU2 .LT. MAX( 1, M - P ) ) THEN
         INFO = -15
      ELSE IF( WANTV1T .AND. LDV1T .LT. MAX( 1, Q ) ) THEN
         INFO = -17
      END IF
*
      R = MIN( P, M-P, Q, M-Q )
*
*     Compute workspace
*
*       WORK layout:
*     |-------------------------------------------------------|
*     | LWORKOPT (1)                                          |
*     |-------------------------------------------------------|
*     | PHI (MAX(1,R-1))                                      |
*     |-------------------------------------------------------|
*     | TAUP1 (MAX(1,P))                        | B11D (R)    |
*     | TAUP2 (MAX(1,M-P))                      | B11E (R-1)  |
*     | TAUQ1 (MAX(1,Q))                        | B12D (R)    |
*     |-----------------------------------------| B12E (R-1)  |
*     | SORBDB WORK | SORGQR WORK | SORGLQ WORK | B21D (R)    |
*     |             |             |             | B21E (R-1)  |
*     |             |             |             | B22D (R)    |
*     |             |             |             | B22E (R-1)  |
*     |             |             |             | SBBCSD WORK |
*     |-------------------------------------------------------|
*
      IF( INFO .EQ. 0 ) THEN
         IPHI = 2
         IB11D = IPHI + MAX( 1, R-1 )
         IB11E = IB11D + MAX( 1, R )
         IB12D = IB11E + MAX( 1, R - 1 )
         IB12E = IB12D + MAX( 1, R )
         IB21D = IB12E + MAX( 1, R - 1 )
         IB21E = IB21D + MAX( 1, R )
         IB22D = IB21E + MAX( 1, R - 1 )
         IB22E = IB22D + MAX( 1, R )
         IBBCSD = IB22E + MAX( 1, R - 1 )
         ITAUP1 = IPHI + MAX( 1, R-1 )
         ITAUP2 = ITAUP1 + MAX( 1, P )
         ITAUQ1 = ITAUP2 + MAX( 1, M-P )
         IORBDB = ITAUQ1 + MAX( 1, Q )
         IORGQR = ITAUQ1 + MAX( 1, Q )
         IORGLQ = ITAUQ1 + MAX( 1, Q )
         LORGQRMIN = 1
         LORGQROPT = 1
         LORGLQMIN = 1
         LORGLQOPT = 1
         IF( R .EQ. Q ) THEN
            CALL SORBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA,
     $                    DUM1, DUM1, DUM1, DUM1, WORK, -1,
     $                    CHILDINFO )
            LORBDB = INT( WORK(1) )
            IF( WANTU1 .AND. P .GT. 0 ) THEN
               CALL SORGQR( P, P, Q, U1, LDU1, DUM1, WORK(1), -1,
     $                      CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            ENDIF
            IF( WANTU2 .AND. M-P .GT. 0 ) THEN
               CALL SORGQR( M-P, M-P, Q, U2, LDU2, DUM1, WORK(1), -1,
     $                      CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, M-P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTV1T .AND. Q .GT. 0 ) THEN
               CALL SORGLQ( Q-1, Q-1, Q-1, V1T, LDV1T,
     $                      DUM1, WORK(1), -1, CHILDINFO )
               LORGLQMIN = MAX( LORGLQMIN, Q-1 )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            END IF
            CALL SBBCSD( JOBU1, JOBU2, JOBV1T, 'N', 'N', M, P, Q, THETA,
     $                   DUM1, U1, LDU1, U2, LDU2, V1T, LDV1T, DUM2,
     $                   1, DUM1, DUM1, DUM1, DUM1, DUM1,
     $                   DUM1, DUM1, DUM1, WORK(1), -1, CHILDINFO
     $                 )
            LBBCSD = INT( WORK(1) )
         ELSE IF( R .EQ. P ) THEN
            CALL SORBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA,
     $                    DUM1, DUM1, DUM1, DUM1, WORK(1), -1,
     $                    CHILDINFO )
            LORBDB = INT( WORK(1) )
            IF( WANTU1 .AND. P .GT. 0 ) THEN
               CALL SORGQR( P-1, P-1, P-1, U1(2,2), LDU1, DUM1,
     $                      WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, P-1 )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTU2 .AND. M-P .GT. 0 ) THEN
               CALL SORGQR( M-P, M-P, Q, U2, LDU2, DUM1, WORK(1), -1,
     $                      CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, M-P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTV1T .AND. Q .GT. 0 ) THEN
               CALL SORGLQ( Q, Q, R, V1T, LDV1T, DUM1, WORK(1), -1,
     $                      CHILDINFO )
               LORGLQMIN = MAX( LORGLQMIN, Q )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            END IF
            CALL SBBCSD( JOBV1T, 'N', JOBU1, JOBU2, 'T', M, Q, P, THETA,
     $                   DUM1, V1T, LDV1T, DUM2, 1, U1, LDU1, U2,
     $                   LDU2, DUM1, DUM1, DUM1, DUM1, DUM1,
     $                   DUM1, DUM1, DUM1, WORK(1), -1, CHILDINFO
     $                 )
            LBBCSD = INT( WORK(1) )
         ELSE IF( R .EQ. M-P ) THEN
            CALL SORBDB3( M, P, Q, X11, LDX11, X21, LDX21, THETA,
     $                    DUM1, DUM1, DUM1, DUM1, WORK(1), -1,
     $                    CHILDINFO )
            LORBDB = INT( WORK(1) )
            IF( WANTU1 .AND. P .GT. 0 ) THEN
               CALL SORGQR( P, P, Q, U1, LDU1, DUM1, WORK(1), -1,
     $                      CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTU2 .AND. M-P .GT. 0 ) THEN
               CALL SORGQR( M-P-1, M-P-1, M-P-1, U2(2,2), LDU2, DUM1,
     $                      WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, M-P-1 )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTV1T .AND. Q .GT. 0 ) THEN
               CALL SORGLQ( Q, Q, R, V1T, LDV1T, DUM1, WORK(1), -1,
     $                      CHILDINFO )
               LORGLQMIN = MAX( LORGLQMIN, Q )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            END IF
            CALL SBBCSD( 'N', JOBV1T, JOBU2, JOBU1, 'T', M, M-Q, M-P,
     $                   THETA, DUM1, DUM2, 1, V1T, LDV1T, U2, LDU2,
     $                   U1, LDU1, DUM1, DUM1, DUM1, DUM1,
     $                   DUM1, DUM1, DUM1, DUM1, WORK(1), -1,
     $                   CHILDINFO )
            LBBCSD = INT( WORK(1) )
         ELSE
            CALL SORBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA,
     $                    DUM1, DUM1, DUM1, DUM1, DUM1,
     $                    WORK(1), -1, CHILDINFO )
            LORBDB = M + INT( WORK(1) )
            IF( WANTU1 .AND. P .GT. 0 ) THEN
               CALL SORGQR( P, P, M-Q, U1, LDU1, DUM1, WORK(1), -1,
     $                      CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTU2 .AND. M-P .GT. 0 ) THEN
               CALL SORGQR( M-P, M-P, M-Q, U2, LDU2, DUM1, WORK(1),
     $                      -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, M-P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTV1T .AND. Q .GT. 0 ) THEN
               CALL SORGLQ( Q, Q, Q, V1T, LDV1T, DUM1, WORK(1), -1,
     $                      CHILDINFO )
               LORGLQMIN = MAX( LORGLQMIN, Q )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            END IF
            CALL SBBCSD( JOBU2, JOBU1, 'N', JOBV1T, 'N', M, M-P, M-Q,
     $                   THETA, DUM1, U2, LDU2, U1, LDU1, DUM2, 1,
     $                   V1T, LDV1T, DUM1, DUM1, DUM1, DUM1,
     $                   DUM1, DUM1, DUM1, DUM1, WORK(1), -1,
     $                   CHILDINFO )
            LBBCSD = INT( WORK(1) )
         END IF
         LWORKMIN = MAX( IORBDB+LORBDB-1,
     $                   IORGQR+LORGQRMIN-1,
     $                   IORGLQ+LORGLQMIN-1,
     $                   IBBCSD+LBBCSD-1 )
         LWORKOPT = MAX( IORBDB+LORBDB-1,
     $                   IORGQR+LORGQROPT-1,
     $                   IORGLQ+LORGLQOPT-1,
     $                   IBBCSD+LBBCSD-1 )
         WORK(1) = LWORKOPT
         IF( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -19
         END IF
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SORCSD2BY1', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
      LORGQR = LWORK-IORGQR+1
      LORGLQ = LWORK-IORGLQ+1
*
*     Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q,
*     in which R = MIN(P,M-P,Q,M-Q)
*
      IF( R .EQ. Q ) THEN
*
*        Case 1: R = Q
*
*        Simultaneously bidiagonalize X11 and X21
*
         CALL SORBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA,
     $                 WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2),
     $                 WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO )
*
*        Accumulate Householder reflectors
*
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            CALL SLACPY( 'L', P, Q, X11, LDX11, U1, LDU1 )
            CALL SORGQR( P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR),
     $                   LORGQR, CHILDINFO )
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            CALL SLACPY( 'L', M-P, Q, X21, LDX21, U2, LDU2 )
            CALL SORGQR( M-P, M-P, Q, U2, LDU2, WORK(ITAUP2),
     $                   WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            V1T(1,1) = ONE
            DO J = 2, Q
               V1T(1,J) = ZERO
               V1T(J,1) = ZERO
            END DO
            IF( Q.GT.1 ) THEN
               CALL SLACPY( 'U', Q-1, Q-1, X21(1,2), LDX21, V1T(2,2),
     $                   LDV1T )
               CALL SORGLQ( Q-1, Q-1, Q-1, V1T(2,2), LDV1T,WORK(ITAUQ1),
     $                   WORK(IORGLQ), LORGLQ, CHILDINFO )
            ENDIF
         END IF
*
*        Simultaneously diagonalize X11 and X21.
*
         CALL SBBCSD( JOBU1, JOBU2, JOBV1T, 'N', 'N', M, P, Q, THETA,
     $                WORK(IPHI), U1, LDU1, U2, LDU2, V1T, LDV1T,
     $                DUM2, 1, WORK(IB11D), WORK(IB11E), WORK(IB12D),
     $                WORK(IB12E), WORK(IB21D), WORK(IB21E),
     $                WORK(IB22D), WORK(IB22E), WORK(IBBCSD), LBBCSD,
     $                CHILDINFO )
*
*        Permute rows and columns to place zero submatrices in
*        preferred positions
*
         IF( Q .GT. 0 .AND. WANTU2 ) THEN
            DO I = 1, Q
               IWORK(I) = M - P - Q + I
            END DO
            DO I = Q + 1, M - P
               IWORK(I) = I - Q
            END DO
            CALL SLAPMT( .FALSE., M-P, M-P, U2, LDU2, IWORK )
         END IF
      ELSE IF( R .EQ. P ) THEN
*
*        Case 2: R = P
*
*        Simultaneously bidiagonalize X11 and X21
*
         CALL SORBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA,
     $                 WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2),
     $                 WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO )
*
*        Accumulate Householder reflectors
*
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            U1(1,1) = ONE
            DO J = 2, P
               U1(1,J) = ZERO
               U1(J,1) = ZERO
            END DO
            CALL SLACPY( 'L', P-1, P-1, X11(2,1), LDX11, U1(2,2), LDU1 )
            CALL SORGQR( P-1, P-1, P-1, U1(2,2), LDU1, WORK(ITAUP1),
     $                   WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            CALL SLACPY( 'L', M-P, Q, X21, LDX21, U2, LDU2 )
            CALL SORGQR( M-P, M-P, Q, U2, LDU2, WORK(ITAUP2),
     $                   WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            CALL SLACPY( 'U', P, Q, X11, LDX11, V1T, LDV1T )
            CALL SORGLQ( Q, Q, R, V1T, LDV1T, WORK(ITAUQ1),
     $                   WORK(IORGLQ), LORGLQ, CHILDINFO )
         END IF
*
*        Simultaneously diagonalize X11 and X21.
*
         CALL SBBCSD( JOBV1T, 'N', JOBU1, JOBU2, 'T', M, Q, P, THETA,
     $                WORK(IPHI), V1T, LDV1T, DUM1, 1, U1, LDU1, U2,
     $                LDU2, WORK(IB11D), WORK(IB11E), WORK(IB12D),
     $                WORK(IB12E), WORK(IB21D), WORK(IB21E),
     $                WORK(IB22D), WORK(IB22E), WORK(IBBCSD), LBBCSD,
     $                CHILDINFO )
*
*        Permute rows and columns to place identity submatrices in
*        preferred positions
*
         IF( Q .GT. 0 .AND. WANTU2 ) THEN
            DO I = 1, Q
               IWORK(I) = M - P - Q + I
            END DO
            DO I = Q + 1, M - P
               IWORK(I) = I - Q
            END DO
            CALL SLAPMT( .FALSE., M-P, M-P, U2, LDU2, IWORK )
         END IF
      ELSE IF( R .EQ. M-P ) THEN
*
*        Case 3: R = M-P
*
*        Simultaneously bidiagonalize X11 and X21
*
         CALL SORBDB3( M, P, Q, X11, LDX11, X21, LDX21, THETA,
     $                 WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2),
     $                 WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO )
*
*        Accumulate Householder reflectors
*
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            CALL SLACPY( 'L', P, Q, X11, LDX11, U1, LDU1 )
            CALL SORGQR( P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR),
     $                   LORGQR, CHILDINFO )
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            U2(1,1) = ONE
            DO J = 2, M-P
               U2(1,J) = ZERO
               U2(J,1) = ZERO
            END DO
            CALL SLACPY( 'L', M-P-1, M-P-1, X21(2,1), LDX21, U2(2,2),
     $                   LDU2 )
            CALL SORGQR( M-P-1, M-P-1, M-P-1, U2(2,2), LDU2,
     $                   WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            CALL SLACPY( 'U', M-P, Q, X21, LDX21, V1T, LDV1T )
            CALL SORGLQ( Q, Q, R, V1T, LDV1T, WORK(ITAUQ1),
     $                   WORK(IORGLQ), LORGLQ, CHILDINFO )
         END IF
*
*        Simultaneously diagonalize X11 and X21.
*
         CALL SBBCSD( 'N', JOBV1T, JOBU2, JOBU1, 'T', M, M-Q, M-P,
     $                THETA, WORK(IPHI), DUM1, 1, V1T, LDV1T, U2,
     $                LDU2, U1, LDU1, WORK(IB11D), WORK(IB11E),
     $                WORK(IB12D), WORK(IB12E), WORK(IB21D),
     $                WORK(IB21E), WORK(IB22D), WORK(IB22E),
     $                WORK(IBBCSD), LBBCSD, CHILDINFO )
*
*        Permute rows and columns to place identity submatrices in
*        preferred positions
*
         IF( Q .GT. R ) THEN
            DO I = 1, R
               IWORK(I) = Q - R + I
            END DO
            DO I = R + 1, Q
               IWORK(I) = I - R
            END DO
            IF( WANTU1 ) THEN
               CALL SLAPMT( .FALSE., P, Q, U1, LDU1, IWORK )
            END IF
            IF( WANTV1T ) THEN
               CALL SLAPMR( .FALSE., Q, Q, V1T, LDV1T, IWORK )
            END IF
         END IF
      ELSE
*
*        Case 4: R = M-Q
*
*        Simultaneously bidiagonalize X11 and X21
*
         CALL SORBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA,
     $                 WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2),
     $                 WORK(ITAUQ1), WORK(IORBDB), WORK(IORBDB+M),
     $                 LORBDB-M, CHILDINFO )
*
*        Accumulate Householder reflectors
*
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            CALL SCOPY( M-P, WORK(IORBDB+P), 1, U2, 1 )
         END IF
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            CALL SCOPY( P, WORK(IORBDB), 1, U1, 1 )
            DO J = 2, P
               U1(1,J) = ZERO
            END DO
            CALL SLACPY( 'L', P-1, M-Q-1, X11(2,1), LDX11, U1(2,2),
     $                   LDU1 )
            CALL SORGQR( P, P, M-Q, U1, LDU1, WORK(ITAUP1),
     $                   WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            DO J = 2, M-P
               U2(1,J) = ZERO
            END DO
            CALL SLACPY( 'L', M-P-1, M-Q-1, X21(2,1), LDX21, U2(2,2),
     $                   LDU2 )
            CALL SORGQR( M-P, M-P, M-Q, U2, LDU2, WORK(ITAUP2),
     $                   WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            CALL SLACPY( 'U', M-Q, Q, X21, LDX21, V1T, LDV1T )
            CALL SLACPY( 'U', P-(M-Q), Q-(M-Q), X11(M-Q+1,M-Q+1), LDX11,
     $                   V1T(M-Q+1,M-Q+1), LDV1T )
            CALL SLACPY( 'U', -P+Q, Q-P, X21(M-Q+1,P+1), LDX21,
     $                   V1T(P+1,P+1), LDV1T )
            CALL SORGLQ( Q, Q, Q, V1T, LDV1T, WORK(ITAUQ1),
     $                   WORK(IORGLQ), LORGLQ, CHILDINFO )
         END IF
*
*        Simultaneously diagonalize X11 and X21.
*
         CALL SBBCSD( JOBU2, JOBU1, 'N', JOBV1T, 'N', M, M-P, M-Q,
     $                THETA, WORK(IPHI), U2, LDU2, U1, LDU1, DUM1, 1,
     $                V1T, LDV1T, WORK(IB11D), WORK(IB11E), WORK(IB12D),
     $                WORK(IB12E), WORK(IB21D), WORK(IB21E),
     $                WORK(IB22D), WORK(IB22E), WORK(IBBCSD), LBBCSD,
     $                CHILDINFO )
*
*        Permute rows and columns to place identity submatrices in
*        preferred positions
*
         IF( P .GT. R ) THEN
            DO I = 1, R
               IWORK(I) = P - R + I
            END DO
            DO I = R + 1, P
               IWORK(I) = I - R
            END DO
            IF( WANTU1 ) THEN
               CALL SLAPMT( .FALSE., P, P, U1, LDU1, IWORK )
            END IF
            IF( WANTV1T ) THEN
               CALL SLAPMR( .FALSE., P, Q, V1T, LDV1T, IWORK )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of SORCSD2BY1
*
      END

