*> \brief \b CUNBDB2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CUNBDB2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunbdb2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunbdb2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunbdb2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CUNBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI,
*                           TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21
*       ..
*       .. Array Arguments ..
*       REAL               PHI(*), THETA(*)
*       COMPLEX            TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*),
*      $                   X11(LDX11,*), X21(LDX21,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*>\verbatim
*>
*> CUNBDB2 simultaneously bidiagonalizes the blocks of a tall and skinny
*> matrix X with orthonormal columns:
*>
*>                            [ B11 ]
*>      [ X11 ]   [ P1 |    ] [  0  ]
*>      [-----] = [---------] [-----] Q1**T .
*>      [ X21 ]   [    | P2 ] [ B21 ]
*>                            [  0  ]
*>
*> X11 is P-by-Q, and X21 is (M-P)-by-Q. P must be no larger than M-P,
*> Q, or M-Q. Routines CUNBDB1, CUNBDB3, and CUNBDB4 handle cases in
*> which P is not the minimum dimension.
*>
*> The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P),
*> and (M-Q)-by-(M-Q), respectively. They are represented implicitly by
*> Householder vectors.
*>
*> B11 and B12 are P-by-P bidiagonal matrices represented implicitly by
*> angles THETA, PHI.
*>
*>\endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           The number of rows X11 plus the number of rows in X21.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>           The number of rows in X11. 0 <= P <= min(M-P,Q,M-Q).
*> \endverbatim
*>
*> \param[in] Q
*> \verbatim
*>          Q is INTEGER
*>           The number of columns in X11 and X21. 0 <= Q <= M.
*> \endverbatim
*>
*> \param[in,out] X11
*> \verbatim
*>          X11 is COMPLEX array, dimension (LDX11,Q)
*>           On entry, the top block of the matrix X to be reduced. On
*>           exit, the columns of tril(X11) specify reflectors for P1 and
*>           the rows of triu(X11,1) specify reflectors for Q1.
*> \endverbatim
*>
*> \param[in] LDX11
*> \verbatim
*>          LDX11 is INTEGER
*>           The leading dimension of X11. LDX11 >= P.
*> \endverbatim
*>
*> \param[in,out] X21
*> \verbatim
*>          X21 is COMPLEX array, dimension (LDX21,Q)
*>           On entry, the bottom block of the matrix X to be reduced. On
*>           exit, the columns of tril(X21) specify reflectors for P2.
*> \endverbatim
*>
*> \param[in] LDX21
*> \verbatim
*>          LDX21 is INTEGER
*>           The leading dimension of X21. LDX21 >= M-P.
*> \endverbatim
*>
*> \param[out] THETA
*> \verbatim
*>          THETA is REAL array, dimension (Q)
*>           The entries of the bidiagonal blocks B11, B21 are defined by
*>           THETA and PHI. See Further Details.
*> \endverbatim
*>
*> \param[out] PHI
*> \verbatim
*>          PHI is REAL array, dimension (Q-1)
*>           The entries of the bidiagonal blocks B11, B21 are defined by
*>           THETA and PHI. See Further Details.
*> \endverbatim
*>
*> \param[out] TAUP1
*> \verbatim
*>          TAUP1 is COMPLEX array, dimension (P-1)
*>           The scalar factors of the elementary reflectors that define
*>           P1.
*> \endverbatim
*>
*> \param[out] TAUP2
*> \verbatim
*>          TAUP2 is COMPLEX array, dimension (Q)
*>           The scalar factors of the elementary reflectors that define
*>           P2.
*> \endverbatim
*>
*> \param[out] TAUQ1
*> \verbatim
*>          TAUQ1 is COMPLEX array, dimension (Q)
*>           The scalar factors of the elementary reflectors that define
*>           Q1.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>           The dimension of the array WORK. LWORK >= M-Q.
*>
*>           If LWORK = -1, then a workspace query is assumed; the routine
*>           only calculates the optimal size of the WORK array, returns
*>           this value as the first entry of the WORK array, and no error
*>           message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           = 0:  successful exit.
*>           < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*>
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup unbdb2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry
*>  in each bidiagonal band is a product of a sine or cosine of a THETA
*>  with a sine or cosine of a PHI. See [1] or CUNCSD for details.
*>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.
*>  See CUNCSD2BY1 for details on generating P1, P2, and Q1 using CUNGQR
*>  and CUNGLQ.
*> \endverbatim
*
*> \par References:
*  ================
*>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
*>      Algorithms, 50(1):33-65, 2009.
*>
*  =====================================================================
      SUBROUTINE CUNBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA,
     $                    PHI,
     $                    TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21
*     ..
*     .. Array Arguments ..
      REAL               PHI(*), THETA(*)
      COMPLEX            TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*),
     $                   X11(LDX11,*), X21(LDX21,*)
*     ..
*
*  ====================================================================
*
*     .. Parameters ..
      COMPLEX            NEGONE
      PARAMETER          ( NEGONE = (-1.0E0,0.0E0) )
*     ..
*     .. Local Scalars ..
      REAL               C, S
      INTEGER            CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5,
     $                   LWORKMIN, LWORKOPT
      LOGICAL            LQUERY
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLARF1F, CLARFGP, CUNBDB5, CSROT, CSCAL,
     $                   CLACGV,
     $                   XERBLA
*     ..
*     .. External Functions ..
      REAL               SCNRM2, SROUNDUP_LWORK
      EXTERNAL           SCNRM2, SROUNDUP_LWORK
*     ..
*     .. Intrinsic Function ..
      INTRINSIC          ATAN2, COS, MAX, SIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test input arguments
*
      INFO = 0
      LQUERY = LWORK .EQ. -1
*
      IF( M .LT. 0 ) THEN
         INFO = -1
      ELSE IF( P .LT. 0 .OR. P .GT. M-P ) THEN
         INFO = -2
      ELSE IF( Q .LT. 0 .OR. Q .LT. P .OR. M-Q .LT. P ) THEN
         INFO = -3
      ELSE IF( LDX11 .LT. MAX( 1, P ) ) THEN
         INFO = -5
      ELSE IF( LDX21 .LT. MAX( 1, M-P ) ) THEN
         INFO = -7
      END IF
*
*     Compute workspace
*
      IF( INFO .EQ. 0 ) THEN
         ILARF = 2
         LLARF = MAX( P-1, M-P, Q-1 )
         IORBDB5 = 2
         LORBDB5 = Q-1
         LWORKOPT = MAX( ILARF+LLARF-1, IORBDB5+LORBDB5-1 )
         LWORKMIN = LWORKOPT
         WORK(1) = SROUNDUP_LWORK(LWORKOPT)
         IF( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) THEN
           INFO = -14
         END IF
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'CUNBDB2', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Reduce rows 1, ..., P of X11 and X21
*
      DO I = 1, P
*
         IF( I .GT. 1 ) THEN
            CALL CSROT( Q-I+1, X11(I,I), LDX11, X21(I-1,I), LDX21, C,
     $                  S )
         END IF
         CALL CLACGV( Q-I+1, X11(I,I), LDX11 )
         CALL CLARFGP( Q-I+1, X11(I,I), X11(I,I+1), LDX11, TAUQ1(I) )
         C = REAL( X11(I,I) )
         CALL CLARF1F( 'R', P-I, Q-I+1, X11(I,I), LDX11, TAUQ1(I),
     $                 X11(I+1,I), LDX11, WORK(ILARF) )
         CALL CLARF1F( 'R', M-P-I+1, Q-I+1, X11(I,I), LDX11,
     $                 TAUQ1(I), X21(I,I), LDX21, WORK(ILARF) )
         CALL CLACGV( Q-I+1, X11(I,I), LDX11 )
         S = SQRT( SCNRM2( P-I, X11(I+1,I), 1 )**2
     $           + SCNRM2( M-P-I+1, X21(I,I), 1 )**2 )
         THETA(I) = ATAN2( S, C )
*
         CALL CUNBDB5( P-I, M-P-I+1, Q-I, X11(I+1,I), 1, X21(I,I), 1,
     $                 X11(I+1,I+1), LDX11, X21(I,I+1), LDX21,
     $                 WORK(IORBDB5), LORBDB5, CHILDINFO )
         CALL CSCAL( P-I, NEGONE, X11(I+1,I), 1 )
         CALL CLARFGP( M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) )
         IF( I .LT. P ) THEN
            CALL CLARFGP( P-I, X11(I+1,I), X11(I+2,I), 1, TAUP1(I) )
            PHI(I) = ATAN2( REAL( X11(I+1,I) ), REAL( X21(I,I) ) )
            C = COS( PHI(I) )
            S = SIN( PHI(I) )
            CALL CLARF1F( 'L', P-I, Q-I, X11(I+1,I), 1,
     $                    CONJG(TAUP1(I)), X11(I+1,I+1), LDX11,
     $                    WORK(ILARF) )
         END IF
         CALL CLARF1F( 'L', M-P-I+1, Q-I, X21(I,I), 1,
     $                 CONJG(TAUP2(I)), X21(I,I+1), LDX21,
     $                 WORK(ILARF) )
*
      END DO
*
*     Reduce the bottom-right portion of X21 to the identity matrix
*
      DO I = P + 1, Q
         CALL CLARFGP( M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) )
         CALL CLARF1F( 'L', M-P-I+1, Q-I, X21(I,I), 1,
     $                 CONJG(TAUP2(I)), X21(I,I+1), LDX21,
     $                 WORK(ILARF) )
      END DO
*
      RETURN
*
*     End of CUNBDB2
*
      END

