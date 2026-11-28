*> \brief \b CERRGG
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CERRGG( PATH, NUNIT )
*
*       .. Scalar Arguments ..
*       CHARACTER*3        PATH
*       INTEGER            NUNIT
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CERRGG tests the error exits for CGGES, CGGESX, CGGEV, CGGEVX,
*> CGGES3, CGGEV3, CGGGLM, CGGHRD, CGGLSE, CGGQRF, CGGRQF,
*> CGGSVD3, CGGSVP3, CHGEQZ, CTGEVC, CTGEXC, CTGSEN, CTGSJA,
*> CTGSNA, CTGSYL, CUNBDB, and CUNCSD.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] PATH
*> \verbatim
*>          PATH is CHARACTER*3
*>          The LAPACK path name for the routines to be tested.
*> \endverbatim
*>
*> \param[in] NUNIT
*> \verbatim
*>          NUNIT is INTEGER
*>          The unit number for output.
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
*> \ingroup complex_eig
*
*  =====================================================================
      SUBROUTINE CERRGG( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX, LW, LWQ, NQ
      PARAMETER          ( NMAX = 3, LW = 6*NMAX, LWQ = 100, NQ = 4 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            DUMMYK, DUMMYL, I, IFST, IHI, ILO, ILST, INFO,
     $                   J, K, M, NCYCLE, NT, SDIM, LWORK
      REAL               ANRM, BNRM, DIF, SCALE, TOLA, TOLB
      REAL               BIG, RNAN, RONE, TINY, TOL
      COMPLEX            EALPHA, EBETA
*     ..
*     .. Local Arrays ..
      LOGICAL            BW( NMAX ), SEL( NMAX )
      INTEGER            IW( LW ), IDUM(NMAX)
      REAL               LS( NMAX ), R1( NMAX ), R2( NMAX ),
     $                   RCE( NMAX ), RCV( NMAX ), RS( NMAX ), RW( LW )
      COMPLEX            A( NMAX, NMAX ), ALPHA( NMAX ),
     $                   B( NMAX, NMAX ), BETA( NMAX ), Q( NMAX, NMAX ),
     $                   TAU( NMAX ), U( NMAX, NMAX ), V( NMAX, NMAX ),
     $                   W( LW ), Z( NMAX, NMAX )
      LOGICAL            BWQ( NQ )
      REAL               RWQ( LWQ )
      COMPLEX            AQ( NQ, NQ ), ALQ( NQ ), BEQ( NQ ),
     $                   BQ( NQ, NQ ), WQ( LWQ )
*     ..
*     .. External Functions ..
      LOGICAL            CLCTES, CLCTSX, LSAMEN
      EXTERNAL           CLCTES, CLCTSX, LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGGES, CGGESX, CGGEV, CGGEVX, CGGGLM, CGGHRD,
     $                   CGGLSE, CGGQRF, CGGRQF, CHGEQZ,
     $                   CHKXER, CTGEVC, CTGEXC, CTGSEN, CTGSJA, CTGSNA,
     $                   CTGSYL, CUNBDB, CUNCSD, CGGES3, CGGEV3, CGGHD3,
     $                   CGGSVD3, CGGSVP3, XLAENV
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
*
*     Set the variables to innocuous values.
*
      DO 20 J = 1, NMAX
         SEL( J ) = .TRUE.
         DO 10 I = 1, NMAX
            A( I, J ) = ZERO
            B( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, NMAX
         A( I, I ) = ONE
         B( I, I ) = ONE
   30 CONTINUE
      OK = .TRUE.
      TOLA = ONE
      TOLB = ONE
      IFST = 1
      ILST = 1
      NT = 0
      LWORK = 1
*
*     Call XLAENV to set the parameters used in CLAQZ0
*
      CALL XLAENV( 12, 10 )
      CALL XLAENV( 13, 12 )
      CALL XLAENV( 14, 13 )
      CALL XLAENV( 15, 2 )
      CALL XLAENV( 17, 10 )
*
*     Test error exits for the GG path.
*
      IF( LSAMEN( 2, C2, 'GG' ) ) THEN
*
*        CGGHRD
*
         SRNAMT = 'CGGHRD'
         INFOT = 1
         CALL CGGHRD( '/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'CGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGHRD( 'N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'CGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGHRD( 'N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'CGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGGHRD( 'N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'CGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGHRD( 'N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'CGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGGHRD( 'N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'CGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CGGHRD( 'N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'CGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGGHRD( 'V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'CGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CGGHRD( 'N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'CGGHRD', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        CGGHD3
*
         SRNAMT = 'CGGHD3'
         INFOT = 1
         CALL CGGHD3( '/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGHD3( 'N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGHD3( 'N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGGHD3( 'N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGHD3( 'N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGGHD3( 'N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CGGHD3( 'N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGGHD3( 'V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CGGHD3( 'N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGHD3', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        CHGEQZ
*
         SRNAMT = 'CHGEQZ'
         INFOT = 1
         CALL CHGEQZ( '/', 'N', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CHGEQZ( 'E', '/', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CHGEQZ( 'E', 'N', '/', 0, 1, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CHGEQZ( 'E', 'N', 'N', -1, 0, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CHGEQZ( 'E', 'N', 'N', 0, 0, 0, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CHGEQZ( 'E', 'N', 'N', 0, 1, 1, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CHGEQZ( 'E', 'N', 'N', 2, 1, 1, A, 1, B, 2, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CHGEQZ( 'E', 'N', 'N', 2, 1, 1, A, 2, B, 1, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL CHGEQZ( 'E', 'V', 'N', 2, 1, 1, A, 2, B, 2, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL CHGEQZ( 'E', 'N', 'V', 2, 1, 1, A, 2, B, 2, ALPHA, BETA,
     $                Q, 1, Z, 1, W, 1, RW, INFO )
         CALL CHKXER( 'CHGEQZ', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        CTGEVC
*
         SRNAMT = 'CTGEVC'
         INFOT = 1
         CALL CTGEVC( '/', 'A', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'CTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CTGEVC( 'R', '/', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'CTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CTGEVC( 'R', 'A', SEL, -1, A, 1, B, 1, Q, 1, Z, 1, 0, M,
     $                W, RW, INFO )
         CALL CHKXER( 'CTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CTGEVC( 'R', 'A', SEL, 2, A, 1, B, 2, Q, 1, Z, 2, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'CTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CTGEVC( 'R', 'A', SEL, 2, A, 2, B, 1, Q, 1, Z, 2, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'CTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CTGEVC( 'L', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'CTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CTGEVC( 'R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W,
     $                RW, INFO )
         CALL CHKXER( 'CTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CTGEVC( 'R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 2, 1, M, W,
     $                RW, INFO )
         CALL CHKXER( 'CTGEVC', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*     Test error exits for the GSV path.
*
      ELSE IF( LSAMEN( 3, PATH, 'GSV' ) ) THEN
*
*        CGGSVD3
*
         SRNAMT = 'CGGSVD3'
         INFOT = 1
         CALL CGGSVD3( '/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGSVD3( 'N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGSVD3( 'N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGGSVD3( 'N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGSVD3( 'N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B,
     $                 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGGSVD3( 'N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B,
     $                 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGGSVD3( 'N', 'N', 'N', 2, 1, 1, DUMMYK, DUMMYL, A, 1, B,
     $                 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CGGSVD3( 'N', 'N', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B,
     $                 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL CGGSVD3( 'U', 'N', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B,
     $                 2, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL CGGSVD3( 'N', 'V', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B,
     $                 2, R1, R2, U, 2, V, 1, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL CGGSVD3( 'N', 'N', 'Q', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B,
     $                 2, R1, R2, U, 2, V, 2, Q, 1, W, LWORK, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL CGGSVD3( 'N', 'N', 'N', 1, 1, 1, DUMMYK, DUMMYL, A, 1, B,
     $                 1, R1, R2, U, 1, V, 1, Q, 1, W, 1, RW, IDUM,
     $                 INFO )
         CALL CHKXER( 'CGGSVD3', INFOT, NOUT, LERR, OK )
         NT = NT + 12
*
*        CGGSVP3
*
         SRNAMT = 'CGGSVP3'
         INFOT = 1
         CALL CGGSVP3( '/', 'N', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGSVP3( 'N', '/', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGSVP3( 'N', 'N', '/', 0, 0, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGGSVP3( 'N', 'N', 'N', -1, 0, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGSVP3( 'N', 'N', 'N', 0, -1, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGGSVP3( 'N', 'N', 'N', 0, 0, -1, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGGSVP3( 'N', 'N', 'N', 2, 1, 1, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGGSVP3( 'N', 'N', 'N', 1, 2, 1, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL CGGSVP3( 'U', 'N', 'N', 2, 2, 2, A, 2, B, 2, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL CGGSVP3( 'N', 'V', 'N', 2, 2, 2, A, 2, B, 2, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 2, V, 1, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL CGGSVP3( 'N', 'N', 'Q', 2, 2, 2, A, 2, B, 2, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 2, V, 2, Q, 1, IW, RW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 25
         CALL CGGSVP3( 'N', 'N', 'N', 1, 1, 1, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W,
     $                 1, INFO )
         CALL CHKXER( 'CGGSVP3', INFOT, NOUT, LERR, OK )
         NT = NT + 12
*
*        CTGSJA
*
         SRNAMT = 'CTGSJA'
         INFOT = 1
         CALL CTGSJA( '/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CTGSJA( 'N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CTGSJA( 'N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CTGSJA( 'N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CTGSJA( 'N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CTGSJA( 'N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CTGSJA( 'N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 0, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CTGSJA( 'N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                0, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL CTGSJA( 'U', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 0, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL CTGSJA( 'N', 'V', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 0, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL CTGSJA( 'N', 'N', 'Q', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 0, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'CTGSJA', INFOT, NOUT, LERR, OK )
         NT = NT + 11
*
*     Test error exits for the GLM path.
*
      ELSE IF( LSAMEN( 3, PATH, 'GLM' ) ) THEN
*
*        CGGGLM
*
         SRNAMT = 'CGGGLM'
         INFOT = 1
         CALL CGGGLM( -1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGGLM( 0, -1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGGLM( 0, 1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGGLM( 0, 0, -1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGGLM( 1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGGLM( 0, 0, 0, A, 0, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGGGLM( 0, 0, 0, A, 1, B, 0, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CGGGLM( 1, 1, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, 1,
     $                INFO )
         CALL CHKXER( 'CGGGLM', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*     Test error exits for the LSE path.
*
      ELSE IF( LSAMEN( 3, PATH, 'LSE' ) ) THEN
*
*        CGGLSE
*
         SRNAMT = 'CGGLSE'
         INFOT = 1
         CALL CGGLSE( -1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGLSE( 0, -1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGLSE( 0, 0, -1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGLSE( 0, 0, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGLSE( 0, 1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGLSE( 0, 0, 0, A, 0, B, 1, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGGLSE( 0, 0, 0, A, 1, B, 0, TAU, ALPHA, BETA, W, LW,
     $                INFO )
         CALL CHKXER( 'CGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CGGLSE( 1, 1, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, 1,
     $                INFO )
         CALL CHKXER( 'CGGLSE', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*     Test error exits for the CSD path.
*
      ELSE IF( LSAMEN( 3, PATH, 'CSD' ) ) THEN
*
*        CUNCSD
*
         SRNAMT = 'CUNCSD'
         INFOT = 7
         CALL CUNCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 -1, 0, 0, A, 1, A,
     $                 1, A, 1, A, 1, RS,
     $                 A, 1, A, 1, A, 1, A,
     $                 1, W, LW, RW, LW, IW, INFO )
         CALL CHKXER( 'CUNCSD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CUNCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, -1, 0, A, 1, A,
     $                 1, A, 1, A, 1, RS,
     $                 A, 1, A, 1, A, 1, A,
     $                 1, W, LW, RW, LW, IW, INFO )
         CALL CHKXER( 'CUNCSD', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CUNCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, -1, A, 1, A,
     $                 1, A, 1, A, 1, RS,
     $                 A, 1, A, 1, A, 1, A,
     $                 1, W, LW, RW, LW, IW, INFO )
         CALL CHKXER( 'CUNCSD', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CUNCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, -1, A,
     $                 1, A, 1, A, 1, RS,
     $                 A, 1, A, 1, A, 1, A,
     $                 1, W, LW, RW, LW, IW, INFO )
         CALL CHKXER( 'CUNCSD', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL CUNCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, 1, A,
     $                 1, A, 1, A, 1, RS,
     $                 A, -1, A, 1, A, 1, A,
     $                 1, W, LW, RW, LW, IW, INFO )
         CALL CHKXER( 'CUNCSD', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL CUNCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, 1, A,
     $                 1, A, 1, A, 1, RS,
     $                 A, 1, A, -1, A, 1, A,
     $                 1, W, LW, RW, LW, IW, INFO )
         CALL CHKXER( 'CUNCSD', INFOT, NOUT, LERR, OK )
         INFOT = 24
         CALL CUNCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, 1, A,
     $                 1, A, 1, A, 1, RS,
     $                 A, 1, A, 1, A, -1, A,
     $                 1, W, LW, RW, LW, IW, INFO )
         CALL CHKXER( 'CUNCSD', INFOT, NOUT, LERR, OK )
         INFOT = 26
         CALL CUNCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, 1, A,
     $                 1, A, 1, A, 1, RS,
     $                 A, 1, A, 1, A, 1, A,
     $                 -1, W, LW, RW, LW, IW, INFO )
         CALL CHKXER( 'CUNCSD', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        CUNBDB
*
         SRNAMT = 'CUNBDB'
         INFOT = 3
         CALL CUNBDB( 'N', 'N', -1, 0, 0, A, 1, A, 1, A, 1, A, 1,
     $                R1, R1, A, A, A, A, W, LW, INFO )
         CALL CHKXER( 'CUNBDB', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CUNBDB( 'N', 'N', 1, -1, 0, A, 1, A, 1, A, 1, A, 1,
     $                R1, R1, A, A, A, A, W, LW, INFO )
         CALL CHKXER( 'CUNBDB', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CUNBDB( 'N', 'N', 1, 1, -1, A, 1, A, 1, A, 1, A, 1,
     $                R1, R1, A, A, A, A, W, LW, INFO )
         CALL CHKXER( 'CUNBDB', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*     Test error exits for the GQR path.
*
      ELSE IF( LSAMEN( 3, PATH, 'GQR' ) ) THEN
*
*        CGGQRF
*
         SRNAMT = 'CGGQRF'
         INFOT = 1
         CALL CGGQRF( -1, 0, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGQRF( 0, -1, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGQRF( 0, 0, -1, A, 1, ALPHA, B, 1, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGQRF( 0, 0, 0, A, 0, ALPHA, B, 1, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGGQRF( 0, 0, 0, A, 1, ALPHA, B, 0, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGGQRF( 1, 1, 2, A, 1, ALPHA, B, 1, BETA, W, 1, INFO )
         CALL CHKXER( 'CGGQRF', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        CGGRQF
*
         SRNAMT = 'CGGRQF'
         INFOT = 1
         CALL CGGRQF( -1, 0, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGRQF( 0, -1, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGRQF( 0, 0, -1, A, 1, ALPHA, B, 1, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGRQF( 0, 0, 0, A, 0, ALPHA, B, 1, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGGRQF( 0, 0, 0, A, 1, ALPHA, B, 0, BETA, W, LW, INFO )
         CALL CHKXER( 'CGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGGRQF( 1, 1, 2, A, 1, ALPHA, B, 1, BETA, W, 1, INFO )
         CALL CHKXER( 'CGGRQF', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*     Test error exits for the CGS, CGV, CGX, and CXV paths.
*
      ELSE IF( LSAMEN( 3, PATH, 'CGS' ) .OR.
     $         LSAMEN( 3, PATH, 'CGV' ) .OR.
     $         LSAMEN( 3, PATH, 'CGX' ) .OR. LSAMEN( 3, PATH, 'CXV' ) )
     $          THEN
*
*        CGGES
*
         SRNAMT = 'CGGES '
         INFOT = 1
         CALL CGGES( '/', 'N', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $               BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGES( 'N', '/', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $               BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGES( 'N', 'V', '/', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $               BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGES( 'N', 'V', 'S', CLCTES, -1, A, 1, B, 1, SDIM, ALPHA,
     $               BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGGES( 'N', 'V', 'S', CLCTES, 1, A, 0, B, 1, SDIM, ALPHA,
     $               BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CGGES( 'N', 'V', 'S', CLCTES, 1, A, 1, B, 0, SDIM, ALPHA,
     $               BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL CGGES( 'N', 'V', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $               BETA, Q, 0, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL CGGES( 'V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA,
     $               BETA, Q, 1, U, 2, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL CGGES( 'N', 'V', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $               BETA, Q, 1, U, 0, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL CGGES( 'V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA,
     $               BETA, Q, 2, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL CGGES( 'V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA,
     $               BETA, Q, 2, U, 2, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES ', INFOT, NOUT, LERR, OK )
         NT = NT + 11
*
*        CGGES3
*
         SRNAMT = 'CGGES3'
         INFOT = 1
         CALL CGGES3( '/', 'N', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $                BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGES3( 'N', '/', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $                BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGES3( 'N', 'V', '/', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $                BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGES3( 'N', 'V', 'S', CLCTES, -1, A, 1, B, 1, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGGES3( 'N', 'V', 'S', CLCTES, 1, A, 0, B, 1, SDIM, ALPHA,
     $                BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CGGES3( 'N', 'V', 'S', CLCTES, 1, A, 1, B, 0, SDIM, ALPHA,
     $                BETA, Q, 1, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL CGGES3( 'N', 'V', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $                BETA, Q, 0, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL CGGES3( 'V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA,
     $                BETA, Q, 1, U, 2, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL CGGES3( 'N', 'V', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA,
     $                BETA, Q, 1, U, 0, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL CGGES3( 'V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA,
     $                BETA, Q, 2, U, 1, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL CGGES3( 'V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA,
     $                BETA, Q, 2, U, 2, W, 1, RW, BW, INFO )
         CALL CHKXER( 'CGGES3', INFOT, NOUT, LERR, OK )
         NT = NT + 11
*
*        CGGESX
*
         SRNAMT = 'CGGESX'
         INFOT = 1
         CALL CGGESX( '/', 'N', 'S', CLCTSX, 'N', 1, A, 1, B, 1, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGESX( 'N', '/', 'S', CLCTSX, 'N', 1, A, 1, B, 1, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGESX( 'V', 'V', '/', CLCTSX, 'N', 1, A, 1, B, 1, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, '/', 1, A, 1, B, 1, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, 'B', -1, A, 1, B, 1, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, 'B', 1, A, 0, B, 1, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, 'B', 1, A, 1, B, 0, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, 'B', 1, A, 1, B, 1, SDIM,
     $                ALPHA, BETA, Q, 0, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, 'B', 2, A, 2, B, 2, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, 'B', 1, A, 1, B, 1, SDIM,
     $                ALPHA, BETA, Q, 1, U, 0, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, 'B', 2, A, 2, B, 2, SDIM,
     $                ALPHA, BETA, Q, 2, U, 1, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 21
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, 'B', 2, A, 2, B, 2, SDIM,
     $                ALPHA, BETA, Q, 2, U, 2, RCE, RCV, W, 1, RW, IW,
     $                1, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 24
         CALL CGGESX( 'V', 'V', 'S', CLCTSX, 'V', 1, A, 1, B, 1, SDIM,
     $                ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 32, RW, IW,
     $                0, BW, INFO )
         CALL CHKXER( 'CGGESX', INFOT, NOUT, LERR, OK )
         NT = NT + 13
*
*        When the QZ iteration fails, CGGES, CGGES3 and CGGESX must
*        still return the eigenvalues that were computed in the units
*        of the input.  The leading 3-by-3 block of A is dense, the
*        entry (4,4) is isolated, and B has a NaN in its 3-by-3 block,
*        so the QZ iteration fails with INFO = 1, 2 or 3 and
*        ALPHA(4) and BETA(4) are correct for the isolated
*        eigenvalue.  The case J = 1 scales A up, J = 2 scales B up
*        (with the NaN in A) and J = 3 scales A down.  A and B
*        themselves are left unchanged.
*
         RONE = ONE
         RNAN = SQRT( -RONE )
         BIG = 1.0E+30
         TINY = 1.0E-30
         TOL = 1.0E-4
         DO 42 J = 1, 3
            DO 41 K = 1, 3
               DO 35 I = 1, NQ
                  DO 34 M = 1, NQ
                     AQ( I, M ) = ZERO
                     BQ( I, M ) = ZERO
   34             CONTINUE
   35          CONTINUE
               DO 37 I = 1, 3
                  DO 36 M = 1, 3
                     AQ( I, M ) = ONE
   36             CONTINUE
                  BQ( I, I ) = ONE
   37          CONTINUE
               AQ( 1, 1 ) = 2.0E+0
               AQ( 2, 2 ) = 3.0E+0
               AQ( 3, 3 ) = 4.0E+0
               AQ( 4, 4 ) = 5.0E+0
               BQ( 4, 4 ) = 2.0E+0
               IF( J.EQ.2 ) THEN
                  AQ( 2, 3 ) = RNAN
                  BQ( 1, 1 ) = BIG
                  BQ( 2, 2 ) = BIG
                  BQ( 3, 3 ) = BIG
                  BQ( 4, 4 ) = 2.0E+0*BIG
               ELSE
                  BQ( 2, 3 ) = RNAN
                  DO 39 I = 1, NQ
                     DO 38 M = 1, NQ
                        IF( J.EQ.1 ) THEN
                           AQ( I, M ) = AQ( I, M )*BIG
                        ELSE
                           AQ( I, M ) = AQ( I, M )*TINY
                        END IF
   38                CONTINUE
   39             CONTINUE
               END IF
               EALPHA = AQ( 4, 4 )
               EBETA = BQ( 4, 4 )
               INFOT = 0
               LERR = .FALSE.
               IF( K.EQ.1 ) THEN
                  SRNAMT = 'CGGES '
                  CALL CGGES( 'N', 'N', 'S', CLCTES, NQ, AQ, NQ, BQ,
     $                        NQ, SDIM, ALQ, BEQ, Q, 1, U, 1,
     $                        WQ, LWQ, RWQ, BWQ, INFO )
               ELSE IF( K.EQ.2 ) THEN
                  SRNAMT = 'CGGES3'
                  CALL CGGES3( 'N', 'N', 'S', CLCTES, NQ, AQ, NQ,
     $                         BQ, NQ, SDIM, ALQ, BEQ, Q, 1, U, 1,
     $                         WQ, LWQ, RWQ, BWQ, INFO )
               ELSE
                  SRNAMT = 'CGGESX'
                  CALL CGGESX( 'N', 'N', 'S', CLCTES, 'N', NQ, AQ,
     $                         NQ, BQ, NQ, SDIM, ALQ, BEQ, Q, 1,
     $                         U, 1, RCE, RCV, WQ, LWQ, RWQ, IW, 1, BWQ,
     $                         INFO )
               END IF
               IF( LERR .OR. INFO.LT.1 .OR. INFO.GT.3 .OR.
     $             .NOT.( ABS( ALQ( 4 )-EALPHA ).LE.
     $                    TOL*ABS( EALPHA ) ) .OR.
     $             .NOT.( ABS( BEQ( 4 )-EBETA ).LE.
     $                    TOL*ABS( EBETA ) ) ) THEN
                  WRITE( NOUT, FMT = 9997 )SRNAMT( 1: 6 ), J, INFO,
     $               REAL( ALQ( 4 ) ), REAL( BEQ( 4 ) )
                  OK = .FALSE.
               END IF
   41       CONTINUE
   42    CONTINUE
         NT = NT + 9
*
*        CGGEV
*
         SRNAMT = 'CGGEV '
         INFOT = 1
         CALL CGGEV( '/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGEV( 'N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGEV( 'V', 'V', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGEV( 'V', 'V', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGGEV( 'V', 'V', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGGEV( 'N', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGGEV( 'V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CGGEV( 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 0,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CGGEV( 'V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CGGEV( 'V', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $               W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        CGGEV3
*
         SRNAMT = 'CGGEV3'
         INFOT = 1
         CALL CGGEV3( '/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGEV3( 'N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGEV3( 'V', 'V', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGEV3( 'V', 'V', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGGEV3( 'V', 'V', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGGEV3( 'N', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGGEV3( 'V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CGGEV3( 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 0,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CGGEV3( 'V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CGGEV3( 'V', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1,
     $                W, 1, RW, INFO )
         CALL CHKXER( 'CGGEV3', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        CGGEVX
*
         SRNAMT = 'CGGEVX'
         INFOT = 1
         CALL CGGEVX( '/', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGGEVX( 'N', '/', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGGEVX( 'N', 'N', '/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGGEVX( 'N', 'N', 'N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGGEVX( 'N', 'N', 'N', 'N', -1, A, 1, B, 1, ALPHA, BETA,
     $                Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE,
     $                RCV, W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGGEVX( 'N', 'N', 'N', 'N', 1, A, 0, B, 1, ALPHA, BETA, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CGGEVX( 'N', 'N', 'N', 'N', 1, A, 1, B, 0, ALPHA, BETA, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CGGEVX( 'N', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q,
     $                0, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CGGEVX( 'N', 'V', 'N', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q,
     $                1, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CGGEVX( 'N', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q,
     $                1, U, 0, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CGGEVX( 'N', 'N', 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q,
     $                2, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 25
         CALL CGGEVX( 'N', 'N', 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q,
     $                2, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 0, RW, IW, BW, INFO )
         CALL CHKXER( 'CGGEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 12
*
*        CTGEXC
*
         SRNAMT = 'CTGEXC'
         INFOT = 3
         CALL CTGEXC( .TRUE., .TRUE., -1, A, 1, B, 1, Q, 1, Z, 1, IFST,
     $                ILST, INFO )
         CALL CHKXER( 'CTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CTGEXC( .TRUE., .TRUE., 1, A, 0, B, 1, Q, 1, Z, 1, IFST,
     $                ILST, INFO )
         CALL CHKXER( 'CTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CTGEXC( .TRUE., .TRUE., 1, A, 1, B, 0, Q, 1, Z, 1, IFST,
     $                ILST, INFO )
         CALL CHKXER( 'CTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CTGEXC( .FALSE., .TRUE., 1, A, 1, B, 1, Q, 0, Z, 1, IFST,
     $                ILST, INFO )
         CALL CHKXER( 'CTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CTGEXC( .TRUE., .TRUE., 1, A, 1, B, 1, Q, 0, Z, 1, IFST,
     $                ILST, INFO )
         CALL CHKXER( 'CTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CTGEXC( .TRUE., .FALSE., 1, A, 1, B, 1, Q, 1, Z, 0, IFST,
     $                ILST, INFO )
         CALL CHKXER( 'CTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CTGEXC( .TRUE., .TRUE., 1, A, 1, B, 1, Q, 1, Z, 0, IFST,
     $                ILST, INFO )
         CALL CHKXER( 'CTGEXC', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
*        CTGSEN
*
         SRNAMT = 'CTGSEN'
         INFOT = 1
         CALL CTGSEN( -1, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, ALPHA,
     $                BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CTGSEN( 1, .TRUE., .TRUE., SEL, -1, A, 1, B, 1, ALPHA,
     $                BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 0, B, 1, ALPHA,
     $                BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 1, B, 0, ALPHA,
     $                BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL CTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, ALPHA,
     $                BETA, Q, 0, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, ALPHA,
     $                BETA, Q, 1, Z, 0, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 21
         CALL CTGSEN( 3, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, ALPHA,
     $                BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, -5, IW,
     $                1, INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 23
         CALL CTGSEN( 0, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, ALPHA,
     $                BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW,
     $                0, INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 23
         CALL CTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, ALPHA,
     $                BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW,
     $                0, INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 23
         CALL CTGSEN( 5, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, ALPHA,
     $                BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW,
     $                1, INFO )
         CALL CHKXER( 'CTGSEN', INFOT, NOUT, LERR, OK )
         NT = NT + 11
*
*        CTGSNA
*
         SRNAMT = 'CTGSNA'
         INFOT = 1
         CALL CTGSNA( '/', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CTGSNA( 'B', '/', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CTGSNA( 'B', 'A', SEL, -1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CTGSNA( 'B', 'A', SEL, 1, A, 0, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CTGSNA( 'B', 'A', SEL, 1, A, 1, B, 0, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CTGSNA( 'E', 'A', SEL, 1, A, 1, B, 1, Q, 0, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CTGSNA( 'E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 0, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CTGSNA( 'E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                0, M, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL CTGSNA( 'E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 0, IW, INFO )
         CALL CHKXER( 'CTGSNA', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        CTGSYL
*
         SRNAMT = 'CTGSYL'
         INFOT = 1
         CALL CTGSYL( '/', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CTGSYL( 'N', -1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CTGSYL( 'N', 0, 0, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CTGSYL( 'N', 0, 1, 0, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CTGSYL( 'N', 0, 1, 1, A, 0, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CTGSYL( 'N', 0, 1, 1, A, 1, B, 0, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CTGSYL( 'N', 0, 1, 1, A, 1, B, 1, Q, 0, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CTGSYL( 'N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 0, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL CTGSYL( 'N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 0, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL CTGSYL( 'N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 0,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL CTGSYL( 'N', 1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL CTGSYL( 'N', 2, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'CTGSYL', INFOT, NOUT, LERR, OK )
         NT = NT + 12
      END IF
*
*     Print a summary line.
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH, NT
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (',
     $      I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )
 9997 FORMAT( ' *** ', A, ' did not return the computed eigenvalue ',
     $      '4 in the units of the input after a QZ failure (case ', I1,
     $      ', INFO = ', I3, ', ALPHA = ', 1P, E12.4, ', BETA = ',
     $      E12.4, ') ***' )
*
      RETURN
*
*     End of CERRGG
*
      END
