*> \brief \b DERRGG
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DERRGG( PATH, NUNIT )
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
*> DERRGG tests the error exits for DGGES, DGGESX, DGGEV,  DGGEVX,
*> DGGGLM, DGGHRD, DGGLSE, DGGQRF, DGGRQF, DGGSVD3,
*> DGGSVP3, DHGEQZ, DORBDB, DORCSD, DTGEVC, DTGEXC, DTGSEN,
*> DTGSJA, DTGSNA, DGGES3, DGGEV3, and DTGSYL.
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
*> \ingroup double_eig
*
*  =====================================================================
      SUBROUTINE DERRGG( PATH, NUNIT )
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
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            DUMMYK, DUMMYL, I, IFST, ILO, IHI, ILST, INFO,
     $                   J, K, M, NCYCLE, NT, SDIM, LWORK
      DOUBLE PRECISION   ANRM, BNRM, DIF, SCALE, TOLA, TOLB
      DOUBLE PRECISION   BIG, EALPHA, EBETA, RNAN, RONE, TINY, TOL
*     ..
*     .. Local Arrays ..
      LOGICAL            BW( NMAX ), SEL( NMAX )
      INTEGER            IW( NMAX ), IDUM(NMAX)
      DOUBLE PRECISION   A( NMAX, NMAX ), B( NMAX, NMAX ), LS( NMAX ),
     $                   Q( NMAX, NMAX ), R1( NMAX ), R2( NMAX ),
     $                   R3( NMAX ), RCE( 2 ), RCV( 2 ), RS( NMAX ),
     $                   TAU( NMAX ), U( NMAX, NMAX ), V( NMAX, NMAX ),
     $                   W( LW ), Z( NMAX, NMAX )
      LOGICAL            BWQ( NQ )
      DOUBLE PRECISION   AQ( NQ, NQ ), ARQ( NQ ), AIQ( NQ ), BEQ( NQ ),
     $                   BQ( NQ, NQ ), WQ( LWQ )
*     ..
*     .. External Functions ..
      LOGICAL            DLCTES, DLCTSX, LSAMEN
      EXTERNAL           DLCTES, DLCTSX, LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHKXER, DGGES, DGGESX, DGGEV, DGGEVX, DGGGLM,
     $                   DGGHRD, DGGLSE, DGGQRF, DGGRQF,
     $                   DHGEQZ, DORBDB, DORCSD, DTGEVC, DTGEXC, DTGSEN,
     $                   DTGSJA, DTGSNA, DTGSYL, DGGHD3, DGGES3, DGGEV3,
     $                   DGGSVD3, DGGSVP3, XLAENV
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
*        DGGHRD
*
         SRNAMT = 'DGGHRD'
         INFOT = 1
         CALL DGGHRD( '/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'DGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGHRD( 'N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'DGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGHRD( 'N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'DGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGGHRD( 'N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'DGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGHRD( 'N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'DGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGGHRD( 'N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'DGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DGGHRD( 'N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'DGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL DGGHRD( 'V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'DGGHRD', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL DGGHRD( 'N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO )
         CALL CHKXER( 'DGGHRD', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        DGGHD3
*
         SRNAMT = 'DGGHD3'
         INFOT = 1
         CALL DGGHD3( '/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'DGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGHD3( 'N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'DGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGHD3( 'N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'DGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGGHD3( 'N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'DGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGHD3( 'N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'DGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGGHD3( 'N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'DGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DGGHD3( 'N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'DGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL DGGHD3( 'V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'DGGHD3', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL DGGHD3( 'N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW,
     $                INFO )
         CALL CHKXER( 'DGGHD3', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        DHGEQZ
*
         SRNAMT = 'DHGEQZ'
         INFOT = 1
         CALL DHGEQZ( '/', 'N', 'N', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q,
     $                1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DHGEQZ( 'E', '/', 'N', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q,
     $                1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DHGEQZ( 'E', 'N', '/', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q,
     $                1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DHGEQZ( 'E', 'N', 'N', -1, 0, 0, A, 1, B, 1, R1, R2, R3,
     $                Q, 1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DHGEQZ( 'E', 'N', 'N', 0, 0, 0, A, 1, B, 1, R1, R2, R3, Q,
     $                1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DHGEQZ( 'E', 'N', 'N', 0, 1, 1, A, 1, B, 1, R1, R2, R3, Q,
     $                1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DHGEQZ( 'E', 'N', 'N', 2, 1, 1, A, 1, B, 2, R1, R2, R3, Q,
     $                1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DHGEQZ( 'E', 'N', 'N', 2, 1, 1, A, 2, B, 1, R1, R2, R3, Q,
     $                1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL DHGEQZ( 'E', 'V', 'N', 2, 1, 1, A, 2, B, 2, R1, R2, R3, Q,
     $                1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL DHGEQZ( 'E', 'N', 'V', 2, 1, 1, A, 2, B, 2, R1, R2, R3, Q,
     $                1, Z, 1, W, LW, INFO )
         CALL CHKXER( 'DHGEQZ', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        DTGEVC
*
         SRNAMT = 'DTGEVC'
         INFOT = 1
         CALL DTGEVC( '/', 'A', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W,
     $                INFO )
         CALL CHKXER( 'DTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DTGEVC( 'R', '/', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W,
     $                INFO )
         CALL CHKXER( 'DTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DTGEVC( 'R', 'A', SEL, -1, A, 1, B, 1, Q, 1, Z, 1, 0, M,
     $                W, INFO )
         CALL CHKXER( 'DTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DTGEVC( 'R', 'A', SEL, 2, A, 1, B, 2, Q, 1, Z, 2, 0, M, W,
     $                INFO )
         CALL CHKXER( 'DTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DTGEVC( 'R', 'A', SEL, 2, A, 2, B, 1, Q, 1, Z, 2, 0, M, W,
     $                INFO )
         CALL CHKXER( 'DTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DTGEVC( 'L', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W,
     $                INFO )
         CALL CHKXER( 'DTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DTGEVC( 'R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W,
     $                INFO )
         CALL CHKXER( 'DTGEVC', INFOT, NOUT, LERR, OK )
         INFOT = 13
         CALL DTGEVC( 'R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 2, 1, M, W,
     $                INFO )
         CALL CHKXER( 'DTGEVC', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*     Test error exits for the GSV path.
*
      ELSE IF( LSAMEN( 3, PATH, 'GSV' ) ) THEN
*
*        DGGSVD3
*
         SRNAMT = 'DGGSVD3'
         INFOT = 1
         CALL DGGSVD3( '/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGSVD3( 'N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGSVD3( 'N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGGSVD3( 'N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGSVD3( 'N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DGGSVD3( 'N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DGGSVD3( 'N', 'N', 'N', 2, 1, 1, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGGSVD3( 'N', 'N', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DGGSVD3( 'U', 'N', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B,
     $               2, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL DGGSVD3( 'N', 'V', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B,
     $               2, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL DGGSVD3( 'N', 'N', 'Q', 1, 2, 1, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL DGGSVD3( 'N', 'N', 'N', 1, 1, 1, DUMMYK, DUMMYL, A, 1, B,
     $               1, R1, R2, U, 1, V, 1, Q, 1, W, 1, IDUM, INFO )
         CALL CHKXER( 'DGGSVD3', INFOT, NOUT, LERR, OK )
         NT = NT + 12
*
*        DGGSVP3
*
         SRNAMT = 'DGGSVP3'
         INFOT = 1
         CALL DGGSVP3( '/', 'N', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGSVP3( 'N', '/', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGSVP3( 'N', 'N', '/', 0, 0, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGGSVP3( 'N', 'N', 'N', -1, 0, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGSVP3( 'N', 'N', 'N', 0, -1, 0, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DGGSVP3( 'N', 'N', 'N', 0, 0, -1, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGGSVP3( 'N', 'N', 'N', 2, 1, 1, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DGGSVP3( 'N', 'N', 'N', 1, 2, 1, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DGGSVP3( 'U', 'N', 'N', 2, 2, 2, A, 2, B, 2, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL DGGSVP3( 'N', 'V', 'N', 1, 2, 1, A, 1, B, 2, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL DGGSVP3( 'N', 'N', 'Q', 1, 1, 2, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 LWORK, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         INFOT = 24
         CALL DGGSVP3( 'N', 'N', 'N', 1, 1, 1, A, 1, B, 1, TOLA, TOLB,
     $                 DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W,
     $                 1, INFO )
         CALL CHKXER( 'DGGSVP3', INFOT, NOUT, LERR, OK )
         NT = NT + 12
*
*        DTGSJA
*
         SRNAMT = 'DTGSJA'
         INFOT = 1
         CALL DTGSJA( '/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DTGSJA( 'N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DTGSJA( 'N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DTGSJA( 'N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DTGSJA( 'N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DTGSJA( 'N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DTGSJA( 'N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 0, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DTGSJA( 'N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                0, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL DTGSJA( 'U', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 0, V, 1, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL DTGSJA( 'N', 'V', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 0, Q, 1, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL DTGSJA( 'N', 'N', 'Q', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B,
     $                1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 0, W,
     $                NCYCLE, INFO )
         CALL CHKXER( 'DTGSJA', INFOT, NOUT, LERR, OK )
         NT = NT + 11
*
*     Test error exits for the GLM path.
*
      ELSE IF( LSAMEN( 3, PATH, 'GLM' ) ) THEN
*
*        DGGGLM
*
         SRNAMT = 'DGGGLM'
         INFOT = 1
         CALL DGGGLM( -1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGGLM( 0, -1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGGLM( 0, 1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGGLM( 0, 0, -1, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGGLM( 1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGGLM( 0, 0, 0, A, 0, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGGGLM( 0, 0, 0, A, 1, B, 0, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGGLM', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGGGLM( 1, 1, 1, A, 1, B, 1, R1, R2, R3, W, 1, INFO )
         CALL CHKXER( 'DGGGLM', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*     Test error exits for the LSE path.
*
      ELSE IF( LSAMEN( 3, PATH, 'LSE' ) ) THEN
*
*        DGGLSE
*
         SRNAMT = 'DGGLSE'
         INFOT = 1
         CALL DGGLSE( -1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGLSE( 0, -1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGLSE( 0, 0, -1, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGLSE( 0, 0, 1, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGLSE( 0, 1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGLSE( 0, 0, 0, A, 0, B, 1, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGGLSE( 0, 0, 0, A, 1, B, 0, R1, R2, R3, W, LW, INFO )
         CALL CHKXER( 'DGGLSE', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGGLSE( 1, 1, 1, A, 1, B, 1, R1, R2, R3, W, 1, INFO )
         CALL CHKXER( 'DGGLSE', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*     Test error exits for the CSD path.
*
      ELSE IF( LSAMEN( 3, PATH, 'CSD' ) ) THEN
*
*        DORCSD
*
         SRNAMT = 'DORCSD'
         INFOT = 7
         CALL DORCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 -1, 0, 0, A, 1, A,
     $                 1, A, 1, A, 1, A,
     $                 A, 1, A, 1, A, 1, A,
     $                 1, W, LW, IW, INFO )
         CALL CHKXER( 'DORCSD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DORCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, -1, 0, A, 1, A,
     $                 1, A, 1, A, 1, A,
     $                 A, 1, A, 1, A, 1, A,
     $                 1, W, LW, IW, INFO )
         CALL CHKXER( 'DORCSD', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DORCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, -1, A, 1, A,
     $                 1, A, 1, A, 1, A,
     $                 A, 1, A, 1, A, 1, A,
     $                 1, W, LW, IW, INFO )
         CALL CHKXER( 'DORCSD', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL DORCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, -1, A,
     $                 1, A, 1, A, 1, A,
     $                 A, 1, A, 1, A, 1, A,
     $                 1, W, LW, IW, INFO )
         CALL CHKXER( 'DORCSD', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL DORCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, 1, A,
     $                 1, A, 1, A, 1, A,
     $                 A, -1, A, 1, A, 1, A,
     $                 1, W, LW, IW, INFO )
         CALL CHKXER( 'DORCSD', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL DORCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, 1, A,
     $                 1, A, 1, A, 1, A,
     $                 A, 1, A, -1, A, 1, A,
     $                 1, W, LW, IW, INFO )
         CALL CHKXER( 'DORCSD', INFOT, NOUT, LERR, OK )
         INFOT = 24
         CALL DORCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, 1, A,
     $                 1, A, 1, A, 1, A,
     $                 A, 1, A, 1, A, -1, A,
     $                 1, W, LW, IW, INFO )
         CALL CHKXER( 'DORCSD', INFOT, NOUT, LERR, OK )
         INFOT = 26
         CALL DORCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'N',
     $                 1, 1, 1, A, 1, A,
     $                 1, A, 1, A, 1, A,
     $                 A, 1, A, 1, A, 1, A,
     $                 -1, W, LW, IW, INFO )
         CALL CHKXER( 'DORCSD', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        DORBDB
*
         SRNAMT = 'DORBDB'
         INFOT = 3
         CALL DORBDB( 'N', 'N', -1, 0, 0, A, 1, A, 1, A, 1, A, 1,
     $                A, A, A, A, A, A, W, LW, INFO )
         CALL CHKXER( 'DORBDB', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DORBDB( 'N', 'N', 1, -1, 0, A, 1, A, 1, A, 1, A, 1,
     $                A, A, A, A, A, A, W, LW, INFO )
         CALL CHKXER( 'DORBDB', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DORBDB( 'N', 'N', 1, 1, -1, A, 1, A, 1, A, 1, A, 1,
     $                A, A, A, A, A, A, W, LW, INFO )
         CALL CHKXER( 'DORBDB', INFOT, NOUT, LERR, OK )
         NT = NT + 3
*
*     Test error exits for the GQR path.
*
      ELSE IF( LSAMEN( 3, PATH, 'GQR' ) ) THEN
*
*        DGGQRF
*
         SRNAMT = 'DGGQRF'
         INFOT = 1
         CALL DGGQRF( -1, 0, 0, A, 1, R1, B, 1, R2, W, LW, INFO )
         CALL CHKXER( 'DGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGQRF( 0, -1, 0, A, 1, R1, B, 1, R2, W, LW, INFO )
         CALL CHKXER( 'DGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGQRF( 0, 0, -1, A, 1, R1, B, 1, R2, W, LW, INFO )
         CALL CHKXER( 'DGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGQRF( 0, 0, 0, A, 0, R1, B, 1, R2, W, LW, INFO )
         CALL CHKXER( 'DGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGGQRF( 0, 0, 0, A, 1, R1, B, 0, R2, W, LW, INFO )
         CALL CHKXER( 'DGGQRF', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL DGGQRF( 1, 1, 2, A, 1, R1, B, 1, R2, W, 1, INFO )
         CALL CHKXER( 'DGGQRF', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*        DGGRQF
*
         SRNAMT = 'DGGRQF'
         INFOT = 1
         CALL DGGRQF( -1, 0, 0, A, 1, R1, B, 1, R2, W, LW, INFO )
         CALL CHKXER( 'DGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGRQF( 0, -1, 0, A, 1, R1, B, 1, R2, W, LW, INFO )
         CALL CHKXER( 'DGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGRQF( 0, 0, -1, A, 1, R1, B, 1, R2, W, LW, INFO )
         CALL CHKXER( 'DGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGRQF( 0, 0, 0, A, 0, R1, B, 1, R2, W, LW, INFO )
         CALL CHKXER( 'DGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGGRQF( 0, 0, 0, A, 1, R1, B, 0, R2, W, LW, INFO )
         CALL CHKXER( 'DGGRQF', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL DGGRQF( 1, 1, 2, A, 1, R1, B, 1, R2, W, 1, INFO )
         CALL CHKXER( 'DGGRQF', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
*     Test error exits for the DGS, DGV, DGX, and DXV paths.
*
      ELSE IF( LSAMEN( 3, PATH, 'DGS' ) .OR.
     $         LSAMEN( 3, PATH, 'DGV' ) .OR.
     $         LSAMEN( 3, PATH, 'DGX' ) .OR. LSAMEN( 3, PATH, 'DXV' ) )
     $          THEN
*
*        DGGES
*
         SRNAMT = 'DGGES '
         INFOT = 1
         CALL DGGES( '/', 'N', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2,
     $               R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGES( 'N', '/', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2,
     $               R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGES( 'N', 'V', '/', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2,
     $               R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGES( 'N', 'V', 'S', DLCTES, -1, A, 1, B, 1, SDIM, R1,
     $               R2, R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGGES( 'N', 'V', 'S', DLCTES, 1, A, 0, B, 1, SDIM, R1, R2,
     $               R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DGGES( 'N', 'V', 'S', DLCTES, 1, A, 1, B, 0, SDIM, R1, R2,
     $               R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL DGGES( 'N', 'V', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2,
     $               R3, Q, 0, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL DGGES( 'V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1, R2,
     $               R3, Q, 1, U, 2, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL DGGES( 'N', 'V', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2,
     $               R3, Q, 1, U, 0, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL DGGES( 'V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1, R2,
     $               R3, Q, 2, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         INFOT = 19
         CALL DGGES( 'V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1, R2,
     $               R3, Q, 2, U, 2, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES ', INFOT, NOUT, LERR, OK )
         NT = NT + 11
*
*        DGGES3
*
         SRNAMT = 'DGGES3 '
         INFOT = 1
         CALL DGGES3( '/', 'N', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1,
     $                R2, R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGES3( 'N', '/', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1,
     $                R2, R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGES3( 'N', 'V', '/', DLCTES, 1, A, 1, B, 1, SDIM, R1,
     $                R2, R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGES3( 'N', 'V', 'S', DLCTES, -1, A, 1, B, 1, SDIM, R1,
     $                R2, R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGGES3( 'N', 'V', 'S', DLCTES, 1, A, 0, B, 1, SDIM, R1,
     $                R2, R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DGGES3( 'N', 'V', 'S', DLCTES, 1, A, 1, B, 0, SDIM, R1,
     $                R2, R3, Q, 1, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL DGGES3( 'N', 'V', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1,
     $                R2, R3, Q, 0, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL DGGES3( 'V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1,
     $                R2, R3, Q, 1, U, 2, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL DGGES3( 'N', 'V', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1,
     $                R2, R3, Q, 1, U, 0, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 17
         CALL DGGES3( 'V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1,
     $                R2, R3, Q, 2, U, 1, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         INFOT = 19
         CALL DGGES3( 'V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1,
     $                R2, R3, Q, 2, U, 2, W, 1, BW, INFO )
         CALL CHKXER( 'DGGES3 ', INFOT, NOUT, LERR, OK )
         NT = NT + 11
*
*        DGGESX
*
         SRNAMT = 'DGGESX'
         INFOT = 1
         CALL DGGESX( '/', 'N', 'S', DLCTSX, 'N', 1, A, 1, B, 1, SDIM,
     $                R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGESX( 'N', '/', 'S', DLCTSX, 'N', 1, A, 1, B, 1, SDIM,
     $                R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGESX( 'V', 'V', '/', DLCTSX, 'N', 1, A, 1, B, 1, SDIM,
     $                R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, '/', 1, A, 1, B, 1, SDIM,
     $                R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, 'B', -1, A, 1, B, 1, SDIM,
     $                R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, 'B', 1, A, 0, B, 1, SDIM,
     $                R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, 'B', 1, A, 1, B, 0, SDIM,
     $                R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, 'B', 1, A, 1, B, 1, SDIM,
     $                R1, R2, R3, Q, 0, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, 'B', 2, A, 2, B, 2, SDIM,
     $                R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, 'B', 1, A, 1, B, 1, SDIM,
     $                R1, R2, R3, Q, 1, U, 0, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, 'B', 2, A, 2, B, 2, SDIM,
     $                R1, R2, R3, Q, 2, U, 1, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, 'B', 2, A, 2, B, 2, SDIM,
     $                R1, R2, R3, Q, 2, U, 2, RCE, RCV, W, 1, IW, 1, BW,
     $                INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         INFOT = 24
         CALL DGGESX( 'V', 'V', 'S', DLCTSX, 'V', 1, A, 1, B, 1, SDIM,
     $                R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 32, IW, 0,
     $                BW, INFO )
         CALL CHKXER( 'DGGESX', INFOT, NOUT, LERR, OK )
         NT = NT + 13
*
*        When the QZ iteration fails, DGGES, DGGES3 and DGGESX must
*        still return the eigenvalues that were computed in the units
*        of the input.  The leading 3-by-3 block of A is dense, the
*        entry (4,4) is isolated, and B has a NaN in its 3-by-3 block,
*        so the QZ iteration fails with INFO = 1, 2 or 3 and
*        ALPHAR(4), ALPHAI(4) and BETA(4) are correct for the isolated
*        eigenvalue.  The case J = 1 scales A up, J = 2 scales B up
*        (with the NaN in A) and J = 3 scales A down.  A and B
*        themselves are left unchanged.
*
         RONE = ONE
         RNAN = SQRT( -RONE )
         BIG = 1.0D+300
         TINY = 1.0D-300
         TOL = 1.0D-10
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
               AQ( 1, 1 ) = 2.0D+0
               AQ( 2, 2 ) = 3.0D+0
               AQ( 3, 3 ) = 4.0D+0
               AQ( 4, 4 ) = 5.0D+0
               BQ( 4, 4 ) = 2.0D+0
               IF( J.EQ.2 ) THEN
                  AQ( 2, 3 ) = RNAN
                  BQ( 1, 1 ) = BIG
                  BQ( 2, 2 ) = BIG
                  BQ( 3, 3 ) = BIG
                  BQ( 4, 4 ) = 2.0D+0*BIG
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
                  SRNAMT = 'DGGES '
                  CALL DGGES( 'N', 'N', 'S', DLCTES, NQ, AQ, NQ, BQ,
     $                        NQ, SDIM, ARQ, AIQ, BEQ, Q, 1, U, 1,
     $                        WQ, LWQ, BWQ, INFO )
               ELSE IF( K.EQ.2 ) THEN
                  SRNAMT = 'DGGES3'
                  CALL DGGES3( 'N', 'N', 'S', DLCTES, NQ, AQ, NQ,
     $                         BQ, NQ, SDIM, ARQ, AIQ, BEQ, Q, 1, U, 1,
     $                         WQ, LWQ, BWQ, INFO )
               ELSE
                  SRNAMT = 'DGGESX'
                  CALL DGGESX( 'N', 'N', 'S', DLCTES, 'N', NQ, AQ,
     $                         NQ, BQ, NQ, SDIM, ARQ, AIQ, BEQ, Q, 1,
     $                         U, 1, RCE, RCV, WQ, LWQ, IW, 1, BWQ,
     $                         INFO )
               END IF
               IF( LERR .OR. INFO.LT.1 .OR. INFO.GT.3 .OR.
     $             AIQ( 4 ).NE.ZERO .OR.
     $             .NOT.( ABS( ARQ( 4 )-EALPHA ).LE.
     $                    TOL*ABS( EALPHA ) ) .OR.
     $             .NOT.( ABS( BEQ( 4 )-EBETA ).LE.
     $                    TOL*ABS( EBETA ) ) ) THEN
                  WRITE( NOUT, FMT = 9997 )SRNAMT( 1: 6 ), J, INFO,
     $               ARQ( 4 ), BEQ( 4 )
                  OK = .FALSE.
               END IF
   41       CONTINUE
   42    CONTINUE
         NT = NT + 9
*
*        DGGEV
*
         SRNAMT = 'DGGEV '
         INFOT = 1
         CALL DGGEV( '/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W,
     $               1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGEV( 'N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W,
     $               1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGEV( 'V', 'V', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1,
     $               W, 1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGEV( 'V', 'V', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, W,
     $               1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGGEV( 'V', 'V', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, W,
     $               1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGGEV( 'N', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, W,
     $               1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGGEV( 'V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, W,
     $               1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL DGGEV( 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 0, W,
     $               1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL DGGEV( 'V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, W,
     $               1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DGGEV( 'V', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W,
     $               1, INFO )
         CALL CHKXER( 'DGGEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        DGGEV3
*
         CALL XLAENV( 12, 20 )
         CALL XLAENV( 13, 4 )
         CALL XLAENV( 14, 13 )
         CALL XLAENV( 15, 2 )
         CALL XLAENV( 17, 10 )
         SRNAMT = 'DGGEV3 '
         INFOT = 1
         CALL DGGEV3( '/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGEV3( 'N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGEV3( 'V', 'V', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGEV3( 'V', 'V', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGGEV3( 'V', 'V', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGGEV3( 'N', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DGGEV3( 'V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL DGGEV3( 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 0,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL DGGEV3( 'V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DGGEV3( 'V', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1,
     $                W, 1, INFO )
         CALL CHKXER( 'DGGEV3 ', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
*        DGGEVX
*
         SRNAMT = 'DGGEVX'
         INFOT = 1
         CALL DGGEVX( '/', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DGGEVX( 'N', '/', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DGGEVX( 'N', 'N', '/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DGGEVX( 'N', 'N', 'N', '/', 1, A, 1, B, 1, R1, R2, R3, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DGGEVX( 'N', 'N', 'N', 'N', -1, A, 1, B, 1, R1, R2, R3, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DGGEVX( 'N', 'N', 'N', 'N', 1, A, 0, B, 1, R1, R2, R3, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DGGEVX( 'N', 'N', 'N', 'N', 1, A, 1, B, 0, R1, R2, R3, Q,
     $                1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL DGGEVX( 'N', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q,
     $                0, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL DGGEVX( 'N', 'V', 'N', 'N', 2, A, 2, B, 2, R1, R2, R3, Q,
     $                1, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DGGEVX( 'N', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q,
     $                1, U, 0, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DGGEVX( 'N', 'N', 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q,
     $                2, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         INFOT = 26
         CALL DGGEVX( 'N', 'N', 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q,
     $                2, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV,
     $                W, 1, IW, BW, INFO )
         CALL CHKXER( 'DGGEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 12
*
*        DTGEXC
*
         SRNAMT = 'DTGEXC'
         INFOT = 3
         CALL DTGEXC( .TRUE., .TRUE., -1, A, 1, B, 1, Q, 1, Z, 1, IFST,
     $                ILST, W, 1, INFO )
         CALL CHKXER( 'DTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DTGEXC( .TRUE., .TRUE., 1, A, 0, B, 1, Q, 1, Z, 1, IFST,
     $                ILST, W, 1, INFO )
         CALL CHKXER( 'DTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DTGEXC( .TRUE., .TRUE., 1, A, 1, B, 0, Q, 1, Z, 1, IFST,
     $                ILST, W, 1, INFO )
         CALL CHKXER( 'DTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DTGEXC( .FALSE., .TRUE., 1, A, 1, B, 1, Q, 0, Z, 1, IFST,
     $                ILST, W, 1, INFO )
         CALL CHKXER( 'DTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DTGEXC( .TRUE., .TRUE., 1, A, 1, B, 1, Q, 0, Z, 1, IFST,
     $                ILST, W, 1, INFO )
         CALL CHKXER( 'DTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL DTGEXC( .TRUE., .FALSE., 1, A, 1, B, 1, Q, 1, Z, 0, IFST,
     $                ILST, W, 1, INFO )
         CALL CHKXER( 'DTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL DTGEXC( .TRUE., .TRUE., 1, A, 1, B, 1, Q, 1, Z, 0, IFST,
     $                ILST, W, 1, INFO )
         CALL CHKXER( 'DTGEXC', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL DTGEXC( .TRUE., .TRUE., 1, A, 1, B, 1, Q, 1, Z, 1, IFST,
     $                ILST, W, 0, INFO )
         CALL CHKXER( 'DTGEXC', INFOT, NOUT, LERR, OK )
         NT = NT + 8
*
*        DTGSEN
*
         SRNAMT = 'DTGSEN'
         INFOT = 1
         CALL DTGSEN( -1, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, R1, R2,
     $                R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DTGSEN( 1, .TRUE., .TRUE., SEL, -1, A, 1, B, 1, R1, R2,
     $                R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 0, B, 1, R1, R2, R3,
     $                Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL DTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 1, B, 0, R1, R2, R3,
     $                Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL DTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, R1, R2, R3,
     $                Q, 0, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, R1, R2, R3,
     $                Q, 1, Z, 0, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL DTGSEN( 0, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, R1, R2, R3,
     $                Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL DTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, R1, R2, R3,
     $                Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 22
         CALL DTGSEN( 2, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, R1, R2, R3,
     $                Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 24
         CALL DTGSEN( 0, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, R1, R2, R3,
     $                Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 24
         CALL DTGSEN( 1, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, R1, R2, R3,
     $                Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         INFOT = 24
         CALL DTGSEN( 2, .TRUE., .TRUE., SEL, 1, A, 1, B, 1, R1, R2, R3,
     $                Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 1,
     $                INFO )
         CALL CHKXER( 'DTGSEN', INFOT, NOUT, LERR, OK )
         NT = NT + 12
*
*        DTGSNA
*
         SRNAMT = 'DTGSNA'
         INFOT = 1
         CALL DTGSNA( '/', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DTGSNA( 'B', '/', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DTGSNA( 'B', 'A', SEL, -1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DTGSNA( 'B', 'A', SEL, 1, A, 0, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DTGSNA( 'B', 'A', SEL, 1, A, 1, B, 0, Q, 1, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DTGSNA( 'E', 'A', SEL, 1, A, 1, B, 1, Q, 0, U, 1, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DTGSNA( 'E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 0, R1, R2,
     $                1, M, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL DTGSNA( 'E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                0, M, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSNA', INFOT, NOUT, LERR, OK )
         INFOT = 18
         CALL DTGSNA( 'E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2,
     $                1, M, W, 0, IW, INFO )
         CALL CHKXER( 'DTGSNA', INFOT, NOUT, LERR, OK )
         NT = NT + 9
*
*        DTGSYL
*
         SRNAMT = 'DTGSYL'
         INFOT = 1
         CALL DTGSYL( '/', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DTGSYL( 'N', -1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DTGSYL( 'N', 0, 0, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DTGSYL( 'N', 0, 1, 0, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL DTGSYL( 'N', 0, 1, 1, A, 0, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DTGSYL( 'N', 0, 1, 1, A, 1, B, 0, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL DTGSYL( 'N', 0, 1, 1, A, 1, B, 1, Q, 0, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL DTGSYL( 'N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 0, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 14
         CALL DTGSYL( 'N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 0, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 16
         CALL DTGSYL( 'N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 0,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL DTGSYL( 'N', 1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL DTGSYL( 'N', 2, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1,
     $                SCALE, DIF, W, 1, IW, INFO )
         CALL CHKXER( 'DTGSYL', INFOT, NOUT, LERR, OK )
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
*     End of DERRGG
*
      END
