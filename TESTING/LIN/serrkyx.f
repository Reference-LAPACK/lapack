*> \brief \b SERRKYX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SERRKY( PATH, NUNIT )
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
*> SERRKY tests the error exits for the REAL routines
*> for symmetric indefinite matrices.
*>
*> Note that this file is used only when the XBLAS are available,
*> otherwise serrsy.f defines this subroutine.
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
*> \ingroup single_lin
*
*  =====================================================================
      SUBROUTINE SERRKY( PATH, NUNIT )
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
      INTEGER            NMAX
      PARAMETER          ( NMAX = 4 )
*     ..
*     .. Local Scalars ..
      CHARACTER          EQ
      CHARACTER*2        C2
      INTEGER            I, INFO, J, N_ERR_BNDS, NPARAMS
      REAL               ANRM, RCOND, BERR
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX ), IW( NMAX )
      REAL               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ),
     $                   E( NMAX ), R1( NMAX ), R2( NMAX ), W( 3*NMAX ),
     $                   X( NMAX ), S( NMAX ), ERR_BNDS_N( NMAX, 3 ),
     $                   ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, SKYTF2, SKYTRF,
     $                   SKYTRI, SKYTRI2, SKYTRI2X, SKYTRS
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
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
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
         DO 10 I = 1, NMAX
            A( I, J ) = 1. / REAL( I+J )
            AF( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
         B( J ) = 0.E+0
         E( J ) = 0.E+0
         R1( J ) = 0.E+0
         R2( J ) = 0.E+0
         W( J ) = 0.E+0
         X( J ) = 0.E+0
         IP( J ) = J
         IW( J ) = J
   20 CONTINUE
      ANRM = 1.0
      RCOND = 1.0
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'KY' ) ) THEN
*
*        Test error exits of the routines that use factorization
*        of a symmetric indefinite matrix with patrial
*        (Bunch-Kaufman) pivoting.
*
*        SKYTRF
*
         SRNAMT = 'SKYTRF'
         INFOT = 1
         CALL SKYTRF( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SKYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKYTRF( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SKYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SKYTRF( 'U', 2, A, 1, IP, W, 4, INFO )
         CALL CHKXER( 'SKYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SKYTRF( 'U', 0, A, 1, IP, W, 0, INFO )
         CALL CHKXER( 'SKYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL SKYTRF( 'U', 0, A, 1, IP, W, -2, INFO )
         CALL CHKXER( 'SKYTRF', INFOT, NOUT, LERR, OK )
*
*        SKYTF2
*
         SRNAMT = 'SKYTF2'
         INFOT = 1
         CALL SKYTF2( '/', 0, A, 1, IP, INFO )
         CALL CHKXER( 'SKYTF2', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKYTF2( 'U', -1, A, 1, IP, INFO )
         CALL CHKXER( 'SKYTF2', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SKYTF2( 'U', 2, A, 1, IP, INFO )
         CALL CHKXER( 'SKYTF2', INFOT, NOUT, LERR, OK )
*
*        SKYTRI
*
         SRNAMT = 'SKYTRI'
         INFOT = 1
         CALL SKYTRI( '/', 0, A, 1, IP, W, INFO )
         CALL CHKXER( 'SKYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKYTRI( 'U', -1, A, 1, IP, W, INFO )
         CALL CHKXER( 'SKYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SKYTRI( 'U', 2, A, 1, IP, W, INFO )
         CALL CHKXER( 'SKYTRI', INFOT, NOUT, LERR, OK )
*
*        SKYTRI2
*
         SRNAMT = 'SKYTRI2'
         INFOT = 1
         CALL SKYTRI2( '/', 0, A, 1, IP, W, IW, INFO )
         CALL CHKXER( 'SKYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKYTRI2( 'U', -1, A, 1, IP, W, IW, INFO )
         CALL CHKXER( 'SKYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SKYTRI2( 'U', 2, A, 1, IP, W, IW, INFO )
         CALL CHKXER( 'SKYTRI', INFOT, NOUT, LERR, OK )
*
*        SKYTRI2X
*
         SRNAMT = 'SKYTRI2X'
         INFOT = 1
         CALL SKYTRI2X( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SKYTRI2X', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKYTRI2X( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SKYTRI2X', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL SKYTRI2X( 'U', 2, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'SKYTRI2X', INFOT, NOUT, LERR, OK )
*
*        SKYTRS
*
         SRNAMT = 'SKYTRS'
         INFOT = 1
         CALL SKYTRS( '/', 0, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'SKYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL SKYTRS( 'U', -1, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'SKYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL SKYTRS( 'U', 0, -1, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'SKYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL SKYTRS( 'U', 2, 1, A, 1, IP, B, 2, INFO )
         CALL CHKXER( 'SKYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL SKYTRS( 'U', 2, 1, A, 2, IP, B, 1, INFO )
         CALL CHKXER( 'SKYTRS', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of SERRKYX
*
      END
