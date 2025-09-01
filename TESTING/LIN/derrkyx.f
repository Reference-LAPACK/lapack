*> \brief \b DERRKYX
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
*> SERRKY tests the error exits for the DOUBLE PRECISION routines
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
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE SERRKY( PATH, NUNIT )
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
      INTEGER            NMAX
      PARAMETER          ( NMAX = 4 )
*     ..
*     .. Local Scalars ..
      CHARACTER          EQ
      CHARACTER*2        C2
      INTEGER            I, INFO, J, N_ERR_BNDS, NPARAMS
      DOUBLE PRECISION   ANRM, RCOND, BERR
*     ..
*     .. Local Arrays ..
      INTEGER            IP( NMAX ), IW( NMAX )
      DOUBLE PRECISION   A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ),
     $                   E( NMAX ), R1( NMAX ), R2( NMAX ), W( 3*NMAX ),
     $                   X( NMAX ), S( NMAX ), ERR_BNDS_N( NMAX, 3 ),
     $                   ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DKYTF2, DKYTRF,
     $                   DKYTRI, DKYTRI2, DKYTRI2X, DKYTRS
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
      INTRINSIC          DOUBLE PRECISION
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
            A( I, J ) = 1. / DOUBLE PRECISION( I+J )
            AF( I, J ) = 1. / DOUBLE PRECISION( I+J )
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
*        DKYTRF
*
         SRNAMT = 'DKYTRF'
         INFOT = 1
         CALL DKYTRF( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DKYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKYTRF( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DKYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DKYTRF( 'U', 2, A, 1, IP, W, 4, INFO )
         CALL CHKXER( 'DKYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DKYTRF( 'U', 0, A, 1, IP, W, 0, INFO )
         CALL CHKXER( 'DKYTRF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DKYTRF( 'U', 0, A, 1, IP, W, -2, INFO )
         CALL CHKXER( 'DKYTRF', INFOT, NOUT, LERR, OK )
*
*        DKYTF2
*
         SRNAMT = 'DKYTF2'
         INFOT = 1
         CALL DKYTF2( '/', 0, A, 1, IP, INFO )
         CALL CHKXER( 'DKYTF2', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKYTF2( 'U', -1, A, 1, IP, INFO )
         CALL CHKXER( 'DKYTF2', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DKYTF2( 'U', 2, A, 1, IP, INFO )
         CALL CHKXER( 'DKYTF2', INFOT, NOUT, LERR, OK )
*
*        DKYTRI
*
         SRNAMT = 'DKYTRI'
         INFOT = 1
         CALL DKYTRI( '/', 0, A, 1, IP, W, INFO )
         CALL CHKXER( 'DKYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKYTRI( 'U', -1, A, 1, IP, W, INFO )
         CALL CHKXER( 'DKYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DKYTRI( 'U', 2, A, 1, IP, W, INFO )
         CALL CHKXER( 'DKYTRI', INFOT, NOUT, LERR, OK )
*
*        DKYTRI2
*
         SRNAMT = 'DKYTRI2'
         INFOT = 1
         CALL DKYTRI2( '/', 0, A, 1, IP, W, IW, INFO )
         CALL CHKXER( 'DKYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKYTRI2( 'U', -1, A, 1, IP, W, IW, INFO )
         CALL CHKXER( 'DKYTRI', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DKYTRI2( 'U', 2, A, 1, IP, W, IW, INFO )
         CALL CHKXER( 'DKYTRI', INFOT, NOUT, LERR, OK )
*
*        DKYTRI2X
*
         SRNAMT = 'DKYTRI2X'
         INFOT = 1
         CALL DKYTRI2X( '/', 0, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DKYTRI2X', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKYTRI2X( 'U', -1, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DKYTRI2X', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DKYTRI2X( 'U', 2, A, 1, IP, W, 1, INFO )
         CALL CHKXER( 'DKYTRI2X', INFOT, NOUT, LERR, OK )
*
*        DKYTRS
*
         SRNAMT = 'DKYTRS'
         INFOT = 1
         CALL DKYTRS( '/', 0, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'DKYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DKYTRS( 'U', -1, 0, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'DKYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL DKYTRS( 'U', 0, -1, A, 1, IP, B, 1, INFO )
         CALL CHKXER( 'DKYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL DKYTRS( 'U', 2, 1, A, 1, IP, B, 2, INFO )
         CALL CHKXER( 'DKYTRS', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL DKYTRS( 'U', 2, 1, A, 2, IP, B, 1, INFO )
         CALL CHKXER( 'DKYTRS', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of DERRKYX
*
      END
