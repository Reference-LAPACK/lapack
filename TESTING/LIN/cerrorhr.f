*> \brief \b CERRORHR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CERRORHR( PATH, NUNIT )
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
*> CERRORHR tests the error exits for high-level CORHR and
*> low-level CLAORHR that does Householder reconstruction from
*> the ouput of tall-skinny factorization CLATSQR.
*>
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
*> \date June 2019
*
*> \ingroup complex_lin
*
*  =====================================================================
      SUBROUTINE CERRORHR( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine (version 3.9.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2019
*
*     .. Scalar Arguments ..
      CHARACTER(LEN=3)   PATH
      INTEGER            NUNIT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J, LWMIN
*     ..
*     .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX ), T1( NMAX, NMAX ),
     $                   T2( NMAX, NMAX ), W( NMAX ), D(NMAX),
     $                   W2(NMAX, NMAX), C( NMAX, NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, CLAORHR, CORHR
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER(LEN=32)  SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, CMPLX
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
*     Set the variables to innocuous values.
*
      DO J = 1, NMAX
         DO I = 1, NMAX
            A( I, J ) = CMPLX( 1.E+0 / REAL( I+J ) )
            C( I, J ) = CMPLX( 1.E+0 / REAL( I+J ) )
            T1( I, J ) = CMPLX( 1.E+0 / REAL( I+J ) )
            T2( I, J ) = CMPLX( 1.E+0 / REAL( I+J ) )
            W2( I, J ) = CMPLX( 1.E+0 / DBLE( I+J ) )
         END DO
         W( J ) = ( 0.E+0, 0.E+0 )
         D( J ) = ( 0.E+0, 0.E+0 )
      END DO
      OK = .TRUE.
*
*     Error exits for Householder reconstruction
*
*     CORHR
*
      SRNAMT = 'CORHR'
*
      INFOT = 1
      CALL CORHR( -1, 0, 2, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 2
      CALL CORHR( 0, -1, 3, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 1, 2, 3, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 3
      CALL CORHR( 2, 2, 2, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 2, 2, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 1, 1, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 4
      CALL CORHR( 0, 0, 1, 0, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 0, 0, 1, -1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 6
      CALL CORHR( 0, 0, 1, 1, A, -1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 0, 0, 1, 1, A, 0, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 2, 2, 3, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 8
      CALL CORHR( 0, 0, 1, 1, A, 1, T1, -1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 0, 0, 1, 1, A, 1, T1, 0, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 2, 2, 3, 2, A, 2, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 9
      CALL CORHR( 0, 0, 1, 1, A, 1, T1, 1, -1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 0, 0, 1, 1, A, 1, T1, 1, 0, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 11
      CALL CORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, -1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 0, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 2, 2, 3, 1, A, 2, T1, 1, 2, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 14
      CALL CORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, -2, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 0, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      CALL CORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 1, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
      LWMIN = 8
      CALL CORHR( 2, 2, 3, 2, A, 2, T1, 2, 1, T2, 1, D,
     $            W, LWMIN-1, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )

      LWMIN = 20
      CALL CORHR( 8, 2, 4, 3, A, 8, T1, 3, 2, T2, 2, D,
     $            W, LWMIN-1, INFO )
      CALL CHKXER( 'CORHR', INFOT, NOUT, LERR, OK )
*
*     CLAORHR
*
      SRNAMT = 'CLAORHR'
*
      INFOT = 1
      CALL CLAORHR( -1, 0, 2, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 2
      CALL CLAORHR( 0, -1, 3, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 1, 2, 3, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 3
      CALL CLAORHR( 2, 2, 2, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 2, 2, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 1, 1, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 4
      CALL CLAORHR( 0, 0, 1, 0, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 0, 0, 1, -1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 6
      CALL CLAORHR( 0, 0, 1, 1, A, -1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 0, 0, 1, 1, A, 0, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 2, 2, 3, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 8
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, -1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, 0, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 2, 2, 3, 2, A, 2, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 9
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, 1, -1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, 1, 0, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 11
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, -1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 0, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 2, 2, 3, 1, A, 2, T1, 1, 2, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 14
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, -1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 0, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 2, 2, 3, 1, A, 2, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 10, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 16
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, -2, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      CALL CLAORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $              W2, 1, W, 0, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      LWMIN = 4
      CALL CLAORHR( 2, 2, 3, 2, A, 2, T1, 2, 1, T2, 1, D,
     $              W2, 2, W, LWMIN-1, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
      LWMIN = 4
      CALL CLAORHR( 8, 2, 4, 3, A, 8, T1, 3, 2, T2, 2, D,
     $              W2, 8, W, LWMIN-1, INFO )
      CALL CHKXER( 'CLAORHR', INFOT, NOUT, LERR, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of CERRORHR
*
      END
