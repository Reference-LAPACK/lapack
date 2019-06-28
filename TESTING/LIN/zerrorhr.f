*> \brief \b ZERRORHR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZERRORHR( PATH, NUNIT )
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
*> ZERRORHR tests the error exits for ZORHR that does Householder
*> reconstruction from the ouput of tall-skinny factorization ZLATSQR.
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
*> \ingroup complex16_lin
*
*  =====================================================================
      SUBROUTINE ZERRORHR( PATH, NUNIT )
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
      COMPLEX*16         A( NMAX, NMAX ), T1( NMAX, NMAX ),
     $                   T2( NMAX, NMAX ), W( NMAX ), D(NMAX),
     $                   C( NMAX, NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, ZORHR
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
      INTRINSIC          DBLE, DCMPLX
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
            A( I, J ) = DCMPLX( 1.D+0 / DBLE( I+J ) )
            C( I, J ) = DCMPLX( 1.D+0 / DBLE( I+J ) )
            T1( I, J ) = DCMPLX( 1.D+0 / DBLE( I+J ) )
            T2( I, J ) = DCMPLX( 1.D+0 / DBLE( I+J ) )
         END DO
         W( J ) = ( 0.D+0, 0.D+0 )
         D( J ) = ( 0.D+0, 0.D+0 )
      END DO
      OK = .TRUE.
*
*     Error exits for Householder reconstruction
*
*     ZORHR
*
      SRNAMT = 'ZORHR'
*
      INFOT = 1
      CALL ZORHR( -1, 0, 2, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 2
      CALL ZORHR( 0, -1, 3, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 1, 2, 3, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 3
      CALL ZORHR( 2, 2, 2, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 2, 2, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 1, 1, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 4
      CALL ZORHR( 0, 0, 1, 0, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 0, 0, 1, -1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 6
      CALL ZORHR( 0, 0, 1, 1, A, -1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 0, 0, 1, 1, A, 0, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 2, 2, 3, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 8
      CALL ZORHR( 0, 0, 1, 1, A, 1, T1, -1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 0, 0, 1, 1, A, 1, T1, 0, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 2, 2, 3, 2, A, 2, T1, 1, 1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 9
      CALL ZORHR( 0, 0, 1, 1, A, 1, T1, 1, -1, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 0, 0, 1, 1, A, 1, T1, 1, 0, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 11
      CALL ZORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, -1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 0, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 2, 2, 3, 1, A, 2, T1, 1, 2, T2, 1, D,
     $            W, 10, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 14
      CALL ZORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, -2, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 0, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      CALL ZORHR( 0, 0, 1, 1, A, 1, T1, 1, 1, T2, 1, D,
     $            W, 1, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
      LWMIN = 8
      CALL ZORHR( 2, 2, 3, 2, A, 2, T1, 2, 1, T2, 1, D,
     $            W, LWMIN-1, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )

      LWMIN = 20
      CALL ZORHR( 8, 2, 4, 3, A, 8, T1, 3, 2, T2, 2, D,
     $            W, LWMIN-1, INFO )
      CALL CHKXER( 'ZORHR', INFOT, NOUT, LERR, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of ZERRORHR
*
      END
