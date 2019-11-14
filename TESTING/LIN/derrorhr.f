*> \brief \b DERRORHR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DERRORHR( PATH, NUNIT )
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
*> DERRORHR tests the error exits for DORHR that does Householder
*> reconstruction from the ouput of tall-skinny factorization DLATSQR.
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
*> \date November 2019
*
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE DERRORHR( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine (version 3.9.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2019
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
      INTEGER            I, INFO, J
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   A( NMAX, NMAX ), T( NMAX, NMAX ),
     $                   R( NMAX, NMAX ), D(NMAX)
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DORHR
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
      INTRINSIC          DBLE
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
            A( I, J ) = 1.D0 / DBLE( I+J )
            R( I, J ) = 1.D0 / DBLE( I+J )
            T( I, J ) = 1.D0 / DBLE( I+J )
         END DO
         D( J ) = 0.D0
      END DO
      OK = .TRUE.
*
*     Error exits for Householder reconstruction
*
*     DORHR
*
      SRNAMT = 'DORHR'
*
      INFOT = 1
      CALL DORHR( -1, 0, A, 1, 1, T, 1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 2
      CALL DORHR( 0, -1, A, 1, 1, T, 1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
      CALL DORHR( 1, 2, A, 1, 1, T, 1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 4
      CALL DORHR( 0, 0, A, -1, 1, T, 1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
      CALL DORHR( 0, 0, A, 0, 1, T, 1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
      CALL DORHR( 2, 0, A, 1, 1, T, 1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 5
      CALL DORHR( 0, 0, A, 1, -1, T, 1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
      CALL DORHR( 0, 0, A, 1, 0, T, 1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
      INFOT = 7
      CALL DORHR( 0, 0, A, 1, 1, T, -1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
      CALL DORHR( 0, 0, A, 1, 1, T, 0, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
      CALL DORHR( 4, 3, A, 4, 2, T, 1, D, INFO )
      CALL CHKXER( 'DORHR', INFOT, NOUT, LERR, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of DERRORHR
*
      END
