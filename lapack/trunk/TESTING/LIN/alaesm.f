      SUBROUTINE ALAESM( PATH, OK, NOUT )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            OK
      CHARACTER*3        PATH
      INTEGER            NOUT
*     ..
*
*  Purpose
*  =======
*
*  ALAESM prints a summary of results from one of the -ERR- routines.
*
*  Arguments
*  =========
*
*  PATH    (input) CHARACTER*3
*          The LAPACK path name.
*
*  OK      (input) LOGICAL
*          The flag from CHKXER that indicates whether or not the tests
*          of error exits passed.
*
*  NOUT    (input) INTEGER
*          The unit number on which results are to be printed.
*          NOUT >= 0.
*
*  =====================================================================
*
*     .. Executable Statements ..
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits'
     $       )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )
      RETURN
*
*     End of ALAESM
*
      END
