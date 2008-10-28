      SUBROUTINE ALASUM( TYPE, NOUT, NFAIL, NRUN, NERRS )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*3        TYPE
      INTEGER            NFAIL, NOUT, NRUN, NERRS
*     ..
*
*  Purpose
*  =======
*
*  ALASUM prints a summary of results from one of the -CHK- routines.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*3
*          The LAPACK path name.
*
*  NOUT    (input) INTEGER
*          The unit number on which results are to be printed.
*          NOUT >= 0.
*
*  NFAIL   (input) INTEGER
*          The number of tests which did not pass the threshold ratio.
*
*  NRUN    (input) INTEGER
*          The total number of tests.
*
*  NERRS   (input) INTEGER
*          The number of error messages recorded.
*
*  =====================================================================
*
*     .. Executable Statements ..
*
      IF( NFAIL.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )TYPE, NFAIL, NRUN
      ELSE
         WRITE( NOUT, FMT = 9998 )TYPE, NRUN
      END IF
      IF( NERRS.GT.0 ) THEN
         WRITE( NOUT, FMT = 9997 )NERRS
      END IF
*
 9999 FORMAT( 1X, A3, ': ', I6, ' out of ', I6,
     $      ' tests failed to pass the threshold' )
 9998 FORMAT( /1X, 'All tests for ', A3,
     $      ' routines passed the threshold (', I6, ' tests run)' )
 9997 FORMAT( 6X, I6, ' error messages recorded' )
      RETURN
*
*     End of ALASUM
*
      END
