*> \brief \b ZCHKEC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZCHKEC( THRESH, TSTERR, NIN, NOUT )
*
*       .. Scalar Arguments ..
*       LOGICAL            TSTERR
*       INTEGER            NIN, NOUT
*       DOUBLE PRECISION   THRESH
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZCHKEC tests eigen- condition estimation routines
*>        ZTRSYL, CTREXC, CTRSNA, CTRSEN
*>
*> In all cases, the routine runs through a fixed set of numerical
*> examples, subjects them to various tests, and compares the test
*> results to a threshold THRESH. In addition, ZTRSNA and CTRSEN are
*> tested by reading in precomputed examples from a file (on input unit
*> NIN).  Output is written to output unit NOUT.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] THRESH
*> \verbatim
*>          THRESH is DOUBLE PRECISION
*>          Threshold for residual tests.  A computed test ratio passes
*>          the threshold if it is less than THRESH.
*> \endverbatim
*>
*> \param[in] TSTERR
*> \verbatim
*>          TSTERR is LOGICAL
*>          Flag that indicates whether error exits are to be tested.
*> \endverbatim
*>
*> \param[in] NIN
*> \verbatim
*>          NIN is INTEGER
*>          The logical unit number for input.
*> \endverbatim
*>
*> \param[in] NOUT
*> \verbatim
*>          NOUT is INTEGER
*>          The logical unit number for output.
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
*> \ingroup complex16_eig
*
*  =====================================================================
      SUBROUTINE ZCHKEC( THRESH, TSTERR, NIN, NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NIN, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            OK
      CHARACTER*3        PATH
      INTEGER            KTREXC, KTRSEN, KTRSNA, KTRSYL, KTRSYL3,
     $                   LTREXC, LTRSYL, NTESTS, NTREXC, NTRSYL
      DOUBLE PRECISION   EPS, RTREXC, SFMIN
*     ..
*     .. Local Arrays ..
      INTEGER            FTRSYL( 3 ), ITRSYL( 2 ), LTRSEN( 3 ),
     $                   LTRSNA( 3 ), NTRSEN( 3 ), NTRSNA( 3 )
      DOUBLE PRECISION   RTRSEN( 3 ), RTRSNA( 3 ), RTRSYL( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZERREC, ZGET35, ZGET36, ZGET37, ZGET38, ZSYL01
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Executable Statements ..
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'EC'
      EPS = DLAMCH( 'P' )
      SFMIN = DLAMCH( 'S' )
      WRITE( NOUT, FMT = 9994 )
      WRITE( NOUT, FMT = 9993 )EPS, SFMIN
      WRITE( NOUT, FMT = 9992 )THRESH
*
*     Test error exits if TSTERR is .TRUE.
*
      IF( TSTERR )
     $   CALL ZERREC( PATH, NOUT )
*
      OK = .TRUE.
      CALL ZGET35( RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL, NIN )
      IF( RTRSYL( 1 ).GT.THRESH ) THEN
         OK = .FALSE.
         WRITE( NOUT, FMT = 9999 )RTRSYL( 1 ), LTRSYL, NTRSYL, KTRSYL
      END IF
*
      CALL ZSYL01( THRESH, FTRSYL, RTRSYL, ITRSYL, KTRSYL3 )
      IF( FTRSYL( 1 ).GT.0 ) THEN
         OK = .FALSE.
         WRITE( NOUT, FMT = 9970 )FTRSYL( 1 ), RTRSYL( 1 ), THRESH
      END IF
      IF( FTRSYL( 2 ).GT.0 ) THEN
         OK = .FALSE.
         WRITE( NOUT, FMT = 9971 )FTRSYL( 2 ), RTRSYL( 2 ), THRESH
      END IF
      IF( FTRSYL( 3 ).GT.0 ) THEN
         OK = .FALSE.
         WRITE( NOUT, FMT = 9972 )FTRSYL( 3 )
      END IF
*
      CALL ZGET36( RTREXC, LTREXC, NTREXC, KTREXC, NIN )
      IF( RTREXC.GT.THRESH .OR. NTREXC.GT.0 ) THEN
         OK = .FALSE.
         WRITE( NOUT, FMT = 9998 )RTREXC, LTREXC, NTREXC, KTREXC
      END IF
*
      CALL ZGET37( RTRSNA, LTRSNA, NTRSNA, KTRSNA, NIN )
      IF( RTRSNA( 1 ).GT.THRESH .OR. RTRSNA( 2 ).GT.THRESH .OR.
     $    NTRSNA( 1 ).NE.0 .OR. NTRSNA( 2 ).NE.0 .OR. NTRSNA( 3 ).NE.0 )
     $     THEN
         OK = .FALSE.
         WRITE( NOUT, FMT = 9997 )RTRSNA, LTRSNA, NTRSNA, KTRSNA
      END IF
*
      CALL ZGET38( RTRSEN, LTRSEN, NTRSEN, KTRSEN, NIN )
      IF( RTRSEN( 1 ).GT.THRESH .OR. RTRSEN( 2 ).GT.THRESH .OR.
     $    NTRSEN( 1 ).NE.0 .OR. NTRSEN( 2 ).NE.0 .OR. NTRSEN( 3 ).NE.0 )
     $     THEN
         OK = .FALSE.
         WRITE( NOUT, FMT = 9996 )RTRSEN, LTRSEN, NTRSEN, KTRSEN
      END IF
*
      NTESTS = KTRSYL + KTRSYL3 + KTREXC + KTRSNA + KTRSEN
      IF( OK )
     $   WRITE( NOUT, FMT = 9995 )PATH, NTESTS
*
 9999 FORMAT( ' Error in ZTRSYL: RMAX =', D12.3, / ' LMAX = ', I8,
     $      ' NINFO=', I8, ' KNT=', I8 )
 9998 FORMAT( ' Error in ZTREXC: RMAX =', D12.3, / ' LMAX = ', I8,
     $      ' NINFO=', I8, ' KNT=', I8 )
 9997 FORMAT( ' Error in ZTRSNA: RMAX =', 3D12.3, / ' LMAX = ', 3I8,
     $      ' NINFO=', 3I8, ' KNT=', I8 )
 9996 FORMAT( ' Error in ZTRSEN: RMAX =', 3D12.3, / ' LMAX = ', 3I8,
     $      ' NINFO=', 3I8, ' KNT=', I8 )
 9995 FORMAT( / 1X, 'All tests for ', A3,
     $      ' routines passed the threshold ( ', I6, ' tests run)' )
 9994 FORMAT( ' Tests of the Nonsymmetric eigenproblem condition',
     $      ' estimation routines', / ' ZTRSYL, ZTREXC, ZTRSNA, ZTRSEN',
     $      / )
 9993 FORMAT( ' Relative machine precision (EPS) = ', D16.6,
     $      / ' Safe minimum (SFMIN)             = ', D16.6, / )
 9992 FORMAT( ' Routines pass computational tests if test ratio is ',
     $      'less than', F8.2, / / )
 9970 FORMAT( 'Error in ZTRSYL: ', I8, ' tests fail the threshold.', /
     $      'Maximum test ratio =', D12.3, ' threshold =', D12.3 )
 9971 FORMAT( 'Error in ZTRSYL3: ', I8, ' tests fail the threshold.', /
     $      'Maximum test ratio =', D12.3, ' threshold =', D12.3 )
 9972 FORMAT( 'ZTRSYL and ZTRSYL3 compute an inconsistent scale ',
     $      'factor in ', I8, ' tests.')
      RETURN
*
*     End of ZCHKEC
*
      END
