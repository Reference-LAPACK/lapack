      SUBROUTINE ATIMIN( PATH, LINE, NSUBS, NAMES, TIMSUB, NOUT, INFO )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      CHARACTER*( * )    PATH
      INTEGER            INFO, NOUT, NSUBS
*     ..
*     .. Array Arguments ..
      LOGICAL            TIMSUB( * )
      CHARACTER*( * )    NAMES( * )
*     ..
*
*  Purpose
*  =======
*
*  ATIMIN interprets the input line for the timing routines.
*  The LOGICAL array TIMSUB returns .true. for each routine to be
*  timed and .false. for the routines which are not to be timed.
*
*  Arguments
*  =========
*
*  PATH    (input) CHARACTER*(*)
*          The LAPACK path name of the calling routine.  The path name
*          may be at most 6 characters long.  If LINE(1:LEN(PATH)) is
*          the same as PATH, then the input line is searched for NSUBS
*          non-blank characters, otherwise, the input line is assumed to
*          specify a single subroutine name.
*
*  LINE    (input) CHARACTER*80
*          The input line to be evaluated.  The path or subroutine name
*          must begin in column 1 and the part of the line after the
*          name is used to indicate the routines to be timed.
*          See below for further details.
*
*  NSUBS   (input) INTEGER
*          The number of subroutines in the LAPACK path name of the
*          calling routine.
*
*  NAMES   (input) CHARACTER*(*) array, dimension (NSUBS)
*          The names of the subroutines in the LAPACK path name of the
*          calling routine.
*
*  TIMSUB  (output) LOGICAL array, dimension (NSUBS)
*          For each I from 1 to NSUBS, TIMSUB( I ) is set to .true. if
*          the subroutine NAMES( I ) is to be timed; otherwise,
*          TIMSUB( I ) is set to .false.
*
*  NOUT    (input) INTEGER
*          The unit number on which error messages will be printed.
*
*  INFO    (output) INTEGER
*          The return status of this routine.
*          = -1:  Unrecognized path or subroutine name
*          =  0:  Normal return
*          =  1:  Name was recognized, but no timing requested
*
*  Further Details
*  ======= =======
*
*  An input line begins with a subroutine or path name, optionally
*  followed by one or more non-blank characters indicating the specific
*  routines to be timed.
*
*  If the character string in PATH appears at the beginning of LINE,
*  up to NSUBS routines may be timed.  If LINE is blank after the path
*  name, all the routines in the path will be timed.  If LINE is not
*  blank after the path name, the rest of the line is searched
*  for NSUBS nonblank characters, and if the i-th such character is
*  't' or 'T', then the i-th subroutine in this path will be timed.
*  For example, the input line
*     SGE    T T T T
*  requests timing of the first 4 subroutines in the SGE path.
*
*  If the character string in PATH does not appear at the beginning of
*  LINE, then LINE is assumed to begin with a subroutine name.  The name
*  is assumed to end in column 6 or in column i if column i+1 is blank
*  and i+1 <= 6.  If LINE is completely blank after the subroutine name,
*  the routine will be timed.  If LINE is not blank after the subroutine
*  name, then the subroutine will be timed if the first non-blank after
*  the name is 't' or 'T'.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            REQ
      CHARACTER(32)      CNAME
      INTEGER            I, ISTART, ISTOP, ISUB, LCNAME, LNAMES, LPATH
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      LOGICAL            LSAME, LSAMEN
      EXTERNAL           LSAME, LSAMEN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN, MIN
*     ..
*     .. Executable Statements ..
*
*
*     Initialize
*
      INFO = 0
      LCNAME = 1
      DO 10 I = 2, 6
         IF( LINE( I: I ).EQ.' ' )
     $      GO TO 20
         LCNAME = I
   10 CONTINUE
   20 CONTINUE
      LPATH = MIN( LCNAME+1, LEN( PATH ) )
      LNAMES = MIN( LCNAME+1, LEN( NAMES( 1 ) ) )
      CNAME = LINE( 1: LCNAME )
*
      DO 30 I = 1, NSUBS
         TIMSUB( I ) = .FALSE.
   30 CONTINUE
      ISTOP = 0
*
*     Check for a valid path or subroutine name.
*
      IF( LCNAME.LE.LEN( PATH ) .AND. LSAMEN( LPATH, CNAME, PATH ) )
     $     THEN
         ISTART = 1
         ISTOP = NSUBS
      ELSE IF( LCNAME.LE.LEN( NAMES( 1 ) ) ) THEN
         DO 40 I = 1, NSUBS
            IF( LSAMEN( LNAMES, CNAME, NAMES( I ) ) ) THEN
               ISTART = I
               ISTOP = I
            END IF
   40    CONTINUE
      END IF
*
      IF( ISTOP.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
 9999    FORMAT( 1X, A, ':  Unrecognized path or subroutine name', / )
         INFO = -1
         GO TO 110
      END IF
*
*     Search the rest of the input line for 1 or NSUBS nonblank
*     characters, where 'T' or 't' means 'Time this routine'.
*
      ISUB = ISTART
      DO 50 I = LCNAME + 1, 80
         IF( LINE( I: I ).NE.' ' ) THEN
            TIMSUB( ISUB ) = LSAME( LINE( I: I ), 'T' )
            ISUB = ISUB + 1
            IF( ISUB.GT.ISTOP )
     $         GO TO 60
         END IF
   50 CONTINUE
   60 CONTINUE
*
*     If no characters appear after the routine or path name, then
*     time the routine or all the routines in the path.
*
      IF( ISUB.EQ.ISTART ) THEN
         DO 70 I = ISTART, ISTOP
            TIMSUB( I ) = .TRUE.
   70    CONTINUE
      ELSE
*
*        Test to see if any timing was requested.
*
         REQ = .FALSE.
         DO 80 I = ISTART, ISUB - 1
            REQ = REQ .OR. TIMSUB( I )
   80    CONTINUE
         IF( .NOT.REQ ) THEN
            WRITE( NOUT, FMT = 9998 )CNAME(1:ILA_LEN_TRIM(CNAME))
 9998       FORMAT( 1X, A, ' was not timed', / )
            INFO = 1
            GO TO 110
         END IF
   90    CONTINUE
*
*       If fewer than NSUBS characters are specified for a path name,
*       the rest are assumed to be 'F'.
*
         DO 100 I = ISUB, ISTOP
            TIMSUB( I ) = .FALSE.
  100    CONTINUE
      END IF
  110 CONTINUE
      RETURN
*
*     End of ATIMIN
*
      END
