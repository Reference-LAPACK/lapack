      SUBROUTINE CHKXER( SRNAMT, INFOT, NOUT, LERR, OK )
*
*  Tests whether XERBLA has detected an error when it should.
*
*  Auxiliary routine for test program for Level 2 Blas.
*
*  -- Written on 10-August-1987.
*     Richard Hanson, Sandia National Labs.
*     Jeremy Du Croz, NAG Central Office.
*
*  =====================================================================
*
*     .. Scalar Arguments ..
      LOGICAL            LERR, OK
      CHARACTER*(*)      SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
*     ..
*     .. Executable Statements ..
      IF( .NOT.LERR ) THEN
         WRITE( NOUT, FMT = 9999 )INFOT,
     $        SRNAMT( 1:LEN_TRIM( SRNAMT ) )
         OK = .FALSE.
      END IF
      LERR = .FALSE.
      RETURN
*
 9999 FORMAT( ' *** Illegal value of parameter number ', I2,
     $      ' not detected by ', A6, ' ***' )
*
*     End of CHKXER.
*
      END
