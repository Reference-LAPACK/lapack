      SUBROUTINE SPRTB2( LAB1, LAB2, LAB3, NN, NVAL, NLDA, RESLTS, LDR1,
     $                   LDR2, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    LAB1, LAB2, LAB3
      INTEGER            LDR1, LDR2, NLDA, NN, NOUT
*     ..
*     .. Array Arguments ..
      INTEGER            NVAL( NN )
      REAL               RESLTS( LDR1, LDR2, * )
*     ..
*
*  Purpose
*  =======
*
*  SPRTB2 prints a table of timing data for the solve routines.
*  There are 4 rows to each table, corresponding to
*  NRHS = 1, 2, N/2, and N,  or  NRHS = 1, 2, K/2, K for the
*  band routines.
*
*  Arguments
*  =========
*
*  LAB1    (input) CHARACTER*(*)
*          The label for the rows.
*
*  LAB2    (input) CHARACTER*(*)
*          The label for the columns.
*
*  LAB3    CHARACTER*(*)
*          The name of the variable used in the row headers (usually
*          N or K).
*
*  NN      (input) INTEGER
*          The number of values of NVAL, and also the number of columns
*          of the table.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of LAB2 used for the data in each column.
*
*  NLDA    (input) INTEGER
*          The number of values of LDA, hence the number of rows for
*          each value of NRHS.
*
*  RESLTS  (input) REAL array, dimension (LDR1, LDR2, NLDA)
*          The timing results for each value of N, K, and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= 4.
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max( 1, NN ).
*
*  NOUT    (input) INTEGER
*          The unit number on which the table is to be printed.
*          NOUT >= 0.
*
*  =====================================================================
*
*     .. Local Scalars ..
      CHARACTER(32)      COLLAB
      INTEGER            I, IC, INB, J, K, LNB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN, MAX
*     ..
*     .. Executable Statements ..
*
      IF( NOUT.LE.0 )
     $   RETURN
      WRITE( NOUT, FMT = 9999 )LAB2, ( NVAL( I ), I = 1, NN )
      WRITE( NOUT, FMT = 9998 )LAB1
*
*     Find the first and last non-blank characters in LAB3.
*
      INB = 0
      DO 10 I = 1, LEN( LAB3 )
         IF( INB.EQ.0 .AND. LAB3( I: I ).NE.' ' )
     $      INB = I
         IF( LAB3( I: I ).NE.' ' )
     $      LNB = I
   10 CONTINUE
      IF( INB.EQ.0 ) THEN
         INB = 1
         LNB = 1
      END IF
*
      DO 50 I = 1, 4
         IF( I.EQ.1 ) THEN
            COLLAB = '     1'
         ELSE IF( I.EQ.2 ) THEN
            COLLAB = '     2'
         ELSE IF( I.EQ.3 ) THEN
            COLLAB = '    /2'
            DO 20 J = LNB, MAX( INB, LNB-3 ), -1
               IC = 4 - ( LNB-J )
               COLLAB( IC: IC ) = LAB3( J: J )
   20       CONTINUE
         ELSE IF( I.EQ.4 ) THEN
            COLLAB = ' '
            DO 30 J = LNB, MAX( INB, LNB-5 ), -1
               IC = 6 - ( LNB-J )
               COLLAB( IC: IC ) = LAB3( J: J )
   30       CONTINUE
         END IF
         WRITE( NOUT, FMT = 9997 )COLLAB,
     $      ( RESLTS( I, J, 1 ), J = 1, NN )
         DO 40 K = 2, NLDA
            WRITE( NOUT, FMT = 9996 )( RESLTS( I, J, K ), J = 1, NN )
   40    CONTINUE
         IF( NLDA.GT.1 )
     $      WRITE( NOUT, FMT = * )
   50 CONTINUE
      IF( NLDA.EQ.1 )
     $   WRITE( NOUT, FMT = * )
*
 9999 FORMAT( 6X, A4, I6, 11I8 )
 9998 FORMAT( 3X, A4 )
 9997 FORMAT( 1X, A, 1X, 12F8.1 )
 9996 FORMAT( 8X, 12F8.1 )
*
      RETURN
*
*     End of SPRTB2
*
      END
