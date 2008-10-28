      SUBROUTINE SPRTBL( LAB1, LAB2, NK, KVAL, NN, NVAL, NLDA, RESLTS,
     $                   LDR1, LDR2, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    LAB1, LAB2
      INTEGER            LDR1, LDR2, NK, NLDA, NN, NOUT
*     ..
*     .. Array Arguments ..
      INTEGER            KVAL( NK ), NVAL( NN )
      REAL               RESLTS( LDR1, LDR2, * )
*     ..
*
*  Purpose
*  =======
*
*  SPRTBL prints a table of timing data for the timing programs.
*  The table has NK block rows and NN columns, with NLDA
*  individual rows in each block row.
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
*  NK      (input) INTEGER
*          The number of values of KVAL, and also the number of block
*          rows of the table.
*
*  KVAL    (input) INTEGER array, dimension (NK)
*          The values of LAB1 used for the data in each block row.
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
*          each value of KVAL.
*
*  RESLTS  (input) REAL array, dimension (LDR1, LDR2, NLDA)
*          The timing results for each value of N, K, and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max( 1, NK ).
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
      INTEGER            I, J, K
*     ..
*     .. Executable Statements ..
*
      IF( NOUT.LE.0 )
     $   RETURN
      WRITE( NOUT, FMT = 9999 )LAB2, ( NVAL( I ), I = 1, NN )
      WRITE( NOUT, FMT = 9998 )LAB1
*
      DO 20 I = 1, NK
         IF( LAB1.EQ.' ' ) THEN
            WRITE( NOUT, FMT = 9996 )( RESLTS( 1, J, 1 ), J = 1, NN )
         ELSE
            WRITE( NOUT, FMT = 9997 )KVAL( I ),
     $         ( RESLTS( I, J, 1 ), J = 1, NN )
         END IF
         DO 10 K = 2, NLDA
            WRITE( NOUT, FMT = 9996 )( RESLTS( I, J, K ), J = 1, NN )
   10    CONTINUE
         IF( NLDA.GT.1 )
     $      WRITE( NOUT, FMT = * )
   20 CONTINUE
      IF( NLDA.EQ.1 )
     $   WRITE( NOUT, FMT = * )
      RETURN
*
 9999 FORMAT( 6X, A4, I6, 11I8 )
 9998 FORMAT( 3X, A4 )
 9997 FORMAT( 1X, I6, 1X, 12F8.1 )
 9996 FORMAT( 8X, 12F8.1 )
*
*     End of SPRTBL
*
      END
