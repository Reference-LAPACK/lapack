      SUBROUTINE DPRTB3( LAB1, LAB2, NK, KVAL, LVAL, NN, NVAL, NLDA,
     $                   RESLTS, LDR1, LDR2, NOUT )
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
      INTEGER            KVAL( NK ), LVAL( NK ), NVAL( NN )
      DOUBLE PRECISION   RESLTS( LDR1, LDR2, * )
*     ..
*
*  Purpose
*  =======
*
*  DPRTB3 prints a table of timing data for the timing programs.
*  The table has NK block rows and NN columns, with NLDA
*  individual rows in each block row.  Each block row depends on two
*  parameters K and L, specified as an ordered pair in the arrays KVAL
*  and LVAL.
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
*          The values of the parameter K.  Each block row depends on
*          the pair of parameters (K, L).
*
*  LVAL    (input) INTEGER array, dimension (NK)
*          The values of the parameter L.  Each block row depends on
*          the pair of parameters (K, L).
*
*  NN      (input) INTEGER
*          The number of values of NVAL, and also the number of columns
*          of the table.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of N used for the data in each column.
*
*  NLDA    (input) INTEGER
*          The number of values of LDA, hence the number of rows for
*          each value of KVAL.
*
*  RESLTS  (input) DOUBLE PRECISION array, dimension (LDR1, LDR2, NLDA)
*          The timing results for each value of N, K, and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(1,NK).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NN).
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
            WRITE( NOUT, FMT = 9997 )KVAL( I ), LVAL( I ),
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
 9999 FORMAT( 10X, A4, I7, 11I8 )
 9998 FORMAT( 1X, A11 )
 9997 FORMAT( 1X, '(', I4, ',', I4, ') ', 12F8.1 )
 9996 FORMAT( 13X, 12F8.1 )
*
*     End of DPRTB3
*
      END
