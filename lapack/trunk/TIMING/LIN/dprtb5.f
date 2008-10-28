      SUBROUTINE DPRTB5( LAB1, LABM, LABN, NK, KVAL, NM, MVAL, NVAL,
     $                   NLDA, RESLTS, LDR1, LDR2, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    LAB1, LABM, LABN
      INTEGER            LDR1, LDR2, NK, NLDA, NM, NOUT
*     ..
*     .. Array Arguments ..
      INTEGER            KVAL( NK ), MVAL( NM ), NVAL( NM )
      DOUBLE PRECISION   RESLTS( LDR1, LDR2, * )
*     ..
*
*  Purpose
*  =======
*
*  DPRTB5 prints a table of timing data for the timing programs.
*  The table has NK block rows and NM columns, with NLDA
*  individual rows in each block row.  Each column depends on two
*  parameters M and N, specified as an ordered pair in the arrays MVAL
*  and NVAL.
*
*  Arguments
*  =========
*
*  LAB1    (input) CHARACTER*(*)
*          The label for the rows.
*
*  LABM    (input) CHARACTER*(*)
*          The first label for the columns.
*
*  LABN    (input) CHARACTER*(*)
*          The second label for the columns.
*
*  NK      (input) INTEGER
*          The number of values of KVAL, and also the number of block
*          rows of the table.
*
*  KVAL    (input) INTEGER array, dimension (NK)
*          The values of LAB1 used for the data in each block row.
*
*  NM      (input) INTEGER
*          The number of values of MVAL and NVAL, and also the number of
*          columns of the table.  Each column depends on the pair of
*          parameters (M,N).
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the parameter M.
*
*  NVAL    (input) INTEGER array, dimension (NM)
*          The values of the parameter N.
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
*          The second dimension of RESLTS.  LDR2 >= max(1,NM).
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
      WRITE( NOUT, FMT = 9999 )LABM, ( MVAL( I ), I = 1, NM )
      WRITE( NOUT, FMT = 9999 )LABN, ( NVAL( I ), I = 1, NM )
      WRITE( NOUT, FMT = 9998 )LAB1
*
      DO 20 I = 1, NK
         IF( LAB1.EQ.' ' ) THEN
            WRITE( NOUT, FMT = 9996 )( RESLTS( 1, J, 1 ), J = 1, NM )
         ELSE
            WRITE( NOUT, FMT = 9997 )KVAL( I ),
     $         ( RESLTS( I, J, 1 ), J = 1, NM )
         END IF
         DO 10 K = 2, NLDA
            WRITE( NOUT, FMT = 9996 )( RESLTS( I, J, K ), J = 1, NM )
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
*     End of DPRTB5
*
      END
