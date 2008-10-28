      SUBROUTINE SPRTBS( LAB1, LAB2, NTYPES, DOTYPE, NSIZES, NN, NPARMS,
     $                   DOLINE, RESLTS, LDR1, LDR2, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    LAB1, LAB2
      INTEGER            LDR1, LDR2, NOUT, NPARMS, NSIZES, NTYPES
*     ..
*     .. Array Arguments ..
      LOGICAL            DOLINE( NPARMS ), DOTYPE( NTYPES )
      INTEGER            NN( NSIZES )
      REAL               RESLTS( LDR1, LDR2, * )
*     ..
*
*  Purpose
*  =======
*
*     SPRTBS prints a table of timing data for the timing programs.
*     The table has NTYPES block rows and NSIZES columns, with NPARMS
*     individual rows in each block row.
*
*  Arguments (none are modified)
*  =========
*
*  LAB1   - CHARACTER*(*)
*           The label for the rows.
*
*  LAB2   - CHARACTER*(*)
*           The label for the columns.
*
*  NTYPES - INTEGER
*           The number of values of DOTYPE, and also the
*           number of sets of rows of the table.
*
*  DOTYPE - LOGICAL array of dimension( NTYPES )
*           If DOTYPE(j) is .TRUE., then block row j (which includes
*           data from RESLTS( i, j, k ), for all i and k) will be
*           printed.  If DOTYPE(j) is .FALSE., then block row j will
*           not be printed.
*
*  NSIZES - INTEGER
*           The number of values of NN, and also the
*           number of columns of the table.
*
*  NN   -   INTEGER array of dimension( NSIZES )
*           The values of N used to label each column.
*
*  NPARMS - INTEGER
*           The number of values of LDA, hence the
*           number of rows for each value of DOTYPE.
*
*  DOLINE - LOGICAL array of dimension( NPARMS )
*           If DOLINE(i) is .TRUE., then row i (which includes data
*           from RESLTS( i, j, k ) for all j and k) will be printed.
*           If DOLINE(i) is .FALSE., then row i will not be printed.
*
*  RESLTS - REAL array of dimension( LDR1, LDR2, NSIZES )
*           The timing results.  The first index indicates the row,
*           the second index indicates the block row, and the last
*           indicates the column.
*
*  LDR1   - INTEGER
*           The first dimension of RESLTS.  It must be at least
*           min( 1, NPARMS ).
*
*  LDR2   - INTEGER
*           The second dimension of RESLTS.  It must be at least
*           min( 1, NTYPES ).
*
*  NOUT   - INTEGER
*           The output unit number on which the table
*           is to be printed.  If NOUT <= 0, no output is printed.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, ILINE, J, K
*     ..
*     .. Executable Statements ..
*
      IF( NOUT.LE.0 )
     $   RETURN
      IF( NPARMS.LE.0 )
     $   RETURN
      WRITE( NOUT, FMT = 9999 )LAB2, ( NN( I ), I = 1, NSIZES )
      WRITE( NOUT, FMT = 9998 )LAB1
*
      DO 20 J = 1, NTYPES
         ILINE = 0
         IF( DOTYPE( J ) ) THEN
            DO 10 I = 1, NPARMS
               IF( DOLINE( I ) ) THEN
                  ILINE = ILINE + 1
                  IF( ILINE.LE.1 ) THEN
                     WRITE( NOUT, FMT = 9997 )J,
     $                  ( RESLTS( I, J, K ), K = 1, NSIZES )
                  ELSE
                     WRITE( NOUT, FMT = 9996 )( RESLTS( I, J, K ),
     $                  K = 1, NSIZES )
                  END IF
               END IF
   10       CONTINUE
            IF( ILINE.GT.1 .AND. J.LT.NTYPES )
     $         WRITE( NOUT, FMT = * )
         END IF
   20 CONTINUE
      RETURN
*
 9999 FORMAT( 6X, A4, I6, 11I9 )
 9998 FORMAT( 3X, A4 )
 9997 FORMAT( 3X, I4, 4X, 1P, 12( 1X, G8.2 ) )
 9996 FORMAT( 11X, 1P, 12( 1X, G8.2 ) )
*
*     End of SPRTBS
*
      END
