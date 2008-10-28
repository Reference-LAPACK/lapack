      SUBROUTINE SPRTLS( ISUB, SUBNAM, NDATA, NM, MVAL, NN, NVAL,
     $                   NNS, NSVAL, NNB, NBVAL, NXVAL, NLDA, LDAVAL, 
     $                   MTYPE, RSLTS, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*(*)       SUBNAM
      INTEGER            ISUB, MTYPE, NDATA, NLDA, NM, NN, NNB,
     $                   NNS, NOUT
*     ..
*     .. Array Arguments ..
      INTEGER            LDAVAL( * ), MVAL( * ), NBVAL( * ),
     $                   NSVAL( * ), NVAL( * ), NXVAL( * )
      REAL               RSLTS( 6, 6, * ) 
*     ..
*
*  Purpose
*  =======
*
*  SPRTLS prints a table of timing data for the least squares routines.
*
*  Arguments
*  =========
*
*  ISUB    (input) INTEGER
*          Subroutine index.
*
*  SUBNAM  (input) CHARACTER*(*)
*          Subroutine name. 
*
*  NDATA   (input) INTEGER
*          Number of components for subroutine SUBNAM.
*
*  NM      (input) INTEGER
*          The number of values of M contained in the vector MVAL.
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix row dimension M.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension N.
*
*  NNS     (input) INTEGER
*          The number of values of NRHS contained in the vector NSVAL.
*
*  NSVAL   (input) INTEGER array, dimension (NNS)
*          The values of the number of right hand sides NRHS.
*
*  NNB     (input) INTEGER
*          The number of values of NB and NX contained in the
*          vectors NBVAL and NXVAL.  The blocking parameters are used
*          in pairs (NB,NX).
*
*  NBVAL   (input) INTEGER array, dimension (NNB)
*          The values of the blocksize NB.
*
*  NXVAL   (input) INTEGER array, dimension (NNB)
*          The values of the crossover point NX.
*
*  NLDA    (input) INTEGER
*          The number of values of LDA contained in the vector LDAVAL.
*
*  LDAVAL  (input) INTEGER array, dimension (NLDA)
*          The values of the leading dimension of the array A.
*
*  MTYPE   (input) INTEGER
*          Number of matrix types.
*
*  RSLTS   (workspace) REAL array
*          dimension( 6, 6, number of runs )
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            ICASE, IDATA, ILDA, IM, IN, INB, INS,
     $                   ITYPE, LDA, M, N, NB, NRHS, NX
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
*     ..
*     .. Executable Statements ..
*
      ICASE = 1
*
      DO 70 IM = 1, NM
         M = MVAL( IM )
         DO 60 IN = 1, NN
            N = NVAL( IN )
            DO 50 INS = 1, NNS
               NRHS = NSVAL( INS )
               DO 40 ILDA = 1, NLDA
                  LDA = MAX( 1, LDAVAL( ILDA ) )
                  IF( ISUB.EQ.2 ) THEN
                     WRITE( NOUT, FMT = 9999 ) M, N, NRHS, LDA
                     WRITE( NOUT, FMT = 9998 )
     $                    SUBNAM(1:ILA_LEN_TRIM( SUBNAM )), ( IDATA,
     $                           IDATA = 1, NDATA-1 )
                     DO 10 ITYPE = 1, MTYPE
                        WRITE( NOUT, FMT = 9997 ) ITYPE,
     $                       ( RSLTS( IDATA, ITYPE, ICASE ),
     $                       IDATA = 1, NDATA )
   10                CONTINUE
                     ICASE = ICASE + 1
                  ELSE
                     DO 30 INB = 1, NNB
                        NB = NBVAL( INB )
                        NX = NXVAL( INB )
                        WRITE( NOUT, FMT = 9996 ) M, N, NRHS, LDA,
     $                       NB, NX               
                        WRITE( NOUT, FMT = 9998 )
     $                       SUBNAM(1:ILA_LEN_TRIM( SUBNAM )), ( IDATA,
     $                              IDATA = 1, NDATA-1 )
                        DO 20 ITYPE = 1, MTYPE
                           WRITE( NOUT, FMT = 9997 ) ITYPE,
     $                          ( RSLTS( IDATA, ITYPE, ICASE ),
     $                          IDATA = 1, NDATA )
   20                   CONTINUE
                        ICASE = ICASE + 1
   30                CONTINUE
                  END IF
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
*   
 9999 FORMAT( / ' M = ', I5, ', N = ', I5, ', NRHS = ', I5,
     $        ', LDA = ', I5 )
 9998 FORMAT( / ' TYPE ', 4X, A, 1X, 8( 4X, 'comp.', I2, : ) )
 9997 FORMAT( I5, 2X, 1P, 6G11.2 )
 9996 FORMAT( / ' M = ', I5, ', N = ', I5, ', NRHS = ', I5,
     $        ', LDA = ', I5, ', NB = ', I3, ', NX = ', I3 )
      RETURN
*
*     End of SPRTLS
*
      END
