      PROGRAM TEST5
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Parameters ..
      INTEGER            NMAX, ITS
      PARAMETER          ( NMAX = 100, ITS = 5000 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ALPHA, AVG, T1, T2, TNOSEC
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   X( NMAX ), Y( NMAX )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DSECND
      EXTERNAL           DSECND
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*
*     Initialize X and Y
*
      DO 10 I = 1, NMAX
         X( I ) = DBLE( 1 ) / DBLE( I )
         Y( I ) = DBLE( NMAX-I ) / DBLE( NMAX )
   10 CONTINUE
      ALPHA = 0.315D0
*
*     Time 1,000,000 DAXPY operations
*
      T1 = DSECND( )
      DO 30 J = 1, ITS
         DO 20 I = 1, NMAX
            Y( I ) = Y( I ) + ALPHA*X( I )
   20    CONTINUE
         ALPHA = -ALPHA
   30 CONTINUE
      T2 = DSECND( )
      WRITE( 6, 9999 )T2 - T1
      IF( T2-T1.GT.0.0D0 ) THEN
         WRITE( 6, 9998 )1.0D0 / ( T2-T1 )
      ELSE
         WRITE( 6, 9994 )
      END IF
      TNOSEC = T2 - T1
*
*     Time 1,000,000 DAXPY operations with DSECND in the outer loop
*
      T1 = DSECND( )
      DO 50 J = 1, ITS
         DO 40 I = 1, NMAX
            Y( I ) = Y( I ) + ALPHA*X( I )
   40    CONTINUE
         ALPHA = -ALPHA
         T2 = DSECND( )
   50 CONTINUE
*
*     Compute the time in milliseconds used by an average call
*     to DSECND.
*
      WRITE( 6, 9997 )T2 - T1
      AVG = ( ( T2-T1 )-TNOSEC )*1000.D0 / DBLE( ITS )
      WRITE( 6, 9996 )AVG
*
*     Compute the equivalent number of floating point operations used
*     by an average call to DSECND.
*
      IF( TNOSEC.GT.0.0D0 )
     $   WRITE( 6, 9995 )1000.D0*AVG / TNOSEC
*
 9999 FORMAT( ' Time for 1,000,000 DAXPY ops  = ', G10.3, ' seconds' )
 9998 FORMAT( ' DAXPY performance rate        = ', G10.3, ' mflops ' )
 9997 FORMAT( ' Including DSECND, time        = ', G10.3, ' seconds' )
 9996 FORMAT( ' Average time for DSECND       = ', G10.3,
     $      ' milliseconds' )
 9995 FORMAT( ' Equivalent floating point ops = ', G10.3, ' ops' )
 9994 FORMAT( ' *** Error:  Time for operations was zero' )
      CALL MYSUB(NMAX,X,Y)
      END
      SUBROUTINE MYSUB(N,X,Y)
      INTEGER N
      DOUBLE PRECISION X(N), Y(N)
      RETURN
      END
