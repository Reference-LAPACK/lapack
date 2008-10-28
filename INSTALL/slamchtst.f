      PROGRAM TEST2
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Local Scalars ..
      REAL               BASE, EMAX, EMIN, EPS, RMAX, RMIN, RND, SFMIN,
     $                   T, PREC
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Executable Statements ..
*
      EPS   = SLAMCH( 'Epsilon' )
      SFMIN = SLAMCH( 'Safe minimum' )
      BASE  = SLAMCH( 'Base' )
      PREC  = SLAMCH( 'Precision' )
      T     = SLAMCH( 'Number of digits in mantissa' )
      RND   = SLAMCH( 'Rounding mode' )
      EMIN  = SLAMCH( 'Minimum exponent' )
      RMIN  = SLAMCH( 'Underflow threshold' )
      EMAX  = SLAMCH( 'Largest exponent' )
      RMAX  = SLAMCH( 'Overflow threshold' )
*
      WRITE( 6, * )' Epsilon                      = ', EPS
      WRITE( 6, * )' Safe minimum                 = ', SFMIN
      WRITE( 6, * )' Base                         = ', BASE
      WRITE( 6, * )' Precision                    = ', PREC
      WRITE( 6, * )' Number of digits in mantissa = ', T
      WRITE( 6, * )' Rounding mode                = ', RND
      WRITE( 6, * )' Minimum exponent             = ', EMIN
      WRITE( 6, * )' Underflow threshold          = ', RMIN
      WRITE( 6, * )' Largest exponent             = ', EMAX
      WRITE( 6, * )' Overflow threshold           = ', RMAX
      WRITE( 6, * )' Reciprocal of safe minimum   = ', 1 / SFMIN
*
      END
