*> \brief \b TSTIEE
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      PROGRAM TSTIEE
*
*  -- LAPACK test routine --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Local Scalars ..
      INTEGER            IEEEOK
*     ..
*     .. Executable Statements ..
*
      WRITE( 6, FMT = * )
     $   'We are about to check whether infinity arithmetic'
      WRITE( 6, FMT = * )'can be trusted.  If this test hangs, set'
      WRITE( 6, FMT = * )
     $   'ILAENV = 0 for ISPEC = 11 in LAPACK/SRC/ilaenv.f'
*
      IEEEOK = ILAENV( 11, 'ILAENV', 'N', 1, 2, 3, 4 )
      WRITE( 6, FMT = * )
*
      IF( IEEEOK.EQ.0 ) THEN
         WRITE( 6, FMT = * )
     $      'Infinity arithmetic did not perform per the ieee spec'
      ELSE
         WRITE( 6, FMT = * )
     $      'Infinity arithmetic performed as per the ieee spec.'
         WRITE( 6, FMT = * )
     $      'However, this is not an exhaustive test and does not'
         WRITE( 6, FMT = * )
     $      'guarantee that infinity arithmetic meets the',
     $      ' ieee spec.'
      END IF
*
      WRITE( 6, FMT = * )
*     ILAENV( 10, ...) checks both infinity and NaN arithmetic
*     infinity has already been checked so checking NaN now
      WRITE( 6, FMT = * )
     $   'We are about to check whether NaN arithmetic'
      WRITE( 6, FMT = * )'can be trusted.  If this test hangs, set'
      WRITE( 6, FMT = * )
     $   'ILAENV = 0 for ISPEC = 10 in LAPACK/SRC/ilaenv.f'
      IEEEOK = ILAENV( 10, 'ILAENV', 'N', 1, 2, 3, 4 )
*
      WRITE( 6, FMT = * )
      IF( IEEEOK.EQ.0 ) THEN
         WRITE( 6, FMT = * )
     $      'NaN arithmetic did not perform per the ieee spec'
      ELSE
         WRITE( 6, FMT = * )'NaN arithmetic performed as per the ieee',
     $      ' spec.'
         WRITE( 6, FMT = * )
     $      'However, this is not an exhaustive test and does not'
         WRITE( 6, FMT = * )'guarantee that NaN arithmetic meets the',
     $      ' ieee spec.'
      END IF
      WRITE( 6, FMT = * )
*
      END
