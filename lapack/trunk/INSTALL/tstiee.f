C> \brief \b TSTIEE
C>\details
C> \b Purpose:
C>\verbatim
C>
C> TEST IEEE
C>
C>\endverbatim
C> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
C> \date November 2011
C> \ingroup auxOTHERauxiliary

      PROGRAM TSTIEE
C
C  -- LAPACK test routine (version 3.2) --
C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
C     November 2006
C
C     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
C     ..
C     .. Local Scalars ..
      INTEGER            IEEEOK
C     ..
C     .. Executable Statements ..
C
      WRITE( 6, FMT = * )
     $   'We are about to check whether infinity arithmetic'
      WRITE( 6, FMT = * )'can be trusted.  If this test hangs, set'
      WRITE( 6, FMT = * )
     $   'ILAENV = 0 for ISPEC = 10 in LAPACK/SRC/ilaenv.f'
C
      IEEEOK = ILAENV( 10, 'ILAENV', 'N', 1, 2, 3, 4 )
      WRITE( 6, FMT = * )
C
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
C
      WRITE( 6, FMT = * )
      WRITE( 6, FMT = * )
     $   'We are about to check whether NaN arithmetic'
      WRITE( 6, FMT = * )'can be trusted.  If this test hangs, set'
      WRITE( 6, FMT = * )
     $   'ILAENV = 0 for ISPEC = 11 in LAPACK/SRC/ilaenv.f'
      IEEEOK = ILAENV( 11, 'ILAENV', 'N', 1, 2, 3, 4 )
C
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
C
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
C
C  -- LAPACK auxiliary routine (version 3.2) --
C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
C     November 2006
C
C     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
C     ..
C
C  Purpose
C  =======
C
C  ILAENV is called from the LAPACK routines to choose problem-dependent
C  parameters for the local environment.  See ISPEC for a description of
C  the parameters.
C
C  This version provides a set of parameters which should give good,
C  but not optimal, performance on many of the currently available
C  computers.  Users are encouraged to modify this subroutine to set
C  the tuning parameters for their particular machine using the option
C  and problem size information in the arguments.
C
C  This routine will not function correctly if it is converted to all
C  lower case.  Converting it to all upper case is allowed.
C
C  Arguments
C  =========
C
C  ISPEC   (input) INTEGER
C          Specifies the parameter to be returned as the value of
C          ILAENV.
C          = 1: the optimal blocksize; if this value is 1, an unblocked
C               algorithm will give the best performance.
C          = 2: the minimum block size for which the block routine
C               should be used; if the usable block size is less than
C               this value, an unblocked routine should be used.
C          = 3: the crossover point (in a block routine, for N less
C               than this value, an unblocked routine should be used)
C          = 4: the number of shifts, used in the nonsymmetric
C               eigenvalue routines
C          = 5: the minimum column dimension for blocking to be used;
C               rectangular blocks must have dimension at least k by m,
C               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
C          = 6: the crossover point for the SVD (when reducing an m by n
C               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
C               this value, a QR factorization is used first to reduce
C               the matrix to a triangular form.)
C          = 7: the number of processors
C          = 8: the crossover point for the multishift QR and QZ methods
C               for nonsymmetric eigenvalue problems.
C          = 9: maximum size of the subproblems at the bottom of the
C               computation tree in the divide-and-conquer algorithm
C               (used by xGELSD and xGESDD)
C          =10: ieee NaN arithmetic can be trusted not to trap
C          =11: infinity arithmetic can be trusted not to trap
C
C  NAME    (input) CHARACTER*(*)
C          The name of the calling subroutine, in either upper case or
C          lower case.
C
C  OPTS    (input) CHARACTER*(*)
C          The character options to the subroutine NAME, concatenated
C          into a single character string.  For example, UPLO = 'U',
C          TRANS = 'T', and DIAG = 'N' for a triangular routine would
C          be specified as OPTS = 'UTN'.
C
C  N1      (input) INTEGER
C  N2      (input) INTEGER
C  N3      (input) INTEGER
C  N4      (input) INTEGER
C          Problem dimensions for the subroutine NAME; these may not all
C          be required.
C
C (ILAENV) (output) INTEGER
C          >= 0: the value of the parameter specified by ISPEC
C          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
C
C  Further Details
C  ===============
C
C  The following conventions have been used when calling ILAENV from the
C  LAPACK routines:
C  1)  OPTS is a concatenation of all of the character options to
C      subroutine NAME, in the same order that they appear in the
C      argument list for NAME, even if they are not used in determining
C      the value of the parameter specified by ISPEC.
C  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
C      that they appear in the argument list for NAME.  N1 is used
C      first, N2 second, and so on, and unused problem dimensions are
C      passed a value of -1.
C  3)  The parameter value returned by ILAENV is checked for validity in
C      the calling subroutine.  For example, ILAENV is used to retrieve
C      the optimal blocksize for STRTRI as follows:
C
C      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
C      IF( NB.LE.1 ) NB = MAX( 1, N )
C
C  =====================================================================
C
C     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
C     ..
C     .. External Functions ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
C     ..
C     .. Executable Statements ..
C
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,
     $        1100 ) ISPEC
C
C     Invalid value for ISPEC
C
      ILAENV = -1
      RETURN
C
  100 CONTINUE
C
C     Convert NAME to upper case if the first character is lower case.
C
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
C
C        ASCII character set
C
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
C
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
C
C        EBCDIC character set
C
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
C
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
C
C        Prime machines:  ASCII+128
C
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
C
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
C
      GO TO ( 110, 200, 300 ) ISPEC
C
  110 CONTINUE
C
C     ISPEC = 1:  block size
C
C     In these examples, separate code is provided for setting NB for
C     real and complex.  We assume that NB will take the same value in
C     single or double precision.
C
      NB = 1
C
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
C
  200 CONTINUE
C
C     ISPEC = 2:  minimum block size
C
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
C
  300 CONTINUE
C
C     ISPEC = 3:  crossover point
C
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
C
  400 CONTINUE
C
C     ISPEC = 4:  number of shifts (used by xHSEQR)
C
      ILAENV = 6
      RETURN
C
  500 CONTINUE
C
C     ISPEC = 5:  minimum column dimension (not used)
C
      ILAENV = 2
      RETURN
C
  600 CONTINUE 
C
C     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
C
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
C
  700 CONTINUE
C
C     ISPEC = 7:  number of processors (not used)
C
      ILAENV = 1
      RETURN
C
  800 CONTINUE
C
C     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
C
      ILAENV = 50
      RETURN
C
  900 CONTINUE
C
C     ISPEC = 9:  maximum size of the subproblems at the bottom of the
C                 computation tree in the divide-and-conquer algorithm
C                 (used by xGELSD and xGESDD)
C
      ILAENV = 25
      RETURN
C
 1000 CONTINUE
C
C     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
C
      ILAENV = 1
      IF (ILAENV .EQ. 1) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 ) 
      ENDIF
      RETURN
C
 1100 CONTINUE
C
C     ISPEC = 11: infinity arithmetic can be trusted not to trap
C
      ILAENV = 1
      IF (ILAENV .EQ. 1) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 ) 
      ENDIF
      RETURN
C
C     End of ILAENV
C
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE ) 
C
C  -- LAPACK auxiliary routine (version 3.2) --
C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
C     November 2006
C
C     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ZERO, ONE
C     ..
C
C  Purpose
C  =======
C
C  IEEECK is called from the ILAENV to verify that Inifinity and 
C  possibly NaN arithmetic is safe (i.e. will not trap).
C
C  Arguments
C  =========
C
C  ISPEC   (input) INTEGER
C          Specifies whether to test just for inifinity arithmetic
C          or whether to test for infinity and NaN arithmetic.
C          = 0: Verify infinity arithmetic only.
C          = 1: Verify infinity and NaN arithmetic.
C
C  ZERO    (input) REAL
C          Must contain the value 0.0
C          This is passed to prevent the compiler from optimizing 
C          away this code.
C
C  ONE     (input) REAL
C          Must contain the value 1.0
C          This is passed to prevent the compiler from optimizing 
C          away this code.
C
C  RETURN VALUE:  INTEGER
C          = 0:  Arithmetic failed to produce the correct answers
C          = 1:  Arithmetic produced the correct answers
C
C     .. Local Scalars ..
      REAL POSINF, NEGINF, NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGZRO,
     $     NEWZRO
C     ..
C     .. Executable Statements ..
      IEEECK = 1

      POSINF = ONE /ZERO
      IF ( POSINF .LE. ONE ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      NEGINF = -ONE / ZERO
      IF ( NEGINF .GE. ZERO ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      NEGZRO = ONE / ( NEGINF + ONE )
      IF ( NEGZRO .NE. ZERO ) THEN
         IEEECK = 0
         RETURN
      ENDIF
         
      NEGINF = ONE / NEGZRO 
      IF ( NEGINF .GE. ZERO ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      NEWZRO = NEGZRO + ZERO
      IF ( NEWZRO .NE. ZERO ) THEN
         IEEECK = 0
         RETURN
      ENDIF
         
      POSINF = ONE / NEWZRO
      IF ( POSINF .LE. ONE ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      NEGINF = NEGINF * POSINF 
      IF ( NEGINF .GE. ZERO ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      POSINF = POSINF * POSINF 
      IF ( POSINF .LE. ONE ) THEN
         IEEECK = 0
         RETURN
      ENDIF



C
C     Return if we were only asked to check infinity arithmetic
C
      IF (ISPEC .EQ. 0 ) RETURN

      NAN1 = POSINF + NEGINF

      NAN2 = POSINF / NEGINF
      
      NAN3 = POSINF / POSINF
      
      NAN4 = POSINF * ZERO
      
      NAN5 = NEGINF * NEGZRO

      NAN6 = NAN5 * 0.0

      IF ( NAN1 .EQ. NAN1 ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      IF ( NAN2 .EQ. NAN2 ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      IF ( NAN3 .EQ. NAN3 ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      IF ( NAN4 .EQ. NAN4 ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      IF ( NAN5 .EQ. NAN5 ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      IF ( NAN6 .EQ. NAN6 ) THEN
         IEEECK = 0
         RETURN
      ENDIF

      RETURN
      END
