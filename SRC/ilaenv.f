*> \brief \b ILAENV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ILAENV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*       .. Scalar Arguments ..
*       CHARACTER*( * )    NAME, OPTS
*       INTEGER            ISPEC, N1, N2, N3, N4
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ILAENV is called from the LAPACK routines to choose problem-dependent
*> parameters for the local environment.  See ISPEC for a description of
*> the parameters.
*>
*> ILAENV returns an INTEGER
*> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
*> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*>
*> This version provides a set of parameters which should give good,
*> but not optimal, performance on many of the currently available
*> computers.  Users are encouraged to modify this subroutine to set
*> the tuning parameters for their particular machine using the option
*> and problem size information in the arguments.
*>
*> Routine names are compared case-insensitively using LSAME and LSAMEN;
*> no conversion of NAME to upper case is needed.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ISPEC
*> \verbatim
*>          ISPEC is INTEGER
*>          Specifies the parameter to be returned as the value of
*>          ILAENV.
*>          = 1: the optimal blocksize; if this value is 1, an unblocked
*>               algorithm will give the best performance.
*>          = 2: the minimum block size for which the block routine
*>               should be used; if the usable block size is less than
*>               this value, an unblocked routine should be used.
*>          = 3: the crossover point (in a block routine, for N less
*>               than this value, an unblocked routine should be used)
*>          = 4: the number of shifts, used in the nonsymmetric
*>               eigenvalue routines (DEPRECATED)
*>          = 5: the minimum column dimension for blocking to be used;
*>               rectangular blocks must have dimension at least k by m,
*>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*>          = 6: the crossover point for the SVD (when reducing an m by n
*>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*>               this value, a QR factorization is used first to reduce
*>               the matrix to a triangular form.)
*>          = 7: the number of processors
*>          = 8: the crossover point for the multishift QR method
*>               for nonsymmetric eigenvalue problems (DEPRECATED)
*>          = 9: maximum size of the subproblems at the bottom of the
*>               computation tree in the divide-and-conquer algorithm
*>               (used by xGELSD and xGESDD)
*>          =10: ieee infinity and NaN arithmetic can be trusted not to trap
*>          =11: infinity arithmetic can be trusted not to trap
*>          12 <= ISPEC <= 17:
*>               xHSEQR or related subroutines,
*>               see IPARMQ for detailed explanation
*> \endverbatim
*>
*> \param[in] NAME
*> \verbatim
*>          NAME is CHARACTER*(*)
*>          The name of the calling subroutine, in upper, lower, or
*>          mixed case.
*> \endverbatim
*>
*> \param[in] OPTS
*> \verbatim
*>          OPTS is CHARACTER*(*)
*>          The character options to the subroutine NAME, concatenated
*>          into a single character string.  For example, UPLO = 'U',
*>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*>          be specified as OPTS = 'UTN'.
*> \endverbatim
*>
*> \param[in] N1
*> \verbatim
*>          N1 is INTEGER
*> \endverbatim
*>
*> \param[in] N2
*> \verbatim
*>          N2 is INTEGER
*> \endverbatim
*>
*> \param[in] N3
*> \verbatim
*>          N3 is INTEGER
*> \endverbatim
*>
*> \param[in] N4
*> \verbatim
*>          N4 is INTEGER
*>          Problem dimensions for the subroutine NAME; these may not all
*>          be required.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup ilaenv
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The following conventions have been used when calling ILAENV from the
*>  LAPACK routines:
*>  1)  OPTS is a concatenation of all of the character options to
*>      subroutine NAME, in the same order that they appear in the
*>      argument list for NAME, even if they are not used in determining
*>      the value of the parameter specified by ISPEC.
*>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*>      that they appear in the argument list for NAME.  N1 is used
*>      first, N2 second, and so on, and unused problem dimensions are
*>      passed a value of -1.
*>  3)  The parameter value returned by ILAENV is checked for validity in
*>      the calling subroutine.  For example, ILAENV is used to retrieve
*>      the optimal blocksize for STRTRI as follows:
*>
*>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*>      IF( NB.LE.1 ) NB = MAX( 1, N )
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            NB, NBMIN, NX
      LOGICAL            CNAME, SNAME, TWOSTAGE
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*16
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      LOGICAL            LSAME, LSAMEN
      EXTERNAL           IEEECK, IPARMQ, LSAME, LSAMEN
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160, 160)ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Compare names without changing their case.
*
      ILAENV = 1
      SUBNAM = NAME
*
      C1 = SUBNAM( 1: 1 )
      SNAME = LSAME( C1, 'S' ) .OR. LSAME( C1, 'D' )
      CNAME = LSAME( C1, 'C' ) .OR. LSAME( C1, 'Z' )
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
      TWOSTAGE = LEN( SUBNAM ).GE.11 .AND. LSAME( SUBNAM( 11: 11 ), '2'
     $    )
*
      GO TO ( 50, 60, 70 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( LSAMEN( 5, SUBNAM(2:6), 'LAORH' ) ) THEN
*
*        This is for *LAORHR_GETRFNP routine
*
         IF( SNAME ) THEN
             NB = 32
         ELSE
             NB = 32
         END IF
      ELSE IF( LSAMEN( 5, SUBNAM(2:6), 'TZRZF' ) ) THEN
         IF( SNAME ) THEN
             NB = 64
         ELSE
             NB = 32
         END IF
      ELSE IF( LSAMEN( 5, SUBNAM(2:6), 'ORMRZ' ) ) THEN
             NB = 64
      ELSE IF( LSAMEN( 5, SUBNAM(2:6), 'UNMRZ' ) ) THEN
             NB = 32
      ELSE IF( LSAMEN( 2, C2, 'GE' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( LSAMEN( 3, C3, 'QRF' ) .OR. LSAMEN( 3, C3, 'RQF' )
     $       .OR. LSAMEN( 3, C3, 'LQF' ) .OR. LSAMEN( 3, C3, 'QLF' ) )
     $       THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( LSAMEN( 3, C3, 'QR ' )) THEN
            IF( N3 .EQ. 1) THEN
               IF( SNAME ) THEN
*     M*N
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               ELSE
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               END IF
            ELSE
               IF( SNAME ) THEN
                  NB = 1
               ELSE
                  NB = 1
               END IF
            END IF
         ELSE IF( LSAMEN( 3, C3, 'LQ ' )) THEN
            IF( N3 .EQ. 2) THEN
               IF( SNAME ) THEN
*     M*N
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               ELSE
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               END IF
            ELSE
               IF( SNAME ) THEN
                  NB = 1
               ELSE
                  NB = 1
               END IF
            END IF
         ELSE IF( LSAMEN( 3, C3, 'HRD' ) ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( LSAMEN( 3, C3, 'BRD' ) ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( LSAMEN( 3, C3, 'TRI' ) ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( LSAMEN( 5, SUBNAM( 4: 7 ), 'QP3RK' ) ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      ELSE IF( LSAMEN( 2, C2, 'PO' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( LSAMEN( 2, C2, 'SY' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            IF( SNAME ) THEN
               IF( TWOSTAGE ) THEN
                  NB = 192
               ELSE
                  NB = 64
               END IF
            ELSE
               IF( TWOSTAGE ) THEN
                  NB = 192
               ELSE
                  NB = 64
               END IF
            END IF
         ELSE IF( SNAME .AND. LSAMEN( 3, C3, 'TRD' ) ) THEN
            NB = 32
         ELSE IF( SNAME .AND. LSAMEN( 3, C3, 'GST' ) ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. LSAMEN( 2, C2, 'HE' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            IF( TWOSTAGE ) THEN
               NB = 192
            ELSE
               NB = 64
            END IF
         ELSE IF( LSAMEN( 3, C3, 'TRD' ) ) THEN
            NB = 32
         ELSE IF( LSAMEN( 3, C3, 'GST' ) ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. LSAMEN( 2, C2, 'OR' ) ) THEN
         IF( LSAME( C3( 1: 1 ), 'G' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NB = 32
            END IF
         ELSE IF( LSAME( C3( 1: 1 ), 'M' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. LSAMEN( 2, C2, 'UN' ) ) THEN
         IF( LSAME( C3( 1: 1 ), 'G' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NB = 32
            END IF
         ELSE IF( LSAME( C3( 1: 1 ), 'M' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
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
      ELSE IF( LSAMEN( 2, C2, 'PB' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
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
      ELSE IF( LSAMEN( 2, C2, 'TR' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRI' ) ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF ( LSAMEN( 3, C3, 'EVC' ) ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( LSAMEN( 3, C3, 'SYL' ) ) THEN
*           The upper bound is to prevent overly aggressive scaling.
            IF( SNAME ) THEN
               NB = MIN( MAX( 48, INT( ( MIN( N1, N2 ) * 16 ) / 100) ),
     $                   240 )
            ELSE
               NB = MIN( MAX( 24, INT( ( MIN( N1, N2 ) * 8 ) / 100) ),
     $                   80 )
            END IF
         END IF
      ELSE IF( LSAMEN( 2, C2, 'LA' ) ) THEN
         IF( LSAMEN( 3, C3, 'UUM' ) ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( LSAMEN( 3, C3, 'LLM' ) ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. LSAMEN( 2, C2, 'ST' ) ) THEN
         IF( LSAMEN( 3, C3, 'EBZ' ) ) THEN
            NB = 1
         END IF
      ELSE IF( LSAMEN( 2, C2, 'GG' ) ) THEN
         NB = 32
         IF( LSAMEN( 3, C3, 'HD3' ) ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      END IF
      ILAENV = NB
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( LSAMEN( 2, C2, 'GE' ) ) THEN
         IF( LSAMEN( 3, C3, 'QRF' ) .OR. LSAMEN( 3, C3, 'RQF' ) .OR.
     $       LSAMEN( 3, C3, 'LQF' ) .OR. LSAMEN( 3, C3, 'QLF' ) ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( LSAMEN( 3, C3, 'HRD' ) ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( LSAMEN( 3, C3, 'BRD' ) ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( LSAMEN( 3, C3, 'TRI' ) ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( LSAMEN( 5, SUBNAM( 4: 7 ), 'QP3RK' ) ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF

      ELSE IF( LSAMEN( 2, C2, 'SY' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. LSAMEN( 3, C3, 'TRD' ) ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. LSAMEN( 2, C2, 'HE' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRD' ) ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. LSAMEN( 2, C2, 'OR' ) ) THEN
         IF( LSAME( C3( 1: 1 ), 'G' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NBMIN = 2
            END IF
         ELSE IF( LSAME( C3( 1: 1 ), 'M' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. LSAMEN( 2, C2, 'UN' ) ) THEN
         IF( LSAME( C3( 1: 1 ), 'G' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NBMIN = 2
            END IF
         ELSE IF( LSAME( C3( 1: 1 ), 'M' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( LSAMEN( 2, C2, 'GG' ) ) THEN
         NBMIN = 2
         IF( LSAMEN( 3, C3, 'HD3' ) ) THEN
            NBMIN = 2
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( LSAMEN( 2, C2, 'GE' ) ) THEN
         IF( LSAMEN( 3, C3, 'QRF' ) .OR. LSAMEN( 3, C3, 'RQF' ) .OR.
     $       LSAMEN( 3, C3, 'LQF' ) .OR. LSAMEN( 3, C3, 'QLF' ) ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( LSAMEN( 3, C3, 'HRD' ) ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( LSAMEN( 3, C3, 'BRD' ) ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( LSAMEN( 5, SUBNAM( 4: 7 ), 'QP3RK' ) ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( LSAMEN( 2, C2, 'SY' ) ) THEN
         IF( SNAME .AND. LSAMEN( 3, C3, 'TRD' ) ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. LSAMEN( 2, C2, 'HE' ) ) THEN
         IF( LSAMEN( 3, C3, 'TRD' ) ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. LSAMEN( 2, C2, 'OR' ) ) THEN
         IF( LSAME( C3( 1: 1 ), 'G' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. LSAMEN( 2, C2, 'UN' ) ) THEN
         IF( LSAME( C3( 1: 1 ), 'G' ) ) THEN
            IF( LSAMEN( 2, C4, 'QR' ) .OR. LSAMEN( 2, C4, 'RQ' ) .OR.
     $          LSAMEN( 2, C4, 'LQ' ) .OR. LSAMEN( 2, C4, 'QL' ) .OR.
     $          LSAMEN( 2, C4, 'HR' ) .OR. LSAMEN( 2, C4, 'TR' ) .OR.
     $          LSAMEN( 2, C4, 'BR' ) ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( LSAMEN( 2, C2, 'GG' ) ) THEN
         NX = 128
         IF( LSAMEN( 3, C3, 'HD3' ) ) THEN
            NX = 128
         END IF
      ELSE IF( LSAMEN( 2, C2, 'LA' ) ) THEN
         IF( LSAMEN( 3, C3, 'RFT' ) ) THEN
            NX = 64
         ELSE IF( LSAMEN( 3, C3, 'LLM' ) ) THEN
*
*           Note that this matches on *lallm*, but we only use it in
*           *lallm_recursive at the moment
         ELSE IF( LSAMEN( 3, C3, 'UUM' ) ) THEN
*
*           Note that this matches on *lauum_recursive, which is the
*           only place we currently use this
*
            NX = 64
         END IF
      END IF
      ILAENV = NX
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: ieee infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 17: xHSEQR or related subroutines.
*
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
*
*     End of ILAENV
*
      END
