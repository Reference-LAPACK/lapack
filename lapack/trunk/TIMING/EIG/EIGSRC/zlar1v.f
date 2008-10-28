      SUBROUTINE ZLAR1V( N, B1, BN, SIGMA, D, L, LD, LLD, GERSCH, Z,
     $                   ZTZ, MINGMA, R, ISUPPZ, WORK )
*
*  -- LAPACK auxiliary routine (instru to count ops, version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            B1, BN, N, R
      DOUBLE PRECISION   MINGMA, SIGMA, ZTZ
*     ..
*     .. Array Arguments ..
      INTEGER            ISUPPZ( * )
      DOUBLE PRECISION   D( * ), GERSCH( * ), L( * ), LD( * ), LLD( * ),
     $                   WORK( * )
      COMPLEX*16         Z( * )
*     ..
*     Common block to return operation count
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
*
*  Purpose
*  =======
*
*  ZLAR1V computes the (scaled) r-th column of the inverse of
*  the sumbmatrix in rows B1 through BN of the tridiagonal matrix
*  L D L^T - sigma I. The following steps accomplish this computation :
*  (a) Stationary qd transform,  L D L^T - sigma I = L(+) D(+) L(+)^T,
*  (b) Progressive qd transform, L D L^T - sigma I = U(-) D(-) U(-)^T,
*  (c) Computation of the diagonal elements of the inverse of
*      L D L^T - sigma I by combining the above transforms, and choosing
*      r as the index where the diagonal of the inverse is (one of the)
*      largest in magnitude.
*  (d) Computation of the (scaled) r-th column of the inverse using the
*      twisted factorization obtained by combining the top part of the
*      the stationary and the bottom part of the progressive transform.
*
*  Arguments
*  =========
*
*  N        (input) INTEGER
*           The order of the matrix L D L^T.
*
*  B1       (input) INTEGER
*           First index of the submatrix of L D L^T.
*
*  BN       (input) INTEGER
*           Last index of the submatrix of L D L^T.
*
*  SIGMA    (input) DOUBLE PRECISION
*           The shift. Initially, when R = 0, SIGMA should be a good
*           approximation to an eigenvalue of L D L^T.
*
*  L        (input) DOUBLE PRECISION array, dimension (N-1)
*           The (n-1) subdiagonal elements of the unit bidiagonal matrix
*           L, in elements 1 to N-1.
*
*  D        (input) DOUBLE PRECISION array, dimension (N)
*           The n diagonal elements of the diagonal matrix D.
*
*  LD       (input) DOUBLE PRECISION array, dimension (N-1)
*           The n-1 elements L(i)*D(i).
*
*  LLD      (input) DOUBLE PRECISION array, dimension (N-1)
*           The n-1 elements L(i)*L(i)*D(i).
*
*  GERSCH   (input) DOUBLE PRECISION array, dimension (2*N)
*           The n Gerschgorin intervals. These are used to restrict
*           the initial search for R, when R is input as 0.
*
*  Z        (output) COMPLEX*16 array, dimension (N)
*           The (scaled) r-th column of the inverse. Z(R) is returned
*           to be 1.
*
*  ZTZ      (output) DOUBLE PRECISION
*           The square of the norm of Z.
*
*  MINGMA   (output) DOUBLE PRECISION
*           The reciprocal of the largest (in magnitude) diagonal
*           element of the inverse of L D L^T - sigma I.
*
*  R        (input/output) INTEGER
*           Initially, R should be input to be 0 and is then output as
*           the index where the diagonal element of the inverse is
*           largest in magnitude. In later iterations, this same value
*           of R should be input.
*
*  ISUPPZ   (output) INTEGER array, dimension (2)
*           The support of the vector in Z, i.e., the vector Z is
*           nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ).
*
*  WORK     (workspace) DOUBLE PRECISION array, dimension (4*N)
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Inderjit Dhillon, IBM Almaden, USA
*     Osni Marques, LBNL/NERSC, USA
*     Ken Stanley, Computer Science Division, University of
*       California at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLKSIZ
      PARAMETER          ( BLKSIZ = 32 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            SAWNAN
      INTEGER            FROM, I, INDP, INDS, INDUMN, J, R1, R2, TO
      DOUBLE PRECISION   DMINUS, DPLUS, EPS, S, TMP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      EPS = DLAMCH( 'Precision' )
      IF( R.EQ.0 ) THEN
*
*        Eliminate the top and bottom indices from the possible values
*        of R where the desired eigenvector is largest in magnitude.
*
         R1 = B1
         DO 10 I = B1, BN
            IF( SIGMA.GE.GERSCH( 2*I-1 ) .OR. SIGMA.LE.GERSCH( 2*I ) )
     $           THEN
               R1 = I
               GO TO 20
            END IF
   10    CONTINUE
   20    CONTINUE
         R2 = BN
         DO 30 I = BN, B1, -1
            IF( SIGMA.GE.GERSCH( 2*I-1 ) .OR. SIGMA.LE.GERSCH( 2*I ) )
     $           THEN
               R2 = I
               GO TO 40
            END IF
   30    CONTINUE
   40    CONTINUE
      ELSE
         R1 = R
         R2 = R
      END IF
*
      INDUMN = N
      INDS = 2*N + 1
      INDP = 3*N + 1
      SAWNAN = .FALSE.
*
*     Compute the stationary transform (using the differential form)
*     untill the index R2
*
      IF( B1.EQ.1 ) THEN
         WORK( INDS ) = ZERO
      ELSE
         WORK( INDS ) = LLD( B1-1 )
      END IF
      OPS = OPS + DBLE( 1 )
      S = WORK( INDS ) - SIGMA
      DO 50 I = B1, R2 - 1
         OPS = OPS + DBLE( 5 )
         DPLUS = D( I ) + S
         WORK( I ) = LD( I ) / DPLUS
         WORK( INDS+I ) = S*WORK( I )*L( I )
         S = WORK( INDS+I ) - SIGMA
   50 CONTINUE
*
      IF( .NOT.( S.GT.ZERO .OR. S.LT.ONE ) ) THEN
*
*        Run a slower version of the above loop if a NaN is detected
*
         SAWNAN = .TRUE.
         J = B1 + 1
   60    CONTINUE
         IF( WORK( INDS+J ).GT.ZERO .OR. WORK( INDS+J ).LT.ONE ) THEN
            J = J + 1
            GO TO 60
         END IF
         WORK( INDS+J ) = LLD( J )
         S = WORK( INDS+J ) - SIGMA
         DO 70 I = J + 1, R2 - 1
            OPS = OPS + DBLE( 3 )
            DPLUS = D( I ) + S
            WORK( I ) = LD( I ) / DPLUS
            IF( WORK( I ).EQ.ZERO ) THEN
               WORK( INDS+I ) = LLD( I )
            ELSE
               OPS = OPS + DBLE( 2 )
               WORK( INDS+I ) = S*WORK( I )*L( I )
            END IF
            S = WORK( INDS+I ) - SIGMA
   70    CONTINUE
      END IF
      OPS = OPS + DBLE( 1 )
      WORK( INDP+BN-1 ) = D( BN ) - SIGMA
      DO 80 I = BN - 1, R1, -1
         OPS = OPS + DBLE( 5 )
         DMINUS = LLD( I ) + WORK( INDP+I )
         TMP = D( I ) / DMINUS
         WORK( INDUMN+I ) = L( I )*TMP
         WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - SIGMA
   80 CONTINUE
      TMP = WORK( INDP+R1-1 )
      IF( .NOT.( TMP.GT.ZERO .OR. TMP.LT.ONE ) ) THEN
*
*        Run a slower version of the above loop if a NaN is detected
*
         SAWNAN = .TRUE.
         J = BN - 3
   90    CONTINUE
         IF( WORK( INDP+J ).GT.ZERO .OR. WORK( INDP+J ).LT.ONE ) THEN
            J = J - 1
            GO TO 90
         END IF
         OPS = OPS + DBLE( 1 )
         WORK( INDP+J ) = D( J+1 ) - SIGMA
         DO 100 I = J, R1, -1
            OPS = OPS + DBLE( 3 )
            DMINUS = LLD( I ) + WORK( INDP+I )
            TMP = D( I ) / DMINUS
            WORK( INDUMN+I ) = L( I )*TMP
            IF( TMP.EQ.ZERO ) THEN
               OPS = OPS + DBLE( 1 )
               WORK( INDP+I-1 ) = D( I ) - SIGMA
            ELSE
               OPS = OPS + DBLE( 2 )
               WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - SIGMA
            END IF
  100    CONTINUE
      END IF
*
*     Find the index (from R1 to R2) of the largest (in magnitude)
*     diagonal element of the inverse
*
      MINGMA = WORK( INDS+R1-1 ) + WORK( INDP+R1-1 )
      IF( MINGMA.EQ.ZERO )
     $   MINGMA = EPS*WORK( INDS+R1-1 )
      R = R1
      DO 110 I = R1, R2 - 1
         OPS = OPS + DBLE( 1 )
         TMP = WORK( INDS+I ) + WORK( INDP+I )
         IF( TMP.EQ.ZERO ) THEN
            OPS = OPS + DBLE( 1 )
            TMP = EPS*WORK( INDS+I )
         END IF
         IF( ABS( TMP ).LT.ABS( MINGMA ) ) THEN
            MINGMA = TMP
            R = I + 1
         END IF
  110 CONTINUE
*
*     Compute the (scaled) r-th column of the inverse
*
      ISUPPZ( 1 ) = B1
      ISUPPZ( 2 ) = BN
      Z( R ) = CONE
      ZTZ = ONE
      IF( .NOT.SAWNAN ) THEN
         FROM = R - 1
         TO = MAX( R-BLKSIZ, B1 )
  120    CONTINUE
         IF( FROM.GE.B1 ) THEN
            DO 130 I = FROM, TO, -1
               OPS = OPS + DBLE( 9 )
               Z( I ) = -( WORK( I )*Z( I+1 ) )
               ZTZ = ZTZ + DBLE( Z( I )*Z( I ) )
  130       CONTINUE
            IF( ABS( Z( TO ) ).LE.EPS .AND. ABS( Z( TO+1 ) ).LE.EPS )
     $           THEN
               ISUPPZ( 1 ) = TO + 2
            ELSE
               FROM = TO - 1
               TO = MAX( TO-BLKSIZ, B1 )
               GO TO 120
            END IF
         END IF
         FROM = R + 1
         TO = MIN( R+BLKSIZ, BN )
  140    CONTINUE
         IF( FROM.LE.BN ) THEN
            DO 150 I = FROM, TO
               OPS = OPS + DBLE( 9 )
               Z( I ) = -( WORK( INDUMN+I-1 )*Z( I-1 ) )
               ZTZ = ZTZ + DBLE( Z( I )*Z( I ) )
  150       CONTINUE
            IF( ABS( Z( TO ) ).LE.EPS .AND. ABS( Z( TO-1 ) ).LE.EPS )
     $           THEN
               ISUPPZ( 2 ) = TO - 2
            ELSE
               FROM = TO + 1
               TO = MIN( TO+BLKSIZ, BN )
               GO TO 140
            END IF
         END IF
      ELSE
         DO 160 I = R - 1, B1, -1
            IF( Z( I+1 ).EQ.ZERO ) THEN
               OPS = OPS + DBLE( 3 )
               Z( I ) = -( LD( I+1 ) / LD( I ) )*Z( I+2 )
            ELSE IF( ABS( Z( I+1 ) ).LE.EPS .AND. ABS( Z( I+2 ) ).LE.
     $               EPS ) THEN
               ISUPPZ( 1 ) = I + 3
               GO TO 170
            ELSE
               OPS = OPS + DBLE( 2 )
               Z( I ) = -( WORK( I )*Z( I+1 ) )
            END IF
            OPS = OPS + DBLE( 7 )
            ZTZ = ZTZ + DBLE( Z( I )*Z( I ) )
  160    CONTINUE
  170    CONTINUE
         DO 180 I = R, BN - 1
            IF( Z( I ).EQ.ZERO ) THEN
               OPS = OPS + DBLE( 3 )
               Z( I+1 ) = -( LD( I-1 ) / LD( I ) )*Z( I-1 )
            ELSE IF( ABS( Z( I ) ).LE.EPS .AND. ABS( Z( I-1 ) ).LE.EPS )
     $                THEN
               ISUPPZ( 2 ) = I - 2
               GO TO 190
            ELSE
               OPS = OPS + DBLE( 2 )
               Z( I+1 ) = -( WORK( INDUMN+I )*Z( I ) )
            END IF
            OPS = OPS + DBLE( 7 )
            ZTZ = ZTZ + DBLE( Z( I+1 )*Z( I+1 ) )
  180    CONTINUE
  190    CONTINUE
      END IF
      DO 200 I = B1, ISUPPZ( 1 ) - 3
         Z( I ) = ZERO
  200 CONTINUE
      DO 210 I = ISUPPZ( 2 ) + 3, BN
         Z( I ) = ZERO
  210 CONTINUE
*
      RETURN
*
*     End of ZLAR1V
*
      END
