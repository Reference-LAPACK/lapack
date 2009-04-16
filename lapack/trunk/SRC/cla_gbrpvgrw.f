      REAL FUNCTION CLA_GBRPVGRW( N, KL, KU, NCOLS, AB, LDAB, AFB,
     $                            LDAFB )
*
*     -- LAPACK routine (version 3.2.1)                                 --
*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and --
*     -- Jason Riedy of Univ. of California Berkeley.                 --
*     -- April 2009                                                   --
*
*     -- LAPACK is a software package provided by Univ. of Tennessee, --
*     -- Univ. of California Berkeley and NAG Ltd.                    --
*
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      INTEGER            N, KL, KU, NCOLS, LDAB, LDAFB
*     ..
*     .. Array Arguments ..
      COMPLEX            AB( LDAB, * ), AFB( LDAFB, * )
*     ..
*
*  Purpose
*  =======
*
*  CLA_GBRPVGRW computes the reciprocal pivot growth factor
*  norm(A)/norm(U). The "max absolute element" norm is used. If this is
*  much less than 1, the stability of the LU factorization of the
*  (equilibrated) matrix A could be poor. This also means that the
*  solution X, estimated condition numbers, and error bounds could be
*  unreliable.
*
*  Arguments
*  =========
*
*     N       (input) INTEGER
*     The number of linear equations, i.e., the order of the
*     matrix A.  N >= 0.
*
*     KL      (input) INTEGER
*     The number of subdiagonals within the band of A.  KL >= 0.
*
*     KU      (input) INTEGER
*     The number of superdiagonals within the band of A.  KU >= 0.
*
*     NCOLS   (input) INTEGER
*     The number of columns of the matrix A.  NCOLS >= 0.
*
*     AB      (input) COMPLEX array, dimension (LDAB,N)
*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
*     The j-th column of A is stored in the j-th column of the
*     array AB as follows:
*     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
*
*     LDAB    (input) INTEGER
*     The leading dimension of the array AB.  LDAB >= KL+KU+1.
*
*     AFB     (input) COMPLEX array, dimension (LDAFB,N)
*     Details of the LU factorization of the band matrix A, as
*     computed by CGBTRF.  U is stored as an upper triangular
*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
*     and the multipliers used during the factorization are stored
*     in rows KL+KU+2 to 2*KL+KU+1.
*
*     LDAFB   (input) INTEGER
*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, KD
      REAL               AMAX, UMAX, RPVGRW
      COMPLEX            ZDUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, REAL, AIMAG
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      RPVGRW = 1.0

      KD = KU + 1
      DO J = 1, NCOLS
         AMAX = 0.0
         UMAX = 0.0
         DO I = MAX( J-KU, 1 ), MIN( J+KL, N )
            AMAX = MAX( CABS1( AB( KD+I-J, J ) ), AMAX )
         END DO
         DO I = MAX( J-KU, 1 ), J
            UMAX = MAX( CABS1( AFB( KD+I-J, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      CLA_GBRPVGRW = RPVGRW
      END FUNCTION
