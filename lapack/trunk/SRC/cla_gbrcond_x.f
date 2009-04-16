      REAL FUNCTION CLA_GBRCOND_X( TRANS, N, KL, KU, AB, LDAB, AFB,
     $                             LDAFB, IPIV, X, INFO, WORK, RWORK )
*
*     -- LAPACK routine (version 3.2.1)                               --
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
      CHARACTER          TRANS
      INTEGER            N, KL, KU, KD, KE, LDAB, LDAFB, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ),
     $                   X( * )
      REAL               RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*     CLA_GBRCOND_X Computes the infinity norm condition number of
*     op(A) * diag(X) where X is a COMPLEX vector.
*
*  Arguments
*  =========
*
*     TRANS   (input) CHARACTER*1
*     Specifies the form of the system of equations:
*       = 'N':  A * X = B     (No transpose)
*       = 'T':  A**T * X = B  (Transpose)
*       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)
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
*     IPIV    (input) INTEGER array, dimension (N)
*     The pivot indices from the factorization A = P*L*U
*     as computed by CGBTRF; row i of the matrix was interchanged
*     with row IPIV(i).
*
*     X       (input) COMPLEX array, dimension (N)
*     The vector X in the formula op(A) * diag(X).
*
*     INFO    (output) INTEGER
*       = 0:  Successful exit.
*     i > 0:  The ith argument is invalid.
*
*     WORK    (input) COMPLEX array, dimension (2*N).
*     Workspace.
*
*     RWORK   (input) REAL array, dimension (N).
*     Workspace.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            NOTRANS
      INTEGER            KASE, I, J
      REAL               AINVNM, ANORM, TMP
      COMPLEX            ZDUM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACN2, CGBTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      CLA_GBRCOND_X = 0.0E+0
*
      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      IF ( .NOT. NOTRANS .AND. .NOT. LSAME(TRANS, 'T') .AND. .NOT.
     $     LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 .OR. KL.GT.N-1 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 .OR. KU.GT.N-1 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KU+1 ) THEN
         INFO = -6
      ELSE IF( LDAFB.LT.2*KL+KU+1 ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLA_GBRCOND_X', -INFO )
         RETURN
      END IF
*
*     Compute norm of op(A)*op2(C).
*
      KD = KU + 1
      KE = KL + 1
      ANORM = 0.0
      IF ( NOTRANS ) THEN
         DO I = 1, N
            TMP = 0.0E+0
            DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
               TMP = TMP + CABS1( AB( KD+I-J, J) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0E+0
            DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
               TMP = TMP + CABS1( AB( KE-I+J, I ) * X( J ) )
            END DO
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         CLA_GBRCOND_X = 1.0E+0
         RETURN
      ELSE IF( ANORM .EQ. 0.0E+0 ) THEN
         RETURN
      END IF
*
*     Estimate the norm of inv(op(A)).
*
      AINVNM = 0.0E+0
*
      KASE = 0
   10 CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.2 ) THEN
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
*
            IF ( NOTRANS ) THEN
               CALL CGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            ELSE
               CALL CGBTRS( 'Conjugate transpose', N, KL, KU, 1, AFB,
     $              LDAFB, IPIV, WORK, N, INFO )
            ENDIF
*
*           Multiply by inv(X).
*
            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO
         ELSE
*
*           Multiply by inv(X').
*
            DO I = 1, N
               WORK( I ) = WORK( I ) / X( I )
            END DO
*
            IF ( NOTRANS ) THEN
               CALL CGBTRS( 'Conjugate transpose', N, KL, KU, 1, AFB,
     $              LDAFB, IPIV, WORK, N, INFO )
            ELSE
               CALL CGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            END IF
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM .NE. 0.0E+0 )
     $   CLA_GBRCOND_X = 1.0E+0 / AINVNM
*
      RETURN
*
      END
