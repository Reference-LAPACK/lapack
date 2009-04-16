      DOUBLE PRECISION FUNCTION DLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB,
     $                                       AFB, LDAFB, IPIV, CMODE, C,
     $                                       INFO, WORK, IWORK )
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
      INTEGER            N, LDAB, LDAFB, INFO, KL, KU, CMODE
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ),
     $                   C( * )
*     ..
*
*  Purpose
*  =======
*
*     DLA_GERCOND Estimates the Skeel condition number of  op(A) * op2(C)
*     where op2 is determined by CMODE as follows
*     CMODE =  1    op2(C) = C
*     CMODE =  0    op2(C) = I
*     CMODE = -1    op2(C) = inv(C)
*     The Skeel condition number  cond(A) = norminf( |inv(A)||A| )
*     is computed by computing scaling factors R such that
*     diag(R)*A*op2(C) is row equilibrated and computing the standard
*     infinity-norm condition number.
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
*     AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
*     The j-th column of A is stored in the j-th column of the
*     array AB as follows:
*     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
*
*     LDAB    (input) INTEGER
*     The leading dimension of the array AB.  LDAB >= KL+KU+1.
*
*     AFB     (input) DOUBLE PRECISION array, dimension (LDAFB,N)
*     Details of the LU factorization of the band matrix A, as
*     computed by DGBTRF.  U is stored as an upper triangular
*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
*     and the multipliers used during the factorization are stored
*     in rows KL+KU+2 to 2*KL+KU+1.
*
*     LDAFB   (input) INTEGER
*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.
*
*     IPIV    (input) INTEGER array, dimension (N)
*     The pivot indices from the factorization A = P*L*U
*     as computed by DGBTRF; row i of the matrix was interchanged
*     with row IPIV(i).
*
*     CMODE   (input) INTEGER
*     Determines op2(C) in the formula op(A) * op2(C) as follows:
*     CMODE =  1    op2(C) = C
*     CMODE =  0    op2(C) = I
*     CMODE = -1    op2(C) = inv(C)
*
*     C       (input) DOUBLE PRECISION array, dimension (N)
*     The vector C in the formula op(A) * op2(C).
*
*     INFO    (output) INTEGER
*       = 0:  Successful exit.
*     i > 0:  The ith argument is invalid.
*
*     WORK    (input) DOUBLE PRECISION array, dimension (5*N).
*     Workspace.
*
*     IWORK   (input) INTEGER array, dimension (N).
*     Workspace.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            NOTRANS
      INTEGER            KASE, I, J, KD, KE
      DOUBLE PRECISION   AINVNM, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACN2, DGBTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      DLA_GBRCOND = 0.0D+0
*
      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      IF ( .NOT. NOTRANS .AND. .NOT. LSAME(TRANS, 'T')
     $     .AND. .NOT. LSAME(TRANS, 'C') ) THEN
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
         CALL XERBLA( 'DLA_GBRCOND', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 ) THEN
         DLA_GBRCOND = 1.0D+0
         RETURN
      END IF
*
*     Compute the equilibration matrix R such that
*     inv(R)*A*C has unit 1-norm.
*
      KD = KU + 1
      KE = KL + 1
      IF ( NOTRANS ) THEN
         DO I = 1, N
            TMP = 0.0D+0
               IF ( CMODE .EQ. 1 ) THEN
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                     TMP = TMP + ABS( AB( KD+I-J, J ) * C( J ) )
                  END DO
               ELSE IF ( CMODE .EQ. 0 ) THEN
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                     TMP = TMP + ABS( AB( KD+I-J, J ) )
                  END DO
               ELSE
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                     TMP = TMP + ABS( AB( KD+I-J, J ) / C( J ) )
                  END DO
               END IF
            WORK( 2*N+I ) = TMP
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0D+0
            IF ( CMODE .EQ. 1 ) THEN
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + ABS( AB( KE-I+J, I ) * C( J ) )
               END DO
            ELSE IF ( CMODE .EQ. 0 ) THEN
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + ABS( AB( KE-I+J, I ) )
               END DO
            ELSE
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + ABS( AB( KE-I+J, I ) / C( J ) )
               END DO
            END IF
            WORK( 2*N+I ) = TMP
         END DO
      END IF
*
*     Estimate the norm of inv(op(A)).
*
      AINVNM = 0.0D+0

      KASE = 0
   10 CONTINUE
      CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.2 ) THEN
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO

            IF ( NOTRANS ) THEN
               CALL DGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            ELSE
               CALL DGBTRS( 'Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV,
     $              WORK, N, INFO )
            END IF
*
*           Multiply by inv(C).
*
            IF ( CMODE .EQ. 1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            ELSE IF ( CMODE .EQ. -1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF
         ELSE
*
*           Multiply by inv(C').
*
            IF ( CMODE .EQ. 1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            ELSE IF ( CMODE .EQ. -1 ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF

            IF ( NOTRANS ) THEN
               CALL DGBTRS( 'Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV,
     $              WORK, N, INFO )
            ELSE
               CALL DGBTRS( 'No transpose', N, KL, KU, 1, AFB, LDAFB,
     $              IPIV, WORK, N, INFO )
            END IF
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM .NE. 0.0D+0 )
     $   DLA_GBRCOND = ( 1.0D+0 / AINVNM )
*
      RETURN
*
      END
