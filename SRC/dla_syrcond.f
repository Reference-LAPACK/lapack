      DOUBLE PRECISION FUNCTION DLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, 
     $                                       IPIV, CMODE, C, INFO, WORK,
     $                                       IWORK )
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
      CHARACTER          UPLO
      INTEGER            N, LDA, LDAF, INFO, CMODE
*     ..
*     .. Array Arguments
      INTEGER            IWORK( * ), IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * )
*     ..
*
*  Purpose
*  =======
*
*     DLA_SYRCOND estimates the Skeel condition number of  op(A) * op2(C)
*     where op2 is determined by CMODE as follows
*     CMODE =  1    op2(C) = C
*     CMODE =  0    op2(C) = I
*     CMODE = -1    op2(C) = inv(C)
*     The Skeel condition number cond(A) = norminf( |inv(A)||A| )
*     is computed by computing scaling factors R such that
*     diag(R)*A*op2(C) is row equilibrated and computing the standard
*     infinity-norm condition number.
*
*  Arguments
*  ==========
*
*     UPLO    (input) CHARACTER*1
*       = 'U':  Upper triangle of A is stored;
*       = 'L':  Lower triangle of A is stored.
*
*     N       (input) INTEGER
*     The number of linear equations, i.e., the order of the
*     matrix A.  N >= 0.
*
*     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*     On entry, the N-by-N matrix A.
*
*     LDA     (input) INTEGER
*     The leading dimension of the array A.  LDA >= max(1,N).
*
*     AF      (input) DOUBLE PRECISION array, dimension (LDAF,N)
*     The block diagonal matrix D and the multipliers used to
*     obtain the factor U or L as computed by DSYTRF.
*
*     LDAF    (input) INTEGER
*     The leading dimension of the array AF.  LDAF >= max(1,N).
*
*     IPIV    (input) INTEGER array, dimension (N)
*     Details of the interchanges and the block structure of D
*     as determined by DSYTRF.
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
*     WORK    (input) DOUBLE PRECISION array, dimension (3*N).
*     Workspace.
*
*     IWORK   (input) INTEGER array, dimension (N).
*     Workspace.
*
*  =====================================================================
*
*     .. Local Scalars ..
      CHARACTER          NORMIN
      INTEGER            KASE, I, J
      DOUBLE PRECISION   AINVNM, SMLNUM, TMP
      LOGICAL            UP
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACN2, DLATRS, DRSCL, XERBLA, DSYTRS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      DLA_SYRCOND = 0.0D+0
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLA_SYRCOND', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 ) THEN
         DLA_SYRCOND = 1.0D+0
         RETURN
      END IF
      UP = .FALSE.
      IF ( LSAME( UPLO, 'U' ) ) UP = .TRUE.
*
*     Compute the equilibration matrix R such that
*     inv(R)*A*C has unit 1-norm.
*
      IF ( UP ) THEN
         DO I = 1, N
            TMP = 0.0D+0
            IF ( CMODE .EQ. 1 ) THEN
               DO J = 1, I
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               END DO
            ELSE IF ( CMODE .EQ. 0 ) THEN
               DO J = 1, I
                  TMP = TMP + ABS( A( J, I ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( I, J ) )
               END DO
            ELSE
               DO J = 1, I
                  TMP = TMP + ABS( A( J, I ) / C( J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               END DO
            END IF
            WORK( 2*N+I ) = TMP
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0D+0
            IF ( CMODE .EQ. 1 ) THEN
               DO J = 1, I
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               END DO
            ELSE IF ( CMODE .EQ. 0 ) THEN
               DO J = 1, I
                  TMP = TMP + ABS( A( I, J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( J, I ) )
               END DO
            ELSE
               DO J = 1, I
                  TMP = TMP + ABS( A( I, J) / C( J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + ABS( A( J, I) / C( J ) )
               END DO
            END IF
            WORK( 2*N+I ) = TMP
         END DO
      ENDIF
*
*     Estimate the norm of inv(op(A)).
*
      SMLNUM = DLAMCH( 'Safe minimum' )
      AINVNM = 0.0D+0
      NORMIN = 'N'

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

            IF ( UP ) THEN
               CALL DSYTRS( 'U', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ELSE
               CALL DSYTRS( 'L', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ENDIF
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

            IF ( UP ) THEN
               CALL DSYTRS( 'U', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ELSE
               CALL DSYTRS( 'L', N, 1, AF, LDAF, IPIV, WORK, N, INFO )
            ENDIF
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO
         END IF
*
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM .NE. 0.0D+0 )
     $   DLA_SYRCOND = ( 1.0D+0 / AINVNM )
*
      RETURN
*
      END
