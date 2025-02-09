*> \brief \b SLA_SYRCOND estimates the Skeel condition number for a symmetric indefinite matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SLA_SYRCOND + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syrcond.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syrcond.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syrcond.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       REAL FUNCTION SLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV, CMODE,
*                                  C, INFO, WORK, IWORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            N, LDA, LDAF, INFO, CMODE
*       ..
*       .. Array Arguments
*       INTEGER            IWORK( * ), IPIV( * )
*       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    SLA_SYRCOND estimates the Skeel condition number of  op(A) * op2(C)
*>    where op2 is determined by CMODE as follows
*>    CMODE =  1    op2(C) = C
*>    CMODE =  0    op2(C) = I
*>    CMODE = -1    op2(C) = inv(C)
*>    The Skeel condition number cond(A) = norminf( |inv(A)||A| )
*>    is computed by computing scaling factors R such that
*>    diag(R)*A*op2(C) is row equilibrated and computing the standard
*>    infinity-norm condition number.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>       = 'U':  Upper triangle of A is stored;
*>       = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>     The number of linear equations, i.e., the order of the
*>     matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>     On entry, the N-by-N matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>     The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] AF
*> \verbatim
*>          AF is REAL array, dimension (LDAF,N)
*>     The block diagonal matrix D and the multipliers used to
*>     obtain the factor U or L as computed by SSYTRF.
*> \endverbatim
*>
*> \param[in] LDAF
*> \verbatim
*>          LDAF is INTEGER
*>     The leading dimension of the array AF.  LDAF >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>     Details of the interchanges and the block structure of D
*>     as determined by SSYTRF.
*> \endverbatim
*>
*> \param[in] CMODE
*> \verbatim
*>          CMODE is INTEGER
*>     Determines op2(C) in the formula op(A) * op2(C) as follows:
*>     CMODE =  1    op2(C) = C
*>     CMODE =  0    op2(C) = I
*>     CMODE = -1    op2(C) = inv(C)
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is REAL array, dimension (N)
*>     The vector C in the formula op(A) * op2(C).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>       = 0:  Successful exit.
*>     i > 0:  The ith argument is invalid.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (3*N).
*>     Workspace.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N).
*>     Workspace.
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
*> \ingroup la_hercond
*
*  =====================================================================
      REAL FUNCTION SLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV,
     $                           CMODE,
     $                           C, INFO, WORK, IWORK )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            N, LDA, LDAF, INFO, CMODE
*     ..
*     .. Array Arguments
      INTEGER            IWORK( * ), IPIV( * )
      REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      CHARACTER          NORMIN
      INTEGER            KASE, I, J
      REAL               AINVNM, SMLNUM, TMP
      LOGICAL            UP
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLACN2, XERBLA, SSYTRS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      SLA_SYRCOND = 0.0
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLA_SYRCOND', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 ) THEN
         SLA_SYRCOND = 1.0
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
            TMP = 0.0
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
            TMP = 0.0
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
      SMLNUM = SLAMCH( 'Safe minimum' )
      AINVNM = 0.0
      NORMIN = 'N'

      KASE = 0
   10 CONTINUE
      CALL SLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.2 ) THEN
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO

            IF ( UP ) THEN
               CALL SSYTRS( 'U', N, 1, AF, LDAF, IPIV, WORK, N,
     $                      INFO )
            ELSE
               CALL SSYTRS( 'L', N, 1, AF, LDAF, IPIV, WORK, N,
     $                      INFO )
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
*           Multiply by inv(C**T).
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
               CALL SSYTRS( 'U', N, 1, AF, LDAF, IPIV, WORK, N,
     $                      INFO )
            ELSE
               CALL SSYTRS( 'L', N, 1, AF, LDAF, IPIV, WORK, N,
     $                      INFO )
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
      IF( AINVNM .NE. 0.0 )
     $   SLA_SYRCOND = ( 1.0 / AINVNM )
*
      RETURN
*
*     End of SLA_SYRCOND
*
      END
