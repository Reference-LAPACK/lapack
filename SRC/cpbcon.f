*> \brief \b CPBCON
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CPBCON + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbcon.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbcon.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbcon.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,
*                          RWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, KD, LDAB, N
*       REAL               ANORM, RCOND
*       ..
*       .. Array Arguments ..
*       REAL               RWORK( * )
*       COMPLEX            AB( LDAB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CPBCON estimates the reciprocal of the condition number (in the
*> 1-norm) of a complex Hermitian positive definite band matrix using
*> the Cholesky factorization A = U**H*U or A = L*L**H computed by
*> CPBTRF.
*>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the
*> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangular factor stored in AB;
*>          = 'L':  Lower triangular factor stored in AB.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] KD
*> \verbatim
*>          KD is INTEGER
*>          The number of superdiagonals of the matrix A if UPLO = 'U',
*>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
*> \endverbatim
*>
*> \param[in] AB
*> \verbatim
*>          AB is COMPLEX array, dimension (LDAB,N)
*>          The triangular factor U or L from the Cholesky factorization
*>          A = U**H*U or A = L*L**H of the band matrix A, stored in the
*>          first KD+1 rows of the array.  The j-th column of U or L is
*>          stored in the j-th column of the array AB as follows:
*>          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
*>          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= KD+1.
*> \endverbatim
*>
*> \param[in] ANORM
*> \verbatim
*>          ANORM is REAL
*>          The 1-norm (or infinity-norm) of the Hermitian band matrix A.
*> \endverbatim
*>
*> \param[out] RCOND
*> \verbatim
*>          RCOND is REAL
*>          The reciprocal of the condition number of the matrix A,
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
*>          estimate of the 1-norm of inv(A) computed in this routine.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
*> \ingroup pbcon
*
*  =====================================================================
      SUBROUTINE CPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,
     $                   RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL               ANORM, RCOND
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            AB( LDAB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      CHARACTER          NORMIN
      INTEGER            IX, KASE
      REAL               AINVNM, SCALE, SCALEL, SCALEU, SMLNUM
      COMPLEX            ZDUM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICAMAX
      REAL               SLAMCH
      EXTERNAL           LSAME, ICAMAX, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACN2, CLATBS, CSRSCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, REAL
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPBCON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      SMLNUM = SLAMCH( 'Safe minimum' )
*
*     Estimate the 1-norm of the inverse.
*
      KASE = 0
      NORMIN = 'N'
   10 CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( UPPER ) THEN
*
*           Multiply by inv(U**H).
*
            CALL CLATBS( 'Upper', 'Conjugate transpose', 'Non-unit',
     $                   NORMIN, N, KD, AB, LDAB, WORK, SCALEL, RWORK,
     $                   INFO )
            NORMIN = 'Y'
*
*           Multiply by inv(U).
*
            CALL CLATBS( 'Upper', 'No transpose', 'Non-unit', NORMIN,
     $                   N,
     $                   KD, AB, LDAB, WORK, SCALEU, RWORK, INFO )
         ELSE
*
*           Multiply by inv(L).
*
            CALL CLATBS( 'Lower', 'No transpose', 'Non-unit', NORMIN,
     $                   N,
     $                   KD, AB, LDAB, WORK, SCALEL, RWORK, INFO )
            NORMIN = 'Y'
*
*           Multiply by inv(L**H).
*
            CALL CLATBS( 'Lower', 'Conjugate transpose', 'Non-unit',
     $                   NORMIN, N, KD, AB, LDAB, WORK, SCALEU, RWORK,
     $                   INFO )
         END IF
*
*        Multiply by 1/SCALE if doing so will not cause overflow.
*
         SCALE = SCALEL*SCALEU
         IF( SCALE.NE.ONE ) THEN
            IX = ICAMAX( N, WORK, 1 )
            IF( SCALE.LT.CABS1( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO )
     $         GO TO 20
            CALL CSRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
   20 CONTINUE
*
      RETURN
*
*     End of CPBCON
*
      END
