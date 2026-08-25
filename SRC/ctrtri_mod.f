*> \brief \b CTRTRI_MOD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CTRTRI_MOD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/DTRTRI_MOD.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/DTRTRI_MOD.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/DTRTRI_MOD.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE CTRTRI_MOD( UPLO, DIAG, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CTRTRI_MOD computes the inverse of a real upper or lower triangular
*> matrix A.
*>
*> This is the Level 3 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  A is upper triangular;
*>          = 'L':  A is lower triangular.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          = 'N':  A is non-unit triangular;
*>          = 'U':  A is unit triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the triangular matrix A.  If UPLO = 'U', the
*>          leading N-by-N upper triangular part of the array A contains
*>          the upper triangular matrix, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of the array A contains
*>          the lower triangular matrix, and the strictly upper
*>          triangular part of A is not referenced.  If DIAG = 'U', the
*>          diagonal elements of A are also not referenced and are
*>          assumed to be 1.
*>          On exit, the (triangular) inverse of the original matrix, in
*>          the same storage format.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*>               matrix is singular and its inverse can not be computed.
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
*> \ingroup trtri
*
*  =====================================================================
      RECURSIVE SUBROUTINE CTRTRI_MOD( UPLO, DIAG, N, A, LDA, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX   ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            K, J, JB, NB, NX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRI_MOD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine when to crossover to level2 version.
*
      !NX = ILAENV( 3, 'DTRTRI_MOD', UPLO // DIAG, N, -1, -1, -1 )
      NX = 1

      IF( N.LE.NX ) THEN
         CALL CTRTI2_MOD(UPLO, DIAG, N, A, LDA, INFO)
         RETURN
      END IF

      K = N/2

      IF( UPPER ) THEN
         ! Compute X_{22}
         CALL CTRTRI_MOD('Upper', DIAG, N-K, A(K+1,K+1), LDA, INFO)
         ! Propogate to A_{12}
         CALL CTRMM('Right', 'Upper', 'No Transpose', DIAG, K, N-K,
     $      -ONE, A(K+1,K+1), LDA, A(1,K+1), LDA)
         ! Solve for X_{12}
         CALL CTRSM_MOD('Left', 'Upper', 'No Transpose', DIAG,
     $      K, N-K, ONE, A, LDA, A(1,K+1), LDA)
         ! Solve for X_{11}
         CALL CTRTRI_MOD('Upper', DIAG, K, A, LDA, INFO)
      ELSE ! A is lower
         ! Compute X_{11}
         CALL CTRTRI_MOD('Lower', DIAG, K, A, LDA, INFO)
         ! Propogate to A_{21}
         CALL CTRMM('Right', 'Lower', 'No Transpose', DIAG, N-K, K,
     $      -ONE, A, LDA, A(K+1,1), LDA)
         ! Solve for X_{21}
         CALL CTRSM_MOD('Left', 'Lower', 'No Transpose', DIAG,
     $      N-K, K, ONE, A(K+1,K+1), LDA, A(K+1,1), LDA)
         ! Compute X_{22}
         CALL CTRTRI_MOD('Lower', DIAG, N-K, A(K+1,K+1), LDA, INFO)
      END IF

      END SUBROUTINE
