*> \brief \b ZTRTI2_MOD computes the inverse of a triangular matrix allowing for
*>    a zero diagonal (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZTRTI2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrti2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrti2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrti2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRTI2_MOD computes the inverse of a real upper or lower triangular
*> matrix. It allows for a zero diagonal. If a diagonal is 0, then
*> that column is set to 0
*>
*> This is the Level 2 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the matrix A is upper or lower triangular.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          Specifies whether or not the matrix A is unit triangular.
*>          = 'N':  Non-unit triangular
*>          = 'U':  Unit triangular
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
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the triangular matrix A.  If UPLO = 'U', the
*>          leading n by n upper triangular part of the array A contains
*>          the upper triangular matrix, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading n by n lower triangular part of the array A contains
*>          the lower triangular matrix, and the strictly upper
*>          triangular part of A is not referenced.  If DIAG = 'U', the
*>          diagonal elements of A are also not referenced and are
*>          assumed to be 1.
*>
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
*>          < 0: if INFO = -k, the k-th argument had an illegal value
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
*> \ingroup trti2
*
*  =====================================================================
      SUBROUTINE ZTRTI2_MOD( UPLO, DIAG, N, A, LDA, INFO )
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
      COMPLEX*16   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UNIT, UPPER
      INTEGER            J
      COMPLEX*16   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      UNIT = LSAME( DIAG, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UNIT .AND. .NOT.LSAME( DIAG, 'N' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO J = N, 1, -1
            ! Determine if we have a 0 diagonal.
            ! If not, we compute this column as normal
            ! Otherwise we set the current column to 0.
            ! Then all 0s should propogate to the right in the
            ! same column
            IF( UNIT ) THEN
               CALL ZSCAL(J-1, -ONE, A(1,J), 1)
            ELSE
               ! This is the only case where we need to check
               ! for zero diagonals
               IF( A(J,J).EQ.ZERO ) THEN
                  ! Set this column to 0
                  CALL ZLASET('a', J, 1, ZERO, ZERO, A(1,J), LDA)
               ELSE
                  ! In this case we must manually invert the diagonal
                  ! element then propogate
                  A(J,J) = ONE / A(J,J)
                  CALL ZSCAL(J-1, -A(J,J), A(1,J), 1)
               END IF
            END IF
            ! Regardless of what we have done above, we want to replace
            ! the current column of A with the solution to
            ! A(1:j-1,1:j-1) x = A(1:j-1,j) and store this
            ! inside A(1:j-1,j) but potentially having some
            ! zeros on the diagonal of A, so we call our modified routine
            CALL ZTRSM_MOD('Left', 'Upper', 'No Transpose', DIAG,
     $         J-1, 1, ONE, A, LDA, A(1,J), LDA)
         END DO
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO J = 1, N
            IF( UNIT ) THEN
               CALL ZSCAL(N-J, -ONE, A(J+1,J), 1)
            ELSE
               IF( A(J,J).EQ.ZERO ) THEN
                  ! Set this column to 0
                  CALL ZLASET('a', N-J, 1, ZERO, ZERO, A(J+1,J), LDA)
                  !call zlaset('a', 1, j-1, zero, zero, a(j+1,1), lda)
               ELSE
                  ! In this case we must manually invert the diagonal
                  ! element then propogate
                  A(J,J) = ONE / A(J,J)
                  CALL ZSCAL(N-J, -A(J,J), A(J+1,J), 1)
               END IF
            END IF
            ! Regardless of what we have done above, we want to replace
            ! the current column of A with the solution to
            ! A(j+1:n,j+1:n) x = A(j+1:n,j) and store this
            ! inside A(j+1:n,j) but potentially having some
            ! zeros on the diagonal of A, so we call our modified routine
            CALL ZTRSM_MOD('Left', 'Lower', 'No Transpose', DIAG,
     $         N-J, 1, ONE, A(J+1,J+1), LDA, A(J+1,J), LDA)
         END DO
      END IF
*
      RETURN
*
*     End of ZTRTI2_MOD
*
      END
