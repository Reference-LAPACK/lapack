*> \brief \b ZTPTTR copies a triangular matrix from the standard packed format (TP) to the standard full format (TR).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZTPTTR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpttr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpttr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpttr.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTPTTR( UPLO, N, AP, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, N, LDA
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), AP( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTPTTR copies a triangular matrix A from standard packed format (TP)
*> to standard full format (TR).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  A is upper triangular.
*>          = 'L':  A is lower triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A. N >= 0.
*> \endverbatim
*>
*> \param[in] AP
*> \verbatim
*>          AP is COMPLEX*16 array, dimension ( N*(N+1)/2 ),
*>          On entry, the upper or lower triangular matrix A, packed
*>          columnwise in a linear array. The j-th column of A is stored
*>          in the array AP as follows:
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, N )
*>          On exit, the triangular matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
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
*> \ingroup tpttr
*
*  =====================================================================
      SUBROUTINE ZTPTTR( UPLO, N, AP, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N, LDA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), AP( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER
      INTEGER            I, J, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTPTTR', -INFO )
         RETURN
      END IF
*
      IF( LOWER ) THEN
         K = 0
         DO J = 1, N
            DO I = J, N
               K = K + 1
               A( I, J ) = AP( K )
            END DO
         END DO
      ELSE
         K = 0
         DO J = 1, N
            DO I = 1, J
               K = K + 1
               A( I, J ) = AP( K )
            END DO
         END DO
      END IF
*
*
      RETURN
*
*     End of ZTPTTR
*
      END
