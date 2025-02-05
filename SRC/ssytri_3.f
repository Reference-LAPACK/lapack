*> \brief \b SSYTRI_3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SSYTRI_3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri_3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri_3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri_3.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SSYTRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK,
*                            INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL               A( LDA, * ), E( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> SSYTRI_3 computes the inverse of a real symmetric indefinite
*> matrix A using the factorization computed by SSYTRF_RK or SSYTRF_BK:
*>
*>     A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
*>
*> where U (or L) is unit upper (or lower) triangular matrix,
*> U**T (or L**T) is the transpose of U (or L), P is a permutation
*> matrix, P**T is the transpose of P, and D is symmetric and block
*> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
*>
*> SSYTRI_3 sets the leading dimension of the workspace  before calling
*> SSYTRI_3X that actually computes the inverse.  This is the blocked
*> version of the algorithm, calling Level 3 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are
*>          stored as an upper or lower triangular matrix.
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
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
*>          A is REAL array, dimension (LDA,N)
*>          On entry, diagonal of the block diagonal matrix D and
*>          factors U or L as computed by SSYTRF_RK and SSYTRF_BK:
*>            a) ONLY diagonal elements of the symmetric block diagonal
*>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*>               (superdiagonal (or subdiagonal) elements of D
*>                should be provided on entry in array E), and
*>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
*>               If UPLO = 'L': factor L in the subdiagonal part of A.
*>
*>          On exit, if INFO = 0, the symmetric inverse of the original
*>          matrix.
*>             If UPLO = 'U': the upper triangular part of the inverse
*>             is formed and the part of A below the diagonal is not
*>             referenced;
*>             If UPLO = 'L': the lower triangular part of the inverse
*>             is formed and the part of A above the diagonal is not
*>             referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is REAL array, dimension (N)
*>          On entry, contains the superdiagonal (or subdiagonal)
*>          elements of the symmetric block diagonal matrix D
*>          with 1-by-1 or 2-by-2 diagonal blocks, where
*>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
*>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
*>
*>          NOTE: For 1-by-1 diagonal block D(k), where
*>          1 <= k <= N, the element E(k) is not referenced in both
*>          UPLO = 'U' or UPLO = 'L' cases.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by SSYTRF_RK or SSYTRF_BK.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK)).
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of WORK.
*>          If N = 0, LWORK >= 1, else LWORK >= (N+NB+1)*(NB+3).
*>
*>          If LWORK = -1, then a workspace query is assumed;
*>          the routine only calculates the optimal size of the optimal
*>          size of the WORK array, returns this value as the first
*>          entry of the WORK array, and no error message related to
*>          LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
*>               inverse could not be computed.
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
*> \ingroup hetri_3
*
*> \par Contributors:
*  ==================
*> \verbatim
*>
*>  November 2017,  Igor Kozachenko,
*>                  Computer Science Division,
*>                  University of California, Berkeley
*>
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE SSYTRI_3( UPLO, N, A, LDA, E, IPIV, WORK, LWORK,
     $                     INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            UPPER, LQUERY
      INTEGER            LWKOPT, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SROUNDUP_LWORK
      EXTERNAL           LSAME, ILAENV, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSYTRI_3X, XERBLA
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
      LQUERY = ( LWORK.EQ.-1 )
*
*     Determine the block size
*
      IF( N.EQ.0 ) THEN
         LWKOPT = 1
      ELSE
         NB = MAX( 1, ILAENV( 1, 'SSYTRI_3', UPLO, N, -1, -1, -1 ) )
         LWKOPT = ( N+NB+1 ) * ( NB+3 )
      END IF
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.LWKOPT .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTRI_3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      CALL SSYTRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO )
*
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
*
      RETURN
*
*     End of SSYTRI_3
*
      END
