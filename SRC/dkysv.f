*> \brief <b> DKYSV computes the solution to system of linear equations A * X = B for KY matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DKYSV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/skysv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/skysv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/skysv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DKYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
*                         LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DKYSV computes the solution to a real system of linear equations
*>    A * X = B,
*> where A is an N-by-N skew-symmetric matrix and X and B are N-by-NRHS
*> matrices.
*>
*> The partial pivoting method is used to factor A as
*>    A = U * D * U**T,  if UPLO = 'U', or
*>    A = L * D * L**T,  if UPLO = 'L',
*> where U (or L) is a product of permutation and unit upper (lower)
*> triangular matrices, and D is skew-symmetric and block diagonal with
*> 1-by-1 and 2-by-2 diagonal blocks. All 2-by-2 diagonal blocks are
*> nonsingular and all 1-by-1 diagonal blocks are 0. If N is odd, there
*> is at least one 1-by-1 diagonal block. The factored form of A is then
*> used to solve the system of equations A * X = B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored;
*>          = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of linear equations, i.e., the order of the
*>          matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the skew-symmetric matrix A.  If UPLO = 'U', the strictly
*>          upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the leading N-by-N lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          strictly lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the leading N-by-N upper
*>          triangular part of A is not referenced.
*>
*>          On exit, if INFO = 0, the block diagonal matrix D and the
*>          multipliers used to obtain the factor U or L from the
*>          factorization A = U*D*U**T or A = L*D*L**T as computed by
*>          DKYTRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges of D, as determined by DKYTRF.
*>
*>          The elements of array IPIV are combined in pair, and the first 
*>          (if UPLO = 'U') or the second (if UPLO = 'L') element in
*>          the pair always keeps the value 0.  If N is odd, the first
*>          (if UPLO = 'U') or the last (if UPLO = 'L') element of IPIV is
*>          0, which is the only element not in pair. So we only use the
*>          first (if UPLO = 'L') or the second (if UPLO = 'U') element in
*>          the pair to determine the interchanges.
*>
*>          If IPIV(k)
*>          = 0: there was no interchange.
*>          > 0: rows and columns k-1 and IPIV(k) were interchanged, if
*>               UPLO = 'U', and rows and columns k+1 and IPIV(k) were
*>               interchanged, if UPLO = 'L'.
*>          < 0: rows and columns k and k-1 were interchanged,
*>               then rows and columns k-1 and -IPIV(k) were interchanged, if
*>               UPLO = 'U', and rows and columns k and k+1 were interchanged,
*>               then rows and columns k+1 and -IPIV(k) were interchanged, if
*>               UPLO = 'L'.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*>          On entry, the N-by-NRHS right hand side matrix B.
*>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of WORK.  LWORK >= 1, and for best performance
*>          LWORK >= max(1,N*NB), where NB is the optimal blocksize for
*>          DKYTRF.
*>          for LWORK < N, TRS will be done with Level BLAS 2
*>          for LWORK >= N, TRS will be done with Level BLAS 3
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, D(i-1,i) (if UPLO = 'U') or D(i+1,i) (if UPLO = 'L')
*>               is exactly zero.  The factorization has been completed,
*>               but the block diagonal matrix D is exactly singular,
*>               so the solution could not be computed.
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
*> \ingroup kysv
*
*  =====================================================================
      SUBROUTINE DKYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
     $                  LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            LWKOPT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DKYTRF, DKYTRS, DKYTRS2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND.
     $    .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            CALL DKYTRF( UPLO, N, A, LDA, IPIV, WORK, -1, INFO )
            LWKOPT = INT( WORK( 1 ) )
         END IF
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DKYSV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Compute the factorization A = U*D*U**T or A = L*D*L**T.
*
      CALL DKYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         IF ( LWORK.LT.N ) THEN
*
*        Solve with TRS ( Use Level BLAS 2)
*
            CALL DKYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
         ELSE
*
*        Solve with TRS2 ( Use Level BLAS 3)
*
            CALL DKYTRS2( UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,INFO )
*
         END IF
*
      END IF
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DKYSV
*
      END
