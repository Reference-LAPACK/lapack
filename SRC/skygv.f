*> \brief \b SKYGV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SKYGV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/skygv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/skygv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/skygv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
*                         LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKYGV computes all the eigenvalues, and optionally, the eigenvectors
*> of a real generalized skew-symmetric-definite eigenproblem, of the form
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
*> Here A is assumed to be skew-symmetric and B is assumed to be symmetric
*> positive definite.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ITYPE
*> \verbatim
*>          ITYPE is INTEGER
*>          Specifies the problem type to be solved:
*>          = 1:  A*x = (lambda)*B*x
*>          = 2:  A*B*x = (lambda)*x
*>          = 3:  B*A*x = (lambda)*x
*> \endverbatim
*>
*> \param[in] JOBZ
*> \verbatim
*>          JOBZ is CHARACTER*1
*>          = 'N':  Compute eigenvalues only;
*>          = 'V':  Compute eigenvalues and eigenvectors.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangles of A and B are stored;
*>          = 'L':  Lower triangles of A and B are stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA, N)
*>          On entry, the skew-symmetric matrix A.  If UPLO = 'U',
*>          the strictly N-by-N upper triangular part of A contains the
*>          upper triangular part of the matrix A.  If UPLO = 'L',
*>          the strictly N-by-N lower triangular part of A contains
*>          the lower triangular part of the matrix A.
*>
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*>          matrix Z, which leads to the block diagonal form in W.
*>          The matrix are normalized as follows:
*>          if ITYPE = 1 or 2, Z**T*B*Z = I;
*>          if ITYPE = 3, Z**T*inv(B)*Z = I.
*>          The eigenvectors of the matrix can be evaluated directly.
*>          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
*>          or the lower triangle (if UPLO='L') of A, including the
*>          diagonal, is destroyed.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB, N)
*>          On entry, the symmetric positive definite matrix B.
*>          If UPLO = 'U', the leading N-by-N upper triangular part of B
*>          contains the upper triangular part of the matrix B.
*>          If UPLO = 'L', the leading N-by-N lower triangular part of B
*>          contains the lower triangular part of the matrix B.
*>
*>          On exit, if INFO <= N, the part of B containing the matrix is
*>          overwritten by the triangular factor U or L from the Cholesky
*>          factorization B = U**T*U or B = L*L**T.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is REAL array, dimension (N)
*>          If INFO = 0, the (N-1) lower subdiagonal elements of the
*>          block diagonal matrix at front, and zero at last.
*>		    The matrix consists of 2-by-2 skew-symmetric blocks, and zeros.
*>          The values in W, which represent blocks, are always
*>          positive, and sorted in descending order.
*>          The eigenvalues of each blocks can be evaluated directly.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of the array WORK.  LWORK >= max(1,3*N-1).
*>          For optimal efficiency, LWORK >= (NB+2)*N,
*>          where NB is the blocksize for SSYTRD returned by ILAENV.
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
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  SPOTRF or SKYEV returned an error code:
*>             <= N:  if INFO = i, SKYEV failed to converge;
*>                    i off-diagonal elements of an intermediate
*>                    tridiagonal form did not converge to zero;
*>             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
*>                    minor of order i of B is not positive definite.
*>                    The factorization of B could not be completed and
*>                    no eigenvalues or eigenvectors were computed.
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
*> \ingroup kygv
*
*  =====================================================================
      SUBROUTINE SKYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
     $                  LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, WANTZ
      CHARACTER          TRANS
      INTEGER            LWKMIN, LWKOPT, NB, NEIG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SROUNDUP_LWORK
      EXTERNAL           ILAENV, LSAME, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SPOTRF, SKYEV, SKYGST, STRMM, STRSM,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
         LWKMIN = MAX( 1, 2*N - 1 )
         NB = ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( LWKMIN, ( NB + 1 )*N )
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -11
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SKYGV ', -INFO )
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
*     Form a Cholesky factorization of B.
*
      CALL SPOTRF( UPLO, N, B, LDB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem and solve.
*
      CALL SKYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      CALL SKYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
      IF( WANTZ ) THEN
*
*        Backtransform eigenvectors to the original problem.
*
         NEIG = N
         IF( INFO.GT.0 )
     $      NEIG = INFO - 1
         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) THEN
*
*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
*           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
*
            IF( UPPER ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'T'
            END IF
*
            CALL STRSM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG,
     $                  ONE,
     $                  B, LDB, A, LDA )
*
         ELSE IF( ITYPE.EQ.3 ) THEN
*
*           For B*A*x=(lambda)*x;
*           backtransform eigenvectors: x = L*y or U**T*y
*
            IF( UPPER ) THEN
               TRANS = 'T'
            ELSE
               TRANS = 'N'
            END IF
*
            CALL STRMM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG,
     $                  ONE,
     $                  B, LDB, A, LDA )
         END IF
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN
*
*     End of SKYGV
*
      END
