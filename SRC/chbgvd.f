*> \brief \b CHBGVD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CHBGVD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chbgvd.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chbgvd.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chbgvd.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CHBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W,
*                          Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK,
*                          LIWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, UPLO
*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LRWORK,
*      $                   LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       REAL               RWORK( * ), W( * )
*       COMPLEX            AB( LDAB, * ), BB( LDBB, * ), WORK( * ),
*      $                   Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CHBGVD computes all the eigenvalues, and optionally, the eigenvectors
*> of a complex generalized Hermitian-definite banded eigenproblem, of
*> the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian
*> and banded, and B is also positive definite.  If eigenvectors are
*> desired, it uses a divide and conquer algorithm.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
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
*> \param[in] KA
*> \verbatim
*>          KA is INTEGER
*>          The number of superdiagonals of the matrix A if UPLO = 'U',
*>          or the number of subdiagonals if UPLO = 'L'. KA >= 0.
*> \endverbatim
*>
*> \param[in] KB
*> \verbatim
*>          KB is INTEGER
*>          The number of superdiagonals of the matrix B if UPLO = 'U',
*>          or the number of subdiagonals if UPLO = 'L'. KB >= 0.
*> \endverbatim
*>
*> \param[in,out] AB
*> \verbatim
*>          AB is COMPLEX array, dimension (LDAB, N)
*>          On entry, the upper or lower triangle of the Hermitian band
*>          matrix A, stored in the first ka+1 rows of the array.  The
*>          j-th column of A is stored in the j-th column of the array AB
*>          as follows:
*>          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).
*>
*>          On exit, the contents of AB are destroyed.
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= KA+1.
*> \endverbatim
*>
*> \param[in,out] BB
*> \verbatim
*>          BB is COMPLEX array, dimension (LDBB, N)
*>          On entry, the upper or lower triangle of the Hermitian band
*>          matrix B, stored in the first kb+1 rows of the array.  The
*>          j-th column of B is stored in the j-th column of the array BB
*>          as follows:
*>          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;
*>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).
*>
*>          On exit, the factor S from the split Cholesky factorization
*>          B = S**H*S, as returned by CPBSTF.
*> \endverbatim
*>
*> \param[in] LDBB
*> \verbatim
*>          LDBB is INTEGER
*>          The leading dimension of the array BB.  LDBB >= KB+1.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is REAL array, dimension (N)
*>          If INFO = 0, the eigenvalues in ascending order.
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is COMPLEX array, dimension (LDZ, N)
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
*>          eigenvectors, with the i-th column of Z holding the
*>          eigenvector associated with W(i). The eigenvectors are
*>          normalized so that Z**H*B*Z = I.
*>          If JOBZ = 'N', then Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1, and if
*>          JOBZ = 'V', LDZ >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*>          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If N <= 1,               LWORK >= 1.
*>          If JOBZ = 'N' and N > 1, LWORK >= N.
*>          If JOBZ = 'V' and N > 1, LWORK >= 2*N**2.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal sizes of the WORK, RWORK and
*>          IWORK arrays, returns these values as the first entries of
*>          the WORK, RWORK and IWORK arrays, and no error message
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (MAX(1,LRWORK))
*>          On exit, if INFO=0, RWORK(1) returns the optimal LRWORK.
*> \endverbatim
*>
*> \param[in] LRWORK
*> \verbatim
*>          LRWORK is INTEGER
*>          The dimension of array RWORK.
*>          If N <= 1,               LRWORK >= 1.
*>          If JOBZ = 'N' and N > 1, LRWORK >= N.
*>          If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.
*>
*>          If LRWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal sizes of the WORK, RWORK
*>          and IWORK arrays, returns these values as the first entries
*>          of the WORK, RWORK and IWORK arrays, and no error message
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*>          On exit, if INFO=0, IWORK(1) returns the optimal LIWORK.
*> \endverbatim
*>
*> \param[in] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          The dimension of array IWORK.
*>          If JOBZ = 'N' or N <= 1, LIWORK >= 1.
*>          If JOBZ = 'V' and N > 1, LIWORK >= 3 + 5*N.
*>
*>          If LIWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal sizes of the WORK, RWORK
*>          and IWORK arrays, returns these values as the first entries
*>          of the WORK, RWORK and IWORK arrays, and no error message
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, and i is:
*>             <= N:  the algorithm failed to converge:
*>                    i off-diagonal elements of an intermediate
*>                    tridiagonal form did not converge to zero;
*>             > N:   if INFO = N + i, for 1 <= i <= N, then CPBSTF
*>                    returned INFO = i: B is not positive definite.
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
*> \ingroup hbgvd
*
*> \par Contributors:
*  ==================
*>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
*
*  =====================================================================
      SUBROUTINE CHBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB,
     $                   W,
     $                   Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK,
     $                   LIWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LRWORK,
     $                   LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL               RWORK( * ), W( * )
      COMPLEX            AB( LDAB, * ), BB( LDBB, * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            CONE, CZERO
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ),
     $                   CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, WANTZ
      CHARACTER          VECT
      INTEGER            IINFO, INDE, INDWK2, INDWRK, LIWMIN, LLRWK,
     $                   LLWK2, LRWMIN, LWMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SROUNDUP_LWORK
      EXTERNAL           LSAME, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSTERF, XERBLA, CGEMM, CHBGST, CHBTRD,
     $                   CLACPY,
     $                   CPBSTF, CSTEDC
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      INFO = 0
      IF( N.LE.1 ) THEN
         LWMIN = 1+N
         LRWMIN = 1+N
         LIWMIN = 1
      ELSE IF( WANTZ ) THEN
         LWMIN = 2*N**2
         LRWMIN = 1 + 5*N + 2*N**2
         LIWMIN = 3 + 5*N
      ELSE
         LWMIN = N
         LRWMIN = N
         LIWMIN = 1
      END IF
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KA.LT.0 ) THEN
         INFO = -4
      ELSE IF( KB.LT.0 .OR. KB.GT.KA ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KA+1 ) THEN
         INFO = -7
      ELSE IF( LDBB.LT.KB+1 ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         RWORK( 1 ) = REAL( LRWMIN )
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -14
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -16
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -18
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHBGVD', -INFO )
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
*     Form a split Cholesky factorization of B.
*
      CALL CPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem.
*
      INDE = 1
      INDWRK = INDE + N
      INDWK2 = 1 + N*N
      LLWK2 = LWORK - INDWK2 + 2
      LLRWK = LRWORK - INDWRK + 2
      CALL CHBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ,
     $             WORK, RWORK, IINFO )
*
*     Reduce Hermitian band matrix to tridiagonal form.
*
      IF( WANTZ ) THEN
         VECT = 'U'
      ELSE
         VECT = 'N'
      END IF
      CALL CHBTRD( VECT, UPLO, N, KA, AB, LDAB, W, RWORK( INDE ), Z,
     $             LDZ, WORK, IINFO )
*
*     For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEDC.
*
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL CSTEDC( 'I', N, W, RWORK( INDE ), WORK, N,
     $                WORK( INDWK2 ),
     $                LLWK2, RWORK( INDWRK ), LLRWK, IWORK, LIWORK,
     $                INFO )
         CALL CGEMM( 'N', 'N', N, N, N, CONE, Z, LDZ, WORK, N, CZERO,
     $               WORK( INDWK2 ), N )
         CALL CLACPY( 'A', N, N, WORK( INDWK2 ), N, Z, LDZ )
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      RWORK( 1 ) = REAL( LRWMIN )
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of CHBGVD
*
      END
