*> \brief \b DKYTRF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DKYTRF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/skytrf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/skytrf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/skytrf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DKYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DKYTRF computes the factorization of a real skew-symmetric matrix A using
*> the Bunch partial pivoting method.  The form of the
*> factorization is
*>
*>    A = U**T*D*U  or  A = L*D*L**T
*>
*> where U (or L) is a product of permutation and unit upper (lower)
*> triangular matrices, and D is skew-symmetric and block diagonal with
*> 1-by-1 and 2-by-2 diagonal blocks.  All 2-by-2 diagonal blocks are
*> nonsingular and all 1-by-1 diagonal blocks are 0. If N is odd, there
*> is at least one 1-by-1 diagonal block.
*>
*> This is the blocked version of the algorithm, calling Level 3 BLAS.
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
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the skew-symmetric matrix A.  If UPLO = 'U', the
*>          strictly upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the leading N-by-N lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          strictly lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the leading N-by-N upper
*>          triangular part of A is not referenced.
*>
*>          On exit, the block diagonal matrix D and the multipliers used
*>          to obtain the factor U or L (see below for further details).
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
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of WORK.  LWORK >=1.  For best performance
*>          LWORK >= N*NB, where NB is the block size returned by ILAENV.
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
*> \ingroup kytrf
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  If UPLO = 'U', then A = U**T*D*U, where
*>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
*>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
*>  1 in steps of 2, and D is a block diagonal matrix with 2-by-2
*>  diagonal blocks D(k).  P(k) is a permutation matrix as defined by
*>  IPIV(k), and U(k) is a unit upper triangular matrix, such that if
*>  the diagonal block D(k) is of order 2, namely s = 2, then
*>
*>             (   I    v    0   )   k-s
*>     U(k) =  (   0    I    0   )   s
*>             (   0    0    I   )   n-k
*>                k-s   s   n-k
*>
*>  The strictly upper triangle of D(k) overwrites A(k-1,k), and v overwrites
*>  A(1:k-2,k-1:k).
*>
*>  If UPLO = 'L', then A = L*D*L**T, where
*>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
*>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
*>  n in steps of 2, and D is a block diagonal matrix with 2-by-2
*>  diagonal blocks D(k).  P(k) is a permutation matrix as defined by
*>  IPIV(k), and L(k) is a unit lower triangular matrix, such that if
*>  the diagonal block D(k) is of order 2, namely s = 2, then
*>
*>             (   I    0     0   )  k-1
*>     L(k) =  (   0    I     0   )  s
*>             (   0    v     I   )  n-k-s+1
*>                k-1   s  n-k-s+1
*>
*>  The strictly lower triangle of D(k) overwrites A(k+1,k), and v overwrites
*>  A(k+2:n,k:k+1).
*>
*>  Remind that if n is odd, A is always singular.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DKYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
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
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAKYF, DKYTF2, XERBLA
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
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size
*
         NB = ILAENV( 1, 'DKYTRF', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, N*NB )
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DKYTRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = LDWORK*NB
         IF( LWORK.LT.IWS ) THEN
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'DKYTRF', UPLO, N, -1, -1,
     $                   -1 ) )
         END IF
      ELSE
         IWS = 1
      END IF
      IF( NB.LT.NBMIN )
     $   NB = N
*
      IF( UPPER ) THEN
*
*        Factorize A as U**T*D*U using the upper triangle of A
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        KB, where KB is the number of columns factorized by DLAKYF;
*        KB is either NB or NB-1, or K for the last block
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop
*
         IF( K.LT.1 )
     $      GO TO 40
*
         IF( K.GT.NB ) THEN
*
*           Factorize columns k-kb+1:k of A and use blocked code to
*           update columns 1:k-kb
*
            CALL DLAKYF( UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK,
     $                   IINFO )
         ELSE
*
*           Use unblocked code to factorize columns 1:k of A
*
            CALL DKYTF2( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO
*
*        Decrease K and return to the start of the main loop
*
         K = K - KB
         GO TO 10
*
      ELSE
*
*        Factorize A as L*D*L**T using the lower triangle of A
*
*        K is the main loop index, increasing from 1 to N in steps of
*        KB, where KB is the number of columns factorized by DLAKYF;
*        KB is either NB or NB-1, or N-K+1 for the last block
*
         K = 1
   20    CONTINUE
*
*        If K > N, exit from loop
*
         IF( K.GT.N )
     $      GO TO 40
*
         IF( K.LE.N-NB ) THEN
*
*           Factorize columns k:k+kb-1 of A and use blocked code to
*           update columns k+kb:n
*
            CALL DLAKYF( UPLO, N-K+1, NB, KB, A( K, K ), LDA,
     $                   IPIV( K ),
     $                   WORK, LDWORK, IINFO )
         ELSE
*
*           Use unblocked code to factorize columns k:n of A
*
            CALL DKYTF2( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ),
     $                   IINFO )
            KB = N - K + 1
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + K - 1
*
*        Adjust IPIV
*
         DO 30 J = K, K + KB - 1
            IF( IPIV( J ).GT.0 ) THEN
               IPIV( J ) = IPIV( J ) + K - 1
            ELSEIF( IPIV( J ).LT.0 ) THEN
               IPIV( J ) = IPIV( J ) - K + 1
            END IF
   30    CONTINUE
*
*        Increase K and return to the start of the main loop
*
         K = K + KB
         GO TO 20
*
      END IF
*
   40 CONTINUE
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DKYTRF
*
      END
