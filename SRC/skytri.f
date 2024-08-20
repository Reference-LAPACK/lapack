*> \brief \b SKYTRI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SKYTRI + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/skytri.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/skytri.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/skytri.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL               A( LDA, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKYTRI computes the inverse of a real skew-symmetric indefinite matrix
*> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
*> SSYTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
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
*>          On entry, the block diagonal matrix D and the multipliers
*>          used to obtain the factor U or L as computed by SKYTRF.
*>
*>          On exit, if INFO = 0, the (skew-symmetric) inverse of the original
*>          matrix.  If UPLO = 'U', the upper triangular part of the
*>          inverse is formed and the part of A below the diagonal is not
*>          referenced; if UPLO = 'L' the lower triangular part of the
*>          inverse is formed and the part of A above the diagonal is
*>          not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by SSYTRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (N)
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
*> \ingroup kytri
*
*  =====================================================================
      SUBROUTINE SKYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), WORK( * )
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
      INTEGER            K, KP, KSTEP
      REAL               TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT
      EXTERNAL           LSAME, SDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SSWAP, SKYMV, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
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
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SKYTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. MOD(N,2).NE.0 )
     $   RETURN
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         DO 10 INFO = N, 2, -2
            IF( A( INFO - 1, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         DO 20 INFO = 1, N-1, 2
            IF( A( INFO + 1, INFO ).EQ.ZERO )
     $         RETURN
   20    CONTINUE
      END IF
      INFO = 0
*
      IF( UPPER ) THEN
*
*        Compute inv(A) from the factorization A = U*D*U**T.
*
*        K is the main loop index, increasing from 1 to N in steps of 2
*
         K = 1
   30    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GE.N )
     $      GO TO 40
*
*
*        2 x 2 diagonal block
*
*        Invert the diagonal block.
*
         A( K, K+1 ) = -ONE / A( K, K+1 )
*
*        Compute columns K and K+1 of the inverse.
*
         IF( K.GT.1 ) THEN
            CALL SCOPY( K-1, A( 1, K ), 1, WORK, 1 )
            CALL SKYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                  A( 1, K ), 1 )
            A( K, K+1 ) = A( K, K+1 ) +
     $                    SDOT( K-1, A( 1, K ), 1, A( 1, K+1 ),
     $                     1 )
            CALL SCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
            CALL SKYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                  A( 1, K+1 ), 1 )
         END IF
         KSTEP = 2
*
         KP = IPIV( K+1 )
*
*        Interchange rows and columns K and KP in the leading
*        submatrix A(1:k+1,1:k+1)
*
         IF( KP.GT.0 ) THEN
            CALL SSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            CALL SSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
            CALL SSCAL( K-KP, -ONE, A( KP, K ), 1)
            CALL SSCAL( K-KP-1, -ONE, A( KP, KP+1 ), LDA )
            TEMP = A( K, K+1 )
            A( K, K+1 ) = A( KP, K+1 )
            A( KP, K+1 ) = TEMP
         ELSEIF( KP.LT.0 ) THEN
            KP = -KP
            CALL SSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            CALL SSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
            CALL SSCAL( K-KP, -ONE, A( KP, K ), 1)
            CALL SSCAL( K-KP-1, -ONE, A( KP, KP+1 ), LDA )
            TEMP = A( K, K+1 )
            A( K, K+1 ) = A( KP, K+1 )
            A( KP, K+1 ) = TEMP
            CALL SSWAP( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
            A( K, K+1 ) = -A( K, K+1 )
         END IF
*
         K = K + KSTEP
         GO TO 30
   40    CONTINUE
*
      ELSE
*
*        Compute inv(A) from the factorization A = L*D*L**T.
*
*        K is the main loop index, increasing from 1 to N in steps of 2
*
         K = N
   50    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LE.1 )
     $      GO TO 60
*
*
*        2 x 2 diagonal block
*
*        Invert the diagonal block.
*
         A( K, K-1 ) = -ONE / A( K, K-1 )
*
*        Compute columns K-1 and K of the inverse.
*
         IF( K.LT.N ) THEN
            CALL SCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
            CALL SKYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK,
     $                  1, ZERO, A( K+1, K ), 1 )
            A( K, K-1 ) = A( K, K-1 ) +
     $                    SDOT( N-K, A( K+1, K ), 1, A( K+1, K-1 ),
     $                    1 )
            CALL SCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
            CALL SKYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK,
     $                  1, ZERO, A( K+1, K-1 ), 1 )
         END IF
         KSTEP = 2
*
         KP = IPIV( K-1 )
*
*        Interchange rows and columns K and KP in the trailing
*        submatrix A(k-1:n,k-1:n)
*
         IF( KP.GT.0 ) THEN
            IF( KP.LT.N )
     $         CALL SSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            CALL SSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
            CALL SSCAL( KP-K, -ONE, A( K+1, K ), 1)
            CALL SSCAL( KP-K-1, -ONE, A( KP, K+1 ), LDA )
            TEMP = A( K, K-1 )
            A( K, K-1 ) = A( KP, K-1 )
            A( KP, K-1 ) = TEMP
         ELSEIF( KP.LT.0 ) THEN
            KP = -KP
            IF( KP.LT.N )
     $         CALL SSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            CALL SSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
            CALL SSCAL( KP-K, -ONE, A( K+1, K ), 1)
            CALL SSCAL( KP-K-1, -ONE, A( KP, K+1 ), LDA )
            TEMP = A( K, K-1 )
            A( K, K-1 ) = A( KP, K-1 )
            A( KP, K-1 ) = TEMP
            CALL SSWAP( N-K, A( K+1, K ), 1, A( K+1, K-1 ), 1 )
            A( K, K-1 ) = -A( K, K-1 )
         END IF
*
         K = K - KSTEP
         GO TO 50
   60    CONTINUE
      END IF
*
      RETURN
*
*     End of SKYTRI
*
      END
