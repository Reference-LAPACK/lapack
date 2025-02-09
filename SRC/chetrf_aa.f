*> \brief \b CHETRF_AA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CHETRF_AA + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrf_aa.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrf_aa.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrf_aa.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CHETRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER    UPLO
*       INTEGER      N, LDA, LWORK, INFO
*       ..
*       .. Array Arguments ..
*       INTEGER      IPIV( * )
*       COMPLEX      A( LDA, * ), WORK( * )
*       ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CHETRF_AA computes the factorization of a complex hermitian matrix A
*> using the Aasen's algorithm.  The form of the factorization is
*>
*>    A = U**H*T*U  or  A = L*T*L**H
*>
*> where U (or L) is a product of permutation and unit upper (lower)
*> triangular matrices, and T is a hermitian tridiagonal matrix.
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
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the hermitian matrix A.  If UPLO = 'U', the leading
*>          N-by-N upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the strictly lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          leading N-by-N lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the strictly upper
*>          triangular part of A is not referenced.
*>
*>          On exit, the tridiagonal matrix is stored in the diagonals
*>          and the subdiagonals of A just below (or above) the diagonals,
*>          and L is stored below (or above) the subdiagonals, when UPLO
*>          is 'L' (or 'U').
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
*>          On exit, it contains the details of the interchanges, i.e.,
*>          the row and column k of A were interchanged with the
*>          row and column IPIV(k).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of WORK.
*>          LWORK >= 1, if N <= 1, and LWORK >= 2*N, otherwise.
*>          For optimum performance LWORK >= N*(1+NB), where NB is
*>          the optimal blocksize, returned by ILAENV.
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
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
*> \ingroup hetrf_aa
*
*  =====================================================================
      SUBROUTINE CHETRF_AA( UPLO, N, A, LDA, IPIV,
     $                      WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER    UPLO
      INTEGER      N, LDA, LWORK, INFO
*     ..
*     .. Array Arguments ..
      INTEGER      IPIV( * )
      COMPLEX      A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*     .. Parameters ..
      COMPLEX      ZERO, ONE
      PARAMETER    ( ZERO = (0.0E+0, 0.0E+0), ONE = (1.0E+0, 0.0E+0) )
*
*     .. Local Scalars ..
      LOGICAL      LQUERY, UPPER
      INTEGER      J, LWKMIN, LWKOPT
      INTEGER      NB, MJ, NJ, K1, K2, J1, J2, J3, JB
      COMPLEX      ALPHA
*     ..
*     .. External Functions ..
      LOGICAL      LSAME
      INTEGER      ILAENV
      REAL         SROUNDUP_LWORK
      EXTERNAL     LSAME, ILAENV, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL     CLAHEF_AA, CGEMM, CCOPY, CSWAP, CSCAL,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC    REAL, CONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Determine the block size
*
      NB = ILAENV( 1, 'CHETRF_AA', UPLO, N, -1, -1, -1 )
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LE.1 ) THEN
         LWKMIN = 1
         LWKOPT = 1
      ELSE
         LWKMIN = 2*N
         LWKOPT = (NB+1)*N
      END IF
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHETRF_AA', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return
*
      IF( N.EQ.0 ) THEN
          RETURN
      ENDIF
      IPIV( 1 ) = 1
      IF( N.EQ.1 ) THEN
         A( 1, 1 ) = REAL( A( 1, 1 ) )
         RETURN
      END IF
*
*     Adjust block size based on the workspace size
*
      IF( LWORK.LT.((1+NB)*N) ) THEN
         NB = ( LWORK-N ) / N
      END IF
*
      IF( UPPER ) THEN
*
*        .....................................................
*        Factorize A as U**H*D*U using the upper triangle of A
*        .....................................................
*
*        copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N))
*
         CALL CCOPY( N, A( 1, 1 ), LDA, WORK( 1 ), 1 )
*
*        J is the main loop index, increasing from 1 to N in steps of
*        JB, where JB is the number of columns factorized by CLAHEF;
*        JB is either NB, or N-J+1 for the last block
*
         J = 0
 10      CONTINUE
         IF( J.GE.N )
     $      GO TO 20
*
*        each step of the main loop
*         J is the last column of the previous panel
*         J1 is the first column of the current panel
*         K1 identifies if the previous column of the panel has been
*          explicitly stored, e.g., K1=1 for the first panel, and
*          K1=0 for the rest
*
         J1 = J + 1
         JB = MIN( N-J1+1, NB )
         K1 = MAX(1, J)-J
*
*        Panel factorization
*
         CALL CLAHEF_AA( UPLO, 2-K1, N-J, JB,
     $                      A( MAX(1, J), J+1 ), LDA,
     $                      IPIV( J+1 ), WORK, N, WORK( N*NB+1 ) )
*
*        Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)
*
         DO J2 = J+2, MIN(N, J+JB+1)
            IPIV( J2 ) = IPIV( J2 ) + J
            IF( (J2.NE.IPIV(J2)) .AND. ((J1-K1).GT.2) ) THEN
               CALL CSWAP( J1-K1-2, A( 1, J2 ), 1,
     $                              A( 1, IPIV(J2) ), 1 )
            END IF
         END DO
         J = J + JB
*
*        Trailing submatrix update, where
*         the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and
*         WORK stores the current block of the auxiriarly matrix H
*
         IF( J.LT.N ) THEN
*
*          if the first panel and JB=1 (NB=1), then nothing to do
*
            IF( J1.GT.1 .OR. JB.GT.1 ) THEN
*
*              Merge rank-1 update with BLAS-3 update
*
               ALPHA = CONJG( A( J, J+1 ) )
               A( J, J+1 ) = ONE
               CALL CCOPY( N-J, A( J-1, J+1 ), LDA,
     $                          WORK( (J+1-J1+1)+JB*N ), 1 )
               CALL CSCAL( N-J, ALPHA, WORK( (J+1-J1+1)+JB*N ), 1 )
*
*              K1 identifies if the previous column of the panel has been
*               explicitly stored, e.g., K1=0 and K2=1 for the first panel,
*               and K1=1 and K2=0 for the rest
*
               IF( J1.GT.1 ) THEN
*
*                 Not first panel
*
                  K2 = 1
               ELSE
*
*                 First panel
*
                  K2 = 0
*
*                 First update skips the first column
*
                  JB = JB - 1
               END IF
*
               DO J2 = J+1, N, NB
                  NJ = MIN( NB, N-J2+1 )
*
*                 Update (J2, J2) diagonal block with CGEMV
*
                  J3 = J2
                  DO MJ = NJ-1, 1, -1
                     CALL CGEMM( 'Conjugate transpose', 'Transpose',
     $                            1, MJ, JB+1,
     $                           -ONE, A( J1-K2, J3 ), LDA,
     $                                 WORK( (J3-J1+1)+K1*N ), N,
     $                            ONE, A( J3, J3 ), LDA )
                     J3 = J3 + 1
                  END DO
*
*                 Update off-diagonal block of J2-th block row with CGEMM
*
                  CALL CGEMM( 'Conjugate transpose', 'Transpose',
     $                        NJ, N-J3+1, JB+1,
     $                       -ONE, A( J1-K2, J2 ), LDA,
     $                             WORK( (J3-J1+1)+K1*N ), N,
     $                        ONE, A( J2, J3 ), LDA )
               END DO
*
*              Recover T( J, J+1 )
*
               A( J, J+1 ) = CONJG( ALPHA )
            END IF
*
*           WORK(J+1, 1) stores H(J+1, 1)
*
            CALL CCOPY( N-J, A( J+1, J+1 ), LDA, WORK( 1 ), 1 )
         END IF
         GO TO 10
      ELSE
*
*        .....................................................
*        Factorize A as L*D*L**H using the lower triangle of A
*        .....................................................
*
*        copy first column A(1:N, 1) into H(1:N, 1)
*         (stored in WORK(1:N))
*
         CALL CCOPY( N, A( 1, 1 ), 1, WORK( 1 ), 1 )
*
*        J is the main loop index, increasing from 1 to N in steps of
*        JB, where JB is the number of columns factorized by CLAHEF;
*        JB is either NB, or N-J+1 for the last block
*
         J = 0
 11      CONTINUE
         IF( J.GE.N )
     $      GO TO 20
*
*        each step of the main loop
*         J is the last column of the previous panel
*         J1 is the first column of the current panel
*         K1 identifies if the previous column of the panel has been
*          explicitly stored, e.g., K1=1 for the first panel, and
*          K1=0 for the rest
*
         J1 = J+1
         JB = MIN( N-J1+1, NB )
         K1 = MAX(1, J)-J
*
*        Panel factorization
*
         CALL CLAHEF_AA( UPLO, 2-K1, N-J, JB,
     $                      A( J+1, MAX(1, J) ), LDA,
     $                      IPIV( J+1 ), WORK, N, WORK( N*NB+1 ) )
*
*        Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)
*
         DO J2 = J+2, MIN(N, J+JB+1)
            IPIV( J2 ) = IPIV( J2 ) + J
            IF( (J2.NE.IPIV(J2)) .AND. ((J1-K1).GT.2) ) THEN
               CALL CSWAP( J1-K1-2, A( J2, 1 ), LDA,
     $                              A( IPIV(J2), 1 ), LDA )
            END IF
         END DO
         J = J + JB
*
*        Trailing submatrix update, where
*          A(J2+1, J1-1) stores L(J2+1, J1) and
*          WORK(J2+1, 1) stores H(J2+1, 1)
*
         IF( J.LT.N ) THEN
*
*          if the first panel and JB=1 (NB=1), then nothing to do
*
            IF( J1.GT.1 .OR. JB.GT.1 ) THEN
*
*              Merge rank-1 update with BLAS-3 update
*
               ALPHA = CONJG( A( J+1, J ) )
               A( J+1, J ) = ONE
               CALL CCOPY( N-J, A( J+1, J-1 ), 1,
     $                          WORK( (J+1-J1+1)+JB*N ), 1 )
               CALL CSCAL( N-J, ALPHA, WORK( (J+1-J1+1)+JB*N ), 1 )
*
*              K1 identifies if the previous column of the panel has been
*               explicitly stored, e.g., K1=0 and K2=1 for the first panel,
*               and K1=1 and K2=0 for the rest
*
               IF( J1.GT.1 ) THEN
*
*                 Not first panel
*
                  K2 = 1
               ELSE
*
*                 First panel
*
                  K2 = 0
*
*                 First update skips the first column
*
                  JB = JB - 1
               END IF
*
               DO J2 = J+1, N, NB
                  NJ = MIN( NB, N-J2+1 )
*
*                 Update (J2, J2) diagonal block with CGEMV
*
                  J3 = J2
                  DO MJ = NJ-1, 1, -1
                     CALL CGEMM( 'No transpose',
     $                           'Conjugate transpose',
     $                           MJ, 1, JB+1,
     $                          -ONE, WORK( (J3-J1+1)+K1*N ), N,
     $                                A( J3, J1-K2 ), LDA,
     $                           ONE, A( J3, J3 ), LDA )
                     J3 = J3 + 1
                  END DO
*
*                 Update off-diagonal block of J2-th block column with CGEMM
*
                  CALL CGEMM( 'No transpose', 'Conjugate transpose',
     $                        N-J3+1, NJ, JB+1,
     $                       -ONE, WORK( (J3-J1+1)+K1*N ), N,
     $                             A( J2, J1-K2 ), LDA,
     $                        ONE, A( J3, J2 ), LDA )
               END DO
*
*              Recover T( J+1, J )
*
               A( J+1, J ) = CONJG( ALPHA )
            END IF
*
*           WORK(J+1, 1) stores H(J+1, 1)
*
            CALL CCOPY( N-J, A( J+1, J+1 ), 1, WORK( 1 ), 1 )
         END IF
         GO TO 11
      END IF
*
   20 CONTINUE
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      RETURN
*
*     End of CHETRF_AA
*
      END
