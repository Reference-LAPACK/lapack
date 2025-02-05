*> \brief \b SPFTRF
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SPFTRF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spftrf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spftrf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spftrf.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SPFTRF( TRANSR, UPLO, N, A, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANSR, UPLO
*       INTEGER            N, INFO
*       ..
*       .. Array Arguments ..
*       REAL               A( 0: * )
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SPFTRF computes the Cholesky factorization of a real symmetric
*> positive definite matrix A.
*>
*> The factorization has the form
*>    A = U**T * U,  if UPLO = 'U', or
*>    A = L  * L**T,  if UPLO = 'L',
*> where U is an upper triangular matrix and L is lower triangular.
*>
*> This is the block version of the algorithm, calling Level 3 BLAS.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANSR
*> \verbatim
*>          TRANSR is CHARACTER*1
*>          = 'N':  The Normal TRANSR of RFP A is stored;
*>          = 'T':  The Transpose TRANSR of RFP A is stored.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of RFP A is stored;
*>          = 'L':  Lower triangle of RFP A is stored.
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
*>          A is REAL array, dimension ( N*(N+1)/2 );
*>          On entry, the symmetric matrix A in RFP format. RFP format is
*>          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N'
*>          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is
*>          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'T' then RFP is
*>          the transpose of RFP A as defined when
*>          TRANSR = 'N'. The contents of RFP A are defined by UPLO as
*>          follows: If UPLO = 'U' the RFP A contains the NT elements of
*>          upper packed A. If UPLO = 'L' the RFP A contains the elements
*>          of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR =
*>          'T'. When TRANSR is 'N' the LDA is N+1 when N is even and N
*>          is odd. See the Note below for more details.
*>
*>          On exit, if INFO = 0, the factor U or L from the Cholesky
*>          factorization RFP A = U**T*U or RFP A = L*L**T.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the leading principal minor of order i
*>                is not positive, and the factorization could not be
*>                completed.
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
*> \ingroup pftrf
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  We first consider Rectangular Full Packed (RFP) Format when N is
*>  even. We give an example where N = 6.
*>
*>      AP is Upper             AP is Lower
*>
*>   00 01 02 03 04 05       00
*>      11 12 13 14 15       10 11
*>         22 23 24 25       20 21 22
*>            33 34 35       30 31 32 33
*>               44 45       40 41 42 43 44
*>                  55       50 51 52 53 54 55
*>
*>
*>  Let TRANSR = 'N'. RFP holds AP as follows:
*>  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last
*>  three columns of AP upper. The lower triangle A(4:6,0:2) consists of
*>  the transpose of the first three columns of AP upper.
*>  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first
*>  three columns of AP lower. The upper triangle A(0:2,0:2) consists of
*>  the transpose of the last three columns of AP lower.
*>  This covers the case N even and TRANSR = 'N'.
*>
*>         RFP A                   RFP A
*>
*>        03 04 05                33 43 53
*>        13 14 15                00 44 54
*>        23 24 25                10 11 55
*>        33 34 35                20 21 22
*>        00 44 45                30 31 32
*>        01 11 55                40 41 42
*>        02 12 22                50 51 52
*>
*>  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the
*>  transpose of RFP A above. One therefore gets:
*>
*>
*>           RFP A                   RFP A
*>
*>     03 13 23 33 00 01 02    33 00 10 20 30 40 50
*>     04 14 24 34 44 11 12    43 44 11 21 31 41 51
*>     05 15 25 35 45 55 22    53 54 55 22 32 42 52
*>
*>
*>  We then consider Rectangular Full Packed (RFP) Format when N is
*>  odd. We give an example where N = 5.
*>
*>     AP is Upper                 AP is Lower
*>
*>   00 01 02 03 04              00
*>      11 12 13 14              10 11
*>         22 23 24              20 21 22
*>            33 34              30 31 32 33
*>               44              40 41 42 43 44
*>
*>
*>  Let TRANSR = 'N'. RFP holds AP as follows:
*>  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last
*>  three columns of AP upper. The lower triangle A(3:4,0:1) consists of
*>  the transpose of the first two columns of AP upper.
*>  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first
*>  three columns of AP lower. The upper triangle A(0:1,1:2) consists of
*>  the transpose of the last two columns of AP lower.
*>  This covers the case N odd and TRANSR = 'N'.
*>
*>         RFP A                   RFP A
*>
*>        02 03 04                00 33 43
*>        12 13 14                10 11 44
*>        22 23 24                20 21 22
*>        00 33 34                30 31 32
*>        01 11 44                40 41 42
*>
*>  Now let TRANSR = 'T'. RFP A in both UPLO cases is just the
*>  transpose of RFP A above. One therefore gets:
*>
*>           RFP A                   RFP A
*>
*>     02 12 22 00 01             00 10 20 30 40 50
*>     03 13 23 33 11             33 11 21 31 41 51
*>     04 14 24 34 44             43 44 22 32 42 52
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE SPFTRF( TRANSR, UPLO, N, A, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANSR, UPLO
      INTEGER            N, INFO
*     ..
*     .. Array Arguments ..
      REAL               A( 0: * )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, NISODD, NORMALTRANSR
      INTEGER            N1, N2, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, SSYRK, SPOTRF, STRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPFTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     If N is odd, set NISODD = .TRUE.
*     If N is even, set K = N/2 and NISODD = .FALSE.
*
      IF( MOD( N, 2 ).EQ.0 ) THEN
         K = N / 2
         NISODD = .FALSE.
      ELSE
         NISODD = .TRUE.
      END IF
*
*     Set N1 and N2 depending on LOWER
*
      IF( LOWER ) THEN
         N2 = N / 2
         N1 = N - N2
      ELSE
         N1 = N / 2
         N2 = N - N1
      END IF
*
*     start execution: there are eight cases
*
      IF( NISODD ) THEN
*
*        N is odd
*
         IF( NORMALTRANSR ) THEN
*
*           N is odd and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
*             T1 -> a(0), T2 -> a(n), S -> a(n1)
*
               CALL SPOTRF( 'L', N1, A( 0 ), N, INFO )
               IF( INFO.GT.0 )
     $            RETURN
               CALL STRSM( 'R', 'L', 'T', 'N', N2, N1, ONE, A( 0 ),
     $                     N,
     $                     A( N1 ), N )
               CALL SSYRK( 'U', 'N', N2, N1, -ONE, A( N1 ), N, ONE,
     $                     A( N ), N )
               CALL SPOTRF( 'U', N2, A( N ), N, INFO )
               IF( INFO.GT.0 )
     $            INFO = INFO + N1
*
            ELSE
*
*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
*             T1 -> a(n2), T2 -> a(n1), S -> a(0)
*
               CALL SPOTRF( 'L', N1, A( N2 ), N, INFO )
               IF( INFO.GT.0 )
     $            RETURN
               CALL STRSM( 'L', 'L', 'N', 'N', N1, N2, ONE, A( N2 ),
     $                     N,
     $                     A( 0 ), N )
               CALL SSYRK( 'U', 'T', N2, N1, -ONE, A( 0 ), N, ONE,
     $                     A( N1 ), N )
               CALL SPOTRF( 'U', N2, A( N1 ), N, INFO )
               IF( INFO.GT.0 )
     $            INFO = INFO + N1
*
            END IF
*
         ELSE
*
*           N is odd and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, TRANSPOSE and N is odd
*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
*              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1
*
               CALL SPOTRF( 'U', N1, A( 0 ), N1, INFO )
               IF( INFO.GT.0 )
     $            RETURN
               CALL STRSM( 'L', 'U', 'T', 'N', N1, N2, ONE, A( 0 ),
     $                     N1,
     $                     A( N1*N1 ), N1 )
               CALL SSYRK( 'L', 'T', N2, N1, -ONE, A( N1*N1 ), N1,
     $                     ONE,
     $                     A( 1 ), N1 )
               CALL SPOTRF( 'L', N2, A( 1 ), N1, INFO )
               IF( INFO.GT.0 )
     $            INFO = INFO + N1
*
            ELSE
*
*              SRPA for UPPER, TRANSPOSE and N is odd
*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
*              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2
*
               CALL SPOTRF( 'U', N1, A( N2*N2 ), N2, INFO )
               IF( INFO.GT.0 )
     $            RETURN
               CALL STRSM( 'R', 'U', 'N', 'N', N2, N1, ONE,
     $                     A( N2*N2 ),
     $                     N2, A( 0 ), N2 )
               CALL SSYRK( 'L', 'N', N2, N1, -ONE, A( 0 ), N2, ONE,
     $                     A( N1*N2 ), N2 )
               CALL SPOTRF( 'L', N2, A( N1*N2 ), N2, INFO )
               IF( INFO.GT.0 )
     $            INFO = INFO + N1
*
            END IF
*
         END IF
*
      ELSE
*
*        N is even
*
         IF( NORMALTRANSR ) THEN
*
*           N is even and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
*              T1 -> a(1), T2 -> a(0), S -> a(k+1)
*
               CALL SPOTRF( 'L', K, A( 1 ), N+1, INFO )
               IF( INFO.GT.0 )
     $            RETURN
               CALL STRSM( 'R', 'L', 'T', 'N', K, K, ONE, A( 1 ),
     $                     N+1,
     $                     A( K+1 ), N+1 )
               CALL SSYRK( 'U', 'N', K, K, -ONE, A( K+1 ), N+1, ONE,
     $                     A( 0 ), N+1 )
               CALL SPOTRF( 'U', K, A( 0 ), N+1, INFO )
               IF( INFO.GT.0 )
     $            INFO = INFO + K
*
            ELSE
*
*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
*              T1 -> a(k+1), T2 -> a(k), S -> a(0)
*
               CALL SPOTRF( 'L', K, A( K+1 ), N+1, INFO )
               IF( INFO.GT.0 )
     $            RETURN
               CALL STRSM( 'L', 'L', 'N', 'N', K, K, ONE, A( K+1 ),
     $                     N+1, A( 0 ), N+1 )
               CALL SSYRK( 'U', 'T', K, K, -ONE, A( 0 ), N+1, ONE,
     $                     A( K ), N+1 )
               CALL SPOTRF( 'U', K, A( K ), N+1, INFO )
               IF( INFO.GT.0 )
     $            INFO = INFO + K
*
            END IF
*
         ELSE
*
*           N is even and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, TRANSPOSE and N is even (see paper)
*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
*
               CALL SPOTRF( 'U', K, A( 0+K ), K, INFO )
               IF( INFO.GT.0 )
     $            RETURN
               CALL STRSM( 'L', 'U', 'T', 'N', K, K, ONE, A( K ), N1,
     $                     A( K*( K+1 ) ), K )
               CALL SSYRK( 'L', 'T', K, K, -ONE, A( K*( K+1 ) ), K,
     $                     ONE,
     $                     A( 0 ), K )
               CALL SPOTRF( 'L', K, A( 0 ), K, INFO )
               IF( INFO.GT.0 )
     $            INFO = INFO + K
*
            ELSE
*
*              SRPA for UPPER, TRANSPOSE and N is even (see paper)
*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
*
               CALL SPOTRF( 'U', K, A( K*( K+1 ) ), K, INFO )
               IF( INFO.GT.0 )
     $            RETURN
               CALL STRSM( 'R', 'U', 'N', 'N', K, K, ONE,
     $                     A( K*( K+1 ) ), K, A( 0 ), K )
               CALL SSYRK( 'L', 'N', K, K, -ONE, A( 0 ), K, ONE,
     $                     A( K*K ), K )
               CALL SPOTRF( 'L', K, A( K*K ), K, INFO )
               IF( INFO.GT.0 )
     $            INFO = INFO + K
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of SPFTRF
*
      END
