*> \brief \b CTFSM solves a matrix equation (one operand is a triangular matrix in RFP format).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CTFSM + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctfsm.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctfsm.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctfsm.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A,
*                         B, LDB )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANSR, DIAG, SIDE, TRANS, UPLO
*       INTEGER            LDB, M, N
*       COMPLEX            ALPHA
*       ..
*       .. Array Arguments ..
*       COMPLEX            A( 0: * ), B( 0: LDB-1, 0: * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> Level 3 BLAS like routine for A in RFP Format.
*>
*> CTFSM solves the matrix equation
*>
*>    op( A )*X = alpha*B  or  X*op( A ) = alpha*B
*>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**H.
*>
*> A is in Rectangular Full Packed (RFP) Format.
*>
*> The matrix X is overwritten on B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANSR
*> \verbatim
*>          TRANSR is CHARACTER*1
*>          = 'N':  The Normal Form of RFP A is stored;
*>          = 'C':  The Conjugate-transpose Form of RFP A is stored.
*> \endverbatim
*>
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether op( A ) appears on the left
*>           or right of X as follows:
*>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*>
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the RFP matrix A came from
*>           an upper or lower triangular matrix as follows:
*>           UPLO = 'U' or 'u' RFP A came from an upper triangular matrix
*>           UPLO = 'L' or 'l' RFP A came from a  lower triangular matrix
*>
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS  specifies the form of op( A ) to be used
*>           in the matrix multiplication as follows:
*>
*>              TRANS  = 'N' or 'n'   op( A ) = A.
*>
*>              TRANS  = 'C' or 'c'   op( A ) = conjg( A' ).
*>
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not RFP A is unit
*>           triangular as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*>
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*>           zero then  A is not referenced and  B need not be set before
*>           entry.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension (NT)
*>           NT = N*(N+1)/2 if SIDE='R' and NT = M*(M+1)/2 otherwise.
*>           On entry, the matrix A in RFP Format. 
*>           RFP Format is described by TRANSR, UPLO and N as follows:
*>           If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;
*>           K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If
*>           TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A as
*>           defined when TRANSR = 'N'. The contents of RFP A are defined
*>           by UPLO as follows: If UPLO = 'U' the RFP A contains the NT
*>           elements of upper packed A either in normal or
*>           conjugate-transpose Format. If UPLO = 'L' the RFP A contains
*>           the NT elements of lower packed A either in normal or
*>           conjugate-transpose Format. The LDA of RFP A is (N+1)/2 when
*>           TRANSR = 'C'. When TRANSR is 'N' the LDA is N+1 when N is
*>           even and is N when is odd.
*>           See the Note below for more details. Unchanged on exit.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,N)
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain  the  right-hand  side  matrix  B,  and  on exit  is
*>           overwritten by the solution matrix  X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
*>           max( 1, m ).
*>           Unchanged on exit.
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
*> \ingroup tfsm
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  We first consider Standard Packed Format when N is even.
*>  We give an example where N = 6.
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
*>  conjugate-transpose of the first three columns of AP upper.
*>  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first
*>  three columns of AP lower. The upper triangle A(0:2,0:2) consists of
*>  conjugate-transpose of the last three columns of AP lower.
*>  To denote conjugate we place -- above the element. This covers the
*>  case N even and TRANSR = 'N'.
*>
*>         RFP A                   RFP A
*>
*>                                -- -- --
*>        03 04 05                33 43 53
*>                                   -- --
*>        13 14 15                00 44 54
*>                                      --
*>        23 24 25                10 11 55
*>
*>        33 34 35                20 21 22
*>        --
*>        00 44 45                30 31 32
*>        -- --
*>        01 11 55                40 41 42
*>        -- -- --
*>        02 12 22                50 51 52
*>
*>  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-
*>  transpose of RFP A above. One therefore gets:
*>
*>
*>           RFP A                   RFP A
*>
*>     -- -- -- --                -- -- -- -- -- --
*>     03 13 23 33 00 01 02    33 00 10 20 30 40 50
*>     -- -- -- -- --                -- -- -- -- --
*>     04 14 24 34 44 11 12    43 44 11 21 31 41 51
*>     -- -- -- -- -- --                -- -- -- --
*>     05 15 25 35 45 55 22    53 54 55 22 32 42 52
*>
*>
*>  We next  consider Standard Packed Format when N is odd.
*>  We give an example where N = 5.
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
*>  conjugate-transpose of the first two   columns of AP upper.
*>  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first
*>  three columns of AP lower. The upper triangle A(0:1,1:2) consists of
*>  conjugate-transpose of the last two   columns of AP lower.
*>  To denote conjugate we place -- above the element. This covers the
*>  case N odd  and TRANSR = 'N'.
*>
*>         RFP A                   RFP A
*>
*>                                   -- --
*>        02 03 04                00 33 43
*>                                      --
*>        12 13 14                10 11 44
*>
*>        22 23 24                20 21 22
*>        --
*>        00 33 34                30 31 32
*>        -- --
*>        01 11 44                40 41 42
*>
*>  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-
*>  transpose of RFP A above. One therefore gets:
*>
*>
*>           RFP A                   RFP A
*>
*>     -- -- --                   -- -- -- -- -- --
*>     02 12 22 00 01             00 10 20 30 40 50
*>     -- -- -- --                   -- -- -- -- --
*>     03 13 23 33 11             33 11 21 31 41 51
*>     -- -- -- -- --                   -- -- -- --
*>     04 14 24 34 44             43 44 22 32 42 52
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA,
     $                  A,
     $                  B, LDB )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANSR, DIAG, SIDE, TRANS, UPLO
      INTEGER            LDB, M, N
      COMPLEX            ALPHA
*     ..
*     .. Array Arguments ..
      COMPLEX            A( 0: * ), B( 0: LDB-1, 0: * )
*     ..
*
*  =====================================================================
*     ..
*     .. Parameters ..
      COMPLEX            CONE, CZERO
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ),
     $                   CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LSIDE, MISODD, NISODD, NORMALTRANSR,
     $                   NOTRANS
      INTEGER            M1, M2, N1, N2, K, INFO, I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CGEMM, CTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LSIDE = LSAME( SIDE, 'L' )
      LOWER = LSAME( UPLO, 'L' )
      NOTRANS = LSAME( TRANS, 'N' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSIDE .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.NOTRANS .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -4
      ELSE IF( .NOT.LSAME( DIAG, 'N' ) .AND.
     $         .NOT.LSAME( DIAG, 'U' ) )
     $         THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTFSM ', -INFO )
         RETURN
      END IF
*
*     Quick return when ( (N.EQ.0).OR.(M.EQ.0) )
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
*
*     Quick return when ALPHA.EQ.(0E+0,0E+0)
*
      IF( ALPHA.EQ.CZERO ) THEN
         DO 20 J = 0, N - 1
            DO 10 I = 0, M - 1
               B( I, J ) = CZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
      IF( LSIDE ) THEN
*
*        SIDE = 'L'
*
*        A is M-by-M.
*        If M is odd, set NISODD = .TRUE., and M1 and M2.
*        If M is even, NISODD = .FALSE., and M.
*
         IF( MOD( M, 2 ).EQ.0 ) THEN
            MISODD = .FALSE.
            K = M / 2
         ELSE
            MISODD = .TRUE.
            IF( LOWER ) THEN
               M2 = M / 2
               M1 = M - M2
            ELSE
               M1 = M / 2
               M2 = M - M1
            END IF
         END IF
*
         IF( MISODD ) THEN
*
*           SIDE = 'L' and N is odd
*
            IF( NORMALTRANSR ) THEN
*
*              SIDE = 'L', N is odd, and TRANSR = 'N'
*
               IF( LOWER ) THEN
*
*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
*                    TRANS = 'N'
*
                     IF( M.EQ.1 ) THEN
                        CALL CTRSM( 'L', 'L', 'N', DIAG, M1, N,
     $                              ALPHA,
     $                              A, M, B, LDB )
                     ELSE
                        CALL CTRSM( 'L', 'L', 'N', DIAG, M1, N,
     $                              ALPHA,
     $                              A( 0 ), M, B, LDB )
                        CALL CGEMM( 'N', 'N', M2, N, M1, -CONE,
     $                              A( M1 ),
     $                              M, B, LDB, ALPHA, B( M1, 0 ), LDB )
                        CALL CTRSM( 'L', 'U', 'C', DIAG, M2, N, CONE,
     $                              A( M ), M, B( M1, 0 ), LDB )
                     END IF
*
                  ELSE
*
*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
*                    TRANS = 'C'
*
                     IF( M.EQ.1 ) THEN
                        CALL CTRSM( 'L', 'L', 'C', DIAG, M1, N,
     $                              ALPHA,
     $                              A( 0 ), M, B, LDB )
                     ELSE
                        CALL CTRSM( 'L', 'U', 'N', DIAG, M2, N,
     $                              ALPHA,
     $                              A( M ), M, B( M1, 0 ), LDB )
                        CALL CGEMM( 'C', 'N', M1, N, M2, -CONE,
     $                              A( M1 ),
     $                              M, B( M1, 0 ), LDB, ALPHA, B, LDB )
                        CALL CTRSM( 'L', 'L', 'C', DIAG, M1, N, CONE,
     $                              A( 0 ), M, B, LDB )
                     END IF
*
                  END IF
*
               ELSE
*
*                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U'
*
                  IF( .NOT.NOTRANS ) THEN
*
*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
*                    TRANS = 'N'
*
                     CALL CTRSM( 'L', 'L', 'N', DIAG, M1, N, ALPHA,
     $                           A( M2 ), M, B, LDB )
                     CALL CGEMM( 'C', 'N', M2, N, M1, -CONE, A( 0 ),
     $                           M,
     $                           B, LDB, ALPHA, B( M1, 0 ), LDB )
                     CALL CTRSM( 'L', 'U', 'C', DIAG, M2, N, CONE,
     $                           A( M1 ), M, B( M1, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
*                    TRANS = 'C'
*
                     CALL CTRSM( 'L', 'U', 'N', DIAG, M2, N, ALPHA,
     $                           A( M1 ), M, B( M1, 0 ), LDB )
                     CALL CGEMM( 'N', 'N', M1, N, M2, -CONE, A( 0 ),
     $                           M,
     $                           B( M1, 0 ), LDB, ALPHA, B, LDB )
                     CALL CTRSM( 'L', 'L', 'C', DIAG, M1, N, CONE,
     $                           A( M2 ), M, B, LDB )
*
                  END IF
*
               END IF
*
            ELSE
*
*              SIDE = 'L', N is odd, and TRANSR = 'C'
*
               IF( LOWER ) THEN
*
*                 SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and
*                    TRANS = 'N'
*
                     IF( M.EQ.1 ) THEN
                        CALL CTRSM( 'L', 'U', 'C', DIAG, M1, N,
     $                              ALPHA,
     $                              A( 0 ), M1, B, LDB )
                     ELSE
                        CALL CTRSM( 'L', 'U', 'C', DIAG, M1, N,
     $                              ALPHA,
     $                              A( 0 ), M1, B, LDB )
                        CALL CGEMM( 'C', 'N', M2, N, M1, -CONE,
     $                              A( M1*M1 ), M1, B, LDB, ALPHA,
     $                              B( M1, 0 ), LDB )
                        CALL CTRSM( 'L', 'L', 'N', DIAG, M2, N, CONE,
     $                              A( 1 ), M1, B( M1, 0 ), LDB )
                     END IF
*
                  ELSE
*
*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and
*                    TRANS = 'C'
*
                     IF( M.EQ.1 ) THEN
                        CALL CTRSM( 'L', 'U', 'N', DIAG, M1, N,
     $                              ALPHA,
     $                              A( 0 ), M1, B, LDB )
                     ELSE
                        CALL CTRSM( 'L', 'L', 'C', DIAG, M2, N,
     $                              ALPHA,
     $                              A( 1 ), M1, B( M1, 0 ), LDB )
                        CALL CGEMM( 'N', 'N', M1, N, M2, -CONE,
     $                              A( M1*M1 ), M1, B( M1, 0 ), LDB,
     $                              ALPHA, B, LDB )
                        CALL CTRSM( 'L', 'U', 'N', DIAG, M1, N, CONE,
     $                              A( 0 ), M1, B, LDB )
                     END IF
*
                  END IF
*
               ELSE
*
*                 SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'U'
*
                  IF( .NOT.NOTRANS ) THEN
*
*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and
*                    TRANS = 'N'
*
                     CALL CTRSM( 'L', 'U', 'C', DIAG, M1, N, ALPHA,
     $                           A( M2*M2 ), M2, B, LDB )
                     CALL CGEMM( 'N', 'N', M2, N, M1, -CONE, A( 0 ),
     $                           M2,
     $                           B, LDB, ALPHA, B( M1, 0 ), LDB )
                     CALL CTRSM( 'L', 'L', 'N', DIAG, M2, N, CONE,
     $                           A( M1*M2 ), M2, B( M1, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and
*                    TRANS = 'C'
*
                     CALL CTRSM( 'L', 'L', 'C', DIAG, M2, N, ALPHA,
     $                           A( M1*M2 ), M2, B( M1, 0 ), LDB )
                     CALL CGEMM( 'C', 'N', M1, N, M2, -CONE, A( 0 ),
     $                           M2,
     $                           B( M1, 0 ), LDB, ALPHA, B, LDB )
                     CALL CTRSM( 'L', 'U', 'N', DIAG, M1, N, CONE,
     $                           A( M2*M2 ), M2, B, LDB )
*
                  END IF
*
               END IF
*
            END IF
*
         ELSE
*
*           SIDE = 'L' and N is even
*
            IF( NORMALTRANSR ) THEN
*
*              SIDE = 'L', N is even, and TRANSR = 'N'
*
               IF( LOWER ) THEN
*
*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
*                    and TRANS = 'N'
*
                     CALL CTRSM( 'L', 'L', 'N', DIAG, K, N, ALPHA,
     $                           A( 1 ), M+1, B, LDB )
                     CALL CGEMM( 'N', 'N', K, N, K, -CONE, A( K+1 ),
     $                           M+1, B, LDB, ALPHA, B( K, 0 ), LDB )
                     CALL CTRSM( 'L', 'U', 'C', DIAG, K, N, CONE,
     $                           A( 0 ), M+1, B( K, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
*                    and TRANS = 'C'
*
                     CALL CTRSM( 'L', 'U', 'N', DIAG, K, N, ALPHA,
     $                           A( 0 ), M+1, B( K, 0 ), LDB )
                     CALL CGEMM( 'C', 'N', K, N, K, -CONE, A( K+1 ),
     $                           M+1, B( K, 0 ), LDB, ALPHA, B, LDB )
                     CALL CTRSM( 'L', 'L', 'C', DIAG, K, N, CONE,
     $                           A( 1 ), M+1, B, LDB )
*
                  END IF
*
               ELSE
*
*                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U'
*
                  IF( .NOT.NOTRANS ) THEN
*
*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
*                    and TRANS = 'N'
*
                     CALL CTRSM( 'L', 'L', 'N', DIAG, K, N, ALPHA,
     $                           A( K+1 ), M+1, B, LDB )
                     CALL CGEMM( 'C', 'N', K, N, K, -CONE, A( 0 ),
     $                           M+1,
     $                           B, LDB, ALPHA, B( K, 0 ), LDB )
                     CALL CTRSM( 'L', 'U', 'C', DIAG, K, N, CONE,
     $                           A( K ), M+1, B( K, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
*                    and TRANS = 'C'
                     CALL CTRSM( 'L', 'U', 'N', DIAG, K, N, ALPHA,
     $                           A( K ), M+1, B( K, 0 ), LDB )
                     CALL CGEMM( 'N', 'N', K, N, K, -CONE, A( 0 ),
     $                           M+1,
     $                           B( K, 0 ), LDB, ALPHA, B, LDB )
                     CALL CTRSM( 'L', 'L', 'C', DIAG, K, N, CONE,
     $                           A( K+1 ), M+1, B, LDB )
*
                  END IF
*
               END IF
*
            ELSE
*
*              SIDE = 'L', N is even, and TRANSR = 'C'
*
               IF( LOWER ) THEN
*
*                 SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L',
*                    and TRANS = 'N'
*
                     CALL CTRSM( 'L', 'U', 'C', DIAG, K, N, ALPHA,
     $                           A( K ), K, B, LDB )
                     CALL CGEMM( 'C', 'N', K, N, K, -CONE,
     $                           A( K*( K+1 ) ), K, B, LDB, ALPHA,
     $                           B( K, 0 ), LDB )
                     CALL CTRSM( 'L', 'L', 'N', DIAG, K, N, CONE,
     $                           A( 0 ), K, B( K, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L',
*                    and TRANS = 'C'
*
                     CALL CTRSM( 'L', 'L', 'C', DIAG, K, N, ALPHA,
     $                           A( 0 ), K, B( K, 0 ), LDB )
                     CALL CGEMM( 'N', 'N', K, N, K, -CONE,
     $                           A( K*( K+1 ) ), K, B( K, 0 ), LDB,
     $                           ALPHA, B, LDB )
                     CALL CTRSM( 'L', 'U', 'N', DIAG, K, N, CONE,
     $                           A( K ), K, B, LDB )
*
                  END IF
*
               ELSE
*
*                 SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'U'
*
                  IF( .NOT.NOTRANS ) THEN
*
*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U',
*                    and TRANS = 'N'
*
                     CALL CTRSM( 'L', 'U', 'C', DIAG, K, N, ALPHA,
     $                           A( K*( K+1 ) ), K, B, LDB )
                     CALL CGEMM( 'N', 'N', K, N, K, -CONE, A( 0 ), K,
     $                           B,
     $                           LDB, ALPHA, B( K, 0 ), LDB )
                     CALL CTRSM( 'L', 'L', 'N', DIAG, K, N, CONE,
     $                           A( K*K ), K, B( K, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U',
*                    and TRANS = 'C'
*
                     CALL CTRSM( 'L', 'L', 'C', DIAG, K, N, ALPHA,
     $                           A( K*K ), K, B( K, 0 ), LDB )
                     CALL CGEMM( 'C', 'N', K, N, K, -CONE, A( 0 ), K,
     $                           B( K, 0 ), LDB, ALPHA, B, LDB )
                     CALL CTRSM( 'L', 'U', 'N', DIAG, K, N, CONE,
     $                           A( K*( K+1 ) ), K, B, LDB )
*
                  END IF
*
               END IF
*
            END IF
*
         END IF
*
      ELSE
*
*        SIDE = 'R'
*
*        A is N-by-N.
*        If N is odd, set NISODD = .TRUE., and N1 and N2.
*        If N is even, NISODD = .FALSE., and K.
*
         IF( MOD( N, 2 ).EQ.0 ) THEN
            NISODD = .FALSE.
            K = N / 2
         ELSE
            NISODD = .TRUE.
            IF( LOWER ) THEN
               N2 = N / 2
               N1 = N - N2
            ELSE
               N1 = N / 2
               N2 = N - N1
            END IF
         END IF
*
         IF( NISODD ) THEN
*
*           SIDE = 'R' and N is odd
*
            IF( NORMALTRANSR ) THEN
*
*              SIDE = 'R', N is odd, and TRANSR = 'N'
*
               IF( LOWER ) THEN
*
*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
*                    TRANS = 'N'
*
                     CALL CTRSM( 'R', 'U', 'C', DIAG, M, N2, ALPHA,
     $                           A( N ), N, B( 0, N1 ), LDB )
                     CALL CGEMM( 'N', 'N', M, N1, N2, -CONE, B( 0,
     $                           N1 ),
     $                           LDB, A( N1 ), N, ALPHA, B( 0, 0 ),
     $                           LDB )
                     CALL CTRSM( 'R', 'L', 'N', DIAG, M, N1, CONE,
     $                           A( 0 ), N, B( 0, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
*                    TRANS = 'C'
*
                     CALL CTRSM( 'R', 'L', 'C', DIAG, M, N1, ALPHA,
     $                           A( 0 ), N, B( 0, 0 ), LDB )
                     CALL CGEMM( 'N', 'C', M, N2, N1, -CONE, B( 0,
     $                           0 ),
     $                           LDB, A( N1 ), N, ALPHA, B( 0, N1 ),
     $                           LDB )
                     CALL CTRSM( 'R', 'U', 'N', DIAG, M, N2, CONE,
     $                           A( N ), N, B( 0, N1 ), LDB )
*
                  END IF
*
               ELSE
*
*                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
*                    TRANS = 'N'
*
                     CALL CTRSM( 'R', 'L', 'C', DIAG, M, N1, ALPHA,
     $                           A( N2 ), N, B( 0, 0 ), LDB )
                     CALL CGEMM( 'N', 'N', M, N2, N1, -CONE, B( 0,
     $                           0 ),
     $                           LDB, A( 0 ), N, ALPHA, B( 0, N1 ),
     $                           LDB )
                     CALL CTRSM( 'R', 'U', 'N', DIAG, M, N2, CONE,
     $                           A( N1 ), N, B( 0, N1 ), LDB )
*
                  ELSE
*
*                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
*                    TRANS = 'C'
*
                     CALL CTRSM( 'R', 'U', 'C', DIAG, M, N2, ALPHA,
     $                           A( N1 ), N, B( 0, N1 ), LDB )
                     CALL CGEMM( 'N', 'C', M, N1, N2, -CONE, B( 0,
     $                           N1 ),
     $                           LDB, A( 0 ), N, ALPHA, B( 0, 0 ), LDB )
                     CALL CTRSM( 'R', 'L', 'N', DIAG, M, N1, CONE,
     $                           A( N2 ), N, B( 0, 0 ), LDB )
*
                  END IF
*
               END IF
*
            ELSE
*
*              SIDE = 'R', N is odd, and TRANSR = 'C'
*
               IF( LOWER ) THEN
*
*                 SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and
*                    TRANS = 'N'
*
                     CALL CTRSM( 'R', 'L', 'N', DIAG, M, N2, ALPHA,
     $                           A( 1 ), N1, B( 0, N1 ), LDB )
                     CALL CGEMM( 'N', 'C', M, N1, N2, -CONE, B( 0,
     $                           N1 ),
     $                           LDB, A( N1*N1 ), N1, ALPHA, B( 0, 0 ),
     $                           LDB )
                     CALL CTRSM( 'R', 'U', 'C', DIAG, M, N1, CONE,
     $                           A( 0 ), N1, B( 0, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and
*                    TRANS = 'C'
*
                     CALL CTRSM( 'R', 'U', 'N', DIAG, M, N1, ALPHA,
     $                           A( 0 ), N1, B( 0, 0 ), LDB )
                     CALL CGEMM( 'N', 'N', M, N2, N1, -CONE, B( 0,
     $                           0 ),
     $                           LDB, A( N1*N1 ), N1, ALPHA, B( 0, N1 ),
     $                           LDB )
                     CALL CTRSM( 'R', 'L', 'C', DIAG, M, N2, CONE,
     $                           A( 1 ), N1, B( 0, N1 ), LDB )
*
                  END IF
*
               ELSE
*
*                 SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'U'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and
*                    TRANS = 'N'
*
                     CALL CTRSM( 'R', 'U', 'N', DIAG, M, N1, ALPHA,
     $                           A( N2*N2 ), N2, B( 0, 0 ), LDB )
                     CALL CGEMM( 'N', 'C', M, N2, N1, -CONE, B( 0,
     $                           0 ),
     $                           LDB, A( 0 ), N2, ALPHA, B( 0, N1 ),
     $                           LDB )
                     CALL CTRSM( 'R', 'L', 'C', DIAG, M, N2, CONE,
     $                           A( N1*N2 ), N2, B( 0, N1 ), LDB )
*
                  ELSE
*
*                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and
*                    TRANS = 'C'
*
                     CALL CTRSM( 'R', 'L', 'N', DIAG, M, N2, ALPHA,
     $                           A( N1*N2 ), N2, B( 0, N1 ), LDB )
                     CALL CGEMM( 'N', 'N', M, N1, N2, -CONE, B( 0,
     $                           N1 ),
     $                           LDB, A( 0 ), N2, ALPHA, B( 0, 0 ),
     $                           LDB )
                     CALL CTRSM( 'R', 'U', 'C', DIAG, M, N1, CONE,
     $                           A( N2*N2 ), N2, B( 0, 0 ), LDB )
*
                  END IF
*
               END IF
*
            END IF
*
         ELSE
*
*           SIDE = 'R' and N is even
*
            IF( NORMALTRANSR ) THEN
*
*              SIDE = 'R', N is even, and TRANSR = 'N'
*
               IF( LOWER ) THEN
*
*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
*                    and TRANS = 'N'
*
                     CALL CTRSM( 'R', 'U', 'C', DIAG, M, K, ALPHA,
     $                           A( 0 ), N+1, B( 0, K ), LDB )
                     CALL CGEMM( 'N', 'N', M, K, K, -CONE, B( 0, K ),
     $                           LDB, A( K+1 ), N+1, ALPHA, B( 0, 0 ),
     $                           LDB )
                     CALL CTRSM( 'R', 'L', 'N', DIAG, M, K, CONE,
     $                           A( 1 ), N+1, B( 0, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
*                    and TRANS = 'C'
*
                     CALL CTRSM( 'R', 'L', 'C', DIAG, M, K, ALPHA,
     $                           A( 1 ), N+1, B( 0, 0 ), LDB )
                     CALL CGEMM( 'N', 'C', M, K, K, -CONE, B( 0, 0 ),
     $                           LDB, A( K+1 ), N+1, ALPHA, B( 0, K ),
     $                           LDB )
                     CALL CTRSM( 'R', 'U', 'N', DIAG, M, K, CONE,
     $                           A( 0 ), N+1, B( 0, K ), LDB )
*
                  END IF
*
               ELSE
*
*                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
*                    and TRANS = 'N'
*
                     CALL CTRSM( 'R', 'L', 'C', DIAG, M, K, ALPHA,
     $                           A( K+1 ), N+1, B( 0, 0 ), LDB )
                     CALL CGEMM( 'N', 'N', M, K, K, -CONE, B( 0, 0 ),
     $                           LDB, A( 0 ), N+1, ALPHA, B( 0, K ),
     $                           LDB )
                     CALL CTRSM( 'R', 'U', 'N', DIAG, M, K, CONE,
     $                           A( K ), N+1, B( 0, K ), LDB )
*
                  ELSE
*
*                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
*                    and TRANS = 'C'
*
                     CALL CTRSM( 'R', 'U', 'C', DIAG, M, K, ALPHA,
     $                           A( K ), N+1, B( 0, K ), LDB )
                     CALL CGEMM( 'N', 'C', M, K, K, -CONE, B( 0, K ),
     $                           LDB, A( 0 ), N+1, ALPHA, B( 0, 0 ),
     $                           LDB )
                     CALL CTRSM( 'R', 'L', 'N', DIAG, M, K, CONE,
     $                           A( K+1 ), N+1, B( 0, 0 ), LDB )
*
                  END IF
*
               END IF
*
            ELSE
*
*              SIDE = 'R', N is even, and TRANSR = 'C'
*
               IF( LOWER ) THEN
*
*                 SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L',
*                    and TRANS = 'N'
*
                     CALL CTRSM( 'R', 'L', 'N', DIAG, M, K, ALPHA,
     $                           A( 0 ), K, B( 0, K ), LDB )
                     CALL CGEMM( 'N', 'C', M, K, K, -CONE, B( 0, K ),
     $                           LDB, A( ( K+1 )*K ), K, ALPHA,
     $                           B( 0, 0 ), LDB )
                     CALL CTRSM( 'R', 'U', 'C', DIAG, M, K, CONE,
     $                           A( K ), K, B( 0, 0 ), LDB )
*
                  ELSE
*
*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L',
*                    and TRANS = 'C'
*
                     CALL CTRSM( 'R', 'U', 'N', DIAG, M, K, ALPHA,
     $                           A( K ), K, B( 0, 0 ), LDB )
                     CALL CGEMM( 'N', 'N', M, K, K, -CONE, B( 0, 0 ),
     $                           LDB, A( ( K+1 )*K ), K, ALPHA,
     $                           B( 0, K ), LDB )
                     CALL CTRSM( 'R', 'L', 'C', DIAG, M, K, CONE,
     $                           A( 0 ), K, B( 0, K ), LDB )
*
                  END IF
*
               ELSE
*
*                 SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'U'
*
                  IF( NOTRANS ) THEN
*
*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U',
*                    and TRANS = 'N'
*
                     CALL CTRSM( 'R', 'U', 'N', DIAG, M, K, ALPHA,
     $                           A( ( K+1 )*K ), K, B( 0, 0 ), LDB )
                     CALL CGEMM( 'N', 'C', M, K, K, -CONE, B( 0, 0 ),
     $                           LDB, A( 0 ), K, ALPHA, B( 0, K ), LDB )
                     CALL CTRSM( 'R', 'L', 'C', DIAG, M, K, CONE,
     $                           A( K*K ), K, B( 0, K ), LDB )
*
                  ELSE
*
*                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U',
*                    and TRANS = 'C'
*
                     CALL CTRSM( 'R', 'L', 'N', DIAG, M, K, ALPHA,
     $                           A( K*K ), K, B( 0, K ), LDB )
                     CALL CGEMM( 'N', 'N', M, K, K, -CONE, B( 0, K ),
     $                           LDB, A( 0 ), K, ALPHA, B( 0, 0 ), LDB )
                     CALL CTRSM( 'R', 'U', 'C', DIAG, M, K, CONE,
     $                           A( ( K+1 )*K ), K, B( 0, 0 ), LDB )
*
                  END IF
*
               END IF
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of CTFSM
*
      END
