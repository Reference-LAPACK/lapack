*> \brief \b SLAORHR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SLAORHR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaorhr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaorhr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaorhr.f">
*> [TXT]</a>
*>
*  Definition:
*  ===========
*
*       SUBROUTINE SLAORHR( M, N, MB1, NB1, A, LDA, T1, LDT1, NB2, T2,
*      $                    LDT2, D, W, LDW, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER           INFO, LDA, LDT1, LDT2, LDW, LWORK, M, N, MB1,
*      $                  NB1, NB2
*       ..
*       .. Array Arguments ..
*       REAL              A( LDA, * ), D( * ), T1( LDT1, * ),
*      $                  T2( LDT2, * ),  W( LDW, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SLAORHR is an auxiliary routine of SORHR. SLAORHR is a low-level
*> routine that does the computation. SLAORHR performs Householder
*> reconstruction from the orthogonal M-by-M matrix Q_in stored in
*> an implicit form in the matrices A and T1 computed by the tall-skinny
*> QR factorization routine DLATSQR. The reconstructed orthogonal M-by-M
*> matrix Q_out is stored in a different implicit form in the matrices
*> A and T2. This implicit output format in A and T2 is the same as the
*> output in A and T from QR factorization routine SGEQRT. SLAORHR uses
*> two work arrays: a two-dimensional array W of dimensions (LDW, N) and
*> a one-dimensional array WORK of dimension LWORK.
*> See Further Details section.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] MB1
*> \verbatim
*>          MB1 is INTEGER
*>          The row block size used by SLATSQR to compute the input
*>          arrays A and T1 (MB1 is called MB in SLATSQR,
*>          see Further Details). MB1 > N.
*>          (Note that if MB1 > M, then M is used instead of MB1
*>          as the row block size).
*> \endverbatim
*>
*> \param[in] NB1
*> \verbatim
*>          NB1 is INTEGER
*>          The column block size used by SLATSQR to compute the input
*>          arrays A and T1 (NB1 is called NB in SLATSQR,
*>          see Further Details). NB1 >= 1.
*>          (Note that if NB1 > N, then N is used instead of NB1
*>          as the column block size).
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>
*>          The elements on and above the diagonal:
*>
*>             On entry, the elements on and above the diagonal
*>             contain the N-by-N upper-triangular matrix R_in
*>             computed by SLATSQR. See Further Details.
*>
*>             On exit, the elements on and above the diagonal
*>             contain the N-by-N upper-triangular matrix R_out
*>             reconstructed by this routine. See Further Details.
*>
*>             NOTE: In general R_in is not equal to R_out element-wise.
*>
*>
*>          The elements below the diagonal:
*>
*>             On entry, the elements below the diagonal represent
*>             the M-by-M matrix Q_in in implicit form by the columns
*>             of blocked matrix V computed by SLATSQR
*>             (same format as the output A in SLATSQR).
*>             See Further Details.
*>
*>             On exit, the elements below the diagonal represent
*>             the unit lower-trapezoidal M-by-N NB2-size column-blocked
*>             matrix Y, where the subdiagonal elements are the column
*>             vectors Y(i). The unit entries along the diagonal of Y
*>             are not stored in A. The matrix Y along with the matrix
*>             T2 defines Q_out. The vectors Y(i) define the elementary
*>             reflectors (same format as the output A in SGEQRT).
*>             See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] T1
*> \verbatim
*>          T1 is REAL array, dimension (LDT1, N * NIRB),
*.
*>          where NIRB = Number_of_input_row_blocks
*>                     = max( 1, CEIL((M-N)/(MB1-N)) )
*>          Let NICB = Number_of_input_col_blocks
*>                   = CEIL(N/NB1)
*>
*>          The upper-triangular block reflectors used to define
*>          Q_in stored in compact form in NIRB block reflector
*>          sequences. Each of NIRB block reflector sequences is stored
*>          in a larger NB1-by-N column block of T1 and consists of NICB
*>          smaller NB1-by-NB1 upper-triangular column blocks.
*>          (same format as the output T in SLATSQR).
*>          The matrix T1 and the matrix V stored in A define Q_in.
*>          See Further Details.
*> \endverbatim
*>
*> \param[in] LDT1
*> \verbatim
*>          LDT1 is INTEGER
*>          The leading dimension of the array T1.
*>          LDT1 >= max(1,min(NB1,N)).
*> \endverbatim
*>
*> \param[in] NB2
*> \verbatim
*>          NB2 is INTEGER
*>          The column block size to be used in the blocked Householder
*>          QR output of this routine in the array A and array T2.
*>          NB2 can be chosen independently of NB1. NB2 >= 1.
*>          (Note that if NB2 > N, then N is used instead of NB2
*>          as the column block size.)
*> \endverbatim
*>
*> \param[out] T2
*> \verbatim
*>          T2 is REAL array, dimension (LDT2, N)
*>
*>          Let NOCB = Number_of_output_col_blocks
*>                   = CEIL(N/NB2)
*>
*>          On exit, T2(1:NB2, 1:N) contains NOCB upper-triangular
*>          block reflectors used to define Q_out stored in compact
*>          form as a sequence of upper-triangular NB2-by-NB2 column
*>          blocks (same format as the output T in SGEQRT).
*>          The matrix T2 and the matrix Y stored in A define Q_out.
*>          See Further Details.
*> \endverbatim
*>
*> \param[in] LDT2
*> \verbatim
*>          LDT2 is INTEGER
*>          The leading dimension of the array T2.
*>          LDT2 >= max(1,min(NB2,N)).
*> \endverbatim
*>
*> \param[out] D
*> \verbatim
*>          D is REAL array, dimension min(M,N)
*>          The elements can be only plus or minus one.
*>
*>          D(i) is constructed as D(i)=-SIGN(Q1_in_i(i,i)), where
*>          1 <= i <= min(M,N), Q1_in is the left M-by-N
*>          submatrix of the input orthogonal M-by-M matrix Q_in
*>          implicitly defined by the matrix A and the matrix T1 on
*>          entry. Q1_in_i is the value after performing i-1 steps
*>          of "modified" Gaussian elimination. See Further Details.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is (workspace) REAL array, dimension (LDW, N)
*>          This workspace is used to store array C in the call to
*>          the SLAMTSQR routine.
*> \endverbatim
*>
*> \param[in] LDW
*> \verbatim
*>          LDW is INTEGER
*>          The leading dimension of the work array W.  LDW >= M.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is (workspace) REAL array,
*>          dimension (max(1,LWORK))
*>          This workspace is used as the workspace W in the call to
*>          the SLAMTSQR routine.
*>          See the documentation of SLAMTSQR.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          The dimension of the array WORK.  LWORK >= NB1*N.
*>          If LWORK = -1, then a workspace query is assumed.
*>          The routine only calculates the optimal size of the WORK
*>          array, returns this value as the first entry of the WORK
*>          array, and no error message related to LWORK is issued
*>          by XERBLA. See the documentation of SLAMTSQR.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*>
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*> =========
*> On entry:
*> =========
*>
*> See the documentation of SLATSQR.
*> On entry in A, the output from the M-by-N matrix A from SLATSQR is
*> used. The subdiagonal elements represent the M-by-M matrix Q_in in
*> the implicit form by the row and column blocked matrix V computed
*> by SLATSQR. The elements on and above the diagonal represent the
*> upper-triangular N-by-N matrix R_in computed by SLATSQR.
*>
*>          Let NIRB = Number_of_input_row_blocks
*>                   = max( 1, CEIL((M-N)/(MB1-N)) )
*>              NICB = Number_of_input_col_blocks
*>                   = CEIL(N/NB1)
*>
*> On entry in T1, the output from the matrix T from SLATSQR. The matrix
*> T1(1:NB1, 1:(N * NIRB)) contains the upper-triangular block
*> reflectors stored in compact form in NIRB block sequences
*> (larger blocks) of upper-triangular blocks. Each of the NIRB block
*> reflector sequences is stored in a larger NB1-by-N column block
*> of T1 and consists of NICB smaller NB1-by-NB1 upper-triangular column
*> blocks. Within each sequence,  the order of the last block is
*> LAST_NB1 = N - (NICB-1)*NB1.
*>
*> Therefore, the input column block size NB1 and input column block
*> size MB1 must be equal to the column block size NB and row block
*> size MB respectively used in SLATSQR.
*>
*> ========
*> On exit:
*> ========
*>
*> The output stored in the matrices A and T2 uses the same format
*> as the output from the SGEQRT routine stored in the matrices A and T.
*> See the documentation for SGEQRT.
*>
*> The Householder QR of a M-by-N tall-skinny matrix, Q1_out
*> where M>N is given by (with matrix dimensions):
*>
*>     Q1_out{M-by-N} = Q_in{M-by-M} * ( R_out{N-by-N} )
*>                                     ( 0{(M-N)-by-N} ),
*> where 0 is a zero matrix.
*>
*>
*> The M-by-M factor Q_out:
*> ========================
*>
*> The computed M-by-M orthogonal factor Q_out is defined implicitly as
*> a product of orthogonal matrices Q_out(i). Each Q_out(i) is stored in
*> the compact WY-representation format in the corresponding blocks of
*> matrices Y (stored in A) and T2.
*>
*> The M-by-N unit lower-trapezoidal matrix Y stored in the M-by-N
*> matrix A contains the column vectors Y(i) in NB2-size column
*> blocks YB(j). For example, YB(1) contains the columns
*> Y(1), Y(2), ... Y(NB2). NOTE: The unit entries on
*> the diagonal of Y are not stored in A.
*>
*> The number of column blocks is
*>
*>     NOCB = Number_of_output_col_blocks = CEIL(N/NB2)
*>
*> where each block is of order NB2 except for the last block, which
*> is of order LAST_NB2 = N - (NOCB-1)*NB2.
*>
*> For example, if M=6,  N=5 and NB2=2, the matrix Y is
*>
*>
*>     Y = (    YB(1),   YB(2), YB(3) ) =
*>
*>       = (   1                      )
*>         ( y21    1                 )
*>         ( y31  y32    1            )
*>         ( y41  y42  y43   1        )
*>         ( y51  y52  y53  y54    1  )
*>         ( y61  y62  y63  y54   y65 )
*>
*>
*> For each of the column blocks YB(i), an upper-triangular block
*> reflector TB2(i) is computed. These blocks are stored as
*> a sequence of upper-triangular column blocks in the NB2-by-N
*> matrix T2. The size of each TB2(i) block is NB2-by-NB2, except
*> for the last block, whose size is LAST_NB2-by-LAST_NB2.
*>
*> For example, if M=6,  N=5 and NB2=2, the matrix T2 is
*>
*>     T2 = (       TB2(1),       TB2(2), TB2(3) ) =
*>
*>        = ( t2_11  t2_12  t2_13  t2_14  t2_15  )
*>          (        t2_22         t2_24         )
*>
*>
*> The M-by-M factor Q_out is given as a product of NOCB
*> orthogonal M-by-M matrices Q_out(i).
*>
*>     Q_out = Q_out(1) * Q_out(2) * ... * Q_out(NOCB),
*>
*> where each matrix Q_out(i) is given by the WY-representation
*> using corresponding blocks from the matrices Y and T2:
*>
*>     Q_out(i) = I - YB(i) * TB2(i) * (YB(i))**T,
*>
*> where I is the identity matrix. Here is the formula with matrix
*> dimensions:
*>
*>  Q(i){M-by-M} = I{M-by-M} -
*>    YB(i){M-by-INB2} * TB2(i){INB2-by-INB2} * (YB(i))**T {INB2-by-M},
*>
*> where INB2 = NB2, except for the last block NOCB
*> for which INB2=LAST_NB2.
*>
*> The N-by-N factor R_out:
*> ========================
*>
*> The computed upper-triangular N-by-N factor R_out is stored
*> explicitly in the elements on and above the diagonal of the
*> matrix A.
*>
*> =====
*> NOTE:
*> =====
*>
*> Note that in general the N-by-N factor R_in is not equal to
*> the N-by-N R_out element-wise.
*>
*>
*> For the details of the algorithm, see [1].
*>
*> [1] "Reconstructing Householder vectors from tall-skinny QR",
*>     G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen,
*>     E. Solomonik, J. Parallel Distrib. Comput.,
*>     vol. 85, pp. 3-31, 2015.
*> \endverbatim
*>
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2019
*
*> \ingroup realGEcomputational
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*> June 2019, Igor Kozachenko,
*>            Computer Science Division,
*>            University of California, Berkeley
*>
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE SLAORHR( M, N, MB1, NB1, A, LDA, T1, LDT1, NB2, T2,
     $                    LDT2, D, W, LDW, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine (version 3.9.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2019
*
*     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDT1, LDT2, LDW, LWORK, M, N, MB1,
     $                  NB1, NB2
*     ..
*     .. Array Arguments ..
      REAL              A( LDA, * ), D( * ), T1( LDT1, * ),
     $                  T2( LDT2, * ),  W( LDW, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, JB, JBTEMP1, JBTEMP2, JNB,
     $                   LW,  NB1LOCAL, NB2LOCAL, NPLUSONE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLAORHR_GETRFNP, SLAMTSQR, SLASET,
     $                   SSCAL, STRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      LQUERY  = LWORK.EQ.-1
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. M.LT.N ) THEN
         INFO = -2
      ELSE IF( MB1.LE.N ) THEN
         INFO = -3
      ELSE IF( NB1.LT.1 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDT1.LT.MAX( 1, MIN( NB1, N ) ) ) THEN
         INFO = -8
      ELSE IF( NB2.LT.1 ) THEN
         INFO = -9
      ELSE IF( LDT2.LT.MAX( 1, MIN( NB2, N ) ) ) THEN
         INFO = -11
      ELSE IF( LDW.LT.MAX( 1, M ) ) THEN
         INFO = -14
      ELSE
*
*        Test the input LWORK, the dimension of the array WORK,
*        and set the optimal size LW for LWORK.
*
*        Set the block size for the input column blocks.
*
         NB1LOCAL = MIN( NB1, N )
*
*        WORK in SLAORHR is only used as a work array in
*        the call to SLAMTSQR, therefore the optimal size LW for WORK
*        is the same as the optimal size for WORK in SLAMTSQR.
*
         LW = N * NB1LOCAL
*
         IF( ( LWORK.LT.MAX( 1, LW ) ) .AND. (.NOT.LQUERY) ) THEN
            INFO = -16
         END IF
*
      END IF
*
*     Handle error in the input parameters and return workspace query.
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAORHR', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
         WORK( 1 ) = REAL( LW )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN( M, N ).EQ.0 ) THEN
         WORK( 1 ) = REAL( LW )
         RETURN
      END IF
*
*     Set the block size for the output column blocks.
*
      NB2LOCAL = MIN( NB2, N )
*
*     (1) Form explicitly the tall-skinny M-by-N left submatrix Q1_in
*     of M-by-M orthogonal matrix Q_in, which is implicitly stored in
*     the input arrays A and T1, by the following operation using
*     the routine SLAMTSQR.
*
*         Q1_in = Q_in * ( I ), where I is a N-by-N identity matrix,
*                        ( 0 )        0 is a (M-N)-by-N zero matrix.
*
*     (1a) Form M-by-N matrix in the array W with ones on the diagonal
*     and zeros elsewhere.
*
      CALL SLASET( 'F', M, N, ZERO, ONE, W, LDW )
*
*     (1b)  On input, W stores ( I ); on output, W stores Q1_in.
*                              ( 0 )
*
      CALL SLAMTSQR( 'L', 'N', M, N, N, MB1, NB1LOCAL, A, LDA, T1, LDT1,
     $               W, LDW, WORK, LWORK, IINFO )
*
*     (2) Perform the modified LU-decomposition C1 = Q1_in - S1 = L1*U
*     to compute the unit lower-trapezoidal Y1, since Y1 = L1. (ones
*     on the diagonal are not stored). Q1_in and S1 are M-by-N matrices,
*     U is a N-by-N matrix.
*
*     Divide C1 into the upper N-by-N submatrix C11 and the lower
*     (M-N)-by-N submatrix C21.
*
*      C1 = ( C11 )
*           ( C21 )
*
*     (2-1) Factor L11 and U11.
*
      CALL SLAORHR_GETRFNP( N, N, W, LDW, D, IINFO )
*
*     (2-2) Solve for L21.
*
      CALL STRSM( 'R', 'U', 'N', 'N', M-N, N, ONE, W, LDW,
     $            W( N+1, 1 ), LDW )
*
*     (3) Copy the result Y1 (i.e. L1) from the work array W into the
*     output array A column-by-column, and only subdiagonal elements.
*     (Y1 is a M-by-N unit lower-trapezoidal where ones on the
*     diagonal are not stored)
*
      DO J = 1, N
         CALL SCOPY( M-J, W( J+1, J ), 1, A( J+1, J ), 1 )
      END DO
*
*     (4) Reconstruct the block reflector T2 stored in T2(1:NB2, 1:N)
*     as a sequence of upper-triangular blocks with NB2-size column
*     blocking.
*
*     Loop over the column blocks of size NB2 of the array W(1:M,1:N)
*     and the array T2(1:NB2,1:N), JB is the column index of a column
*     block, JNB is the column block size at each step JB.
*
      NPLUSONE = N + 1
      DO JB = 1, N, NB2LOCAL
*
*        (4-0) Determine the column block size JNB.
*
         JNB = MIN( NPLUSONE-JB, NB2LOCAL )
*
*        (4-1) Copy the upper-triangular part of the current JNB-by-JNB
*        diagonal block U(JB) (of the N-by-N matrix U) stored
*        in W(JB:JB+JNB-1,JB:JB+JNB-1) into the upper-triangular part
*        of the current JNB-by-JNB block T2(1:JNB,JB:JB+JNB-1)
*        column-by-column, total JNB*(JNB+1)/2 elements.
*
         JBTEMP1 = JB - 1
         DO J = JB, JB+JNB-1
            CALL SCOPY( J-JBTEMP1, W( JB, J ), 1, T2( 1, J ), 1 )
         END DO
*
*        (4-2) Perform on the upper-triangular part of the current
*        JNB-by-JNB diagonal block U(JB) (of the N-by-N matrix U) stored
*        in T2(1:JNB,JB:JB+JNB-1) the following operation in place:
*        (-1)*U(JB)*S(JB), i.e the result will be stored in the upper-
*        triangular part of T2(1:JNB,JB:JB+JNB-1). This multiplication
*        of the JNB-by-JNB diagonal block U(JB) by the JNB-by-JNB
*        diagonal block S(JB) of the N-by-N sign matrix S from the
*        right means changing the sign of each J-th column of the block
*        U(JB) according to the sign of the diagonal element of the block
*        S(JB), i.e. S(J,J) that is stored in the array element D(J).
*        (NOTE: The N-by-N sign matrix S is the upper submatrix of the
*        M-by-N sign matrix S1, mentioned above).
*
         DO J = JB, JB+JNB-1
            IF( D( J ).EQ.ONE ) THEN
               CALL SSCAL( J-JBTEMP1, -ONE, T2( 1, J ), 1 )
            END IF
         END DO
*
*        (4-3) Perform the triangular solve for the current block
*        matrix X(JB):
*
*               X(JB) * (A(JB)**T) = B(JB), where:
*
*               A(JB)**T  is a JNB-by-JNB unit upper-triangular
*                         coefficient block, and A(JB)=L1(JB), which
*                         is a JNB-by-JNB unit lower-triangular block
*                         stored in W(JB:JB+JNB-1,JB:JB+JNB-1).
*                         The N-by-N matrix L is the upper part
*                         of the M-by-N lower-trapezoidal matrix L1
*                         stored in W(1:M,1:N);
*
*               B(JB)     is a JNB-by-JNB  upper-triangular right-hand
*                         side block, B(JB) = (-1)*U(JB)*S(JB), and
*                         B(JB) is stored in T2(1:JNB,JB:JB+JNB-1);
*
*               X(JB)     is a JNB-by-JNB upper-triangular solution
*                         block, X(JB) is the upper-triangular block
*                         reflector T(JB), and X(JB) is stored
*                         in T2(1:JNB,JB:JB+JNB-1).
*
*             In other words, we perform the triangular solve for the
*             upper-triangular block T(JB):
*
*               T(JB) * (L1(JB)**T) = (-1)*U(JB)*S(JB).
*
*             Even though the blocks X(JB) and B(JB) are upper-
*             triangular, the routine STRSM will access all JNB**2
*             elements of the square T2(1:JNB,JB:JB+JNB-1). Therefore,
*             we need to set to zero the elements of the block
*             T2(1:JNB,JB:JB+JNB-1) below the diagonal before the call
*             to STRSM.
*
*        (4-3a) Set the elements to zero.
*
         JBTEMP2 = JB - 2
         DO J = JB, JB+JNB-2
            DO I = J-JBTEMP2, NB2LOCAL
               T2( I, J ) = ZERO
            END DO
         END DO
*
*        (4-3b) Perform the triangular solve.
*
         CALL STRSM( 'R', 'L', 'T', 'U', JNB, JNB, ONE,
     $               W( JB, JB ), LDW, T2( 1, JB ), LDT2 )
*
      END DO
*
*     (5) Compute in place, i.e. in the array A, the upper-triangular
*         N-by-N output R_out by performing R_out = S * R_in. (S is the
*         upper N-by-N submatrix of the M-by-N sign matrix S1). This
*         multiplication by the sign matrix S on the left means
*         changing the sign of I-ths row in the matrix R_in according
*         to the sign of I-th diagonal element D(I) of the matrix S.
*
      DO I = 1, N
         IF( D( I ).EQ.-ONE ) THEN
            CALL SSCAL( NPLUSONE-I, -ONE, A( I, I ), LDA )
         END IF
      END DO
*
      WORK( 1 ) = DCMPLX( LW )
      RETURN
*
*     End of SLAORHR
*
      END