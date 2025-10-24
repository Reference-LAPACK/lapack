*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGECX computes a CX factorization of a real M-by-N matrix A:
*>
*>   A * P(K) = C*X + A_resid, where
*> 
*>   C is an M-by-K matrix which is a subset of K columns selected
*>     from the original matrix A,
*>
*>   X is a K-by-N matrix that minimizes the Frobenius norm of the
*>     residual matrix A_resid, X = pseudoinv(C) * A,
*>
*>   P(K) is an N-by-N permutation matrix chosen so that the first
*>     K columns of A*P(K) equal C,
*> 
*>   A_resid is an M-by-N residual matrix.
*>
*> The column selection for the matrix C has two stages.
*> 
*> Column selection stage 1.
*> =========================
*> 
*> The user can select N_sel columns and deselect N_desel columns
*> of the matrix A that MUST be included and excluded respectively
*> from the matrix C a priori, before running the column selection
*> algorithm. This is controlled by the flags in the array 
*> SEL_DESEL_COLS. The deselected columns are permuted to the right 
*> side of the array A and selected columns are permuted to the left
*> side of the array A. The details of the column permutation 
*> (i.e. the column permutation matrix P(K)) are stored in the 
*> array JPIV. This feature can be used when the goal is to approximate
*> the deselected columns by linear combinations of K selected columns,
*> where the K columns MUST include the N_sel selected columns.
*>
*> Column selection stage 2.
*> =========================
*>
*> The routine runs the column selection algorithm that can
*> be controlled with three stopping criteria described below.
*> For the column selection, the routine uses a truncated (rank K) in 
*> Householder QR factorization with column pivoting algorithm
*> DGEQP3RK routine. Note, that before running the column selection
*> algorithm, the user can deselect M_desel rows of the matrix A that
*> should NOT be considered by the column selection algorithm (i.e.
*> during the factorization). This is controlled by the flags in 
*> the array DESEL_ROWS. The deselected rows are permuted to the
*> bottom of the array A. The details of the row permutation (i.e. the
*> row permutation matrix) are stored in the array IPIV. This feature
*> can be used when the goal is to use the deselected rows as test data,
*> and the selected rows as training data.
*>
*> This means that the column selection factorization algorithm is 
*> effectively running on the submatrix A_sub=A(1:M_sub,1:N_sub) of
*> the matrix A after the permutations described above. Here M_sub is 
*> the number of rows of the matrix A minus the number of deselected 
*> rows M_desel, i.e. M_sub = M - M_desel, and  N_sub is the number 
*> of columns of the matrix A minus the number of deselected columns 
*> N_desel, i.e. N_sub = N - N_desel.
*>
*> Column selection criteria.
*> ==========================
*>
*> The column selection criteria (i.e. when to stop the factorization)
*> can be any of the following:
*>
*>   1) The input parameter KMAXFREE, the maximum number of columns
*>      to factorize outside of the N_sel preselected columns,
*>      i.e. the factorization rank is limited to N_sel + KMAXFREE.
*>      If N_sel + KMAXFREE >= min(M_sub, N_sub), the criterion 
*>      is not used.
*>
*>   2) The input parameter ABSTOL, the absolute tolerance for
*>      the maximum column 2-norm of the submatrix residual
*>      A_sub_resid = A(K+1:M_sub, K+1:N_sub).
*>      This means that the factorization stops if this norm is less
*>      or equal to ABSTOL. If ABSTOL < 0.0, the criterion is not used.
*>
*>   3) The input parameter RELTOL, the tolerance for the maximum
*>      column 2-norm matrix of the submatrix residual 
*>      A_sub_resid(K) = A(K+1:M_sub, K+1:N_sub) divided
*>      by the maximum column 2-norm of the submatrix 
*>      A_sub = A(1:M_sub, 1:N_sub).
*>      This means that the factorization stops when the ratio of the
*>      maximum column 2-norm of A_sub_resid to the maximum column
*>      2-norm of A_sub is less than or equal to RELTOL.
*>      If RELTOL < 0.0, the criterion is not used.
*>
*> The algorithm stops when any of these conditions is first
*> satisfied, otherwise the whole submatrix A_sub is factorized.
*>
*> For a full rank factorization of the matrix A_sub, use selection 
*> criteria that satisfy N_sel + KMAXFREE >= min(M_sub,N_sub) and 
*> ABSTOL < 0.0 and RELTOL < 0.0.
*>
*> If the user wants to verify whether the columns of the matrix C are
*> sufficiently linearly independent for their intended use, the user 
*> can compute the condition number of its R factor by calling DTRCON 
*> on the upper-triangular part of A(1:K,1:K) of the output array A.
*> 
*> How N_sel affects the column selection algorithm.
*> =================================================
*>
*> As mentioned above, the N_sel selected columns are permuted to the 
*> right side of the array A, and will be included in the column 
*> selection. Then the routine runs the factorization of that block 
*> A(1:M_sub,1:N_sel), and if any of the three stopping criteria is met
*> immediately after factoring the first N_sel columns the routine exits
*> (i.e. there is no requirement to select extra columns,
*> if the absolute or relative tolerance of the maximum column 2-norm of
*> the residual is satisfied). In this case, the number 
*> of selected columns would be K = N_sel. Otherwise, the factorization
*> routine finds a new column to select with the maximum column 2-norm 
*> in the residual A(N_sel+1:M_sub,N_sel+1:N_sub), and permutes that
*> column to the right side of A(1:M,N_sel+1:N_sub). Then the routine 
*> checks if the stopping criteria are met in the next residual
*> A(N_sel+2:M_sub,N_sel+2:N_sub), and so on.
*>
*> Computation of the matrix factors.
*> ==================================
*>
*> When the columns are selected for the factor C, and:
*>  a) If the flag RET_C is set, then the routine explicitly returns
*>     the matrix C, otherwise the routine returns only the indices of 
*>     the selected columns from the original matrix A stored in the JPIV
*>     array as the first K elements.
*>  b) If the flag COMP_X is set, then the routine also explicitly 
*>     computes and returns the factor X = pseudoinv(C) * A. 
*> 
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] FACT
*> \verbatim
*>          FACT is CHARACTER*1
*>          Specifies how the factors of CX factorization
*>          are returned.
*>
*>          = 'P' or 'p' : return only the column permutaion matrix P
*>                         in the array JPIV. The first K elements
*>                         of the array JPIV contain indeces of 
*>                         the factor C colums that were selected
*>                         from the matrix A.
*>                         (fastest, smallest memory space)
*>
*>          = 'C' or 'c' : return the column permutaion matrix P
*>                         in the array JPIV and the factor C
*>                         explicitly in the array C
*>                         (slower, more memory space)
*>
*>          = 'X' or 'x' : return the column permutaion matrix P
*>                         in the array JPIV, and both factors
*>                         C and X exlplicitly in the arrays
*>                         C and X respectively.
*>                         (slowest, largest memory space)
*> \endverbatim
*>
*> \param[in] USESD
*> \verbatim
*>          USESD is CHARACTER*1
*>          Specifies if row deselection and column 
*>          preselection-deselection functionality is turned ON or OFF.
*>
*>          = 'N' or 'n' : Both row deselection and column
*>                         preselection-deselection are OFF. 
*>                         Both arrays DESEL_ROWS and
*>                         SEL_DESEL_COLS are not used.
*>
*>          = 'R' or 'r' : Only row deselection is ON. 
*>                         Column preselection-deselection is OFF.
*>                         The array SEL_DESEL_COLS is not used.
*>
*>          = 'C' or 'c' : Only column preselection-deselection is ON.
*>                         Row deselection is OFF. 
*>                         The array DESEL_ROWS is not used.
*>
*>          = 'A' or 'a' : Means "All".
*>                         Both row deselection and column
*>                         preselection-deselection are ON. 
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] DESEL_ROWS
*> \verbatim
*>          DESEL_ROWS is INTEGER array, dimension (M)
*>          This is a row deselection mask array that separates
*.          the matrix A rows into 2 sets.
*>
*>          a) If DESEL_ROWS(I) = -1, the I-th row of the matrix A is
*>             deselected by the user, i.e. chosen to be excluded from 
*.             the algorithm and will be permuted to the bottom of A.
*>             The number of deselected rows is denoted by M_desel.         
*>         
*>          b) If DESEL_ROWS(I) not equal -1,
*>             the I-th row of A is a free row and will be used by the
*>             algorithm. This defines a set of M_sub = M - M_desel
*>             rows that the algorithm will work on. After permutation,
*>             this set will be in the top of the matrix A.
*> \endverbatim
*>
*> \param[in] SEL_DESEL_COLS
*> \verbatim
*>          SEL_DESEL_COLS is INTEGER array, dimension (N)
*>          This is a column preselection/deselection mask array that
*.          separates the matrix A columns into 3 sets.
*>
*>          a) If SEL_DESEL_COLS(J) = +1, the J-th column of the matrix 
*>             A is selected by the user to be included in the factor C
*>             and will be permuted to the left side of the array A.
*>             The number of selected columns is denoted by N_sel.
*>
*>          b) If SEL_DESEL_COLS(J) = -1, the J-th column of the matrix
*>             A is deselected by the user, i.e. chosen to be excluded
*>             from the factor C and will be permuted to the right side
*>             of the array A. The number of deselected columns is
*>             denoted by N_desel.          
*>         
*>          c) If SEL_DESEL_COLS(J) not equal 1, and not equal -1,
*>             the J-th column of A is a free column and will be used by
*>             the algorithm to determine if this column has to be 
*>             selected. This defines a set of 
*>             N_free = N - N_sel - N_desel.
*>         
*>          NOTE: Error returned as INFO = -6 means that the number of
*>          preselected N_sel colunms is larger than M_sub. 
*>          Therefore, the factoriaztion of all N_sel preselected
*>          columns cannot be completed.           
*> \endverbatim
*>
*> \param[in] KMAXFREE
*> \verbatim
*>          KMAXFREE is INTEGER
*>
*>          The first column selection stopping criterion in the 
*>          column selection stage 2.
*>
*>          The maximum number of columns of the matrix A_sub to select
*>          during the factorization stage, KMAXFREE >= 0 
*>          
*>          KMAXFREE does not include the preselected columns.
*>          N_sel + KMAXFREE is the maximum factorization rank of
*>          the matrix A_sub = A(1:M_sub, 1:N_sub).
*>
*>          a) If N_sel + KMAXFREE >= min(M_sub, N_sub), then this 
*>                stopping criterion is not used, i.e. columns are selected
*>                in the factorization stage depending on 
*>                ABSTOL and RELTOL.
*>
*>          b) If KMAXFREE = 0, then this stopping criterion is
*>                satisfied on input and the routine exits without 
*>                performing column selection stage 2 on the submatrix 
*>                A_sub. This means that the matrix
*>                A_free = A(N_sel+1:M_sub, N_sel+1:N_sub) is not modified.
*>                and A_free is itself the residual for the factorization.
*> \endverbatim
*>
*> \param[in] ABSTOL
*> \verbatim
*>          ABSTOL is DOUBLE PRECISION, cannot be NaN.
*>
*>          The second column selection stopping criterion in the 
*>          column selection stage 2.
*>
*>          Here, SAFMIN = DLAMCH('S').
*>
*>          The absolute tolerance (stopping threshold) for
*>          maximum column 2-norm of the residual matrix
*>          A_sub_resid(K) = A_sub(K+1:M_sub, K+1:N_sub),
*>          when K columns were factorized.
*>          The algorithm converges (stops the factorization) when
*>          the maximum column 2-norm of the residual matrix
*>          A_sub_resid is less than or equal to ABSTOL.
*>
*>          a) If ABSTOL is NaN, then no computation is performed
*>                and an error message ( INFO = -8 ) is issued
*>                by XERBLA.
*>
*>          b) If ABSTOL < 0.0, then this stopping criterion is not
*>                used, factorize columns depending on KMAXFREE
*>                and RELTOL.
*>                This includes the case ABSTOL = -Inf.
*>      
*>          c) If 0.0 <= ABSTOL < 2*SAFMIN, then ABSTOL = 2*SAFMIN
*>                is used. This includes the case ABSTOL = -0.0.
*>
*>          d) If 2*SAFMIN <= ABSTOL then the input value
*>                of ABSTOL is used.
*>
*>          Here, maxcol2norm(A_free) is the maximum column 2-norm
*>          of the matrix A_free = A(N_sel+1:M_sub, N_sel+1:N_sub).
*>
*>          If ABSTOL chosen above is >= maxcol2norm(A_free), then
*>          this stopping criterion is satisfied after the matrix
*>          A_sel = A(1:M_sub, 1:N_sel) is factorized and the 
*>          routine exits immediately after maxcol2norm(A_free) is
*>          computed to return it in MAXC2NORMK. This means that 
*>          the factorization residual
*>          A_sub_resid = A_free = A(N_sel+1:M_sub, N_sel+1:N_sub)
*>          is not modified.
*>          Also RELMAXC2NORMK of A_free is returned.
*>          This includes the case ABSTOL = +Inf.          
*> \endverbatim
*>
*> \param[in] RELTOL
*> \verbatim
*>          RELTOL is DOUBLE PRECISION, cannot be NaN.
*>
*>          The third column selection stopping criterion in the 
*>          column selection stage 2.
*>
*>          Here, EPS = DLAMCH('E').
*>
*>          The tolerance (stopping threshold) for the ratio
*>          maxcol2norm(A_sub_resid(K))/maxcol2norm(A_sub) of
*>          the maximum column 2-norm of the residual matrix 
*>          A_sub_resid(K) = A_sub(K+1:M_sub, K+1:N_sub) and 
*>          the maximum column 2-norm of the original submatrix
*>          A_sub = A(1:M_sub, 1:N_sub). The algorithm
*>          converges (stops the factorization), when
*>          maxcol2norm(A_sub_resid(K))/maxcol2norm(A_sub) is 
*>          less than or equal to RELTOL.
*>
*>          a) If RELTOL is NaN, then no computation is performed
*>                and an error message ( INFO = -9 ) is issued
*>                by XERBLA.
*>
*>          b) If RELTOL < 0.0, then this stopping criterion is not
*>                used, factorize columns depending on KMAXFREE 
*>                and ABSTOL.
*>                This includes the case RELTOL = -Inf.
*>
*>          c) If 0.0 <= RELTOL < EPS, then RELTOL = EPS is used.
*>                This includes the case RELTOL = -0.0.
*>
*>          d) If EPS <= RELTOL then the input value of RELTOL
*>                is used.
*>
*>          If RELTOL chosen above is >= 1.0, then this stopping
*>          criterion is satisfied on input and routine exits
*>          immediately after A_sel = A(1:M_sub, 1:N_sel))
*>          is factorized and maxcol2norm(A_free) is computed to
*>          return it in MAXC2NORMK. This means that 
*>          the factorization residual 
*>          A_sub_resid = A_free = A(N_sel+1:M_sub, N_sel+1:N_sub)
*>          is not modified.
*>          Also RELMAXC2NORMK is returned as 1.0.
*>          This includes the case RELTOL = +Inf.
*>
*>          NOTE: We recommend RELTOL to satisfy
*>                min(max(M_sub,N_sub)*EPS, sqrt(EPS)) <= RELTOL
*>
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>
*>          On entry:
*>            the M-by-N matrix A.
*>
*>          On exit:
*>            NOTE DEFINITIONS: M_sub = M_free,
*>                              N_sub = N_sel + N_free
*>
*>            The output parameter K, the number of selected columns,
*>            is described later. 
*>
*>            1) If K = 0, A(1:M,1:N) contains the original matrix A.
*>
*>            2) If K > 0, A(1:M,1:N): contains the following parts:
*>                  
*>            (a) If M_sub < M (which is the same as M_desel > 0),
*>               the subarray A(M_sub+1:M,1:N) contains the deselected
*>               rows.
*> 
*>            (b) If N_sub < N ( which is the same as N_desel > 1 ).
*>               the subarray A(1:M,N_sub+1:N) contains the
*>               deselected columns.
*>
*>            (c) If N_sel > 0,
*>               the union of the subarray A(1:M_sub, 1:N_sel)
*>               and the subarray A(1:N_sel, 1:N_sub) contains parts
*>               of the factors obtained by computing Householder QR 
*>               factorization WITHOUT column pivoting of N_sel
*>               preselected columns using DGEQRF routine.
*>               
*>            (d) The subarray A(N_sel:M_sub, N_sel:N_sub) contains 
*>               parts of the factors obtained by computing a truncated
*>               (rank K) Householder QR factorization with
*>               column pivoting using DGEQP3RK on the matrix
*>               A_free = A(N_sel+1:M_sub, N_sel+1:N_sub) which
*>               is the result of applying selection and deselection
*>               of columns, applying deselection of rows to the 
*>               original matrix A, and applying orthogonal
*>               transformation from the factorization of the first
*>               N_sel columns as described in part (c).
*>
*>              1. The elements below the diagonal of the subarray
*>                 A_sub(1:M_sub,1:K) together with TAU(1:K)
*>                 represent the orthogonal matrix Q(K) as a
*>                 product of K Householder elementary reflectors.
*>
*>              2. The elements on and above the diagonal of
*>                 the subarray A_sub(1:K,1:N_sub) contain 
*>                 K-by-N_sub upper-trapezoidal matrix
*>                 R_sub_approx(K) = ( R_sub11(K), R_sub12(K) ).
*>                 NOTE: If K=min(M_sub,N_sub), i.e. full rank 
*>                       factorization, then R_sub_approx(K) is the
*>                       full factor R which is upper-trapezoidal.
*>                       If, in addition, M_sub>=N_sub, then R is
*>                       upper-triangular.
*>
*>              3. The subarray A_sub(K+1:M_sub,K+1:N_sub) contains
*>                 (M_sub-K)-by-(N_sub-K) rectangular matrix
*>                 A_sub_resid(K).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] K
*> \verbatim
*>          K is INTEGER
*>          The number of columns that were selected.
*>          (K is the factorization rank)
*>          0 <= K <= min( M_sub, min(N_sel+KMAXFREE, N_sub) ).
*>
*>          If K = 0, the arrays A, TAU were not modified.
*> \endverbatim
*>
*> \param[out] MAXC2NRMK
*> \verbatim
*>          MAXC2NRMK is DOUBLE PRECISION
*>          The maximum column 2-norm of the residual matrix
*>          A_sub_resid(K) = A_sub(K+1:M_sub, K+1:N_sub),
*>          when factorization stopped at rank K. MAXC2NRMK >= 0.
*>
*>          a) If K = 0, i.e. the factorization was not performed,
*>             the matrix A_sub = A(1:M_sub, 1:N_sub) was not modified
*>             and is itself a residual matrix, then MAXC2NRMK equals
*>             the maximum column 2-norm of the original matrix A_sub.
*>
*>          b) If 0 < K < min(M_sub, N_sub), then MAXC2NRMK is returned.
*>
*>          c) If K = min(M_sub, N_sub), i.e. the whole matrix A_sub was
*>             factorized and there is no factorization residual matrix,
*>             then MAXC2NRMK = 0.0.
*>
*>          NOTE: MAXC2NRMK at the factorization step K would equal
*>                to the diagonal element R_sub(K+1,K+1) of the factor
*>                R_sub in the next factorization step K+1.
*> \endverbatim
*>
*> \param[out] RELMAXC2NRMK
*> \verbatim
*>          RELMAXC2NRMK is DOUBLE PRECISION
*>          The ratio MAXC2NRMK / MAXC2NRM of the maximum column
*>          2-norm of the residual matrix 
*>          A_sub_resid(K) = A_sub(K+1:M_sub, K+1:N_sub) (when
*>          factorization stopped at rank K) and maximum column 2-norm 
*>          MAXC2NRM of the matrix A_sub = A(1:M_sub, 1:N_sub).
*>          RELMAXC2NRMK >= 0.
*>
*>          a) If K = 0, i.e. the factorization was not performed,
*>             the matrix  A_sub was not modified
*>             and is itself a residual matrix,
*>             then RELMAXC2NRMK = 1.0.
*>
*>          b) If 0 < K < min(M_sub,N_sub), then
*>                RELMAXC2NRMK = MAXC2NRMK / MAXC2NRM is returned.
*>
*>          c) If K = min(M_sub,N_sub), i.e. the whole matrix A_sub was
*>             factorized and there is no residual matrix
*>             A_sub_resid(K), then RELMAXC2NRMK = 0.0.
*>
*>         NOTE: RELMAXC2NRMK at the factorization step K would equal
*>               abs(R_sub(K+1,K+1))/MAXC2NRM in the next
*>               factorization step K+1, where R_sub(K+1,K+1) is the 
*>               diaginal element of the factor R_sub in the next
*>               factorization step K+1.
*> \endverbatim
*>
*> \param[out] FNRMK
*> \verbatim
*>          FNRMK is DOUBLE PRECISION
*>          Frobenius norm of the factorization residual matrix 
*>          A_sub_resid(K) = A_sub(K+1:M_sub, K+1:N_sub).
*>          FNRMK >= 0.0
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (M)
*>          Row permutation indices due to row
*>          deselection, for 1 <= i <= M.
*>          If IPIV(i)= k, then the row i of A_sub was the
*>          the row k of A.
*> \endverbatim
*>
*> \param[out] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (N)
*>          Column permutation indices, for 1 <= j <= N.
*>          If JPIV(j)= k, then the column j of A*P was the
*>          the column k of A.
*> 
*>          The first K elements of the array JPIV contain
*>          indeces of the factor C colums that were selected
*>          from the matrix A.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors.
*>
*>          If 0 < K <= MIN(M_sub,N_sub), only elements TAU(1:K) of
*>          the array TAU may be modified. The elements
*>          TAU(K+1:min(M_sub,N_sub)) are set to zero.
*>          The elements of TAU(min(M_sub,N_sub)+1:N) are not
*>          modified.
*> \endverbatim
*>
*> \param[out] C
*> \verbatim
*>          C is DOUBLE PRECISION array.
*>          If FACT = 'C' or 'X':
*>             If K > 0, C is the M-by-K factor C
*>                       and array has dimension (LDC,N),
*>          If FACT = 'N':
*>             array is not used and can have linear dimension >=1.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C.
*>          If FACT = 'C' or 'X', LDC >= max(1,M).
*>          If FACT = 'P', LDC >= 1.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is DOUBLE PRECISION array.
*>          If FACT = 'X':
*>             If K > 0, C is the K-by-N factor X
*>                       and array has dimension (LDX,N).
*>          If FACT = 'P':
*>             array is not used and can have linear dimension >=1.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X.
*>          If FACT = 'X', LDC >= max(1,M).
*>          If FACT = 'P', LDC >= 1.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO>=0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= 3*N_sub+1.
*>          For optimal performance LWORK >= 2*N_sub+( N_sub+1 )*NB,
*>          where NB is the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N).
*>          Is a work array. ( IWORK is used by DGEQP3RK to store indices
*>          of "bad" columns for norm downdating in the residual
*>          matrix in the blocked step auxiliary subroutine DLAQP3RK ).
*> \endverbatim
*>
*> \param[out] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          The dimension of the array LIWORK. LIWORK >= N
*>
*>          If LIWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the IWORK array, and no error
*>          message related to LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit.
*>          < 0: if INFO = -i, the i-th argument had an illegal value.
*>          > 0: if INFO =  i, the i-th diagonal element of the
*>                triangular factor of in the matrix C is zero,
*>                so that C does not have full rank; X cannot be
*>                computed as  the least squares solution to C*X = A.
*> \endverbatim
*  =====================================================================      
      SUBROUTINE         DGECX( FACT, USESD, M, N,
     $                          DESEL_ROWS, SEL_DESEL_COLS,
     $                          KMAXFREE, ABSTOL, RELTOL, A, LDA,
     $                          K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,       
     $                          IPIV, JPIV, TAU, C, LDC, X, LDX,
     $                          WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--     
*
*     .. Scalar Arguments ..
      CHARACTER           FACT, USESD
      INTEGER             INFO, K, KMAXFREE, LDA, LDC, LDX, LIWORK,
     $                    LWORK, M, N
      DOUBLE PRECISION    ABSTOL, ABSTOLFREE, MAXC2NRMK, RELTOL, 
     $                    RELTOLFREE, RELMAXC2NRMK, FNRMK
*     ..
*     .. Array Arguments ..      
      INTEGER             DESEL_ROWS( * ), IPIV( * ), JPIV( * ),
     $                    SEL_DESEL_COLS( * ), IWORK( * )       
      DOUBLE PRECISION    A( LDA, * ), C( LDC, * ), TAU( * ),
     $                    X( LDX, *), WORK( * )
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            INB
      PARAMETER          ( INB = 1 )
      DOUBLE PRECISION   ZERO, TWO, MINUSONE
      PARAMETER          ( ZERO = 0.0D+0, TWO = 2.0D+0,
     $                   MINUSONE = -1.0D+0 ) 
*     ..
*     .. Local Scalars ..
      LOGICAL            LIQUERY, LQUERY, RETURNC, RETURNX,
     $                   USE_DESEL_ROWS, USE_SEL_DESEL_COLS, USETOL
      INTEGER            I, J, NSUB, MFREE, MSUB, MNSUB, NSEL, JDESEL,
     $                   ITEMP, IINFO, KP, KP0, KFREE, MRESID, NRESID,
     $                   NRHS, LWKMIN, LWKOPT, JP, JJ, JPW, MINMN,
     $                   NBOPT
      DOUBLE PRECISION   EPS, MAXC2NRM, SAFMIN, RELMAXC2NRMKFREE

*     .. External Subroutines ..
      EXTERNAL           DLACPY, DGELS, DGEQP3RK, DGEQRF, DORMQR,
     $                   DSWAP, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN, LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE, DNRM2
      EXTERNAL           DISNAN, DLAMCH, DLANGE, DNRM2, IDAMAX,
     $                   ILAENV, LSAME    
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN      
*     ..
*     .. Executable Statements ..      
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      LIQUERY = ( LIWORK.EQ.-1 )
*          
      RETURNX = LSAME( FACT, 'X' ) 
      RETURNC = LSAME( FACT, 'C' ) .OR. RETURNX
*
      USE_DESEL_ROWS = LSAME( USESD, 'R' ) .OR. LSAME( USESD, 'A' ) 
      USE_SEL_DESEL_COLS = LSAME( USESD, 'C') .OR. LSAME( USESD, 'A' )
*      
      IF ( .NOT.(RETURNC .OR. LSAME( FACT, 'P') ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( USE_DESEL_ROWS .OR. USE_SEL_DESEL_COLS 
     $       .OR. LSAME( USESD, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KMAXFREE.LT.0 ) THEN
         INFO = -7
      ELSE IF( DISNAN( ABSTOL ) ) THEN
         INFO = -8
      ELSE IF( DISNAN( RELTOL ) ) THEN
         INFO = -9
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -11 
      ELSE IF( ( RETURNC .AND. LDC.LT.MAX( 1, M ) ) .OR.
     $    ( .NOT.RETURNC .AND. LDC.LT.1 )) THEN
         INFO = -20
      ELSE IF( ( RETURNX .AND. LDX.LT.MAX( 1, M ) ) .OR.
     $    ( .NOT.RETURNX .AND. LDX.LT.1 )) THEN
         INFO = -22
      END IF
*      
*     If the above input parameters are valid:       
*       a) Test the input workspace size LWORK for the minimum
*          size requirement LWKMIN. 
*       b) Determine the optimal block size NB and optimal
*          workspace size LWKOPT to be returned in WORK(1)
*          in case of: (1) LWORK < LWKMIN, (2) LQUERY = .TRUE.,
*          (3) when routine exits.
*     Here, LWKMIN is the miminum workspace required for unblocked
*     code.
*
      IF( INFO.EQ.0 ) THEN
         MINMN = MIN( M, N )
         IF( MINMN.EQ.0 ) THEN
            LWKMIN = 1
            LWKOPT = 1
         ELSE
            LWKMIN = 3*N + 1
*
*           Assign to NBOPT optimal block size.
*
            NBOPT = ILAENV( INB, 'DGEQRF', ' ', M, N, -1, -1 )
            LWKOPT =  1000
         END IF
         WORK( 1 ) = DBLE( LWKOPT )
*
         IF( ( LWORK.LT.LWKMIN ) .AND. .NOT.LQUERY ) THEN
            INFO = -24
         END IF
      END IF    
*
*     ==================================================================
*    
      K = 0
*
      EPS = DLAMCH('Epsilon')
*
      USETOL = .FALSE.
*
*     Adjust ABSTOL only if nonnegative. Negative value means disabled.
*     We need to keep negtive value for later use in criterion
*     check.
*
      IF( ABSTOL.GE.ZERO ) THEN
         SAFMIN = DLAMCH('Safe minimum')
         ABSTOL = MAX( ABSTOL, TWO*SAFMIN )
         USETOL = .TRUE.
      END IF
*
*     Ajust RELTOL only if nonnegative. Negative value means disabled.
*     We need to keep negtive value for later use in criterion
*     check.
*
      IF( RELTOL.GE.ZERO ) THEN
         RELTOL = MAX( RELTOL, EPS )
         USETOL = .TRUE.
      END IF
*
*     ==================================================================
*
*     If we need to return factor C, copy the original unctouched matrix
*     A into the array C.
*
      IF( RETURNC ) THEN
         CALL DLACPY( 'F', M, N, A, LDA, C, LDC )
      END IF
*
*     If we need to return factor X, copy the original unctouched matrix
*     A into the array X.
*
      IF( RETURNX ) THEN
         CALL DLACPY( 'F', M, N, A, LDA, X, LDX )
      END IF
*      
*     ==================================================================
*     Permute the deselected rows to the bottom of the matrix A.
*     1) Order of free rows is preserved.
*     2) Order of deselected rows is not preserved.
*     ==================================================================
*
*     I is the index of DESEL_ROWS array and row I
*     of the matrix A.
*     MFREE is the number of free rows, also the pointer to the last
*     free row.
*     (For each position I, we check if this position is a FREE row.
*     If it is a FREE row we increment the MFREE pointer, otherwise we 
*     do not change the MFREE pointer. Also, if it is a FREE row, we move
*     this row from the larger (or same) I index into samaller (or same)
*     MFREE index. This way we move all the FREE rows to the lower index
*     block preserving FREE row order. Deselected rows will be )
*      
      IF( USE_DESEL_ROWS ) THEN
*      
         MFREE = 0 
         DO I = 1, M, 1
*
*           Initialize row pivot array IPIV.          
            IPIV( I ) = I
*
            IF( DESEL_ROWS(I).NE.-1 ) THEN
               MFREE = MFREE + 1
*
*              This is the check whether the deselected row is
*              on the deselected place already.
*            
               IF( I.NE.MFREE ) THEN
*
*                 Here, we swap A(I,1:N) into A(MFREE,1:N)
*
                  CALL DSWAP( N, A( I, 1 ), LDA, A( MFREE, 1 ), LDA )                   
                  IPIV( I ) = IPIV( MFREE )
                  IPIV( MFREE ) = I
                  ITEMP = DESEL_ROWS( I )
                  DESEL_ROWS( I ) = DESEL_ROWS( MFREE )
                  DESEL_ROWS( MFREE ) = ITEMP
               END IF                       
            END IF
*            
         END DO      
*         
      ELSE
*
*        We do not row deselection DESEL_ROWS array.
*        Initialize row pivot array IPIV.
*
         DO I = 1, M, 1
            IPIV( I ) = I
         END DO
*
         MFREE = M
      END IF
      MSUB = M
*            
*     ==================================================================
*     Permute the pseselected columns to the left and deselected 
*     columns to the right of the matrix A.
*     1) Order of preselected columns is preserved.
*     2) Order of free columns is not preserved.
*     3) Order of deselected columns is not preserved.
*     ==================================================================
*
*     J is the index of SEL_DESEL_COLS array and column J
*     of the matrix A.
*
*     Column selection.
*     NSEL is the number of selected columns, also the pointer to the last
*     selected column.
*
      NSEL = 0
      IF( USE_SEL_DESEL_COLS ) THEN
*                        
         DO J = 1, N, 1
*         
*           Initialize column pivot array JPIV.         
            JPIV( J ) = J
*             
            IF( SEL_DESEL_COLS(J).EQ.1 ) THEN
               NSEL = NSEL + 1
*
*              This is the check whether the selected column is
*              on the selected place already.
*            
               IF( J.NE.NSEL ) THEN
*
*                 Here, we swap the column A(1:M,J) into A(1:M,NSEL)
*
                  CALL DSWAP( M, A( 1, J ), 1, A( 1, NSEL ), 1 )                     
                  JPIV( J ) = JPIV( NSEL )
                  JPIV( NSEL ) = J
                  SEL_DESEL_COLS( J ) = SEL_DESEL_COLS( NSEL )
                  SEL_DESEL_COLS( NSEL ) = 1
               END IF
            END IF
         END DO
*
*        Column deselection.
*
         JDESEL = N+1
         DO J = N, NSEL+1, -1
            IF( SEL_DESEL_COLS( J ).EQ.-1 ) THEN
               JDESEL = JDESEL - 1
*
*              This is the check whether the deselected column is
*              on the deselected place already.
*            
               IF( J.NE.JDESEL ) THEN
*
*                 Here, we swap the column A(1:M,J) into A(1:M,JDESEL)
*
                  CALL DSWAP( M, A( 1, J ), 1, A( 1, JDESEL ), 1 )
                  ITEMP = JPIV( J )                   
                  JPIV( J ) = JPIV( JDESEL )
                  JPIV( JDESEL ) = ITEMP
                  SEL_DESEL_COLS( J ) = SEL_DESEL_COLS( JDESEL )
                  SEL_DESEL_COLS( JDESEL ) = -1
               END IF        
            END IF
         END DO
*         
         NSUB = JDESEL - 1
*
      ELSE
*
*        We do not column selection deselection SEL_DESEL_COLS array.
*        Initialize column pivot array JPIV.
*
         DO J = 1, N, 1
            JPIV( J ) = J
         END DO
*
         NSUB = N            
      END IF
*    
*     ==================================================================
*     Compute the complete column 2-norms of the submatrix 
*     A_sub=A(1:MSUB, 1:NSUB) and store them in WORK(NSUB+1:2*NSUB).
*
      DO J = 1, NSUB
         WORK( NSUB+J ) = DNRM2( MSUB, A( 1, J ), 1 )
      END DO
*      
*     Compute the column index and the maximum column 2-norm
*     for the submatrix A_sub=A(1:MSUB, 1:NSUB).
*
      KP0 = IDAMAX( NSUB, WORK( NSUB+1 ), 1 )
      MAXC2NRM = WORK( NSUB + KP0 )
*      
*     ==================================================================
*     Process preselected columns
*
*     Compute the QR factorization of NSEL preselected columns (1:NSEL)
*     the submatrix A_sub=(1:MSUB, 1:NSUB) and update
*     remaining NFEE free columns (NSEL+1:NSUB).
*     MSUB = MFREE, NSUB = MSEL + NFREE 
*
      MNSUB = MIN(MSUB, NSUB)
      MRESID = MSUB-NSEL
      NRESID = NSUB-NSEL
      IF( NSEL.GT.0 ) THEN
         IF( MSUB.LT.NSEL ) THEN
*           TODO: Move this part to the top of the routine. 
*           a) Case MSUB < NSEL. 
*           When the number of preselected columns
*           is larger than MSUB, hence the factorization of all NSEL 
*           columns cannot be completed. Return from the routine with the
*           error of COL_SEL_DESEL parameter. NSEL cannot be larger than
*           MSUB.
*          
            INFO = -6
            WORK( 1 ) = DBLE( LWKOPT )
            RETURN
         ELSE IF( MSUB.EQ.NSEL.OR. 
     $                ( MSUB.GT.NSEL.AND.NSEL.EQ.NSUB )) THEN
*
*           b) Case MSUB = NSEL.
*           c-1) Case MSUB > NSEL and NSEL = NSUB.
*
*           There will be no residual submatrix after factorization 
*           of NSEL columns at step K = NSEL:
*           A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB).
*           Therefore, ther is no need to do the factorization of NSEL
*           columns. Set norms to ZERO and return from the routine.
*
            K = NSEL
            MAXC2NRMK = ZERO
            RELMAXC2NRMK = ZERO
            FNRMK = ZERO
*            
            DO J = K + 1, MNSUB
               TAU( J ) = ZERO
            END DO
*
*           Factorization is done. Go to computation of the factor C.
* 
            GO TO 10  
         ELSE
*
*           (c-2) Case MSUB > NSEL and NSEL < NSUB.
*
*           There is a submatrix residual at step K=NSEL
*           A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB)
*
            CALL DGEQRF( MSUB, NSEL, A, LDA, TAU, WORK, LWORK, IINFO )
*
*           Apply Q**T from the left to A(NSEL+1:MSUB, NSEL+1:NSUB)            
*
            CALL DORMQR( 'Left', 'Transpose', MSUB, NSUB-NSEL, NSEL,
     $        A, LDA, TAU, A( 1, NSEL+1 ), LDA, WORK, LWORK, IINFO )
*
*           Compute the complete column 2-norms of the submatrix
*           residual at step NSEL
*           A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB) and
*           store them in WORK(NSUB+NSEL+1:2*NSUB).
*
            DO J = NSEL+1, NSUB
               WORK( NSUB+J ) = DNRM2( MRESID, A( NSEL+1, J ), 1 )
            END DO
*      
*           Compute the column index and the maximum column 2-norm
*           and the relative maximum column 2-norm for the submatrix
*           residual.
*
            KP =  IDAMAX( NRESID, WORK( NSUB+NSEL+1 ), 1 )
*
            K = NSEL
            MAXC2NRMK = WORK( NSUB + NSEL + KP )
            RELMAXC2NRMK = MAXC2NRMK / MAXC2NRM
*
*           Test for the first, second and third tolerance stopping
*           criteria after factorizarion of preselected columns.
*           If any of them is met, return. Otherwise,
*           proceed with factorization of the NFREE free columns.
*           NOTE: There is no need to test for ABSTOL.GE.ZERO, since
*           MAXC2NRMK is non-negative. Similarly, there is no need
*           to test for RELTOL.GE.ZERO, since RELMAXC2NRMK is
*           non-negative.            
*                                        
            IF( KMAXFREE.EQ.0 
     $        .OR. MAXC2NRMK.LE.ABSTOL
     $        .OR. RELMAXC2NRMK.LE.RELTOL ) THEN
*
*              NOTE: In this (c-2) case. There is a submatrix
*              residual A_sub_resid(NSEL). We do not need to have a check
*              for  MIN(MRESID, NRESID) = 0 to call DLANGE.
*
               FNRMK = DLANGE( 'F', MRESID, NRESID, A(NSEL+1,NSEL+1),
     $                         LDA, WORK )
*
               DO J = K + 1, MNSUB
                  TAU( J ) = ZERO
               END DO
*
*              Factorization is done. Go to computation of the factor C.
*                
               GO TO 10              
            END IF
*
*
*
         END IF   
      END IF
*      
*     ==================================================================
*
*     Factorize NFREE free columns of
*     A_free = A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB),
*     KFREE is the number of columns that were actually factorized among
*     NFREE columns.
*
*     Disable RELTOLFREE when calling DGEQP3RK for free columns 
*     factorization, since it expects RELTOLFREE with respect to 
*     the residual matrix A_sub_resid(NSEL), not the whole original
*     marix A. We can use RELTOL criterion by passing it to
*     ABSTOLFREE as RELTOL*MAXC2NRM. We need to make sure that
*     the negative values of ABSTOL and RELTOL are propagated
*     to ABSTOLFREE and RELTOLFREE, since negative vaslues means
*     that the criterionis is disabled.
*
      IF( USETOL ) THEN
         ABSTOLFREE = MAX( ABSTOL, RELTOL * MAXC2NRM )
      ELSE
         ABSTOLFREE = MINUSONE
      END IF
      RELTOLFREE = MINUSONE
      NRHS = 0
*
      CALL DGEQP3RK( MRESID, NRESID, NRHS, KMAXFREE,
     $               ABSTOLFREE, RELTOLFREE,
     $               A( K+1, K+1 ), LDA, KFREE, MAXC2NRMK, 
     $               RELMAXC2NRMKFREE, JPIV( K+1 ), TAU( K+1 ),
     $               WORK, LWORK, IWORK, IINFO )
*
*     1) Adjust the return value for the number of factorized 
*        columns K for the whole submatrix A_sub.
*     2) MAXC2NRMK is returned transparently without change 
*        as it is returned from DGEQP3RK.
*     3) Adjust the return value RELMAXC2NRMK for the whole
*        submatrix A_sub. We do not use RELMAXC2NRMKFREE
*        returned from DGEQP3RK.
*
      K = K + KFREE
      RELMAXC2NRMK =  MAXC2NRMK / MAXC2NRM 
*
*     Now, MRESID and NRESID is the number of rows and columns
*     respectively in  A_free_resid = A(K+1:MSUB,K+1:NSUB).
*      
      MRESID = MRESID-KFREE
      NRESID = NRESID-KFREE
      IF( MIN( MRESID, NRESID ).NE.0 ) THEN
         FNRMK = DLANGE( 'F',  MRESID, NRESID, A( K+1, K+1 ),
     $                   LDA, WORK )
      ELSE
         FNRMK = ZERO
      END IF 
*
*     Compute the factor C.
*
   10 CONTINUE 
*   
      IF( RETURNC .AND. K.GT.0 ) THEN
*
*        Apply interchanges to columns 1:K in the matrix C in place,
*        which stores the original matrix A.
*        IWORK is used to keep track of original column indices,
*        when swaping.

         DO J = 1, N, 1
            IWORK( J ) = J
         END DO
         DO J = 1, K, 1
            JP = JPIV( J )
            IF( J.NE.JP ) THEN
               DO JJ = J, N, 1
                  IF( JP.EQ.IWORK( JJ ) ) THEN
                     JPW = JJ
                  END IF
               END DO
               IF( J.NE.JPW ) THEN
                  CALL DSWAP( M, C( 1, J ), 1, C( 1, JPW ), 1 )
                  ITEMP = IWORK( J )
                  IWORK( J ) = IWORK( JPW )
                  IWORK( JPW ) = ITEMP
               END IF
            END IF  
         END DO 
*
      END IF
*
*     Return matrix X.
*
      IF( RETURNX .AND. K.GT.0 ) THEN
*
*        We need to use C and A to compute X = pseudoinv(C) * A, as
*        the Linear Least Squares problem C*X = A. We use LLS routine
*        that uses QR factorization.  For that purpose, we store C into
*        WORK array WORK(1:M*K), and the matrix A was copied into
*        the array X at the begining of the routine.  
*
*        Copy matrix C into work array WORK.
*
         CALL DLACPY( 'F', M, K, C, LDC, WORK, M )
*
         CALL DGELS( 'N', M, K, N, WORK, M, X, LDX,
     $               WORK( M*K+1 ), LWORK,
     $               IINFO )
         INFO = IINFO
*         
      END IF
*
      WORK( 1 ) = DBLE( LWKOPT ) 
*
*     DGECX
*
      END 