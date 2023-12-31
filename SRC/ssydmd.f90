      SUBROUTINE SSYDMD( JOBS, JOBZ, JOBR, JOBF,           &
                         WHTSVD,  WHTSYM, WHTEIG,          &
                         M, N, X, LDX, Y, LDY, NRNK, TOL,  &
                         K, EIGS,  Z, LDZ,  RES,           &
                         B, LDB, W,  LDW,   S, LDS,        &
                         WORK, LWORK, IWORK, LIWORK, INFO ) 
!.....
      USE                   iso_fortran_env
      IMPLICIT NONE
      INTEGER, PARAMETER :: WP = real32
!.....
!     Scalar arguments
      CHARACTER, INTENT(IN)   :: JOBS,   JOBZ,  JOBR,  JOBF
      INTEGER,   INTENT(IN)   :: WHTSVD, WHTSYM, WHTEIG,    &
                                 M,   N, LDX,    LDY,       &
                                 NRNK,   LDZ, LDB, LDW, LDS,&
                                 LWORK,  LIWORK
      INTEGER,   INTENT(OUT)  :: K, INFO
      REAL(KIND=WP), INTENT(IN)  :: TOL
!     Array arguments      
      REAL(KIND=WP), INTENT(INOUT) :: X(LDX,*), Y(LDY,*)
      REAL(KIND=WP), INTENT(OUT)   :: Z(LDZ,*), B(LDB,*), & 
                                      W(LDW,*), S(LDS,*)
      REAL(KIND=WP), INTENT(OUT)   :: EIGS(*),  RES(*)
      REAL(KIND=WP), INTENT(OUT)   :: WORK(*)  
      INTEGER,       INTENT(OUT)   :: IWORK(*)
!............................................................
!     Purpose
!     =======
!     SSYDMD computes the Dynamic Mode Decomposition (DMD) for
!     a pair of data snapshot matrices. For the input matrices
!     X and Y such that Y = A*X with an unaccessible symmetric 
!     matrix A, SSYDMD computes a certain number of Ritz pairs 
!     of A using the standard Rayleigh-Ritz extraction from a 
!     subspace of range(X) that is determined using the leading  
!     left singular vectors of X. Optionally, SSYDMD returns  
!     the residuals of the computed Ritz pairs, the information 
!     needed for a refinement of the Ritz vectors, or the 
!     eigenvectors of the Exact DMD. 
!     For furter details see the references listed below.
!     For more details of the implementation see [3], [4].
!
!     References
!     ==========
!     [1] P. Schmid: Dynamic mode decomposition of numerical
!         and experimental data,
!         Journal of Fluid Mechanics 656, 5-28, 2010.
!     [2] Z. Drmac, I. Mezic, R. Mohr: Data driven modal
!         decompositions: analysis and enhancements,
!         SIAM J. on Sci. Comp. 40 (4), A2253-A2285, 2018.
!     [3] Z. Drmac: A LAPACK implementation of the Dynamic
!         Mode Decomposition I. Technical report. AIMDyn Inc.
!         October 2022, and LAPACK Working Note 298.
!     [4] Z. Drmac: A LAPACK implementation of the Dynamic
!         Mode Decomposition II. The symmetric/Hermitian DMD
!         (xSYDMD/xHEDMD) Technical report. AIMDyn Inc. 
!         November 2022.  LAPACK Working Note 300.    
!     [5] J. Tu, C. W. Rowley, D. M. Luchtenburg, S. L. 
!         Brunton, N. Kutz: On Dynamic Mode Decomposition:
!         Theory and Applications, Journal of Computational
!         Dynamics 1(2), 391 -421, 2014.
!     [6] P. J. Baddoo, B. Herrmann, B. J. McKeon, 
!         J. N.  Kutz, S. L. Brunton: Physics-informed 
!         dynamic mode decomposition (piDMD), arXiv:2112.04307.
!         
!......................................................................
!     Developed and supported by:
!     ===========================
!     Developed and coded by Zlatko Drmac, Faculty of Science,
!     University of Zagreb;  drmac@math.hr
!     In cooperation with
!     AIMdyn Inc., Santa Barbara, CA.
!     and supported by
!     - DARPA SBIR project "Koopman Operator-Based Forecasting
!     for Nonstationary Processes from Near-Term, Limited
!     Observational Data" Contract No: W31P4Q-21-C-0007
!     - DARPA PAI project "Physics-Informed Machine Learning
!     Methodologies" Contract No: HR0011-18-9-0033
!     - DARPA MoDyL project "A Data-Driven, Operator-Theoretic
!     Framework for Space-Time Analysis of Process Dynamics"
!     Contract No: HR0011-16-C-0116
!     Any opinions, findings and conclusions or recommendations 
!     expressed in this material are those of the author and 
!     do not necessarily reflect the views of the DARPA SBIR 
!     Program Office
!============================================================
!     Distribution Statement A:                       
!     Approved for Public Release, Distribution Unlimited.      
!     
!============================================================      
!............................................................
!     Arguments
!     =========
!     JOBS (input) CHARACTER*1
!     Determines whether the initial data snapshots are scaled
!     by a diagonal matrix.
!     'S' :: The data snapshots matrices X and Y are multiplied
!            with a diagonal matrix D so that X*D has unit
!            nonzero columns (in the Euclidean 2-norm)
!     'C' :: The snapshots are scaled as with the 'S' option.
!            If it is found that an i-th column of X is zero
!            vector and the corresponding i-th column of Y is
!            non-zero, then the i-th column of Y is set to
!            zero and a warning flag is raised.
!     'Y' :: The data snapshots matrices X and Y are multiplied
!            by a diagonal matrix D so that Y*D has unit
!            nonzero columns (in the Euclidean 2-norm)
!     'N' :: No data scaling.
!.....
!     JOBZ (input) CHARACTER*1
!     Determines whether the eigenvectors (Koopman modes) will
!     be computed.
!     'V' :: The eigenvectors (Koopman modes) will be computed
!            and returned in the matrix Z.
!            See the description of Z.
!     'F' :: The eigenvectors (Koopman modes) will be returned
!            in factored form as the product X(:,1:K)*W, where X
!            contains a POD basis (leading left singular vectors
!            of the data matrix X) and W contains the eigenvectors
!            of the corresponding Rayleigh quotient.
!            See the descriptions of K, X, W, Z.
!     'N' :: The eigenvectors are not computed.
!.....
!     JOBR (input) CHARACTER*1 
!     Determines whether to compute the residuals.
!     'R' :: The residuals for the computed eigenpairs will be
!            computed and stored in the array RES.
!            See the description of RES.
!            For this option to be legal, JOBZ must be 'V'.
!     'N' :: The residuals are not computed.
!.....
!     JOBF (input) CHARACTER*1
!     Specifies whether to store information needed for post-
!     processing (e.g. computing refined Ritz vectors)
!     'R' :: The matrix needed for the refinement of the Ritz
!            vectors is computed and stored in the array B.
!            See the description of B.
!     'E' :: The unscaled eigenvectors of the Exact DMD are 
!            computed and returned in the array B. See the
!            description of B.
!     'X' :: The Exact DMD vectors are orthogonalized and
!            returned in the array B. To preserve the 
!            residuals of the orthogonalized EDMD vectors
!            they are reordered and the reordering permutation
!            is stored and returned in the array IWORK.
!            See the descriptions of B and IWORK, and [4].
!     'N' :: No eigenvector refinement data is computed.
!.....
!     WHTSVD (input) INTEGER, WHSTVD in { 1, 2, 3, 4 }
!     Allows for a selection of the SVD algorithm from the
!     LAPACK library.
!     1 :: SGESVD (the QR SVD algorithm)
!     2 :: SGESDD (the Divide and Conquer algortihm; if enough
!          workspace available, this is the fastest option)
!     3 :: SGESVDQ (the preconditioned QR SVD  ; this and 4
!          are the most accurate options)
!     4 :: SGEJSV (the preconditioned Jacobi SVD; this and 3
!          are the most accurate options)
!     For the four methods above, a significant difference in
!     the accuracy of small singular values is possible if
!     the snapshots vary in norm so that X is severely
!     ill-conditioned. If small (smaller than EPS*||X||)
!     singular values are of interest and JOBS=='N',  then
!     the options (3, 4) give the most accurate results, where
!     the option 4 is slightly better and with stronger 
!     theoretical background.
!     If JOBS=='S', i.e. the columns of X will be normalized,
!     then all methods give nearly equally accurate results.
!.....
!     WHTSYM (input) INTEGER
!     Specifies the method for restoring the symmetry of the
!     Rayleigh quotient.
!     1 :: The lower triangle of the computed Rayleigh
!          quotient is used to symmetrize the matrix,
!     2 :: The formulas for the lower triangle of a 
!          truncated solution of the symmetric Procrustes
!          problem are used to symmetrize the computed
!          Rayleigh quotient.
!.....
!     WHTEIG (input) INTEGER
!     Specifies the symmetric eigensolver to compute the 
!     eigenvalues and eigenvectors of the symmetric Rayleigh 
!     quotient.
!     1 :: SSYEV  (the QR algorithm)
!     2 :: SSYEVD (the divide and conquer algorithm)
!.....
!     M (input) INTEGER, M>= 0
!     The state space dimension (the row dimension of X, Y).
!.....
!     N (input) INTEGER, 0 <= N <= M
!     The number of data snapshot pairs
!     (the number of columns of X and Y).
!.....
!     X (input/output) REAL(KIND=WP) M-by-N array
!     > On entry, X contains the data snapshot matrix X. It is 
!     assumed that the column norms of X are in the range of
!     the normalized floating point numbers. 
!     < On exit, the leading K columns of X contain a POD basis,
!     i.e. the leading K left singular vectors of the input
!     data matrix X, U(:,1:K). All N columns of X contain all
!     left singular vectors of the input matrix X.
!     See the descriptions of K, Z and W.  
!.....
!     LDX (input) INTEGER, LDX >= M
!     The leading dimension of the array X.
!.....
!     Y (input/workspace/output) REAL(KIND=WP) M-by-N array
!     > On entry, Y contains the data snapshot matrix Y
!     < On exit,
!     If JOBR == 'R', the leading K columns of Y  contain
!     the residual vectors for the computed Ritz pairs.
!     See the description of RES.
!     If JOBR == 'N', Y contains the original input data.
!.....
!     LDY (input) INTEGER , LDY >= M
!     The leading dimension of the array Y.
!.....
!     NRNK (input) INTEGER
!     Determines the mode how to compute the numerical rank,
!     i.e. how to truncate small singular values of the input
!     matrix X. On input, if
!     NRNK = -1 :: i-th singular value sigma(i) is truncated
!                  if sigma(i) <= TOL*sigma(1)
!     NRNK = -2 :: i-th singular value sigma(i) is truncated
!                  if sigma(i) <= TOL*sigma(i-1)
!     The numerical rank can be enforced by using positive 
!     value of NRNK as follows:      
!     0 < NRNK <= N :: at most NRNK largest singular values
!     will be used. If the number of the computed nonzero
!     singular values is less than NRNK, then only those
!     nonzero values will be used and the actually used
!     dimension is less than NRNK. The actual number of
!     the nonzero singular values is returned in the variable
!     K. See the descriptions of TOL and  K.
!.....
!     TOL (input) REAL(KIND=WP), 0 <= TOL < 1
!     The tolerance for truncating small singular values.
!     See the description of NRNK.
!.....
!     K (output) INTEGER,  0 <= K <= N
!     The dimension of the POD basis for the data snapshot
!     matrix X and the number of the computed Ritz pairs.
!     The value of K is determinet according to the rule set
!     by the parameters NRNK and TOL.
!     See the descriptions of NRNK and TOL.
!.....
!     EIGS (output) REAL(KIND=WP) N-by-1 array
!     The leading K (K<=N) entries of EIGS contain
!     the computed eigenvalues in ascending order.
!     If the eigenvectors are requested, then Z(:,i)
!     corresponds to EIGS(i). If JOBF == 'X', then
!     orthonormalised Exact DMD vectors are stored
!     in the array B and to the eigenvector B(:,i)
!     the corresponding eigenvalue is EIGS(IWORK(i)).
!     See the descriptions of K, Z, B and IWORK.
!.....      
!     Z (workspace/output) REAL(KIND=WP)  M-by-N array
!     If JOBZ =='V' then Z contains Ritz vectors.
!     If JOBZ == 'F', then the above descriptions hold for
!     the columns of X(:,1:K)*W(1:K,1:K), where the columns 
!     of W(1:k,1:K) are the computed eigenvectors of the 
!     K-by-K Rayleigh quotient. 
!     See the descriptions of EIGS, X and W.
!.....
!     LDZ (input) INTEGER , LDZ >= M
!     The leading dimension of the array Z.
!.....
!     RES (output) REAL(KIND=WP) N-by-1 array
!     RES(1:K) contains the residuals for the K computed 
!     Ritz pairs. 
!     RES(i) = || A * Z(:,i) - EIGS(i)*Z(:,i))||_2.    
!     If JOBF == 'X', the array IWORK on exit 
!     contains the permutation that sorts RES in
!     ascending order. 
!     See the description of JOBF, EIGS, Z and IWORK.      
!.....
!     B (output) REAL(KIND=WP)  M-by-N array.
!     IF JOBF =='R', B(1:M,1:K) contains A*U(:,1:K), and can
!     be used for computing the refined vectors; see further 
!     details in the provided references. 
!     If JOBF == 'E', B(1:M,1;K) contains 
!     A*U(:,1:K)*W(1:K,1:K), which are the vectors from the
!     Exact DMD, up to scaling by the inverse eigenvalues. 
!     Note that the EDMD vectors may not be even numerically
!     orthogonal and that the non-orthogonality may be 
!     substantial.
!     If JOBF == 'X', then the EDMD vectors 
!     A*U(:,1:K)*W(1:K,1:K) are orthonormalized. To preserve
!     information on the residuals, they are reordered and
!     the reordering permutation is stored in the array IWORK.
!     If JOBF =='N', then B is not referenced.
!     See the descriptions of JOBF, X, W, K, IWORK.      
!.....
!     LDB (input) INTEGER, LDB >= M
!     The leading dimension of the array B.
!.....
!     W (workspace/output) REAL(KIND=WP) N-by-N array
!     On exit, W(1:K,1:K) contains the K computed 
!     eigenvectors of the matrix Rayleigh quotient. 
!     The Ritz vectors (returned in Z) are the
!     product of X (containing a POD basis for the input
!     matrix X) and W. See the descriptions of K, S, X and Z.
!     W is also used as a workspace to temporarily store the
!     left singular vectors of X.
!.....
!     LDW (input) INTEGER, LDW >= N
!     The leading dimension of the array W.
!.....      
!     S (workspace/output) REAL(KIND=WP) N-by-N array
!     The array S(1:K,1:K) is used for the matrix Rayleigh
!     quotient. This content is overwritten during
!     the eigenvalue decomposition.
!     See the description of K.
!.....
!     LDS (input) INTEGER, LDS >= N
!     The leading dimension of the array S.
!.....
!     WORK (workspace/output) REAL(KIND=WP) LWORK-by-1 array
!     On exit, WORK(1:N) contains the singular values of 
!     X (for JOBS=='N') or column scaled X (JOBS=='S', 'C').
!     If WHTSVD==4, then WORK(N+1) and WORK(N+2) contain
!     scaling factor WORK(N+2)/WORK(N+1) used to scale X
!     and Y to avoid overflow in the SVD of X.
!     This may be of interest if the scaling option is off
!     and as many as possible smallest eigenvalues are 
!     desired to the highest feasible accuracy. 
!     If the call to SSYDMD is only workspace query, then
!     WORK(1) contains the minimal workspace length and
!     WORK(2) is the optimal workspace length. Hence, the
!     length of WORK is at least 2.
!     See the description of LWORK.
!.....
!     LWORK (input) INTEGER
!     The minimal length of the workspace vector WORK.
!     LWORK is calculated as follows:
!     If WHTSVD == 1 ::
!        If JOBZ == 'V', then
!        LWORK >= MAX(2, N + LWORK_SVD, N+LWORK_EIG),
!        where LWORK_EIG is the work length for the 
!        symmetric eigenvalue solver. If
!        -  WHTEIG = 1(SSYEV), LWORK_EIG = MAX(1,3*N-1)
!        -  WHTEIG = 2(SSYEVD), LWORK_EIG = 1+6*N+2*N**2  
!        If JOBZ == 'N'  then
!        LWORK >= MAX(2, N + LWORK_SVD, N+LWORK_EIG),
!        where 
!        - if WHTEIG = 1, LWORK_EIG = MAX(1,3*N-1)
!        - if WHTEIG = 2, LWORK_EIG = 2*N+1 
!        Here LWORK_SVD = MAX(1,3*N+M,5*N) is the minimal
!        workspace length of SGESVD     
!     If WHTSVD == 2 ::
!        If JOBZ == 'V', then
!        LWORK >= MAX(2, N + LWORK_SVD, N+LWORK_EIG)
!        LWORK_EIG is the work length for the 
!        symmetric eigenvalue solver. If 
!        -  WHTEIG = 1(SSYEV), LWORK_EIG = MAX(1,3*N-1)
!        -  WHTEIG = 2(SSYEVD), LWORK_EIG = 1+6*N+2*N**2      
!        If JOBZ == 'N', then
!        LWORK >= MAX(2, N + LWORK_SVD, N+LWORK_EIG),
!        where
!        - if WHTEIG = 1, LWORK_EIG = MAX(1,3*N-1)
!        - if WHTEIG = 2, LWORK_EIG = 2*N+1          
!        Here LWORK_SVD = MAX(M, 5*N*N+4*N)+3*N*N is the
!        minimal workspace length of SGESDD.  
!     If WHTSVD == 3 ::
!        If JOBZ == 'V', then
!        LWORK >= MAX(2, N+LWORK_SVD,N+LOWRK_EIG)
!        LWORK_EIG is the work length for the 
!        symmetric eigenvalue solver. If 
!        -  WHTEIG = 1(SSYEV), LWORK_EIG = MAX(1,3*N-1)
!        -  WHTEIG = 2(SSYEVD), LWORK_EIG = 1+6*N+2*N**2     
!        If JOBZ == 'N', then
!        LWORK >= MAX(2, N+LWORK_SVD,N+LWORK_EIG),
!        where
!        - if WHTEIG = 1, LWORK_EIG = MAX(1,3*N-1)
!        - if WHTEIG = 2, LWORK_EIG = 2*N+1          
!        Here LWORK_SVD = N+M+MAX(3*N+1,
!                        MAX(1,3*N+M,5*N),MAX(1,N))
!        is the minimal workspace length of SGESVDQ.
!     If WHTSVD == 4 ::
!        If JOBZ == 'V', then
!        LWORK >= MAX(2, N+LWORK_SVD,N+LWORK_EIG)
!        LWORK_EIG is the work length for the 
!        symmetric eigenvalue solver. If 
!        -  WHTEIG = 1(SSYEV), LWORK_EIG = MAX(1,3*N-1)
!        -  WHTEIG = 2(SSYEVD), LWORK_EIG = 1+6*N+2*N**2       
!        If JOBZ == 'N', then
!        LWORK >= MAX(2, N+LWORK_SVD,N+LWORK_EIG),
!        where
!        - if WHTEIG = 1, LWORK_EIG = MAX(1,3*N-1)
!        - if WHTEIG = 2, LWORK_EIG = 2*N+1       
!        Here LWORK_SVD = MAX(7,2*M+N,6*N+2*N*N) is the
!        minimal workspace length of SGEJSV.
!     Further, if JOBF=='X', then in addition
!     LWORK = MAX(LWORK,2*N+N) (for SGEQRF and SORGQR).
!     The above expressions are not simplified in order to
!     make the usage of WORK more transparent, and for 
!     easier checking. In any case, LWORK >= 2.     
!     If on entry LWORK = -1, then a workspace query is
!     assumed and the procedure only computes the minimal
!     and the optimal workspace lengths for both WORK and
!     IWORK. See the descriptions of WORK and IWORK.     
!.....          
!     IWORK (workspace/output) INTEGER LIWORK-by-1 array
!     Workspace that is required if WHTSVD equals
!     2 , 3 or 4. Further, if  JOBF=='X', it is used to return 
!     ordering of the orthonormalized Exact DMD eigenvectors,
!     so that EIGS(IWORK(i)) is the eigenvalue that corresponds to
!     the i-th EDMD vector. See the descriptions of JOBF and B.
!     If on entry LWORK =-1 or LIWORK=-1, then the
!     minimal length of IWORK is computed and returned in
!     IWORK(1). See the description of LIWORK.
!.....
!     LIWORK (input) INTEGER
!     The minimal length of the workspace vector IWORK.
!     LIWORK is determined in two steps. First:      
!     If WHTSVD == 1, then only IWORK(1) is used; LIWORK >=1
!     If WHTSVD == 2, then LIWORK >= MAX(1,8*MIN(M,N)).
!     If WHTSVD == 3, then LIWORK >= MAX(1,M+N-1).
!     If WHTSVD == 4, then LIWORK >= MAX(3,M+3*N).
!     If WHTEIG == 2 and JOBZ == 'V', then LIWORK >= MAX(1,3+5*N)   
!     Then, if JOBF == 'X', then LIWORK = MAX(LIWORK,N).
!     If on entry LIWORK = -1, then a workspace query is
!     assumed and the procedure only computes the minimal
!     and the optimal workspace lengths for both WORK and
!     IWORK. See the descriptions of WORK and IWORK.
!.....
!     INFO (output) INTEGER
!     -i < 0 :: On entry, the i-th argument had an
!               illegal value
!        = 0 :: Successful return.
!        = 1 :: Void input. Quick exit (M=0 or N=0).
!        = 2 :: The SVD computation of X did not converge.
!               Suggestion: Check the input data and/or
!               repeat with different WHTSVD.
!        = 3 :: The computation of the eigenvalues did not
!               converge.
!        = 4 :: If data scaling was requested on input and
!               the procedure found inconsistency in the data
!               such that for some column index i,
!               X(:,i) = 0 but Y(:,i) /= 0, then Y(:,i) is set
!               to zero if JOBS=='C'. The computation proceeds
!               with original or modified data and warning
!               flag is set with INFO=4.
!.............................................................
!.............................................................
!     Parameters
!     ~~~~~~~~~~
      REAL(KIND=WP), PARAMETER ::  ONE = 1.0_WP
      REAL(KIND=WP), PARAMETER :: ZERO = 0.0_WP
      
!     Local scalars
!     ~~~~~~~~~~~~~
      REAL(KIND=WP) :: OFL,    ROOTSC, SCALE,  SMALL,  & 
                       SSUM,   XSCL1,  XSCL2
      INTEGER       :: i,   j, IMINWR,  INFO1, INFO2,  &  
                       IWRSDD, LWRKEV, LWRSDD, LWRSVD, &
                       LWRSVQ, MLWORK, MWRKEV, MWRSDD, &
                       MWRSVD, MWRSVJ, MWRSVQ, NUMRNK, &
                       OLWORK
      LOGICAL       :: BADXY,  FORWRD, LQUERY, SCCOLX, &
                       SCCOLY, WNTEX,  WNTREF, WNTRES, &
                       WNTVEC
      CHARACTER     :: JOBZL,  T_OR_N
      CHARACTER     :: JSVOPT
      
!     Local arrays      
!     ~~~~~~~~~~~~      
      REAL(KIND=WP) :: RDUMMY(2), RDUMMY2(2)
!     External funcions (BLAS and LAPACK)
!     ~~~~~~~~~~~~~~~~~
      REAL(KIND=WP) SLANGE, SLAMCH, SNRM2
      EXTERNAL      SLANGE, SLAMCH, SNRM2, ISAMAX, ISAMIN
      INTEGER       ISAMAX, ISAMIN
      LOGICAL       SISNAN, LSAME
      EXTERNAL      SISNAN, LSAME 

!     External subroutines (BLAS and LAPACK)
!     ~~~~~~~~~~~~~~~~~~~~
      EXTERNAL      SAXPY,  SGEMM,   SSCAL 
      EXTERNAL      SSYEV,  SSYEVD,  SGEJSV, SGEQRF, SGESDD, &
                    SGESVD, SGESVDQ, SLACPY, SLAPMT, SLASCL, &
                    SLASSQ, SORGQR,  XERBLA
      
!     Intrinsic functions
!     ~~~~~~~~~~~~~~~~~~~
      INTRINSIC     INT, FLOAT, MAX, SQRT       
!............................................................  
!
!    Test the input arguments
!
      WNTRES = LSAME(JOBR,'R')
      SCCOLX = LSAME(JOBS,'S') .OR. LSAME(JOBS,'C')
      SCCOLY = LSAME(JOBS,'Y')
      WNTVEC = LSAME(JOBZ,'V')
      WNTREF = LSAME(JOBF,'R') 
      WNTEX  = LSAME(JOBF,'E') .OR. LSAME(JOBF,'X')
      INFO   = 0
      LQUERY = ( ( LWORK == -1 ) .OR. ( LIWORK == -1 ) )
!
      IF ( .NOT. (SCCOLX .OR. SCCOLY .OR. &
                                  LSAME(JOBS,'N')) )     THEN
          INFO = -1
      ELSE IF ( .NOT. (WNTVEC .OR. LSAME(JOBZ,'N')        & 
                              .OR. LSAME(JOBZ,'F')) )    THEN
          INFO = -2
      ELSE IF ( .NOT. (WNTRES .OR. LSAME(JOBR,'N')) .OR.  & 
                ( WNTRES .AND. (.NOT.WNTVEC) ) )         THEN
          INFO = -3
      ELSE IF ( .NOT. (WNTREF .OR. WNTEX .OR.             & 
                LSAME(JOBF,'N') ) )                      THEN
          INFO = -4
      ELSE IF ( .NOT.((WHTSVD == 1) .OR. (WHTSVD == 2) .OR.  & 
                      (WHTSVD == 3) .OR. (WHTSVD == 4) ))THEN
          INFO = -5
      ELSE IF ( .NOT.((WHTSYM == 1) .OR. (WHTSYM == 2))) THEN
          INFO = -6
      ELSE IF ( .NOT.((WHTEIG == 1) .OR. (WHTEIG == 2))) THEN
          INFO = -7
      ELSE IF ( M < 0 )   THEN
          INFO = -8
      ELSE IF ( ( N < 0 ) .OR. ( N > M ) ) THEN
          INFO = -9
      ELSE IF ( LDX < M ) THEN
          INFO = -11
      ELSE IF ( LDY < M ) THEN
          INFO = -13
      ELSE IF ( .NOT. (( NRNK == -2).OR.(NRNK == -1).OR. & 
                ((NRNK >= 1).AND.(NRNK <=N ))) )         THEN
          INFO = -14
      ELSE IF ( ( TOL < ZERO ) .OR. ( TOL >= ONE ) )     THEN
          INFO = -15
      ELSE IF ( LDZ < M ) THEN
          INFO = -19
      ELSE IF ( (WNTREF .OR. WNTEX ) .AND. ( LDB < M ) ) THEN
          INFO = -22
      ELSE IF ( LDW < N ) THEN
          INFO = -24
      ELSE IF ( LDS < N ) THEN
          INFO = -26
      END IF
!      
      IF ( INFO == 0 ) THEN  
          ! Compute the minimal and the optimal workspace 
          ! requirements. Simulate running the code and 
          ! determine minimal and optimal sizes of the
          ! workspace at any moment of the run.
         IF ( N == 0 ) THEN
             ! Quick return. All output except K is void. 
             ! INFO=1 signals the void input.
             ! In case of a workspace query, the default 
             ! minimal workspace lengths are returned.             
            IF ( LQUERY ) THEN  
                IWORK(1) = 1
                WORK(1)  = 2
                WORK(2)  = 2
            ELSE                
               K = 0
            END IF             
            INFO = 1  
            RETURN
         END IF
         MLWORK = MAX(2,N)
         OLWORK = MAX(2,N)
         IMINWR = 1
         SELECT CASE ( WHTSVD )
         CASE (1)
             ! The following is specified as the minimal 
             ! length of WORK in the definition of SGESVD:
             ! MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
             MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
             MLWORK = MAX(MLWORK,N + MWRSVD)
             IF ( LQUERY ) THEN 
                CALL SGESVD( 'O', 'S', M, N, X, LDX, WORK, & 
                           B, LDB, W, LDW, RDUMMY, -1, INFO1 )   
                LWRSVD = MAX( MWRSVD, INT( RDUMMY(1) ) )
                OLWORK = MAX(OLWORK,N + LWRSVD)           
             END IF  
         CASE (2)
             ! The following is specified as the minimal 
             ! length of WORK in the definition of SGESDD:
             ! MWRSDD = 3*MIN(M,N)*MIN(M,N) + 
             ! MAX( MAX(M,N),5*MIN(M,N)*MIN(M,N)+4*MIN(M,N) )
             ! IMINWR = 8*MIN(M,N)
             MWRSDD = 3*MIN(M,N)*MIN(M,N) +                &  
              MAX( MAX(M,N),5*MIN(M,N)*MIN(M,N)+4*MIN(M,N) )
             MLWORK = MAX(MLWORK,N + MWRSDD) 
             IMINWR = 8*MIN(M,N) 
             IF ( LQUERY ) THEN
                CALL SGESDD( 'O', M, N, X, LDX, WORK, B,     &
                     LDB, W, LDW, RDUMMY, -1, IWORK, INFO1 ) 
                LWRSDD = MAX( MWRSDD, INT( RDUMMY(1) ) )
                OLWORK = MAX(OLWORK,N + LWRSDD)
             END IF                     
         CASE (3)
             !LWQP3 = 3*N+1
             !LWORQ = MAX(N, 1)
             !MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
             !MWRSVQ = N + MAX( LWQP3, MWRSVD, LWORQ ) + MAX(M,2)
             !MLWORK = N +  MWRSVQ 
             !IMINWR = M+N-1            
             CALL SGESVDQ( 'H', 'P', 'N', 'R', 'R', M, N, &
                             X, LDX, WORK, Z, LDZ, W, LDW,   & 
                             NUMRNK, IWORK, LIWORK, RDUMMY,  &
                             -1, RDUMMY2, -1, INFO1 )
             IMINWR = IWORK(1)
             MWRSVQ = INT(RDUMMY(2)) + INT(RDUMMY2(1))
             MLWORK = MAX(MLWORK,N+MWRSVQ) 
             IF ( LQUERY ) THEN 
                LWRSVQ = MAX( MWRSVQ, INT(RDUMMY(1)) ) 
                OLWORK = MAX(OLWORK,N+LWRSVQ+INT(RDUMMY2(1)))
             END IF  
         CASE (4) 
             JSVOPT = 'J'
             !MWRSVJ = MAX( 7, 2*M+N, 6*N+2*N*N ) ! for JSVOPT='V'
             MWRSVJ = MAX( 7, 2*M+N, 4*N+N*N, 2*N+N*N+6 )
             MLWORK = MAX(MLWORK,N+MWRSVJ) 
             IMINWR = MAX( 3, M+3*N )
             IF ( LQUERY ) THEN 
                OLWORK =  MAX(OLWORK,N+MWRSVJ) 
             END IF
         END SELECT
         IF ( WNTVEC .OR. WNTEX .OR. LSAME(JOBZ,'F') ) THEN 
             JOBZL = 'V'
         ELSE
             JOBZL = 'N'
         END IF
         SELECT CASE ( WHTEIG )
         CASE (1)
         ! Workspace calculation to the SSYEV call
         MWRKEV = MAX( 1, 3*N-1 ) 
         MLWORK = MAX(MLWORK,N+MWRKEV) 
         IF ( LQUERY ) THEN 
            CALL SSYEV( JOBZL, 'U', N, S, LDS, EIGS, RDUMMY, &
                  -1, INFO1 )   ! LAPACK CALL  
                LWRKEV = MAX( MWRKEV, INT(RDUMMY(1)) )
                OLWORK = MAX( OLWORK, N+LWRKEV )
         END IF 
         CASE (2)
          IF ( LSAME(JOBZL,'V') ) THEN 
             MWRKEV = 1 + 6*N + 2*N*N
             IWRSDD = 3+5*N
         ELSE
             MWRKEV = MAX( 1, 2*N+1)
             IWRSDD = 1
         END IF 
         MLWORK = MAX(MLWORK,N+MWRKEV)         
         IF ( LQUERY ) THEN 
            CALL SSYEVD( JOBZL, 'U', N, S, LDS, EIGS, RDUMMY, &
                  -1, IWORK, -1, INFO1 )   ! LAPACK CALL  
                LWRKEV = MAX( MWRKEV, INT(RDUMMY(1)) )
                OLWORK = MAX( OLWORK, N+LWRKEV )
                IWRSDD = IWORK(1)
         END IF 
         IMINWR = MAX(IMINWR,IWRSDD)  
         END SELECT
         IF ( LSAME(JOBF,'X') ) THEN 
         MLWORK = MAX(MLWORK,2*N+N) 
         ! SGEQRF and SORGQR need >= 2*N locations
         IF ( LQUERY ) THEN
             CALL SGEQRF( M, N,   B, LDB, RDUMMY, RDUMMY, &
                         -1, INFO1 )
            OLWORK = MAX( OLWORK, 2*N+INT(RDUMMY(1)) )
            CALL SORGQR( M, N, N, B, LDB, RDUMMY, RDUMMY, &
                        -1, INFO1 ) 
            OLWORK = MAX( OLWORK, 2*N+INT(RDUMMY(1)) )
         END IF     
          IMINWR = MAX( IMINWR, N )
          END IF 
         IF ( LIWORK < IMINWR .AND. (.NOT.LQUERY) ) INFO = -30
         IF (  LWORK < MLWORK .AND. (.NOT.LQUERY) ) INFO = -28   
      END IF
!      
      IF( INFO /= 0 ) THEN
         CALL XERBLA( 'SSYDMD', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
!     Return minimal and optimal workspace sizes
          IWORK(1) = IMINWR
          WORK(1)  = MLWORK
          WORK(2)  = OLWORK
          RETURN
      END IF   
!............................................................
!      
      OFL   = SLAMCH('O')
      SMALL = SLAMCH('S')
      BADXY = .FALSE.
!
!     <1> Optional scaling of the snapshots (columns of X, Y)
!     ==========================================================
      IF ( SCCOLX ) THEN 
          ! The columns of X will be normalized. 
          ! To prevent overflows, the column norms of X are 
          ! carefully computed using SLASSQ.    
          K = 0 
          DO i = 1, N
            !WORK(i) = DNRM2( M, X(1,i), 1 )  
            SCALE  = ZERO
            CALL SLASSQ( M, X(1,i), 1, SCALE, SSUM )
            IF ( SISNAN(SCALE) .OR. SISNAN(SSUM) ) THEN
                K    =  0 
                INFO = -10 
                CALL XERBLA('SSYDMD',-INFO)
            END IF 
            IF ( (SCALE /= ZERO) .AND. (SSUM /= ZERO) ) THEN 
               ROOTSC = SQRT(SSUM)
               IF ( SCALE .GE. (OFL / ROOTSC) ) THEN 
!                 Norm of X(:,i) overflows. First, X(:,i) 
!                 is scaled by
!                 ( ONE / ROOTSC ) / SCALE = 1/||X(:,i)||_2. 
!                 Next, the norm of X(:,i) is stored without
!                 overflow as WORK(i) = - SCALE * (ROOTSC/M), 
!                 the minus sign indicating the 1/M factor.  
!                 Scaling is performed without overflow, and
!                 underflow may occur in the smallest entries
!                 of X(:,i). The relative backward and forward
!                 errors are small in the ell_2 norm.                      
                  CALL SLASCL( 'G', 0, 0, SCALE, ONE/ROOTSC, &
                               M, 1, X(1,i), M, INFO1 )                 
                  WORK(i) = - SCALE * ( ROOTSC / FLOAT(M) )
               ELSE 
!                 X(:,i) will be scaled to unit 2-norm
                  WORK(i) =   SCALE * ROOTSC 
                  CALL SLASCL( 'G',0, 0, WORK(i), ONE, M, 1, &
                               X(1,i), M, INFO1 )              ! LAPACK CALL
!                 X(1:M,i) = (ONE/WORK(i)) * X(1:M,i)          ! INTRINSIC
               END IF 
            ELSE
               WORK(i) = ZERO
               K = K + 1 
            END IF          
          END DO
          IF ( K == N ) THEN
          ! All columns of X are zero. Return error code -8.
          ! (the 8th input variable had an illegal value)
          K = 0 
          INFO = -8
          CALL XERBLA('SSYDMD',-INFO)
          RETURN       
          END IF
          DO i = 1, N
!           Now, apply the same scaling to the columns of Y.        
            IF ( WORK(i) >  ZERO ) THEN 
                CALL SSCAL( M, ONE/WORK(i), Y(1,i), 1 )  ! BLAS CALL
!               Y(1:M,i) = (ONE/WORK(i)) * Y(1:M,i)      ! INTRINSIC                 
            ELSE IF ( WORK(i) < ZERO ) THEN 
                CALL SLASCL( 'G', 0, 0, -WORK(i),          & 
                     ONE/FLOAT(M), M, 1, Y(1,i), M, INFO1 ) ! LAPACK CALL
            ELSE IF ( Y(ISAMAX(M, Y(1,i),1),i )  & 
                                            /= ZERO ) THEN 
!               X(:,i) is zero vector. For consistency, 
!               Y(:,i) should also be zero. If Y(:,i) is not
!               zero, then the data might be inconsistent or
!               corrupted. If JOBS == 'C', Y(:,i) is set to 
!               zero and a warning flag is raised. 
!               The computation continues but the
!               situation will be reported in the output.
                BADXY = .TRUE.
                IF ( LSAME(JOBS,'C')) & 
                CALL SSCAL( M, ZERO, Y(1,i), 1 )  ! BLAS CALL                  
            END IF         
          END DO
      END IF  
  ! 
      IF ( SCCOLY ) THEN 
          ! The columns of Y will be normalized. 
          ! To prevent overflows, the column norms of Y are 
          ! carefully computed using SLASSQ.         
          DO i = 1, N
            !WORK(i) = DNRM2( M, Y(1,i), 1 )  
            SCALE  = ZERO
            CALL SLASSQ( M, Y(1,i), 1, SCALE, SSUM )
            IF ( SISNAN(SCALE) .OR. SISNAN(SSUM) ) THEN
                K    =  0 
                INFO = -12
                CALL XERBLA('SSYDMD',-INFO)
            END IF
            IF ( SCALE /= ZERO  .AND. (SSUM /= ZERO) ) THEN 
               ROOTSC = SQRT(SSUM)
               IF ( SCALE .GE. (OFL / ROOTSC) ) THEN 
!                 Norm of Y(:,i) overflows. First, Y(:,i) 
!                 is scaled by
!                 ( ONE / ROOTSC ) / SCALE = 1/||Y(:,i)||_2. 
!                 Next, the norm of Y(:,i) is stored without
!                 overflow as WORK(i) = - SCALE * (ROOTSC/M), 
!                 the minus sign indicating the 1/M factor.  
!                 Scaling is performed without overflow, and
!                 underflow may occur in the smallest entries
!                 of Y(:,i). The relative backward and forward
!                 errors are small in the ell_2 norm.                      
                  CALL SLASCL( 'G', 0, 0, SCALE, ONE/ROOTSC, &
                               M, 1, Y(1,i), M, INFO1 )                 
                  WORK(i) = - SCALE * ( ROOTSC / FLOAT(M) )
               ELSE 
!                 X(:,i) will be scaled to unit 2-norm
                  WORK(i) =   SCALE * ROOTSC 
                  CALL SLASCL( 'G',0, 0, WORK(i), ONE, M, 1, &
                               Y(1,i), M, INFO1 )              ! LAPACK CALL
!                 Y(1:M,i) = (ONE/WORK(i)) * Y(1:M,i)          ! INTRINSIC
               END IF 
            ELSE
               WORK(i) = ZERO
            END IF          
         END DO
         DO i = 1, N
!           Now, apply the same scaling to the columns of X.              
            IF ( WORK(i) >  ZERO ) THEN 
                CALL SSCAL( M, ONE/WORK(i), X(1,i), 1 )  ! BLAS CALL
!               X(1:M,i) = (ONE/WORK(i)) * X(1:M,i)      ! INTRINSIC                 
            ELSE IF ( WORK(i) < ZERO ) THEN 
                CALL SLASCL( 'G', 0, 0, -WORK(i),          & 
                     ONE/FLOAT(M), M, 1, X(1,i), M, INFO1 ) ! LAPACK CALL
            ELSE IF ( X(ISAMAX(M, X(1,i),1),i )  & 
                                           /= ZERO ) THEN 
!               Y(:,i) is zero vector.  If X(:,i) is not
!               zero, then a warning flag is raised. 
!               The computation continues but the
!               situation will be reported in the output.
                BADXY = .TRUE.
            END IF         
         END DO
       END IF  
!      
!     <2> SVD of the data snapshot matrix X.
!     ===================================== 
!     The left singular vectors are stored in the array X.
!     The right singular vectors are in the array W. 
!     The array W will later on contain the eigenvectors 
!     of a Rayleigh quotient.
      NUMRNK = N 
      SELECT CASE ( WHTSVD )
         CASE (1) 
             CALL SGESVD( 'O', 'S', M, N, X, LDX, WORK, B,     & 
                  LDB, W, LDW, WORK(N+1), LWORK-N, INFO1 ) ! LAPACK CALL
             T_OR_N = 'T'
         CASE (2) 
            CALL SGESDD( 'O', M, N, X, LDX, WORK, B, LDB, W,   &
                 LDW, WORK(N+1), LWORK-N, IWORK, INFO1 )   ! LAPACK CALL
            T_OR_N = 'T'
         CASE (3)  
              CALL SGESVDQ( 'H', 'P', 'N', 'R', 'R', M, N,     &
                   X, LDX, WORK, Z, LDZ, W, LDW,               & 
                   NUMRNK, IWORK, LIWORK, WORK(N+MAX(2,M)+1),  & 
                   LWORK-N-MAX(2,M), WORK(N+1), MAX(2,M), INFO1)     ! LAPACK CALL  
              CALL SLACPY( 'A', M, NUMRNK, Z, LDZ, X, LDX )   ! LAPACK CALL 
         T_OR_N = 'T'
         CASE (4) 
              CALL SGEJSV( 'F', 'U', JSVOPT, 'N', 'N', 'P', M, & 
                   N, X, LDX, WORK, Z, LDZ, W, LDW, &
                   WORK(N+1), LWORK-N, IWORK, INFO1 )    ! LAPACK CALL 
              CALL SLACPY( 'A', M, N, Z, LDZ, X, LDX )   ! LAPACK CALL
              T_OR_N = 'N'
              XSCL1 = WORK(N+1)
              XSCL2 = WORK(N+2)
              IF ( XSCL1 /=  XSCL2 ) THEN 
                 ! This is an exceptional situation. If the
                 ! data matrices are not scaled and the 
                 ! largest singular value of X overflows. 
                 ! In that case SGEJSV can return the SVD
                 ! in scaled form. The scaling factor can be used 
                 ! to rescale the data (X and Y).              
                 CALL SLASCL( 'G', 0, 0, XSCL1, XSCL2, M, N, Y, LDY, INFO2  )    
              END IF       
      END SELECT  
!         
      IF ( INFO1 > 0 ) THEN 
         ! The SVD selected subroutine did not converge. 
         ! Return with an error code.  
         INFO = 2 
         RETURN
      END IF       
! 
      IF ( WORK(1) == ZERO ) THEN
          ! The largest computed singular value of (scaled)
          ! X is zero. Return error code -8 
          ! (the 8th input variable had an illegal value).
          K = 0 
          INFO = -8 
          CALL XERBLA('SSYDMD',-INFO)
          RETURN       
      END IF
!      
      !<3> Determine the numerical rank of the data 
      !    snapshots matrix X. This depends on the 
      !    parameters NRNK and TOL.
                
      SELECT CASE ( NRNK )
          CASE ( -1 )
               K = 1 
               DO i = 2, NUMRNK 
                 IF ( ( WORK(i) <= WORK(1)*TOL ) .OR. &
                      ( WORK(i) <= SMALL ) ) EXIT  
                 K = K + 1     
               END DO
          CASE ( -2 )
               K = 1 
               DO i = 1, NUMRNK-1 
                 IF ( ( WORK(i+1) <= WORK(i)*TOL  ) .OR. & 
                      ( WORK(i) <= SMALL ) ) EXIT  
                 K = K + 1     
               END DO
          CASE DEFAULT
               K = 1
               DO i = 2, NRNK
                  IF ( WORK(i) <= SMALL ) EXIT
                  K = K + 1
               END DO
          END SELECT       
      !   Now, U = X(1:M,1:K) is the SVD/POD basis for the  
      !   snapshot data in the input matrix X.
      !<4> Compute the Rayleigh quotient S = U^T * A * U.
      !    Depending on the requsted outputs, the computation
      !    is organized to compute additional auxiliary 
      !    matrices (for the residuals and refinements).
      !    
      !    In all formulas below, we need V_k*Sigma_k^(-1)
      !    where either V_k is in W(1:N,1:K), or V_k^T is in
      !    W(1:K,1:N). Here Sigma_k=diag(WORK(1:K)).
      IF ( LSAME(T_OR_N, 'N') ) THEN
          DO i = 1, K 
           CALL SSCAL( N, ONE/WORK(i), W(1,i), 1 )    ! BLAS CALL   
           ! W(1:N,i) = (ONE/WORK(i)) * W(1:N,i)      ! INTRINSIC
          END DO
      ELSE
          ! This non-unit stride access is due to the fact
          ! that SGESVD, SGESVDQ and SGESDD return the 
          ! transposed matrix of the right singular vectors.
          !DO i = 1, K 
          ! CALL SSCAL( N, ONE/WORK(i), W(i,1), LDW )  ! BLAS CALL   
          ! ! W(i,1:N) = (ONE/WORK(i)) * W(i,1:N)      ! INTRINSIC
          !END DO
          DO i = 1, K
              WORK(N+i) = ONE/WORK(i)
          END DO      
          DO j = 1, N 
             DO i = 1, K 
                 W(i,j) = (WORK(N+i))*W(i,j)
             END DO
          END DO       
      END IF
!      
      IF ( WNTREF ) THEN      
         !
         ! Need A*U(:,1:K)=Y*V_k*inv(diag(WORK(1:K))) 
         ! for computing the refined Ritz vectors 
         ! (optionally, outside SSYDMD).
          CALL SGEMM( 'N', T_OR_N, M, K, N, ONE, Y, LDY, W, & 
                      LDW, ZERO, Z, LDZ )                        ! BLAS CALL  
          ! Z(1:M,1:K)=MATMUL(Y(1:M,1:N),TRANSPOSE(W(1:K,1:N)))  ! INTRINSIC, for T_OR_N=='T'
          ! Z(1:M,1:K)=MATMUL(Y(1:M,1:N),W(1:N,1:K))             ! INTRINSIC, for T_OR_N=='N'                      
          !          
          ! At this point Z contains
          ! A * U(:,1:K) = Y * V_k * Sigma_k^(-1), and 
          ! this is needed for computing the residuals. 
          ! This matrix is  returned in the array B and
          ! it can be used to compute refined Ritz vectors.
          CALL SLACPY( 'A', M, K, Z, LDZ, B, LDB )   ! BLAS CALL
          ! B(1:M,1:K) = Z(1:M,1:K)                  ! INTRINSIC 
          
          CALL SGEMM( 'T', 'N', K, K, M, ONE, X, LDX, Z, & 
                      LDZ, ZERO, S, LDS )                        ! BLAS CALL
          ! S(1:K,1:K) = MATMUL(TANSPOSE(X(1:M,1:K)),Z(1:M,1:K)) ! INTRINSIC
          ! At this point S = U^T * A * U is the Rayleigh quotient.     
      ELSE         
        ! A * U(:,1:K) is not explicitly needed and the 
        ! computation is organized differently. The Rayleigh
        ! quotient is computed more efficiently.   
        CALL SGEMM( 'T', 'N', K, N, M, ONE, X, LDX, Y, LDY, & 
                   ZERO, Z, LDZ )                                   ! BLAS CALL
        ! Z(1:K,1:N) = MATMUL( TRANSPOSE(X(1:M,1:K)), Y(1:M,1:N) )  ! INTRINSIC
        ! In the two SGEMM calls here, can use K for LDZ.
        CALL SGEMM( 'N', T_OR_N, K, K, N, ONE, Z, LDZ, W, & 
                    LDW, ZERO, S, LDS )                         ! BLAS CALL       
        ! S(1:K,1:K) = MATMUL(Z(1:K,1:N),TRANSPOSE(W(1:K,1:N))) ! INTRINSIC, for T_OR_N=='T'
        ! S(1:K,1:K) = MATMUL(Z(1:K,1:N),(W(1:N,1:K)))          ! INTRINSIC, for T_OR_N=='N'
        ! At this point S = U^T * A * U is the Rayleigh quotient.
        ! If the residuals are requested, save scaled V_k into Z. 
        ! Recal that V_k or V_k^T is stored in W.
        IF ( WNTRES .OR. WNTEX ) THEN
          IF ( LSAME(T_OR_N, 'N') ) THEN 
              CALL SLACPY( 'A', N, K, W, LDW, Z, LDZ )
          ELSE
              CALL SLACPY( 'A', K, N, W, LDW, Z, LDZ )           
          END IF       
        END IF
      END IF
      
      SELECT CASE ( WHTSYM )
      CASE (1)    
           CALL SLACPY( 'L', K, K, S, LDS, W, LDW ) 
      CASE (2)
          ! This is the symmetrizer from the piDMD [6],
          ! based on a solution of the symmetric Procrustes
          ! problem. Here included  for comparisons/study and
          ! for the sake of completeness.
          DO i = 1, K-1 
              W(i,i) = S(i,i)
              DO j = i+1, K
                 W(j,i) = ( WORK(i)*(S(j,i)*WORK(i))   +  &
                            WORK(j)*(S(i,j)*WORK(j)) ) /  &
                           ( WORK(i)**2 + WORK(j)**2 ) 
              END DO
          END DO
          W(k,k) = S(k,k)
      END SELECT
      !
      !<5> Compute the Ritz values and (if requested) the 
      !   right eigenvectors of the Rayleigh quotient.
      !
      !    The LAPACK eigensolvers SSYEV and SSYEVD return the
      !    eigenvectors in the array that contains upper or 
      !    lower triangle of the symmetric Rayleigh quotient. 
      ! 
      SELECT CASE ( WHTEIG )
      CASE (1) 
          CALL SSYEV( JOBZL, 'L', K, W, LDW, EIGS, WORK(N+1),   &
                  LWORK-N, INFO1 )   ! LAPACK CALL            
      CASE (2)
           CALL SSYEVD( JOBZL, 'L', K, W, LDW, EIGS, WORK(N+1), &
                  LWORK-N, IWORK, LIWORK, INFO1 )   ! LAPACK CALL       
      END SELECT
                  
      !
      ! W(1:K,1:K) contains the eigenvectors of the Rayleigh 
      ! quotient. 
      IF ( INFO1 > 0 ) THEN
         ! SSYEV/SSYEVD failed to compute the eigenvalues and 
         ! eigenvectors of the Rayleigh quotient. 
         INFO = 3 
         RETURN
      END IF
!      
      ! <6> Compute the eigenvectors (if requested) and, 
      ! the residuals (if requested).
      !
      IF ( WNTVEC .OR. WNTEX ) THEN        
      IF ( WNTRES ) THEN
          IF ( WNTREF ) THEN 
            ! Here, if the refinement is requested, we have 
            ! A*U(:,1:K) already computed and stored in Z. 
            ! For the residuals, need Y = A * U(:,1;K) * W.
            ! W is stored in S.  ? copy in Q 
            CALL SGEMM( 'N', 'N', M, K, K, ONE, Z, LDZ, W, & 
                       LDW, ZERO, Y, LDY )               ! BLAS CALL 
            ! Y(1:M,1:K) = Z(1:M,1:K) * W(1:K,1:K)       ! INTRINSIC
            ! This frees Z; Y contains A * U(:,1:K) * W.
          ELSE
            ! Compute S = V_k * Sigma_k^(-1) * W, where  
            ! V_k * Sigma_k^(-1) is stored in Z 
            CALL SGEMM( T_OR_N, 'N', N, K, K, ONE, Z, LDZ, &
                       W, LDW, ZERO, S, LDS) 
            ! Then, compute Z = Y * S = 
            ! = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
            ! = A * U(:,1:K) * W(1:K,1:K)
            CALL SGEMM( 'N', 'N', M, K, N, ONE, Y, LDY, S, &
                       LDS, ZERO, Z, LDZ)
            ! Save a copy of Z into Y and free Z for holding
            ! the Ritz vectors.
            CALL SLACPY( 'A', M, K, Z, LDZ, Y, LDY )
            IF ( WNTEX ) CALL SLACPY( 'A', M, K, Z, LDZ, B, LDB )
          END IF   
      ELSE IF ( WNTEX ) THEN
          ! Compute S = V_k * Sigma_k^(-1) * W, where  
            ! V_k * Sigma_k^(-1) is stored in Z 
            CALL SGEMM( T_OR_N, 'N', N, K, K, ONE, Z, LDZ, &
                       W, LDW, ZERO, S, LDS) 
            ! Then, compute Z = Y * S = 
            ! = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
            ! = A * U(:,1:K) * W(1:K,1:K)
            CALL SGEMM( 'N', 'N', M, K, N, ONE, Y, LDY, S, &
                       LDS, ZERO, B, LDB )
            ! The above call replaces the following two calls
            ! that were used in the developing-testing phase.
            ! CALL SGEMM( 'N', 'N', M, K, N, ONE, Y, LDY, S, &
            !           LDS, ZERO, Z, LDZ)
            ! Save a copy of Z into B and free Z for holding
            ! the Ritz vectors.
            ! CALL SLACPY( 'A', M, K, Z, LDZ, B, LDB )
      END IF
!      
      ! Compute the Ritz vectors
      IF ( WNTVEC ) CALL SGEMM( 'N', 'N', M, K, K, ONE, X, LDX, W, LDW, & 
                   ZERO, Z, LDZ )                           ! BLAS CALL
      ! Z(1:M,1:K) = MATMUL(X(1:M,1:K), W(1:K,1:K))         ! INTRINSIC 
!      
      IF ( WNTRES ) THEN
         DO i = 1, K 
                ! Compute the residuals
                CALL SAXPY( M, -EIGS(i), Z(1,i), 1, Y(1,i), 1 )       ! BLAS CALL
                ! Y(1:M,i) = Y(1:M,i) - EIGS(i) * Z(1:M,i)            ! INTRINSIC
                RES(i) = SNRM2( M, Y(1,i), 1)                         ! BLAS CALL
         END DO
      END IF
      END IF 
!      
      IF ( LSAME(JOBF,'X') ) THEN 
          ! If the Exact DMD eigenvectors are requested, the
          ! original EDMD vectors must be orthogonalized. 
          ! Orthogonalization may change the vector so that the
          ! corresponding residuals may increase. (Data driven
          ! setting does not allow recomputing the Razleigh 
          ! quotients.) To preseve the quality of the best EDMD
          ! vectors, orthogonalization is prformed in order of
          ! increasing residuals. For more details see [4].
          DO i = 1, K
              IWORK(i) = i
          END DO
          CALL SCOPY( K, RES, 1, WORK(N+1), 1 )
          DO i = 1, K-1
             j = ISAMIN( K-i+1, WORK(N+i), 1 ) + i - 1
             IF ( j /= i ) THEN
                 INFO1    = IWORK(i)
                 IWORK(i) = IWORK(j)
                 IWORK(j) = INFO1
                 SCALE     = WORK(N+i)
                 WORK(N+i) = WORK(N+j)
                 WORK(N+j) = SCALE
             END IF     
          END DO
          FORWRD = .TRUE.
          CALL SLAPMT( FORWRD, M, K, B, LDB, IWORK )
          ! Here we need the Gram-Schmidt orthogonalization 
          ! of the columns of B. The following two lines 
          ! use the QR factorization subroutine SGEQRF. This
          ! can be replaced with a more efficient Gram-Schmidt
          ! implementation. The matrix B is not expected to
          ! be ill-conditioned, so Gram-Schmid will be OK.
          CALL SGEQRF( M, K, B, LDB, WORK(N+1), WORK(N+K+1), &
                      LWORK-(N+K), INFO1 )
          CALL SORGQR( M, K, K, B, LDB, WORK(N+1),           &
                      WORK(N+K+1), LWORK-(N+K), INFO1 )   
      END IF
      
      IF ( WHTSVD == 4 ) THEN 
          WORK(N+1) = XSCL1
          WORK(N+2) = XSCL2
      END IF      
!      
!     Successful exit.     
      IF ( .NOT. BADXY ) THEN 
         INFO = 0
      ELSE
         ! A warning on possible data inconsistency. 
         ! This shouild be a rare event. 
         INFO = 4
      END IF
!............................................................      
      RETURN
!     ......
      END SUBROUTINE SSYDMD 
    