* 
*  Definition:
*  ===========
*
*       SUBROUTINE SGEQR( M, N, A, LDA, WORK1, LWORK1, WORK2, LWORK2,
*                        INFO)
* 
*       .. Scalar Arguments ..
*       INTEGER           INFO, LDA, M, N, LWORK1, LWORK2
*       ..
*       .. Array Arguments ..
*       REAL              A( LDA, * ), WORK1( * ), WORK2( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> 
*> SGEQR computes a QR factorization of an M-by-N matrix A, 
*> using SLATSQR when A is tall and skinny 
*> (M sufficiently greater than N), and otherwise SGEQRT:          
*> A = Q * R .
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
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix R 
*>          (R is upper triangular if M >= N);
*>          the elements below the diagonal represent Q (see Further Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK1
*> \verbatim
*>          WORK1 is REAL array, dimension (MAX(1,LWORK1))
*>          WORK1 contains part of the data structure used to store Q.
*>          WORK1(1): algorithm type = 1, to indicate output from 
*>                    DLATSQR or DGEQRT
*>          WORK1(2): optimum size of WORK1
*>          WORK1(3): minimum size of WORK1
*>          WORK1(4): row block size
*>          WORK1(5): column block size
*>          WORK1(6:LWORK1): data structure needed for Q, computed by 
*>                           SLATSQR or SGEQRT
*> \endverbatim
*>
*> \param[in] LWORK1
*> \verbatim
*>          LWORK1 is INTEGER
*>          The dimension of the array WORK1.
*>          If LWORK1 = -1, then a query is assumed. In this case the 
*>          routine calculates the optimal size of WORK1 and
*>          returns this value in WORK1(2), and calculates the minimum 
*>          size of WORK1 and returns this value in WORK1(3). 
*>          No error message related to LWORK1 is issued by XERBLA when 
*>          LWORK1 = -1.
*> \endverbatim
*>
*> \param[out] WORK2
*> \verbatim
*>         (workspace) REAL array, dimension (MAX(1,LWORK2))      
*> \endverbatim
*>
*> \param[in] LWORK2
*> \verbatim
*>          LWORK2 is INTEGER
*>          The dimension of the array WORK2.
*>          If LWORK2 = -1, then a query is assumed. In this case the 
*>          routine calculates the optimal size of WORK2 and 
*>          returns this value in WORK2(1), and calculates the minimum
*>          size of WORK2 and returns this value in WORK2(2).
*>          No error message related to LWORK2 is issued by XERBLA when
*>          LWORK2 = -1.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>  Depending on the matrix dimensions M and N, and row and column
*>  block sizes MB and NB returned by ILAENV, GEQR will use either
*>  LATSQR (if the matrix is tall-and-skinny) or GEQRT to compute
*>  the QR decomposition. 
*>  The output of LATSQR or GEQRT representing Q is stored in A and in
*>  array WORK1(6:LWORK1) for later use. 
*>  WORK1(2:5) contains the matrix dimensions M,N and block sizes MB,NB 
*>  which are needed to interpret A and WORK1(6:LWORK1) for later use. 
*>  WORK1(1)=1 indicates that the code needed to take WORK1(2:5) and 
*>  decide whether LATSQR or GEQRT was used is the same as used below in 
*>  GEQR. For a detailed description of A and WORK1(6:LWORK1), see 
*>  Further Details in LATSQR or GEQRT.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE SGEQR( M, N, A, LDA, WORK1, LWORK1, WORK2, LWORK2, 
     $   INFO)
*
*  -- LAPACK computational routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
*     November 2013
*
*     .. Scalar Arguments ..
      INTEGER           INFO, LDA, M, N, LWORK1, LWORK2
*     ..
*     .. Array Arguments ..
      REAL              A( LDA, * ), WORK1( * ), WORK2( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL    LQUERY, LMINWS
      INTEGER    MB, NB, I, II, KK, MINLW1, NBLCKS
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           SLATSQR, SGEQRT, XERBLA
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
*
      LQUERY = ( LWORK1.EQ.-1 .OR. LWORK2.EQ.-1 )
*
*     Determine the block size 
*    
      IF ( MIN(M,N).GT.0 ) THEN
        MB = ILAENV( 1, 'SGEQR ', ' ', M, N, 1, -1)
        NB = ILAENV( 1, 'SGEQR ', ' ', M, N, 2, -1)
      ELSE
        MB = M
        NB = 1
      END IF
      IF( MB.GT.M.OR.MB.LE.N) MB = M
      IF( NB.GT.MIN(M,N).OR.NB.LT.1) NB = 1
      MINLW1 = N + 5
      IF ((MB.GT.N).AND.(M.GT.N)) THEN
        IF(MOD(M-N, MB-N).EQ.0) THEN
          NBLCKS = (M-N)/(MB-N)
        ELSE
          NBLCKS = (M-N)/(MB-N) + 1
        END IF
      ELSE
        NBLCKS = 1
      END IF
*
*     Determine if the workspace size satisfies minimum size
*  
      LMINWS = .FALSE.
      IF((LWORK1.LT.MAX(1, NB*N*NBLCKS+5) 
     $    .OR.(LWORK2.LT.NB*N)).AND.(LWORK2.GE.N).AND.(LWORK1.GT.N+5) 
     $    .AND.(.NOT.LQUERY)) THEN
        IF (LWORK1.LT.MAX(1, NB * N * NBLCKS+5)) THEN
            LMINWS = .TRUE.
            NB = 1
        END IF  
        IF (LWORK1.LT.MAX(1, N * NBLCKS+5)) THEN
            LMINWS = .TRUE.
            MB = M 
        END IF
        IF (LWORK2.LT.NB*N) THEN
            LMINWS = .TRUE.
            NB = 1
        END IF
      END IF
*
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      ELSE IF( LWORK1.LT.MAX( 1, NB * N * NBLCKS + 5 ) 
     $   .AND.(.NOT.LQUERY).AND.(.NOT.LMINWS)) THEN
        INFO = -6
      ELSE IF( (LWORK2.LT.MAX(1,N*NB)).AND.(.NOT.LQUERY) 
     $   .AND.(.NOT.LMINWS)) THEN
        INFO = -8 
      END IF    

      IF( INFO.EQ.0)  THEN
        WORK1(1) = 1
        WORK1(2) = NB * N * NBLCKS + 5
        WORK1(3) = MINLW1
        WORK1(4) = MB
        WORK1(5) = NB
        WORK2(1) = NB * N
        WORK2(2) = N
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'SGEQR', -INFO )
        RETURN
      ELSE IF (LQUERY) THEN
       RETURN
      END IF
*
*     Quick return if possible
*
      IF( MIN(M,N).EQ.0 ) THEN
          RETURN
      END IF
*
*     The QR Decomposition
*
      IF((M.LE.N).OR.(MB.LE.N).OR.(MB.GE.M)) THEN
         CALL SGEQRT( M, N, NB, A, LDA, WORK1(6), NB, WORK2, INFO)
      ELSE 
         CALL SLATSQR( M, N, MB, NB, A, LDA, WORK1(6), NB, WORK2, 
     $                    LWORK2, INFO)
      END IF
      RETURN
*     
*     End of SGEQR
*
      END