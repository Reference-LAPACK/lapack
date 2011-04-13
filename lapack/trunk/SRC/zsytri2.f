      SUBROUTINE ZSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*  -- Written by Julie Langou of the Univ. of TN    --
*
* @precisions normal z -> s d c
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZSYTRI2 computes the inverse of a COMPLEX*16 hermitian indefinite matrix
*  A using the factorization A = U*D*U**T or A = L*D*L**T computed by
*  ZSYTRF. ZSYTRI2 sets the LEADING DIMENSION of the workspace
*  before calling ZSYTRI2X that actually computes the inverse.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the details of the factorization are stored
*          as an upper or lower triangular matrix.
*          = 'U':  Upper triangular, form is A = U*D*U**T;
*          = 'L':  Lower triangular, form is A = L*D*L**T.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the NB diagonal matrix D and the multipliers
*          used to obtain the factor U or L as computed by ZSYTRF.
*
*          On exit, if INFO = 0, the (symmetric) inverse of the original
*          matrix.  If UPLO = 'U', the upper triangular part of the
*          inverse is formed and the part of A below the diagonal is not
*          referenced; if UPLO = 'L' the lower triangular part of the
*          inverse is formed and the part of A above the diagonal is
*          not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          Details of the interchanges and the NB structure of D
*          as determined by ZSYTRF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N+NB+1)*(NB+3)
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          WORK is size >= (N+NB+1)*(NB+3)
*          If LDWORK = -1, then a workspace query is assumed; the routine
*           calculates:
*              - the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array,
*              - and no error message related to LDWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
*               inverse could not be computed.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            UPPER, LQUERY
      INTEGER            MINSIZE, NBMAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZSYTRI2X
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
*     Get blocksize
      NBMAX = ILAENV( 1, 'ZSYTRF', UPLO, N, -1, -1, -1 )
      IF ( NBMAX .GE. N ) THEN
         MINSIZE = N
      ELSE
         MINSIZE = (N+NBMAX+1)*(NBMAX+3)
      END IF
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF (LWORK .LT. MINSIZE .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
*
*     Quick return if possible
*
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYTRI2', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         WORK(1)=MINSIZE
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      
      IF( NBMAX .GE. N ) THEN
         CALL ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
      ELSE
         CALL ZSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NBMAX, INFO )
      END IF
      RETURN
*
*     End of ZSYTRI2
*
      END
