*> \brief \b SLATRN
*
* =========== DOCUMENTATION ===========
*
*> \par Purpose:
* =============
*>
*> \details \b Purpose:
*> \verbatim
*>
*> SLATRN performs an in-place transpose of an N-by-N single precision
*> matrix A. It utilizes OpenMP parallelization and a blocked iteration 
*> scheme to optimize cache usage during the transposition process.
*> \endverbatim
*
* Arguments:
* ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA, N)
*>          On entry, the N-by-N matrix to be transposed.
*>          On exit, A is overwritten by its transpose.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
      SUBROUTINE SLATRN( N, A, LDA )

#if defined(_OPENMP)
      use omp_lib
#endif
*
* -- LAPACK auxiliary routine --
*
* .. Scalar Arguments ..
      INTEGER            LDA, N
* ..
* .. Array Arguments ..
      REAL   A( LDA, * )
* ..
*
* =====================================================================
*
* .. Parameters ..
      INTEGER            SIZE
      PARAMETER          ( SIZE = 16 )
* ..
* .. Local Scalars ..
      INTEGER            BLOCK, I, J, NUM_THREADS, TAIL_START
      REAL   TEMP
* ..
* .. External Functions ..
#if defined(_OPENMP)
      INTEGER            OMP_GET_MAX_THREADS
      EXTERNAL           OMP_GET_MAX_THREADS
#endif
* ..
* .. Intrinsic Functions ..
      INTRINSIC          MIN
* ..
* .. Executable Statements ..
*
#if defined(_OPENMP)
      IF( N.GE.256 ) THEN
         NUM_THREADS = N / 64
      ELSE
         NUM_THREADS = 1
      END IF
      NUM_THREADS = MIN( NUM_THREADS, OMP_GET_MAX_THREADS() )
#else
      NUM_THREADS = 1
#endif
   
!$OMP PARALLEL DO PRIVATE( BLOCK, I, J, TEMP )
!$OMP+            NUM_THREADS( NUM_THREADS ) SCHEDULE( DYNAMIC, 1 )
      DO BLOCK = 1, N - SIZE + 1, SIZE
*
* This pair takes care of the main diagonal block
         DO I = BLOCK, BLOCK + SIZE - 1
            DO J = I + 1, BLOCK + SIZE - 1
               TEMP      = A( I, J )
               A( I, J ) = A( J, I )
               A( J, I ) = TEMP
            END DO
         END DO
*
* Transpose sub-matrices not on the main diagonal
         DO I = BLOCK + SIZE, N
            DO J = BLOCK, BLOCK + SIZE - 1
               TEMP      = A( I, J )
               A( I, J ) = A( J, I )
               A( J, I ) = TEMP
            END DO
         END DO
*
      END DO
!$OMP END PARALLEL DO
*
* Transpose remaining elements along the main diagonal
* Note: TAIL_START is explicitly calculated because BLOCK is 
* undefined outside the OpenMP parallel construct.
      TAIL_START = ( N / SIZE ) * SIZE + 1
      DO I = TAIL_START, N
         DO J = I + 1, N
            TEMP      = A( I, J )
            A( I, J ) = A( J, I )
            A( J, I ) = TEMP
         END DO
      END DO
*
      RETURN
*
* End of SLATRN
*
      END