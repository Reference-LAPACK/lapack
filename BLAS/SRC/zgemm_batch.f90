!> \brief \b ZGEMM_BATCH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEMM_BATCH(TRANSA_ARRAY, TRANSB_ARRAY,
!                              M_ARRAY, N_ARRAY, K_ARRAY,
!                              ALPHA_ARRAY,
!                              A_ARRAY, LDA_ARRAY,
!                              B_ARRAY, LDB_ARRAY,
!                              BETA_ARRAY,
!                              C_ARRAY, LDC_ARRAY,
!                              GROUP_COUNT, GROUP_SIZE)
!
!       .. Scalar Arguments ..
!       INTEGER GROUP_COUNT
!       ..
!       .. Array Arguments ..
!       CHARACTER TRANSA_ARRAY(GROUP_COUNT), TRANSB_ARRAY(GROUP_COUNT)
!       INTEGER M_ARRAY(GROUP_COUNT), N_ARRAY(GROUP_COUNT), K_ARRAY(GROUP_COUNT)
!       COMPLEX*16 ALPHA_ARRAY(GROUP_COUNT),BETA_ARRAY(GROUP_COUNT)
!       INTEGER LDA_ARRAY(GROUP_COUNT), LDB_ARRAY(GROUP_COUNT), LDC_ARRAY(GROUP_COUNT)
!       INTEGER GROUP_SIZE(GROUP_COUNT)
!       ..
!       .. Pointer Arguments ..
!       TYPE(C_PTR) A_ARRAY(*), B_ARRAY(*), C_ARRAY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEMM_BATCH  performs a series of the matrix-matrix operations with each ji'th matrix:
!>
!>    C_ji := alpha_i*op_i( A_ji )*op( B_ji ) + beta_i*C_ji,
!>
!> where  op_i( X ) is one of
!>
!>    op_i( X_ji ) = X_ji   or   op_i( X_ji ) = X_ji**T,
!>
!> alpha_i and beta_i are scalars, and A_ji, B_ji and C_ji are matrices, with op_i( A_ji )
!> an m_i by k_i matrix,  op_i( B_ji )  a  k_i by n_i matrix and  C_ji an m_i by n_i matrix.
!> Group count defines i and group_size(i) defines j.
!>
!> More generally,
!>
!>    idx = 1
!>    for i in 1..group_count
!>      alpha, beta = alpha(i), beta(i)
!>      for j in 1..group_size(i)
!>        A, B, C = A_ARRAY(idx), B_ARRAY(idx), C_ARRAY(idx)
!>        C := alpha*op(A)*op(B) + beta*C
!>        idx = idx + 1
!>
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA_ARRAY
!> \verbatim
!>          TRANSA_ARRAY is CHARACTER*1 array
!>           On entry, TRANSA_ARRAY(i) specifies the form of op_i( A_ji ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA_ARRAY(i) = 'N' or 'n',  op_i( A_ji ) = A_ji.
!>
!>              TRANSA_ARRAY(i) = 'T' or 't',  op_i( A_ji ) = A_ji**T.
!>
!>              TRANSA_ARRAY(i) = 'C' or 'c',  op_i( A_ji ) = A_ji**H.
!> \endverbatim
!>
!> \param[in] TRANSB_ARRAY
!> \verbatim
!>          TRANSB_ARRAY is CHARACTER*1 array
!>           On entry, TRANSB_ARRAY(i) specifies the form of op_i( B_ji ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSB_ARRAY(i) = 'N' or 'n',  op_i( B_ji ) = B_ji.
!>
!>              TRANSB_ARRAY(i) = 'T' or 't',  op_i( B_ji ) = B_ji**T.
!>
!>              TRANSB_ARRAY(i) = 'C' or 'c',  op_i( B_ji ) = B_ji**H.
!> \endverbatim
!>
!> \param[in] M_ARRAY
!> \verbatim
!>          M_ARRAY is INTEGER array
!>           On entry,  M_ARRAY(i)  specifies  the number  of rows  of the  matrixes
!>           op_i( A_ji )  and the number of rows of the matrixes  C_ji.
!>           Each M_ARRAY(i)  must  be at least  zero.
!> \endverbatim
!>
!> \param[in] N_ARRAY
!> \verbatim
!>          N_ARRAY is INTEGER array
!>           On entry,  N_ARRAY(i)  specifies the number  of columns of the matrixes
!>           op_i( B_ji ) and the number of columns of the matrixes C_ji.
!>           Each N_ARRAY(i) must be at least zero.
!> \endverbatim
!>
!> \param[in] K_ARRAY
!> \verbatim
!>          K_ARRAY is INTEGER array
!>           On entry,  K_ARRAY(i)  specifies  the number of columns of the matrixes
!>           op_i( A_ji ) and the number of rows of the matrixes op_i( B_ji ).
!>           Each K_ARRAY(i) must be at least  zero.
!> \endverbatim
!>
!> \param[in] ALPHA_ARRAY
!> \verbatim
!>          ALPHA_ARRAY is COMPLEX*16 array.
!>           On entry, ALPHA_ARRAY(i) specifies the scalar alpha_i.
!> \endverbatim
!>
!> \param[in] A_ARRAY
!> \verbatim
!>          A_ARRAY is POINTER array, dimension ( sum( GROUP_SIZE ) ),
!>           to COMPLEX*16 arrays, dimension ( LDA_i, ka_i ),
!>           where ka_i is k_i  when  TRANSA(i) = 'N' or 'n',  and is  m_i  otherwise.
!>           Before entry with  TRANSA = 'N' or 'n',  the leading  m_i by k_i elements
!>           at address A(ji)  must contain the matrix  A_ji,  otherwise
!>           the leading  k_i by m_i elements at address A(ji)  must contain  the
!>           matrix A_ji.
!> \endverbatim
!>
!> \param[in] LDA_ARRAY
!> \verbatim
!>          LDA_ARRAY is INTEGER array
!>           On entry, LDA_ARRAY(i) specifies the first dimension of A_ji as declared
!>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!>           LDA_ARRAY(i) must be at least  max( 1, m_i ), otherwise  LDA must be at
!>           least  max( 1, k_i ).
!> \endverbatim
!>
!> \param[in] B_ARRAY
!> \verbatim
!>          B_ARRAY is POINTER array, dimension ( sum( GROUP_SIZE ) ),
!>           to COMPLEX*16 arrays, dimension ( LDB_i, kb_i ),
!>           where kb_i is n_i  when  TRANSB(i) = 'N' or 'n',  and is  k_i  otherwise.
!>           Before entry with  TRANSB = 'N' or 'n',  the leading  k_i by n_i elements
!>           at address B(ji)  must contain the matrix  B_ji,  otherwise
!>           the leading  n_i by k_i elements at address B(ji)  must contain  the
!>           matrix B_ji.
!> \endverbatim
!>
!> \param[in] LDB_ARRAY
!> \verbatim
!>          LDB_ARRAY is INTEGER array
!>           On entry, LDB_ARRAY(i) specifies the first dimension of B_ji as declared
!>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!>           LDB must be at least  max( 1, k_i ), otherwise  LDB must be at
!>           least  max( 1, n_i ).
!> \endverbatim
!>
!> \param[in] BETA_ARRAY
!> \verbatim
!>          BETA_ARRAY is COMPLEX*16 array.
!>           On entry,  BETA_ARRAY(i)  specifies the scalar  beta.  When  BETA_ARRAY(i)  is
!>           supplied as zero then C_ji need not be set on input.
!> \endverbatim
!>
!> \param[in,out] C_ARRAY
!> \verbatim
!>          C_ARRAY is POINTER array, dimension ( sum( GROUP_SIZE ) ),
!>           to COMPLEX*16 arrays, dimension ( LDC_i, n_i ).
!>           Before entry,  the leading  m_i by n_i elements
!>           at address C(ji)  must contain the matrix  C_ji,  except when BETA_ARRAY(i)
!>           is zero, in which case C_ji need not be set on entry.
!>           On exit, the array  C_ji  is overwritten by the  m_i by n_i  matrix
!>           ( alpha_i*op_i( A_ji )*op_i( B_ji ) + beta_i*C_ji ).
!> \endverbatim
!>
!> \param[in] LDC_ARRAY
!> \verbatim
!>          LDC_ARRAY is INTEGER array
!>           On entry, LDC_ARRAY(i) specifies the first dimension of C_ji as declared
!>           in  the  calling  (sub)  program.   LDC_ARRAY(i)  must  be  at  least
!>           max( 1, m_i ).
!> \endverbatim
!>
!> \param[in] GROUP_COUNT
!> \verbatim
!>          GROUP_COUNT is INTEGER
!>           On entry, GROUP_COUNT specifies the number of groups that determines index i.
!> \endverbatim
!>
!> \param[in] GROUP_SIZE
!> \verbatim
!>          GROUP_SIZE is INTEGER array
!>           On entry, GROUP_SIZE specifies the number of elements in each groups that determines index j.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Igor S. Gerasimov
!
!> \ingroup gemm_batch
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>  Original API is taken from:
!>     https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-2/gemm-batch.html
!>
!>  -- Written on 23-October-2023.
!>
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGEMM_BATCH(TRANSA_ARRAY, TRANSB_ARRAY, &
                             M_ARRAY, N_ARRAY, K_ARRAY, &
                             ALPHA_ARRAY, &
                             A_ARRAY, LDA_ARRAY, &
                             B_ARRAY, LDB_ARRAY, &
                             BETA_ARRAY, &
                             C_ARRAY, LDC_ARRAY, &
                             GROUP_COUNT, GROUP_SIZE)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_ASSOCIATED
!
!  -- Reference BLAS level3 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER GROUP_COUNT
!     ..
!     .. Array Arguments ..
      CHARACTER TRANSA_ARRAY(GROUP_COUNT), TRANSB_ARRAY(GROUP_COUNT)
      INTEGER M_ARRAY(GROUP_COUNT), N_ARRAY(GROUP_COUNT), K_ARRAY(GROUP_COUNT)
      COMPLEX*16 ALPHA_ARRAY(GROUP_COUNT), BETA_ARRAY(GROUP_COUNT)
      INTEGER LDA_ARRAY(GROUP_COUNT), LDB_ARRAY(GROUP_COUNT), LDC_ARRAY(GROUP_COUNT)
      INTEGER GROUP_SIZE(GROUP_COUNT)
!     ..
!     .. Pointer Arguments ..
      TYPE(C_PTR) A_ARRAY(*), B_ARRAY(*), C_ARRAY(*)
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
      EXTERNAL XERBLAI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      INTEGER I, J, IDX, INFO
      LOGICAL NOTA, NOTB
      INTEGER NROWA, NROWB
!     ..
!     .. Local Addresses ..
      COMPLEX*16, POINTER :: A, B, C
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (GROUP_COUNT.LT.0) THEN
        INFO = 15
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('ZGEMM_BATCH ', INFO)
        RETURN
      END IF
      DO I = 1, GROUP_COUNT
        INFO = 0
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA and NROWB  as the number of rows of  A
!     and  B  respectively.
!
        NOTA = LSAME(TRANSA_ARRAY(I),'N')
        NOTB = LSAME(TRANSB_ARRAY(I),'N')
        IF (NOTA) THEN
          NROWA = M_ARRAY(I)
        ELSE
          NROWA = K_ARRAY(I)
        END IF
        IF (NOTB) THEN
          NROWB = K_ARRAY(I)
        ELSE
          NROWB = N_ARRAY(I)
        END IF
        IF ((.NOT.NOTA) .AND. &
            (.NOT.LSAME(TRANSA_ARRAY(I),'C')) .AND. &
            (.NOT.LSAME(TRANSA_ARRAY(I),'T'))) THEN
          INFO = 1
        ELSE IF ((.NOT.NOTB) .AND. &
               (.NOT.LSAME(TRANSB_ARRAY(I),'C')) .AND. &
               (.NOT.LSAME(TRANSB_ARRAY(I),'T'))) THEN
          INFO = 2
        ELSE IF (M_ARRAY(I).LT.0) THEN
          INFO = 3
        ELSE IF (N_ARRAY(I).LT.0) THEN
          INFO = 4
        ELSE IF (K_ARRAY(I).LT.0) THEN
          INFO = 5
        ELSE IF (LDA_ARRAY(I).LT.MAX(1,NROWA)) THEN
          INFO = 8
        ELSE IF (LDB_ARRAY(I).LT.MAX(1,NROWB)) THEN
          INFO = 10
        ELSE IF (LDC_ARRAY(I).LT.MAX(1,M_ARRAY(I))) THEN
          INFO = 13
        ELSE IF (GROUP_SIZE(I).LT.0) THEN
          INFO = 15
        END IF
        IF (INFO.NE.0) THEN
          CALL XERBLAI('ZGEMM_BATCH ',INFO,I)
          RETURN
        END IF
      END DO
      IDX = 1
      DO I = 1, GROUP_COUNT
        DO J = 1, GROUP_SIZE(I)
          INFO = 0
          IF (.NOT.C_ASSOCIATED(A_ARRAY(IDX))) THEN
            INFO = 7
          ELSE IF (.NOT.C_ASSOCIATED(B_ARRAY(IDX))) THEN
            INFO = 9
          ELSE IF (.NOT.C_ASSOCIATED(C_ARRAY(IDX))) THEN
            INFO = 12
          END IF
          IF (INFO.NE.0) THEN
            CALL XERBLAI('ZGEMM_BATCH ',INFO,IDX)
            RETURN
          END IF
          IDX = IDX + 1
        END DO
      END DO
!
!     Do computations.
!
      IDX = 1
      DO I = 1, GROUP_COUNT
        DO J = 1, GROUP_SIZE(I)
          CALL C_F_POINTER(A_ARRAY(IDX), A)
          CALL C_F_POINTER(B_ARRAY(IDX), B)
          CALL C_F_POINTER(C_ARRAY(IDX), C)
          CALL ZGEMM(TRANSA_ARRAY(I), TRANSB_ARRAY(I), &
                     M_ARRAY(I), N_ARRAY(I), K_ARRAY(I), &
                     ALPHA_ARRAY(I), &
                     A, LDA_ARRAY(I), &
                     B, LDB_ARRAY(I), &
                     BETA_ARRAY(I), &
                     C, LDC_ARRAY(I))
          IDX = IDX + 1
        END DO
      END DO
      RETURN
!
!     End of ZGEMM_BATCH.
!
      END
