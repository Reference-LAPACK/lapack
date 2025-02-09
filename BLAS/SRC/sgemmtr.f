*> \brief \b SGEMMTR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGEMMTR(UPLO,TRANSA,TRANSB,N,K,ALPHA,A,LDA,B,LDB,BETA,
*                         C,LDC)
*
*       .. Scalar Arguments ..
*       REAL ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,N
*       CHARACTER TRANSA,TRANSB, UPLO
*       ..
*       .. Array Arguments ..
*       REAL A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGEMMTR  performs one of the matrix-matrix operations
*>
*>    C := alpha*op( A )*op( B ) + beta*C,
*>
*> where  op( X ) is one of
*>
*>    op( X ) = X   or   op( X ) = X**T,
*>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
*> an n by k matrix,  op( B )  a  k by n matrix and  C an n by n matrix.
*> Thereby, the routine only accesses and updates the upper or lower
*> triangular part of the result matrix C. This behaviour can be used if
*> the resulting matrix C is known to be symmetric.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the lower or the upper
*>           triangular part of C is access and updated.
*>
*>              UPLO = 'L' or 'l', the lower triangular part of C is used.
*>
*>              UPLO = 'U' or 'u', the upper triangular part of C is used.
*> \endverbatim
*
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n',  op( A ) = A.
*>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c',  op( A ) = A**T.
*> \endverbatim
*>
*> \param[in] TRANSB
*> \verbatim
*>          TRANSB is CHARACTER*1
*>           On entry, TRANSB specifies the form of op( B ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSB = 'N' or 'n',  op( B ) = B.
*>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.
*>
*>              TRANSB = 'C' or 'c',  op( B ) = B**T.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N specifies the number of rows and columns of
*>           the matrix C, the number of columns of op(B) and the number
*>           of rows of op(A).  N must be at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry,  K  specifies  the number of columns of the matrix
*>           op( A ) and the number of rows of the matrix op( B ). K must
*>           be at least  zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is REAL.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension ( LDA, ka ), where ka is
*>           k  when  TRANSA = 'N' or 'n',  and is  n  otherwise.
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  n by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by m  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*>           LDA must be at least  max( 1, n ), otherwise  LDA must be at
*>           least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is REAL array, dimension ( LDB, kb ), where kb is
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*>           part of the array  B  must contain the matrix  B,  otherwise
*>           the leading  n by k  part of the array  B  must contain  the
*>           matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*>           least  max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is REAL.
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*>           supplied as zero then C need not be set on input.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is REAL array, dimension ( LDC, N )
*>           Before entry, the leading  n by n  part of the array  C must
*>           contain the matrix  C,  except when  beta  is zero, in which
*>           case C need not be set on entry.
*>           On exit, the upper or lower triangular part of the matrix
*>           C  is overwritten by the n by n matrix
*>           ( alpha*op( A )*op( B ) + beta*C ).
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, n ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Martin Koehler
*
*> \ingroup gemmtr
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 19-July-2023.
*>     Martin Koehler, MPI Magdeburg
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE SGEMMTR(UPLO,TRANSA,TRANSB,N,K,ALPHA,A,LDA,B,LDB,
     +         BETA,C,LDC)
      IMPLICIT NONE
*
*  -- Reference BLAS level3 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANSA,TRANSB,UPLO
*     ..
*     .. Array Arguments ..
      REAL A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      REAL TEMP
      INTEGER I,INFO,J,L,NROWA,NROWB, ISTART, ISTOP
      LOGICAL NOTA,NOTB, UPPER
*     ..
*     .. Parameters ..
      REAL ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA and NROWB  as the number of rows of  A
*     and  B  respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      IF (NOTB) THEN
          NROWB = K
      ELSE
          NROWB = N
      END IF
      UPPER = LSAME(UPLO, 'U')
*
*     Test the input parameters.
*
      INFO = 0
      IF ((.NOT. UPPER) .AND. (.NOT. LSAME(UPLO, 'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND.
     +    (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 2
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.
     +         (.NOT.LSAME(TRANSB,'T'))) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
          INFO = 10
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 13
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('SGEMMTR',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
*     And if  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  IF (UPPER) THEN
                      ISTART = 1
                      ISTOP  = J
                  ELSE
                      ISTART = J
                      ISTOP  = N
                  END IF

                  DO 10 I = ISTART, ISTOP
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  IF (UPPER) THEN
                      ISTART = 1
                      ISTOP  = J
                  ELSE
                      ISTART = J
                      ISTOP  = N
                  END IF

                  DO 30 I = ISTART, ISTOP
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (NOTB) THEN
          IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 J = 1,N
                  IF (UPPER) THEN
                      ISTART = 1
                      ISTOP  = J
                  ELSE
                      ISTART = J
                      ISTOP  = N
                  END IF
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = ISTART, ISTOP
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = ISTART, ISTOP
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      TEMP = ALPHA*B(L,J)
                      DO 70 I = ISTART, ISTOP
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B + beta*C
*
              DO 120 J = 1,N
                  IF (UPPER) THEN
                      ISTART = 1
                      ISTOP  = J
                  ELSE
                      ISTART = J
                      ISTOP  = N
                  END IF

                  DO 110 I = ISTART, ISTOP
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (NOTA) THEN
*
*           Form  C := alpha*A*B**T + beta*C
*
              DO 170 J = 1,N
                  IF (UPPER) THEN
                      ISTART = 1
                      ISTOP  = J
                  ELSE
                      ISTART = J
                      ISTOP  = N
                  END IF

                  IF (BETA.EQ.ZERO) THEN
                      DO 130 I = ISTART,ISTOP
                          C(I,J) = ZERO
  130                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 140 I = ISTART,ISTOP
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 150 I = ISTART,ISTOP
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B**T + beta*C
*
              DO 200 J = 1,N
                  IF (UPPER) THEN
                      ISTART = 1
                      ISTOP  = J
                  ELSE
                      ISTART = J
                      ISTOP  = N
                  END IF

                  DO 190 I = ISTART, ISTOP
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of SGEMMTR
*
      END
