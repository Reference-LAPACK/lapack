*> \brief \b ZLAQHP scales a Hermitian matrix stored in packed form.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZLAQHP + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqhp.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqhp.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqhp.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAQHP( UPLO, N, AP, S, SCOND, AMAX, EQUED )
*
*       .. Scalar Arguments ..
*       CHARACTER          EQUED, UPLO
*       INTEGER            N
*       DOUBLE PRECISION   AMAX, SCOND
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   S( * )
*       COMPLEX*16         AP( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLAQHP equilibrates a Hermitian matrix A using the scaling factors
*> in the vector S.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          Hermitian matrix A is stored.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] AP
*> \verbatim
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
*>          On entry, the upper or lower triangle of the Hermitian matrix
*>          A, packed columnwise in a linear array.  The j-th column of A
*>          is stored in the array AP as follows:
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
*>
*>          On exit, the equilibrated matrix:  diag(S) * A * diag(S), in
*>          the same storage format as A.
*> \endverbatim
*>
*> \param[in] S
*> \verbatim
*>          S is DOUBLE PRECISION array, dimension (N)
*>          The scale factors for A.
*> \endverbatim
*>
*> \param[in] SCOND
*> \verbatim
*>          SCOND is DOUBLE PRECISION
*>          Ratio of the smallest S(i) to the largest S(i).
*> \endverbatim
*>
*> \param[in] AMAX
*> \verbatim
*>          AMAX is DOUBLE PRECISION
*>          Absolute value of largest matrix entry.
*> \endverbatim
*>
*> \param[out] EQUED
*> \verbatim
*>          EQUED is CHARACTER*1
*>          Specifies whether or not equilibration was done.
*>          = 'N':  No equilibration.
*>          = 'Y':  Equilibration was done, i.e., A has been replaced by
*>                  diag(S) * A * diag(S).
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  THRESH is a threshold value used to decide if scaling should be done
*>  based on the ratio of the scaling factors.  If SCOND < THRESH,
*>  scaling is done.
*>
*>  LARGE and SMALL are threshold values used to decide if scaling should
*>  be done based on the absolute size of the largest matrix element.
*>  If AMAX > LARGE or AMAX < SMALL, scaling is done.
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
*> \ingroup laqhp
*
*  =====================================================================
      SUBROUTINE ZLAQHP( UPLO, N, AP, S, SCOND, AMAX, EQUED )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, UPLO
      INTEGER            N
      DOUBLE PRECISION   AMAX, SCOND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   S( * )
      COMPLEX*16         AP( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, THRESH
      PARAMETER          ( ONE = 1.0D+0, THRESH = 0.1D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JC
      DOUBLE PRECISION   CJ, LARGE, SMALL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Initialize LARGE and SMALL.
*
      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL
*
      IF( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) THEN
*
*        No equilibration
*
         EQUED = 'N'
      ELSE
*
*        Replace A by diag(S) * A * diag(S).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Upper triangle of A is stored.
*
            JC = 1
            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = 1, J - 1
                  AP( JC+I-1 ) = CJ*S( I )*AP( JC+I-1 )
   10          CONTINUE
               AP( JC+J-1 ) = CJ*CJ*DBLE( AP( JC+J-1 ) )
               JC = JC + J
   20       CONTINUE
         ELSE
*
*           Lower triangle of A is stored.
*
            JC = 1
            DO 40 J = 1, N
               CJ = S( J )
               AP( JC ) = CJ*CJ*DBLE( AP( JC ) )
               DO 30 I = J + 1, N
                  AP( JC+I-J ) = CJ*S( I )*AP( JC+I-J )
   30          CONTINUE
               JC = JC + N - J + 1
   40       CONTINUE
         END IF
         EQUED = 'Y'
      END IF
*
      RETURN
*
*     End of ZLAQHP
*
      END
