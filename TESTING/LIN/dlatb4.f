*> \brief \b DLATB4
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE,
*                          CNDNUM, DIST )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIST, TYPE
*       CHARACTER*3        PATH
*       INTEGER            IMAT, KL, KU, M, MODE, N
*       DOUBLE PRECISION   ANORM, CNDNUM
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLATB4 sets parameters for the matrix generator based on the type of
*> matrix to be generated.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] PATH
*> \verbatim
*>          PATH is CHARACTER*3
*>          The LAPACK path name.
*> \endverbatim
*>
*> \param[in] IMAT
*> \verbatim
*>          IMAT is INTEGER
*>          An integer key describing which matrix to generate for this
*>          path.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows in the matrix to be generated.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns in the matrix to be generated.
*> \endverbatim
*>
*> \param[out] TYPE
*> \verbatim
*>          TYPE is CHARACTER*1
*>          The type of the matrix to be generated:
*>          = 'S':  symmetric matrix
*>          = 'P':  symmetric positive (semi)definite matrix
*>          = 'N':  nonsymmetric matrix
*> \endverbatim
*>
*> \param[out] KL
*> \verbatim
*>          KL is INTEGER
*>          The lower band width of the matrix to be generated.
*> \endverbatim
*>
*> \param[out] KU
*> \verbatim
*>          KU is INTEGER
*>          The upper band width of the matrix to be generated.
*> \endverbatim
*>
*> \param[out] ANORM
*> \verbatim
*>          ANORM is DOUBLE PRECISION
*>          The desired norm of the matrix to be generated.  The diagonal
*>          matrix of singular values or eigenvalues is scaled by this
*>          value.
*> \endverbatim
*>
*> \param[out] MODE
*> \verbatim
*>          MODE is INTEGER
*>          A key indicating how to choose the vector of eigenvalues.
*> \endverbatim
*>
*> \param[out] CNDNUM
*> \verbatim
*>          CNDNUM is DOUBLE PRECISION
*>          The desired condition number.
*> \endverbatim
*>
*> \param[out] DIST
*> \verbatim
*>          DIST is CHARACTER*1
*>          The type of distribution to be used by the random number
*>          generator.
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
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE DLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE,
     $                   CNDNUM, DIST )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIST, TYPE
      CHARACTER*3        PATH
      INTEGER            IMAT, KL, KU, M, MODE, N
      DOUBLE PRECISION   ANORM, CNDNUM
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   SHRINK, TENTH
      PARAMETER          ( SHRINK = 0.25D+0, TENTH = 0.1D+0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST
      CHARACTER*2        C2
      INTEGER            MAT
      DOUBLE PRECISION   BADC1, BADC2, EPS, LARGE, SMALL
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAMEN, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Save statement ..
      SAVE               EPS, SMALL, LARGE, BADC1, BADC2, FIRST
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     Set some constants for use in the subroutine.
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         EPS = DLAMCH( 'Precision' )
         BADC2 = TENTH / EPS
         BADC1 = SQRT( BADC2 )
         SMALL = DLAMCH( 'Safe minimum' )
         LARGE = ONE / SMALL
         SMALL = SHRINK*( SMALL / EPS )
         LARGE = ONE / SMALL
      END IF
*
      C2 = PATH( 2: 3 )
*
*     Set some parameters we don't plan to change.
*
      DIST = 'S'
      MODE = 3
*
      IF( LSAMEN( 2, C2, 'QR' ) .OR. LSAMEN( 2, C2, 'LQ' ) .OR.
     $    LSAMEN( 2, C2, 'QL' ) .OR. LSAMEN( 2, C2, 'RQ' ) ) THEN
*
*        xQR, xLQ, xQL, xRQ:  Set parameters to generate a general
*                             M x N matrix.
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
*        Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
            KU = 0
         ELSE IF( IMAT.EQ.2 ) THEN
            KL = 0
            KU = MAX( N-1, 0 )
         ELSE IF( IMAT.EQ.3 ) THEN
            KL = MAX( M-1, 0 )
            KU = 0
         ELSE
            KL = MAX( M-1, 0 )
            KU = MAX( N-1, 0 )
         END IF
*
*        Set the condition number and norm.
*
         IF( IMAT.EQ.5 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.6 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.7 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.8 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'QK' ) ) THEN
*
*        xQK: truncated QR with pivoting.
*             Set parameters to generate a general
*             M x N matrix.
*
*        Set TYPE, the type of matrix to be generated.  'N' is nonsymmetric.
*
         TYPE = 'N'
*
*        Set DIST, the type of distribution for the random
*        number generator. 'S' is
*
         DIST = 'S'
*
*        Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.2 ) THEN
*
*           2. Random, Diagonal, CNDNUM = 2
*
            KL = 0
            KU = 0
            CNDNUM = TWO
            ANORM = ONE
            MODE = 3
         ELSE IF( IMAT.EQ.3 ) THEN
*
*           3. Random, Upper triangular,  CNDNUM = 2
*
            KL = 0
            KU = MAX( N-1, 0 )
            CNDNUM = TWO
            ANORM = ONE
            MODE = 3
         ELSE IF( IMAT.EQ.4 ) THEN
*
*          4. Random, Lower triangular,  CNDNUM = 2
*
            KL = MAX( M-1, 0 )
            KU = 0
            CNDNUM = TWO
            ANORM = ONE
            MODE = 3
         ELSE
*
*           5.-19. Rectangular matrix
*
            KL = MAX( M-1, 0 )
            KU = MAX( N-1, 0 )
*
            IF( IMAT.GE.5 .AND. IMAT.LE.14 ) THEN
*
*              5.-14. Random, CNDNUM = 2.
*
               CNDNUM = TWO
               ANORM = ONE
               MODE = 3
*
            ELSE IF( IMAT.EQ.15 ) THEN
*
*              15. Random, CNDNUM = sqrt(0.1/EPS)
*
               CNDNUM = BADC1
               ANORM = ONE
               MODE = 3
*
            ELSE IF( IMAT.EQ.16 ) THEN
*
*              16. Random, CNDNUM = 0.1/EPS
*
               CNDNUM = BADC2
               ANORM = ONE
               MODE = 3
*
            ELSE IF( IMAT.EQ.17 ) THEN
*
*              17. Random, CNDNUM = 0.1/EPS,
*                  one small singular value S(N)=1/CNDNUM
*
               CNDNUM = BADC2
               ANORM = ONE
               MODE = 2
*
            ELSE IF( IMAT.EQ.18 ) THEN
*
*              18. Random, scaled near underflow
*
               CNDNUM = TWO
               ANORM = SMALL
               MODE = 3
*
            ELSE IF( IMAT.EQ.19 ) THEN
*
*              19. Random, scaled near overflow
*
               CNDNUM = TWO
               ANORM = LARGE
               MODE = 3
*
            END IF
*
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'GE' ) ) THEN
*
*        xGE:  Set parameters to generate a general M x N matrix.
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
*        Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
            KU = 0
         ELSE IF( IMAT.EQ.2 ) THEN
            KL = 0
            KU = MAX( N-1, 0 )
         ELSE IF( IMAT.EQ.3 ) THEN
            KL = MAX( M-1, 0 )
            KU = 0
         ELSE
            KL = MAX( M-1, 0 )
            KU = MAX( N-1, 0 )
         END IF
*
*        Set the condition number and norm.
*
         IF( IMAT.EQ.8 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.9 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.10 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.11 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN
*
*        xGB:  Set parameters to generate a general banded matrix.
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
*        Set the condition number and norm.
*
         IF( IMAT.EQ.5 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.6 ) THEN
            CNDNUM = TENTH*BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.7 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.8 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'GT' ) ) THEN
*
*        xGT:  Set parameters to generate a general tridiagonal matrix.
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
*        Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
         ELSE
            KL = 1
         END IF
         KU = KL
*
*        Set the condition number and norm.
*
         IF( IMAT.EQ.3 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.4 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.5 .OR. IMAT.EQ.11 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.6 .OR. IMAT.EQ.12 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'PP' ) ) THEN
*
*        xPO, xPP: Set parameters to generate a
*        symmetric positive definite matrix.
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = C2( 1: 1 )
*
*        Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
         ELSE
            KL = MAX( N-1, 0 )
         END IF
         KU = KL
*
*        Set the condition number and norm.
*
         IF( IMAT.EQ.6 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.7 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.8 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.9 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
*
      ELSE IF( LSAMEN( 2, C2, 'SY' ) .OR. LSAMEN( 2, C2, 'SP' ) ) THEN
*
*        xSY, xSP: Set parameters to generate a
*        symmetric matrix.
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = C2( 1: 1 )
*
*        Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
         ELSE
            KL = MAX( N-1, 0 )
         END IF
         KU = KL
*
*        Set the condition number and norm.
*
         IF( IMAT.EQ.7 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.8 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.9 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.10 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'KY' ) ) THEN
*
*        xKY: Set parameters to generate a
*        skew-symmetric matrix.
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = C2( 1: 1 )
*
*        Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
         ELSE
            KL = MAX( N-1, 0 )
         END IF
         KU = KL
*
*        Set the condition number and norm.
*
         IF( IMAT.EQ.7 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.8 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.9 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.10 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'PB' ) ) THEN
*
*        xPB:  Set parameters to generate a symmetric band matrix.
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = 'P'
*
*        Set the norm and condition number.
*
         IF( IMAT.EQ.5 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.6 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.7 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.8 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'PT' ) ) THEN
*
*        xPT:  Set parameters to generate a symmetric positive definite
*        tridiagonal matrix.
*
         TYPE = 'P'
         IF( IMAT.EQ.1 ) THEN
            KL = 0
         ELSE
            KL = 1
         END IF
         KU = KL
*
*        Set the condition number and norm.
*
         IF( IMAT.EQ.3 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.4 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.5 .OR. IMAT.EQ.11 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.6 .OR. IMAT.EQ.12 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'TR' ) .OR. LSAMEN( 2, C2, 'TP' ) ) THEN
*
*        xTR, xTP:  Set parameters to generate a triangular matrix
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
*        Set the lower and upper bandwidths.
*
         MAT = ABS( IMAT )
         IF( MAT.EQ.1 .OR. MAT.EQ.7 ) THEN
            KL = 0
            KU = 0
         ELSE IF( IMAT.LT.0 ) THEN
            KL = MAX( N-1, 0 )
            KU = 0
         ELSE
            KL = 0
            KU = MAX( N-1, 0 )
         END IF
*
*        Set the condition number and norm.
*
         IF( MAT.EQ.3 .OR. MAT.EQ.9 ) THEN
            CNDNUM = BADC1
         ELSE IF( MAT.EQ.4 ) THEN
            CNDNUM = BADC2
         ELSE IF( MAT.EQ.10 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( MAT.EQ.5 ) THEN
            ANORM = SMALL
         ELSE IF( MAT.EQ.6 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'TB' ) ) THEN
*
*        xTB:  Set parameters to generate a triangular band matrix.
*
*        Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
*        Set the norm and condition number.
*
         MAT = ABS( IMAT )
         IF( MAT.EQ.2 .OR. MAT.EQ.8 ) THEN
            CNDNUM = BADC1
         ELSE IF( MAT.EQ.3 .OR. MAT.EQ.9 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( MAT.EQ.4 ) THEN
            ANORM = SMALL
         ELSE IF( MAT.EQ.5 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
      END IF
      IF( N.LE.1 )
     $   CNDNUM = ONE
*
      RETURN
*
*     End of DLATB4
*
      END
