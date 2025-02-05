*> \brief \b DLA_GERPVGRW
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DLA_GERPVGRW + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_gerpvgrw.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_gerpvgrw.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_gerpvgrw.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLA_GERPVGRW( N, NCOLS, A, LDA, AF,
*                                               LDAF )
*
*       .. Scalar Arguments ..
*       INTEGER            N, NCOLS, LDA, LDAF
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>
*> DLA_GERPVGRW computes the reciprocal pivot growth factor
*> norm(A)/norm(U). The "max absolute element" norm is used. If this is
*> much less than 1, the stability of the LU factorization of the
*> (equilibrated) matrix A could be poor. This also means that the
*> solution X, estimated condition numbers, and error bounds could be
*> unreliable.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>     The number of linear equations, i.e., the order of the
*>     matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NCOLS
*> \verbatim
*>          NCOLS is INTEGER
*>     The number of columns of the matrix A. NCOLS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>     On entry, the N-by-N matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>     The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] AF
*> \verbatim
*>          AF is DOUBLE PRECISION array, dimension (LDAF,N)
*>     The factors L and U from the factorization
*>     A = P*L*U as computed by DGETRF.
*> \endverbatim
*>
*> \param[in] LDAF
*> \verbatim
*>          LDAF is INTEGER
*>     The leading dimension of the array AF.  LDAF >= max(1,N).
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
*> \ingroup la_gerpvgrw
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DLA_GERPVGRW( N, NCOLS, A, LDA, AF,
     $         LDAF )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            N, NCOLS, LDA, LDAF
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   AMAX, UMAX, RPVGRW
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      RPVGRW = 1.0D+0

      DO J = 1, NCOLS
         AMAX = 0.0D+0
         UMAX = 0.0D+0
         DO I = 1, N
            AMAX = MAX( ABS( A( I, J ) ), AMAX )
         END DO
         DO I = 1, J
            UMAX = MAX( ABS( AF( I, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0D+0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      DLA_GERPVGRW = RPVGRW
*
*     End of DLA_GERPVGRW
*
      END
