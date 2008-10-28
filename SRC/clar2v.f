      SUBROUTINE CLAR2V( N, X, Y, Z, INCX, C, S, INCC )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCC, INCX, N
*     ..
*     .. Array Arguments ..
      REAL               C( * )
      COMPLEX            S( * ), X( * ), Y( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  CLAR2V applies a vector of complex plane rotations with real cosines
*  from both sides to a sequence of 2-by-2 complex Hermitian matrices,
*  defined by the elements of the vectors x, y and z. For i = 1,2,...,n
*
*     (       x(i)  z(i) ) :=
*     ( conjg(z(i)) y(i) )
*
*       (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )
*       ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of plane rotations to be applied.
*
*  X       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)
*          The vector x; the elements of x are assumed to be real.
*
*  Y       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)
*          The vector y; the elements of y are assumed to be real.
*
*  Z       (input/output) COMPLEX array, dimension (1+(N-1)*INCX)
*          The vector z.
*
*  INCX    (input) INTEGER
*          The increment between elements of X, Y and Z. INCX > 0.
*
*  C       (input) REAL array, dimension (1+(N-1)*INCC)
*          The cosines of the plane rotations.
*
*  S       (input) COMPLEX array, dimension (1+(N-1)*INCC)
*          The sines of the plane rotations.
*
*  INCC    (input) INTEGER
*          The increment between elements of C and S. INCC > 0.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IX
      REAL               CI, SII, SIR, T1I, T1R, T5, T6, XI, YI, ZII,
     $                   ZIR
      COMPLEX            SI, T2, T3, T4, ZI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, CMPLX, CONJG, REAL
*     ..
*     .. Executable Statements ..
*
      IX = 1
      IC = 1
      DO 10 I = 1, N
         XI = REAL( X( IX ) )
         YI = REAL( Y( IX ) )
         ZI = Z( IX )
         ZIR = REAL( ZI )
         ZII = AIMAG( ZI )
         CI = C( IC )
         SI = S( IC )
         SIR = REAL( SI )
         SII = AIMAG( SI )
         T1R = SIR*ZIR - SII*ZII
         T1I = SIR*ZII + SII*ZIR
         T2 = CI*ZI
         T3 = T2 - CONJG( SI )*XI
         T4 = CONJG( T2 ) + SI*YI
         T5 = CI*XI + T1R
         T6 = CI*YI - T1R
         X( IX ) = CI*T5 + ( SIR*REAL( T4 )+SII*AIMAG( T4 ) )
         Y( IX ) = CI*T6 - ( SIR*REAL( T3 )-SII*AIMAG( T3 ) )
         Z( IX ) = CI*T3 + CONJG( SI )*CMPLX( T6, T1I )
         IX = IX + INCX
         IC = IC + INCC
   10 CONTINUE
      RETURN
*
*     End of CLAR2V
*
      END
