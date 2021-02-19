subroutine SLARTG( f, g, c, s, r )
   use LA_CONSTANTS32, only: wp, zero, half, one, rtmin, rtmax, safmin, safmax
!
!  LAPACK auxiliary routine
!  E. Anderson
!  July 30, 2016
!
!  .. Scalar Arguments ..
   real(wp) :: c, f, g, r, s
!  ..
!
!  Purpose
!  =======
!
!  SLARTG generates a plane rotation so that
!
!     [  C  S  ]  .  [ F ]  =  [ R ]
!     [ -S  C  ]     [ G ]     [ 0 ]
!
!  where C**2 + S**2 = 1.
!
!  The mathematical formulas used for C and S are
!     R = sign(F) * sqrt(F**2 + G**2)
!     C = F / R
!     S = G / R
!  Hence C >= 0.  The algorithm used to compute these quantities
!  incorporates scaling to avoid overflow or underflow in computing the
!  square root of the sum of squares.
!
!  This version is discontinuous in R at F = 0 but it returns the same
!  C and S as CLARTG for complex inputs (F,0) and (G,0).
!
!  Arguments
!  =========
!
!  F       (input) REAL
!          The first component of vector to be rotated.
!
!  G       (input) REAL
!          The second component of vector to be rotated.
!
!  C       (output) REAL
!          The cosine of the rotation.
!
!  S       (output) REAL
!          The sine of the rotation.
!
!  R       (output) REAL
!          The nonzero component of the rotated vector.
!
!  =====================================================================
!
!  .. Local Scalars ..
   real(wp) :: d, f1, fs, g1, gs, p, u, uu
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, sign, sqrt
!  ..
!  .. Executable Statements ..
!
   f1 = abs( f )
   g1 = abs( g )
   if( g == zero ) then
      c = one
      s = zero
      r = f
   else if( f == zero ) then
      c = zero
      s = sign( one, g )
      r = g1
   else if( f1 > rtmin .and. f1 < rtmax .and. &
            g1 > rtmin .and. g1 < rtmax ) then
      d = sqrt( f*f + g*g )
      p = one / d
      c = f1*p
      s = g*sign( p, f )
      r = sign( d, f )
   else
      u = min( safmax, max( safmin, f1, g1 ) )
      uu = one / u
      fs = f*uu
      gs = g*uu
      d = sqrt( fs*fs + gs*gs )
      p = one / d
      c = abs( fs )*p
      s = gs*sign( p, f )
      r = sign( d, f )*u
   end if
   return
end subroutine
