subroutine ZLARTG( f, g, c, s, r )
   use LA_CONSTANTS, only: wp, zero, one, two, czero, rtmin, rtmax, &
                             safmin, safmax
!
!  LAPACK auxiliary routine
!  E. Anderson
!  August 4, 2016
!
!  .. Scalar Arguments ..
   real(wp)           c
   complex(wp)        f, g, r, s
!  ..
!
!  Purpose
!  =======
!
!  ZLARTG generates a plane rotation so that
!
!     [  C         S  ] . [ F ]  =  [ R ]
!     [ -conjg(S)  C  ]   [ G ]     [ 0 ]
!
!  where C is real and C**2 + |S|**2 = 1.
!
!  The mathematical formulas used for C and S are
!
!     sgn(x) = {  x / |x|,   x != 0
!              {  1,         x = 0
!
!     R = sgn(F) * sqrt(|F|**2 + |G|**2)
!
!     C = |F| / sqrt(|F|**2 + |G|**2)
!
!     S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2)
!
!  When F and G are real, the formulas simplify to C = F/R and
!  S = G/R, and the returned values of C, S, and R should be
!  identical to those returned by DLARTG.
!
!  The algorithm used to compute these quantities incorporates scaling
!  to avoid overflow or underflow in computing the square root of the
!  sum of squares.
!
!  Arguments
!  =========
!
!  F       (input) COMPLEX
!          The first component of vector to be rotated.
!
!  G       (input) COMPLEX
!          The second component of vector to be rotated.
!
!  C       (output) REAL
!          The cosine of the rotation.
!
!  S       (output) COMPLEX
!          The sine of the rotation.
!
!  R       (output) COMPLEX
!          The nonzero component of the rotated vector.
!
!  =====================================================================
!
!  .. Local Scalars ..
   real(wp) :: d, f1, f2, g1, g2, h2, p, u, uu, v, vv, w
   complex(wp) :: fs, gs, t
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, aimag, conjg, max, min, real, sqrt
!  ..
!  .. Statement Functions ..
   real(wp) :: ABSSQ
!  ..
!  .. Statement Function definitions ..
   ABSSQ( t ) = real( t )**2 + aimag( t )**2
!  ..
!  .. Executable Statements ..
!
   if( g == czero ) then
      c = one
      s = czero
      r = f
   else if( f == czero ) then
      c = zero
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         g2 = ABSSQ( g )
         d = sqrt( g2 )
         s = conjg( g ) / d
         r = d
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, g1 ) )
         uu = one / u
         gs = g*uu
         g2 = ABSSQ( gs )
         d = sqrt( g2 )
         s = conjg( gs ) / d
         r = d*u
      end if
   else
      f1 = max( abs(real(f)), abs(aimag(f)) )
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( f1 > rtmin .and. f1 < rtmax .and. &
          g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         f2 = ABSSQ( f )
         g2 = ABSSQ( g )
         h2 = f2 + g2
         if( f2 > rtmin .and. h2 < rtmax ) then
            d = sqrt( f2*h2 )
         else
            d = sqrt( f2 )*sqrt( h2 )
         end if
         p = 1 / d
         c = f2*p
         s = conjg( g )*( f*p )
         r = f*( h2*p )
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, f1, g1 ) )
         uu = one / u
         gs = g*uu
         g2 = ABSSQ( gs )
         if( f1*uu < rtmin ) then
!
!           f is not well-scaled when scaled by g1.
!           Use a different scaling for f.
!
            v = min( safmax, max( safmin, f1 ) )
            vv = one / v
            w = v * uu
            fs = f*vv
            f2 = ABSSQ( fs )
            h2 = f2*w**2 + g2
         else
!
!           Otherwise use the same scaling for f and g.
!
            w = one
            fs = f*uu
            f2 = ABSSQ( fs )
            h2 = f2 + g2
         end if
         if( f2 > rtmin .and. h2 < rtmax ) then
            d = sqrt( f2*h2 )
         else
            d = sqrt( f2 )*sqrt( h2 )
         end if
         p = 1 / d
         c = ( f2*p )*w
         s = conjg( gs )*( fs*p )
         r = ( fs*( h2*p ) )*u
      end if
   end if
   return
end subroutine
