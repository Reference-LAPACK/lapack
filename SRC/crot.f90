subroutine CROT( n, x, incx, y, incy, c, s )
   use LA_CONSTANTS32, only: wp
!
!  Updated Level 1 BLAS
!  E. Anderson
!  August 11, 2016
!
!  .. Scalar Arguments ..
   integer :: incx, incy, n
   real(wp) :: c
   complex(wp) :: s
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*), y(*)
!  ..
!
!  Purpose
!  =======
!
!  CROT applies a plane rotation to vectors x and y:
!
!     [  c         s ] [ x(1)  x(2)  ...  x(n) ]
!     [ -conjg(s)  c ] [ y(1)  y(2)  ...  y(n) ]
!
!  where c is real, s is complex, and c**2 + conjg(s)*s = 1.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements of the vectors x and y.
!
!  X       (input/output) COMPLEX array, dimension (1+(N-1)*abs(INCX))
!          On entry, the n-element vector x.
!          On exit, the vector sum c*x + s*y.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector x.
!          If INCX >= 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!
!  Y       (input/output) COMPLEX array, dimension (1+(N-1)*abs(INCY))
!          On entry, the n-element vector y.
!          On exit, the vector sum -conjg(s)*x + c*y.
!
!  INCY    (input) INTEGER
!          The increment between successive values of the vector y.
!          If INCY >= 0, Y(1+(i-1)*INCY) = y(i) for 1 <= i <= n
!          If INCY < 0, Y(1-(n-i)*INCY) = y(i) for 1 <= i <= n
!
!  C       (input) REAL
!          The scalar c.
!
!  S       (input) COMPLEX
!          The scalar s.
!
! =====================================================================
!
!  .. Local Scalars ..
   integer :: i, ix, iy
   complex(wp) :: stmp
!  ..
!
!  Quick return if possible
!
   if( n <= 0 ) return
!
   if( incx == 1 .and. incy == 1 ) then
!
!     Special case for incx = incy = 1
!
      do i = 1, n
         stmp = c*x(i) + s*y(i)
         y(i) = c*y(i) - conjg(s)*x(i)
         x(i) = stmp
      end do
   else
!
!     General increment case
!
      ix = 1
      iy = 1
      if( incx < 0 ) ix = 1 - (n-1)*incx
      if( incy < 0 ) iy = 1 - (n-1)*incy
      do i = 1, n
         stmp = c*x(ix) + s*y(iy)
         y(iy) = c*y(iy) - conjg(s)*x(ix)
         x(ix) = stmp
         ix = ix + incx
         iy = iy + incy
      end do
   end if
   return
end subroutine
