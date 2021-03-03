!> \brief \b ZROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZROT( N, X, INCX, Y, INCY, C, S )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       REAL               C
!       COMPLEX            S
!       ..
!       .. Array Arguments ..
!       COMPLEX            X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZROT applies a plane rotation to vectors x and y:
!>
!>    [  c         s ] [ x(1)  x(2)  ...  x(n) ]
!>    [ -conjg(s)  c ] [ y(1)  y(2)  ...  y(n) ]
!>
!> where c is real, s is complex, and c**2 + conjg(s)*s = 1.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements in the vectors X and Y.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX(dp) array, dimension (N)
!>          On input, the vector X.
!>          On output, X is overwritten with C*X + S*Y.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of X.  INCX <> 0.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX(dp) array, dimension (N)
!>          On input, the vector Y.
!>          On output, Y is overwritten with -CONJG(S)*X + C*Y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          The increment between successive values of Y.  INCX <> 0.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL(dp)
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is COMPLEX(dp)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!> \endverbatim
!
!  =====================================================================
subroutine ZROT( n, x, incx, y, incy, c, s )
   use LA_CONSTANTS, only: wp
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
