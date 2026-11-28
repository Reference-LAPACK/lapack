! Check inverse iteration on a tridiagonal with a tightly clustered spectrum.
program test_zstein_cluster
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: n = 13, wp = kind(1.0d0)
  integer :: iblock(n), isplit(n), ifail(n), iwork(3*n), m, nsplit, info, j
  real(wp) :: d(n), e(n-1), w(n), work(5*n), result(2), rwork(n)
  complex(wp) :: z(n,n), testwork(n*n)
  external :: dstebz, zstein, zstt21

  d = 0
  do j = 1, n-1, 4
    e(j:j+3) = [1e-10_wp, 1e-10_wp, 1.0_wp, 1e-120_wp]
  end do
  call dstebz('A','B',n,0.0_wp,0.0_wp,0,0,-1.0_wp,d,e,m,nsplit,w,iblock,isplit,work,iwork,info)
  if (info /= 0 .or. m /= n .or. nsplit /= 1) stop 1
  call zstein(n,d,e,m,w,iblock,isplit,z,n,work,iwork,ifail,info)
  if (info /= 0 .or. any(ifail /= 0)) stop 2
  if (.not.all(ieee_is_finite(abs(z))) .or. .not.all(ieee_is_finite(w))) stop 3
  call zstt21(n,0,d,e,w,e,z,n,testwork,rwork,result)
  if (.not.all(result <= 32)) then
    print *, 'ZSTEIN reconstruction and orthogonality ratios:', result
    stop 4
  end if
end program
