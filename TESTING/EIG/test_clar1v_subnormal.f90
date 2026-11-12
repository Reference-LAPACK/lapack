! A subnormal nonzero progressive quotient must not lose the stored product.
! Exercise the public CSTEMR route, both TRYRAC settings, scale/sign
! neighbors, ordinary coupling and a split matrix; compare independent QR.
program test_clar1v_subnormal
  use, intrinsic :: iso_fortran_env, only: wp => real32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: n=4
  real(wp) :: d(n), e(n), w(n), work(18*n), d0(n), e0(n)
  real(wp) :: wr(n), er(n), zr(n,n), qrwork(2*n), a(n,n)
  real(wp) :: edge, factor, anorm, tol, residual, orthogonality, eigenvalues
  complex(wp) :: z(n,n), gram(n,n)
  real(wp), external :: slamch
  integer :: i, mode, case_id, m, info, isuppz(2*n), iwork(10*n), checks
  logical :: tryrac
  external :: cstemr, ssteqr

  tol=100*n*slamch('Precision')
  checks=0
  do case_id=1,8
    edge=1.073741824e-6_wp
    factor=1.0_wp
    select case(case_id)
    case(2)
      edge=scale(edge,-1)
    case(3)
      edge=scale(edge,1)
    case(4)
      edge=1.0e-3_wp
    case(5)
      factor=scale(factor,-16)
    case(6)
      factor=scale(factor,16)
    case(8)
      edge=0.0_wp
    end select
    d0=0.0_wp
    e0=[edge,edge,1.0_wp,0.0_wp]*factor
    if(case_id==7) e0(2)=-e0(2)
    a=0.0_wp
    do i=1,n-1
      a(i,i+1)=e0(i)
      a(i+1,i)=e0(i)
    end do
    anorm=maxval(sum(abs(a),dim=1))
    wr=d0; er=e0
    call ssteqr('I',n,wr,er,zr,n,qrwork,info)
    if(info/=0 .or. .not.all(ieee_is_finite(wr))) stop 1
    do mode=0,1
      d=d0; e=e0; z=0; w=0; m=0
      tryrac=mode==1
      call cstemr('V','A',n,d,e,0.0_wp,0.0_wp,0,0,m,w,z,n,n, &
                  isuppz,tryrac,work,size(work),iwork,size(iwork),info)
      if(info/=0 .or. m/=n) stop 1
      if(.not.all(ieee_is_finite(w)) .or. .not.(all(ieee_is_finite(real(z,wp))) .and. all(ieee_is_finite(aimag(z))))) stop 1
      gram=matmul(transpose(conjg(z)),z)
      do i=1,n
        gram(i,i)=gram(i,i)-1.0_wp
      end do
      orthogonality=maxval(abs(gram))
      residual=0.0_wp
      do i=1,n
        residual=max(residual,maxval(abs(matmul(a/anorm,z(:,i))-(w(i)/anorm)*z(:,i))))
      end do
      eigenvalues=maxval(abs(w/anorm-wr/anorm))
      print '(a,3i3,3es14.5)', 'CSTEMR case/TRYRAC/INFO, residual/Gram/eigenvalues: ', &
            case_id,mode,info,residual,orthogonality,eigenvalues
      if(.not.(residual<=tol .and. orthogonality<=tol .and. eigenvalues<=tol)) stop 1
      if(any(isuppz<1) .or. any(isuppz>n)) stop 1
      checks=checks+1
    end do
  end do
  if(checks/=16) stop 1
  print '(a)', 'CSTEMR subnormal regression: 16 checks passed'
end program
