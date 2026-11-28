! Clustered eigenvalues in a weakly connected tridiagonal can trigger
! overflow during inverse iteration (Reference-LAPACK/lapack#1474).
program test_dstein_cluster
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer :: failures
  real(dp), external :: dlamch
  external :: dstebz, dstein, dstevx, dsyevx
  failures = 0
  call clustered(3, 1.0d-120, 1.0d0, 4)
  call clustered(11, 1.0d-30, 1.0d0, 4)
  call clustered(11, 1.0d-30, 1.0d-80, 1)
  ! The near-zero shift and the large second diagonal stress forward substitution.
  call exercise([0.0d0,1.0d75], [1.0d-145], 2*dlamch('Safe minimum'), 1, 1, [0.0d0,1.0d75])
  call coarse()
  if (failures /= 0) stop 1

contains

  subroutine clustered(k, link, factor, methods)
    integer, intent(in) :: k, methods
    real(dp), intent(in) :: link, factor
    real(dp) :: d(4*k+1), e(4*k)
    integer :: j
    d = 0
    do j = 1, 4*k, 4
      e(j:j+3) = [1.0d-10, 1.0d-10, 1.0d0, link]
    end do
    call exercise(d, factor*e, -1.0d0, methods, size(d))
  end subroutine

  subroutine coarse()
    ! The five eigenvalues of this Toeplitz matrix are known analytically.
    ! Positive ABSTOL intentionally supplies approximate shifts to DSTEIN.
    integer, parameter :: n = 5
    real(dp) :: d(n), e(n-1), exact(n)
    integer :: j
    d = 2
    e = 0.25d0
    do j = 1, n
      exact(j) = 2.0d0 - 0.5d0 * cos(real(j,dp) * acos(-1.0d0) / (n+1))
    end do
    call exercise(d, e, 1.0d-2, 4, n, exact)
  end subroutine

  subroutine exercise(d, e, abstol, methods, last, exact)
    integer, intent(in) :: methods, last
    real(dp), intent(in) :: d(:), e(:), abstol
    real(dp), intent(in), optional :: exact(:)
    integer :: n, m, nsplit, info, j, method, first
    real(dp) :: a(size(d),size(d)), dc(size(d)), ec(size(e)), w(size(d)), z(size(d),size(d))
    real(dp) :: work(8*size(d)), anorm, resid, orth, r(size(d)), g(size(d)), tol, otol, gap
    integer :: iblock(size(d)), isplit(size(d)), iwork(5*size(d)), ifail(size(d))
    real(dp) :: ac(size(d),size(d)), norm_error
    character(1) :: uplo, range
    character(8), parameter :: names(4) = ['DSTEIN  ', 'DSTEVX  ', 'DSYEVX U', 'DSYEVX L']

    n = size(d)
    a = 0
    do j = 1, n
      a(j,j) = d(j)
      if (j < n) then
        a(j,j+1) = e(j)
        a(j+1,j) = e(j)
      end if
    end do
    anorm = maxval(sum(abs(a),dim=1))
    tol = 32*n*epsilon(1.0d0)
    otol = tol
    if (present(exact)) then
      ! A vector belonging to the exact eigenvalue can have residual ABSTOL
      ! against an approximate shift. The separated eigenvalues bound the
      ! loss of orthogonality attributable to that uncertainty.
      gap = minval(exact(2:) - exact(:n-1))
      tol = max(tol, abstol / anorm)
      otol = max(otol, 2*abstol / (gap-2*abstol))
    end if
    do method = 1, methods
      ifail = -1
      if (method == 1) then
        first = 1
        range = 'A'
        if (last < n) range = 'I'
        call dstebz(range,'B',n,0.0d0,0.0d0,first,last,abstol,d,e,m,nsplit,w,iblock,isplit,work,iwork,info)
        if (info /= 0 .or. m /= last .or. nsplit /= 1) then
          print *, 'DSTEBZ N, M, NSPLIT, INFO:', n, m, nsplit, info
          stop 2
        end if
        call dstein(n,d,e,m,w,iblock,isplit,z,n,work,iwork,ifail,info)
      else
        ! A subset forces inverse iteration instead of the full-spectrum QR shortcut.
        first = 2
        if (method == 2) then
          dc = d
          ec = e
          call dstevx('V','I',n,dc,ec,0.0d0,0.0d0,first,last,abstol,m,w,z,n,work,iwork,ifail,info)
        else
          uplo = 'U'
          if (method == 4) uplo = 'L'
          ac = a
          call dsyevx('V','I',uplo,n,ac,n,0.0d0,0.0d0,first,last,abstol,m,w,z,n,work,size(work),iwork,ifail,info)
        end if
      end if
      if (info /= 0 .or. m /= last-first+1) then
        print *, names(method), ' N, M, INFO:', n, m, info
        failures = failures + 1
        cycle
      end if
      if (any(ifail(1:m) /= 0) .or. .not. all(ieee_is_finite(z(:,1:m))) &
          .or. .not. all(ieee_is_finite(w(1:m)))) then
        print *, names(method), ' invalid output, N=', n
        failures = failures + 1
        cycle
      end if
      resid = 0
      orth = 0
      norm_error = 0
      do j = 1, m
        r = matmul(a,z(:,j)) - w(j)*z(:,j)
        g(1:m) = matmul(transpose(z(:,1:m)),z(:,j))
        norm_error = max(norm_error, abs(g(j)-1))
        g(j) = 0
        resid = max(resid, norm2(r) / anorm)
        orth = max(orth, maxval(abs(g(1:m))))
      end do
      if (.not.(resid <= tol .and. orth <= otol .and. norm_error <= 32*n*epsilon(1.0d0))) then
        print *, names(method), ' N, residual/norm, orthogonality, norm error:', n, resid, orth, norm_error
        failures = failures + 1
      end if
      if (present(exact)) then
        if (any(abs(w(1:m)-exact(first:last)) > abstol)) then
          print *, names(method), ' eigenvalue outside ABSTOL, N=', n
          failures = failures + 1
        end if
      end if
    end do
  end subroutine
end program
