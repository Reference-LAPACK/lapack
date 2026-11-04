!> \brief Tests DSTEDS using the existing DSTT21/DSTT22 ratios.
!>
!> RESULT(1): eigenvalues, maximum over all selections and M=0.
!> RESULT(2:3): full decomposition (DSTT21).
!> RESULT(4:5): lowest, highest, middle and nonconsecutive subsets
!>              (DSTT22).
!> RESULT(6): workspace-query consistency and preservation, sorting,
!>             and that DSTEDS writes no element of Z outside
!>             Z(1:N,1:M) (0 or 1/ULP).
!> Every selection is computed with LDZ = N and with LDZ = N+3.
!> N, AD, AE describe the original tridiagonal; REF contains its sorted
!> eigenvalues from the caller's reference solver.  Inputs are preserved.
!> Workspace is allocated from a query, independent of DCHKST workspace.
subroutine dchksteds(n, ad, ae, ref, result, info)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: ad(*), ae(*), ref(*)
  double precision, intent(out) :: result(6)
  integer, intent(out) :: info
  integer :: nn, lw, liw, q, mode, m, i, j, ierr, lwc, liwc, ldz, jl
  double precision :: ulp, unfl, wq(1), pair(2), den, diff
  double precision, allocatable :: d(:), e(:), z(:), work(:), &
                                   testwork(:,:), selected(:)
  integer, allocatable :: isel(:), iw(:)
  integer :: iq(1)
  double precision, external :: dlamch
  external :: dsteds, dstt21, dstt22
  result = 0d0
  info = 0
  ulp = dlamch('Precision')
  unfl = dlamch('Safe minimum')
  nn = max(1,n)
  allocate(d(nn), e(nn), z((nn+3)*nn), isel(nn), selected(nn), &
           testwork(nn,nn+1), stat=ierr)
  if(ierr /= 0) then
    info = ierr
    result = 1d0/ulp
    return
  endif
  do i=1,nn
    isel(i) = i
  enddo
  ! Independent real, integer and joint queries must agree and preserve
  ! input/output arrays.  Queries may use only one-element work arrays.
  d = -7d0
  e = -7d0
  z = -7d0
  lw = 0
  liw = 0
  do q=1,3
    i = 1
    j = 1
    if(q /= 2) i = -1
    if(q /= 1) j = -1
    call dsteds(n, d, e, n, isel, z, nn, wq, i, iq, j, info)
    if(info /= 0) goto 900
    if(q == 1) then
      lw = int(wq(1))
      liw = iq(1)
    endif
    if(wq(1) /= dble(lw) .or. iq(1) /= liw .or. &
       lw < 1 .or. liw < 1 .or. any(d /= -7d0) .or. &
       any(e /= -7d0) .or. any(z /= -7d0)) result(6)=1d0/ulp
    do i=1,nn
      if(isel(i) /= i) result(6)=1d0/ulp
    enddo
  enddo
  allocate(work(lw), iw(liw), stat=ierr)
  if(ierr /= 0) then
    info = ierr
    goto 900
  endif
  ! Full spectrum and M=0, then four types of selected column sets,
  ! each with a tight and a padded leading dimension.
  do jl=0,1
  ldz=nn+3*jl
  do mode=0,5
    select case(mode)
    case(0)
      m = n
      do i=1,m
        isel(i)=i
      enddo
    case(1)
      m = 0
    case default
      m = (n+2)/3
      do i=1,m
        select case(mode)
        case(2)
          isel(i)=i
        case(3)
          isel(i)=n-m+i
        case(4)
          isel(i)=(n-m)/2+i
        case(5)
          if(m==1) then
            isel(i)=1
          else
            isel(i)=1+((i-1)*(n-1))/(m-1)
          endif
        end select
      enddo
    end select
    d(1:n)=ad(1:n)
    if(n>1) e(1:n-1)=ae(1:n-1)
    z=-7d0
    call dsteds(n, d, e, m, isel, z, ldz, wq, -1, iq, -1, info)
    if(info/=0) goto 900
    lwc=int(wq(1))
    liwc=iq(1)
    work=-7d0
    iw=-7
    call dsteds(n, d, e, m, isel, z, ldz, work, lwc, iw, liwc, info)
    if(info/=0) goto 900
    if(work(1)/=dble(lwc) .or. iw(1)/=liwc) &
      result(6)=1d0/ulp
    if(any(work(lwc+1:)/=-7d0) .or. any(iw(liwc+1:)/=-7)) &
      result(6)=1d0/ulp
    do j=1,m
      if(any(z((j-1)*ldz+n+1:j*ldz)/=-7d0)) result(6)=1d0/ulp
    enddo
    if(any(z(m*ldz+1:)/=-7d0)) result(6)=1d0/ulp
    do i=2,n
      if(d(i)<d(i-1)) result(6)=1d0/ulp
    enddo
    den=0d0
    diff=0d0
    do i=1,n
      den=max(den,abs(ref(i)),abs(d(i)))
      diff=max(diff,abs(ref(i)-d(i)))
    enddo
    result(1)=max(result(1),diff/max(unfl,ulp*max(den,diff)))
    if(mode==0) then
      call dstt21(n,0,ad,ae,d,e,z,ldz,testwork,pair)
      result(2:3)=max(result(2:3),pair)
    else if(mode>=2) then
      do i=1,m
        selected(i)=d(isel(i))
      enddo
      call dstt22(n,m,0,ad,ae,selected,e,z,ldz,testwork,max(1,m),pair)
      result(4:5)=max(result(4:5),pair)
    endif
  enddo
  enddo
  return
900 continue
  result=1d0/ulp
end subroutine dchksteds
