! Exercise the INFO encoding of a failing merge (secular equation) in
! DSTEDS.  DSTEQR and DLAED4 below are deliberate mocks: DSTEQR returns
! an orthogonal leaf basis, so that the first merge has two
! nondeflated poles, and DLAED4 then fails.  DSTEDS, DLAES0, DLAES1
! and DLAED9_BR are compiled into this program by the build files, so
! that both mocks are reached even when LAPACK is a shared library.
! No eigenvalue accuracy or performance is measured by this program.
program test_dsteds_merge_failures
  implicit none
  integer :: n, i, mode, info, expected, nerr, ntest, lw, liw
  integer :: rows_done, block_start, leaf1, leaf2, nsolve
  integer :: iq(1), isel(1)
  common /dsteds_merge_test/ rows_done, block_start, leaf1, leaf2, &
                             nsolve
  double precision :: wq(1), dummy(1)
  double precision, allocatable :: d(:),e(:),w(:)
  integer, allocatable :: iw(:)
  external :: dsteds,dlasum
  n=50000
  allocate(d(n),e(n))
  isel=0
  dummy=0d0
  call dsteds(n,d,e,0,isel,dummy,n,wq,-1,iq,-1,info)
  if(info/=0) then
    deallocate(d,e)
    stop 1
  endif
  lw=int(wq(1))
  liw=iq(1)
  allocate(w(lw),iw(liw))
  nerr=0
  ntest=0
  do mode=1,2
    rows_done=0
    leaf1=0
    leaf2=0
    nsolve=0
    d=2d0
    e=0.1d0
    ! Mode 1: the first merge of the matrix fails.
    ! Mode 2: rows 1:45000 split into single-leaf blocks, so the first
    ! merge, rows 45001 onward, fails where its global row/column
    ! encoding does not fit INT32.
    block_start=1
    if(mode==2) then
      block_start=45001
      do i=25,45000,25
        e(i)=0d0
      enddo
    endif
    call dsteds(n,d,e,0,isel,dummy,n,w,lw,iw,liw,info)
    ntest=ntest+1
    if(nsolve/=1 .or. leaf1==0 .or. leaf2==0) then
      nerr=nerr+1
      print *, 'Mock did not reach the requested merge: ',mode,nsolve
    endif
    if(block_start>(huge(0)-(block_start+leaf1+leaf2-1))/(n+1)) then
      expected=1
    else
      expected=block_start*(n+1)+block_start+leaf1+leaf2-1
    endif
    ntest=ntest+1
    if(info/=expected) then
      nerr=nerr+1
      print *, 'Wrong failure encoding: ',mode,info,expected
    endif
  enddo
  call dlasum('DST',6,nerr,ntest)
  deallocate(d,e,w,iw)
  if(nerr/=0) stop 1
end program test_dsteds_merge_failures

subroutine dsteqr(compz,n,d,e,z,ldz,work,info)
  implicit none
  character :: compz
  integer :: n,ldz,info,i,j
  double precision :: d(*),e(*),z(ldz,*),work(*)
  integer :: rows_done, block_start, leaf1, leaf2, nsolve
  common /dsteds_merge_test/ rows_done, block_start, leaf1, leaf2, &
                             nsolve
  ! Record the sizes of the first two leaves of the failing block.
  if(rows_done>=block_start-1) then
    if(leaf1==0) then
      leaf1=n
    else if(leaf2==0) then
      leaf2=n
    endif
  endif
  rows_done=rows_done+n
  info=0
  ! Sorted poles and an identity basis: the rank-one update of the
  ! first merge then couples only one pole of each child.
  do i=1,n
    d(i)=dble(i)/dble(n+1)
    do j=1,n
      z(i,j)=0d0
    enddo
    z(i,i)=1d0
  enddo
end subroutine dsteqr

subroutine dlaed4(n,i,d,z,delta,rho,dlam,info)
  implicit none
  integer :: n,i,info
  double precision :: d(*),z(*),delta(*),rho,dlam
  integer :: rows_done, block_start, leaf1, leaf2, nsolve
  common /dsteds_merge_test/ rows_done, block_start, leaf1, leaf2, &
                             nsolve
  nsolve=nsolve+1
  info=1
end subroutine dlaed4
