! Boundary and query tests for DSTEDS. Numerical checks use DCHKSTEDS,
! DSTT21 and DSTT22 with the threshold read from steds.in.
program test_dsteds_interfaces
  implicit none
  integer :: sizes(8), leaves(2), n, nn, js, jl, kind, sc, i, j, m
  integer :: info, nerr, ntest, iq(1), isq(1), big, levels, sz, expected
  integer :: infot, nout
  logical :: lerr, ok
  character(32) :: srnamt
  common /infoc/ infot, nout, ok, lerr
  common /srnamc/ srnamt
  double precision :: thresh, result(6), scales(3), fac, wq(1)
  double precision, allocatable :: d(:), e(:), ref(:), ec(:), rw(:)
  double precision, external :: dlamch
  external :: dchksteds, dsteqr, dsteds, xlaenv, dlasum, chkxer
  external :: set_lapack_xerbla, xer_replace
  call set_lapack_xerbla(xer_replace)
  read(*,*) thresh
  if(.not.(thresh>0d0)) stop 1
  infot=0
  nout=6
  ok=.true.
  lerr=.false.
  srnamt='DSTEDS'
  nerr=0
  ntest=0
  sizes=[0,1,2,3,25,26,51,129]
  leaves=[2,25]
  scales=[1d0,1d-280,1d280]
  call xlaenv(1,1)
  do jl=1,size(leaves)
    call xlaenv(9,leaves(jl))
    do js=1,size(sizes)
      n=sizes(js)
      nn=max(1,n)
      allocate(d(nn), e(nn), ref(nn), ec(nn), rw(max(1,2*n)))
      do kind=1,8
        do sc=1,size(scales)
          d=0d0
          e=0d0
          do i=1,n
            select case(kind)
            case(1)                   ! Zero matrix.
              d(i)=0d0
            case(2)                   ! Identity and repeated eigenvalues.
              d(i)=1d0
            case(3)                   ! Unsorted diagonal, repeated entries.
              d(i)=dble(mod(7*i,5)-2)
            case(4)                   ! Negative offdiagonal Toeplitz.
              d(i)=2d0
              e(i)=-1d0
            case(5)                   ! Split blocks, interleaving spectra.
              d(i)=dble(mod(17*i,11))
              if(mod(i,11)/=0) e(i)=0.1d0
            case(6)                   ! Tight clusters with signed couplings.
              d(i)=1d0+dble(mod(i,5))*1d-12
              e(i)=(-1d0)**i*1d-14
            case(7)                   ! Graded spectrum.
              d(i)=10d0**(-12d0*dble(i-1)/dble(max(1,n-1)))
              e(i)=0.01d0*d(i)
            case(8)                   ! Zero diagonal and signed couplings.
              e(i)=(-1d0)**i
            end select
          enddo
          fac=scales(sc)
          d=d*fac
          e=e*fac
          ref=d
          ec=e
          call dsteqr('N',n,ref,ec,rw,nn,rw,info)
          call check(info==0,'DSTEQR reference')
          if(info/=0) cycle
          call dchksteds(n,d,e,ref,result,info)
          call check(info==0,'DSTEDS computation')
          do j=1,6
            ntest=ntest+1
            if(.not.(result(j)>=0d0 .and. result(j)<thresh)) then
              nerr=nerr+1
              write(*,'(A,5I6,A,ES12.4)') 'DSTEDS n/type/scale/leaf/test:', &
                n,kind,sc,leaves(jl),j,' ratio=',result(j)
            endif
          enddo
        enddo
      enddo
      deallocate(d,e,ref,ec,rw)
    enddo
  enddo
  ! N=0/1 quick returns need one element in each workspace. The query
  ! must not inspect ISEL or data, even if the arrays have only one entry.
  do n=0,1
    isq=0
    wq=-7d0
    call dsteds(n, wq, wq, n, isq, wq, 1, wq, -1, iq, -1, info)
    call check(info==0 .and. wq(1)==1d0 .and. iq(1)==1, &
               'N=0/1 minimal query')
  enddo
  ! These are dimension-only queries: no large matrix is allocated.
  ! They also run with -ftrapv in the bounds-checked test build.
  do j=1,3
    big=huge(0)
    big=big/4
    m=0
    if(j==2) m=big
    if(j==3) big=huge(0)
    infot=9
    call dsteds(big, wq, wq, m, isq, wq, big, wq, -1, iq, -1, info)
    call chkxer('DSTEDS',infot,nout,lerr,ok)
    call check(info==-9,'Unrepresentable workspace query')
  enddo
  infot=0
  ! Above 2**53 a double cannot represent every INTEGER in the _64
  ! interface.  The returned floating size must round upward.  This
  ! query allocates only scalar arrays, not a matrix of this dimension.
  if(bit_size(0)>32) then
    big=134217731
    sz=big
    levels=0
    do while(sz>25)
      sz=sz/2+mod(sz,2)
      levels=levels+1
    enddo
    expected=big*(25+7*levels+19+64+2*big)
    call dsteds(big, wq, wq, big, isq, wq, big, wq, -1, iq, -1, info)
    call check(info==0 .and. int(wq(1))>=expected, &
               '64-bit workspace rounds upward')
  endif
  call check(ok .and. .not.lerr,'No unexpected XERBLA')
  call dlasum('DST',6,nerr,ntest)
  if(nerr/=0) stop 1
contains
  subroutine check(valid,label)
    logical, intent(in) :: valid
    character(*), intent(in) :: label
    ntest=ntest+1
    if(.not.valid) then
      nerr=nerr+1
      print *, 'DSTEDS interface failure: ',label
    endif
  end subroutine check
end program test_dsteds_interfaces
