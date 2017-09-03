    program spmmtest

    implicit none

    doubleprecision A(6),B(4,3),C(3,3)
    integer M,N,LDB,LDC
    doubleprecision ALPHA,BETA

    M=3
    N=3
    LDB=4
    LDC=3
    ALPHA=1.0
    BETA=0.0

    A(1)=1
    A(2)=2
    A(3)=3
    A(4)=4
    A(5)=5
    A(6)=6

    B(1,1)=5
    B(2,1)=4
    B(3,1)=3

    !  [ 1 2 3 ]   [ 5 0 0 ]   [ 22 0 0 ]
    !  [ 2 4 5 ] X [ 4 0 0 ] = [ 38 0 0 ]
    !  [ 3 5 6 ]   [ 3 0 0 ]   [ 53 0 0 ]
    call DSPMM('R','L','N',M,N,A,ALPHA,B,LDB,BETA,C,LDC)
    print *, C
    print *
    !  [ 5 0 0 ]   [ 1 2 3 ]   [ 5 10 15]
    !  [ 4 0 0 ] X [ 2 4 5 ] = [ 4  8 12]
    !  [ 3 0 0 ]   [ 3 5 6 ]   [ 3  6  9]
    call DSPMM('L','L','N',M,N,A,ALPHA,B,LDB,BETA,C,LDC)
    print *, C
    print *
    !  [ 1 2 4 ]   [ 5 0 0 ]   [ 25 0 0 ]
    !  [ 2 3 5 ] X [ 4 0 0 ] = [ 37 0 0 ]
    !  [ 4 5 6 ]   [ 3 0 0 ]   [ 58 0 0 ]
    call DSPMM('R','U','N',M,N,A,ALPHA,B,LDB,BETA,C,LDC)
    print *, C
    print *
    !  [ 5 0 0 ]   [ 1 2 4 ]   [ 5 10 20]
    !  [ 4 0 0 ] X [ 2 3 5 ] = [ 4  8 16]
    !  [ 3 0 0 ]   [ 4 5 6 ]   [ 3  6 12]
    call DSPMM('L','U','N',M,N,A,ALPHA,B,LDB,BETA,C,LDC)
    print *, C
    print *
    !  [ 1 2 4 ]   [ 5 4 3 ]   [ 5  4  3]
    !  [ 2 3 5 ] X [ 0 0 0 ] = [10  8  6]
    !  [ 4 5 6 ]   [ 0 0 0 ]   [20 16 12]
    call DSPMM('R','U','T',M,N,A,ALPHA,B,LDB,BETA,C,LDC)
    print *, C
    print *
    !  [ 5 4 3 ]   [ 1 2 4 ]   [ 25 37 58 ]
    !  [ 0 0 0 ] X [ 2 3 5 ] = [  0  0  0 ]
    !  [ 0 0 0 ]   [ 4 5 6 ]   [  0  0  0 ]
    call DSPMM('L','U','T',M,N,A,ALPHA,B,LDB,BETA,C,LDC)
    print *, C
    print *

    end

