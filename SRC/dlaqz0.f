*> \brief \b DLAQZ0
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAQZ0 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqz0.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqz0.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqz0.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAQZ0( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB,
*                          ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK,
*                          LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          WANTS, WANTQ, WANTZ
*       INTEGER            IHI, ILO, INFO, LDA, LDQ, LDB, LDZ, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * ), BETA( * ),
*      $                   A( LDA, * ), Q( LDQ, * ), B( LDB, * ),
*      $                   WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAQZ0 computes the eigenvalues of a real matrix pair (H,T),
*> where H is an upper Hessenberg matrix and T is upper triangular,
*> using the double-shift QZ method.
*> Matrix pairs of this type are produced by the reduction to
*> generalized upper Hessenberg form of a real matrix pair (A,B):
*>
*>    A = Q1*H*Z1**T,  B = Q1*T*Z1**T,
*>
*> as computed by DGGHRD.
*>
*> If JOB='S', then the Hessenberg-triangular pair (H,T) is
*> also reduced to generalized Schur form,
*>
*>    H = Q*S*Z**T,  T = Q*P*Z**T,
*>
*> where Q and Z are orthogonal matrices, P is an upper triangular
*> matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2
*> diagonal blocks.
*>
*> The 1-by-1 blocks correspond to real eigenvalues of the matrix pair
*> (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of
*> eigenvalues.
*>
*> Additionally, the 2-by-2 upper triangular diagonal blocks of P
*> corresponding to 2-by-2 blocks of S are reduced to positive diagonal
*> form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0,
*> P(j,j) > 0, and P(j+1,j+1) > 0.
*>
*> Optionally, the orthogonal matrix Q from the generalized Schur
*> factorization may be postmultiplied into an input matrix Q1, and the
*> orthogonal matrix Z may be postmultiplied into an input matrix Z1.
*> If Q1 and Z1 are the orthogonal matrices from DGGHRD that reduced
*> the matrix pair (A,B) to generalized upper Hessenberg form, then the
*> output matrices Q1*Q and Z1*Z are the orthogonal factors from the
*> generalized Schur factorization of (A,B):
*>
*>    A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T.
*>
*> To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently,
*> of (A,B)) are computed as a pair of values (alpha,beta), where alpha is
*> complex and beta real.
*> If beta is nonzero, lambda = alpha / beta is an eigenvalue of the
*> generalized nonsymmetric eigenvalue problem (GNEP)
*>    A*x = lambda*B*x
*> and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
*> alternate form of the GNEP
*>    mu*A*y = B*y.
*> Real eigenvalues can be read directly from the generalized Schur
*> form:
*>   alpha = S(i,i), beta = P(i,i).
*>
*> Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*>      Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*>      pp. 241--256.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTS
*> \verbatim
*>          WANTS is CHARACTER*1
*>          = 'E': Compute eigenvalues only;
*>          = 'S': Compute eigenvalues and the Schur form.
*> \endverbatim
*>
*> \param[in] WANTQ
*> \verbatim
*>          WANTQ is CHARACTER*1
*>          = 'N': Left Schur vectors (Q) are not computed;
*>          = 'I': Q is initialized to the unit matrix and the matrix Q
*>                 of left Schur vectors of (A,B) is returned;
*>          = 'V': Q must contain an orthogonal matrix Q1 on entry and
*>                 the product Q1*Q is returned.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is CHARACTER*1
*>          = 'N': Right Schur vectors (Z) are not computed;
*>          = 'I': Z is initialized to the unit matrix and the matrix Z
*>                 of right Schur vectors of (A,B) is returned;
*>          = 'V': Z must contain an orthogonal matrix Z1 on entry and
*>                 the product Z1*Z is returned.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, Q, and Z.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          ILO and IHI mark the rows and columns of A which are in
*>          Hessenberg form.  It is assumed that A is already upper
*>          triangular in rows and columns 1:ILO-1 and IHI+1:N.
*>          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          On entry, the N-by-N upper Hessenberg matrix A.
*>          On exit, if JOB = 'S', A contains the upper quasi-triangular
*>          matrix S from the generalized Schur factorization.
*>          If JOB = 'E', the diagonal blocks of A match those of S, but
*>          the rest of A is unspecified.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max( 1, N ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB, N)
*>          On entry, the N-by-N upper triangular matrix B.
*>          On exit, if JOB = 'S', B contains the upper triangular
*>          matrix P from the generalized Schur factorization;
*>          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S
*>          are reduced to positive diagonal form, i.e., if A(j+1,j) is
*>          non-zero, then B(j+1,j) = B(j,j+1) = 0, B(j,j) > 0, and
*>          B(j+1,j+1) > 0.
*>          If JOB = 'E', the diagonal blocks of B match those of P, but
*>          the rest of B is unspecified.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max( 1, N ).
*> \endverbatim
*>
*> \param[out] ALPHAR
*> \verbatim
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)
*>          The real parts of each scalar alpha defining an eigenvalue
*>          of GNEP.
*> \endverbatim
*>
*> \param[out] ALPHAI
*> \verbatim
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)
*>          The imaginary parts of each scalar alpha defining an
*>          eigenvalue of GNEP.
*>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
*>          positive, then the j-th and (j+1)-st eigenvalues are a
*>          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j).
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>          The scalars beta that define the eigenvalues of GNEP.
*>          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
*>          beta = BETA(j) represent the j-th eigenvalue of the matrix
*>          pair (A,B), in one of the forms lambda = alpha/beta or
*>          mu = beta/alpha.  Since either lambda or mu may overflow,
*>          they should not, in general, be computed.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
*>          On entry, if COMPQ = 'V', the orthogonal matrix Q1 used in
*>          the reduction of (A,B) to generalized Hessenberg form.
*>          On exit, if COMPQ = 'I', the orthogonal matrix of left Schur
*>          vectors of (A,B), and if COMPQ = 'V', the orthogonal matrix
*>          of left Schur vectors of (A,B).
*>          Not referenced if COMPQ = 'N'.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.  LDQ >= 1.
*>          If COMPQ='V' or 'I', then LDQ >= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*>          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in
*>          the reduction of (A,B) to generalized Hessenberg form.
*>          On exit, if COMPZ = 'I', the orthogonal matrix of
*>          right Schur vectors of (H,T), and if COMPZ = 'V', the
*>          orthogonal matrix of right Schur vectors of (A,B).
*>          Not referenced if COMPZ = 'N'.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1.
*>          If COMPZ='V' or 'I', then LDZ >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          = 1,...,N: the QZ iteration did not converge.  (A,B) is not
*>                     in Schur form, but ALPHAR(i), ALPHAI(i), and
*>                     BETA(i), i=INFO+1,...,N should be correct.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Thijs Steel
*
*> \date May 2020
*
*> \ingroup doubleGEcomputational
*>
*  =====================================================================
      subroutine dlaqz0(wantS,wantQ,wantZ,n,ilo,ihi,A,ldA,B,ldB,alphar,
     $   alphai,beta,Q,ldQ,Z,ldZ,work,lwork,info)
      implicit none

*     Arguments
      character,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ,lwork

      integer,intent(out) :: info

      double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $   Z(ldZ,*),alphar(*),alphai(*),beta(*),work(*)

*     Parameters
      double precision :: zero,one,half
      parameter(zero=0.0d0,one=1.0d0,half=0.5d0)

*     Local scalars
      double precision :: smlnum,ulp,eshift,safmin,safmax,c1,s1,temp
      integer :: istart,istop,iiter,maxit,istart2,k,ld,nshifts,nblock,
     $   nw,nmin,nibble,n_undeflated,n_deflated,ns,sweep_info,shiftpos,
     $   lworkreq,k2,istartm,istopm,iwants,iwantq,iwantz,norm_info,
     $   aed_info,nwr,nbr,nsr,itemp1,itemp2
      logical :: ilschur,ilq,ilz,lquery
      character :: jbcmpz*3

*     External Functions
      double precision,external :: dlamch
      logical,external :: lsame
      integer,external :: ilaenv

*
*     Decode wantS,wantQ,wantZ
*      
      if(lsame(wantS,'E')) then
         ilschur = .false.
         iwants = 1
      else if(lsame(wantS,'S')) then
         ilschur = .true.
         iwants = 2
      else
         iwants = 0
      end if

      if(lsame(wantQ,'N')) then
         ilq = .false.
         iwantq = 1
      else if(lsame(wantQ,'V')) then
         ilq = .true.
         iwantq = 2
      else if(lsame(wantQ,'I')) then
         ilq = .true.
         iwantq = 3
      else
         iwantq = 0
      end if

      if(lsame(wantZ,'N')) then
         ilz = .false.
         iwantz = 1
      else if(lsame(wantZ,'V')) then
         ilz = .true.
         iwantz = 2
      else if(lsame(wantZ,'I')) then
         ilz = .true.
         iwantz = 3
      else
         iwantz = 0
      end if
*
*     Check Argument Values
*
      info = 0
      if( iwants.EQ.0 ) then
         info =-1
      else if( iwantq.EQ.0 ) then
         info =-2
      else if( iwantz.EQ.0 ) then
         info =-3
      else if( n.lt.0 ) then
         info =-4
      else if( ilo.lt.1 ) then
         info =-5
      else if( ihi.GT.n .OR. ihi.lt.ilo-1 ) then
         info =-6
      else if( lda.lt.n ) then
         info =-8
      else if( ldb.lt.n ) then
         info =-10
      else if( ldq.lt.1 .OR. ( ilq .and. ldq.lt.n ) ) then
         info =-15
      else if( ldz.lt.1 .OR. ( ilz .and. ldz.lt.n ) ) then
         info =-17
      end if
      if( info.NE.0 ) then
         call xerbla( 'DLAQZ0',-info )
         return
      end if
   
*
*     Quick return if possible
*
      if( n.le.0 ) then
         work( 1 ) = dble( 1 )
         return
      end if

*
*     Get the parameters
*
      jbcmpz(1:1)=wantS
      jbcmpz(2:2)=wantQ
      jbcmpz(3:3)=wantZ

      nwr = ilaenv( 13,'DLAQR0',jbcmpz,n,ilo,ihi,lwork )
      nwr = max( 2,nwr )
      nwr = min( ihi-ilo+1,( n-1 ) / 3,nwr )
      
      nsr = ilaenv( 15,'DLAQR0',jbcmpz,n,ilo,ihi,lwork )
      nsr = min( nsr,( n+6 ) / 9,ihi-ilo )
      nsr = max( 2,nsr-mod( nsr,2 ) )

*     
*     I don't really know how to change ilaenv, the above code
*     doesn't contain all the parameters and the parameters
*     it does set are overwritten in the following line.
*  
      call getparameters(n,nsr,nbr,nwr,nmin,nibble)

      if( n .lt. nmin ) then
         call dhgeqz(wantS,wantQ,wantZ,n,ilo,ihi,A,ldA,B,ldB,alphar,
     $      alphai,beta,Q,ldQ,Z,ldZ,work,lwork,info)
         return
      end if

*
*     Find out required workspace
*

*     Workspace query to dlaqz4
      nw = max(nwr,nmin)
      call dlaqz4(ilschur,ilq,ilz,n,ilo,ihi,nw,A,ldA,B,ldB,Q,ldQ,Z,ldZ,
     $   n_undeflated,n_deflated,alphar,alphai,beta,work,nw,work,nw,
     $   work,-1,aed_info)
      itemp1 = int(work(1))
*     Workspace query to dlaqz5
      call dlaqz5(ilschur,ilq,ilz,n,ilo,ihi,nsr,nbr,alphar,alphai,beta,
     $   A,ldA,B,ldB,Q,ldQ,Z,ldZ,work,nbr,work,nbr,work,-1,sweep_info)
      itemp2 = int(work(1))

      lworkreq = max(itemp1+2*nw**2,itemp2+2*nbr**2)
      if (lwork .eq.-1) then
         work(1) = dble(lworkreq)
         return
      else if (lwork .lt. lworkreq) then
         info =-19
      end if
      if( info.NE.0 ) then
         call xerbla( 'DLAQZ0',-info )
         return
      end if
*
*     Initialize Q and Z
*
      if( iwantq.eq.3 ) call dlaset( 'Full',n,n,zero,one,Q,ldq )
      if( iwantz.eq.3 ) call dlaset( 'Full',n,n,zero,one,Z,ldz )

*     Get machine constants
      safmin = dlamch('SAFE MINIMUM')
      safmax = one/safmin
      call dlabad(safmin,safmax)
      ulp = dlamch('precision')
      smlnum = safmin*(dble(n)/ulp)

      istart = ilo
      istop = ihi
      maxit = 3*(ihi-ilo+1)
      ld = 0

      do iiter = 1,maxit
         if(iiter .ge. maxit) then
            info = istop
            goto 80
         end if
         if (istart+1 .ge. istop) then
            exit
         end if

*        Check deflations at the end
         if (abs(A(istop-1,istop-2)) .le. ulp*(abs(A(istop-1,
     $      istop-1))+abs(A(istop-2,istop-2)))) then
            A(istop-1,istop-2) = zero
            istop = istop-2
            ld = 0
            eshift = zero
         else if (abs(A(istop,istop-1)) .le. ulp*(abs(A(istop,
     $      istop))+abs(A(istop-1,istop-1)))) then
            A(istop,istop-1) = zero
            istop = istop-1
            ld = 0
            eshift = zero
         end if
*        Check deflations at the start
         if (abs(A(istart+2,istart+1)) .le. ulp*(abs(A(istart+1,
     $      istart+1))+abs(A(istart+2,istart+2)))) then
            A(istart+2,istart+1) = zero
            istart = istart+2
            ld = 0
            eshift = zero
         else if (abs(A(istart+1,istart)) .le. ulp*(abs(A(istart,
     $      istart))+abs(A(istart+1,istart+1)))) then
            A(istart+1,istart) = zero
            istart = istart+1
            ld = 0
            eshift = zero
         end if

         if (istart+1 .ge. istop) then
            exit
         end if

*        Check interior deflations
         istart2 = istart
         do k = istop,istart+1,-1
            if (abs(A(k,k-1)) .le. ulp*(abs(A(k,k))+abs(A(k-1,
     $         k-1)))) then
               A(k,k-1) = zero
               istart2 = k
               exit
            end if
         end do

*        Get range to apply rotations to
         if (ilschur) then
            istartm = 1
            istopm = n
         else
            istartm = istart2
            istopm = istop
         end if

*        Check infinite eigenvalues, this is done without blocking so might
*        slow down the method when many infinite eigenvalues are present
         k = istop
         do while (k.ge.istart2)
            temp = zero
            if(k .lt. istop) then
               temp = temp+abs(B(k,k+1))
            end if
            if(k .gt. istart2) then
               temp = temp+abs(B(k-1,k))
            end if

            if(abs(B(k,k)) .lt. max(smlnum,ulp*temp)) then
*              A diagonal element of B is negligable, move it
*              to the top and deflate it
               
               do k2 = k,istart2+1,-1
                  call dlartg(B(k2-1,k2),B(k2-1,k2-1),c1,s1,temp)
                  B(k2-1,k2) = temp
                  B(k2-1,k2-1) = zero

                  call drot(k2-2-istartm+1,B(istartm,k2),1,B(istartm,
     $               k2-1),1,c1,s1)
                  call drot(min(k2+1,istop)-istartm+1,A(istartm,k2),1,
     $               A(istartm,k2-1),1,c1,s1)
                  if (ilz) then
                     call drot(n,Z(1,k2),1,Z(1,k2-1),1,c1,s1)
                  end if

                  if(k2.lt.istop) then
                     call dlartg(A(k2,k2-1),A(k2+1,k2-1),c1,s1,temp)
                     A(k2,k2-1) = temp
                     A(k2+1,k2-1) = zero

                     call drot(istopm-k2+1,A(k2,k2),ldA,A(k2+1,k2),ldA,
     $                  c1,s1)
                     call drot(istopm-k2+1,B(k2,k2),ldB,B(k2+1,k2),ldB,
     $                  c1,s1)
                     if(ilq) then
                        call drot(n,Q(1,k2),1,Q(1,k2+1),1,c1,s1)
                     end if
                  end if

               end do

               if(istart2.lt.istop)then
                  call dlartg(A(istart2,istart2),A(istart2+1,istart2),
     $               c1,s1,temp)
                  A(istart2,istart2) = temp
                  A(istart2+1,istart2) = zero

                  call drot(istopm-(istart2+1)+1,A(istart2,istart2+1),
     $               ldA,A(istart2+1,istart2+1),ldA,c1,s1)
                  call drot(istopm-(istart2+1)+1,B(istart2,istart2+1),
     $               ldB,B(istart2+1,istart2+1),ldB,c1,s1)
                  if(ilq) then
                     call drot(n,Q(1,istart2),1,Q(1,istart2+1),1,c1,s1)
                  end if
               end if

               istart2 = istart2+1
   
            end if
            k = k-1
         end do

*        istart2 now points to the top of the bottom right
*        unreduced Hessenberg block
         if (istart2 .ge. istop) then
            istop = istart2-1
            ld = 0
            eshift = zero
            cycle
         end if

         nw = nwr
         nshifts = nsr
         nblock = nbr

         if (istop-istart2+1 .lt. nmin) then
*           Setting nw to the size of the subblock will make AED deflate
*           all the eigenvalues. This is slightly more efficient than just
*           using qz_small because the off diagonal part gets updated via BLAS.
            if (istop-istart+1 .lt. nmin) then
               nw = istop-istart+1
               istart2 = istart
            else
               nw = istop-istart2+1
            end if
         end if

*
*        Time for AED
*
         call dlaqz4(ilschur,ilq,ilz,n,istart2,istop,nw,A,ldA,B,ldB,Q,
     $      ldQ,Z,ldZ,n_undeflated,n_deflated,alphar,alphai,beta,work,
     $      nw,work(nw**2+1),nw,work(2*nw**2+1),lwork-2*nw**2,aed_info)

         if (n_deflated > 0) then
            istop = istop-n_deflated
            ld = 0
            eshift = zero
         end if

         if (100*n_deflated > nibble*(n_deflated+
     $      n_undeflated) .or.istop-istart2+1 .lt. nmin) then
*           AED has uncovered many eigenvalues. Skip a QZ sweep and run
*           AED again.
            cycle
         end if

         ld = ld+1

         ns = min(nshifts,istop-istart2)
         ns = min(ns,n_undeflated)
         shiftpos = istop-n_deflated-n_undeflated+1

         if (mod(ld,6) .eq. 0) then
* 
*           Exceptional shift.  Chosen for no particularly good reason.
*
            if((dble(maxit)*safmin)*abs(A(istop,
     $         istop-1)).lt.abs(A(istop-1,istop-1))) then
               eshift = A(istop,istop-1)/B(istop-1,istop-1)
            else
               eshift = eshift+one/(safmin*dble(maxit))
            end if
            alphar(shiftpos) = one
            alphar(shiftpos+1) = zero
            alphai(shiftpos) = zero
            alphai(shiftpos+1) = zero
            beta(shiftpos) = eshift
            beta(shiftpos+1) = eshift
            ns = 2
         end if

*
*        Time for a QZ sweep
*
         call dlaqz5(ilschur,ilq,ilz,n,istart2,istop,ns,nblock,
     $      alphar(shiftpos),alphai(shiftpos),beta(shiftpos),A,ldA,B,
     $      ldB,Q,ldQ,Z,ldZ,work,nblock,work(nblock**2+1),nblock,
     $      work(2*nblock**2+1),lwork-2*nblock**2,sweep_info)

      end do

   80 call dlaqz3(ilschur,ilq,ilz,n,ilo,ihi,A,ldA,B,ldB,Q,ldQ,Z,ldZ,
     $alphar,alphai,beta,norm_info)

      end subroutine

* Subroutine: getparameters
*
*------------------------------------------------------------------------------
* DESCRIPTION:
*>  Returns some tuning parameters for a given pencil size
*
*   This subroutine currently only accounts for the size of the unreduced
*   subblock, higher performance might be attained by also accounting for the
*   size of the full matrix.
*
*------------------------------------------------------------------------------
* ARGUMENTS
*>  n          integer [IN]
*                 The size of the subblock
*>  nshifts    integer [OUT]
*                 The desired number of simultanious shifts.
*                 This should always be a multiple of two.
*>  nblock     integer [OUT]
*                 The desired block size to use during the qz sweep.
*                 Good values for this are a multiple of the blas kernel size.
*>  nw         integer [OUT]
*                 The desired deflation window size
*>  nmin       integer [OUT]
*                 Threshold value to switch between double and multishift code.
*>  nibble     integer [OUT]
*                 Threshold value to skip a QZ sweep, expressed as a percentage.
*
*------------------------------------------------------------------------------
      subroutine getparameters(n,nshifts,nblock,nw,nmin,nibble)

      integer,intent(in) :: n
      integer,intent(out) :: nshifts,nblock,nw,nmin,nibble

      real,parameter :: c = 0.1

      integer :: k

      nmin = 20
      nibble = 8

      if (n .lt. 30) then
         nshifts = 2
         nw = 4
      else if (n .lt. 150) then
         nshifts = 4
         nw = 8
      else if (n .lt. 590) then
         nshifts = 32
         nw = 48
      else if (n .lt. 3000) then
         nshifts = 40
         nw = 96
      else if (n .lt. 6000) then
         nshifts = 64
         nw = 96
      else
         nshifts = 64
         nw = 96
      end if

      k = int(nshifts/sqrt( 1+2*nshifts/(c*n) ))
      k = ((k-1)/4)*4+4

      nblock = nshifts+k

      end subroutine getparameters
