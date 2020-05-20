*> \brief \b DLAQZ6
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAQZ6 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqz6.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqz6.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqz6.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAQZ6( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB,
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
*> DLAQZ6 computes the eigenvalues of a real matrix pair (H,T),
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
      subroutine dlaqz6(wantS,wantQ,wantZ,n,ilo,ihi,A,ldA,B,ldB,Q,ldQ,Z,
     $   ldZ,alphar,alphai,beta,info,n_shifts)
      implicit none

*     Arguments
      logical,intent(in) :: wantS,wantQ,wantZ
      integer,intent(in) :: n,ilo,ihi,ldA,ldB,ldQ,ldZ
      integer,intent(inout) :: n_shifts,info

      double precision,intent(inout) :: A(ldA,*),B(ldB,*),Q(ldQ,*),
     $   Z(ldZ,*),alphar(*),alphai(*),beta(*)

*     Parameters
      double precision :: zero,one,half
      parameter(zero=0.0d0,one=1.0d0,half=0.5d0)

*     Local scalars
      double precision :: smlnum,ulp,safmin,safmax,sr1,sr2,si,sb1,sb2,
     $   temp,v(3),c1,s1,c2,s2,H(2,3),eshift
      integer :: istartm,istopm,istart,istop,iiter,maxit,istart2,k,k2,
     $   ld,i
      logical :: bulge

*     External Functions
      double precision,external :: dlamch

*     External Subroutines
      external            :: dlartg

      info = 0

*     Get machine constants
      safmin = dlamch('SAFE MINIMUM')
      safmax = one/safmin
      call dlabad(safmin,safmax)
      ulp = dlamch('precision')
      smlnum = safmin*(dble(n)/ulp)

      istart = ilo
      istop = ihi
      maxit = 30*(ihi-ilo+1)
      ld = 0
      n_shifts = 0

      do iiter = 1,maxit+1
         if(iiter .ge. maxit) then
            info = istop
            goto 80
         end if
         if (istart+1 .ge. istop) then
            exit
         end if

*        Check deflations at the end
         if (abs(A(istop-1,istop-2)) .le. max(smlnum,ulp*(abs(A(istop-1,
     $      istop-1))+abs(A(istop-2,istop-2))) )) then
            A(istop-1,istop-2) = zero
            istop = istop-2
         else if (abs(A(istop,istop-1)) .le. max(smlnum,
     $      ulp*(abs(A(istop,istop))+abs(A(istop-1,istop-1))) )) then
            A(istop,istop-1) = zero
            istop = istop-1
            ld = 0
            eshift = zero
         end if

         if (istart+1 .ge. istop) then
            exit
         end if

*        Check interior deflations
         istart2 = istart
         do k = istop-2,istart+1,-1
            if (abs(A(k,k-1)) .le. max(smlnum,ulp*(abs(A(k,k))+abs(A(k-
     $         1,k-1))) )) then
               A(k,k-1) = zero
               istart2 = k
               ld = 0
               eshift = zero
               exit
            end if
         end do

*        Get range to apply rotations to
         if (wantS) then
            istartm = 1
            istopm = n
         else
            istartm = istart2
            istopm = istop
         end if

*        Check infinite eigenvalues
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
                  if (wantZ) then
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
                     if(wantQ) then
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
                  if(wantQ) then
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

         ld = ld+1
*        Calculate shift
         if (mod(ld,6) .eq. 0) then
* 
*           Exceptional shift.  Chosen for no particularly good reason.
*           (Single shift only.)
*
            if((dble(maxit)*safmin)*abs(A(istop,
     $         istop-1)).lt.abs(A(istop-1,istop-1))) then
               eshift = A(istop,istop-1)/B(istop-1,istop-1)
            else
               eshift = eshift+one/(safmin*dble(maxit))
            end if
            sr1 = one
            sr2 =-one
            si = zero
            sb1 = eshift
            sb2 = eshift
         else
            call dlag2(A(istop-1,istop-1),ldA,B(istop-1,istop-1),ldB,
     $         safmin,sb1,sb2,sr1,sr2,si)
         end if

*        Get range to apply rotations to
         if (wantS) then
            istartm = 1
            istopm = n
         else
            istartm = istart2
            istopm = istop
         end if

*        Introduce double shift
         call dlaqz1(A(istart2,istart2),ldA,B(istart2,istart2),ldB,sr1,
     $      sr2,si,sb1,sb2,v)

         temp = v(2)
         call dlartg(temp,v(3),c1,s1,v(2))
         call dlartg(v(1),v(2),c2,s2,temp)

*        Apply rotations from the left
         call drot(istopm-istart2+1,A(istart2+1,istart2),ldA,
     $      A(istart2+2,istart2),ldA,c1,s1)
         call drot(istopm-istart2+1,A(istart2,istart2),ldA,A(istart2+1,
     $      istart2),ldA,c2,s2)
         call drot(istopm-istart2+1,B(istart2+1,istart2),ldB,
     $      B(istart2+2,istart2),ldB,c1,s1)
         call drot(istopm-istart2+1,B(istart2,istart2),ldB,B(istart2+1,
     $      istart2),ldB,c2,s2)
         if (wantQ) then
         call drot(n,Q(1,istart2+1),1,Q(1,istart2+2),1,c1,s1)
            call drot(n,Q(1,istart2),1,Q(1,istart2+1),1,c2,s2)
         end if

*        The actual QZ sweep
         do k = istart2,istop-2

            call dlaqz2(wantQ,wantZ,k,istartm,istopm,istop,A,ldA,B,ldB,
     $         n,1,Q,ldQ,n,1,Z,ldZ)

         end do

         n_shifts = n_shifts+2

      end do

*     Store eigenvalues and standardize the blocks
   80 call dlaqz3(wantS,wantQ,wantZ,n,istart,ihi,A,ldA,B,ldB,Q,ldQ,Z,
     $ldZ,alphar,alphai,beta,info)

      end subroutine