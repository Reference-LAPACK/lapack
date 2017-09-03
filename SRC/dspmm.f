* Purpose:  To compute alpha*op(A,B)+beta*C, op(A,B) could be either
*           A*B,
*           B*A,
*           A*B**T,
*           B**T*A.
*
* \param[in] SIDE
*            SIDE is CHARACTER*1
*            Select which problem to compute:
*            'L'    C=alpha*B*A+beta*C,
*            'R'    C=alpha*A*B+beta*C.
*
* \param[in] UPLO
*            UPLO is CHARACTER*1
*            Select which part of A is stored:
*            'U'    Upper Triangle,
*            'L'    Lower Triangle.
*
* \param[in] TRAN
*            TRAN is CHARACTER*1
*            Select if B is transverse:
*            'N'    No transverse,
*            'T'    Transverse.
*
* \param[in] M
*            M is INTEGER
*            The size of square matrix A.
*
* \param[in] N
*            N is INTEGER
*            Another dimension of matrix B.
*            For SIDE='L', B=>(N,M),
*            For SIDE='R', B=>(M,N).
*
* \param[in] ALPHA
*            ALPHA is DOUBLEPRECISION
*            The factor.
*
* \param[in] A
*            A is DOUBLEPRECISION(*) array of DIMENSION ((M+1)*M/2)
*
* \param[in] B
*            B is DOUBLEPRECISION(*,*) array of DIMENSION (M,N) or (N,M)
*
* \param[in] LDB
*            LDB is INTEGER
*            The leading dimension of matrix B, should be at least max(1,M) or max(1,N).
*
* \param[in] BETA
*            BETA is DOUBLEPRECISION
*            The factor.
*
* \param[in/out] C
*                C is DOUBLEPRECISION(*,*) array of DIMENSION (M,N) or (N,M)
*
* \param[in] LDC
*            LDC is INTEGER
*            The leading dimension of matrix C, should be at least max(1,M) or max(1,N) based on SIDE.
*
      SUBROUTINE DSPMM(SIDE,UPLO,TRAN,M,N,A,ALPHA,B,LDB,BETA,C,LDC)

      !...INPUT ARGUMENTS...
      CHARACTER SIDE,UPLO,TRAN
      INTEGER M,N,LDB,LDC
      DOUBLEPRECISION ALPHA,BETA,A(*),B(LDB,*),C(LDC,*)

      !...TEMP VARIABLES...
      INTEGER I,J,K,X,Y,Z,DIMA,DIMB,PTYPE,TEMPB,INFO
      DOUBLEPRECISION TEMPA
      LOGICAL S,U,T

      !...TWO CONSTANTS...
      DOUBLEPRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)

      !...EXTERNAL SUBROUTINES...
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA

      !...FLAGS...
      S=LSAME(SIDE,'R')
      U=LSAME(UPLO,'L')
      T=LSAME(TRAN,'N')

      !...CHECK IF ACCEPTABLE ARGUMENTS ARE GIVEN...
      INFO=0
      IF((.NOT.S).AND.(.NOT.LSAME(SIDE,'L')))THEN
          INFO=1
      ELSEIF((.NOT.U).AND.(.NOT.LSAME(UPLO,'U')))THEN
          INFO=2
      ELSEIF((.NOT.T).AND.(.NOT.LSAME(TRAN,'T')))THEN
          INFO=3
      ENDIF

      IF(INFO.NE.0)THEN
          CALL XERBLA('DSPMM ',INFO)
          RETURN
      ENDIF

      !...SWITCH TO PROPER DIMENSION...
      !...THE DIMENSION OF C IS AWAYS (DIMA,DIMB)...
      IF(S)THEN
          DIMA=M
          DIMB=N
      ELSE
          DIMA=N
          DIMB=M
      ENDIF

      !...QUICK RETURN...
      IF(ALPHA.EQ.ZERO)THEN
          IF(BETA.EQ.ZERO)THEN
              DO 20 J=1,DIMB
                  DO 10 I=1,DIMA
                      C(I,J)=ZERO
   10             CONTINUE
   20         CONTINUE
          ELSEIF(BETA.NE.ONE)THEN
              DO 40 J=1,DIMB
                  DO 30 I=1,DIMA
                      C(I,J)=BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          ENDIF
          RETURN
      ENDIF

      !...ALPHA.NE.ZERO...
      !...CHECK beta*C FIRST...
      IF(BETA.EQ.ZERO)THEN
          DO 60 J=1,DIMB
              DO 50 I=1,DIMA
                  C(I,J)=ZERO
   50         CONTINUE
   60     CONTINUE
      ELSEIF(BETA.NE.ONE)THEN
          DO 80 J=1,DIMB
              DO 70 I=1,DIMA
                  C(I,J)=BETA*C(I,J)
   70         CONTINUE
   80     CONTINUE
      ENDIF

      !...ASSIGN PROBLEM TYPE ACCORDING TO GIVEN FLAGS...
      PTYPE=0000
      IF(S)PTYPE=PTYPE+1000
      IF(U)PTYPE=PTYPE+100
      IF(T)PTYPE=PTYPE+10
      IF(ALPHA.EQ.ONE)PTYPE=PTYPE+1

      !...NOW DO MULTIPLICATION FOR DIFFERENT CASES...
      !...USE RELATIVE INCREMENT AS INDICES...
      IF(PTYPE.EQ.1111)THEN
          DO 120 J=1,DIMB
              DO 110 I=1,DIMA
                  TEMPA=ZERO
                  X=I
                  Y=1
                  DO 90 K=M-1,M-I+1,-1
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,J)
                      X=X+K
                      Y=Y+1
   90             CONTINUE
                  DO 100 Z=X,M-I+X
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(Z)*B(Y,J)
                      Y=Y+1
  100             CONTINUE
                  C(I,J)=C(I,J)+TEMPA
  110         CONTINUE
  120     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.1101)THEN
          DO 160 J=1,DIMB
              DO 150 I=1,DIMA
                  TEMPA=ZERO
                  X=I
                  Y=1
                  DO 130 K=M-1,M-I+1,-1
                      IF(B(J,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(J,Y)
                      X=X+K
                      Y=Y+1
  130             CONTINUE
                  DO 140 Z=X,M-I+X
                      IF(B(J,Y).NE.ZERO)TEMPA=TEMPA+A(Z)*B(J,Y)
                      Y=Y+1
  140             CONTINUE
                  C(I,J)=C(I,J)+TEMPA
  150         CONTINUE
  160     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.1110)THEN
          DO 200 J=1,DIMB
              DO 190 I=1,DIMA
                  TEMPA=ZERO
                  X=I
                  Y=1
                  DO 170 K=M-1,M-I+1,-1
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,J)
                      X=X+K
                      Y=Y+1
  170             CONTINUE
                  DO 180 Z=X,M-I+X
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(Z)*B(Y,J)
                      Y=Y+1
  180             CONTINUE
                  C(I,J)=C(I,J)+TEMPA*ALPHA
  190         CONTINUE
  200     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.1100)THEN
          DO 240 J=1,DIMB
              DO 230 I=1,DIMA
                  TEMPA=ZERO
                  X=I
                  Y=1
                  DO 210 K=M-1,M-I+1,-1
                      IF(B(J,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(J,Y)
                      X=X+K
                      Y=Y+1
  210             CONTINUE
                  DO 220 Z=X,M-I+X
                      IF(B(J,Y).NE.ZERO)TEMPA=TEMPA+A(Z)*B(J,Y)
                      Y=Y+1
  220             CONTINUE
                  C(I,J)=C(I,J)+TEMPA*ALPHA
  230         CONTINUE
  240     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.1011)THEN
          DO 280 J=1,DIMB
              DO 270 I=1,DIMA
                  TEMPA=ZERO
                  TEMPB=(I-1)*I/2
                  Y=1
                  DO 250 X=TEMPB+1,TEMPB+I
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,J)
                      Y=Y+1
  250             CONTINUE
                  X=X-1
                  DO 260 Z=I,M-1
                      X=X+Z
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,J)
                      Y=Y+1
  260             CONTINUE
                  C(I,J)=C(I,J)+TEMPA
  270         CONTINUE
  280     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.1001)THEN
          DO 320 J=1,DIMB
              DO 310 I=1,DIMA
                  TEMPA=ZERO
                  TEMPB=(I-1)*I/2
                  Y=1
                  DO 290 X=TEMPB+1,TEMPB+I
                      IF(B(J,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(J,Y)
                      Y=Y+1
  290             CONTINUE
                  X=X-1
                  DO 300 Z=I,M-1
                      X=X+Z
                      IF(B(J,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(J,Y)
                      Y=Y+1
  300             CONTINUE
                  C(I,J)=C(I,J)+TEMPA
  310         CONTINUE
  320     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.1010)THEN
          DO 360 J=1,DIMB
              DO 350 I=1,DIMA
                  TEMPA=ZERO
                  TEMPB=(I-1)*I/2
                  Y=1
                  DO 330 X=TEMPB+1,TEMPB+I
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,J)
                      Y=Y+1
  330             CONTINUE
                  X=X-1
                  DO 340 Z=I,M-1
                      X=X+Z
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,J)
                      Y=Y+1
  340             CONTINUE
                  C(I,J)=C(I,J)+TEMPA*ALPHA
  350         CONTINUE
  360     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.1000)THEN
          DO 400 J=1,DIMB
              DO 390 I=1,DIMA
                  TEMPA=ZERO
                  TEMPB=(I-1)*I/2
                  Y=1
                  DO 370 X=TEMPB+1,TEMPB+I
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,J)
                      Y=Y+1
  370             CONTINUE
                  X=X-1
                  DO 380 Z=I,M-1
                      X=X+Z
                      IF(B(Y,J).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,J)
                      Y=Y+1
  380             CONTINUE
                  C(I,J)=C(I,J)+TEMPA*ALPHA
  390         CONTINUE
  400     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.0111)THEN
          DO 440 J=1,DIMB
              DO 430 I=1,DIMA
                  TEMPA=ZERO
                  X=J
                  Y=1
                  DO 410 K=M-1,M-J+1,-1
                      IF(B(I,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(I,Y)
                      X=X+K
                      Y=Y+1
  410             CONTINUE
                  DO 420 Z=X,M-J+X
                      IF(B(I,Y).NE.ZERO)TEMPA=TEMPA+A(Z)*B(I,Y)
                      Y=Y+1
  420             CONTINUE
                  C(I,J)=C(I,J)+TEMPA
  430         CONTINUE
  440     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.0101)THEN
          DO 480 J=1,DIMB
              DO 470 I=1,DIMA
                  TEMPA=ZERO
                  X=J
                  Y=1
                  DO 450 K=M-1,M-J+1,-1
                      IF(B(Y,I).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,I)
                      X=X+K
                      Y=Y+1
  450             CONTINUE
                  DO 460 Z=X,M-J+X
                      IF(B(Y,I).NE.ZERO)TEMPA=TEMPA+A(Z)*B(Y,I)
                      Y=Y+1
  460             CONTINUE
                  C(I,J)=C(I,J)+TEMPA
  470         CONTINUE
  480     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.0110)THEN
          DO 520 J=1,DIMB
              DO 510 I=1,DIMA
                  TEMPA=ZERO
                  X=J
                  Y=1
                  DO 490 K=M-1,M-J+1,-1
                      IF(B(I,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(I,Y)
                      X=X+K
                      Y=Y+1
  490             CONTINUE
                  DO 500 Z=X,M-J+X
                      IF(B(I,Y).NE.ZERO)TEMPA=TEMPA+A(Z)*B(I,Y)
                      Y=Y+1
  500             CONTINUE
                  C(I,J)=C(I,J)+TEMPA*ALPHA
  510         CONTINUE
  520     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.0100)THEN
          DO 560 J=1,DIMB
              DO 550 I=1,DIMA
                  TEMPA=ZERO
                  X=J
                  Y=1
                  DO 530 K=M-1,M-J+1,-1
                      IF(B(Y,I).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,I)
                      X=X+K
                      Y=Y+1
  530             CONTINUE
                  DO 540 Z=X,M-J+X
                      IF(B(Y,I).NE.ZERO)TEMPA=TEMPA+A(Z)*B(Y,I)
                      Y=Y+1
  540             CONTINUE
                  C(I,J)=C(I,J)+TEMPA*ALPHA
  550         CONTINUE
  560     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.0011)THEN
          DO 600 J=1,DIMB
              DO 590 I=1,DIMA
                  TEMPA=ZERO
                  TEMPB=(J-1)*J/2
                  Y=1
                  DO 570 X=TEMPB+1,TEMPB+J
                      IF(B(I,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(I,Y)
                      Y=Y+1
  570             CONTINUE
                  X=X-1
                  DO 580 Z=J,M-1
                      X=X+Z
                      IF(B(I,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(I,Y)
                      Y=Y+1
  580             CONTINUE
                  C(I,J)=C(I,J)+TEMPA
  590         CONTINUE
  600     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.0001)THEN
          DO 640 J=1,DIMB
              DO 630 I=1,DIMA
                  TEMPA=ZERO
                  TEMPB=(J-1)*J/2
                  Y=1
                  DO 610 X=TEMPB+1,TEMPB+J
                      IF(B(Y,I).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,I)
                      Y=Y+1
  610             CONTINUE
                  X=X-1
                  DO 620 Z=J,M-1
                      X=X+Z
                      IF(B(Y,I).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,I)
                      Y=Y+1
  620             CONTINUE
                  C(I,J)=C(I,J)+TEMPA
  630         CONTINUE
  640     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.0010)THEN
          DO 680 J=1,DIMB
              DO 670 I=1,DIMA
                  TEMPA=ZERO
                  TEMPB=(J-1)*J/2
                  Y=1
                  DO 650 X=TEMPB+1,TEMPB+J
                      IF(B(I,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(I,Y)
                      Y=Y+1
  650             CONTINUE
                  X=X-1
                  DO 660 Z=J,M-1
                      X=X+Z
                      IF(B(I,Y).NE.ZERO)TEMPA=TEMPA+A(X)*B(I,Y)
                      Y=Y+1
  660             CONTINUE
                  C(I,J)=C(I,J)+TEMPA*ALPHA
  670         CONTINUE
  680     CONTINUE
          RETURN
      ENDIF

      IF(PTYPE.EQ.0000)THEN
          DO 720 J=1,DIMB
              DO 710 I=1,DIMA
                  TEMPA=ZERO
                  TEMPB=(J-1)*J/2
                  Y=1
                  DO 690 X=TEMPB+1,TEMPB+J
                      IF(B(Y,I).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,I)
                      Y=Y+1
  690             CONTINUE
                  X=X-1
                  DO 700 Z=J,M-1
                      X=X+Z
                      IF(B(Y,I).NE.ZERO)TEMPA=TEMPA+A(X)*B(Y,I)
                      Y=Y+1
  700             CONTINUE
                  C(I,J)=C(I,J)+TEMPA*ALPHA
  710         CONTINUE
  720     CONTINUE
          RETURN
      ENDIF

      END
