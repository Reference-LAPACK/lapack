      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
C
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C

C     .. Scalar Arguments ..
      REAL AI,AR,BI,BR,CI,CR
C     ..
C     .. Local Scalars ..
      REAL AIS,ARS,BIS,BRS,S
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS+AIS*BIS)/S
      CI = (AIS*BRS-ARS*BIS)/S
      END
      SUBROUTINE CINVIT(NM,N,AR,AI,WR,WI,SELECT,MM,M,ZR,ZI,IERR,RM1,RM2,
     +                  RV1,RV2)
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE CX INVIT
C     BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP. VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A COMPLEX UPPER
C     HESSENBERG MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVALUES OF THE MATRIX.  THE EIGENVALUES MUST BE
C          STORED IN A MANNER IDENTICAL TO THAT OF SUBROUTINE  COMLR,
C          WHICH RECOGNIZES POSSIBLE SPLITTING OF THE MATRIX.
C
C        SELECT SPECIFIES THE EIGENVECTORS TO BE FOUND.  THE
C          EIGENVECTOR CORRESPONDING TO THE J-TH EIGENVALUE IS
C          SPECIFIED BY SETTING SELECT(J) TO .TRUE..
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          EIGENVECTORS TO BE FOUND.
C
C     ON OUTPUT
C
C        AR, AI, WI, AND SELECT ARE UNALTERED.
C
C        WR MAY HAVE BEEN ALTERED SINCE CLOSE EIGENVALUES ARE PERTURBED
C          SLIGHTLY IN SEARCHING FOR INDEPENDENT EIGENVECTORS.
C
C        M IS THE NUMBER OF EIGENVECTORS ACTUALLY FOUND.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVECTORS.  THE EIGENVECTORS ARE NORMALIZED
C          SO THAT THE COMPONENT OF LARGEST MAGNITUDE IS 1.
C          ANY VECTOR WHICH FAILS THE ACCEPTANCE TEST IS SET TO ZERO.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -(2*N+1)   IF MORE THAN MM EIGENVECTORS HAVE BEEN SPECIFIED,
C          -K         IF THE ITERATION CORRESPONDING TO THE K-TH
C                     VALUE FAILS,
C          -(N+K)     IF BOTH ERROR SITUATIONS OCCUR.
C
C        RM1, RM2, RV1, AND RV2 ARE TEMPORARY STORAGE ARRAYS.
C
C     THE ALGOL PROCEDURE GUESSVEC APPEARS IN CINVIT IN LINE.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
*
*     GET ULP FROM SLAMCH FOR NEW SMALL PERTURBATION AS IN LAPACK

C     .. Scalar Arguments ..
      INTEGER IERR,M,MM,N,NM
C     ..
C     .. Array Arguments ..
      REAL AI(NM,N),AR(NM,N),RM1(N,N),RM2(N,N),RV1(N),RV2(N),WI(N),
     +     WR(N),ZI(NM,MM),ZR(NM,MM)
      LOGICAL SELECT(N)
C     ..
C     .. Local Scalars ..
      REAL EPS3,GROWTO,ILAMBD,NORM,NORMV,OPST,RLAMBD,UKROOT,ULP,X,Y
      INTEGER I,II,IP1,ITS,J,K,KM1,MP,S,UK
C     ..
C     .. External Functions ..
      REAL PYTHAG,SLAMCH
      EXTERNAL PYTHAG,SLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL CDIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
      IF (N.GT.0) THEN
          ULP = SLAMCH('EPSILON')
C
*
*     INITIALIZE
          OPST = 0
          IERR = 0
          UK = 0
          S = 1
C
          DO 240 K = 1,N
              IF (SELECT(K)) THEN
                  IF (S.GT.MM) THEN
                      GO TO 250
                  ELSE
                      IF (UK.LT.K) THEN
C     .......... CHECK FOR POSSIBLE SPLITTING ..........
                          DO 10 UK = K,N
                              IF (UK.EQ.N) THEN
                                  GO TO 20
                              ELSE IF (AR(UK+1,UK).EQ.0.0E0 .AND.
     +                                 AI(UK+1,UK).EQ.0.0E0) THEN
                                  GO TO 20
                              END IF
   10                     CONTINUE
C     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK
C                (HESSENBERG) MATRIX ..........
   20                     NORM = 0.0E0
                          MP = 1
C
*
*        INCREMENT OPCOUNT FOR LOOP 180
                          OPS = OPS + 6*UK* (UK-1)
                          DO 40 I = 1,UK
                              X = 0.0E0
C
                              DO 30 J = MP,UK
                                  X = X + PYTHAG(AR(I,J),AI(I,J))
   30                         CONTINUE
C
                              IF (X.GT.NORM) NORM = X
                              MP = I
   40                     CONTINUE
C     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
C                AND CLOSE ROOTS ARE MODIFIED BY EPS3 ..........
                          IF (NORM.EQ.0.0E0) NORM = 1.0E0
*         EPS3 = EPSLON(NORM)
*
*        INCREMENT OPCOUNT FOR EPS3, UKROOT
                          OPST = OPST + 3
                          EPS3 = NORM*ULP
C     .......... GROWTO IS THE CRITERION FOR GROWTH ..........
                          UKROOT = UK
                          UKROOT = SQRT(UKROOT)
                          GROWTO = 0.1E0/UKROOT
                      END IF
                      RLAMBD = WR(K)
                      ILAMBD = WI(K)
                      IF (K.NE.1) THEN
                          KM1 = K - 1
   50                     CONTINUE
C     .......... FOR I=K-1 STEP -1 UNTIL 1 DO -- ..........
                          DO 60 II = 1,KM1
                              I = K - II
                              IF (SELECT(I) .AND.
     +                            ABS(WR(I)-RLAMBD).LT.EPS3 .AND.
     +                            ABS(WI(I)-ILAMBD).LT.EPS3) GO TO 70
   60                     CONTINUE
                          GO TO 80
C     .......... PERTURB EIGENVALUE IF IT IS CLOSE
C                TO ANY PREVIOUS EIGENVALUE ..........
   70                     RLAMBD = RLAMBD + EPS3
                          GO TO 50
C
*
*        INCREMENT OPCOUNT FOR LOOP 260.
   80                     OPST = OPST + 2* (K-1)
                          WR(K) = RLAMBD
                      END IF
C     .......... FORM UPPER HESSENBERG (AR,AI)-(RLAMBD,ILAMBD)*I
C                AND INITIAL COMPLEX VECTOR ..........
                      MP = 1
C
*
*        INCREMENT OP COUNT FOR LOOP 320
                      OPS = OPS + 2*UK
                      DO 100 I = 1,UK
C
                          DO 90 J = MP,UK
                              RM1(I,J) = AR(I,J)
                              RM2(I,J) = AI(I,J)
   90                     CONTINUE
C
                          RM1(I,I) = RM1(I,I) - RLAMBD
                          RM2(I,I) = RM2(I,I) - ILAMBD
                          MP = I
                          RV1(I) = EPS3
  100                 CONTINUE
C     .......... TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3 ..........
                      IF (UK.NE.1) THEN
C
*
*        INCREMENT OP COUNT FOR LOOP 400
                          OPS = OPS + (52+4*UK)* (UK-1)
                          DO 130 I = 2,UK
                              MP = I - 1
                              IF (PYTHAG(RM1(I,MP),RM2(I,MP)).GT.
     +                            PYTHAG(RM1(MP,MP),RM2(MP,MP))) THEN
C
                                  DO 110 J = MP,UK
                                      Y = RM1(I,J)
                                      RM1(I,J) = RM1(MP,J)
                                      RM1(MP,J) = Y
                                      Y = RM2(I,J)
                                      RM2(I,J) = RM2(MP,J)
                                      RM2(MP,J) = Y
  110                             CONTINUE
                              END IF
C
                              IF (RM1(MP,MP).EQ.0.0E0 .AND.
     +                            RM2(MP,MP).EQ.0.0E0) RM1(MP,MP) = EPS3
                              CALL CDIV(RM1(I,MP),RM2(I,MP),RM1(MP,MP),
     +                                  RM2(MP,MP),X,Y)
                              IF (X.NE.0.0E0 .OR. Y.NE.0.0E0) THEN
C
                                  DO 120 J = I,UK
                                      RM1(I,J) = RM1(I,J) -
     +                                           X*RM1(MP,J) +
     +                                           Y*RM2(MP,J)
                                      RM2(I,J) = RM2(I,J) -
     +                                           X*RM2(MP,J) -
     +                                           Y*RM1(MP,J)
  120                             CONTINUE
                              END IF
  130                     CONTINUE
C
                      END IF
C
                      IF (RM1(UK,UK).EQ.0.0E0 .AND.
     +                    RM2(UK,UK).EQ.0.0E0) RM1(UK,UK) = EPS3
                      ITS = 0
  140                 CONTINUE
C     .......... BACK SUBSTITUTION
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
                      DO 160 II = 1,UK
                          I = UK + 1 - II
                          X = RV1(I)
                          Y = 0.0E0
                          IF (I.NE.UK) THEN
                              IP1 = I + 1
C
                              DO 150 J = IP1,UK
                                  X = X - RM1(I,J)*RV1(J) +
     +                                RM2(I,J)*RV2(J)
                                  Y = Y - RM1(I,J)*RV2(J) -
     +                                RM2(I,J)*RV1(J)
  150                         CONTINUE
                          END IF
C
                          CALL CDIV(X,Y,RM1(I,I),RM2(I,I),RV1(I),RV2(I))
  160                 CONTINUE
*
*        INCREMENT OP COUNT FOR BACK SUBSTITUTION LOOP 720
                      OPS = OPS + 4*UK* (UK+3)
C     .......... ACCEPTANCE TEST FOR EIGENVECTOR
C                AND NORMALIZATION ..........
                      ITS = ITS + 1
                      NORM = 0.0E0
                      NORMV = 0.0E0
C
*
*        INCREMENT OP COUNT ACCEPTANCE TEST
                      OPS = OPS + 19*UK
                      DO 170 I = 1,UK
                          X = PYTHAG(RV1(I),RV2(I))
                          IF (NORMV.LT.X) THEN
                              NORMV = X
                              J = I
                          END IF
                          NORM = NORM + X
  170                 CONTINUE
C
                      IF (NORM.LT.GROWTO) THEN
C     .......... IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR ..........
                          IF (ITS.GE.UK) THEN
                              GO TO 200
                          ELSE
                              X = UKROOT
                              Y = EPS3/ (X+1.0E0)
                              RV1(1) = EPS3
C
                              DO 180 I = 2,UK
                                  RV1(I) = Y
  180                         CONTINUE
C
                              J = UK - ITS + 1
                              RV1(J) = RV1(J) - EPS3*X
                              GO TO 140
                          END IF
                      END IF
C     .......... ACCEPT VECTOR ..........
                      X = RV1(J)
                      Y = RV2(J)
C
*
*        INCREMENT OP COUNT ACCEPT VECTOR LOOP 820
                      OPS = OPS + 16*UK
                      DO 190 I = 1,UK
                          CALL CDIV(RV1(I),RV2(I),X,Y,ZR(I,S),ZI(I,S))
  190                 CONTINUE
C
                      IF (UK.EQ.N) THEN
                          GO TO 230
                      ELSE
                          J = UK + 1
                          GO TO 210
                      END IF
C     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  200                 J = 1
                      IERR = -K
C*PL*ERROR* Embedded comment after label moved
C     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  210                 DO 220 I = J,N
                          ZR(I,S) = 0.0E0
                          ZI(I,S) = 0.0E0
  220                 CONTINUE
C
  230                 S = S + 1
                  END IF
              END IF
  240     CONTINUE
C
          GO TO 260
C     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
C                SPACE REQUIRED ..........
  250     IF (IERR.NE.0) IERR = IERR - N
          IF (IERR.EQ.0) IERR = - (2*N+1)
  260     M = S - 1
*
*     COMPUTE FINAL OP COUNT
          OPS = OPS + OPST
      END IF
      END
      SUBROUTINE COMQR(NM,N,LOW,IGH,HR,HI,WR,WI,IERR)
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR, NUM. MATH. 12, 369-376(1968) BY MARTIN
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A COMPLEX
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN
C          INFORMATION ABOUT THE UNITARY TRANSFORMATIONS USED IN
C          THE REDUCTION BY  CORTH, IF PERFORMED.
C
C     ON OUTPUT
C
C        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
C          DESTROYED.  THEREFORE, THEY MUST BE SAVED BEFORE
C          CALLING  COMQR  IF SUBSEQUENT CALCULATION OF
C          EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
*
*
C     .. Scalar Arguments ..
      INTEGER IERR,IGH,LOW,N,NM
C     ..
C     .. Array Arguments ..
      REAL HI(NM,N),HR(NM,N),WI(N),WR(N)
C     ..
C     .. Local Scalars ..
      REAL NORM,OPST,OVFL,SI,SMALL,SMLNUM,SR,TI,TR,TST1,TST2,ULP,UNFL,
     +     XI,XR,YI,YR,ZZI,ZZR
      INTEGER EN,ENM1,I,ITN,ITS,J,L,LL,LP1
C     ..
C     .. External Functions ..
      REAL PYTHAG,SLAMCH
      EXTERNAL PYTHAG,SLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL CDIV,CSROOT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,MIN0
C     ..
      IF (N.GT.0) THEN
*
*     COMPUTE THE 1-NORM OF MATRIX H
*
          NORM = 0.0E0
          DO 20 J = LOW,IGH
              SR = 0.0E0
              DO 10 I = LOW,MIN(IGH,J+1)
                  SR = SR + PYTHAG(HR(I,J),HI(I,J))
   10         CONTINUE
              NORM = MAX(NORM,SR)
   20     CONTINUE
*
*     GET SMALL FOR NEW CONVERGENCE CRITERION AS IN LAPACK
*
          UNFL = SLAMCH('SAFE MINIMUM')
          OVFL = SLAMCH('OVERFLOW')
          ULP = SLAMCH('EPSILON')*SLAMCH('BASE')
          SMLNUM = MAX(UNFL* (N/ULP),N/ (ULP*OVFL))
          SMALL = MAX(SMLNUM,ULP*NORM)
*
*
*     INITIALIZE
          ITCNT = 0
          OPST = 0
          IERR = 0
          IF (LOW.NE.IGH) THEN
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
              L = LOW + 1
C
*
*        INCREMENT OP COUNT FOR LOOP 170
              OPS = OPS + (6* (IGH-LOW+1)+32)* (IGH-L+1)
              DO 50 I = L,IGH
                  LL = MIN0(I+1,IGH)
                  IF (HI(I,I-1).NE.0.0E0) THEN
                      NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
                      YR = HR(I,I-1)/NORM
                      YI = HI(I,I-1)/NORM
                      HR(I,I-1) = NORM
                      HI(I,I-1) = 0.0E0
C
                      DO 30 J = I,IGH
                          SI = YR*HI(I,J) - YI*HR(I,J)
                          HR(I,J) = YR*HR(I,J) + YI*HI(I,J)
                          HI(I,J) = SI
   30                 CONTINUE
C
                      DO 40 J = LOW,LL
                          SI = YR*HI(J,I) + YI*HR(J,I)
                          HR(J,I) = YR*HR(J,I) - YI*HI(J,I)
                          HI(J,I) = SI
   40                 CONTINUE
                  END IF
   50         CONTINUE
C
          END IF
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
          DO 60 I = 1,N
              IF (I.LT.LOW .OR. I.GT.IGH) THEN
                  WR(I) = HR(I,I)
                  WI(I) = HI(I,I)
              END IF
   60     CONTINUE
C
          EN = IGH
          TR = 0.0E0
          TI = 0.0E0
          ITN = 30*N
   70     CONTINUE
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
          IF (EN.LT.LOW) THEN
              GO TO 180
          ELSE
              ITS = 0
              ENM1 = EN - 1
   80         CONTINUE
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW E0 -- ..........
              DO 90 LL = LOW,EN
                  L = EN + LOW - LL
                  IF (L.EQ.LOW) THEN
                      GO TO 100
                  ELSE
                      TST1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1)) +
     +                       ABS(HR(L,L)) + ABS(HI(L,L))
*         TST2 = TST1 + ABS(HR(L,L-1))
*         IF (TST2 .EQ. TST1) GO TO 300
                      TST2 = ABS(HR(L,L-1))
                      IF (TST2.LE.MIN(ULP*TST1,SMALL)) GO TO 100
                  END IF
   90         CONTINUE
C     .......... FORM SHIFT ..........
  100         CONTINUE
*
*        INCREMENT OP COUNT FOR CONVERGENCE TEST
              OPS = OPS + 4* (EN-L+1)
              IF (L.NE.EN) THEN
                  IF (ITN.EQ.0) THEN
                      GO TO 170
                  ELSE
                      IF (ITS.EQ.10 .OR. ITS.EQ.20) THEN
C     .......... FORM EXCEPTIONAL SHIFT ..........
                          SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
                          SI = 0.0E0
                      ELSE
*
*        INCREMENT OPCOUNT FOR FOMING SHIFT
                          OPST = OPST + 58
                          SR = HR(EN,EN)
                          SI = HI(EN,EN)
                          XR = HR(ENM1,EN)*HR(EN,ENM1)
                          XI = HI(ENM1,EN)*HR(EN,ENM1)
                          IF (XR.NE.0.0E0 .OR. XI.NE.0.0E0) THEN
                              YR = (HR(ENM1,ENM1)-SR)/2.0E0
                              YI = (HI(ENM1,ENM1)-SI)/2.0E0
                              CALL CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,
     +                                    ZZR,ZZI)
                              IF (YR*ZZR+YI*ZZI.LT.0.0E0) THEN
                                  ZZR = -ZZR
                                  ZZI = -ZZI
                              END IF
                              CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
                              SR = SR - XR
                              SI = SI - XI
                          END IF
                      END IF
C
                      DO 110 I = LOW,EN
                          HR(I,I) = HR(I,I) - SR
                          HI(I,I) = HI(I,I) - SI
  110                 CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 360
                      OPS = OPS + 2*EN
C
                      TR = TR + SR
                      TI = TI + SI
                      ITS = ITS + 1
                      ITN = ITN - 1
*
*       UPDATE ITERATION NUMBER
                      ITCNT = 30*N - ITN
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
                      LP1 = L + 1
C
*
*        INCREMENT OPCOUNT FOR REDUCING TO TRIANGULAR, LOOP 500
                      OPS = OPS + (EN-LP1+1)* (61+10* (EN-LP1))
                      DO 130 I = LP1,EN
                          SR = HR(I,I-1)
                          HR(I,I-1) = 0.0E0
                          NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),
     +                           SR)
                          XR = HR(I-1,I-1)/NORM
                          WR(I-1) = XR
                          XI = HI(I-1,I-1)/NORM
                          WI(I-1) = XI
                          HR(I-1,I-1) = NORM
                          HI(I-1,I-1) = 0.0E0
                          HI(I,I-1) = SR/NORM
C
                          DO 120 J = I,EN
                              YR = HR(I-1,J)
                              YI = HI(I-1,J)
                              ZZR = HR(I,J)
                              ZZI = HI(I,J)
                              HR(I-1,J) = XR*YR + XI*YI + HI(I,I-1)*ZZR
                              HI(I-1,J) = XR*YI - XI*YR + HI(I,I-1)*ZZI
                              HR(I,J) = XR*ZZR - XI*ZZI - HI(I,I-1)*YR
                              HI(I,J) = XR*ZZI + XI*ZZR - HI(I,I-1)*YI
  120                     CONTINUE
  130                 CONTINUE
C
C
                      SI = HI(EN,EN)
                      IF (SI.NE.0.0E0) THEN
                          NORM = PYTHAG(HR(EN,EN),SI)
                          SR = HR(EN,EN)/NORM
                          SI = SI/NORM
                          HR(EN,EN) = NORM
                          HI(EN,EN) = 0.0E0
*
*        INCREMENT OPCOUNT
                          OPST = OPST + 20
                      END IF
C     .......... INVERSE OPERATION (COLUMNS) ..........
                      DO 150 J = LP1,EN
                          XR = WR(J-1)
                          XI = WI(J-1)
C
                          DO 140 I = L,J
                              YR = HR(I,J-1)
                              YI = 0.0E0
                              ZZR = HR(I,J)
                              ZZI = HI(I,J)
                              IF (I.NE.J) THEN
                                  YI = HI(I,J-1)
                                  HI(I,J-1) = XR*YI + XI*YR +
     +                                        HI(J,J-1)*ZZI
                              END IF
                              HR(I,J-1) = XR*YR - XI*YI + HI(J,J-1)*ZZR
                              HR(I,J) = XR*ZZR + XI*ZZI - HI(J,J-1)*YR
                              HI(I,J) = XR*ZZI - XI*ZZR - HI(J,J-1)*YI
  140                     CONTINUE
  150                 CONTINUE
C
*
*        INCREMENT OPCOUNT FOR INVERSE OPERATION LOOP 600
                      OPS = OPS + 10* (EN-LP1+1)* (EN+LP1)
C
                      IF (SI.NE.0.0E0) THEN
C
*
*        INCREMENT OP COUNT FOR LOOP 630
                          OPS = OPS + 6* (EN-L+1)
                          DO 160 I = L,EN
                              YR = HR(I,EN)
                              YI = HI(I,EN)
                              HR(I,EN) = SR*YR - SI*YI
                              HI(I,EN) = SR*YI + SI*YR
  160                     CONTINUE
C
                      END IF
                      GO TO 80
                  END IF
              END IF
C     .......... A ROOT FOUND ..........
              WR(EN) = HR(EN,EN) + TR
              WI(EN) = HI(EN,EN) + TI
              EN = ENM1
              GO TO 70
          END IF
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
  170     IERR = EN
*
*     COMPUTE FINAL OP COUNT
  180     OPS = OPS + OPST
      END IF
      END
      SUBROUTINE COMQR2(NM,N,LOW,IGH,ORTR,ORTI,HR,HI,WR,WI,ZR,ZI,IERR)
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR2, NUM. MATH. 16, 181-204(1970) BY PETERS
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE QR
C     METHOD.  THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX
C     CAN ALSO BE FOUND IF  CORTH  HAS BEEN USED TO REDUCE
C     THIS GENERAL MATRIX TO HESSENBERG FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        ORTR AND ORTI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  CORTH, IF PERFORMED.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, SET ORTR(J) AND
C          ORTI(J) TO 0.0E0 FOR THESE ELEMENTS.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN FURTHER
C          INFORMATION ABOUT THE TRANSFORMATIONS WHICH WERE USED IN THE
C          REDUCTION BY  CORTH, IF PERFORMED.  IF THE EIGENVECTORS OF
C          THE HESSENBERG MATRIX ARE DESIRED, THESE ELEMENTS MAY BE
C          ARBITRARY.
C
C     ON OUTPUT
C
C        ORTR, ORTI, AND THE UPPER HESSENBERG PORTIONS OF HR AND HI
C          HAVE BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS.  THE EIGENVECTORS
C          ARE UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE OF
C          THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
*     THE ORIGINAL DO STATEMENTS
*
*         DO 840 I = 1, ENM1
*         DO 820 J = IP1, N
*         DO 880 JJ = LOW, ENM1
*
*     HAVE BEEN CHANGED TO
*
*         DO 840 I = 1, N
*         DO 820 J = I, N
*         DO 880 JJ = LOW, N
*
*     ACCORDING TO BURT GARBOW'S SUGGESTION ON NA-NET.
*     ZHAOJUN BAI, NOV.28, 1989
C     ------------------------------------------------------------------
C
*
*
C     .. Scalar Arguments ..
      INTEGER IERR,IGH,LOW,N,NM
C     ..
C     .. Array Arguments ..
      REAL HI(NM,N),HR(NM,N),ORTI(IGH),ORTR(IGH),WI(N),WR(N),ZI(NM,N),
     +     ZR(NM,N)
C     ..
C     .. Local Scalars ..
      REAL NORM,OPST,OVFL,SI,SMALL,SMLNUM,SR,TI,TR,TST1,TST2,ULP,UNFL,
     +     XI,XR,YI,YR,ZZI,ZZR
      INTEGER EN,ENM1,I,IEND,II,IP1,ITN,ITS,J,JJ,K,L,LL,LP1,M,NN
C     ..
C     .. External Functions ..
      REAL PYTHAG,SLAMCH
      EXTERNAL PYTHAG,SLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL CDIV,CSROOT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,MIN0
C     ..
      IF (N.GT.0) THEN
*
*     COMPUTE THE 1-NORM OF MATRIX H
*
          NORM = 0.0E0
          DO 20 J = 1,N
              SR = 0.0E0
              DO 10 I = 1,MIN(N,J+1)
                  SR = SR + PYTHAG(HR(I,J),HI(I,J))
   10         CONTINUE
              NORM = MAX(NORM,SR)
   20     CONTINUE
*
*     GET SMALL FOR NEW CONVERGENCE CRITERION AS IN LAPACK
*
          UNFL = SLAMCH('SAFE MINIMUM')
          OVFL = SLAMCH('OVERFLOW')
          ULP = SLAMCH('EPSILON')*SLAMCH('BASE')
          SMLNUM = MAX(UNFL* (N/ULP),N/ (ULP*OVFL))
          SMALL = MAX(SMLNUM,ULP*NORM)
*
*
*     INITIALIZE
          ITCNT = 0
          OPST = 0
          IERR = 0
C     .......... INITIALIZE EIGENVECTOR MATRIX ..........
          DO 40 J = 1,N
C
              DO 30 I = 1,N
                  ZR(I,J) = 0.0E0
                  ZI(I,J) = 0.0E0
   30         CONTINUE
              ZR(J,J) = 1.0E0
   40     CONTINUE
C     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY CORTH ..........
          IEND = IGH - LOW - 1
          IF (IEND) 160,110,50
C*PL*ERROR* Embedded comment after label moved
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
   50     DO 100 II = 1,IEND
              I = IGH - II
              IF (ORTR(I).NE.0.0E0 .OR. ORTI(I).NE.0.0E0) THEN
                  IF (HR(I,I-1).NE.0.0E0 .OR. HI(I,I-1).NE.0.0E0) THEN
C     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
                      NORM = HR(I,I-1)*ORTR(I) + HI(I,I-1)*ORTI(I)
                      IP1 = I + 1
C
                      DO 60 K = IP1,IGH
                          ORTR(K) = HR(K,I-1)
                          ORTI(K) = HI(K,I-1)
   60                 CONTINUE
C
*
*        INCREMENT OP COUNT FOR LOOP 130
                      OPS = OPS + (16* (IGH-I+1)+2)* (IGH-I+1)
                      DO 90 J = I,IGH
                          SR = 0.0E0
                          SI = 0.0E0
C
                          DO 70 K = I,IGH
                              SR = SR + ORTR(K)*ZR(K,J) +
     +                             ORTI(K)*ZI(K,J)
                              SI = SI + ORTR(K)*ZI(K,J) -
     +                             ORTI(K)*ZR(K,J)
   70                     CONTINUE
C
                          SR = SR/NORM
                          SI = SI/NORM
C
                          DO 80 K = I,IGH
                              ZR(K,J) = ZR(K,J) + SR*ORTR(K) -
     +                                  SI*ORTI(K)
                              ZI(K,J) = ZI(K,J) + SR*ORTI(K) +
     +                                  SI*ORTR(K)
   80                     CONTINUE
   90                 CONTINUE
C
                  END IF
              END IF
  100     CONTINUE
C
*
*        INCREMENT OP COUNT FOR COMPUTING NORM IN LOOP 140
          OPS = OPS + 3*IEND
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  110     L = LOW + 1
C
*
*        INCREMENT OP COUNT FOR LOOP 170
          OPS = OPS + (12* (IGH-LOW+1)+42)* (IGH-L+1)
          DO 150 I = L,IGH
              LL = MIN0(I+1,IGH)
              IF (HI(I,I-1).NE.0.0E0) THEN
                  NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
                  YR = HR(I,I-1)/NORM
                  YI = HI(I,I-1)/NORM
                  HR(I,I-1) = NORM
                  HI(I,I-1) = 0.0E0
C
                  DO 120 J = I,N
                      SI = YR*HI(I,J) - YI*HR(I,J)
                      HR(I,J) = YR*HR(I,J) + YI*HI(I,J)
                      HI(I,J) = SI
  120             CONTINUE
C
                  DO 130 J = 1,LL
                      SI = YR*HI(J,I) + YI*HR(J,I)
                      HR(J,I) = YR*HR(J,I) - YI*HI(J,I)
                      HI(J,I) = SI
  130             CONTINUE
C
                  DO 140 J = LOW,IGH
                      SI = YR*ZI(J,I) + YI*ZR(J,I)
                      ZR(J,I) = YR*ZR(J,I) - YI*ZI(J,I)
                      ZI(J,I) = SI
  140             CONTINUE
              END IF
  150     CONTINUE
C
C*PL*ERROR* Embedded comment after label moved
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  160     DO 170 I = 1,N
              IF (I.LT.LOW .OR. I.GT.IGH) THEN
                  WR(I) = HR(I,I)
                  WI(I) = HI(I,I)
              END IF
  170     CONTINUE
C
          EN = IGH
          TR = 0.0E0
          TI = 0.0E0
          ITN = 30*N
  180     CONTINUE
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
          IF (EN.LT.LOW) THEN
              GO TO 320
          ELSE
              ITS = 0
              ENM1 = EN - 1
  190         CONTINUE
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
              DO 200 LL = LOW,EN
                  L = EN + LOW - LL
                  IF (L.EQ.LOW) THEN
                      GO TO 210
                  ELSE
                      TST1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1)) +
     +                       ABS(HR(L,L)) + ABS(HI(L,L))
*         TST2 = TST1 + ABS(HR(L,L-1))
*         IF (TST2 .EQ. TST1) GO TO 300
                      TST2 = ABS(HR(L,L-1))
                      IF (TST2.LE.MIN(ULP*TST1,SMALL)) GO TO 210
                  END IF
  200         CONTINUE
C     .......... FORM SHIFT ..........
  210         CONTINUE
*
*        INCREMENT OP COUNT FOR CONVERGENCE TEST
              OPS = OPS + 4* (EN-L+1)
              IF (L.NE.EN) THEN
                  IF (ITN.EQ.0) THEN
                      GO TO 310
                  ELSE
                      IF (ITS.EQ.10 .OR. ITS.EQ.20) THEN
C     .......... FORM EXCEPTIONAL SHIFT ..........
                          SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
                          SI = 0.0E0
                      ELSE
*
*        INCREMENT OPCOUNT FOR FOMING SHIFT
                          OPST = OPST + 58
                          SR = HR(EN,EN)
                          SI = HI(EN,EN)
                          XR = HR(ENM1,EN)*HR(EN,ENM1)
                          XI = HI(ENM1,EN)*HR(EN,ENM1)
                          IF (XR.NE.0.0E0 .OR. XI.NE.0.0E0) THEN
                              YR = (HR(ENM1,ENM1)-SR)/2.0E0
                              YI = (HI(ENM1,ENM1)-SI)/2.0E0
                              CALL CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,
     +                                    ZZR,ZZI)
                              IF (YR*ZZR+YI*ZZI.LT.0.0E0) THEN
                                  ZZR = -ZZR
                                  ZZI = -ZZI
                              END IF
                              CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
                              SR = SR - XR
                              SI = SI - XI
                          END IF
                      END IF
C
                      DO 220 I = LOW,EN
                          HR(I,I) = HR(I,I) - SR
                          HI(I,I) = HI(I,I) - SI
  220                 CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 360
                      OPS = OPS + 2* (EN-LOW+1)
C
                      TR = TR + SR
                      TI = TI + SI
                      ITS = ITS + 1
                      ITN = ITN - 1
*
*       UPDATE ITERATION NUMBER
                      ITCNT = 30*N - ITN
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
                      LP1 = L + 1
C
*
*        INCREMENT OPCOUNT FOR REDUCING TO TRIANGULAR, LOOP 500
                      OPS = OPS + (EN-LP1+1)* (61+10* (EN-LP1))
                      DO 240 I = LP1,EN
                          SR = HR(I,I-1)
                          HR(I,I-1) = 0.0E0
                          NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),
     +                           SR)
                          XR = HR(I-1,I-1)/NORM
                          WR(I-1) = XR
                          XI = HI(I-1,I-1)/NORM
                          WI(I-1) = XI
                          HR(I-1,I-1) = NORM
                          HI(I-1,I-1) = 0.0E0
                          HI(I,I-1) = SR/NORM
C
                          DO 230 J = I,N
                              YR = HR(I-1,J)
                              YI = HI(I-1,J)
                              ZZR = HR(I,J)
                              ZZI = HI(I,J)
                              HR(I-1,J) = XR*YR + XI*YI + HI(I,I-1)*ZZR
                              HI(I-1,J) = XR*YI - XI*YR + HI(I,I-1)*ZZI
                              HR(I,J) = XR*ZZR - XI*ZZI - HI(I,I-1)*YR
                              HI(I,J) = XR*ZZI + XI*ZZR - HI(I,I-1)*YI
  230                     CONTINUE
  240                 CONTINUE
C
C
                      SI = HI(EN,EN)
                      IF (SI.NE.0.0E0) THEN
                          NORM = PYTHAG(HR(EN,EN),SI)
                          SR = HR(EN,EN)/NORM
                          SI = SI/NORM
                          HR(EN,EN) = NORM
                          HI(EN,EN) = 0.0E0
*
*        INCREMENT OP COUNT
                          OPST = OPST + 20
                          IF (EN.NE.N) THEN
                              IP1 = EN + 1
C
*
*        INCREMENT OP COUNT FOR LOOP 520
                              OPST = OPST + 6* (N-IP1+1)
                              DO 250 J = IP1,N
                                  YR = HR(EN,J)
                                  YI = HI(EN,J)
                                  HR(EN,J) = SR*YR + SI*YI
                                  HI(EN,J) = SR*YI - SI*YR
  250                         CONTINUE
                          END IF
                      END IF
C     .......... INVERSE OPERATION (COLUMNS) ..........
                      DO 280 J = LP1,EN
                          XR = WR(J-1)
                          XI = WI(J-1)
C
                          DO 260 I = 1,J
                              YR = HR(I,J-1)
                              YI = 0.0E0
                              ZZR = HR(I,J)
                              ZZI = HI(I,J)
                              IF (I.NE.J) THEN
                                  YI = HI(I,J-1)
                                  HI(I,J-1) = XR*YI + XI*YR +
     +                                        HI(J,J-1)*ZZI
                              END IF
                              HR(I,J-1) = XR*YR - XI*YI + HI(J,J-1)*ZZR
                              HR(I,J) = XR*ZZR + XI*ZZI - HI(J,J-1)*YR
                              HI(I,J) = XR*ZZI - XI*ZZR - HI(J,J-1)*YI
  260                     CONTINUE
C
                          DO 270 I = LOW,IGH
                              YR = ZR(I,J-1)
                              YI = ZI(I,J-1)
                              ZZR = ZR(I,J)
                              ZZI = ZI(I,J)
                              ZR(I,J-1) = XR*YR - XI*YI + HI(J,J-1)*ZZR
                              ZI(I,J-1) = XR*YI + XI*YR + HI(J,J-1)*ZZI
                              ZR(I,J) = XR*ZZR + XI*ZZI - HI(J,J-1)*YR
                              ZI(I,J) = XR*ZZI - XI*ZZR - HI(J,J-1)*YI
  270                     CONTINUE
  280                 CONTINUE
C
*
*        INCREMENT OPCOUNT FOR INVERSE OPERATION LOOP 600
                      OPS = OPS + (10* (EN+LP1)+20* (IGH-LOW+1))*
     +                      (EN-LP1+1)
C
                      IF (SI.NE.0.0E0) THEN
C
*
*        INCREMENT OPCOUNT FOR LOOP 630 AND 640
                          OPS = OPS + 6*EN + 6* (IGH-LOW+1)
                          DO 290 I = 1,EN
                              YR = HR(I,EN)
                              YI = HI(I,EN)
                              HR(I,EN) = SR*YR - SI*YI
                              HI(I,EN) = SR*YI + SI*YR
  290                     CONTINUE
C
                          DO 300 I = LOW,IGH
                              YR = ZR(I,EN)
                              YI = ZI(I,EN)
                              ZR(I,EN) = SR*YR - SI*YI
                              ZI(I,EN) = SR*YI + SI*YR
  300                     CONTINUE
C
                      END IF
                      GO TO 190
                  END IF
              END IF
C     .......... A ROOT FOUND ..........
              HR(EN,EN) = HR(EN,EN) + TR
              WR(EN) = HR(EN,EN)
              HI(EN,EN) = HI(EN,EN) + TI
              WI(EN) = HI(EN,EN)
              EN = ENM1
              GO TO 180
          END IF
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
  310     IERR = EN
          GO TO 450
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  320     NORM = 0.0E0
C
*
*        INCREMENT OP COUNT FOR LOOP 720
          OPS = OPS + N* (N+1)/2
          DO 340 I = 1,N
C
              DO 330 J = I,N
                  TR = ABS(HR(I,J)) + ABS(HI(I,J))
                  IF (TR.GT.NORM) NORM = TR
  330         CONTINUE
  340     CONTINUE
C
          IF (N.NE.1 .AND. NORM.NE.0.0E0) THEN
C     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
              DO 390 NN = 2,N
                  EN = N + 2 - NN
                  XR = WR(EN)
                  XI = WI(EN)
                  HR(EN,EN) = 1.0E0
                  HI(EN,EN) = 0.0E0
                  ENM1 = EN - 1
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
*
*        INCREMENT OP COUNT FOR COMPUT YR, .. IN LOOP 780
                  OPS = OPS + 22*ENM1
                  DO 380 II = 1,ENM1
                      I = EN - II
                      ZZR = 0.0E0
                      ZZI = 0.0E0
                      IP1 = I + 1
C
*
*        INCREMENT OP COUNT FOR LOOP 740
                      OPS = OPS + 7* (EN-IP1+1)
                      DO 350 J = IP1,EN
                          ZZR = ZZR + HR(I,J)*HR(J,EN) -
     +                          HI(I,J)*HI(J,EN)
                          ZZI = ZZI + HR(I,J)*HI(J,EN) +
     +                          HI(I,J)*HR(J,EN)
  350                 CONTINUE
C
                      YR = XR - WR(I)
                      YI = XI - WI(I)
                      IF (YR.EQ.0.0E0 .AND. YI.EQ.0.0E0) THEN
                          TST1 = NORM
                          YR = TST1
  360                     CONTINUE
                          YR = 0.01E0*YR
                          TST2 = NORM + YR
                          IF (TST2.GT.TST1) GO TO 360
                      END IF
                      CALL CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
*
*        INCREMENT OP COUNT FOR CDIV
                      OPST = OPST + 16
C     .......... OVERFLOW CONTROL ..........
                      TR = ABS(HR(I,EN)) + ABS(HI(I,EN))
                      IF (TR.NE.0.0E0) THEN
                          TST1 = TR
                          TST2 = TST1 + 1.0E0/TST1
                          IF (TST2.LE.TST1) THEN
*
*        INCREMENT OP COUNT FOR LOOP 770
                              OPS = OPS + 2* (EN-I+1)
                              DO 370 J = I,EN
                                  HR(J,EN) = HR(J,EN)/TR
                                  HI(J,EN) = HI(J,EN)/TR
  370                         CONTINUE
                          END IF
                      END IF
  380             CONTINUE
C
  390         CONTINUE
C
C     .......... END BACKSUBSTITUTION ..........
              ENM1 = N - 1
C     .......... VECTORS OF ISOLATED ROOTS ..........
              DO 410 I = 1,N
                  IF (I.LT.LOW .OR. I.GT.IGH) THEN
                      IP1 = I + 1
C
                      DO 400 J = I,N
                          ZR(I,J) = HR(I,J)
                          ZI(I,J) = HI(I,J)
  400                 CONTINUE
                  END IF
  410         CONTINUE
C
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
              DO 440 JJ = LOW,N
                  J = N + LOW - JJ
                  M = MIN0(J,IGH)
C
*
*        INCREMENT OP COUNT FOR LOOP 880
                  OPS = OPS + 8* (M-LOW+1)* (IGH-LOW+1)
                  DO 430 I = LOW,IGH
                      ZZR = 0.0E0
                      ZZI = 0.0E0
C
                      DO 420 K = LOW,M
                          ZZR = ZZR + ZR(I,K)*HR(K,J) - ZI(I,K)*HI(K,J)
                          ZZI = ZZI + ZR(I,K)*HI(K,J) + ZI(I,K)*HR(K,J)
  420                 CONTINUE
C
                      ZR(I,J) = ZZR
                      ZI(I,J) = ZZI
  430             CONTINUE
  440         CONTINUE
C
          END IF
*
*     COMPUTE FINAL OP COUNT
  450     OPS = OPS + OPST
      END IF
      END
      SUBROUTINE CORTH(NM,N,LOW,IGH,AR,AI,ORTR,ORTI)
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
      COMMON /PYTHOP/OPST
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS,OPST
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE ORTHES, NUM. MATH. 12, 349-368(1968)
C     BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.  INFORMATION
C          ABOUT THE UNITARY TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLES UNDER THE
C          HESSENBERG MATRIX.
C
C        ORTR AND ORTI CONTAIN FURTHER INFORMATION ABOUT THE
C          TRANSFORMATIONS.  ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER IGH,LOW,N,NM
C     ..
C     .. Array Arguments ..
      REAL AI(NM,N),AR(NM,N),ORTI(IGH),ORTR(IGH)
C     ..
C     .. Local Scalars ..
      REAL F,FI,FR,G,H,SCALE
      INTEGER I,II,J,JJ,KP1,LA,M,MP
C     ..
C     .. External Functions ..
      REAL PYTHAG
      EXTERNAL PYTHAG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
      IF (N.GT.0) THEN
***
*     INITIALIZE
          OPST = 0
***
          LA = IGH - 1
          KP1 = LOW + 1
          IF (LA.GE.KP1) THEN
C
              DO 90 M = KP1,LA
                  H = 0.0E0
                  ORTR(M) = 0.0E0
                  ORTI(M) = 0.0E0
                  SCALE = 0.0E0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
                  DO 10 I = M,IGH
                      SCALE = SCALE + ABS(AR(I,M-1)) + ABS(AI(I,M-1))
   10             CONTINUE
***
*        INCREMENT OPCOUNT FOR LOOP 90
                  OPS = OPS + 2* (IGH-M+1)
***
C
                  IF (SCALE.NE.0.0E0) THEN
                      MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
                      DO 20 II = M,IGH
                          I = MP - II
                          ORTR(I) = AR(I,M-1)/SCALE
                          ORTI(I) = AI(I,M-1)/SCALE
                          H = H + ORTR(I)*ORTR(I) + ORTI(I)*ORTI(I)
   20                 CONTINUE
***
*        INCREMENT OP COUNT FOR LOOP 100 AND SQRT
                      OPS = OPS + 6* (IGH-M+1) + 1
***
C
                      G = SQRT(H)
                      F = PYTHAG(ORTR(M),ORTI(M))
                      IF (F.EQ.0.0E0) THEN
C
                          ORTR(M) = G
                          AR(M,M-1) = SCALE
                      ELSE
                          H = H + F*G
                          G = G/F
                          ORTR(M) = (1.0E0+G)*ORTR(M)
                          ORTI(M) = (1.0E0+G)*ORTI(M)
                          OPST = OPST + 7
                      END IF
C     .......... FORM (I-(U*UT)/H) * A ..........
                      DO 50 J = M,N
                          FR = 0.0E0
                          FI = 0.0E0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
                          DO 30 II = M,IGH
                              I = MP - II
                              FR = FR + ORTR(I)*AR(I,J) +
     +                             ORTI(I)*AI(I,J)
                              FI = FI + ORTR(I)*AI(I,J) -
     +                             ORTI(I)*AR(I,J)
   30                     CONTINUE
C
                          FR = FR/H
                          FI = FI/H
C
                          DO 40 I = M,IGH
                              AR(I,J) = AR(I,J) - FR*ORTR(I) +
     +                                  FI*ORTI(I)
                              AI(I,J) = AI(I,J) - FR*ORTI(I) -
     +                                  FI*ORTR(I)
   40                     CONTINUE
   50                 CONTINUE
C
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
                      DO 80 I = 1,IGH
                          FR = 0.0E0
                          FI = 0.0E0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
                          DO 60 JJ = M,IGH
                              J = MP - JJ
                              FR = FR + ORTR(J)*AR(I,J) -
     +                             ORTI(J)*AI(I,J)
                              FI = FI + ORTR(J)*AI(I,J) +
     +                             ORTI(J)*AR(I,J)
   60                     CONTINUE
C
                          FR = FR/H
                          FI = FI/H
C
                          DO 70 J = M,IGH
                              AR(I,J) = AR(I,J) - FR*ORTR(J) -
     +                                  FI*ORTI(J)
                              AI(I,J) = AI(I,J) + FR*ORTI(J) -
     +                                  FI*ORTR(J)
   70                     CONTINUE
   80                 CONTINUE
C
***
*        INCREMENT OP COUNT FOR LOOPS 130 AND 160
                      OPS = OPS + (IGH+N-M+1)* ((IGH-M+1)*16+2)
                      OPST = OPST + 4
***
C
                      ORTR(M) = SCALE*ORTR(M)
                      ORTI(M) = SCALE*ORTI(M)
                      AR(M,M-1) = -G*AR(M,M-1)
                      AI(M,M-1) = -G*AI(M,M-1)
                  END IF
   90         CONTINUE
              OPS = OPS + OPST
          END IF
          RETURN
      END IF
*$st$ Unreachable comments ...
C
      END
      SUBROUTINE CSROOT(XR,XI,YR,YI)
C
C     (YR,YI) = COMPLEX SQRT(XR,XI)
C     BRANCH CHOSEN SO THAT YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI)
C

C     .. Scalar Arguments ..
      REAL XI,XR,YI,YR
C     ..
C     .. Local Scalars ..
      REAL S,TI,TR
C     ..
C     .. External Functions ..
      REAL PYTHAG
      EXTERNAL PYTHAG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
      TR = XR
      TI = XI
      S = SQRT(0.5E0* (PYTHAG(TR,TI)+ABS(TR)))
      IF (TR.GE.0.0E0) YR = S
      IF (TI.LT.0.0E0) S = -S
      IF (TR.LE.0.0E0) YI = S
      IF (TR.LT.0.0E0) YR = 0.5E0* (TI/YI)
      IF (TR.GT.0.0E0) YI = 0.5E0* (TI/YR)
      END
      SUBROUTINE HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT.
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED.
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRBAK1, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX HERMITIAN
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     REAL SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  HTRIDI.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  HTRIDI  IN THEIR
C          FULL LOWER TRIANGLES EXCEPT FOR THE DIAGONAL OF AR.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        ZR CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     NOTE THAT THE LAST COMPONENT OF EACH RETURNED VECTOR
C     IS REAL AND THAT VECTOR EUCLIDEAN NORMS ARE PRESERVED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER M,N,NM
C     ..
C     .. Array Arguments ..
      REAL AI(NM,N),AR(NM,N),TAU(2,N),ZI(NM,M),ZR(NM,M)
C     ..
C     .. Local Scalars ..
      REAL H,S,SI
      INTEGER I,J,K,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,REAL
C     ..
      IF (M.NE.0) THEN
*
          OPS = OPS + MAX(0.0E0,8*M*REAL(N)**2-2*M*REAL(N)-4*M)
*
C     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
C                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
C                TRIDIAGONAL MATRIX. ..........
          DO 20 K = 1,N
C
              DO 10 J = 1,M
                  ZI(K,J) = -ZR(K,J)*TAU(2,K)
                  ZR(K,J) = ZR(K,J)*TAU(1,K)
   10         CONTINUE
   20     CONTINUE
C
          IF (N.NE.1) THEN
C     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
              DO 60 I = 2,N
                  L = I - 1
                  H = AI(I,I)
                  IF (H.NE.0.0E0) THEN
C
                      DO 50 J = 1,M
                          S = 0.0E0
                          SI = 0.0E0
C
                          DO 30 K = 1,L
                              S = S + AR(I,K)*ZR(K,J) - AI(I,K)*ZI(K,J)
                              SI = SI + AR(I,K)*ZI(K,J) +
     +                             AI(I,K)*ZR(K,J)
   30                     CONTINUE
C     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
                          S = (S/H)/H
                          SI = (SI/H)/H
C
                          DO 40 K = 1,L
                              ZR(K,J) = ZR(K,J) - S*AR(I,K) - SI*AI(I,K)
                              ZI(K,J) = ZI(K,J) - SI*AR(I,K) + S*AI(I,K)
   40                     CONTINUE
   50                 CONTINUE
C
                  END IF
   60         CONTINUE
C
          END IF
      END IF
      RETURN
*$st$ Unreachable comments ...
C
      END
      SUBROUTINE HTRIDI(NM,N,AR,AI,D,E,E2,TAU)
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT.
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED.
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRED1, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A COMPLEX HERMITIAN MATRIX
C     TO A REAL SYMMETRIC TRIDIAGONAL MATRIX USING
C     UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX HERMITIAN INPUT MATRIX.
C          ONLY THE LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION IN THEIR FULL LOWER
C          TRIANGLES.  THEIR STRICT UPPER TRIANGLES AND THE
C          DIAGONAL OF AR ARE UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
*
C     .. Scalar Arguments ..
      INTEGER N,NM
C     ..
C     .. Array Arguments ..
      REAL AI(NM,N),AR(NM,N),D(N),E(N),E2(N),TAU(2,N)
C     ..
C     .. Local Scalars ..
      REAL F,FI,G,GI,H,HH,SCALE,SI
      INTEGER I,II,J,JP1,K,L
C     ..
C     .. External Functions ..
      REAL PYTHAG
      EXTERNAL PYTHAG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,REAL,SQRT
C     ..
      OPS = OPS + MAX(0.0E0, (16.0E0/3.0E0)*REAL(N)**3+3*REAL(N)**2+
     +      (56.0E0/3.0E0)*N-61)
*
      TAU(1,N) = 1.0E0
      TAU(2,N) = 0.0E0
C
      DO 10 I = 1,N
          D(I) = AR(I,I)
   10 CONTINUE
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 120 II = 1,N
          I = N + 1 - II
          L = I - 1
          H = 0.0E0
          SCALE = 0.0E0
          IF (L.GE.1) THEN
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
              DO 20 K = 1,L
                  SCALE = SCALE + ABS(AR(I,K)) + ABS(AI(I,K))
   20         CONTINUE
C
              IF (SCALE.NE.0.0E0) THEN
C
                  DO 30 K = 1,L
                      AR(I,K) = AR(I,K)/SCALE
                      AI(I,K) = AI(I,K)/SCALE
                      H = H + AR(I,K)*AR(I,K) + AI(I,K)*AI(I,K)
   30             CONTINUE
C
                  E2(I) = SCALE*SCALE*H
                  G = SQRT(H)
                  E(I) = SCALE*G
                  F = PYTHAG(AR(I,L),AI(I,L))
C     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
                  IF (F.EQ.0.0E0) THEN
                      TAU(1,L) = -TAU(1,I)
                      SI = TAU(2,I)
                      AR(I,L) = G
                  ELSE
                      TAU(1,L) = (AI(I,L)*TAU(2,I)-AR(I,L)*TAU(1,I))/F
                      SI = (AR(I,L)*TAU(2,I)+AI(I,L)*TAU(1,I))/F
                      H = H + F*G
                      G = 1.0E0 + G/F
                      AR(I,L) = G*AR(I,L)
                      AI(I,L) = G*AI(I,L)
                      IF (L.EQ.1) GO TO 90
                  END IF
                  F = 0.0E0
C
                  DO 60 J = 1,L
                      G = 0.0E0
                      GI = 0.0E0
C     .......... FORM ELEMENT OF A*U ..........
                      DO 40 K = 1,J
                          G = G + AR(J,K)*AR(I,K) + AI(J,K)*AI(I,K)
                          GI = GI - AR(J,K)*AI(I,K) + AI(J,K)*AR(I,K)
   40                 CONTINUE
C
                      JP1 = J + 1
                      IF (L.GE.JP1) THEN
C
                          DO 50 K = JP1,L
                              G = G + AR(K,J)*AR(I,K) - AI(K,J)*AI(I,K)
                              GI = GI - AR(K,J)*AI(I,K) -
     +                             AI(K,J)*AR(I,K)
   50                     CONTINUE
                      END IF
C     .......... FORM ELEMENT OF P ..........
                      E(J) = G/H
                      TAU(2,J) = GI/H
                      F = F + E(J)*AR(I,J) - TAU(2,J)*AI(I,J)
   60             CONTINUE
C
                  HH = F/ (H+H)
C     .......... FORM REDUCED A ..........
                  DO 80 J = 1,L
                      F = AR(I,J)
                      G = E(J) - HH*F
                      E(J) = G
                      FI = -AI(I,J)
                      GI = TAU(2,J) - HH*FI
                      TAU(2,J) = -GI
C
                      DO 70 K = 1,J
                          AR(J,K) = AR(J,K) - F*E(K) - G*AR(I,K) +
     +                              FI*TAU(2,K) + GI*AI(I,K)
                          AI(J,K) = AI(J,K) - F*TAU(2,K) - G*AI(I,K) -
     +                              FI*E(K) - GI*AR(I,K)
   70                 CONTINUE
   80             CONTINUE
C*PL*ERROR* Embedded comment after label moved
C
   90             DO 100 K = 1,L
                      AR(I,K) = SCALE*AR(I,K)
                      AI(I,K) = SCALE*AI(I,K)
  100             CONTINUE
C
                  TAU(2,L) = -SI
                  GO TO 110
              ELSE
                  TAU(1,L) = 1.0E0
                  TAU(2,L) = 0.0E0
              END IF
          END IF
          E(I) = 0.0E0
          E2(I) = 0.0E0
  110     HH = D(I)
          D(I) = AR(I,I)
          AR(I,I) = HH
          AI(I,I) = SCALE*SQRT(H)
  120 CONTINUE
C
      END
      SUBROUTINE IMTQL1(N,D,E,IERR)
*
*     EISPACK ROUTINE
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN SSTEQR.
*
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE CONTRIBUTIONS TO OPS FROM
*     FUNCTION PYTHAG.  IT IS PASSED TO AND FROM PYTHAG
*     THROUGH COMMON BLOCK PYTHOP.
C     .. Common blocks ..
*
      COMMON /LATIME/OPS,ITCNT
      COMMON /PYTHOP/OPST
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS,OPST
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL1,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 40 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER IERR,N
C     ..
C     .. Array Arguments ..
      REAL D(N),E(N)
C     ..
C     .. Local Scalars ..
      REAL B,C,EPS,F,G,P,R,S,TST
      INTEGER I,II,J,L,M,MML
C     ..
C     .. External Functions ..
      REAL PYTHAG,SLAMCH
      EXTERNAL PYTHAG,SLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SIGN
C     ..
      IERR = 0
      IF (N.NE.1) THEN
*
*        INITIALIZE ITERATION COUNT AND OPST
          ITCNT = 0
          OPST = 0
*
*     DETERMINE THE UNIT ROUNDOFF FOR THIS ENVIRONMENT.
*
          EPS = SLAMCH('EPSILON')
C
          DO 10 I = 2,N
              E(I-1) = E(I)
   10     CONTINUE
C
          E(N) = 0.0E0
C
          DO 90 L = 1,N
              J = 0
   20         CONTINUE
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
              DO 30 M = L,N
                  IF (M.EQ.N) THEN
                      GO TO 40
                  ELSE
                      TST = ABS(E(M))
                      IF (TST.LE.EPS* (ABS(D(M))+ABS(D(M+1))))
     +                    GO TO 40
                  END IF
   30         CONTINUE
*            TST1 = ABS(D(M)) + ABS(D(M+1))
*            TST2 = TST1 + ABS(E(M))
*            IF (TST2 .EQ. TST1) GO TO 120
C
   40         P = D(L)
*
*        INCREMENT OPCOUNT FOR FINDING SMALL SUBDIAGONAL ELEMENT.
              OPS = OPS + 2* (MIN(M,N-1)-L+1)
              IF (M.NE.L) THEN
                  IF (J.EQ.40) THEN
                      GO TO 100
                  ELSE
                      J = J + 1
C     .......... FORM SHIFT ..........
                      G = (D(L+1)-P)/ (2.0E0*E(L))
                      R = PYTHAG(G,1.0E0)
                      G = D(M) - P + E(L)/ (G+SIGN(R,G))
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT.
                      OPS = OPS + 7
                      S = 1.0E0
                      C = 1.0E0
                      P = 0.0E0
                      MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
                      DO 50 II = 1,MML
                          I = M - II
                          F = S*E(I)
                          B = C*E(I)
                          R = PYTHAG(F,G)
                          E(I+1) = R
                          IF (R.EQ.0.0E0) THEN
                              GO TO 60
                          ELSE
                              S = F/R
                              C = G/R
                              G = D(I+1) - P
                              R = (D(I)-G)*S + 2.0E0*C*B
                              P = S*R
                              D(I+1) = G + P
                              G = C*R - B
                          END IF
   50                 CONTINUE
C
                      D(L) = D(L) - P
                      E(L) = G
                      E(M) = 0.0E0
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
                      OPS = OPS + MML*14 + 1
*
*        INCREMENT ITERATION COUNTER
                      ITCNT = ITCNT + 1
                      GO TO 20
C     .......... RECOVER FROM UNDERFLOW ..........
   60                 D(I+1) = D(I+1) - P
                      E(M) = 0.0E0
*
*        INCREMENT OPCOUNT FOR INNER LOOP, WHEN UNDERFLOW OCCURS.
                      OPS = OPS + 2 + (II-1)*14 + 1
                      GO TO 20
                  END IF
              END IF
C     .......... ORDER EIGENVALUES ..........
              IF (L.NE.1) THEN
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
                  DO 70 II = 2,L
                      I = L + 2 - II
                      IF (P.GE.D(I-1)) THEN
                          GO TO 80
                      ELSE
                          D(I) = D(I-1)
                      END IF
   70             CONTINUE
              END IF
C
              I = 1
   80         D(I) = P
   90     CONTINUE
C
          GO TO 110
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 40 ITERATIONS ..........
  100     IERR = L
      END IF
*
*     COMPUTE FINAL OP COUNT
  110 OPS = OPS + OPST
      END
      SUBROUTINE IMTQL2(NM,N,D,E,Z,IERR)
*
*     EISPACK ROUTINE.  MODIFIED FOR COMPARISON WITH LAPACK.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN SSTEQR.
*
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE CONTRIBUTIONS TO OPS FROM
*     FUNCTION PYTHAG.  IT IS PASSED TO AND FROM PYTHAG
*     THROUGH COMMON BLOCK PYTHOP.
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
      COMMON /PYTHOP/OPST
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS,OPST
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C
C        E HAS BEEN DESTROYED.
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 40 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER IERR,N,NM
C     ..
C     .. Array Arguments ..
      REAL D(N),E(N),Z(NM,N)
C     ..
C     .. Local Scalars ..
      REAL B,C,EPS,F,G,P,R,S,TST
      INTEGER I,II,J,K,L,M,MML
C     ..
C     .. External Functions ..
      REAL PYTHAG,SLAMCH
      EXTERNAL PYTHAG,SLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SIGN
C     ..
      IERR = 0
      IF (N.NE.1) THEN
*
*        INITIALIZE ITERATION COUNT AND OPST
          ITCNT = 0
          OPST = 0
*
*     DETERMINE UNIT ROUNDOFF FOR THIS MACHINE.
          EPS = SLAMCH('EPSILON')
C
          DO 10 I = 2,N
              E(I-1) = E(I)
   10     CONTINUE
C
          E(N) = 0.0E0
C
          DO 80 L = 1,N
              J = 0
   20         CONTINUE
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
              DO 30 M = L,N
                  IF (M.EQ.N) THEN
                      GO TO 40
                  ELSE
*            TST1 = ABS(D(M)) + ABS(D(M+1))
*            TST2 = TST1 + ABS(E(M))
*            IF (TST2 .EQ. TST1) GO TO 120
                      TST = ABS(E(M))
                      IF (TST.LE.EPS* (ABS(D(M))+ABS(D(M+1))))
     +                    GO TO 40
                  END IF
   30         CONTINUE
C
   40         P = D(L)
*
*        INCREMENT OPCOUNT FOR FINDING SMALL SUBDIAGONAL ELEMENT.
              OPS = OPS + 2* (MIN(M,N)-L+1)
              IF (M.NE.L) THEN
                  IF (J.EQ.40) THEN
                      GO TO 120
                  ELSE
                      J = J + 1
C     .......... FORM SHIFT ..........
                      G = (D(L+1)-P)/ (2.0E0*E(L))
                      R = PYTHAG(G,1.0E0)
                      G = D(M) - P + E(L)/ (G+SIGN(R,G))
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT.
                      OPS = OPS + 7
                      S = 1.0E0
                      C = 1.0E0
                      P = 0.0E0
                      MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
                      DO 60 II = 1,MML
                          I = M - II
                          F = S*E(I)
                          B = C*E(I)
                          R = PYTHAG(F,G)
                          E(I+1) = R
                          IF (R.EQ.0.0E0) THEN
                              GO TO 70
                          ELSE
                              S = F/R
                              C = G/R
                              G = D(I+1) - P
                              R = (D(I)-G)*S + 2.0E0*C*B
                              P = S*R
                              D(I+1) = G + P
                              G = C*R - B
C     .......... FORM VECTOR ..........
                              DO 50 K = 1,N
                                  F = Z(K,I+1)
                                  Z(K,I+1) = S*Z(K,I) + C*F
                                  Z(K,I) = C*Z(K,I) - S*F
   50                         CONTINUE
                          END IF
   60                 CONTINUE
C
C
                      D(L) = D(L) - P
                      E(L) = G
                      E(M) = 0.0E0
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
                      OPS = OPS + MML* (14+6*N) + 1
*
*        INCREMENT ITERATION COUNTER
                      ITCNT = ITCNT + 1
                      GO TO 20
C     .......... RECOVER FROM UNDERFLOW ..........
   70                 D(I+1) = D(I+1) - P
                      E(M) = 0.0E0
*
*        INCREMENT OPCOUNT FOR INNER LOOP, WHEN UNDERFLOW OCCURS.
                      OPS = OPS + 2 + (II-1)* (14+6*N) + 1
                      GO TO 20
                  END IF
              END IF
   80     CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
          DO 110 II = 2,N
              I = II - 1
              K = I
              P = D(I)
C
              DO 90 J = II,N
                  IF (D(J).LT.P) THEN
                      K = J
                      P = D(J)
                  END IF
   90         CONTINUE
C
              IF (K.NE.I) THEN
                  D(K) = D(I)
                  D(I) = P
C
                  DO 100 J = 1,N
                      P = Z(J,I)
                      Z(J,I) = Z(J,K)
                      Z(J,K) = P
  100             CONTINUE
              END IF
  110     CONTINUE
C
C
          GO TO 130
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 40 ITERATIONS ..........
  120     IERR = L
      END IF
*
*     COMPUTE FINAL OP COUNT
  130 OPS = OPS + OPST
      END
      REAL FUNCTION PYTHAG(A,B)
C
C     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT
*     OPST IS ONLY INCREMENTED HERE
C     .. Common blocks ..
      COMMON /PYTHOP/OPST
C     ..
*     .. SCALARS IN COMMON

C     .. Scalar Arguments ..
      REAL A,B
C     ..
C     .. Scalars in Common ..
      REAL OPST
C     ..
C     .. Local Scalars ..
      REAL P,R,S,T,U
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AMAX1,AMIN1
C     ..
      P = AMAX1(ABS(A),ABS(B))
      IF (P.NE.0.0E0) THEN
          R = (AMIN1(ABS(A),ABS(B))/P)**2
*
*     INCREMENT OPST
          OPST = OPST + 2
   10     CONTINUE
          T = 4.0E0 + R
          IF (T.NE.4.0E0) THEN
              S = R/T
              U = 1.0E0 + 2.0E0*S
              P = U*P
              R = (S/U)**2*R
*
*        INCREMENT OPST
              OPST = OPST + 8
              GO TO 10
          END IF
      END IF
      PYTHAG = P
      END
      REAL FUNCTION EPSLON(X)
C
C     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C
C
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT
C            NUMBERS IS NOT A POWER OF THREE.
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
C            THE ACCURACY USED IN FLOATING POINT VARIABLES
C            THAT ARE STORED IN MEMORY.
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
C     ASSUMPTION 2.
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
C            C  IS NOT EXACTLY EQUAL TO ONE,
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM
C                 THE NEXT LARGER FLOATING POINT NUMBER.
C     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
C     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
C
C     THIS VERSION DATED 4/6/83.
C
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Local Scalars ..
      REAL A,B,C,EPS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      A = 4.0E0/3.0E0
   10 CONTINUE
      B = A - 1.0E0
      C = B + B + B
      EPS = ABS(C-1.0E0)
      IF (EPS.EQ.0.0E0) GO TO 10
      EPSLON = EPS*ABS(X)
      END
      SUBROUTINE BISECT(N,EPS1,D,E,E2,LB,UB,MM,M,W,IND,IERR,RV4,RV5)
*
*     EISPACK ROUTINE.
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN SSTEBZ.
*
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE BISECTION TECHNIQUE
C     IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX WHICH LIE IN A SPECIFIED INTERVAL,
C     USING BISECTION.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C        LB AND UB DEFINE THE INTERVAL TO BE SEARCHED FOR EIGENVALUES.
C          IF LB IS NOT LESS THAN UB, NO EIGENVALUES WILL BE FOUND.
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          EIGENVALUES IN THE INTERVAL.  WARNING. IF MORE THAN
C          MM EIGENVALUES ARE DETERMINED TO LIE IN THE INTERVAL,
C          AN ERROR RETURN IS MADE WITH NO EIGENVALUES FOUND.
C
C     ON OUTPUT
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE.
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        M IS THE NUMBER OF EIGENVALUES DETERMINED TO LIE IN (LB,UB).
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING ORDER.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF M EXCEEDS MM.
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     THE ALGOL PROCEDURE STURMCNT CONTAINED IN TRISTURM
C     APPEARS IN BISECT IN-LINE.
C
C     NOTE THAT SUBROUTINE TQL1 OR IMTQL1 IS GENERALLY FASTER THAN
C     BISECT, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
*        INITIALIZE ITERATION COUNT.
C     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E0)
      REAL RELFAC
      PARAMETER (RELFAC=2.0E0)
C     ..
C     .. Scalar Arguments ..
      REAL EPS1,LB,UB
      INTEGER IERR,M,MM,N
C     ..
C     .. Array Arguments ..
      REAL D(N),E(N),E2(N),RV4(N),RV5(N),W(MM)
      INTEGER IND(MM)
C     ..
C     .. Local Scalars ..
      REAL ATOLI,PIVMIN,RTOLI,SAFEMN,T1,T2,TMP1,TMP2,TNORM,U,ULP,V,X0,
     +     X1,XU
      INTEGER I,II,ISTURM,J,K,L,M1,M2,P,Q,R,S,TAG
C     ..
C     .. External Functions ..
      REAL EPSLON,SLAMCH
      EXTERNAL EPSLON,SLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AMAX1,AMIN1,MAX,MIN
C     ..
      ITCNT = 0
      SAFEMN = SLAMCH('S')
      ULP = SLAMCH('E')*SLAMCH('B')
      RTOLI = ULP*RELFAC
      IERR = 0
      TAG = 0
      T1 = LB
      T2 = UB
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
      DO 30 I = 1,N
          IF (I.EQ.1) GO TO 10
CCC         TST1 = ABS(D(I)) + ABS(D(I-1))
CCC         TST2 = TST1 + ABS(E(I))
CCC         IF (TST2 .GT. TST1) GO TO 40
          TMP1 = E(I)**2
          IF (ABS(D(I)*D(I-1))*ULP**2+SAFEMN.LE.TMP1) GO TO 20
   10     E2(I) = 0.0E0
   20     CONTINUE
   30 CONTINUE
*           INCREMENT OPCOUNT FOR DETERMINING IF MATRIX SPLITS.
      OPS = OPS + 5* (N-1)
C
C                COMPUTE QUANTITIES NEEDED FOR CONVERGENCE TEST.
      TMP1 = D(1) - ABS(E(2))
      TMP2 = D(1) + ABS(E(2))
      PIVMIN = ONE
      DO 40 I = 2,N - 1
          TMP1 = MIN(TMP1,D(I)-ABS(E(I))-ABS(E(I+1)))
          TMP2 = MAX(TMP2,D(I)+ABS(E(I))+ABS(E(I+1)))
          PIVMIN = MAX(PIVMIN,E(I)**2)
   40 CONTINUE
      TMP1 = MIN(TMP1,D(N)-ABS(E(N)))
      TMP2 = MAX(TMP2,D(N)+ABS(E(N)))
      PIVMIN = MAX(PIVMIN,E(N)**2)
      PIVMIN = PIVMIN*SAFEMN
      TNORM = MAX(ABS(TMP1),ABS(TMP2))
      ATOLI = ULP*TNORM
*        INCREMENT OPCOUNT FOR COMPUTING THESE QUANTITIES.
      OPS = OPS + 4* (N-1)
C
C     .......... DETERMINE THE NUMBER OF EIGENVALUES
C                IN THE INTERVAL ..........
      P = 1
      Q = N
      X1 = UB
      ISTURM = 1
      GO TO 200
   50 M = S
      X1 = LB
      ISTURM = 2
      GO TO 200
   60 M = M - S
      IF (M.GT.MM) GO TO 350
      Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
   70 IF (R.EQ.M) GO TO 360
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0E0
C
      DO 90 Q = P,N
          X1 = U
          U = 0.0E0
          V = 0.0E0
          IF (Q.EQ.N) GO TO 80
          U = ABS(E(Q+1))
          V = E2(Q+1)
   80     XU = AMIN1(D(Q)- (X1+U),XU)
          X0 = AMAX1(D(Q)+ (X1+U),X0)
          IF (V.EQ.0.0E0) GO TO 100
   90 CONTINUE
*        INCREMENT OPCOUNT FOR REFINING INTERVAL.
      OPS = OPS + (N-P+1)*2
C
  100 X1 = EPSLON(AMAX1(ABS(XU),ABS(X0)))
      IF (EPS1.LE.0.0E0) EPS1 = -X1
      IF (P.NE.Q) GO TO 110
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1.GT.D(P) .OR. D(P).GE.T2) GO TO 340
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 280
  110 X1 = X1* (Q-P+1)
      LB = AMAX1(T1,XU-X1)
      UB = AMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 200
  120 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 200
  130 M2 = S
      IF (M1.GT.M2) GO TO 340
C     .......... FIND ROOTS BY BISECTION ..........
      X0 = UB
      ISTURM = 5
C
      DO 140 I = M1,M2
          RV5(I) = UB
          RV4(I) = LB
  140 CONTINUE
C     .......... LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
      K = M2
  150 XU = LB
C     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
      DO 170 II = M1,K
          I = M1 + K - II
          IF (XU.GE.RV4(I)) GO TO 160
          XU = RV4(I)
          GO TO 180
  160     CONTINUE
  170 CONTINUE
C
  180 IF (X0.GT.RV5(K)) X0 = RV5(K)
C     .......... NEXT BISECTION STEP ..........
  190 X1 = (XU+X0)*0.5E0
CCC         IF ((X0 - XU) .LE. ABS(EPS1)) GO TO 420
CCC         TST1 = 2.0E0 * (ABS(XU) + ABS(X0))
CCC         TST2 = TST1 + (X0 - XU)
CCC         IF (TST2 .EQ. TST1) GO TO 420
      TMP1 = ABS(X0-XU)
      TMP2 = MAX(ABS(X0),ABS(XU))
      IF (TMP1.LT.MAX(ATOLI,PIVMIN,RTOLI*TMP2)) GO TO 270
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  200 S = P - 1
      U = 1.0E0
C
      DO 230 I = P,Q
          IF (U.NE.0.0E0) GO TO 210
          V = ABS(E(I))/EPSLON(1.0E0)
          IF (E2(I).EQ.0.0E0) V = 0.0E0
          GO TO 220
  210     V = E2(I)/U
  220     U = D(I) - X1 - V
          IF (U.LT.0.0E0) S = S + 1
  230 CONTINUE
*           INCREMENT OPCOUNT FOR STURM SEQUENCE.
      OPS = OPS + (Q-P+1)*3
*           INCREMENT ITERATION COUNTER.
      ITCNT = ITCNT + 1
C
      GO TO (50,60,120,130,240) ISTURM
C     .......... REFINE INTERVALS ..........
  240 IF (S.GE.K) GO TO 260
      XU = X1
      IF (S.GE.M1) GO TO 250
      RV4(M1) = X1
      GO TO 190
  250 RV4(S+1) = X1
      IF (RV5(S).GT.X1) RV5(S) = X1
      GO TO 190
  260 X0 = X1
      GO TO 190
C     .......... K-TH EIGENVALUE FOUND ..........
  270 RV5(K) = X1
      K = K - 1
      IF (K.GE.M1) GO TO 150
C     .......... ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ..........
  280 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 330 L = 1,R
          IF (J.GT.S) GO TO 300
          IF (K.GT.M2) GO TO 340
          IF (RV5(K).GE.W(L)) GO TO 310
C
          DO 290 II = J,S
              I = L + S - II
              W(I+1) = W(I)
              IND(I+1) = IND(I)
  290     CONTINUE
C
  300     W(L) = RV5(K)
          IND(L) = TAG
          K = K + 1
          GO TO 320
  310     J = J + 1
  320     CONTINUE
  330 CONTINUE
C
  340 IF (Q.LT.N) GO TO 70
      GO TO 360
C     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF
C                EIGENVALUES IN INTERVAL ..........
  350 IERR = 3*N + 1
  360 LB = T1
      UB = T2
      RETURN
      END
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,IERR,RV1,RV2,RV3,RV4,RV6)
*
*     EISPACK ROUTINE.
*
*     CONVERGENCE TEST WAS NOT MODIFIED, SINCE IT SHOULD GIVE
*     APPROXIMATELY THE SAME LEVEL OF ACCURACY AS LAPACK ROUTINE,
*     ALTHOUGH THE EIGENVECTORS MAY NOT BE AS CLOSE TO ORTHOGONAL.
*
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
      COMMON /PYTHOP/OPST
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS,OPST
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0E0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0E0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE.
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES.
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT
C
C        ALL INPUT ARRAYS ARE UNALTERED.
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS.
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
*        INITIALIZE ITERATION COUNT.
C     .. Scalar Arguments ..
      INTEGER IERR,M,N,NM
C     ..
C     .. Array Arguments ..
      REAL D(N),E(N),E2(N),RV1(N),RV2(N),RV3(N),RV4(N),RV6(N),W(M),
     +     Z(NM,M)
      INTEGER IND(M)
C     ..
C     .. Local Scalars ..
      REAL EPS2,EPS3,EPS4,NORM,ORDER,U,UK,V,X0,X1,XU
      INTEGER GROUP,I,II,IP,ITS,J,JJ,P,Q,R,S,TAG
C     ..
C     .. External Functions ..
      REAL EPSLON,PYTHAG
      EXTERNAL EPSLON,PYTHAG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AMAX1,SQRT
C     ..
      ITCNT = 0
      IERR = 0
      IF (M.NE.0) THEN
          TAG = 0
          ORDER = 1.0E0 - E2(1)
          Q = 0
   10     CONTINUE
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
          P = Q + 1
C
          DO 20 Q = P,N
              IF (Q.EQ.N) THEN
                  GO TO 30
              ELSE IF (E2(Q+1).EQ.0.0E0) THEN
                  GO TO 30
              END IF
   20     CONTINUE
C     .......... FIND VECTORS BY INVERSE ITERATION ..........
   30     TAG = TAG + 1
          S = 0
C
          DO 210 R = 1,M
              IF (IND(R).EQ.TAG) THEN
                  ITS = 1
                  X1 = W(R)
                  IF (S.EQ.0) THEN
C     .......... CHECK FOR ISOLATED ROOT ..........
                      XU = 1.0E0
                      IF (P.NE.Q) THEN
                          NORM = ABS(D(P))
                          IP = P + 1
C
                          DO 40 I = IP,Q
                              NORM = AMAX1(NORM,ABS(D(I))+ABS(E(I)))
   40                     CONTINUE
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
                          EPS2 = 1.0E-3*NORM
                          EPS3 = EPSLON(NORM)
                          UK = Q - P + 1
                          EPS4 = UK*EPS3
                          UK = EPS4/SQRT(UK)
*           INCREMENT OPCOUNT FOR COMPUTING CRITERIA.
                          OPS = OPS + (Q-IP+4)
                          S = P
                      ELSE
                          RV6(P) = 1.0E0
                          GO TO 180
                      END IF
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
                  ELSE IF (ABS(X1-X0).LT.EPS2) THEN
                      GROUP = GROUP + 1
                      IF (ORDER* (X1-X0).LE.0.0E0) X1 = X0 + ORDER*EPS3
                      GO TO 50
                  END IF
                  GROUP = 0
C     .......... ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ..........
   50             V = 0.0E0
C
                  DO 60 I = P,Q
                      RV6(I) = UK
                      IF (I.NE.P) THEN
                          IF (ABS(E(I)).LT.ABS(U)) THEN
                              XU = E(I)/U
                              RV4(I) = XU
                              RV1(I-1) = U
                              RV2(I-1) = V
                              RV3(I-1) = 0.0E0
                          ELSE
C     .......... WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ..........
                              XU = U/E(I)
                              RV4(I) = XU
                              RV1(I-1) = E(I)
                              RV2(I-1) = D(I) - X1
                              RV3(I-1) = 0.0E0
                              IF (I.NE.Q) RV3(I-1) = E(I+1)
                              U = V - XU*RV2(I-1)
                              V = -XU*RV3(I-1)
                              GO TO 60
                          END IF
                      END IF
                      U = D(I) - X1 - XU*V
                      IF (I.NE.Q) V = E(I+1)
   60             CONTINUE
*           INCREMENT OPCOUNT FOR ELIMINATION.
                  OPS = OPS + (Q-P+1)*5
C
                  IF (U.EQ.0.0E0) U = EPS3
                  RV1(Q) = U
                  RV2(Q) = 0.0E0
                  RV3(Q) = 0.0E0
   70             CONTINUE
C     .......... BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- ..........
                  DO 80 II = P,Q
                      I = P + Q - II
                      RV6(I) = (RV6(I)-U*RV2(I)-V*RV3(I))/RV1(I)
                      V = U
                      U = RV6(I)
   80             CONTINUE
*           INCREMENT OPCOUNT FOR BACK SUBSTITUTION.
                  OPS = OPS + (Q-P+1)*5
C     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ..........
                  IF (GROUP.NE.0) THEN
                      J = R
C
                      DO 120 JJ = 1,GROUP
   90                     CONTINUE
                          J = J - 1
                          IF (IND(J).NE.TAG) GO TO 90
                          XU = 0.0E0
C
                          DO 100 I = P,Q
                              XU = XU + RV6(I)*Z(I,J)
  100                     CONTINUE
C
                          DO 110 I = P,Q
                              RV6(I) = RV6(I) - XU*Z(I,J)
  110                     CONTINUE
C
*              INCREMENT OPCOUNT FOR ORTHOGONALIZING.
                          OPS = OPS + (Q-P+1)*4
  120                 CONTINUE
                  END IF
C
                  NORM = 0.0E0
C
                  DO 130 I = P,Q
                      NORM = NORM + ABS(RV6(I))
  130             CONTINUE
*           INCREMENT OPCOUNT FOR COMPUTING NORM.
                  OPS = OPS + (Q-P+1)
C
                  IF (NORM.GE.1.0E0) THEN
                      GO TO 160
C     .......... FORWARD SUBSTITUTION ..........
                  ELSE IF (ITS.NE.5) THEN
                      IF (NORM.NE.0.0E0) THEN
                          XU = EPS4/NORM
C
                          DO 140 I = P,Q
                              RV6(I) = RV6(I)*XU
  140                     CONTINUE
                      ELSE
                          RV6(S) = EPS4
                          S = S + 1
                          IF (S.GT.Q) S = P
                      END IF
C     .......... ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE ..........
                      DO 150 I = IP,Q
                          U = RV6(I)
C     .......... IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS ..........
                          IF (RV1(I-1).EQ.E(I)) THEN
                              U = RV6(I-1)
                              RV6(I-1) = RV6(I)
                          END IF
                          RV6(I) = U - RV4(I)*RV6(I-1)
  150                 CONTINUE
*           INCREMENT OPCOUNT FOR FORWARD SUBSTITUTION.
                      OPS = OPS + (Q-P+1) + (Q-IP+1)*2
C
                      ITS = ITS + 1
*           INCREMENT ITERATION COUNTER.
                      ITCNT = ITCNT + 1
                      GO TO 70
                  END IF
C     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
                  IERR = -R
                  XU = 0.0E0
                  GO TO 180
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  160             U = 0.0E0
C
                  DO 170 I = P,Q
                      U = PYTHAG(U,RV6(I))
  170             CONTINUE
C
                  XU = 1.0E0/U
C*PL*ERROR* Embedded comment after label moved
C
  180             DO 190 I = 1,N
                      Z(I,R) = 0.0E0
  190             CONTINUE
C
                  DO 200 I = P,Q
                      Z(I,R) = RV6(I)*XU
  200             CONTINUE
*           INCREMENT OPCOUNT FOR NORMALIZING.
                  OPS = OPS + (Q-P+1)
C
                  X0 = X1
              END IF
  210     CONTINUE
C
          IF (Q.LT.N) GO TO 10
*        INCREMENT OPCOUNT FOR USE OF FUNCTION PYTHAG.
          OPS = OPS + OPST
      END IF
      RETURN
      END
      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
*
*     EISPACK ROUTINE.
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN SSTEBZ.
*
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BISECT,
C     NUM. MATH. 9, 386-393(1967) BY BARTH, MARTIN, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX BETWEEN SPECIFIED BOUNDARY INDICES,
C     USING BISECTION.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C        M11 SPECIFIES THE LOWER BOUNDARY INDEX FOR THE DESIRED
C          EIGENVALUES.
C
C        M SPECIFIES THE NUMBER OF EIGENVALUES DESIRED.  THE UPPER
C          BOUNDARY INDEX M22 IS THEN OBTAINED AS M22=M11+M-1.
C
C     ON OUTPUT
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE.
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        LB AND UB DEFINE AN INTERVAL CONTAINING EXACTLY THE DESIRED
C          EIGENVALUES.
C
C        W CONTAINS, IN ITS FIRST M POSITIONS, THE EIGENVALUES
C          BETWEEN INDICES M11 AND M22 IN ASCENDING ORDER.
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF MULTIPLE EIGENVALUES AT INDEX M11 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE,
C          3*N+2      IF MULTIPLE EIGENVALUES AT INDEX M22 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE.
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     NOTE THAT SUBROUTINE TQL1, IMTQL1, OR TQLRAT IS GENERALLY FASTER
C     THAN TRIDIB, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
*        INITIALIZE ITERATION COUNT.
C     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E0)
      REAL RELFAC
      PARAMETER (RELFAC=2.0E0)
C     ..
C     .. Scalar Arguments ..
      REAL EPS1,LB,UB
      INTEGER IERR,M,M11,N
C     ..
C     .. Array Arguments ..
      REAL D(N),E(N),E2(N),RV4(N),RV5(N),W(M)
      INTEGER IND(M)
C     ..
C     .. Local Scalars ..
      REAL ATOLI,PIVMIN,RTOLI,SAFEMN,T1,T2,TMP1,TMP2,TNORM,U,ULP,V,X0,
     +     X1,XU
      INTEGER I,II,ISTURM,J,K,L,M1,M2,M22,P,Q,R,S,TAG
C     ..
C     .. External Functions ..
      REAL EPSLON,SLAMCH
      EXTERNAL EPSLON,SLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AMAX1,AMIN1,MAX
C     ..
      ITCNT = 0
      SAFEMN = SLAMCH('S')
      ULP = SLAMCH('E')*SLAMCH('B')
      RTOLI = ULP*RELFAC
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0E0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES ..........
      PIVMIN = ONE
      DO 30 I = 1,N
          X1 = U
          U = 0.0E0
          IF (I.NE.N) U = ABS(E(I+1))
          XU = AMIN1(D(I)- (X1+U),XU)
          X0 = AMAX1(D(I)+ (X1+U),X0)
          IF (I.EQ.1) GO TO 10
CCC         TST1 = ABS(D(I)) + ABS(D(I-1))
CCC         TST2 = TST1 + ABS(E(I))
CCC         IF (TST2 .GT. TST1) GO TO 40
          TMP1 = E(I)**2
          IF (ABS(D(I)*D(I-1))*ULP**2+SAFEMN.LE.TMP1) THEN
              PIVMIN = MAX(PIVMIN,TMP1)
              GO TO 20
          END IF
   10     E2(I) = 0.0E0
   20     CONTINUE
   30 CONTINUE
      PIVMIN = PIVMIN*SAFEMN
      TNORM = MAX(ABS(XU),ABS(X0))
      ATOLI = ULP*TNORM
*        INCREMENT OPCOUNT FOR DETERMINING IF MATRIX SPLITS.
      OPS = OPS + 9* (N-1)
C
      X1 = N
      X1 = X1*EPSLON(AMAX1(ABS(XU),ABS(X0)))
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
C     .......... DETERMINE AN INTERVAL CONTAINING EXACTLY
C                THE DESIRED EIGENVALUES ..........
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1.EQ.0) GO TO 90
      ISTURM = 1
   40 V = X1
      X1 = XU + (X0-XU)*0.5E0
      IF (X1.EQ.V) GO TO 410
      GO TO 260
   50 IF (S-M1) 60,80,70
   60 XU = X1
      GO TO 40
   70 X0 = X1
      GO TO 40
   80 XU = X1
      T1 = X1
   90 M22 = M1 + M
      IF (M22.EQ.N) GO TO 120
      X0 = T2
      ISTURM = 2
      GO TO 40
  100 IF (S-M22) 60,110,70
  110 T2 = X1
  120 Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  130 IF (R.EQ.M) GO TO 420
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0E0
C
      DO 150 Q = P,N
          X1 = U
          U = 0.0E0
          V = 0.0E0
          IF (Q.EQ.N) GO TO 140
          U = ABS(E(Q+1))
          V = E2(Q+1)
  140     XU = AMIN1(D(Q)- (X1+U),XU)
          X0 = AMAX1(D(Q)+ (X1+U),X0)
          IF (V.EQ.0.0E0) GO TO 160
  150 CONTINUE
*        INCREMENT OPCOUNT FOR REFINING INTERVAL.
      OPS = OPS + (N-P+1)*2
C
  160 X1 = EPSLON(AMAX1(ABS(XU),ABS(X0)))
      IF (EPS1.LE.0.0E0) EPS1 = -X1
      IF (P.NE.Q) GO TO 170
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1.GT.D(P) .OR. D(P).GE.T2) GO TO 400
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 340
  170 X1 = X1* (Q-P+1)
      LB = AMAX1(T1,XU-X1)
      UB = AMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 260
  180 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 260
  190 M2 = S
      IF (M1.GT.M2) GO TO 400
C     .......... FIND ROOTS BY BISECTION ..........
      X0 = UB
      ISTURM = 5
C
      DO 200 I = M1,M2
          RV5(I) = UB
          RV4(I) = LB
  200 CONTINUE
C     .......... LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
      K = M2
  210 XU = LB
C     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
      DO 230 II = M1,K
          I = M1 + K - II
          IF (XU.GE.RV4(I)) GO TO 220
          XU = RV4(I)
          GO TO 240
  220     CONTINUE
  230 CONTINUE
C
  240 IF (X0.GT.RV5(K)) X0 = RV5(K)
C     .......... NEXT BISECTION STEP ..........
  250 X1 = (XU+X0)*0.5E0
CCC         IF ((X0 - XU) .LE. ABS(EPS1)) GO TO 420
CCC         TST1 = 2.0E0 * (ABS(XU) + ABS(X0))
CCC         TST2 = TST1 + (X0 - XU)
CCC         IF (TST2 .EQ. TST1) GO TO 420
      TMP1 = ABS(X0-XU)
      TMP2 = MAX(ABS(X0),ABS(XU))
      IF (TMP1.LT.MAX(ATOLI,PIVMIN,RTOLI*TMP2)) GO TO 330
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  260 S = P - 1
      U = 1.0E0
C
      DO 290 I = P,Q
          IF (U.NE.0.0E0) GO TO 270
          V = ABS(E(I))/EPSLON(1.0E0)
          IF (E2(I).EQ.0.0E0) V = 0.0E0
          GO TO 280
  270     V = E2(I)/U
  280     U = D(I) - X1 - V
          IF (U.LT.0.0E0) S = S + 1
  290 CONTINUE
*           INCREMENT OPCOUNT FOR STURM SEQUENCE.
      OPS = OPS + (Q-P+1)*3
*           INCREMENT ITERATION COUNTER.
      ITCNT = ITCNT + 1
C
      GO TO (50,100,180,190,300) ISTURM
C     .......... REFINE INTERVALS ..........
  300 IF (S.GE.K) GO TO 320
      XU = X1
      IF (S.GE.M1) GO TO 310
      RV4(M1) = X1
      GO TO 250
  310 RV4(S+1) = X1
      IF (RV5(S).GT.X1) RV5(S) = X1
      GO TO 250
  320 X0 = X1
      GO TO 250
C     .......... K-TH EIGENVALUE FOUND ..........
  330 RV5(K) = X1
      K = K - 1
      IF (K.GE.M1) GO TO 210
C     .......... ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ..........
  340 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 390 L = 1,R
          IF (J.GT.S) GO TO 360
          IF (K.GT.M2) GO TO 400
          IF (RV5(K).GE.W(L)) GO TO 370
C
          DO 350 II = J,S
              I = L + S - II
              W(I+1) = W(I)
              IND(I+1) = IND(I)
  350     CONTINUE
C
  360     W(L) = RV5(K)
          IND(L) = TAG
          K = K + 1
          GO TO 380
  370     J = J + 1
  380     CONTINUE
  390 CONTINUE
C
  400 IF (Q.LT.N) GO TO 130
      GO TO 420
C     .......... SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
C                EXACTLY THE DESIRED EIGENVALUES ..........
  410 IERR = 3*N + ISTURM
  420 LB = T1
      UB = T2
      RETURN
      END
      SUBROUTINE CSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, IOPS IS ONLY INCREMENTED
*     IOPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO IOPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/IOPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL IOPS,ITCNT
C     ..
C
C
C     CSVDC IS A SUBROUTINE TO REDUCE A COMPLEX NXP MATRIX X BY
C     UNITARY TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
C     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
C     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
C     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
C
C     ON ENTRY
C
C         X         COMPLEX(LDX,P), WHERE LDX.GE.N.
C                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
C                   DECOMPOSITION IS TO BE COMPUTED.  X IS
C                   DESTROYED BY CSVDC.
C
C         LDX       INTEGER.
C                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C         N         INTEGER.
C                   N IS THE NUMBER OF ROWS OF THE MATRIX X.
C
C         P         INTEGER.
C                   P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
C
C         LDU       INTEGER.
C                   LDU IS THE LEADING DIMENSION OF THE ARRAY U
C                   (SEE BELOW).
C
C         LDV       INTEGER.
C                   LDV IS THE LEADING DIMENSION OF THE ARRAY V
C                   (SEE BELOW).
C
C         WORK      COMPLEX(N).
C                   WORK IS A SCRATCH ARRAY.
C
C         JOB       INTEGER.
C                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR
C                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB
C                   WITH THE FOLLOWING MEANING
C
C                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR
C                                  VECTORS.
C                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS
C                                  IN U.
C                        A.GE.2    RETURNS THE FIRST MIN(N,P)
C                                  LEFT SINGULAR VECTORS IN U.
C                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
C                                  VECTORS.
C                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
C                                  IN V.
C
C     ON RETURN
C
C         S         COMPLEX(MM), WHERE MM=MIN(N+1,P).
C                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
C                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
C                   ORDER OF MAGNITUDE.
C
C         E         COMPLEX(P).
C                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
C                   DISCUSSION OF INFO FOR EXCEPTIONS.
C
C         U         COMPLEX(LDU,K), WHERE LDU.GE.N.  IF JOBA.EQ.1 THEN
C                                   K.EQ.N, IF JOBA.GE.2 THEN
C                                   K.EQ.MIN(N,P).
C                   U CONTAINS THE MATRIX OF LEFT SINGULAR VECTORS.
C                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
C                   OR IF JOBA.GT.2, THEN U MAY BE IDENTIFIED WITH X
C                   IN THE SUBROUTINE CALL.
C
C         V         COMPLEX(LDV,P), WHERE LDV.GE.P.
C                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   V IS NOT REFERENCED IF JOBB.EQ.0.  IF P.LE.N,
C                   THEN V MAY BE IDENTIFIED WHTH X IN THE
C                   SUBROUTINE CALL.
C
C         INFO      INTEGER.
C                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
C                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
C                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
C                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
C                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
C                   B = CTRANS(U)*X*V IS THE BIDIAGONAL MATRIX
C                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
C                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (CTRANS(U)
C                   IS THE CONJUGATE-TRANSPOSE OF U).  THUS THE
C                   SINGULAR VALUES OF X AND B ARE THE SAME.
C
C     LINPACK. THIS VERSION DATED 03/19/79 .
C              CORRECTION TO SHIFT CALCULATION MADE 2/85.
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     CSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     EXTERNAL CSROT
C     BLAS CAXPY,CDOTC,CSCAL,CSWAP,SCNRM2,SROTG
C     FORTRAN ABS,AIMAG,AMAX1,CABS,CMPLX
C     FORTRAN CONJG,MAX0,MIN0,MOD,REAL,SQRT
C
C     INTERNAL VARIABLES
C
*     REAL ZTEST
C
*
*     DECLARE EPS AND SLAMCH FOR NEW STOPPING CRITERION
*
*
*     GET EPS FROM SLAMCH FOR NEW STOPPING CRITERION
C     .. Scalar Arguments ..
      INTEGER INFO,JOB,LDU,LDV,LDX,N,P
C     ..
C     .. Array Arguments ..
      COMPLEX E(*),S(*),U(LDU,*),V(LDV,*),WORK(*),X(LDX,*)
C     ..
C     .. Local Scalars ..
      COMPLEX R,T,ZDUM,ZDUM1,ZDUM2
      REAL B,C,CS,EL,EMM1,EPS,F,G,IOPST,SCALE,SHIFT,SL,SM,SMM1,SN,T1,
     +     TEST
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,MM,
     +        MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      LOGICAL WANTU,WANTV
C     ..
C     .. External Functions ..
      COMPLEX CDOTC
      REAL SCNRM2,SLAMCH
      EXTERNAL CDOTC,SCNRM2,SLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL CAXPY,CSCAL,CSROT,CSWAP,SROTG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,AMAX1,CABS,CMPLX,CONJG,FLOAT,MAX,MAX0,MIN,
     +          MIN0,MOD,REAL,SQRT
C     ..
C     .. Statement Functions ..
      COMPLEX CSIGN
      REAL CABS1
C     ..
C     .. Statement Function definitions ..
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN(ZDUM1,ZDUM2) = CABS(ZDUM1)* (ZDUM2/CABS(ZDUM2))
C     ..
      IF (N.GT.0 .AND. P.GT.0) THEN
          EPS = SLAMCH('EPSILON')
*
C
C     SET THE MAXIMUM NUMBER OF ITERATIONS.
C
          MAXIT = 50
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
          WANTU = .FALSE.
          WANTV = .FALSE.
          JOBU = MOD(JOB,100)/10
          NCU = N
          IF (JOBU.GT.1) NCU = MIN0(N,P)
          IF (JOBU.NE.0) WANTU = .TRUE.
          IF (MOD(JOB,10).NE.0) WANTV = .TRUE.
C
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
C
*
*     INITIALIZE OP COUNT
          IOPST = 0
          INFO = 0
          NCT = MIN0(N-1,P)
          NRT = MAX0(0,MIN0(P-2,N))
          LU = MAX0(NCT,NRT)
          IF (LU.GE.1) THEN
              DO 70 L = 1,LU
                  LP1 = L + 1
                  IF (L.LE.NCT) THEN
C
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
C           PLACE THE L-TH DIAGONAL IN S(L).
C
*
*           INCREMENT OP COUNT
                      IOPS = IOPS + (4* (N-L+1)+2)
                      S(L) = CMPLX(SCNRM2(N-L+1,X(L,L),1),0.0E0)
                      IF (CABS1(S(L)).NE.0.0E0) THEN
                          IF (CABS1(X(L,L)).NE.0.0E0) S(L) = CSIGN(S(L),
     +                        X(L,L))
*
*              INCREMENT OP COUNT
                          IOPS = IOPS + (6* (N-L+1)+23)
                          CALL CSCAL(N-L+1,1.0E0/S(L),X(L,L),1)
                          X(L,L) = (1.0E0,0.0E0) + X(L,L)
                      END IF
                      S(L) = -S(L)
                  END IF
                  IF (P.GE.LP1) THEN
                      DO 10 J = LP1,P
                          IF (L.LE.NCT) THEN
                              IF (CABS1(S(L)).NE.0.0E0) THEN
C
C              APPLY THE TRANSFORMATION.
C
*
*              INCREMENT OP COUNT
                                  IOPS = IOPS + (16* (N-L)+26)
                                  T = -CDOTC(N-L+1,X(L,L),1,X(L,J),1)/
     +                                X(L,L)
                                  CALL CAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                              END IF
                          END IF
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
                          E(J) = CONJG(X(L,J))
   10                 CONTINUE
                  END IF
                  IF (WANTU .AND. L.LE.NCT) THEN
C
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
C           MULTIPLICATION.
C
                      DO 20 I = L,N
                          U(I,L) = X(I,L)
   20                 CONTINUE
                  END IF
                  IF (L.LE.NRT) THEN
C
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
C           L-TH SUPER-DIAGONAL IN E(L).
C
*
*           INCREMENT OP COUNT
                      IOPS = IOPS + (4* (P-L)+3)
                      E(L) = CMPLX(SCNRM2(P-L,E(LP1),1),0.0E0)
                      IF (CABS1(E(L)).NE.0.0E0) THEN
                          IF (CABS1(E(LP1)).NE.0.0E0) E(L) = CSIGN(E(L),
     +                        E(LP1))
*
*              INCREMENT OP COUNT
                          IOPS = IOPS + (6* (P-L)+23)
                          CALL CSCAL(P-L,1.0E0/E(L),E(LP1),1)
                          E(LP1) = (1.0E0,0.0E0) + E(LP1)
                      END IF
                      E(L) = -CONJG(E(L))
                      IF (LP1.LE.N .AND. CABS1(E(L)).NE.0.0E0) THEN
C
C              APPLY THE TRANSFORMATION.
C
                          DO 30 I = LP1,N
                              WORK(I) = (0.0E0,0.0E0)
   30                     CONTINUE
*
*              INCREMENT OP COUNT
                          IOPS = IOPS + FLOAT(16* (N-L)+9)* (P-L)
                          DO 40 J = LP1,P
                              CALL CAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),
     +                                   1)
   40                     CONTINUE
                          DO 50 J = LP1,P
                              CALL CAXPY(N-L,CONJG(-E(J)/E(LP1)),
     +                                   WORK(LP1),1,X(LP1,J),1)
   50                     CONTINUE
                      END IF
                      IF (WANTV) THEN
C
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
C              BACK MULTIPLICATION.
C
                          DO 60 I = LP1,P
                              V(I,L) = E(I)
   60                     CONTINUE
                      END IF
                  END IF
   70         CONTINUE
          END IF
C
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
C
          M = MIN0(P,N+1)
          NCTP1 = NCT + 1
          NRTP1 = NRT + 1
          IF (NCT.LT.P) S(NCTP1) = X(NCTP1,NCTP1)
          IF (N.LT.M) S(M) = (0.0E0,0.0E0)
          IF (NRTP1.LT.M) E(NRTP1) = X(NRTP1,M)
          E(M) = (0.0E0,0.0E0)
C
C     IF REQUIRED, GENERATE U.
C
          IF (WANTU) THEN
              IF (NCU.GE.NCTP1) THEN
                  DO 90 J = NCTP1,NCU
                      DO 80 I = 1,N
                          U(I,J) = (0.0E0,0.0E0)
   80                 CONTINUE
                      U(J,J) = (1.0E0,0.0E0)
   90             CONTINUE
              END IF
              IF (NCT.GE.1) THEN
                  DO 130 LL = 1,NCT
                      L = NCT - LL + 1
                      IF (CABS1(S(L)).EQ.0.0E0) THEN
                          DO 100 I = 1,N
                              U(I,L) = (0.0E0,0.0E0)
  100                     CONTINUE
                          U(L,L) = (1.0E0,0.0E0)
                      ELSE
                          LP1 = L + 1
                          IF (NCU.GE.LP1) THEN
*
*              INCREMENT OP COUNT
                              IOPS = IOPS + (FLOAT(16* (N-L)+25)*
     +                               (NCU-L)+6* (N-L)+9)
                              DO 110 J = LP1,NCU
                                  T = -CDOTC(N-L+1,U(L,L),1,U(L,J),1)/
     +                                U(L,L)
                                  CALL CAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  110                         CONTINUE
                          END IF
                          CALL CSCAL(N-L+1, (-1.0E0,0.0E0),U(L,L),1)
                          U(L,L) = (1.0E0,0.0E0) + U(L,L)
                          LM1 = L - 1
                          IF (LM1.GE.1) THEN
                              DO 120 I = 1,LM1
                                  U(I,L) = (0.0E0,0.0E0)
  120                         CONTINUE
                          END IF
                      END IF
  130             CONTINUE
              END IF
          END IF
C
C     IF IT IS REQUIRED, GENERATE V.
C
          IF (WANTV) THEN
              DO 160 LL = 1,P
                  L = P - LL + 1
                  LP1 = L + 1
                  IF (L.LE.NRT) THEN
                      IF (CABS1(E(L)).NE.0.0E0) THEN
*
*              INCREMENT OP COUNT
                          IOPS = IOPS + (FLOAT(16* (P-L)+9)* (P-L)+1)
                          DO 140 J = LP1,P
                              T = -CDOTC(P-L,V(LP1,L),1,V(LP1,J),1)/
     +                            V(LP1,L)
                              CALL CAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  140                     CONTINUE
                      END IF
                  END IF
                  DO 150 I = 1,P
                      V(I,L) = (0.0E0,0.0E0)
  150             CONTINUE
                  V(L,L) = (1.0E0,0.0E0)
  160         CONTINUE
          END IF
C
C     TRANSFORM S AND E SO THAT THEY ARE REAL.
C
*
*     INCREMENT OP COUNT
          IOPS = IOPS + (2*M-1)
          DO 170 I = 1,M
              IF (CABS1(S(I)).NE.0.0E0) THEN
*
*           INCREMENT OP COUNT
                  IOPS = IOPS + 23
                  IF (WANTU) IOPS = IOPS + 6*N
                  T = CMPLX(CABS(S(I)),0.0E0)
                  R = S(I)/T
                  S(I) = T
                  IF (I.LT.M) E(I) = E(I)/R
                  IF (WANTU) CALL CSCAL(N,R,U(1,I),1)
              END IF
C     ...EXIT
              IF (I.EQ.M) THEN
                  GO TO 180
              ELSE IF (CABS1(E(I)).NE.0.0E0) THEN
*
*           INCREMENT OP COUNT
                  IOPS = IOPS + 20
                  IF (WANTV) IOPS = IOPS + 6*P
                  T = CMPLX(CABS(E(I)),0.0E0)
                  R = T/E(I)
                  E(I) = T
                  S(I+1) = S(I+1)*R
                  IF (WANTV) CALL CSCAL(P,R,V(1,I+1),1)
              END IF
  170     CONTINUE
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
  180     MM = M
*
*     INITIALIZE ITERATION COUNTER
          ITCNT = 0
          ITER = 0
  190     CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
C     ...EXIT
          IF (M.EQ.0) THEN
              GO TO 340
          ELSE
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
*
*        UPDATE ITERATION COUNTER
              ITCNT = ITER
              IF (ITER.LT.MAXIT) THEN
C
C        THIS SECTION OF THE PROGRAM INSPECTS FOR
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
C
C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
C
                  DO 200 LL = 1,M
                      L = M - LL
C        ...EXIT
                      IF (L.EQ.0) THEN
                          GO TO 220
                      ELSE
*
*           INCREMENT OP COUNT
                          IOPST = IOPST + 17
                          TEST = CABS(S(L)) + CABS(S(L+1))
*
*           REPLACE STOPPING CRITERION WITH NEW ONE
*
*           ZTEST = TEST + CABS(E(L))
*           IF (ZTEST .NE. TEST) GO TO 420
                          IF (CABS(E(L)).LE.EPS*TEST) GO TO 210
                      END IF
  200             CONTINUE
                  GO TO 220
*
  210             E(L) = (0.0E0,0.0E0)
C        ......EXIT
  220             IF (L.NE.M-1) THEN
                      LP1 = L + 1
                      MP1 = M + 1
                      DO 230 LLS = LP1,MP1
                          LS = M - LLS + LP1
C           ...EXIT
                          IF (LS.EQ.L) THEN
                              GO TO 250
                          ELSE
                              TEST = 0.0E0
*
*              INCREMENT OP COUNT
                              IOPST = IOPST + 18
                              IF (LS.NE.M) TEST = TEST + CABS(E(LS))
                              IF (LS.NE.L+1) TEST = TEST + CABS(E(LS-1))
*
*              REPLACE STOPPING CRITERION WITH NEW ONE AS IN LAPACK
*
*              ZTEST = TEST + CABS(S(LS))
*              IF (ZTEST .NE. TEST) GO TO 460
                              IF (CABS(S(LS)).LE.EPS*TEST) GO TO 240
                          END IF
  230                 CONTINUE
                      GO TO 250
*
  240                 S(LS) = (0.0E0,0.0E0)
C           ......EXIT
  250                 IF (LS.EQ.L) THEN
                          KASE = 3
                      ELSE IF (LS.NE.M) THEN
                          KASE = 2
                          L = LS
                      ELSE
                          KASE = 1
                      END IF
                  ELSE
                      KASE = 4
                  END IF
                  L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C
                  GO TO (260,280,300,320) KASE
C
C        DEFLATE NEGLIGIBLE S(M).
C
  260             CONTINUE
                  MM1 = M - 1
                  F = REAL(E(M-1))
                  E(M-1) = (0.0E0,0.0E0)
*
*           INCREMENT OP COUNT
                  IOPS = IOPS + ((MM1-L+1)*14-3)
                  IF (WANTV) IOPS = IOPS + FLOAT(MM1-L+1)*12*P
                  DO 270 KK = L,MM1
                      K = MM1 - KK + L
                      T1 = REAL(S(K))
                      CALL SROTG(T1,F,CS,SN)
                      S(K) = CMPLX(T1,0.0E0)
                      IF (K.NE.L) THEN
                          F = -SN*REAL(E(K-1))
                          E(K-1) = CS*E(K-1)
                      END IF
                      IF (WANTV) CALL CSROT(P,V(1,K),1,V(1,M),1,CS,SN)
  270             CONTINUE
                  GO TO 190
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  280             CONTINUE
                  F = REAL(E(L-1))
                  E(L-1) = (0.0E0,0.0E0)
*
*           INCREMENT OP COUNT
                  IOPS = IOPS + (M-L+1)*14
                  IF (WANTU) IOPS = IOPS + FLOAT(M-L+1)*12*N
                  DO 290 K = L,M
                      T1 = REAL(S(K))
                      CALL SROTG(T1,F,CS,SN)
                      S(K) = CMPLX(T1,0.0E0)
                      F = -SN*REAL(E(K))
                      E(K) = CS*E(K)
                      IF (WANTU) CALL CSROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  290             CONTINUE
                  GO TO 190
C
C        PERFORM ONE QR STEP.
C
  300             CONTINUE
C
C           CALCULATE THE SHIFT.
C
*
*           INCREMENT OP COUNT
                  IOPST = IOPST + 48
                  SCALE = AMAX1(CABS(S(M)),CABS(S(M-1)),CABS(E(M-1)),
     +                    CABS(S(L)),CABS(E(L)))
                  SM = REAL(S(M))/SCALE
                  SMM1 = REAL(S(M-1))/SCALE
                  EMM1 = REAL(E(M-1))/SCALE
                  SL = REAL(S(L))/SCALE
                  EL = REAL(E(L))/SCALE
                  B = ((SMM1+SM)* (SMM1-SM)+EMM1**2)/2.0E0
                  C = (SM*EMM1)**2
                  SHIFT = 0.0E0
                  IF (B.NE.0.0E0 .OR. C.NE.0.0E0) THEN
                      SHIFT = SQRT(B**2+C)
                      IF (B.LT.0.0E0) SHIFT = -SHIFT
                      SHIFT = C/ (B+SHIFT)
                  END IF
                  F = (SL+SM)* (SL-SM) + SHIFT
                  G = SL*EL
C
C           CHASE ZEROS.
C
                  MM1 = M - 1
*
*           INCREMENT OP COUNT
                  IOPS = IOPS + (MM1-L+1)*46
                  IF (WANTV) IOPS = IOPS + FLOAT(MM1-L+1)*12*P
                  IF (WANTU) IOPS = IOPS + FLOAT(MAX((MIN(MM1,N-1)-L+1),
     +                              0))*12*N
                  DO 310 K = L,MM1
                      CALL SROTG(F,G,CS,SN)
                      IF (K.NE.L) E(K-1) = CMPLX(F,0.0E0)
                      F = CS*REAL(S(K)) + SN*REAL(E(K))
                      E(K) = CS*E(K) - SN*S(K)
                      G = SN*REAL(S(K+1))
                      S(K+1) = CS*S(K+1)
                      IF (WANTV) CALL CSROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
                      CALL SROTG(F,G,CS,SN)
                      S(K) = CMPLX(F,0.0E0)
                      F = CS*REAL(E(K)) + SN*REAL(S(K+1))
                      S(K+1) = -SN*E(K) + CS*S(K+1)
                      G = SN*REAL(E(K+1))
                      E(K+1) = CS*E(K+1)
                      IF (WANTU .AND. K.LT.N) CALL CSROT(N,U(1,K),1,
     +                    U(1,K+1),1,CS,SN)
  310             CONTINUE
                  E(M-1) = CMPLX(F,0.0E0)
                  ITER = ITER + 1
                  GO TO 190
C
C        CONVERGENCE.
C
  320             CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE
C
                  IF (REAL(S(L)).LT.0.0E0) THEN
                      S(L) = -S(L)
*
*              INCREMENT OP COUNT
                      IF (WANTV) IOPS = IOPS + 6*P
                      IF (WANTV) CALL CSCAL(P, (-1.0E0,0.0E0),V(1,L),1)
                  END IF
  330             CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
                  IF (L.NE.MM) THEN
C           ...EXIT
                      IF (REAL(S(L)).LT.REAL(S(L+1))) THEN
                          T = S(L)
                          S(L) = S(L+1)
                          S(L+1) = T
                          IF (WANTV .AND. L.LT.P) CALL CSWAP(P,V(1,L),1,
     +                        V(1,L+1),1)
                          IF (WANTU .AND. L.LT.N) CALL CSWAP(N,U(1,L),1,
     +                        U(1,L+1),1)
                          L = L + 1
                          GO TO 330
                      END IF
                  END IF
                  ITER = 0
                  M = M - 1
                  GO TO 190
              END IF
          END IF
          INFO = M
C     ......EXIT
*
*     COMPUTE FINAL OPCOUNT
  340     IOPS = IOPS + IOPST
      END IF
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE CQZHES(NM,N,AR,AI,BR,BI,MATZ,ZR,ZI)
C
*
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
*     ------------------------ END TIMING CODE -------------------------
*
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-
C     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR
C     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY
C     CQZVAL  AND POSSIBLY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE
C          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE),
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
C          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
C          OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** INITIALIZE Z **********
C     .. Scalar Arguments ..
      INTEGER N,NM
      LOGICAL MATZ
C     ..
C     .. Array Arguments ..
      REAL AI(NM,N),AR(NM,N),BI(NM,N),BR(NM,N),ZI(NM,N),ZR(NM,N)
C     ..
C     .. Local Scalars ..
      REAL OPST,R,RHO,S,T,TI,U1,U1I,U2,XI,XR,YI,YR
      INTEGER I,IOPST,J,K,K1,L,L1,LB,NK1,NM1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,CABS,CMPLX,REAL,SQRT
C     ..
      IF (MATZ) THEN
C
          DO 20 I = 1,N
C
              DO 10 J = 1,N
                  ZR(I,J) = 0.0
                  ZI(I,J) = 0.0
   10         CONTINUE
C
              ZR(I,I) = 1.0
   20     CONTINUE
      END IF
C     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH
C                TEMPORARILY REAL DIAGONAL ELEMENTS **********
      IF (N.GT.1) THEN
          NM1 = N - 1
C
          DO 120 L = 1,NM1
*        ---------------------- BEGIN TIMING CODE ----------------------
              IOPST = 0
*        ----------------------- END TIMING CODE -----------------------
              L1 = L + 1
              S = 0.0
C
              DO 30 I = L,N
                  S = S + ABS(BR(I,L)) + ABS(BI(I,L))
   30         CONTINUE
*        ---------------------- BEGIN TIMING CODE ----------------------
              IOPST = IOPST + 2* (N+1-L)
*        ----------------------- END TIMING CODE -----------------------
C
              IF (S.NE.0.0) THEN
                  RHO = 0.0
C
                  DO 40 I = L,N
                      BR(I,L) = BR(I,L)/S
                      BI(I,L) = BI(I,L)/S
                      RHO = RHO + BR(I,L)**2 + BI(I,L)**2
   40             CONTINUE
C
                  R = SQRT(RHO)
                  XR = CABS(CMPLX(BR(L,L),BI(L,L)))
                  IF (XR.EQ.0.0) THEN
C
                      BR(L,L) = R
                      U1 = -1.0
                      U1I = 0.0
                  ELSE
*        ---------------------- BEGIN TIMING CODE ----------------------
                      IOPST = IOPST + 8
*        ----------------------- END TIMING CODE -----------------------
                      RHO = RHO + XR*R
                      U1 = -BR(L,L)/XR
                      U1I = -BI(L,L)/XR
                      YR = R/XR + 1.0
                      BR(L,L) = YR*BR(L,L)
                      BI(L,L) = YR*BI(L,L)
                  END IF
C
                  DO 70 J = L1,N
                      T = 0.0
                      TI = 0.0
C
                      DO 50 I = L,N
                          T = T + BR(I,L)*BR(I,J) + BI(I,L)*BI(I,J)
                          TI = TI + BR(I,L)*BI(I,J) - BI(I,L)*BR(I,J)
   50                 CONTINUE
C
                      T = T/RHO
                      TI = TI/RHO
C
                      DO 60 I = L,N
                          BR(I,J) = BR(I,J) - T*BR(I,L) + TI*BI(I,L)
                          BI(I,J) = BI(I,J) - T*BI(I,L) - TI*BR(I,L)
   60                 CONTINUE
C
                      XI = U1*BI(L,J) - U1I*BR(L,J)
                      BR(L,J) = U1*BR(L,J) + U1I*BI(L,J)
                      BI(L,J) = XI
   70             CONTINUE
C
                  DO 100 J = 1,N
                      T = 0.0
                      TI = 0.0
C
                      DO 80 I = L,N
                          T = T + BR(I,L)*AR(I,J) + BI(I,L)*AI(I,J)
                          TI = TI + BR(I,L)*AI(I,J) - BI(I,L)*AR(I,J)
   80                 CONTINUE
C
                      T = T/RHO
                      TI = TI/RHO
C
                      DO 90 I = L,N
                          AR(I,J) = AR(I,J) - T*BR(I,L) + TI*BI(I,L)
                          AI(I,J) = AI(I,J) - T*BI(I,L) - TI*BR(I,L)
   90                 CONTINUE
C
                      XI = U1*AI(L,J) - U1I*AR(L,J)
                      AR(L,J) = U1*AR(L,J) + U1I*AI(L,J)
                      AI(L,J) = XI
  100             CONTINUE
C
                  BR(L,L) = R*S
                  BI(L,L) = 0.0
C
                  DO 110 I = L1,N
                      BR(I,L) = 0.0
                      BI(I,L) = 0.0
  110             CONTINUE
*        ---------------------- BEGIN TIMING CODE ----------------------
                  OPS = OPS + (REAL(16* (N-L)+16*N+30)*REAL(N-L)+
     +                  REAL(24*N+13+IOPST))
              END IF
  120     CONTINUE
*        ----------------------- END TIMING CODE -----------------------
C
C     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
C                ELEMENTS, WHILE KEEPING B TRIANGULAR **********
          DO 200 K = 1,NM1
*        ---------------------- BEGIN TIMING CODE ----------------------
              OPST = 0.0
*        ----------------------- END TIMING CODE -----------------------
              K1 = K + 1
C     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL **********
              IF (AI(N,K).NE.0.0) THEN
                  R = CABS(CMPLX(AR(N,K),AI(N,K)))
                  U1 = AR(N,K)/R
                  U1I = AI(N,K)/R
                  AR(N,K) = R
                  AI(N,K) = 0.0
C
                  DO 130 J = K1,N
                      XI = U1*AI(N,J) - U1I*AR(N,J)
                      AR(N,J) = U1*AR(N,J) + U1I*AI(N,J)
                      AI(N,J) = XI
  130             CONTINUE
C
                  XI = U1*BI(N,N) - U1I*BR(N,N)
                  BR(N,N) = U1*BR(N,N) + U1I*BI(N,N)
                  BI(N,N) = XI
*        ---------------------- BEGIN TIMING CODE ----------------------
                  OPST = OPST + REAL(18+6* (N-K))
              END IF
*        ----------------------- END TIMING CODE -----------------------
              IF (K.EQ.NM1) THEN
                  GO TO 210
              ELSE
                  NK1 = NM1 - K
C     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********
                  DO 190 LB = 1,NK1
                      L = N - LB
                      L1 = L + 1
C     ********** ZERO A(L+1,K) **********
                      S = ABS(AR(L,K)) + ABS(AI(L,K)) + AR(L1,K)
                      IF (S.NE.0.0) THEN
*           -------------------- BEGIN TIMING CODE ---------------------
                          OPST = OPST + REAL(18+20* (2*N-K-L))
*           --------------------- END TIMING CODE ----------------------
                          U1 = AR(L,K)/S
                          U1I = AI(L,K)/S
                          U2 = AR(L1,K)/S
                          R = SQRT(U1*U1+U1I*U1I+U2*U2)
                          U1 = U1/R
                          U1I = U1I/R
                          U2 = U2/R
                          AR(L,K) = R*S
                          AI(L,K) = 0.0
                          AR(L1,K) = 0.0
C
                          DO 140 J = K1,N
                              XR = AR(L,J)
                              XI = AI(L,J)
                              YR = AR(L1,J)
                              YI = AI(L1,J)
                              AR(L,J) = U1*XR + U1I*XI + U2*YR
                              AI(L,J) = U1*XI - U1I*XR + U2*YI
                              AR(L1,J) = U1*YR - U1I*YI - U2*XR
                              AI(L1,J) = U1*YI + U1I*YR - U2*XI
  140                     CONTINUE
C
                          XR = BR(L,L)
                          BR(L,L) = U1*XR
                          BI(L,L) = -U1I*XR
                          BR(L1,L) = -U2*XR
C
                          DO 150 J = L1,N
                              XR = BR(L,J)
                              XI = BI(L,J)
                              YR = BR(L1,J)
                              YI = BI(L1,J)
                              BR(L,J) = U1*XR + U1I*XI + U2*YR
                              BI(L,J) = U1*XI - U1I*XR + U2*YI
                              BR(L1,J) = U1*YR - U1I*YI - U2*XR
                              BI(L1,J) = U1*YI + U1I*YR - U2*XI
  150                     CONTINUE
C     ********** ZERO B(L+1,L) **********
                          S = ABS(BR(L1,L1)) + ABS(BI(L1,L1)) +
     +                        ABS(BR(L1,L))
                          IF (S.NE.0.0) THEN
*           -------------------- BEGIN TIMING CODE ---------------------
                              OPST = OPST + REAL(13+20* (N+L))
*           --------------------- END TIMING CODE ----------------------
                              U1 = BR(L1,L1)/S
                              U1I = BI(L1,L1)/S
                              U2 = BR(L1,L)/S
                              R = SQRT(U1*U1+U1I*U1I+U2*U2)
                              U1 = U1/R
                              U1I = U1I/R
                              U2 = U2/R
                              BR(L1,L1) = R*S
                              BI(L1,L1) = 0.0
                              BR(L1,L) = 0.0
C
                              DO 160 I = 1,L
                                  XR = BR(I,L1)
                                  XI = BI(I,L1)
                                  YR = BR(I,L)
                                  YI = BI(I,L)
                                  BR(I,L1) = U1*XR + U1I*XI + U2*YR
                                  BI(I,L1) = U1*XI - U1I*XR + U2*YI
                                  BR(I,L) = U1*YR - U1I*YI - U2*XR
                                  BI(I,L) = U1*YI + U1I*YR - U2*XI
  160                         CONTINUE
C
                              DO 170 I = 1,N
                                  XR = AR(I,L1)
                                  XI = AI(I,L1)
                                  YR = AR(I,L)
                                  YI = AI(I,L)
                                  AR(I,L1) = U1*XR + U1I*XI + U2*YR
                                  AI(I,L1) = U1*XI - U1I*XR + U2*YI
                                  AR(I,L) = U1*YR - U1I*YI - U2*XR
                                  AI(I,L) = U1*YI + U1I*YR - U2*XI
  170                         CONTINUE
C
                              IF (MATZ) THEN
*           -------------------- BEGIN TIMING CODE ---------------------
                                  OPST = OPST + 20*N
*           --------------------- END TIMING CODE ----------------------
C
                                  DO 180 I = 1,N
                                      XR = ZR(I,L1)
                                      XI = ZI(I,L1)
                                      YR = ZR(I,L)
                                      YI = ZI(I,L)
                                      ZR(I,L1) = U1*XR + U1I*XI + U2*YR
                                      ZI(I,L1) = U1*XI - U1I*XR + U2*YI
                                      ZR(I,L) = U1*YR - U1I*YI - U2*XR
                                      ZI(I,L) = U1*YI + U1I*YR - U2*XI
  180                             CONTINUE
                              END IF
                          END IF
                      END IF
  190             CONTINUE
C
*        ---------------------- BEGIN TIMING CODE ----------------------
                  OPS = OPS + (OPST+REAL(2* (N-1-K)))
              END IF
  200     CONTINUE
*        ----------------------- END TIMING CODE -----------------------
C
      END IF
  210 RETURN
C     ********** LAST CARD OF CQZHES **********
*$st$ Unreachable comments ...
C
      END
      SUBROUTINE CQZVAL(NM,N,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,MATZ,ZR,ZI,
     +                  IERR)
C
*
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
*     ------------------------ END TIMING CODE -------------------------
*
C
C
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,
C     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
C     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING
C     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS
C     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS
C     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY
C     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX
C          WITH REAL SUBDIAGONAL ELEMENTS,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE,
C
C        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET
C          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO
C          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,
C
C        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE
C          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN
C          THE RATIOS ((ALFR+I*ALFI)/BETA),
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF AR(J,J-1) HAS NOT BECOME
C                     ZERO AFTER 50 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      REAL EPS1
      INTEGER IERR,N,NM
      LOGICAL MATZ
C     ..
C     .. Array Arguments ..
      REAL AI(NM,N),ALFI(N),ALFR(N),AR(NM,N),BETA(N),BI(NM,N),BR(NM,N),
     +     ZI(NM,N),ZR(NM,N)
C     ..
C     .. Local Scalars ..
      COMPLEX Z3
      REAL A1,A1I,A2,A33,A33I,A34,A34I,A43,A43I,A44,A44I,ANI,ANORM,B11,
     +     B33,B3344,B3344I,B33I,B44,B44I,BNI,BNORM,EP,EPSA,EPSB,OPST,R,
     +     S,SH,SHI,U1,U1I,U2,XI,XR,YI,YR
      INTEGER EN,ENM2,ENORN,I,IOPST,ITS,J,K,K1,K2,KM1,L,L1,LL,LM1,LOR1,
     +        NA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,CABS,CMPLX,CSQRT,MAX0,REAL,SQRT
C     ..
      IERR = 0
C     ********** COMPUTE EPSA,EPSB **********
      ANORM = 0.0
      BNORM = 0.0
C
      DO 20 I = 1,N
          ANI = 0.0
          IF (I.NE.1) ANI = ABS(AR(I,I-1))
          BNI = 0.0
C
          DO 10 J = I,N
              ANI = ANI + ABS(AR(I,J)) + ABS(AI(I,J))
              BNI = BNI + ABS(BR(I,J)) + ABS(BI(I,J))
   10     CONTINUE
C
          IF (ANI.GT.ANORM) ANORM = ANI
          IF (BNI.GT.BNORM) BNORM = BNI
   20 CONTINUE
C
      IF (ANORM.EQ.0.0) ANORM = 1.0
      IF (BNORM.EQ.0.0) BNORM = 1.0
      EP = EPS1
      IF (EP.GT.0.0) GO TO 40
C     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
      EP = 1.0
   30 EP = EP/2.0
      IF (1.0+EP.GT.1.0) GO TO 30
   40 EPSA = EP*ANORM
      EPSB = EP*BNORM
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COUNT OPS FOR NORMS, BUT NOT FOR CALCULATION OF "EP"
      OPS = OPS + REAL(2*N* (N+1)+2)
      OPST = 0.0
      ITCNT = 0.0
*     ------------------------ END TIMING CODE -------------------------
C     ********** REDUCE A TO TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR **********
      LOR1 = 1
      ENORN = N
      EN = N
C     ********** BEGIN QZ STEP **********
   50 IF (EN.EQ.0) GO TO 300
      IF (.NOT.MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
   60 CONTINUE
C     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- **********
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPS = OPS + OPST
      OPST = 0.0
*     ------------------------ END TIMING CODE -------------------------
      DO 70 LL = 1,EN
          LM1 = EN - LL
          L = LM1 + 1
          IF (L.EQ.1) GO TO 90
          IF (ABS(AR(L,LM1)).LE.EPSA) GO TO 80
   70 CONTINUE
C
   80 AR(L,LM1) = 0.0
C     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********
   90 B11 = CABS(CMPLX(BR(L,L),BI(L,L)))
      IF (B11.EQ.0.0) GO TO 110
      U1 = BR(L,L)/B11
      U1I = BI(L,L)/B11
C
      DO 100 J = L,ENORN
          XI = U1*AI(L,J) - U1I*AR(L,J)
          AR(L,J) = U1*AR(L,J) + U1I*AI(L,J)
          AI(L,J) = XI
          XI = U1*BI(L,J) - U1I*BR(L,J)
          BR(L,J) = U1*BR(L,J) + U1I*BI(L,J)
          BI(L,J) = XI
  100 CONTINUE
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + REAL(7+12* (ENORN+1-L))
*     ------------------------ END TIMING CODE -------------------------
C
      BI(L,L) = 0.0
  110 IF (L.NE.EN) GO TO 120
C     ********** 1-BY-1 BLOCK ISOLATED **********
      ALFR(EN) = AR(EN,EN)
      ALFI(EN) = AI(EN,EN)
      BETA(EN) = B11
      EN = NA
      GO TO 50
C     ********** CHECK FOR SMALL TOP OF B **********
  120 L1 = L + 1
      IF (B11.GT.EPSB) GO TO 140
      BR(L,L) = 0.0
      S = ABS(AR(L,L)) + ABS(AI(L,L)) + ABS(AR(L1,L))
      U1 = AR(L,L)/S
      U1I = AI(L,L)/S
      U2 = AR(L1,L)/S
      R = SQRT(U1*U1+U1I*U1I+U2*U2)
      U1 = U1/R
      U1I = U1I/R
      U2 = U2/R
      AR(L,L) = R*S
      AI(L,L) = 0.0
C
      DO 130 J = L1,ENORN
          XR = AR(L,J)
          XI = AI(L,J)
          YR = AR(L1,J)
          YI = AI(L1,J)
          AR(L,J) = U1*XR + U1I*XI + U2*YR
          AI(L,J) = U1*XI - U1I*XR + U2*YI
          AR(L1,J) = U1*YR - U1I*YI - U2*XR
          AI(L1,J) = U1*YI + U1I*YR - U2*XI
          XR = BR(L,J)
          XI = BI(L,J)
          YR = BR(L1,J)
          YI = BI(L1,J)
          BR(L1,J) = U1*YR - U1I*YI - U2*XR
          BR(L,J) = U1*XR + U1I*XI + U2*YR
          BI(L,J) = U1*XI - U1I*XR + U2*YI
          BI(L1,J) = U1*YI + U1I*YR - U2*XI
  130 CONTINUE
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + REAL(15+40* (ENORN-L))
*     ------------------------ END TIMING CODE -------------------------
C
      LM1 = L
      L = L1
      GO TO 80
C     ********** ITERATION STRATEGY **********
  140 IF (ITS.EQ.50) GO TO 290
      IF (ITS.EQ.10) GO TO 180
C     ********** DETERMINE SHIFT **********
      B33 = BR(NA,NA)
      B33I = BI(NA,NA)
      IF (CABS(CMPLX(B33,B33I)).GE.EPSB) GO TO 150
      B33 = EPSB
      B33I = 0.0
  150 B44 = BR(EN,EN)
      B44I = BI(EN,EN)
      IF (CABS(CMPLX(B44,B44I)).GE.EPSB) GO TO 160
      B44 = EPSB
      B44I = 0.0
  160 B3344 = B33*B44 - B33I*B44I
      B3344I = B33*B44I + B33I*B44
      A33 = AR(NA,NA)*B44 - AI(NA,NA)*B44I
      A33I = AR(NA,NA)*B44I + AI(NA,NA)*B44
      A34 = AR(NA,EN)*B33 - AI(NA,EN)*B33I - AR(NA,NA)*BR(NA,EN) +
     +      AI(NA,NA)*BI(NA,EN)
      A34I = AR(NA,EN)*B33I + AI(NA,EN)*B33 - AR(NA,NA)*BI(NA,EN) -
     +       AI(NA,NA)*BR(NA,EN)
      A43 = AR(EN,NA)*B44
      A43I = AR(EN,NA)*B44I
      A44 = AR(EN,EN)*B33 - AI(EN,EN)*B33I - AR(EN,NA)*BR(NA,EN)
      A44I = AR(EN,EN)*B33I + AI(EN,EN)*B33 - AR(EN,NA)*BI(NA,EN)
      SH = A44
      SHI = A44I
      XR = A34*A43 - A34I*A43I
      XI = A34*A43I + A34I*A43
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + REAL(54)
*     ------------------------ END TIMING CODE -------------------------
      IF (XR.EQ.0.0 .AND. XI.EQ.0.0) GO TO 190
      YR = (A33-SH)/2.0
      YI = (A33I-SHI)/2.0
      Z3 = CSQRT(CMPLX(YR**2-YI**2+XR,2.0*YR*YI+XI))
      U1 = REAL(Z3)
      U1I = AIMAG(Z3)
      IF (YR*U1+YI*U1I.GE.0.0) GO TO 170
      U1 = -U1
      U1I = -U1I
  170 Z3 = (CMPLX(SH,SHI)-CMPLX(XR,XI)/CMPLX(YR+U1,YI+U1I))/
     +     CMPLX(B3344,B3344I)
      SH = REAL(Z3)
      SHI = AIMAG(Z3)
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + REAL(66)
*     ------------------------ END TIMING CODE -------------------------
      GO TO 190
C     ********** AD HOC SHIFT **********
  180 SH = AR(EN,NA) + AR(NA,ENM2)
      SHI = 0.0
C     ********** DETERMINE ZEROTH COLUMN OF A **********
  190 A1 = AR(L,L)/B11 - SH
      A1I = AI(L,L)/B11 - SHI
      A2 = AR(L1,L)/B11
      ITS = ITS + 1
*     ----------------------- BEGIN TIMING CODE ------------------------
      ITCNT = ITCNT + 1.0
*     ------------------------ END TIMING CODE -------------------------
      IF (.NOT.MATZ) LOR1 = L
C     ********** MAIN LOOP **********
      DO 270 K = L,NA
          K1 = K + 1
          K2 = K + 2
          KM1 = MAX0(K-1,L)
C     ********** ZERO A(K+1,K-1) **********
          IF (K.EQ.L) GO TO 200
          A1 = AR(K,KM1)
          A1I = AI(K,KM1)
          A2 = AR(K1,KM1)
  200     S = ABS(A1) + ABS(A1I) + ABS(A2)
          U1 = A1/S
          U1I = A1I/S
          U2 = A2/S
          R = SQRT(U1*U1+U1I*U1I+U2*U2)
          U1 = U1/R
          U1I = U1I/R
          U2 = U2/R
C
          DO 210 J = KM1,ENORN
              XR = AR(K,J)
              XI = AI(K,J)
              YR = AR(K1,J)
              YI = AI(K1,J)
              AR(K,J) = U1*XR + U1I*XI + U2*YR
              AI(K,J) = U1*XI - U1I*XR + U2*YI
              AR(K1,J) = U1*YR - U1I*YI - U2*XR
              AI(K1,J) = U1*YI + U1I*YR - U2*XI
              XR = BR(K,J)
              XI = BI(K,J)
              YR = BR(K1,J)
              YI = BI(K1,J)
              BR(K,J) = U1*XR + U1I*XI + U2*YR
              BI(K,J) = U1*XI - U1I*XR + U2*YI
              BR(K1,J) = U1*YR - U1I*YI - U2*XR
              BI(K1,J) = U1*YI + U1I*YR - U2*XI
  210     CONTINUE
C
          IF (K.EQ.L) GO TO 220
          AI(K,KM1) = 0.0
          AR(K1,KM1) = 0.0
          AI(K1,KM1) = 0.0
C     ********** ZERO B(K+1,K) **********
  220     S = ABS(BR(K1,K1)) + ABS(BI(K1,K1)) + ABS(BR(K1,K))
          U1 = BR(K1,K1)/S
          U1I = BI(K1,K1)/S
          U2 = BR(K1,K)/S
          R = SQRT(U1*U1+U1I*U1I+U2*U2)
          U1 = U1/R
          U1I = U1I/R
          U2 = U2/R
          IF (K.EQ.NA) GO TO 230
          XR = AR(K2,K1)
          AR(K2,K1) = U1*XR
          AI(K2,K1) = -U1I*XR
          AR(K2,K) = -U2*XR
C*PL*ERROR* Embedded comment after label moved
C
  230     DO 240 I = LOR1,K1
              XR = AR(I,K1)
              XI = AI(I,K1)
              YR = AR(I,K)
              YI = AI(I,K)
              AR(I,K1) = U1*XR + U1I*XI + U2*YR
              AI(I,K1) = U1*XI - U1I*XR + U2*YI
              AR(I,K) = U1*YR - U1I*YI - U2*XR
              AI(I,K) = U1*YI + U1I*YR - U2*XI
              XR = BR(I,K1)
              XI = BI(I,K1)
              YR = BR(I,K)
              YI = BI(I,K)
              BR(I,K1) = U1*XR + U1I*XI + U2*YR
              BI(I,K1) = U1*XI - U1I*XR + U2*YI
              BR(I,K) = U1*YR - U1I*YI - U2*XR
              BI(I,K) = U1*YI + U1I*YR - U2*XI
  240     CONTINUE
C
          BI(K1,K1) = 0.0
          BR(K1,K) = 0.0
          BI(K1,K) = 0.0
          IF (.NOT.MATZ) GO TO 260
C
          DO 250 I = 1,N
              XR = ZR(I,K1)
              XI = ZI(I,K1)
              YR = ZR(I,K)
              YI = ZI(I,K)
              ZR(I,K1) = U1*XR + U1I*XI + U2*YR
              ZI(I,K1) = U1*XI - U1I*XR + U2*YI
              ZR(I,K) = U1*YR - U1I*YI - U2*XR
              ZI(I,K) = U1*YI + U1I*YR - U2*XI
  250     CONTINUE
  260     CONTINUE
  270 CONTINUE
C
*
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COUNT OPS FOR STATEMENTS 140 -- 260
      IOPST = 29 + 40* (ENORN-LOR1+4)
      IF (MATZ) IOPST = IOPST + 20*N
      OPST = OPST + (REAL(N-L)*REAL(IOPST)+2)
      IF (L.LE.1) OPST = OPST - 40
*     ------------------------ END TIMING CODE -------------------------
*
C     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP **********
      IF (AI(EN,NA).EQ.0.0) GO TO 60
      R = CABS(CMPLX(AR(EN,NA),AI(EN,NA)))
      U1 = AR(EN,NA)/R
      U1I = AI(EN,NA)/R
      AR(EN,NA) = R
      AI(EN,NA) = 0.0
C
      DO 280 J = EN,ENORN
          XI = U1*AI(EN,J) - U1I*AR(EN,J)
          AR(EN,J) = U1*AR(EN,J) + U1I*AI(EN,J)
          AI(EN,J) = XI
          XI = U1*BI(EN,J) - U1I*BR(EN,J)
          BR(EN,J) = U1*BR(EN,J) + U1I*BI(EN,J)
          BI(EN,J) = XI
  280 CONTINUE
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPST = OPST + REAL(7+12* (EN+1-ENORN))
*     ------------------------ END TIMING CODE -------------------------
C
      GO TO 60
C     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT
C                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
  290 IERR = EN
C     ********** SAVE EPSB FOR USE BY CQZVEC **********
  300 IF (N.GT.1) BR(N,1) = EPSB
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPS = OPS + OPST
      OPST = 0.0
*     ------------------------ END TIMING CODE -------------------------
      RETURN
C     ********** LAST CARD OF CQZVAL **********
      END
      SUBROUTINE CQZVEC(NM,N,AR,AI,BR,BI,ALFR,ALFI,BETA,ZR,ZI)
C
C
C
*
*     ----------------------- BEGIN TIMING CODE ------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      REAL ITCNT,OPS
C     ..
*     ------------------------ END TIMING CODE -------------------------
*
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER
C     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
C     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
C     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL
C          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS
C          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL,
C
C        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,
C
C        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT-
C
C        A IS UNALTERED,
C
C        B HAS BEEN DESTROYED,
C
C        ALFR, ALFI, AND BETA ARE UNALTERED,
C
C        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED
C          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER N,NM
C     ..
C     .. Array Arguments ..
      REAL AI(NM,N),ALFI(N),ALFR(N),AR(NM,N),BETA(N),BI(NM,N),BR(NM,N),
     +     ZI(NM,N),ZR(NM,N)
C     ..
C     .. Local Scalars ..
      COMPLEX Z3
      REAL ALMI,ALMR,BETM,EPSB,R,RI,T,TI,XI
      INTEGER EN,I,II,J,JJ,K,M,NA,NN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC AIMAG,CABS,CMPLX,REAL
C     ..
      IF (N.GT.1) THEN
          EPSB = BR(N,1)
C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
          DO 30 NN = 2,N
              EN = N + 2 - NN
              NA = EN - 1
              ALMR = ALFR(EN)
              ALMI = ALFI(EN)
              BETM = BETA(EN)
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
              DO 20 II = 1,NA
                  I = EN - II
                  R = 0.0
                  RI = 0.0
                  M = I + 1
C
                  DO 10 J = M,EN
                      T = BETM*AR(I,J) - ALMR*BR(I,J) + ALMI*BI(I,J)
                      TI = BETM*AI(I,J) - ALMR*BI(I,J) - ALMI*BR(I,J)
                      IF (J.NE.EN) THEN
                          XI = T*BI(J,EN) + TI*BR(J,EN)
                          T = T*BR(J,EN) - TI*BI(J,EN)
                          TI = XI
                      END IF
                      R = R + T
                      RI = RI + TI
   10             CONTINUE
C
                  T = ALMR*BETA(I) - BETM*ALFR(I)
                  TI = ALMI*BETA(I) - BETM*ALFI(I)
                  IF (T.EQ.0.0 .AND. TI.EQ.0.0) T = EPSB
                  Z3 = CMPLX(R,RI)/CMPLX(T,TI)
                  BR(I,EN) = REAL(Z3)
                  BI(I,EN) = AIMAG(Z3)
   20         CONTINUE
   30     CONTINUE
C
C     ********** END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 2 DO -- **********
          DO 60 JJ = 2,N
              J = N + 2 - JJ
              M = J - 1
C
              DO 50 I = 1,N
C
                  DO 40 K = 1,M
                      ZR(I,J) = ZR(I,J) + ZR(I,K)*BR(K,J) -
     +                          ZI(I,K)*BI(K,J)
                      ZI(I,J) = ZI(I,J) + ZR(I,K)*BI(K,J) +
     +                          ZI(I,K)*BR(K,J)
   40             CONTINUE
   50         CONTINUE
   60     CONTINUE
C
C     ********** NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1 **********
          DO 90 J = 1,N
              T = 0.0
C
              DO 70 I = 1,N
                  R = CABS(CMPLX(ZR(I,J),ZI(I,J)))
                  IF (R.GT.T) T = R
   70         CONTINUE
C
              DO 80 I = 1,N
                  ZR(I,J) = ZR(I,J)/T
                  ZI(I,J) = ZI(I,J)/T
   80         CONTINUE
   90     CONTINUE
C
      END IF
C
*
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPS = OPS + REAL(N)*REAL(14*N**2+15*N-15)/REAL(2)
*     ------------------------ END TIMING CODE -------------------------
*
C     ********** LAST CARD OF CQZVEC **********
      END
