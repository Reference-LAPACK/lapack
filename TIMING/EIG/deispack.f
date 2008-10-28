      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
C
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C

C     .. Scalar Arguments ..
      DOUBLE PRECISION AI,AR,BI,BR,CI,CR
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AIS,ARS,BIS,BRS,S
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
      S = DABS(BR) + DABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS+AIS*BIS)/S
      CI = (AIS*BRS-ARS*BIS)/S
      END
      DOUBLE PRECISION FUNCTION EPSLON(X)
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
      DOUBLE PRECISION X
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,C,EPS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
      A = 4.0D0/3.0D0
   10 CONTINUE
      B = A - 1.0D0
      C = B + B + B
      EPS = DABS(C-1.0D0)
      IF (EPS.EQ.0.0D0) GO TO 10
      EPSLON = EPS*DABS(X)
      END
      SUBROUTINE HQR(NM,N,LOW,IGH,H,WR,WI,IERR)
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
      DOUBLE PRECISION ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR,
C     NUM. MATH. 14, 219-231(1970) BY MARTIN, PETERS, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A REAL
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
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE TRANSFORMATIONS USED IN THE REDUCTION TO HESSENBERG
C          FORM BY  ELMHES  OR  ORTHES, IF PERFORMED, IS STORED
C          IN THE REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.  THEREFORE, IT MUST BE SAVED
C          BEFORE CALLING  HQR  IF SUBSEQUENT CALCULATION AND
C          BACK TRANSFORMATION OF EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C     MODIFIED ON 11/1/89; ADJUSTING INDICES OF LOOPS
C       200, 210, 230, AND 240 TO INCREASE PERFORMANCE. JACK DONGARRA
C
C     ------------------------------------------------------------------
C
*

C     .. Scalar Arguments ..
      INTEGER IERR,IGH,LOW,N,NM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION H(NM,N),WI(N),WR(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION NORM,OPST,OVFL,P,Q,R,S,SMALL,SMLNUM,T,TST1,TST2,
     +                 ULP,UNFL,W,X,Y,ZZ
      INTEGER EN,ENM2,I,ITN,ITS,J,K,L,LL,M,MM,MP2,NA
      LOGICAL NOTLAS
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DSIGN,DSQRT,MAX,MIN,MIN0
C     ..
      IF (N.GT.0) THEN
*
*
*     INITIALIZE
          ITCNT = 0
          OPST = 0
          IERR = 0
          K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
          DO 10 I = 1,N
              K = I
              IF (I.LT.LOW .OR. I.GT.IGH) THEN
                  WR(I) = H(I,I)
                  WI(I) = 0.0D0
              END IF
   10     CONTINUE
*
*        INCREMENT OPCOUNT FOR COMPUTING MATRIX NORM
          OPS = OPS + (IGH-LOW+1)* (IGH-LOW+2)/2
*
*     COMPUTE THE 1-NORM OF MATRIX H
*
          NORM = 0.0D0
          DO 30 J = LOW,IGH
              S = 0.0D0
              DO 20 I = LOW,MIN(IGH,J+1)
                  S = S + DABS(H(I,J))
   20         CONTINUE
              NORM = MAX(NORM,S)
   30     CONTINUE
*
          UNFL = DLAMCH('SAFE MINIMUM')
          OVFL = DLAMCH('OVERFLOW')
          ULP = DLAMCH('EPSILON')*DLAMCH('BASE')
          SMLNUM = MAX(UNFL* (N/ULP),N/ (ULP*OVFL))
          SMALL = MAX(SMLNUM,ULP*NORM)
C
          EN = IGH
          T = 0.0D0
          ITN = 30*N
   40     CONTINUE
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
          IF (EN.LT.LOW) THEN
              GO TO 190
          ELSE
              ITS = 0
              NA = EN - 1
              ENM2 = NA - 1
   50         CONTINUE
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
*     REPLACE SPLITTING CRITERION WITH NEW ONE AS IN LAPACK
*
              DO 60 LL = LOW,EN
                  L = EN + LOW - LL
                  IF (L.EQ.LOW) THEN
                      GO TO 70
                  ELSE
                      S = DABS(H(L-1,L-1)) + DABS(H(L,L))
                      IF (S.EQ.0.0D0) S = NORM
                      IF (DABS(H(L,L-1)).LE.MAX(ULP*S,SMALL))
     +                    GO TO 70
                  END IF
   60         CONTINUE
C     .......... FORM SHIFT ..........
   70         CONTINUE
*
*        INCREMENT OP COUNT FOR CONVERGENCE TEST
              OPS = OPS + 2* (EN-L+1)
              X = H(EN,EN)
              IF (L.EQ.EN) THEN
                  GO TO 170
              ELSE
                  Y = H(NA,NA)
                  W = H(EN,NA)*H(NA,EN)
                  IF (L.NE.NA) THEN
                      IF (ITN.EQ.0) THEN
                          GO TO 180
                      ELSE
                          IF (ITS.EQ.10 .OR. ITS.EQ.20) THEN
C     .......... FORM EXCEPTIONAL SHIFT ..........
*
*        INCREMENT OP COUNT FOR FORMING EXCEPTIONAL SHIFT
                              OPS = OPS + (EN-LOW+6)
                              T = T + X
C
                              DO 80 I = LOW,EN
                                  H(I,I) = H(I,I) - X
   80                         CONTINUE
C
                              S = DABS(H(EN,NA)) + DABS(H(NA,ENM2))
                              X = 0.75D0*S
                              Y = X
                              W = -0.4375D0*S*S
                          END IF
                          ITS = ITS + 1
                          ITN = ITN - 1
*
*       UPDATE ITERATION NUMBER
                          ITCNT = 30*N - ITN
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
*     REPLACE SPLITTING CRITERION WITH NEW ONE AS IN LAPACK
                          DO 90 MM = L,ENM2
                              M = ENM2 + L - MM
                              ZZ = H(M,M)
                              R = X - ZZ
                              S = Y - ZZ
                              P = (R*S-W)/H(M+1,M) + H(M,M+1)
                              Q = H(M+1,M+1) - ZZ - R - S
                              R = H(M+2,M+1)
                              S = DABS(P) + DABS(Q) + DABS(R)
                              P = P/S
                              Q = Q/S
                              R = R/S
                              IF (M.EQ.L) THEN
                                  GO TO 100
                              ELSE
                                  TST1 = DABS(P)* (DABS(H(M-1,M-1))+
     +                                   DABS(ZZ)+DABS(H(M+1,M+1)))
                                  TST2 = DABS(H(M,M-1))*
     +                                   (DABS(Q)+DABS(R))
                                  IF (TST2.LE.MAX(ULP*TST1,SMALL))
     +                                GO TO 100
                              END IF
   90                     CONTINUE
C
  100                     CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 140
                          OPST = OPST + 20* (ENM2-M+1)
                          MP2 = M + 2
C
                          DO 110 I = MP2,EN
                              H(I,I-2) = 0.0D0
                              IF (I.NE.MP2) H(I,I-3) = 0.0D0
  110                     CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
*
*        INCREMENT OPCOUNT FOR LOOP 260
                          OPST = OPST + 18* (NA-M+1)
                          DO 160 K = M,NA
                              NOTLAS = K .NE. NA
                              IF (K.NE.M) THEN
                                  P = H(K,K-1)
                                  Q = H(K+1,K-1)
                                  R = 0.0D0
                                  IF (NOTLAS) R = H(K+2,K-1)
                                  X = DABS(P) + DABS(Q) + DABS(R)
                                  IF (X.EQ.0.0D0) THEN
                                      GO TO 160
                                  ELSE
                                      P = P/X
                                      Q = Q/X
                                      R = R/X
                                  END IF
                              END IF
                              S = DSIGN(DSQRT(P*P+Q*Q+R*R),P)
                              IF (K.NE.M) THEN
                                  H(K,K-1) = -S*X
                              ELSE IF (L.NE.M) THEN
                                  H(K,K-1) = -H(K,K-1)
                              END IF
                              P = P + S
                              X = P/S
                              Y = Q/S
                              ZZ = R/S
                              Q = Q/P
                              R = R/P
                              IF (NOTLAS) THEN
C     .......... ROW MODIFICATION ..........
*
*        INCREMENT OPCOUNT
                                  OPS = OPS + 10* (EN-K+1)
                                  DO 120 J = K,EN
                                      P = H(K,J) + Q*H(K+1,J) +
     +                                    R*H(K+2,J)
                                      H(K,J) = H(K,J) - P*X
                                      H(K+1,J) = H(K+1,J) - P*Y
                                      H(K+2,J) = H(K+2,J) - P*ZZ
  120                             CONTINUE
C
                                  J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
*
*        INCREMENT OPCOUNT
                                  OPS = OPS + 10* (J-L+1)
                                  DO 130 I = L,J
                                      P = X*H(I,K) + Y*H(I,K+1) +
     +                                    ZZ*H(I,K+2)
                                      H(I,K) = H(I,K) - P
                                      H(I,K+1) = H(I,K+1) - P*Q
                                      H(I,K+2) = H(I,K+2) - P*R
  130                             CONTINUE
                              ELSE
C     .......... ROW MODIFICATION ..........
*
*        INCREMENT OPCOUNT
                                  OPS = OPS + 6* (EN-K+1)
                                  DO 140 J = K,EN
                                      P = H(K,J) + Q*H(K+1,J)
                                      H(K,J) = H(K,J) - P*X
                                      H(K+1,J) = H(K+1,J) - P*Y
  140                             CONTINUE
C
                                  J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
*
*        INCREMENT OPCOUNT
                                  OPS = OPS + 6* (J-L+1)
                                  DO 150 I = L,J
                                      P = X*H(I,K) + Y*H(I,K+1)
                                      H(I,K) = H(I,K) - P
                                      H(I,K+1) = H(I,K+1) - P*Q
  150                             CONTINUE
                              END IF
  160                     CONTINUE
C
C
                          GO TO 50
                      END IF
                  END IF
              END IF
C     .......... TWO ROOTS FOUND ..........
              P = (Y-X)/2.0D0
              Q = P*P + W
              ZZ = DSQRT(DABS(Q))
              X = X + T
*
*        INCREMENT OP COUNT FOR FINDING TWO ROOTS.
              OPST = OPST + 8
              IF (Q.LT.0.0D0) THEN
C     .......... COMPLEX PAIR ..........
                  WR(NA) = X + P
                  WR(EN) = X + P
                  WI(NA) = ZZ
                  WI(EN) = -ZZ
              ELSE
C     .......... REAL PAIR ..........
                  ZZ = P + DSIGN(ZZ,P)
                  WR(NA) = X + ZZ
                  WR(EN) = WR(NA)
                  IF (ZZ.NE.0.0D0) WR(EN) = X - W/ZZ
                  WI(NA) = 0.0D0
                  WI(EN) = 0.0D0
              END IF
              EN = ENM2
              GO TO 40
C     .......... ONE ROOT FOUND ..........
  170         WR(EN) = X + T
              WI(EN) = 0.0D0
              EN = NA
              GO TO 40
          END IF
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
  180     IERR = EN
*
*     COMPUTE FINAL OP COUNT
  190     OPS = OPS + OPST
      END IF
      END
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
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
      DOUBLE PRECISION ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
C     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
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
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
C          IDENTITY MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
*

C     .. Scalar Arguments ..
      INTEGER IERR,IGH,LOW,N,NM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION H(NM,N),WI(N),WR(N),Z(NM,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION NORM,OPST,OVFL,P,Q,R,RA,S,SA,SMALL,SMLNUM,T,TST1,
     +                 TST2,ULP,UNFL,VI,VR,W,X,Y,ZZ
      INTEGER EN,ENM2,I,II,ITN,ITS,J,JJ,K,L,LL,M,MM,MP2,NA,NN
      LOGICAL NOTLAS
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
C     ..
C     .. External Subroutines ..
      EXTERNAL CDIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DABS,DMAX1,DSIGN,DSQRT,MAX,MIN,MIN0
C     ..
      IF (N.GT.0) THEN
*
*     INITIALIZE
*
          ITCNT = 0
          OPST = 0
C
          IERR = 0
          K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
          DO 10 I = 1,N
              IF (I.LT.LOW .OR. I.GT.IGH) THEN
                  WR(I) = H(I,I)
                  WI(I) = 0.0D0
              END IF
   10     CONTINUE

*        INCREMENT OPCOUNT FOR COMPUTING MATRIX NORM
          OPS = OPS + (IGH-LOW+1)* (IGH-LOW+2)/2
*
*     COMPUTE THE 1-NORM OF MATRIX H
*
          NORM = 0.0D0
          DO 30 J = LOW,IGH
              S = 0.0D0
              DO 20 I = LOW,MIN(IGH,J+1)
                  S = S + DABS(H(I,J))
   20         CONTINUE
              NORM = MAX(NORM,S)
   30     CONTINUE
C
          UNFL = DLAMCH('SAFE MINIMUM')
          OVFL = DLAMCH('OVERFLOW')
          ULP = DLAMCH('EPSILON')*DLAMCH('BASE')
          SMLNUM = MAX(UNFL* (N/ULP),N/ (ULP*OVFL))
          SMALL = MAX(SMLNUM,MIN((NORM*SMLNUM)*NORM,ULP*NORM))
C
          EN = IGH
          T = 0.0D0
          ITN = 30*N
   40     CONTINUE
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
          IF (EN.LT.LOW) THEN
              GO TO 240
          ELSE
              ITS = 0
              NA = EN - 1
              ENM2 = NA - 1
   50         CONTINUE
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
*     REPLACE SPLITTING CRITERION WITH NEW ONE AS IN LAPACK
*
              DO 60 LL = LOW,EN
                  L = EN + LOW - LL
                  IF (L.EQ.LOW) THEN
                      GO TO 70
                  ELSE
                      S = DABS(H(L-1,L-1)) + DABS(H(L,L))
                      IF (S.EQ.0.0D0) S = NORM
                      IF (ABS(H(L,L-1)).LE.MAX(ULP*S,SMALL)) GO TO 70
                  END IF
   60         CONTINUE
C     .......... FORM SHIFT ..........
   70         CONTINUE
*
*        INCREMENT OP COUNT FOR CONVERGENCE TEST
              OPS = OPS + 2* (EN-L+1)
              X = H(EN,EN)
              IF (L.EQ.EN) THEN
                  GO TO 220
              ELSE
                  Y = H(NA,NA)
                  W = H(EN,NA)*H(NA,EN)
                  IF (L.NE.NA) THEN
                      IF (ITN.EQ.0) THEN
                          GO TO 230
                      ELSE
                          IF (ITS.EQ.10 .OR. ITS.EQ.20) THEN
C     .......... FORM EXCEPTIONAL SHIFT ..........
*
*        INCREMENT OP COUNT
                              OPS = OPS + (EN-LOW+6)
                              T = T + X
C
                              DO 80 I = LOW,EN
                                  H(I,I) = H(I,I) - X
   80                         CONTINUE
C
                              S = DABS(H(EN,NA)) + DABS(H(NA,ENM2))
                              X = 0.75D0*S
                              Y = X
                              W = -0.4375D0*S*S
                          END IF
                          ITS = ITS + 1
                          ITN = ITN - 1
*
*       UPDATE ITERATION NUMBER
                          ITCNT = 30*N - ITN
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
                          DO 90 MM = L,ENM2
                              M = ENM2 + L - MM
                              ZZ = H(M,M)
                              R = X - ZZ
                              S = Y - ZZ
                              P = (R*S-W)/H(M+1,M) + H(M,M+1)
                              Q = H(M+1,M+1) - ZZ - R - S
                              R = H(M+2,M+1)
                              S = DABS(P) + DABS(Q) + DABS(R)
                              P = P/S
                              Q = Q/S
                              R = R/S
                              IF (M.EQ.L) THEN
                                  GO TO 100
                              ELSE
                                  TST1 = DABS(P)* (DABS(H(M-1,M-1))+
     +                                   DABS(ZZ)+DABS(H(M+1,M+1)))
                                  TST2 = DABS(H(M,M-1))*
     +                                   (DABS(Q)+DABS(R))
                                  IF (TST2.LE.MAX(ULP*TST1,SMALL))
     +                                GO TO 100
                              END IF
   90                     CONTINUE
C
  100                     CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 140
                          OPST = OPST + 20* (ENM2-M+1)
                          MP2 = M + 2
C
                          DO 110 I = MP2,EN
                              H(I,I-2) = 0.0D0
                              IF (I.NE.MP2) H(I,I-3) = 0.0D0
  110                     CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
*
*        INCREMENT OPCOUNT FOR LOOP 260
                          OPST = OPST + 18* (NA-M+1)
                          DO 180 K = M,NA
                              NOTLAS = K .NE. NA
                              IF (K.NE.M) THEN
                                  P = H(K,K-1)
                                  Q = H(K+1,K-1)
                                  R = 0.0D0
                                  IF (NOTLAS) R = H(K+2,K-1)
                                  X = DABS(P) + DABS(Q) + DABS(R)
                                  IF (X.EQ.0.0D0) THEN
                                      GO TO 180
                                  ELSE
                                      P = P/X
                                      Q = Q/X
                                      R = R/X
                                  END IF
                              END IF
                              S = DSIGN(DSQRT(P*P+Q*Q+R*R),P)
                              IF (K.NE.M) THEN
                                  H(K,K-1) = -S*X
                              ELSE IF (L.NE.M) THEN
                                  H(K,K-1) = -H(K,K-1)
                              END IF
                              P = P + S
                              X = P/S
                              Y = Q/S
                              ZZ = R/S
                              Q = Q/P
                              R = R/P
                              IF (NOTLAS) THEN
C     .......... ROW MODIFICATION ..........
*
*        INCREMENT OPCOUNT FOR LOOP 230
                                  OPS = OPS + 10* (N-K+1)
                                  DO 120 J = K,N
                                      P = H(K,J) + Q*H(K+1,J) +
     +                                    R*H(K+2,J)
                                      H(K,J) = H(K,J) - P*X
                                      H(K+1,J) = H(K+1,J) - P*Y
                                      H(K+2,J) = H(K+2,J) - P*ZZ
  120                             CONTINUE
C
                                  J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
*
*        INCREMENT OPCOUNT FOR LOOP 240
                                  OPS = OPS + 10*J
                                  DO 130 I = 1,J
                                      P = X*H(I,K) + Y*H(I,K+1) +
     +                                    ZZ*H(I,K+2)
                                      H(I,K) = H(I,K) - P
                                      H(I,K+1) = H(I,K+1) - P*Q
                                      H(I,K+2) = H(I,K+2) - P*R
  130                             CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
*
*        INCREMENT OPCOUNT FOR LOOP 250
                                  OPS = OPS + 10* (IGH-LOW+1)
                                  DO 140 I = LOW,IGH
                                      P = X*Z(I,K) + Y*Z(I,K+1) +
     +                                    ZZ*Z(I,K+2)
                                      Z(I,K) = Z(I,K) - P
                                      Z(I,K+1) = Z(I,K+1) - P*Q
                                      Z(I,K+2) = Z(I,K+2) - P*R
  140                             CONTINUE
                              ELSE
C     .......... ROW MODIFICATION ..........
*
*        INCREMENT OP COUNT FOR LOOP 200
                                  OPS = OPS + 6* (N-K+1)
                                  DO 150 J = K,N
                                      P = H(K,J) + Q*H(K+1,J)
                                      H(K,J) = H(K,J) - P*X
                                      H(K+1,J) = H(K+1,J) - P*Y
  150                             CONTINUE
C
                                  J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
*
*        INCREMENT OPCOUNT FOR LOOP 210
                                  OPS = OPS + 6*J
                                  DO 160 I = 1,J
                                      P = X*H(I,K) + Y*H(I,K+1)
                                      H(I,K) = H(I,K) - P
                                      H(I,K+1) = H(I,K+1) - P*Q
  160                             CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
*
*        INCREMENT OPCOUNT FOR LOOP 220
                                  OPS = OPS + 6* (IGH-LOW+1)
                                  DO 170 I = LOW,IGH
                                      P = X*Z(I,K) + Y*Z(I,K+1)
                                      Z(I,K) = Z(I,K) - P
                                      Z(I,K+1) = Z(I,K+1) - P*Q
  170                             CONTINUE
                              END IF
  180                     CONTINUE
C
C
                          GO TO 50
                      END IF
                  END IF
              END IF
C     .......... TWO ROOTS FOUND ..........
              P = (Y-X)/2.0D0
              Q = P*P + W
              ZZ = DSQRT(DABS(Q))
              H(EN,EN) = X + T
              X = H(EN,EN)
              H(NA,NA) = Y + T
              IF (Q.LT.0.0D0) THEN
C     .......... COMPLEX PAIR ..........
                  WR(NA) = X + P
                  WR(EN) = X + P
                  WI(NA) = ZZ
                  WI(EN) = -ZZ
*
*        INCREMENT OP COUNT FOR FINDING COMPLEX PAIR.
                  OPST = OPST + 9
              ELSE
C     .......... REAL PAIR ..........
                  ZZ = P + DSIGN(ZZ,P)
                  WR(NA) = X + ZZ
                  WR(EN) = WR(NA)
                  IF (ZZ.NE.0.0D0) WR(EN) = X - W/ZZ
                  WI(NA) = 0.0D0
                  WI(EN) = 0.0D0
                  X = H(EN,NA)
                  S = DABS(X) + DABS(ZZ)
                  P = X/S
                  Q = ZZ/S
                  R = DSQRT(P*P+Q*Q)
                  P = P/R
                  Q = Q/R
*
*        INCREMENT OP COUNT FOR FINDING TWO ROOTS.
                  OPST = OPST + 18
*
*        INCREMENT OP COUNT FOR MODIFICATION AND ACCUMULATION
*        IN LOOP 290, 300, 310
                  OPS = OPS + 6* (N-NA+1) + 6*EN + 6* (IGH-LOW+1)
C     .......... ROW MODIFICATION ..........
                  DO 190 J = NA,N
                      ZZ = H(NA,J)
                      H(NA,J) = Q*ZZ + P*H(EN,J)
                      H(EN,J) = Q*H(EN,J) - P*ZZ
  190             CONTINUE
C     .......... COLUMN MODIFICATION ..........
                  DO 200 I = 1,EN
                      ZZ = H(I,NA)
                      H(I,NA) = Q*ZZ + P*H(I,EN)
                      H(I,EN) = Q*H(I,EN) - P*ZZ
  200             CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
                  DO 210 I = LOW,IGH
                      ZZ = Z(I,NA)
                      Z(I,NA) = Q*ZZ + P*Z(I,EN)
                      Z(I,EN) = Q*Z(I,EN) - P*ZZ
  210             CONTINUE
C
              END IF
              EN = ENM2
              GO TO 40
C     .......... ONE ROOT FOUND ..........
  220         H(EN,EN) = X + T
              WR(EN) = H(EN,EN)
              WI(EN) = 0.0D0
              EN = NA
              GO TO 40
          END IF
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
  230     IERR = EN
          GO TO 410
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  240     IF (NORM.NE.0.0D0) THEN
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
              DO 350 NN = 1,N
                  EN = N + 1 - NN
                  P = WR(EN)
                  Q = WI(EN)
                  NA = EN - 1
                  IF (Q) 250,300,350
C     .......... COMPLEX VECTOR ..........
  250             M = NA
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
                  IF (DABS(H(EN,NA)).LE.DABS(H(NA,EN))) THEN
                      CALL CDIV(0.0D0,-H(NA,EN),H(NA,NA)-P,Q,H(NA,NA),
     +                          H(NA,EN))
*
*        INCREMENT OP COUNT IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN)))
                      OPST = OPST + 16
                  ELSE
                      H(NA,NA) = Q/H(EN,NA)
                      H(NA,EN) = - (H(EN,EN)-P)/H(EN,NA)
*
*        INCREMENT OP COUNT.
                      OPST = OPST + 3
                  END IF
                  H(EN,NA) = 0.0D0
                  H(EN,EN) = 1.0D0
                  ENM2 = NA - 1
                  IF (ENM2.NE.0) THEN
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
                      DO 290 II = 1,ENM2
                          I = NA - II
                          W = H(I,I) - P
                          RA = 0.0D0
                          SA = 0.0D0
C
*
*        INCREMENT OP COUNT FOR LOOP 760
                          OPST = OPST + 4* (EN-M+1)
                          DO 260 J = M,EN
                              RA = RA + H(I,J)*H(J,NA)
                              SA = SA + H(I,J)*H(J,EN)
  260                     CONTINUE
C
                          IF (WI(I).GE.0.0D0) THEN
                              M = I
                              IF (WI(I).NE.0.0D0) THEN
C     .......... SOLVE COMPLEX EQUATIONS ..........
                                  X = H(I,I+1)
                                  Y = H(I+1,I)
                                  VR = (WR(I)-P)* (WR(I)-P) +
     +                                 WI(I)*WI(I) - Q*Q
                                  VI = (WR(I)-P)*2.0D0*Q
*
*        INCREMENT OPCOUNT (AVERAGE) FOR SOLVING COMPLEX EQUATIONS
                                  OPST = OPST + 42
                                  IF (VR.EQ.0.0D0 .AND.
     +                                VI.EQ.0.0D0) THEN
                                      TST1 = NORM* (DABS(W)+DABS(Q)+
     +                                       DABS(X)+DABS(Y)+DABS(ZZ))
                                      VR = TST1
  270                                 CONTINUE
                                      VR = 0.01D0*VR
                                      TST2 = TST1 + VR
                                      IF (TST2.GT.TST1) GO TO 270
                                  END IF
                                  CALL CDIV(X*R-ZZ*RA+Q*SA,
     +                                      X*S-ZZ*SA-Q*RA,VR,VI,
     +                                      H(I,NA),H(I,EN))
                                  IF (DABS(X).LE.DABS(ZZ)+DABS(Q)) THEN
                                      CALL CDIV(-R-Y*H(I,NA),
     +                                          -S-Y*H(I,EN),ZZ,Q,
     +                                          H(I+1,NA),H(I+1,EN))
                                  ELSE
                                      H(I+1,NA) = (-RA-W*H(I,NA)+
     +                                            Q*H(I,EN))/X
                                      H(I+1,EN) = (-SA-W*H(I,EN)-
     +                                            Q*H(I,NA))/X
                                  END IF
                              ELSE
                                  CALL CDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
*
*        INCREMENT OP COUNT FOR CDIV
                                  OPST = OPST + 16
                              END IF
C
C     .......... OVERFLOW CONTROL ..........
                              T = DMAX1(DABS(H(I,NA)),DABS(H(I,EN)))
                              IF (T.NE.0.0D0) THEN
                                  TST1 = T
                                  TST2 = TST1 + 1.0D0/TST1
                                  IF (TST2.LE.TST1) THEN
*
*        INCREMENT OP COUNT.
                                      OPST = OPST + 2* (EN-I+1)
                                      DO 280 J = I,EN
                                          H(J,NA) = H(J,NA)/T
                                          H(J,EN) = H(J,EN)/T
  280                                 CONTINUE
                                  END IF
                              END IF
                          ELSE
                              ZZ = W
                              R = RA
                              S = SA
                          END IF
  290                 CONTINUE
C
                  END IF
                  GO TO 350
C     .......... REAL VECTOR ..........
  300             M = EN
                  H(EN,EN) = 1.0D0
                  IF (NA.NE.0) THEN
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
                      DO 340 II = 1,NA
                          I = EN - II
                          W = H(I,I) - P
                          R = 0.0D0
C
*
*        INCREMENT OP COUNT FOR LOOP 610
                          OPST = OPST + 2* (EN-M+1)
                          DO 310 J = M,EN
                              R = R + H(I,J)*H(J,EN)
  310                     CONTINUE
C
                          IF (WI(I).GE.0.0D0) THEN
                              M = I
                              IF (WI(I).NE.0.0D0) THEN
C     .......... SOLVE REAL EQUATIONS ..........
                                  X = H(I,I+1)
                                  Y = H(I+1,I)
                                  Q = (WR(I)-P)* (WR(I)-P) + WI(I)*WI(I)
                                  T = (X*S-ZZ*R)/Q
*
*        INCREMENT OP COUNT FOR SOLVING REAL EQUATION.
                                  OPST = OPST + 13
                                  H(I,EN) = T
                                  IF (DABS(X).LE.DABS(ZZ)) THEN
                                      H(I+1,EN) = (-S-Y*T)/ZZ
                                  ELSE
                                      H(I+1,EN) = (-R-W*T)/X
                                  END IF
                              ELSE
                                  T = W
                                  IF (T.EQ.0.0D0) THEN
                                      TST1 = NORM
                                      T = TST1
  320                                 CONTINUE
                                      T = 0.01D0*T
                                      TST2 = NORM + T
                                      IF (TST2.GT.TST1) GO TO 320
                                  END IF
                                  H(I,EN) = -R/T
                              END IF
C
C     .......... OVERFLOW CONTROL ..........
                              T = DABS(H(I,EN))
                              IF (T.NE.0.0D0) THEN
                                  TST1 = T
                                  TST2 = TST1 + 1.0D0/TST1
                                  IF (TST2.LE.TST1) THEN
*
*        INCREMENT OP COUNT.
                                      OPST = OPST + (EN-I+1)
                                      DO 330 J = I,EN
                                          H(J,EN) = H(J,EN)/T
  330                                 CONTINUE
                                  END IF
                              END IF
                          ELSE
                              ZZ = W
                              S = R
                          END IF
  340                 CONTINUE
C
C     .......... END REAL VECTOR ..........
                  END IF
  350         CONTINUE
C     .......... END COMPLEX VECTOR ..........
C     .......... END BACK SUBSTITUTION.
C                VECTORS OF ISOLATED ROOTS ..........
              DO 370 I = 1,N
                  IF (I.LT.LOW .OR. I.GT.IGH) THEN
C
                      DO 360 J = I,N
                          Z(I,J) = H(I,J)
  360                 CONTINUE
                  END IF
  370         CONTINUE
C
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW DO -- ..........
              DO 400 JJ = LOW,N
                  J = N + LOW - JJ
                  M = MIN0(J,IGH)
C
*
*        INCREMENT OP COUNT.
                  OPS = OPS + 2* (IGH-LOW+1)* (M-LOW+1)
                  DO 390 I = LOW,IGH
                      ZZ = 0.0D0
C
                      DO 380 K = LOW,M
                          ZZ = ZZ + Z(I,K)*H(K,J)
  380                 CONTINUE
C
                      Z(I,J) = ZZ
  390             CONTINUE
  400         CONTINUE
C
          END IF
*
*     COMPUTE FINAL OP COUNT
  410     OPS = OPS + OPST
      END IF
      END
      SUBROUTINE IMTQL1(N,D,E,IERR)
*
*     EISPACK ROUTINE
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN DSTEQR.
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
      DOUBLE PRECISION ITCNT,OPS,OPST
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
      DOUBLE PRECISION D(N),E(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION B,C,EPS,F,G,P,R,S,TST
      INTEGER I,II,J,L,M,MML
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAMCH,PYTHAG
      EXTERNAL DLAMCH,PYTHAG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DSIGN,MIN
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
          EPS = DLAMCH('EPSILON')
C
          DO 10 I = 2,N
              E(I-1) = E(I)
   10     CONTINUE
C
          E(N) = 0.0D0
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
                      G = (D(L+1)-P)/ (2.0D0*E(L))
                      R = PYTHAG(G,1.0D0)
                      G = D(M) - P + E(L)/ (G+DSIGN(R,G))
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT.
                      OPS = OPS + 7
                      S = 1.0D0
                      C = 1.0D0
                      P = 0.0D0
                      MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
                      DO 50 II = 1,MML
                          I = M - II
                          F = S*E(I)
                          B = C*E(I)
                          R = PYTHAG(F,G)
                          E(I+1) = R
                          IF (R.EQ.0.0D0) THEN
                              GO TO 60
                          ELSE
                              S = F/R
                              C = G/R
                              G = D(I+1) - P
                              R = (D(I)-G)*S + 2.0D0*C*B
                              P = S*R
                              D(I+1) = G + P
                              G = C*R - B
                          END IF
   50                 CONTINUE
C
                      D(L) = D(L) - P
                      E(L) = G
                      E(M) = 0.0D0
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
                      OPS = OPS + MML*14 + 1
*
*        INCREMENT ITERATION COUNTER
                      ITCNT = ITCNT + 1
                      GO TO 20
C     .......... RECOVER FROM UNDERFLOW ..........
   60                 D(I+1) = D(I+1) - P
                      E(M) = 0.0D0
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
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN DSTEQR.
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
      DOUBLE PRECISION ITCNT,OPS,OPST
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
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION B,C,EPS,F,G,P,R,S,TST
      INTEGER I,II,J,K,L,M,MML
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAMCH,PYTHAG
      EXTERNAL DLAMCH,PYTHAG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DSIGN,MIN
C     ..
      IERR = 0
      IF (N.NE.1) THEN
*
*        INITIALIZE ITERATION COUNT AND OPST
          ITCNT = 0
          OPST = 0
*
*     DETERMINE UNIT ROUNDOFF FOR THIS MACHINE.
          EPS = DLAMCH('EPSILON')
C
          DO 10 I = 2,N
              E(I-1) = E(I)
   10     CONTINUE
C
          E(N) = 0.0D0
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
                      G = (D(L+1)-P)/ (2.0D0*E(L))
                      R = PYTHAG(G,1.0D0)
                      G = D(M) - P + E(L)/ (G+DSIGN(R,G))
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT.
                      OPS = OPS + 7
                      S = 1.0D0
                      C = 1.0D0
                      P = 0.0D0
                      MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
                      DO 60 II = 1,MML
                          I = M - II
                          F = S*E(I)
                          B = C*E(I)
                          R = PYTHAG(F,G)
                          E(I+1) = R
                          IF (R.EQ.0.0D0) THEN
                              GO TO 70
                          ELSE
                              S = F/R
                              C = G/R
                              G = D(I+1) - P
                              R = (D(I)-G)*S + 2.0D0*C*B
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
                      E(M) = 0.0D0
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
                      OPS = OPS + MML* (14+6*N) + 1
*
*        INCREMENT ITERATION COUNTER
                      ITCNT = ITCNT + 1
                      GO TO 20
C     .......... RECOVER FROM UNDERFLOW ..........
   70                 D(I+1) = D(I+1) - P
                      E(M) = 0.0D0
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
      SUBROUTINE INVIT(NM,N,A,WR,WI,SELECT,MM,M,Z,IERR,RM1,RV1,RV2)
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
      DOUBLE PRECISION ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE INVIT
C     BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A REAL UPPER
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
C        A CONTAINS THE HESSENBERG MATRIX.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVALUES OF THE MATRIX.  THE EIGENVALUES MUST BE
C          STORED IN A MANNER IDENTICAL TO THAT OF SUBROUTINE  HQR,
C          WHICH RECOGNIZES POSSIBLE SPLITTING OF THE MATRIX.
C
C        SELECT SPECIFIES THE EIGENVECTORS TO BE FOUND. THE
C          EIGENVECTOR CORRESPONDING TO THE J-TH EIGENVALUE IS
C          SPECIFIED BY SETTING SELECT(J) TO .TRUE..
C
C        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
C          COLUMNS REQUIRED TO STORE THE EIGENVECTORS TO BE FOUND.
C          NOTE THAT TWO COLUMNS ARE REQUIRED TO STORE THE
C          EIGENVECTOR CORRESPONDING TO A COMPLEX EIGENVALUE.
C
C     ON OUTPUT
C
C        A AND WI ARE UNALTERED.
C
C        WR MAY HAVE BEEN ALTERED SINCE CLOSE EIGENVALUES ARE PERTURBED
C          SLIGHTLY IN SEARCHING FOR INDEPENDENT EIGENVECTORS.
C
C        SELECT MAY HAVE BEEN ALTERED.  IF THE ELEMENTS CORRESPONDING
C          TO A PAIR OF CONJUGATE COMPLEX EIGENVALUES WERE EACH
C          INITIALLY SET TO .TRUE., THE PROGRAM RESETS THE SECOND OF
C          THE TWO ELEMENTS TO .FALSE..
C
C        M IS THE NUMBER OF COLUMNS ACTUALLY USED TO STORE
C          THE EIGENVECTORS.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE NEXT SELECTED EIGENVALUE IS REAL, THE NEXT COLUMN
C          OF Z CONTAINS ITS EIGENVECTOR.  IF THE EIGENVALUE IS
C          COMPLEX, THE NEXT TWO COLUMNS OF Z CONTAIN THE REAL AND
C          IMAGINARY PARTS OF ITS EIGENVECTOR.  THE EIGENVECTORS ARE
C          NORMALIZED SO THAT THE COMPONENT OF LARGEST MAGNITUDE IS 1.
C          ANY VECTOR WHICH FAILS THE ACCEPTANCE TEST IS SET TO ZERO.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -(2*N+1)   IF MORE THAN MM COLUMNS OF Z ARE NECESSARY
C                     TO STORE THE EIGENVECTORS CORRESPONDING TO
C                     THE SPECIFIED EIGENVALUES.
C          -K         IF THE ITERATION CORRESPONDING TO THE K-TH
C                     VALUE FAILS,
C          -(N+K)     IF BOTH ERROR SITUATIONS OCCUR.
C
C        RM1, RV1, AND RV2 ARE TEMPORARY STORAGE ARRAYS.  NOTE THAT RM1
C          IS SQUARE OF DIMENSION N BY N AND, AUGMENTED BY TWO COLUMNS
C          OF Z, IS THE TRANSPOSE OF THE CORRESPONDING ALGOL B ARRAY.
C
C     THE ALGOL PROCEDURE GUESSVEC APPEARS IN INVIT IN LINE.
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
*
*     GET ULP FROM DLAMCH FOR NEW SMALL PERTURBATION AS IN LAPACK

C     .. Scalar Arguments ..
      INTEGER IERR,M,MM,N,NM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NM,N),RM1(N,N),RV1(N),RV2(N),WI(N),WR(N),
     +                 Z(NM,MM)
      LOGICAL SELECT(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION EPS3,GROWTO,ILAMBD,NORM,NORMV,OPST,RLAMBD,T,
     +                 UKROOT,ULP,W,X,Y
      INTEGER I,II,IP,IP1,ITS,J,K,KM1,L,MP,N1,NS,S,UK
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAMCH,PYTHAG
      EXTERNAL DLAMCH,PYTHAG
C     ..
C     .. External Subroutines ..
      EXTERNAL CDIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DSQRT,IABS
C     ..
      IF (N.LE.0) RETURN
      ULP = DLAMCH('EPSILON')
C
*
*     INITIALIZE
      OPST = 0
      IERR = 0
      UK = 0
      S = 1
C     .......... IP = 0, REAL EIGENVALUE
C                     1, FIRST OF CONJUGATE COMPLEX PAIR
C                    -1, SECOND OF CONJUGATE COMPLEX PAIR ..........
      IP = 0
      N1 = N - 1
C
      DO 610 K = 1,N
          IF (WI(K).EQ.0.0D0 .OR. IP.LT.0) GO TO 10
          IP = 1
          IF (SELECT(K) .AND. SELECT(K+1)) SELECT(K+1) = .FALSE.
   10     IF (.NOT.SELECT(K)) GO TO 600
          IF (WI(K).NE.0.0D0) S = S + 1
          IF (S.GT.MM) GO TO 620
          IF (UK.GE.K) GO TO 60
C     .......... CHECK FOR POSSIBLE SPLITTING ..........
          DO 20 UK = K,N
              IF (UK.EQ.N) GO TO 30
              IF (A(UK+1,UK).EQ.0.0D0) GO TO 30
   20     CONTINUE
C     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK
C                (HESSENBERG) MATRIX ..........
   30     NORM = 0.0D0
          MP = 1
C
*
*        INCREMENT OPCOUNT FOR COMPUTING MATRIX NORM
          OPS = OPS + UK* (UK-1)/2
          DO 50 I = 1,UK
              X = 0.0D0
C
              DO 40 J = MP,UK
                  X = X + DABS(A(I,J))
   40         CONTINUE
C
              IF (X.GT.NORM) NORM = X
              MP = I
   50     CONTINUE
C     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
C                AND CLOSE ROOTS ARE MODIFIED BY EPS3 ..........
          IF (NORM.EQ.0.0D0) NORM = 1.0D0
*        EPS3 = EPSLON(NORM)
*
*        INCREMENT OPCOUNT
          OPST = OPST + 3
          EPS3 = NORM*ULP
C     .......... GROWTO IS THE CRITERION FOR THE GROWTH ..........
          UKROOT = UK
          UKROOT = DSQRT(UKROOT)
          GROWTO = 0.1D0/UKROOT
   60     RLAMBD = WR(K)
          ILAMBD = WI(K)
          IF (K.EQ.1) GO TO 100
          KM1 = K - 1
          GO TO 80
C     .......... PERTURB EIGENVALUE IF IT IS CLOSE
C                TO ANY PREVIOUS EIGENVALUE ..........
   70     RLAMBD = RLAMBD + EPS3
C*PL*ERROR* Embedded comment after label moved
C     .......... FOR I=K-1 STEP -1 UNTIL 1 DO -- ..........
   80     DO 90 II = 1,KM1
              I = K - II
              IF (SELECT(I) .AND. DABS(WR(I)-RLAMBD).LT.EPS3 .AND.
     +            DABS(WI(I)-ILAMBD).LT.EPS3) GO TO 70
   90     CONTINUE
*
*        INCREMENT OPCOUNT FOR LOOP 260 (ASSUME THAT ALL EIGENVALUES
*        ARE DIFFERENT)
          OPST = OPST + 2* (K-1)
C
          WR(K) = RLAMBD
C     .......... PERTURB CONJUGATE EIGENVALUE TO MATCH ..........
          IP1 = K + IP
          WR(IP1) = RLAMBD
C     .......... FORM UPPER HESSENBERG A-RLAMBD*I (TRANSPOSED)
C                AND INITIAL REAL VECTOR ..........
  100     MP = 1
C
*
*        INCREMENT OP COUNT FOR LOOP 320
          OPS = OPS + UK
          DO 120 I = 1,UK
C
              DO 110 J = MP,UK
                  RM1(J,I) = A(I,J)
  110         CONTINUE
C
              RM1(I,I) = RM1(I,I) - RLAMBD
              MP = I
              RV1(I) = EPS3
  120     CONTINUE
C
          ITS = 0
          IF (ILAMBD.NE.0.0D0) GO TO 230
C     .......... REAL EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3 ..........
          IF (UK.EQ.1) GO TO 180
C
*
*        INCREMENT OPCOUNT LU DECOMPOSITION
          OPS = OPS + (UK-1)* (UK+2)
          DO 170 I = 2,UK
              MP = I - 1
              IF (DABS(RM1(MP,I)).LE.DABS(RM1(MP,MP))) GO TO 140
C
              DO 130 J = MP,UK
                  Y = RM1(J,I)
                  RM1(J,I) = RM1(J,MP)
                  RM1(J,MP) = Y
  130         CONTINUE
C
  140         IF (RM1(MP,MP).EQ.0.0D0) RM1(MP,MP) = EPS3
              X = RM1(MP,I)/RM1(MP,MP)
              IF (X.EQ.0.0D0) GO TO 160
C
              DO 150 J = I,UK
                  RM1(J,I) = RM1(J,I) - X*RM1(J,MP)
  150         CONTINUE
  160         CONTINUE
  170     CONTINUE
C
C
  180     IF (RM1(UK,UK).EQ.0.0D0) RM1(UK,UK) = EPS3
C*PL*ERROR* Embedded comment after label moved
C     .......... BACK SUBSTITUTION FOR REAL VECTOR
C*PL*ERROR* Embedded comment after label moved
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  190     DO 220 II = 1,UK
              I = UK + 1 - II
              Y = RV1(I)
              IF (I.EQ.UK) GO TO 210
              IP1 = I + 1
C
              DO 200 J = IP1,UK
                  Y = Y - RM1(J,I)*RV1(J)
  200         CONTINUE
C
  210         RV1(I) = Y/RM1(I,I)
  220     CONTINUE
*
*        INCREMENT OP COUNT FOR BACK SUBSTITUTION LOOP 500
          OPS = OPS + UK* (UK+1)
C
          GO TO 480
C     .......... COMPLEX EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3.  STORE IMAGINARY
C                PARTS IN UPPER TRIANGLE STARTING AT (1,3) ..........
  230     NS = N - S
          Z(1,S-1) = -ILAMBD
          Z(1,S) = 0.0D0
          IF (N.EQ.2) GO TO 250
          RM1(1,3) = -ILAMBD
          Z(1,S-1) = 0.0D0
          IF (N.EQ.3) GO TO 250
C
          DO 240 I = 4,N
              RM1(1,I) = 0.0D0
  240     CONTINUE
C*PL*ERROR* Embedded comment after label moved
C
  250     DO 370 I = 2,UK
              MP = I - 1
              W = RM1(MP,I)
              IF (I.LT.N) T = RM1(MP,I+1)
              IF (I.EQ.N) T = Z(MP,S-1)
              X = RM1(MP,MP)*RM1(MP,MP) + T*T
              IF (W*W.LE.X) GO TO 300
              X = RM1(MP,MP)/W
              Y = T/W
              RM1(MP,MP) = W
              IF (I.LT.N) RM1(MP,I+1) = 0.0D0
              IF (I.EQ.N) Z(MP,S-1) = 0.0D0
C
*
*        INCREMENT OPCOUNT FOR LOOP 560
              OPS = OPS + 4* (UK-I+1)
              DO 280 J = I,UK
                  W = RM1(J,I)
                  RM1(J,I) = RM1(J,MP) - X*W
                  RM1(J,MP) = W
                  IF (J.LT.N1) GO TO 260
                  L = J - NS
                  Z(I,L) = Z(MP,L) - Y*W
                  Z(MP,L) = 0.0D0
                  GO TO 270
  260             RM1(I,J+2) = RM1(MP,J+2) - Y*W
                  RM1(MP,J+2) = 0.0D0
  270             CONTINUE
  280         CONTINUE
C
              RM1(I,I) = RM1(I,I) - Y*ILAMBD
              IF (I.LT.N1) GO TO 290
              L = I - NS
              Z(MP,L) = -ILAMBD
              Z(I,L) = Z(I,L) + X*ILAMBD
              GO TO 360
  290         RM1(MP,I+2) = -ILAMBD
              RM1(I,I+2) = RM1(I,I+2) + X*ILAMBD
              GO TO 360
  300         IF (X.NE.0.0D0) GO TO 310
              RM1(MP,MP) = EPS3
              IF (I.LT.N) RM1(MP,I+1) = 0.0D0
              IF (I.EQ.N) Z(MP,S-1) = 0.0D0
              T = 0.0D0
              X = EPS3*EPS3
  310         W = W/X
              X = RM1(MP,MP)*W
              Y = -T*W
C
*
*        INCREMENT OPCOUNT FOR LOOP 620
              OPS = OPS + 6* (UK-I+1)
              DO 340 J = I,UK
                  IF (J.LT.N1) GO TO 320
                  L = J - NS
                  T = Z(MP,L)
                  Z(I,L) = -X*T - Y*RM1(J,MP)
                  GO TO 330
  320             T = RM1(MP,J+2)
                  RM1(I,J+2) = -X*T - Y*RM1(J,MP)
  330             RM1(J,I) = RM1(J,I) - X*RM1(J,MP) + Y*T
  340         CONTINUE
C
              IF (I.LT.N1) GO TO 350
              L = I - NS
              Z(I,L) = Z(I,L) - ILAMBD
              GO TO 360
  350         RM1(I,I+2) = RM1(I,I+2) - ILAMBD
  360         CONTINUE
  370     CONTINUE
*
*        INCREMENT OP COUNT (AVERAGE) FOR COMPUTING
*        THE SCALARS IN LOOP 640
          OPS = OPS + 10* (UK-1)
C
          IF (UK.LT.N1) GO TO 380
          L = UK - NS
          T = Z(UK,L)
          GO TO 390
  380     T = RM1(UK,UK+2)
  390     IF (RM1(UK,UK).EQ.0.0D0 .AND. T.EQ.0.0D0) RM1(UK,UK) = EPS3
C*PL*ERROR* Embedded comment after label moved
C     .......... BACK SUBSTITUTION FOR COMPLEX VECTOR
C*PL*ERROR* Embedded comment after label moved
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  400     DO 470 II = 1,UK
              I = UK + 1 - II
              X = RV1(I)
              Y = 0.0D0
              IF (I.EQ.UK) GO TO 440
              IP1 = I + 1
C
              DO 430 J = IP1,UK
                  IF (J.LT.N1) GO TO 410
                  L = J - NS
                  T = Z(I,L)
                  GO TO 420
  410             T = RM1(I,J+2)
  420             X = X - RM1(J,I)*RV1(J) + T*RV2(J)
                  Y = Y - RM1(J,I)*RV2(J) - T*RV1(J)
  430         CONTINUE
C
  440         IF (I.LT.N1) GO TO 450
              L = I - NS
              T = Z(I,L)
              GO TO 460
  450         T = RM1(I,I+2)
  460         CALL CDIV(X,Y,RM1(I,I),T,RV1(I),RV2(I))
  470     CONTINUE
*
*        INCREMENT OP COUNT FOR LOOP 720.
          OPS = OPS + 4*UK* (UK+3)
C     .......... ACCEPTANCE TEST FOR REAL OR COMPLEX
C                EIGENVECTOR AND NORMALIZATION ..........
  480     ITS = ITS + 1
          NORM = 0.0D0
          NORMV = 0.0D0
C
          DO 500 I = 1,UK
              IF (ILAMBD.EQ.0.0D0) X = DABS(RV1(I))
              IF (ILAMBD.NE.0.0D0) X = PYTHAG(RV1(I),RV2(I))
              IF (NORMV.GE.X) GO TO 490
              NORMV = X
              J = I
  490         NORM = NORM + X
  500     CONTINUE
*
*        INCREMENT OP COUNT ACCEPTANCE TEST
          IF (ILAMBD.EQ.0.0D0) OPS = OPS + UK
          IF (ILAMBD.NE.0.0D0) OPS = OPS + 16*UK
C
          IF (NORM.LT.GROWTO) GO TO 540
C     .......... ACCEPT VECTOR ..........
          X = RV1(J)
          IF (ILAMBD.EQ.0.0D0) X = 1.0D0/X
          IF (ILAMBD.NE.0.0D0) Y = RV2(J)
C
*
*        INCREMENT OPCOUNT FOR LOOP 820
          IF (ILAMBD.EQ.0.0D0) OPS = OPS + UK
          IF (ILAMBD.NE.0.0D0) OPS = OPS + 16*UK
          DO 530 I = 1,UK
              IF (ILAMBD.NE.0.0D0) GO TO 510
              Z(I,S) = RV1(I)*X
              GO TO 520
  510         CALL CDIV(RV1(I),RV2(I),X,Y,Z(I,S-1),Z(I,S))
  520         CONTINUE
  530     CONTINUE
C
          IF (UK.EQ.N) GO TO 590
          J = UK + 1
          GO TO 570
C     .......... IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR ..........
  540     IF (ITS.GE.UK) GO TO 560
          X = UKROOT
          Y = EPS3/ (X+1.0D0)
          RV1(1) = EPS3
C
          DO 550 I = 2,UK
              RV1(I) = Y
  550     CONTINUE
C
          J = UK - ITS + 1
          RV1(J) = RV1(J) - EPS3*X
          IF (ILAMBD.EQ.0.0D0) GO TO 190
          GO TO 400
C     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  560     J = 1
          IERR = -K
C*PL*ERROR* Embedded comment after label moved
C     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  570     DO 580 I = J,N
              Z(I,S) = 0.0D0
              IF (ILAMBD.NE.0.0D0) Z(I,S-1) = 0.0D0
  580     CONTINUE
C
  590     S = S + 1
  600     IF (IP.EQ. (-1)) IP = 0
          IF (IP.EQ.1) IP = -1
  610 CONTINUE
C
      GO TO 630
C     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
C                SPACE REQUIRED ..........
  620 IF (IERR.NE.0) IERR = IERR - N
      IF (IERR.EQ.0) IERR = - (2*N+1)
  630 M = S - 1 - IABS(IP)
*
*     COMPUTE FINAL OP COUNT
      OPS = OPS + OPST
      RETURN
      END
      SUBROUTINE ORTHES(NM,N,LOW,IGH,A,ORT)
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
      DOUBLE PRECISION ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
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
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT
C
C        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLE UNDER THE
C          HESSENBERG MATRIX.
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
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
      DOUBLE PRECISION A(NM,N),ORT(IGH)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,G,H,SCALE
      INTEGER I,II,J,JJ,KP1,LA,M,MP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DSIGN,DSQRT
C     ..
      IF (N.GT.0) THEN
          LA = IGH - 1
          KP1 = LOW + 1
          IF (LA.GE.KP1) THEN
C
*
*     INCREMENT OP COUNR FOR COMPUTING G,H,ORT(M),.. IN LOOP 180
              OPS = OPS + 6* (LA-KP1+1)
              DO 90 M = KP1,LA
                  H = 0.0D0
                  ORT(M) = 0.0D0
                  SCALE = 0.0D0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
*
*     INCREMENT OP COUNT FOR LOOP 90
                  OPS = OPS + (IGH-M+1)
                  DO 10 I = M,IGH
                      SCALE = SCALE + DABS(A(I,M-1))
   10             CONTINUE
C
                  IF (SCALE.NE.0.0D0) THEN
                      MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
*
*     INCREMENT OP COUNT FOR LOOP 100
                      OPS = OPS + 3* (IGH-M+1)
                      DO 20 II = M,IGH
                          I = MP - II
                          ORT(I) = A(I,M-1)/SCALE
                          H = H + ORT(I)*ORT(I)
   20                 CONTINUE
C
                      G = -DSIGN(DSQRT(H),ORT(M))
                      H = H - ORT(M)*G
                      ORT(M) = ORT(M) - G
C     .......... FORM (I-(U*UT)/H) * A ..........
*
*     INCREMENT OP COUNT FOR LOOP 130 AND 160
                      OPS = OPS + (N-M+1+IGH)* (4* (IGH-M+1)+1)
                      DO 50 J = M,N
                          F = 0.0D0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
                          DO 30 II = M,IGH
                              I = MP - II
                              F = F + ORT(I)*A(I,J)
   30                     CONTINUE
C
                          F = F/H
C
                          DO 40 I = M,IGH
                              A(I,J) = A(I,J) - F*ORT(I)
   40                     CONTINUE
   50                 CONTINUE
C
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
                      DO 80 I = 1,IGH
                          F = 0.0D0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
                          DO 60 JJ = M,IGH
                              J = MP - JJ
                              F = F + ORT(J)*A(I,J)
   60                     CONTINUE
C
                          F = F/H
C
                          DO 70 J = M,IGH
                              A(I,J) = A(I,J) - F*ORT(J)
   70                     CONTINUE
   80                 CONTINUE
C
C
                      ORT(M) = SCALE*ORT(M)
                      A(M,M-1) = SCALE*G
                  END IF
   90         CONTINUE
          END IF
          RETURN
      END IF
*$st$ Unreachable comments ...
C
      END
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
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
      DOUBLE PRECISION A,B
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION OPST
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION P,R,S,T,U
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DMAX1,DMIN1
C     ..
      P = DMAX1(DABS(A),DABS(B))
      IF (P.NE.0.0D0) THEN
          R = (DMIN1(DABS(A),DABS(B))/P)**2
*
*     INCREMENT OPST
          OPST = OPST + 2
   10     CONTINUE
          T = 4.0D0 + R
          IF (T.NE.4.0D0) THEN
              S = R/T
              U = 1.0D0 + 2.0D0*S
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
      SUBROUTINE TQLRAT(N,D,E2,IERR)
*
*     EISPACK ROUTINE.
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN DSTEQR.
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
      DOUBLE PRECISION ITCNT,OPS,OPST
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E2 HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
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
      DOUBLE PRECISION D(N),E2(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION B,C,EPS,F,G,H,P,R,S,T,TST
      INTEGER I,II,J,L,L1,M,MML
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAMCH,EPSLON,PYTHAG
      EXTERNAL DLAMCH,EPSLON,PYTHAG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DABS,DSIGN,DSQRT,MIN,SQRT
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
          EPS = DLAMCH('EPSILON')
C
          DO 10 I = 2,N
              E2(I-1) = E2(I)
   10     CONTINUE
C
          F = 0.0D0
          T = 0.0D0
          E2(N) = 0.0D0
C
          DO 90 L = 1,N
              J = 0
              H = DABS(D(L)) + DSQRT(E2(L))
              IF (T.LE.H) THEN
                  T = H
                  B = EPSLON(T)
                  C = B*B
*
*     INCREMENT OPCOUNT FOR THIS SECTION.
*     (FUNCTION EPSLON IS COUNTED AS 6 FLOPS.  THIS IS THE MINIMUM
*     NUMBER REQUIRED, BUT COUNTING THEM EXACTLY WOULD AFFECT
*     THE TIMING.)
                  OPS = OPS + 9
              END IF
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
              DO 20 M = L,N
                  IF (M.EQ.N) THEN
                      GO TO 30
                  ELSE
                      TST = SQRT(ABS(E2(M)))
                      IF (TST.LE.EPS* (ABS(D(M))+ABS(D(M+1))))
     +                    GO TO 30
                  END IF
   20         CONTINUE
*            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
C
   30         CONTINUE
*
*        INCREMENT OPCOUNT FOR FINDING SMALL SUBDIAGONAL ELEMENT.
              OPS = OPS + 3* (MIN(M,N-1)-L+1)
              IF (M.NE.L) THEN
   40             CONTINUE
                  IF (J.EQ.30) THEN
                      GO TO 100
                  ELSE
                      J = J + 1
C     .......... FORM SHIFT ..........
                      L1 = L + 1
                      S = DSQRT(E2(L))
                      G = D(L)
                      P = (D(L1)-G)/ (2.0D0*S)
                      R = PYTHAG(P,1.0D0)
                      D(L) = S/ (P+DSIGN(R,P))
                      H = G - D(L)
C
                      DO 50 I = L1,N
                          D(I) = D(I) - H
   50                 CONTINUE
C
                      F = F + H
*
*        INCREMENT OPCOUNT FOR FORMING SHIFT AND SUBTRACTING.
                      OPS = OPS + 8 + (I-L1+1)
C     .......... RATIONAL QL TRANSFORMATION ..........
                      G = D(M)
                      IF (G.EQ.0.0D0) G = B
                      H = G
                      S = 0.0D0
                      MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
                      DO 60 II = 1,MML
                          I = M - II
                          P = G*H
                          R = P + E2(I)
                          E2(I+1) = S*R
                          S = E2(I)/R
                          D(I+1) = H + S* (H+D(I))
                          G = D(I) - E2(I)/G
                          IF (G.EQ.0.0D0) G = B
                          H = G*P/R
   60                 CONTINUE
C
                      E2(L) = S*G
                      D(L) = H
*
*        INCREMENT OPCOUNT FOR INNER LOOP.
                      OPS = OPS + MML*11 + 1
*
*        INCREMENT ITERATION COUNTER
                      ITCNT = ITCNT + 1
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
                      IF (H.NE.0.0D0) THEN
                          IF (DABS(E2(L)).GT.DABS(C/H)) THEN
                              E2(L) = H*E2(L)
                              IF (E2(L).NE.0.0D0) GO TO 40
                          END IF
                      END IF
                  END IF
              END IF
              P = D(L) + F
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
C                EIGENVALUE AFTER 30 ITERATIONS ..........
  100     IERR = L
      END IF
*
*     COMPUTE FINAL OP COUNT
  110 OPS = OPS + OPST
      END
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT.
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED.
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ITCNT,OPS
C     ..
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
C     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
C          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
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
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,G,H,SCALE
      INTEGER I,II,J,JP1,K,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DBLE,DSIGN,DSQRT,MAX
C     ..
      OPS = OPS + MAX(0.0D0, (4.0D0/3.0D0)*DBLE(N)**3+12.0D0*DBLE(N)**2+
     +      (11.0D0/3.0D0)*N-22)
*
      DO 10 I = 1,N
          D(I) = A(N,I)
          A(N,I) = A(I,I)
   10 CONTINUE
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 130 II = 1,N
          I = N + 1 - II
          L = I - 1
          H = 0.0D0
          SCALE = 0.0D0
          IF (L.GE.1) THEN
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
              DO 20 K = 1,L
                  SCALE = SCALE + DABS(D(K))
   20         CONTINUE
C
              IF (SCALE.NE.0.0D0) THEN
C
                  DO 30 K = 1,L
                      D(K) = D(K)/SCALE
                      H = H + D(K)*D(K)
   30             CONTINUE
C
                  E2(I) = SCALE*SCALE*H
                  F = D(L)
                  G = -DSIGN(DSQRT(H),F)
                  E(I) = SCALE*G
                  H = H - F*G
                  D(L) = F - G
                  IF (L.NE.1) THEN
C     .......... FORM A*U ..........
                      DO 40 J = 1,L
                          E(J) = 0.0D0
   40                 CONTINUE
C
                      DO 60 J = 1,L
                          F = D(J)
                          G = E(J) + A(J,J)*F
                          JP1 = J + 1
                          IF (L.GE.JP1) THEN
C
                              DO 50 K = JP1,L
                                  G = G + A(K,J)*D(K)
                                  E(K) = E(K) + A(K,J)*F
   50                         CONTINUE
                          END IF
C
                          E(J) = G
   60                 CONTINUE
C     .......... FORM P ..........
                      F = 0.0D0
C
                      DO 70 J = 1,L
                          E(J) = E(J)/H
                          F = F + E(J)*D(J)
   70                 CONTINUE
C
                      H = F/ (H+H)
C     .......... FORM Q ..........
                      DO 80 J = 1,L
                          E(J) = E(J) - H*D(J)
   80                 CONTINUE
C     .......... FORM REDUCED A ..........
                      DO 100 J = 1,L
                          F = D(J)
                          G = E(J)
C
                          DO 90 K = J,L
                              A(K,J) = A(K,J) - F*E(K) - G*D(K)
   90                     CONTINUE
  100                 CONTINUE
C
                  END IF
C
                  DO 110 J = 1,L
                      F = D(J)
                      D(J) = A(L,J)
                      A(L,J) = A(I,J)
                      A(I,J) = F*SCALE
  110             CONTINUE
                  GO TO 130
              ELSE
C
                  DO 120 J = 1,L
                      D(J) = A(L,J)
                      A(L,J) = A(I,J)
                      A(I,J) = 0.0D0
  120             CONTINUE
              END IF
          END IF
C
          E(I) = 0.0D0
          E2(I) = 0.0D0
  130 CONTINUE
C
C
      END
      SUBROUTINE BISECT(N,EPS1,D,E,E2,LB,UB,MM,M,W,IND,IERR,RV4,RV5)
*
*     EISPACK ROUTINE.
*     MODIFIED FOR COMPARISON WITH LAPACK ROUTINES.
*
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN DSTEBZ.
*
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ITCNT,OPS
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
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION RELFAC
      PARAMETER (RELFAC=2.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS1,LB,UB
      INTEGER IERR,M,MM,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION D(N),E(N),E2(N),RV4(N),RV5(N),W(MM)
      INTEGER IND(MM)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ATOLI,PIVMIN,RTOLI,SAFEMN,T1,T2,TMP1,TMP2,TNORM,
     +                 U,ULP,V,X0,X1,XU
      INTEGER I,II,ISTURM,J,K,L,M1,M2,P,Q,R,S,TAG
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAMCH,EPSLON
      EXTERNAL DLAMCH,EPSLON
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DABS,DMAX1,DMIN1,MAX,MIN
C     ..
      ITCNT = 0
      SAFEMN = DLAMCH('S')
      ULP = DLAMCH('E')*DLAMCH('B')
      RTOLI = ULP*RELFAC
      IERR = 0
      TAG = 0
      T1 = LB
      T2 = UB
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
      DO 30 I = 1,N
          IF (I.EQ.1) GO TO 10
CCC         TST1 = DABS(D(I)) + DABS(D(I-1))
CCC         TST2 = TST1 + DABS(E(I))
CCC         IF (TST2 .GT. TST1) GO TO 40
          TMP1 = E(I)**2
          IF (ABS(D(I)*D(I-1))*ULP**2+SAFEMN.LE.TMP1) GO TO 20
   10     E2(I) = 0.0D0
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
      U = 0.0D0
C
      DO 90 Q = P,N
          X1 = U
          U = 0.0D0
          V = 0.0D0
          IF (Q.EQ.N) GO TO 80
          U = DABS(E(Q+1))
          V = E2(Q+1)
   80     XU = DMIN1(D(Q)- (X1+U),XU)
          X0 = DMAX1(D(Q)+ (X1+U),X0)
          IF (V.EQ.0.0D0) GO TO 100
   90 CONTINUE
*        INCREMENT OPCOUNT FOR REFINING INTERVAL.
      OPS = OPS + (N-P+1)*2
C
  100 X1 = EPSLON(DMAX1(DABS(XU),DABS(X0)))
      IF (EPS1.LE.0.0D0) EPS1 = -X1
      IF (P.NE.Q) GO TO 110
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1.GT.D(P) .OR. D(P).GE.T2) GO TO 340
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 280
  110 X1 = X1* (Q-P+1)
      LB = DMAX1(T1,XU-X1)
      UB = DMIN1(T2,X0+X1)
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
  190 X1 = (XU+X0)*0.5D0
CCC         IF ((X0 - XU) .LE. DABS(EPS1)) GO TO 420
CCC         TST1 = 2.0D0 * (DABS(XU) + DABS(X0))
CCC         TST2 = TST1 + (X0 - XU)
CCC         IF (TST2 .EQ. TST1) GO TO 420
      TMP1 = ABS(X0-XU)
      TMP2 = MAX(ABS(X0),ABS(XU))
      IF (TMP1.LT.MAX(ATOLI,PIVMIN,RTOLI*TMP2)) GO TO 270
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  200 S = P - 1
      U = 1.0D0
C
      DO 230 I = P,Q
          IF (U.NE.0.0D0) GO TO 210
          V = DABS(E(I))/EPSLON(1.0D0)
          IF (E2(I).EQ.0.0D0) V = 0.0D0
          GO TO 220
  210     V = E2(I)/U
  220     U = D(I) - X1 - V
          IF (U.LT.0.0D0) S = S + 1
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
      DOUBLE PRECISION ITCNT,OPS,OPST
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
C          0.0D0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0D0
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
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
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
      DOUBLE PRECISION D(N),E(N),E2(N),RV1(N),RV2(N),RV3(N),RV4(N),
     +                 RV6(N),W(M),Z(NM,M)
      INTEGER IND(M)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION EPS2,EPS3,EPS4,NORM,ORDER,U,UK,V,X0,X1,XU
      INTEGER GROUP,I,II,IP,ITS,J,JJ,P,Q,R,S,TAG
C     ..
C     .. External Functions ..
      DOUBLE PRECISION EPSLON,PYTHAG
      EXTERNAL EPSLON,PYTHAG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DMAX1,DSQRT
C     ..
      ITCNT = 0
      IERR = 0
      IF (M.NE.0) THEN
          TAG = 0
          ORDER = 1.0D0 - E2(1)
          Q = 0
   10     CONTINUE
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
          P = Q + 1
C
          DO 20 Q = P,N
              IF (Q.EQ.N) THEN
                  GO TO 30
              ELSE IF (E2(Q+1).EQ.0.0D0) THEN
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
                      XU = 1.0D0
                      IF (P.NE.Q) THEN
                          NORM = DABS(D(P))
                          IP = P + 1
C
                          DO 40 I = IP,Q
                              NORM = DMAX1(NORM,DABS(D(I))+DABS(E(I)))
   40                     CONTINUE
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
                          EPS2 = 1.0D-3*NORM
                          EPS3 = EPSLON(NORM)
                          UK = Q - P + 1
                          EPS4 = UK*EPS3
                          UK = EPS4/DSQRT(UK)
*           INCREMENT OPCOUNT FOR COMPUTING CRITERIA.
                          OPS = OPS + (Q-IP+4)
                          S = P
                      ELSE
                          RV6(P) = 1.0D0
                          GO TO 180
                      END IF
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
                  ELSE IF (DABS(X1-X0).LT.EPS2) THEN
                      GROUP = GROUP + 1
                      IF (ORDER* (X1-X0).LE.0.0D0) X1 = X0 + ORDER*EPS3
                      GO TO 50
                  END IF
                  GROUP = 0
C     .......... ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ..........
   50             V = 0.0D0
C
                  DO 60 I = P,Q
                      RV6(I) = UK
                      IF (I.NE.P) THEN
                          IF (DABS(E(I)).LT.DABS(U)) THEN
                              XU = E(I)/U
                              RV4(I) = XU
                              RV1(I-1) = U
                              RV2(I-1) = V
                              RV3(I-1) = 0.0D0
                          ELSE
C     .......... WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ..........
                              XU = U/E(I)
                              RV4(I) = XU
                              RV1(I-1) = E(I)
                              RV2(I-1) = D(I) - X1
                              RV3(I-1) = 0.0D0
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
                  IF (U.EQ.0.0D0) U = EPS3
                  RV1(Q) = U
                  RV2(Q) = 0.0D0
                  RV3(Q) = 0.0D0
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
                          XU = 0.0D0
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
                  NORM = 0.0D0
C
                  DO 130 I = P,Q
                      NORM = NORM + DABS(RV6(I))
  130             CONTINUE
*           INCREMENT OPCOUNT FOR COMPUTING NORM.
                  OPS = OPS + (Q-P+1)
C
                  IF (NORM.GE.1.0D0) THEN
                      GO TO 160
C     .......... FORWARD SUBSTITUTION ..........
                  ELSE IF (ITS.NE.5) THEN
                      IF (NORM.NE.0.0D0) THEN
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
                  XU = 0.0D0
                  GO TO 180
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  160             U = 0.0D0
C
                  DO 170 I = P,Q
                      U = PYTHAG(U,RV6(I))
  170             CONTINUE
C
                  XU = 1.0D0/U
C*PL*ERROR* Embedded comment after label moved
C
  180             DO 190 I = 1,N
                      Z(I,R) = 0.0D0
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
*     CONVERGENCE TEST WAS MODIFIED TO BE THE SAME AS IN DSTEBZ.
*
C
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ITCNT,OPS
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
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION RELFAC
      PARAMETER (RELFAC=2.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS1,LB,UB
      INTEGER IERR,M,M11,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION D(N),E(N),E2(N),RV4(N),RV5(N),W(M)
      INTEGER IND(M)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ATOLI,PIVMIN,RTOLI,SAFEMN,T1,T2,TMP1,TMP2,TNORM,
     +                 U,ULP,V,X0,X1,XU
      INTEGER I,II,ISTURM,J,K,L,M1,M2,M22,P,Q,R,S,TAG
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DLAMCH,EPSLON
      EXTERNAL DLAMCH,EPSLON
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DABS,DMAX1,DMIN1,MAX
C     ..
      ITCNT = 0
      SAFEMN = DLAMCH('S')
      ULP = DLAMCH('E')*DLAMCH('B')
      RTOLI = ULP*RELFAC
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0D0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES ..........
      PIVMIN = ONE
      DO 30 I = 1,N
          X1 = U
          U = 0.0D0
          IF (I.NE.N) U = DABS(E(I+1))
          XU = DMIN1(D(I)- (X1+U),XU)
          X0 = DMAX1(D(I)+ (X1+U),X0)
          IF (I.EQ.1) GO TO 10
CCC         TST1 = DABS(D(I)) + DABS(D(I-1))
CCC         TST2 = TST1 + DABS(E(I))
CCC         IF (TST2 .GT. TST1) GO TO 40
          TMP1 = E(I)**2
          IF (ABS(D(I)*D(I-1))*ULP**2+SAFEMN.LE.TMP1) THEN
              PIVMIN = MAX(PIVMIN,TMP1)
              GO TO 20
          END IF
   10     E2(I) = 0.0D0
   20     CONTINUE
   30 CONTINUE
      PIVMIN = PIVMIN*SAFEMN
      TNORM = MAX(ABS(XU),ABS(X0))
      ATOLI = ULP*TNORM
*        INCREMENT OPCOUNT FOR DETERMINING IF MATRIX SPLITS.
      OPS = OPS + 9* (N-1)
C
      X1 = N
      X1 = X1*EPSLON(DMAX1(DABS(XU),DABS(X0)))
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
      X1 = XU + (X0-XU)*0.5D0
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
      U = 0.0D0
C
      DO 150 Q = P,N
          X1 = U
          U = 0.0D0
          V = 0.0D0
          IF (Q.EQ.N) GO TO 140
          U = DABS(E(Q+1))
          V = E2(Q+1)
  140     XU = DMIN1(D(Q)- (X1+U),XU)
          X0 = DMAX1(D(Q)+ (X1+U),X0)
          IF (V.EQ.0.0D0) GO TO 160
  150 CONTINUE
*        INCREMENT OPCOUNT FOR REFINING INTERVAL.
      OPS = OPS + (N-P+1)*2
C
  160 X1 = EPSLON(DMAX1(DABS(XU),DABS(X0)))
      IF (EPS1.LE.0.0D0) EPS1 = -X1
      IF (P.NE.Q) GO TO 170
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1.GT.D(P) .OR. D(P).GE.T2) GO TO 400
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 340
  170 X1 = X1* (Q-P+1)
      LB = DMAX1(T1,XU-X1)
      UB = DMIN1(T2,X0+X1)
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
  250 X1 = (XU+X0)*0.5D0
CCC         IF ((X0 - XU) .LE. DABS(EPS1)) GO TO 420
CCC         TST1 = 2.0D0 * (DABS(XU) + DABS(X0))
CCC         TST2 = TST1 + (X0 - XU)
CCC         IF (TST2 .EQ. TST1) GO TO 420
      TMP1 = ABS(X0-XU)
      TMP2 = MAX(ABS(X0),ABS(XU))
      IF (TMP1.LT.MAX(ATOLI,PIVMIN,RTOLI*TMP2)) GO TO 330
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  260 S = P - 1
      U = 1.0D0
C
      DO 290 I = P,Q
          IF (U.NE.0.0D0) GO TO 270
          V = DABS(E(I))/EPSLON(1.0D0)
          IF (E2(I).EQ.0.0D0) V = 0.0D0
          GO TO 280
  270     V = E2(I)/U
  280     U = D(I) - X1 - V
          IF (U.LT.0.0D0) S = S + 1
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
      SUBROUTINE DSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
*
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, IOPS IS ONLY INCREMENTED
*     IOPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO IOPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/IOPS,ITCNT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION IOPS,ITCNT
C     ..
C
C
C     DSVDC IS A SUBROUTINE TO REDUCE A DOUBLE PRECISION NXP MATRIX X
C     BY ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
C     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
C     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
C     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
C
C     ON ENTRY
C
C         X         DOUBLE PRECISION(LDX,P), WHERE LDX.GE.N.
C                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
C                   DECOMPOSITION IS TO BE COMPUTED.  X IS
C                   DESTROYED BY DSVDC.
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
C                   LDU IS THE LEADING DIMENSION OF THE ARRAY U.
C                   (SEE BELOW).
C
C         LDV       INTEGER.
C                   LDV IS THE LEADING DIMENSION OF THE ARRAY V.
C                   (SEE BELOW).
C
C         WORK      DOUBLE PRECISION(N).
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
C                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR
C                                  VECTORS IN U.
C                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
C                                  VECTORS.
C                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
C                                  IN V.
C
C     ON RETURN
C
C         S         DOUBLE PRECISION(MM), WHERE MM=MIN(N+1,P).
C                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
C                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
C                   ORDER OF MAGNITUDE.
C
C         E         DOUBLE PRECISION(P),
C                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
C                   DISCUSSION OF INFO FOR EXCEPTIONS.
C
C         U         DOUBLE PRECISION(LDU,K), WHERE LDU.GE.N.  IF
C                                   JOBA.EQ.1 THEN K.EQ.N, IF JOBA.GE.2
C                                   THEN K.EQ.MIN(N,P).
C                   U CONTAINS THE MATRIX OF LEFT SINGULAR VECTORS.
C                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
C                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X
C                   IN THE SUBROUTINE CALL.
C
C         V         DOUBLE PRECISION(LDV,P), WHERE LDV.GE.P.
C                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
C                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N,
C                   THEN V MAY BE IDENTIFIED WITH X IN THE
C                   SUBROUTINE CALL.
C
C         INFO      INTEGER.
C                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
C                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
C                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
C                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
C                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
C                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX
C                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
C                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U)
C                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR
C                   VALUES OF X AND B ARE THE SAME.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C              CORRECTION MADE TO SHIFT 2/84.
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     DSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     EXTERNAL DROT
C     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2,DROTG
C     FORTRAN DABS,DMAX1,MAX0,MIN0,MOD,DSQRT
C
C     INTERNAL VARIABLES
C
*     DOUBLE PRECISION ZTEST,R
*
*     GET EPS FROM DLAMCH FOR NEW STOPPING CRITERION

C     .. Scalar Arguments ..
      INTEGER INFO,JOB,LDU,LDV,LDX,N,P
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION E(*),S(*),U(LDU,*),V(LDV,*),WORK(*),X(LDX,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION B,C,CS,EL,EMM1,EPS,F,G,IOPST,SCALE,SHIFT,SL,SM,
     +                 SMM1,SN,T,T1,TEST
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,MM,
     +        MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      LOGICAL WANTU,WANTV
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DLAMCH,DNRM2
      EXTERNAL DDOT,DLAMCH,DNRM2
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DROT,DROTG,DSCAL,DSWAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DBLE,DMAX1,DSIGN,DSQRT,MAX,MAX0,MIN,MIN0,MOD
C     ..
      IF (N.GT.0 .AND. P.GT.0) THEN
          EPS = DLAMCH('EPSILON')
*
C
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
                      IOPS = IOPS + (2* (N-L+1)+1)
                      S(L) = DNRM2(N-L+1,X(L,L),1)
                      IF (S(L).NE.0.0D0) THEN
                          IF (X(L,L).NE.0.0D0) S(L) = DSIGN(S(L),X(L,L))
*
*              INCREMENT OP COUNT
                          IOPS = IOPS + (N-L+3)
                          CALL DSCAL(N-L+1,1.0D0/S(L),X(L,L),1)
                          X(L,L) = 1.0D0 + X(L,L)
                      END IF
                      S(L) = -S(L)
                  END IF
                  IF (P.GE.LP1) THEN
                      DO 10 J = LP1,P
                          IF (L.LE.NCT) THEN
                              IF (S(L).NE.0.0D0) THEN
C
C              APPLY THE TRANSFORMATION.
C
*
*              INCREMENT OP COUNT
                                  IOPS = IOPS + (4* (N-L)+5)
                                  T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/
     +                                X(L,L)
                                  CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                              END IF
                          END IF
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
                          E(J) = X(L,J)
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
                      IOPS = IOPS + (2* (P-L)+1)
                      E(L) = DNRM2(P-L,E(LP1),1)
                      IF (E(L).NE.0.0D0) THEN
                          IF (E(LP1).NE.0.0D0) E(L) = DSIGN(E(L),E(LP1))
*
*              INCREMENT OP COUNT
                          IOPS = IOPS + (P-L+2)
                          CALL DSCAL(P-L,1.0D0/E(L),E(LP1),1)
                          E(LP1) = 1.0D0 + E(LP1)
                      END IF
                      E(L) = -E(L)
                      IF (LP1.LE.N .AND. E(L).NE.0.0D0) THEN
C
C              APPLY THE TRANSFORMATION.
C
                          DO 30 I = LP1,N
                              WORK(I) = 0.0D0
   30                     CONTINUE
*
*              INCREMENT OP COUNT
                          IOPS = IOPS + DBLE(4* (N-L)+1)* (P-L)
                          DO 40 J = LP1,P
                              CALL DAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),
     +                                   1)
   40                     CONTINUE
                          DO 50 J = LP1,P
                              CALL DAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,
     +                                   X(LP1,J),1)
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
          IF (N.LT.M) S(M) = 0.0D0
          IF (NRTP1.LT.M) E(NRTP1) = X(NRTP1,M)
          E(M) = 0.0D0
C
C     IF REQUIRED, GENERATE U.
C
          IF (WANTU) THEN
              IF (NCU.GE.NCTP1) THEN
                  DO 90 J = NCTP1,NCU
                      DO 80 I = 1,N
                          U(I,J) = 0.0D0
   80                 CONTINUE
                      U(J,J) = 1.0D0
   90             CONTINUE
              END IF
              IF (NCT.GE.1) THEN
                  DO 130 LL = 1,NCT
                      L = NCT - LL + 1
                      IF (S(L).EQ.0.0D0) THEN
                          DO 100 I = 1,N
                              U(I,L) = 0.0D0
  100                     CONTINUE
                          U(L,L) = 1.0D0
                      ELSE
                          LP1 = L + 1
                          IF (NCU.GE.LP1) THEN
*
*              INCREMENT OP COUNT
                              IOPS = IOPS + (DBLE(4* (N-L)+5)* (NCU-L)+
     +                               (N-L+2))
                              DO 110 J = LP1,NCU
                                  T = -DDOT(N-L+1,U(L,L),1,U(L,J),1)/
     +                                U(L,L)
                                  CALL DAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  110                         CONTINUE
                          END IF
                          CALL DSCAL(N-L+1,-1.0D0,U(L,L),1)
                          U(L,L) = 1.0D0 + U(L,L)
                          LM1 = L - 1
                          IF (LM1.GE.1) THEN
                              DO 120 I = 1,LM1
                                  U(I,L) = 0.0D0
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
                      IF (E(L).NE.0.0D0) THEN
*
*              INCREMENT OP COUNT
                          IOPS = IOPS + DBLE(4* (P-L)+1)* (P-L)
                          DO 140 J = LP1,P
                              T = -DDOT(P-L,V(LP1,L),1,V(LP1,J),1)/
     +                            V(LP1,L)
                              CALL DAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  140                     CONTINUE
                      END IF
                  END IF
                  DO 150 I = 1,P
                      V(I,L) = 0.0D0
  150             CONTINUE
                  V(L,L) = 1.0D0
  160         CONTINUE
          END IF
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
          MM = M
*
*     INITIALIZE ITERATION COUNTER
          ITCNT = 0
          ITER = 0
  170     CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
C     ...EXIT
          IF (M.EQ.0) THEN
              GO TO 320
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
                  DO 180 LL = 1,M
                      L = M - LL
C        ...EXIT
                      IF (L.EQ.0) THEN
                          GO TO 200
                      ELSE
*
*           INCREMENT OP COUNT
                          IOPST = IOPST + 2
                          TEST = DABS(S(L)) + DABS(S(L+1))
*
*           REPLACE STOPPING CRITERION WITH NEW ONE AS IN LAPACK
*
*           ZTEST = TEST + DABS(E(L))
*           IF (ZTEST .NE. TEST) GO TO 380
                          IF (DABS(E(L)).LE.EPS*TEST) GO TO 190
                      END IF
  180             CONTINUE
                  GO TO 200
*
  190             E(L) = 0.0D0
C        ......EXIT
  200             IF (L.NE.M-1) THEN
                      LP1 = L + 1
                      MP1 = M + 1
                      DO 210 LLS = LP1,MP1
                          LS = M - LLS + LP1
C           ...EXIT
                          IF (LS.EQ.L) THEN
                              GO TO 230
                          ELSE
                              TEST = 0.0D0
*
*              INCREMENT OP COUNT
                              IOPST = IOPST + 3
                              IF (LS.NE.M) TEST = TEST + DABS(E(LS))
                              IF (LS.NE.L+1) TEST = TEST + DABS(E(LS-1))
*
*              REPLACE STOPPING CRITERION WITH NEW ONE AS IN LAPACK
*
*              ZTEST = TEST + DABS(S(LS))
*              IF (ZTEST .NE. TEST) GO TO 420
                              IF (DABS(S(LS)).LE.EPS*TEST) GO TO 220
                          END IF
  210                 CONTINUE
                      GO TO 230
*
  220                 S(LS) = 0.0D0
C           ......EXIT
  230                 IF (LS.EQ.L) THEN
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
                  GO TO (240,260,280,300) KASE
C
C        DEFLATE NEGLIGIBLE S(M).
C
  240             CONTINUE
                  MM1 = M - 1
                  F = E(M-1)
                  E(M-1) = 0.0D0
*
*           INCREMENT OP COUNT
                  IOPS = IOPS + ((MM1-L+1)*13-2)
                  IF (WANTV) IOPS = IOPS + DBLE(MM1-L+1)*6*P
                  DO 250 KK = L,MM1
                      K = MM1 - KK + L
                      T1 = S(K)
                      CALL DROTG(T1,F,CS,SN)
                      S(K) = T1
                      IF (K.NE.L) THEN
                          F = -SN*E(K-1)
                          E(K-1) = CS*E(K-1)
                      END IF
                      IF (WANTV) CALL DROT(P,V(1,K),1,V(1,M),1,CS,SN)
  250             CONTINUE
                  GO TO 170
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  260             CONTINUE
                  F = E(L-1)
                  E(L-1) = 0.0D0
*
*           INCREMENT OP COUNT
                  IOPS = IOPS + (M-L+1)*13
                  IF (WANTU) IOPS = IOPS + DBLE(M-L+1)*6*N
                  DO 270 K = L,M
                      T1 = S(K)
                      CALL DROTG(T1,F,CS,SN)
                      S(K) = T1
                      F = -SN*E(K)
                      E(K) = CS*E(K)
                      IF (WANTU) CALL DROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  270             CONTINUE
                  GO TO 170
C
C        PERFORM ONE QR STEP.
C
  280             CONTINUE
C
C           CALCULATE THE SHIFT.
C
*
*           INCREMENT OP COUNT
                  IOPST = IOPST + 23
                  SCALE = DMAX1(DABS(S(M)),DABS(S(M-1)),DABS(E(M-1)),
     +                    DABS(S(L)),DABS(E(L)))
                  SM = S(M)/SCALE
                  SMM1 = S(M-1)/SCALE
                  EMM1 = E(M-1)/SCALE
                  SL = S(L)/SCALE
                  EL = E(L)/SCALE
                  B = ((SMM1+SM)* (SMM1-SM)+EMM1**2)/2.0D0
                  C = (SM*EMM1)**2
                  SHIFT = 0.0D0
                  IF (B.NE.0.0D0 .OR. C.NE.0.0D0) THEN
                      SHIFT = DSQRT(B**2+C)
                      IF (B.LT.0.0D0) SHIFT = -SHIFT
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
                  IOPS = IOPS + (MM1-L+1)*38
                  IF (WANTV) IOPS = IOPS + DBLE(MM1-L+1)*6*P
                  IF (WANTU) IOPS = IOPS + DBLE(MAX((MIN(MM1,N-1)-L+1),
     +                              0))*6*N
                  DO 290 K = L,MM1
                      CALL DROTG(F,G,CS,SN)
                      IF (K.NE.L) E(K-1) = F
                      F = CS*S(K) + SN*E(K)
                      E(K) = CS*E(K) - SN*S(K)
                      G = SN*S(K+1)
                      S(K+1) = CS*S(K+1)
                      IF (WANTV) CALL DROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
                      CALL DROTG(F,G,CS,SN)
                      S(K) = F
                      F = CS*E(K) + SN*S(K+1)
                      S(K+1) = -SN*E(K) + CS*S(K+1)
                      G = SN*E(K+1)
                      E(K+1) = CS*E(K+1)
                      IF (WANTU .AND. K.LT.N) CALL DROT(N,U(1,K),1,
     +                    U(1,K+1),1,CS,SN)
  290             CONTINUE
                  E(M-1) = F
                  ITER = ITER + 1
                  GO TO 170
C
C        CONVERGENCE.
C
  300             CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE.
C
                  IF (S(L).LT.0.0D0) THEN
                      S(L) = -S(L)
*
*              INCREMENT OP COUNT
                      IF (WANTV) IOPS = IOPS + P
                      IF (WANTV) CALL DSCAL(P,-1.0D0,V(1,L),1)
                  END IF
  310             CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
                  IF (L.NE.MM) THEN
C           ...EXIT
                      IF (S(L).LT.S(L+1)) THEN
                          T = S(L)
                          S(L) = S(L+1)
                          S(L+1) = T
                          IF (WANTV .AND. L.LT.P) CALL DSWAP(P,V(1,L),1,
     +                        V(1,L+1),1)
                          IF (WANTU .AND. L.LT.N) CALL DSWAP(N,U(1,L),1,
     +                        U(1,L+1),1)
                          L = L + 1
                          GO TO 310
                      END IF
                  END IF
                  ITER = 0
                  M = M - 1
                  GO TO 170
              END IF
          END IF
          INFO = M
C     ......EXIT
*
*     COMPUTE FINAL OPCOUNT
  320     IOPS = IOPS + IOPST
      END IF
      END
      SUBROUTINE QZHES(NM,N,A,B,MATZ,Z)
C
*
*     ---------------------- BEGIN TIMING CODE -------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ITCNT,OPS
C     ..
*     ----------------------- END TIMING CODE --------------------------
*
C
C     THIS SUBROUTINE IS THE FIRST STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM AND THE OTHER
C     TO UPPER TRIANGULAR FORM USING ORTHOGONAL TRANSFORMATIONS.
C     IT IS USUALLY FOLLOWED BY  QZIT,  QZVAL  AND, POSSIBLY,  QZVEC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL GENERAL MATRIX.
C
C        B CONTAINS A REAL GENERAL MATRIX.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO.
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO.
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS IF
C          MATZ HAS BEEN SET TO .TRUE.  OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .......... INITIALIZE Z ..........
C     .. Scalar Arguments ..
      INTEGER N,NM
      LOGICAL MATZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NM,N),B(NM,N),Z(NM,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION R,RHO,S,T,U1,U2,V1,V2
      INTEGER I,J,K,L,L1,LB,NK1,NM1,NM2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DBLE,DSIGN,DSQRT
C     ..
      IF (MATZ) THEN
C
          DO 20 J = 1,N
C
              DO 10 I = 1,N
                  Z(I,J) = 0.0D0
   10         CONTINUE
C
              Z(J,J) = 1.0D0
   20     CONTINUE
      END IF
C     .......... REDUCE B TO UPPER TRIANGULAR FORM ..........
      IF (N.GT.1) THEN
          NM1 = N - 1
C
          DO 120 L = 1,NM1
              L1 = L + 1
              S = 0.0D0
C
              DO 30 I = L1,N
                  S = S + DABS(B(I,L))
   30         CONTINUE
C
              IF (S.NE.0.0D0) THEN
                  S = S + DABS(B(L,L))
                  R = 0.0D0
C
                  DO 40 I = L,N
                      B(I,L) = B(I,L)/S
                      R = R + B(I,L)**2
   40             CONTINUE
C
                  R = DSIGN(DSQRT(R),B(L,L))
                  B(L,L) = B(L,L) + R
                  RHO = R*B(L,L)
C
                  DO 70 J = L1,N
                      T = 0.0D0
C
                      DO 50 I = L,N
                          T = T + B(I,L)*B(I,J)
   50                 CONTINUE
C
                      T = -T/RHO
C
                      DO 60 I = L,N
                          B(I,J) = B(I,J) + T*B(I,L)
   60                 CONTINUE
   70             CONTINUE
C
C
                  DO 100 J = 1,N
                      T = 0.0D0
C
                      DO 80 I = L,N
                          T = T + B(I,L)*A(I,J)
   80                 CONTINUE
C
                      T = -T/RHO
C
                      DO 90 I = L,N
                          A(I,J) = A(I,J) + T*B(I,L)
   90                 CONTINUE
  100             CONTINUE
C
C
                  B(L,L) = -S*R
C
                  DO 110 I = L1,N
                      B(I,L) = 0.0D0
  110             CONTINUE
              END IF
  120     CONTINUE
C
*
*     ---------------------- BEGIN TIMING CODE -------------------------
          OPS = OPS + DBLE(8*N**2+17*N+24)*DBLE(N-1)/3.0D0
*     ----------------------- END TIMING CODE --------------------------
*
C     .......... REDUCE A TO UPPER HESSENBERG FORM, WHILE
C                KEEPING B TRIANGULAR ..........
          IF (N.NE.2) THEN
              NM2 = N - 2
C
              DO 190 K = 1,NM2
                  NK1 = NM1 - K
C     .......... FOR L=N-1 STEP -1 UNTIL K+1 DO -- ..........
                  DO 180 LB = 1,NK1
                      L = N - LB
                      L1 = L + 1
C     .......... ZERO A(L+1,K) ..........
                      S = DABS(A(L,K)) + DABS(A(L1,K))
                      IF (S.NE.0.0D0) THEN
                          U1 = A(L,K)/S
                          U2 = A(L1,K)/S
                          R = DSIGN(DSQRT(U1*U1+U2*U2),U1)
                          V1 = - (U1+R)/R
                          V2 = -U2/R
                          U2 = V2/V1
C
                          DO 130 J = K,N
                              T = A(L,J) + U2*A(L1,J)
                              A(L,J) = A(L,J) + T*V1
                              A(L1,J) = A(L1,J) + T*V2
  130                     CONTINUE
C
                          A(L1,K) = 0.0D0
C
                          DO 140 J = L,N
                              T = B(L,J) + U2*B(L1,J)
                              B(L,J) = B(L,J) + T*V1
                              B(L1,J) = B(L1,J) + T*V2
  140                     CONTINUE
C     .......... ZERO B(L+1,L) ..........
                          S = DABS(B(L1,L1)) + DABS(B(L1,L))
                          IF (S.NE.0.0D0) THEN
                              U1 = B(L1,L1)/S
                              U2 = B(L1,L)/S
                              R = DSIGN(DSQRT(U1*U1+U2*U2),U1)
                              V1 = - (U1+R)/R
                              V2 = -U2/R
                              U2 = V2/V1
C
                              DO 150 I = 1,L1
                                  T = B(I,L1) + U2*B(I,L)
                                  B(I,L1) = B(I,L1) + T*V1
                                  B(I,L) = B(I,L) + T*V2
  150                         CONTINUE
C
                              B(L1,L) = 0.0D0
C
                              DO 160 I = 1,N
                                  T = A(I,L1) + U2*A(I,L)
                                  A(I,L1) = A(I,L1) + T*V1
                                  A(I,L) = A(I,L) + T*V2
  160                         CONTINUE
C
                              IF (MATZ) THEN
C
                                  DO 170 I = 1,N
                                      T = Z(I,L1) + U2*Z(I,L)
                                      Z(I,L1) = Z(I,L1) + T*V1
                                      Z(I,L) = Z(I,L) + T*V2
  170                             CONTINUE
                              END IF
                          END IF
                      END IF
  180             CONTINUE
C
  190         CONTINUE
C
C
*
*     ---------------------- BEGIN TIMING CODE -------------------------
              IF (MATZ) THEN
                  OPS = OPS + DBLE(11*N+20)*DBLE(N-1)*DBLE(N-2)
              ELSE
                  OPS = OPS + DBLE(8*N+20)*DBLE(N-1)*DBLE(N-2)
              END IF
          END IF
      END IF
      RETURN
*$st$ Unreachable comments ...
*     ----------------------- END TIMING CODE --------------------------
*
      END
      SUBROUTINE QZIT(NM,N,A,B,EPS1,MATZ,Z,IERR)
C
*
*     ---------------------- BEGIN TIMING CODE -------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ITCNT,OPS
C     ..
*     ----------------------- END TIMING CODE --------------------------
*
C
C     THIS SUBROUTINE IS THE SECOND STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN D-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM.
C     IT REDUCES THE HESSENBERG MATRIX TO QUASI-TRIANGULAR FORM USING
C     ORTHOGONAL TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX.  IT IS USUALLY PRECEDED BY  QZHES  AND
C     FOLLOWED BY  QZVAL  AND, POSSIBLY,  QZVEC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL UPPER HESSENBERG MATRIX.
C
C        B CONTAINS A REAL UPPER TRIANGULAR MATRIX.
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C        Z CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  QZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT
C
C        A HAS BEEN REDUCED TO QUASI-TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL ARE STILL ZERO AND NO TWO
C          CONSECUTIVE SUBDIAGONAL ELEMENTS ARE NONZERO.
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  THE LOCATION B(N,1) IS USED TO STORE
C          EPS1 TIMES THE NORM OF B FOR LATER USE BY  QZVAL  AND  QZVEC.
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION EPS1
      INTEGER IERR,N,NM
      LOGICAL MATZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NM,N),B(NM,N),Z(NM,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A11,A12,A2,A21,A22,A3,A33,A34,A43,A44,ANI,
     +                 ANORM,B11,B12,B22,B33,B34,B44,BNI,BNORM,EP,EPSA,
     +                 EPSB,OPST,R,S,SH,T,U1,U2,U3,V1,V2,V3
      INTEGER EN,ENM2,ENORN,I,ISH,ITN,ITS,J,K,K1,K2,KM1,L,L1,LD,LL,LM1,
     +        LOR1,NA
      LOGICAL NOTLAS
C     ..
C     .. External Functions ..
      DOUBLE PRECISION EPSLON
      EXTERNAL EPSLON
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DBLE,DSIGN,DSQRT,MAX0,MIN0
C     ..
      IERR = 0
C     .......... COMPUTE EPSA,EPSB ..........
      ANORM = 0.0D0
      BNORM = 0.0D0
C
      DO 20 I = 1,N
          ANI = 0.0D0
          IF (I.NE.1) ANI = DABS(A(I,I-1))
          BNI = 0.0D0
C
          DO 10 J = I,N
              ANI = ANI + DABS(A(I,J))
              BNI = BNI + DABS(B(I,J))
   10     CONTINUE
C
          IF (ANI.GT.ANORM) ANORM = ANI
          IF (BNI.GT.BNORM) BNORM = BNI
   20 CONTINUE
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPS = OPS + DBLE(N* (N+1))
      OPST = 0.0D0
      ITCNT = 0
*     ----------------------- END TIMING CODE --------------------------
*
C
      IF (ANORM.EQ.0.0D0) ANORM = 1.0D0
      IF (BNORM.EQ.0.0D0) BNORM = 1.0D0
      EP = EPS1
      IF (EP.GT.0.0D0) GO TO 30
C     .......... USE ROUNDOFF LEVEL IF EPS1 IS ZERO ..........
      EP = EPSLON(1.0D0)
   30 EPSA = EP*ANORM
      EPSB = EP*BNORM
C     .......... REDUCE A TO QUASI-TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR ..........
      LOR1 = 1
      ENORN = N
      EN = N
      ITN = 30*N
C     .......... BEGIN QZ STEP ..........
   40 IF (EN.LE.2) GO TO 310
      IF (.NOT.MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
   50 ISH = 2
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPS = OPS + OPST
      OPST = 0.0D0
      ITCNT = ITCNT + 1
*     ----------------------- END TIMING CODE --------------------------
*
C     .......... CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- ..........
      DO 60 LL = 1,EN
          LM1 = EN - LL
          L = LM1 + 1
          IF (L.EQ.1) GO TO 80
          IF (DABS(A(L,LM1)).LE.EPSA) GO TO 70
   60 CONTINUE
C
   70 A(L,LM1) = 0.0D0
      IF (L.LT.NA) GO TO 80
C     .......... 1-BY-1 OR 2-BY-2 BLOCK ISOLATED ..........
      EN = LM1
      GO TO 40
C     .......... CHECK FOR SMALL TOP OF B ..........
   80 LD = L
   90 L1 = L + 1
      B11 = B(L,L)
      IF (DABS(B11).GT.EPSB) GO TO 110
      B(L,L) = 0.0D0
      S = DABS(A(L,L)) + DABS(A(L1,L))
      U1 = A(L,L)/S
      U2 = A(L1,L)/S
      R = DSIGN(DSQRT(U1*U1+U2*U2),U1)
      V1 = - (U1+R)/R
      V2 = -U2/R
      U2 = V2/V1
C
      DO 100 J = L,ENORN
          T = A(L,J) + U2*A(L1,J)
          A(L,J) = A(L,J) + T*V1
          A(L1,J) = A(L1,J) + T*V2
          T = B(L,J) + U2*B(L1,J)
          B(L,J) = B(L,J) + T*V1
          B(L1,J) = B(L1,J) + T*V2
  100 CONTINUE
C
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = OPST + DBLE(12* (ENORN+1-L)+11)
*     ----------------------- END TIMING CODE --------------------------
      IF (L.NE.1) A(L,LM1) = -A(L,LM1)
      LM1 = L
      L = L1
      GO TO 70
  110 A11 = A(L,L)/B11
      A21 = A(L1,L)/B11
      IF (ISH.EQ.1) GO TO 130
C     .......... ITERATION STRATEGY ..........
      IF (ITN.EQ.0) GO TO 300
      IF (ITS.EQ.10) GO TO 150
C     .......... DETERMINE TYPE OF SHIFT ..........
      B22 = B(L1,L1)
      IF (DABS(B22).LT.EPSB) B22 = EPSB
      B33 = B(NA,NA)
      IF (DABS(B33).LT.EPSB) B33 = EPSB
      B44 = B(EN,EN)
      IF (DABS(B44).LT.EPSB) B44 = EPSB
      A33 = A(NA,NA)/B33
      A34 = A(NA,EN)/B44
      A43 = A(EN,NA)/B33
      A44 = A(EN,EN)/B44
      B34 = B(NA,EN)/B44
      T = 0.5D0* (A43*B34-A33-A44)
      R = T*T + A34*A43 - A33*A44
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = OPST + DBLE(16)
*     ----------------------- END TIMING CODE --------------------------
      IF (R.LT.0.0D0) GO TO 140
C     .......... DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A ..........
      ISH = 1
      R = DSQRT(R)
      SH = -T + R
      S = -T - R
      IF (DABS(S-A44).LT.DABS(SH-A44)) SH = S
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS OF A.
C                FOR L=EN-2 STEP -1 UNTIL LD DO -- ..........
      DO 120 LL = LD,ENM2
          L = ENM2 + LD - LL
          IF (L.EQ.LD) GO TO 130
          LM1 = L - 1
          L1 = L + 1
          T = A(L,L)
          IF (DABS(B(L,L)).GT.EPSB) T = T - SH*B(L,L)
*        --------------------- BEGIN TIMING CODE -----------------------
          IF (DABS(A(L,LM1)).LE.DABS(T/A(L1,L))*EPSA) THEN
              OPST = OPST + DBLE(5+4* (LL+1-LD))
              GO TO 90
          END IF
  120 CONTINUE
*        ---------------------- END TIMING CODE ------------------------
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = OPST + DBLE(5+4* (ENM2+1-LD))
*     ----------------------- END TIMING CODE --------------------------
C
  130 A1 = A11 - SH
      A2 = A21
      IF (L.NE.LD) A(L,LM1) = -A(L,LM1)
      GO TO 160
C     .......... DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A ..........
  140 A12 = A(L,L1)/B22
      A22 = A(L1,L1)/B22
      B12 = B(L,L1)/B22
      A1 = ((A33-A11)* (A44-A11)-A34*A43+A43*B34*A11)/A21 + A12 -
     +     A11*B12
      A2 = (A22-A11) - A21*B12 - (A33-A11) - (A44-A11) + A43*B34
      A3 = A(L1+1,L1)/B22
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = OPST + DBLE(25)
*     ----------------------- END TIMING CODE --------------------------
      GO TO 160
C     .......... AD HOC SHIFT ..........
  150 A1 = 0.0D0
      A2 = 1.0D0
      A3 = 1.1605D0
  160 ITS = ITS + 1
      ITN = ITN - 1
      IF (.NOT.MATZ) LOR1 = LD
C     .......... MAIN LOOP ..........
      DO 290 K = L,NA
          NOTLAS = K .NE. NA .AND. ISH .EQ. 2
          K1 = K + 1
          K2 = K + 2
          KM1 = MAX0(K-1,L)
          LL = MIN0(EN,K1+ISH)
          IF (NOTLAS) GO TO 190
C     .......... ZERO A(K+1,K-1) ..........
          IF (K.EQ.L) GO TO 170
          A1 = A(K,KM1)
          A2 = A(K1,KM1)
  170     S = DABS(A1) + DABS(A2)
          IF (S.EQ.0.0D0) GO TO 50
          U1 = A1/S
          U2 = A2/S
          R = DSIGN(DSQRT(U1*U1+U2*U2),U1)
          V1 = - (U1+R)/R
          V2 = -U2/R
          U2 = V2/V1
C
          DO 180 J = KM1,ENORN
              T = A(K,J) + U2*A(K1,J)
              A(K,J) = A(K,J) + T*V1
              A(K1,J) = A(K1,J) + T*V2
              T = B(K,J) + U2*B(K1,J)
              B(K,J) = B(K,J) + T*V1
              B(K1,J) = B(K1,J) + T*V2
  180     CONTINUE
C
*        --------------------- BEGIN TIMING CODE -----------------------
          OPST = OPST + DBLE(11+12* (ENORN+1-KM1))
*        ---------------------- END TIMING CODE ------------------------
          IF (K.NE.L) A(K1,KM1) = 0.0D0
          GO TO 250
C     .......... ZERO A(K+1,K-1) AND A(K+2,K-1) ..........
  190     IF (K.EQ.L) GO TO 200
          A1 = A(K,KM1)
          A2 = A(K1,KM1)
          A3 = A(K2,KM1)
  200     S = DABS(A1) + DABS(A2) + DABS(A3)
          IF (S.EQ.0.0D0) GO TO 280
          U1 = A1/S
          U2 = A2/S
          U3 = A3/S
          R = DSIGN(DSQRT(U1*U1+U2*U2+U3*U3),U1)
          V1 = - (U1+R)/R
          V2 = -U2/R
          V3 = -U3/R
          U2 = V2/V1
          U3 = V3/V1
C
          DO 210 J = KM1,ENORN
              T = A(K,J) + U2*A(K1,J) + U3*A(K2,J)
              A(K,J) = A(K,J) + T*V1
              A(K1,J) = A(K1,J) + T*V2
              A(K2,J) = A(K2,J) + T*V3
              T = B(K,J) + U2*B(K1,J) + U3*B(K2,J)
              B(K,J) = B(K,J) + T*V1
              B(K1,J) = B(K1,J) + T*V2
              B(K2,J) = B(K2,J) + T*V3
  210     CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
          OPST = OPST + DBLE(17+20* (ENORN+1-KM1))
*        ---------------------- END TIMING CODE ------------------------
C
          IF (K.EQ.L) GO TO 220
          A(K1,KM1) = 0.0D0
          A(K2,KM1) = 0.0D0
C     .......... ZERO B(K+2,K+1) AND B(K+2,K) ..........
  220     S = DABS(B(K2,K2)) + DABS(B(K2,K1)) + DABS(B(K2,K))
          IF (S.EQ.0.0D0) GO TO 250
          U1 = B(K2,K2)/S
          U2 = B(K2,K1)/S
          U3 = B(K2,K)/S
          R = DSIGN(DSQRT(U1*U1+U2*U2+U3*U3),U1)
          V1 = - (U1+R)/R
          V2 = -U2/R
          V3 = -U3/R
          U2 = V2/V1
          U3 = V3/V1
C
          DO 230 I = LOR1,LL
              T = A(I,K2) + U2*A(I,K1) + U3*A(I,K)
              A(I,K2) = A(I,K2) + T*V1
              A(I,K1) = A(I,K1) + T*V2
              A(I,K) = A(I,K) + T*V3
              T = B(I,K2) + U2*B(I,K1) + U3*B(I,K)
              B(I,K2) = B(I,K2) + T*V1
              B(I,K1) = B(I,K1) + T*V2
              B(I,K) = B(I,K) + T*V3
  230     CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
          OPST = OPST + DBLE(17+20* (LL+1-LOR1))
*        ---------------------- END TIMING CODE ------------------------
C
          B(K2,K) = 0.0D0
          B(K2,K1) = 0.0D0
          IF (.NOT.MATZ) GO TO 250
C
          DO 240 I = 1,N
              T = Z(I,K2) + U2*Z(I,K1) + U3*Z(I,K)
              Z(I,K2) = Z(I,K2) + T*V1
              Z(I,K1) = Z(I,K1) + T*V2
              Z(I,K) = Z(I,K) + T*V3
  240     CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
          OPST = OPST + DBLE(10*N)
*        ---------------------- END TIMING CODE ------------------------
C     .......... ZERO B(K+1,K) ..........
  250     S = DABS(B(K1,K1)) + DABS(B(K1,K))
          IF (S.EQ.0.0D0) GO TO 280
          U1 = B(K1,K1)/S
          U2 = B(K1,K)/S
          R = DSIGN(DSQRT(U1*U1+U2*U2),U1)
          V1 = - (U1+R)/R
          V2 = -U2/R
          U2 = V2/V1
C
          DO 260 I = LOR1,LL
              T = A(I,K1) + U2*A(I,K)
              A(I,K1) = A(I,K1) + T*V1
              A(I,K) = A(I,K) + T*V2
              T = B(I,K1) + U2*B(I,K)
              B(I,K1) = B(I,K1) + T*V1
              B(I,K) = B(I,K) + T*V2
  260     CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
          OPST = OPST + DBLE(11+12* (LL+1-LOR1))
*        ---------------------- END TIMING CODE ------------------------
C
          B(K1,K) = 0.0D0
          IF (.NOT.MATZ) GO TO 280
C
          DO 270 I = 1,N
              T = Z(I,K1) + U2*Z(I,K)
              Z(I,K1) = Z(I,K1) + T*V1
              Z(I,K) = Z(I,K) + T*V2
  270     CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
          OPST = OPST + DBLE(6*N)
  280     CONTINUE
  290 CONTINUE
*        ---------------------- END TIMING CODE ------------------------
C
C     .......... END QZ STEP ..........
      GO TO 50
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
  300 IERR = EN
C     .......... SAVE EPSB FOR USE BY QZVAL AND QZVEC ..........
  310 IF (N.GT.1) B(N,1) = EPSB
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPS = OPS + OPST
      OPST = 0.0D0
*     ----------------------- END TIMING CODE --------------------------
*
      RETURN
      END
      SUBROUTINE QZVAL(NM,N,A,B,ALFR,ALFI,BETA,MATZ,Z)
C
*
*     ---------------------- BEGIN TIMING CODE -------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ITCNT,OPS
C     ..
*     ----------------------- END TIMING CODE --------------------------
*
C
C     THIS SUBROUTINE IS THE THIRD STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM
C     IN QUASI-TRIANGULAR FORM AND THE OTHER IN UPPER TRIANGULAR FORM.
C     IT REDUCES THE QUASI-TRIANGULAR MATRIX FURTHER, SO THAT ANY
C     REMAINING 2-BY-2 BLOCKS CORRESPOND TO PAIRS OF COMPLEX
C     EIGENVALUES, AND RETURNS QUANTITIES WHOSE RATIOS GIVE THE
C     GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY  QZHES
C     AND  QZIT  AND MAY BE FOLLOWED BY  QZVEC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL UPPER QUASI-TRIANGULAR MATRIX.
C
C        B CONTAINS A REAL UPPER TRIANGULAR MATRIX.  IN ADDITION,
C          LOCATION B(N,1) CONTAINS THE TOLERANCE QUANTITY (EPSB)
C          COMPUTED AND SAVED IN  QZIT.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C        Z CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTIONS BY QZHES
C          AND QZIT, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT
C
C        A HAS BEEN REDUCED FURTHER TO A QUASI-TRIANGULAR MATRIX
C          IN WHICH ALL NONZERO SUBDIAGONAL ELEMENTS CORRESPOND TO
C          PAIRS OF COMPLEX EIGENVALUES.
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  B(N,1) IS UNALTERED.
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULAR MATRIX THAT WOULD BE
C          OBTAINED IF A WERE REDUCED COMPLETELY TO TRIANGULAR FORM
C          BY UNITARY TRANSFORMATIONS.  NON-ZERO VALUES OF ALFI OCCUR
C          IN PAIRS, THE FIRST MEMBER POSITIVE AND THE SECOND NEGATIVE.
C
C        BETA CONTAINS THE DIAGONAL ELEMENTS OF THE CORRESPONDING B,
C          NORMALIZED TO BE REAL AND NON-NEGATIVE.  THE GENERALIZED
C          EIGENVALUES ARE THEN THE RATIOS ((ALFR+I*ALFI)/BETA).
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR ALL THREE STEPS) IF MATZ HAS BEEN SET TO .TRUE.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER N,NM
      LOGICAL MATZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NM,N),ALFI(N),ALFR(N),B(NM,N),BETA(N),Z(NM,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A11,A11I,A11R,A12,A12I,A12R,A1I,A2,A21,A22,
     +                 A22I,A22R,A2I,AN,B11,B12,B22,BN,C,CQ,CZ,D,DI,DR,
     +                 E,EI,EPSB,OPST,OPST2,R,S,SQI,SQR,SSI,SSR,SZI,SZR,
     +                 T,TI,TR,U1,U2,V1,V2
      INTEGER EN,I,ISW,J,NA,NN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DBLE,DSIGN,DSQRT
C     ..
      EPSB = B(N,1)
      ISW = 1
C     .......... FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES.
C                FOR EN=N STEP -1 UNTIL 1 DO -- ..........
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPST = 0.0D0
      OPST2 = 0.0D0
*     ----------------------- END TIMING CODE --------------------------
*
      DO 80 NN = 1,N
*
*        --------------------- BEGIN TIMING CODE -----------------------
          OPST = OPST + OPST2
          OPST2 = 0.0D0
*        ---------------------- END TIMING CODE ------------------------
*
          EN = N + 1 - NN
          NA = EN - 1
          IF (ISW.NE.2) THEN
              IF (EN.NE.1) THEN
                  IF (A(EN,NA).NE.0.0D0) THEN
C     .......... 2-BY-2 BLOCK ..........
                      IF (DABS(B(NA,NA)).GT.EPSB) THEN
                          IF (DABS(B(EN,EN)).GT.EPSB) THEN
                              AN = DABS(A(NA,NA)) + DABS(A(NA,EN)) +
     +                             DABS(A(EN,NA)) + DABS(A(EN,EN))
                              BN = DABS(B(NA,NA)) + DABS(B(NA,EN)) +
     +                             DABS(B(EN,EN))
                              A11 = A(NA,NA)/AN
                              A12 = A(NA,EN)/AN
                              A21 = A(EN,NA)/AN
                              A22 = A(EN,EN)/AN
                              B11 = B(NA,NA)/BN
                              B12 = B(NA,EN)/BN
                              B22 = B(EN,EN)/BN
                              E = A11/B11
                              EI = A22/B22
                              S = A21/ (B11*B22)
                              T = (A22-E*B22)/B22
                              IF (DABS(E).GT.DABS(EI)) THEN
                                  E = EI
                                  T = (A11-E*B11)/B11
                              END IF
                              C = 0.5D0* (T-S*B12)
                              D = C*C + S* (A12-E*B12)
*        --------------------- BEGIN TIMING CODE -----------------------
                              OPST2 = OPST2 + DBLE(28)
*        ---------------------- END TIMING CODE ------------------------
                              IF (D.LT.0.0D0) THEN
C     .......... TWO COMPLEX ROOTS ..........
                                  E = E + C
                                  EI = DSQRT(-D)
                                  A11R = A11 - E*B11
                                  A11I = EI*B11
                                  A12R = A12 - E*B12
                                  A12I = EI*B12
                                  A22R = A22 - E*B22
                                  A22I = EI*B22
                                  IF (DABS(A11R)+DABS(A11I)+DABS(A12R)+
     +                                DABS(A12I).LT.DABS(A21)+
     +                                DABS(A22R)+DABS(A22I)) THEN
                                      A1 = A22R
                                      A1I = A22I
                                      A2 = -A21
                                      A2I = 0.0D0
                                  ELSE
                                      A1 = A12R
                                      A1I = A12I
                                      A2 = -A11R
                                      A2I = -A11I
                                  END IF
C     .......... CHOOSE COMPLEX Z ..........
                                  CZ = DSQRT(A1*A1+A1I*A1I)
                                  IF (CZ.EQ.0.0D0) THEN
                                      SZR = 1.0D0
                                      SZI = 0.0D0
                                  ELSE
                                      SZR = (A1*A2+A1I*A2I)/CZ
                                      SZI = (A1*A2I-A1I*A2)/CZ
                                      R = DSQRT(CZ*CZ+SZR*SZR+SZI*SZI)
                                      CZ = CZ/R
                                      SZR = SZR/R
                                      SZI = SZI/R
                                  END IF
                                  IF (AN.LT. (DABS(E)+EI)*BN) THEN
                                      A1 = CZ*A11 + SZR*A12
                                      A1I = SZI*A12
                                      A2 = CZ*A21 + SZR*A22
                                      A2I = SZI*A22
                                  ELSE
                                      A1 = CZ*B11 + SZR*B12
                                      A1I = SZI*B12
                                      A2 = SZR*B22
                                      A2I = SZI*B22
                                  END IF
C     .......... CHOOSE COMPLEX Q ..........
                                  CQ = DSQRT(A1*A1+A1I*A1I)
                                  IF (CQ.EQ.0.0D0) THEN
                                      SQR = 1.0D0
                                      SQI = 0.0D0
                                  ELSE
                                      SQR = (A1*A2+A1I*A2I)/CQ
                                      SQI = (A1*A2I-A1I*A2)/CQ
                                      R = DSQRT(CQ*CQ+SQR*SQR+SQI*SQI)
                                      CQ = CQ/R
                                      SQR = SQR/R
                                      SQI = SQI/R
                                  END IF
C     .......... COMPUTE DIAGONAL ELEMENTS THAT WOULD RESULT
C                IF TRANSFORMATIONS WERE APPLIED ..........
                                  SSR = SQR*SZR + SQI*SZI
                                  SSI = SQR*SZI - SQI*SZR
                                  I = 1
                                  TR = CQ*CZ*A11 + CQ*SZR*A12 +
     +                                 SQR*CZ*A21 + SSR*A22
                                  TI = CQ*SZI*A12 - SQI*CZ*A21 + SSI*A22
                                  DR = CQ*CZ*B11 + CQ*SZR*B12 + SSR*B22
                                  DI = CQ*SZI*B12 + SSI*B22
   10                             CONTINUE
                                  T = TI*DR - TR*DI
                                  J = NA
                                  IF (T.LT.0.0D0) J = EN
                                  R = DSQRT(DR*DR+DI*DI)
                                  BETA(J) = BN*R
                                  ALFR(J) = AN* (TR*DR+TI*DI)/R
                                  ALFI(J) = AN*T/R
                                  IF (I.EQ.1) THEN
                                      I = 2
                                      TR = SSR*A11 - SQR*CZ*A12 -
     +                                     CQ*SZR*A21 + CQ*CZ*A22
                                      TI = -SSI*A11 - SQI*CZ*A12 +
     +                                     CQ*SZI*A21
                                      DR = SSR*B11 - SQR*CZ*B12 +
     +                                     CQ*CZ*B22
                                      DI = -SSI*B11 - SQI*CZ*B12
                                      GO TO 10
                                  END IF
*        --------------------- BEGIN TIMING CODE -----------------------
                                  OPST2 = OPST2 + DBLE(151)
                                  GO TO 70
                              ELSE
C     .......... TWO REAL ROOTS.
C                ZERO BOTH A(EN,NA) AND B(EN,NA) ..........
                                  E = E + (C+DSIGN(DSQRT(D),C))
                                  A11 = A11 - E*B11
                                  A12 = A12 - E*B12
                                  A22 = A22 - E*B22
*        --------------------- BEGIN TIMING CODE -----------------------
                                  OPST2 = OPST2 + DBLE(11)
*        ---------------------- END TIMING CODE ------------------------
                                  IF (DABS(A11)+DABS(A12).LT.
     +                                DABS(A21)+DABS(A22)) THEN
                                      A1 = A22
                                      A2 = A21
                                  ELSE
                                      A1 = A12
                                      A2 = A11
                                  END IF
                              END IF
                          ELSE
                              A1 = A(EN,EN)
                              A2 = A(EN,NA)
                              BN = 0.0D0
                          END IF
C     .......... CHOOSE AND APPLY REAL Z ..........
                          S = DABS(A1) + DABS(A2)
                          U1 = A1/S
                          U2 = A2/S
                          R = DSIGN(DSQRT(U1*U1+U2*U2),U1)
                          V1 = - (U1+R)/R
                          V2 = -U2/R
                          U2 = V2/V1
C
                          DO 20 I = 1,EN
                              T = A(I,EN) + U2*A(I,NA)
                              A(I,EN) = A(I,EN) + T*V1
                              A(I,NA) = A(I,NA) + T*V2
                              T = B(I,EN) + U2*B(I,NA)
                              B(I,EN) = B(I,EN) + T*V1
                              B(I,NA) = B(I,NA) + T*V2
   20                     CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
                          OPST2 = OPST2 + DBLE(11+12*EN)
*        ---------------------- END TIMING CODE ------------------------
C
                          IF (MATZ) THEN
C
                              DO 30 I = 1,N
                                  T = Z(I,EN) + U2*Z(I,NA)
                                  Z(I,EN) = Z(I,EN) + T*V1
                                  Z(I,NA) = Z(I,NA) + T*V2
   30                         CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
                              OPST2 = OPST2 + DBLE(6*N)
                          END IF
*        ---------------------- END TIMING CODE ------------------------
C
                          IF (BN.EQ.0.0D0) THEN
                              GO TO 60
                          ELSE IF (AN.GE.DABS(E)*BN) THEN
                              A1 = B(NA,NA)
                              A2 = B(EN,NA)
                              GO TO 40
                          END IF
                      END IF
                      A1 = A(NA,NA)
                      A2 = A(EN,NA)
C     .......... CHOOSE AND APPLY REAL Q ..........
   40                 S = DABS(A1) + DABS(A2)
                      IF (S.NE.0.0D0) THEN
                          U1 = A1/S
                          U2 = A2/S
                          R = DSIGN(DSQRT(U1*U1+U2*U2),U1)
                          V1 = - (U1+R)/R
                          V2 = -U2/R
                          U2 = V2/V1
C
                          DO 50 J = NA,N
                              T = A(NA,J) + U2*A(EN,J)
                              A(NA,J) = A(NA,J) + T*V1
                              A(EN,J) = A(EN,J) + T*V2
                              T = B(NA,J) + U2*B(EN,J)
                              B(NA,J) = B(NA,J) + T*V1
                              B(EN,J) = B(EN,J) + T*V2
   50                     CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
                          OPST2 = OPST2 + DBLE(11+12* (N+1-NA))
                      END IF
*        ---------------------- END TIMING CODE ------------------------
C
   60                 A(EN,NA) = 0.0D0
                      B(EN,NA) = 0.0D0
                      ALFR(NA) = A(NA,NA)
                      ALFR(EN) = A(EN,EN)
                      IF (B(NA,NA).LT.0.0D0) ALFR(NA) = -ALFR(NA)
                      IF (B(EN,EN).LT.0.0D0) ALFR(EN) = -ALFR(EN)
                      BETA(NA) = DABS(B(NA,NA))
                      BETA(EN) = DABS(B(EN,EN))
                      ALFI(EN) = 0.0D0
                      ALFI(NA) = 0.0D0
                      GO TO 70
                  END IF
              END IF
C     .......... 1-BY-1 BLOCK, ONE REAL ROOT ..........
              ALFR(EN) = A(EN,EN)
              IF (B(EN,EN).LT.0.0D0) ALFR(EN) = -ALFR(EN)
              BETA(EN) = DABS(B(EN,EN))
              ALFI(EN) = 0.0D0
              GO TO 80
          END IF
*        ---------------------- END TIMING CODE ------------------------
   70     ISW = 3 - ISW
   80 CONTINUE
*
*     ---------------------- BEGIN TIMING CODE -------------------------
      OPS = OPS + (OPST+OPST2)
*     ----------------------- END TIMING CODE --------------------------
*
      B(N,1) = EPSB
C
      END
      SUBROUTINE QZVEC(NM,N,A,B,ALFR,ALFI,BETA,Z)
C
*
*     ---------------------- BEGIN TIMING CODE -------------------------
*     COMMON BLOCK TO RETURN OPERATION COUNT AND ITERATION COUNT
*     ITCNT IS INITIALIZED TO 0, OPS IS ONLY INCREMENTED
*     OPST IS USED TO ACCUMULATE SMALL CONTRIBUTIONS TO OPS
*     TO AVOID ROUNDOFF ERROR
C     .. Common blocks ..
      COMMON /LATIME/OPS,ITCNT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ITCNT,OPS
C     ..
*     ----------------------- END TIMING CODE --------------------------
*
C
C     THIS SUBROUTINE IS THE OPTIONAL FOURTH STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM IN
C     QUASI-TRIANGULAR FORM (IN WHICH EACH 2-BY-2 BLOCK CORRESPONDS TO
C     A PAIR OF COMPLEX EIGENVALUES) AND THE OTHER IN UPPER TRIANGULAR
C     FORM.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM AND
C     TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  QZHES,  QZIT, AND  QZVAL.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRICES.
C
C        A CONTAINS A REAL UPPER QUASI-TRIANGULAR MATRIX.
C
C        B CONTAINS A REAL UPPER TRIANGULAR MATRIX.  IN ADDITION,
C          LOCATION B(N,1) CONTAINS THE TOLERANCE QUANTITY (EPSB)
C          COMPUTED AND SAVED IN  QZIT.
C
C        ALFR, ALFI, AND BETA  ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  QZVAL.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  QZHES,  QZIT, AND  QZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT
C
C        A IS UNALTERED.  ITS SUBDIAGONAL ELEMENTS PROVIDE INFORMATION
C           ABOUT THE STORAGE OF THE COMPLEX EIGENVECTORS.
C
C        B HAS BEEN DESTROYED.
C
C        ALFR, ALFI, AND BETA ARE UNALTERED.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF ALFI(I) .EQ. 0.0, THE I-TH EIGENVALUE IS REAL AND
C            THE I-TH COLUMN OF Z CONTAINS ITS EIGENVECTOR.
C          IF ALFI(I) .NE. 0.0, THE I-TH EIGENVALUE IS COMPLEX.
C            IF ALFI(I) .GT. 0.0, THE EIGENVALUE IS THE FIRST OF
C              A COMPLEX PAIR AND THE I-TH AND (I+1)-TH COLUMNS
C              OF Z CONTAIN ITS EIGENVECTOR.
C            IF ALFI(I) .LT. 0.0, THE EIGENVALUE IS THE SECOND OF
C              A COMPLEX PAIR AND THE (I-1)-TH AND I-TH COLUMNS
C              OF Z CONTAIN THE CONJUGATE OF ITS EIGENVECTOR.
C          EACH EIGENVECTOR IS NORMALIZED SO THAT THE MODULUS
C          OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER N,NM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NM,N),ALFI(N),ALFR(N),B(NM,N),BETA(N),Z(NM,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALFM,ALMI,ALMR,BETM,D,DI,DR,EPSB,Q,R,RA,RR,S,SA,
     +                 T,T1,T2,TI,TR,W,W1,X,X1,Y,Z1,ZZ
      INTEGER EN,ENM2,I,II,IN2BY2,ISW,J,JJ,K,M,NA,NN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DBLE,DSQRT
C     ..
      EPSB = B(N,1)
      ISW = 1
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 220 NN = 1,N
*        --------------------- BEGIN TIMING CODE -----------------------
          IN2BY2 = 0
*        ---------------------- END TIMING CODE ------------------------
          EN = N + 1 - NN
          NA = EN - 1
          IF (ISW.EQ.2) GO TO 200
          IF (ALFI(EN).NE.0.0D0) GO TO 80
C     .......... REAL VECTOR ..........
          M = EN
          B(EN,EN) = 1.0D0
          IF (NA.EQ.0) GO TO 210
          ALFM = ALFR(M)
          BETM = BETA(M)
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
          DO 70 II = 1,NA
              I = EN - II
              W = BETM*A(I,I) - ALFM*B(I,I)
              R = 0.0D0
C
              DO 10 J = M,EN
                  R = R + (BETM*A(I,J)-ALFM*B(I,J))*B(J,EN)
   10         CONTINUE
C
              IF (I.EQ.1 .OR. ISW.EQ.2) GO TO 20
              IF (BETM*A(I,I-1).EQ.0.0D0) GO TO 20
              ZZ = W
              S = R
              GO TO 50
   20         M = I
              IF (ISW.EQ.2) GO TO 30
C     .......... REAL 1-BY-1 BLOCK ..........
              T = W
              IF (W.EQ.0.0D0) T = EPSB
              B(I,EN) = -R/T
              GO TO 60
C     .......... REAL 2-BY-2 BLOCK ..........
   30         X = BETM*A(I,I+1) - ALFM*B(I,I+1)
              Y = BETM*A(I+1,I)
              Q = W*ZZ - X*Y
              T = (X*S-ZZ*R)/Q
              B(I,EN) = T
*           ------------------- BEGIN TIMING CODE ----------------------
              IN2BY2 = IN2BY2 + 1
*           -------------------- END TIMING CODE -----------------------
              IF (DABS(X).LE.DABS(ZZ)) GO TO 40
              B(I+1,EN) = (-R-W*T)/X
              GO TO 50
   40         B(I+1,EN) = (-S-Y*T)/ZZ
   50         ISW = 3 - ISW
   60         CONTINUE
   70     CONTINUE
C     .......... END REAL VECTOR ..........
*        --------------------- BEGIN TIMING CODE -----------------------
          OPS = OPS + (5.0D0/2.0D0)*DBLE((EN+2)* (EN-1)+IN2BY2)
*        ---------------------- END TIMING CODE ------------------------
          GO TO 210
C     .......... COMPLEX VECTOR ..........
   80     M = NA
          ALMR = ALFR(M)
          ALMI = ALFI(M)
          BETM = BETA(M)
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
          Y = BETM*A(EN,NA)
          B(NA,NA) = -ALMI*B(EN,EN)/Y
          B(NA,EN) = (ALMR*B(EN,EN)-BETM*A(EN,EN))/Y
          B(EN,NA) = 0.0D0
          B(EN,EN) = 1.0D0
          ENM2 = NA - 1
          IF (ENM2.EQ.0) GO TO 200
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
          DO 190 II = 1,ENM2
              I = NA - II
              W = BETM*A(I,I) - ALMR*B(I,I)
              W1 = -ALMI*B(I,I)
              RA = 0.0D0
              SA = 0.0D0
C
              DO 90 J = M,EN
                  X = BETM*A(I,J) - ALMR*B(I,J)
                  X1 = -ALMI*B(I,J)
                  RA = RA + X*B(J,NA) - X1*B(J,EN)
                  SA = SA + X*B(J,EN) + X1*B(J,NA)
   90         CONTINUE
C
              IF (I.EQ.1 .OR. ISW.EQ.2) GO TO 100
              IF (BETM*A(I,I-1).EQ.0.0D0) GO TO 100
              ZZ = W
              Z1 = W1
              R = RA
              S = SA
              ISW = 2
              GO TO 180
  100         M = I
              IF (ISW.EQ.2) GO TO 140
C     .......... COMPLEX 1-BY-1 BLOCK ..........
              TR = -RA
              TI = -SA
  110         DR = W
              DI = W1
C     .......... COMPLEX DIVIDE (T1,T2) = (TR,TI) / (DR,DI) ..........
  120         IF (DABS(DI).GT.DABS(DR)) GO TO 130
              RR = DI/DR
              D = DR + DI*RR
              T1 = (TR+TI*RR)/D
              T2 = (TI-TR*RR)/D
              GO TO (170,150) ISW
  130         RR = DR/DI
              D = DR*RR + DI
              T1 = (TR*RR+TI)/D
              T2 = (TI*RR-TR)/D
              GO TO (170,150) ISW
C     .......... COMPLEX 2-BY-2 BLOCK ..........
  140         X = BETM*A(I,I+1) - ALMR*B(I,I+1)
              X1 = -ALMI*B(I,I+1)
              Y = BETM*A(I+1,I)
              TR = Y*RA - W*R + W1*S
              TI = Y*SA - W*S - W1*R
              DR = W*ZZ - W1*Z1 - X*Y
              DI = W*Z1 + W1*ZZ - X1*Y
*           ------------------- BEGIN TIMING CODE ----------------------
              IN2BY2 = IN2BY2 + 1
*           -------------------- END TIMING CODE -----------------------
              IF (DR.EQ.0.0D0 .AND. DI.EQ.0.0D0) DR = EPSB
              GO TO 120
  150         B(I+1,NA) = T1
              B(I+1,EN) = T2
              ISW = 1
              IF (DABS(Y).GT.DABS(W)+DABS(W1)) GO TO 160
              TR = -RA - X*B(I+1,NA) + X1*B(I+1,EN)
              TI = -SA - X*B(I+1,EN) - X1*B(I+1,NA)
              GO TO 110
  160         T1 = (-R-ZZ*B(I+1,NA)+Z1*B(I+1,EN))/Y
              T2 = (-S-ZZ*B(I+1,EN)-Z1*B(I+1,NA))/Y
  170         B(I,NA) = T1
              B(I,EN) = T2
  180         CONTINUE
  190     CONTINUE
*        --------------------- BEGIN TIMING CODE -----------------------
          OPS = OPS + DBLE((6*EN-7)* (EN-2)+31*IN2BY2)
*        ---------------------- END TIMING CODE ------------------------
C     .......... END COMPLEX VECTOR ..........
  200     ISW = 3 - ISW
  210     CONTINUE
  220 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 1 DO -- ..........
      DO 250 JJ = 1,N
          J = N + 1 - JJ
C
          DO 240 I = 1,N
              ZZ = 0.0D0
C
              DO 230 K = 1,J
                  ZZ = ZZ + Z(I,K)*B(K,J)
  230         CONTINUE
C
              Z(I,J) = ZZ
  240     CONTINUE
  250 CONTINUE
*     ----------------------- BEGIN TIMING CODE ------------------------
      OPS = OPS + DBLE(N**2)*DBLE(N+1)
*     ------------------------ END TIMING CODE -------------------------
C     .......... NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1.
C                (ISW IS 1 INITIALLY FROM BEFORE) ..........
*     ------------------------ BEGIN TIMING CODE -----------------------
      IN2BY2 = 0
*     ------------------------- END TIMING CODE ------------------------
      DO 330 J = 1,N
          D = 0.0D0
          IF (ISW.EQ.2) GO TO 280
          IF (ALFI(J).NE.0.0D0) GO TO 310
C
          DO 260 I = 1,N
              IF (DABS(Z(I,J)).GT.D) D = DABS(Z(I,J))
  260     CONTINUE
C
          DO 270 I = 1,N
              Z(I,J) = Z(I,J)/D
  270     CONTINUE
C
          GO TO 320
C*PL*ERROR* Embedded comment after label moved
C
  280     DO 290 I = 1,N
              R = DABS(Z(I,J-1)) + DABS(Z(I,J))
              IF (R.NE.0.0D0) R = R*DSQRT((Z(I,J-1)/R)**2+
     +                            (Z(I,J)/R)**2)
              IF (R.GT.D) D = R
  290     CONTINUE
C
          DO 300 I = 1,N
              Z(I,J-1) = Z(I,J-1)/D
              Z(I,J) = Z(I,J)/D
  300     CONTINUE
*        ---------------------- BEGIN TIMING CODE ----------------------
          IN2BY2 = IN2BY2 + 1
*        ----------------------- END TIMING CODE -----------------------
C
  310     ISW = 3 - ISW
  320     CONTINUE
  330 CONTINUE
*     ------------------------ BEGIN TIMING CODE -----------------------
      OPS = OPS + DBLE(N* (N+5*IN2BY2))
*     ------------------------- END TIMING CODE ------------------------
C
      RETURN
      END
