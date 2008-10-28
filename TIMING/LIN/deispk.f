      SUBROUTINE ORTHES(NM,N,LOW,IGH,A,ORT)
C
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
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA.GE.KP1) THEN
C
          DO 90 M = KP1,LA
              H = 0.0D0
              ORT(M) = 0.0D0
              SCALE = 0.0D0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
              DO 10 I = M,IGH
                  SCALE = SCALE + DABS(A(I,M-1))
   10         CONTINUE
C
              IF (SCALE.NE.0.0D0) THEN
                  MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
                  DO 20 II = M,IGH
                      I = MP - II
                      ORT(I) = A(I,M-1)/SCALE
                      H = H + ORT(I)*ORT(I)
   20             CONTINUE
C
                  G = -DSIGN(DSQRT(H),ORT(M))
                  H = H - ORT(M)*G
                  ORT(M) = ORT(M) - G
C     .......... FORM (I-(U*UT)/H) * A ..........
                  DO 50 J = M,N
                      F = 0.0D0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
                      DO 30 II = M,IGH
                          I = MP - II
                          F = F + ORT(I)*A(I,J)
   30                 CONTINUE
C
                      F = F/H
C
                      DO 40 I = M,IGH
                          A(I,J) = A(I,J) - F*ORT(I)
   40                 CONTINUE
   50             CONTINUE
C
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
                  DO 80 I = 1,IGH
                      F = 0.0D0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
                      DO 60 JJ = M,IGH
                          J = MP - JJ
                          F = F + ORT(J)*A(I,J)
   60                 CONTINUE
C
                      F = F/H
C
                      DO 70 J = M,IGH
                          A(I,J) = A(I,J) - F*ORT(J)
   70                 CONTINUE
   80             CONTINUE
C
C
                  ORT(M) = SCALE*ORT(M)
                  A(M,M-1) = SCALE*G
              END IF
   90     CONTINUE
      END IF
      RETURN
*$st$ Unreachable comments ...
C
      END
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
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
      INTRINSIC DABS,DSIGN,DSQRT
C     ..
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
