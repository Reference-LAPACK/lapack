      SUBROUTINE CGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(*),INFO
      COMPLEX A(LDA,*)
C
C     CGEFA FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION.
C
C     CGEFA IS USUALLY CALLED BY CGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR CGECO) = (1 + 9/N)*(TIME FOR CGEFA) .
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT CGESL OR CGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN CGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CAXPY,CSCAL,ICAMAX
C     FORTRAN ABS,AIMAG,REAL
C
C     INTERNAL VARIABLES
C
      COMPLEX T
      INTEGER ICAMAX,J,K,KP1,L,NM1
C
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = ICAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (CABS1(A(L,K)) .EQ. 0.0E0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -(1.0E0,0.0E0)/A(K,K)
            CALL CSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL CAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (CABS1(A(N,N)) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE CPOFA(A,LDA,N,INFO)
      INTEGER LDA,N,INFO
      COMPLEX A(LDA,*)
C
C     CPOFA FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX.
C
C     CPOFA IS USUALLY CALLED BY CPOCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR CPOCO) = (1 + 18/N)*(TIME FOR CPOFA) .
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE HERMITIAN MATRIX TO BE FACTORED.  ONLY THE
C                DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A =
C                CTRANS(R)*R WHERE  CTRANS(R)  IS THE CONJUGATE
C                TRANSPOSE.  THE STRICT LOWER TRIANGLE IS UNALTERED.
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS CDOTC
C     FORTRAN AIMAG,CMPLX,CONJG,REAL,SQRT
C
C     INTERNAL VARIABLES
C
      COMPLEX CDOTC,T
      REAL S
      INTEGER J,JM1,K
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - CDOTC(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + REAL(T*CONJG(T))
   10       CONTINUE
   20       CONTINUE
            S = REAL(A(J,J)) - S
C     ......EXIT
            IF (S .LE. 0.0E0 .OR. AIMAG(A(J,J)) .NE. 0.0E0) GO TO 40
            A(J,J) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
      SUBROUTINE CGTSL(N,C,D,E,B,INFO)
      INTEGER N,INFO
      COMPLEX C(*),D(*),E(*),B(*)
C
C     CGTSL GIVEN A GENERAL TRIDIAGONAL MATRIX AND A RIGHT HAND
C     SIDE WILL FIND THE SOLUTION.
C
C     ON ENTRY
C
C        N       INTEGER
C                IS THE ORDER OF THE TRIDIAGONAL MATRIX.
C
C        C       COMPLEX(N)
C                IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX.
C                C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL.
C                ON OUTPUT C IS DESTROYED.
C
C        D       COMPLEX(N)
C                IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
C                ON OUTPUT D IS DESTROYED.
C
C        E       COMPLEX(N)
C                IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX.
C                E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL.
C                ON OUTPUT E IS DESTROYED.
C
C        B       COMPLEX(N)
C                IS THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       IS THE SOLUTION VECTOR.
C
C        INFO    INTEGER
C                = 0 NORMAL VALUE.
C                = K IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES
C                    EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN
C                    THIS IS DETECTED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
C
C     NO EXTERNALS
C     FORTRAN ABS,AIMAG,REAL
C
C     INTERNAL VARIABLES
C
      INTEGER K,KB,KP1,NM1,NM2
      COMPLEX T
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C     BEGIN BLOCK PERMITTING ...EXITS TO 100
C
         INFO = 0
         C(1) = D(1)
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 40
            D(1) = E(1)
            E(1) = (0.0E0,0.0E0)
            E(N) = (0.0E0,0.0E0)
C
            DO 30 K = 1, NM1
               KP1 = K + 1
C
C              FIND THE LARGEST OF THE TWO ROWS
C
               IF (CABS1(C(KP1)) .LT. CABS1(C(K))) GO TO 10
C
C                 INTERCHANGE ROW
C
                  T = C(KP1)
                  C(KP1) = C(K)
                  C(K) = T
                  T = D(KP1)
                  D(KP1) = D(K)
                  D(K) = T
                  T = E(KP1)
                  E(KP1) = E(K)
                  E(K) = T
                  T = B(KP1)
                  B(KP1) = B(K)
                  B(K) = T
   10          CONTINUE
C
C              ZERO ELEMENTS
C
               IF (CABS1(C(K)) .NE. 0.0E0) GO TO 20
                  INFO = K
C     ............EXIT
                  GO TO 100
   20          CONTINUE
               T = -C(KP1)/C(K)
               C(KP1) = D(KP1) + T*D(K)
               D(KP1) = E(KP1) + T*E(K)
               E(KP1) = (0.0E0,0.0E0)
               B(KP1) = B(KP1) + T*B(K)
   30       CONTINUE
   40    CONTINUE
         IF (CABS1(C(N)) .NE. 0.0E0) GO TO 50
            INFO = N
         GO TO 90
   50    CONTINUE
C
C           BACK SOLVE
C
            NM2 = N - 2
            B(N) = B(N)/C(N)
            IF (N .EQ. 1) GO TO 80
               B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)
               IF (NM2 .LT. 1) GO TO 70
               DO 60 KB = 1, NM2
                  K = NM2 - KB + 1
                  B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)
   60          CONTINUE
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CPTSL(N,D,E,B)
      INTEGER N
      COMPLEX D(*),E(*),B(*)
C
C     CPTSL GIVEN A POSITIVE DEFINITE TRIDIAGONAL MATRIX AND A RIGHT
C     HAND SIDE WILL FIND THE SOLUTION.
C
C     ON ENTRY
C
C        N        INTEGER
C                 IS THE ORDER OF THE TRIDIAGONAL MATRIX.
C
C        D        COMPLEX(N)
C                 IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX.
C                 ON OUTPUT D IS DESTROYED.
C
C        E        COMPLEX(N)
C                 IS THE OFFDIAGONAL OF THE TRIDIAGONAL MATRIX.
C                 E(1) THROUGH E(N-1) SHOULD CONTAIN THE
C                 OFFDIAGONAL.
C
C        B        COMPLEX(N)
C                 IS THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B        CONTAINS THE SOULTION.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
C
C     NO EXTERNALS
C     FORTRAN CONJG,MOD
C
C     INTERNAL VARIABLES
C
      INTEGER K,KBM1,KE,KF,KP1,NM1,NM1D2
      COMPLEX T1,T2
C
C     CHECK FOR 1 X 1 CASE
C
      IF (N .NE. 1) GO TO 10
         B(1) = B(1)/D(1)
      GO TO 70
   10 CONTINUE
         NM1 = N - 1
         NM1D2 = NM1/2
         IF (N .EQ. 2) GO TO 30
            KBM1 = N - 1
C
C           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF
C           SUPERDIAGONAL
C
            DO 20 K = 1, NM1D2
               T1 = CONJG(E(K))/D(K)
               D(K+1) = D(K+1) - T1*E(K)
               B(K+1) = B(K+1) - T1*B(K)
               T2 = E(KBM1)/D(KBM1+1)
               D(KBM1) = D(KBM1) - T2*CONJG(E(KBM1))
               B(KBM1) = B(KBM1) - T2*B(KBM1+1)
               KBM1 = KBM1 - 1
   20       CONTINUE
   30    CONTINUE
         KP1 = NM1D2 + 1
C
C        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER
C
         IF (MOD(N,2) .NE. 0) GO TO 40
            T1 = CONJG(E(KP1))/D(KP1)
            D(KP1+1) = D(KP1+1) - T1*E(KP1)
            B(KP1+1) = B(KP1+1) - T1*B(KP1)
            KP1 = KP1 + 1
   40    CONTINUE
C
C        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP
C        AND BOTTOM
C
         B(KP1) = B(KP1)/D(KP1)
         IF (N .EQ. 2) GO TO 60
            K = KP1 - 1
            KE = KP1 + NM1D2 - 1
            DO 50 KF = KP1, KE
               B(K) = (B(K) - E(K)*B(K+1))/D(K)
               B(KF+1) = (B(KF+1) - CONJG(E(KF))*B(KF))/D(KF+1)
               K = K - 1
   50       CONTINUE
   60    CONTINUE
         IF (MOD(N,2) .EQ. 0) B(1) = (B(1) - E(1)*B(2))/D(1)
   70 CONTINUE
      RETURN
      END
