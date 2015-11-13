#ifdef OMP
MODULE MI_OMP
    IMPLICIT NONE
    INTEGER :: NUM_THREADS = 8
#endif
#ifndef OMP
MODULE MI
    IMPLICIT NONE
#endif
    INTERFACE UNIRNK
      MODULE PROCEDURE R_UNIRNK, I_UNIRNK
    END INTERFACE UNIRNK
    INTERFACE NEARLESS
      MODULE PROCEDURE R_NEARLESS, I_NEARLESS
    END INTERFACE NEARLESS
    INTERFACE PROB2D
        MODULE PROCEDURE I_PROB2D, R_PROB2D
    END INTERFACE PROB2D
    INTERFACE PROBDEF2D
        MODULE PROCEDURE II_PROBDEF2D, RR_PROBDEF2D, RI_PROBDEF2D, IR_PROBDEF2D
    END INTERFACE PROBDEF2D
    INTERFACE MUTUALINFO
        MODULE PROCEDURE R_MUTUALINFO, I_MUTUALINFO
    END INTERFACE MUTUALINFO
    INTERFACE MUTUALINFO_WEIGHT
        MODULE PROCEDURE R_MUTUALINFO_WEIGHT, I_MUTUALINFO_WEIGHT
    END INTERFACE MUTUALINFO_WEIGHT
    INTERFACE DIGITALIZE
        MODULE PROCEDURE R_DIGITALIZE, I_DIGITALIZE
    END INTERFACE DIGITALIZE
    INTERFACE TO1DIM
        MODULE PROCEDURE R_TO1DIM, I_TO1DIM
    END INTERFACE TO1DIM
    INTERFACE MUTUALINFO_PROB
        MODULE PROCEDURE R_MUTUALINFO_PROB, I_MUTUALINFO_PROB
    END INTERFACE MUTUALINFO_PROB
    INTERFACE MUTUALINFO_OTHER
        MODULE PROCEDURE R_MUTUALINFO_OTHER, I_MUTUALINFO_OTHER, RI_MUTUALINFO_OTHER, IR_MUTUALINFO_OTHER
    END INTERFACE MUTUALINFO_OTHER
    INTERFACE MUTUALINFO_SIMP
        MODULE PROCEDURE I_MUTUALINFO_SIMP, R_MUTUALINFO_SIMP
    END INTERFACE MUTUALINFO_SIMP
    INTERFACE MUTUALINFO_SIMP_PROB
        MODULE PROCEDURE I_MUTUALINFO_SIMP_PROB, R_MUTUALINFO_SIMP_PROB
    END INTERFACE MUTUALINFO_SIMP_PROB

    CONTAINS
    
    SUBROUTINE R_UNIRNK (XVALT, IRNGT, NUNI)
    ! __________________________________________________________
    !   UNIRNK = MERGE-SORT RANKING OF AN ARRAY, WITH REMOVAL OF
    !   DUPLICATE ENTRIES.
    !   THE ROUTINE IS SIMILAR TO PURE MERGE-SORT RANKING, BUT ON
    !   THE LAST PASS, IT DISCARDS INDICES THAT CORRESPOND TO
    !   DUPLICATE ENTRIES.
    !   FOR PERFORMANCE REASONS, THE FIRST 2 PASSES ARE TAKEN
    !   OUT OF THE STANDARD LOOP, AND USE DEDICATED CODING.
    ! __________________________________________________________
    ! __________________________________________________________
          REAL(KIND=8), DIMENSION (:), INTENT (IN) :: XVALT
          INTEGER, DIMENSION (:), INTENT (OUT) :: IRNGT
          INTEGER, INTENT (OUT) :: NUNI
    ! __________________________________________________________
          INTEGER, DIMENSION (SIZE(IRNGT)) :: JWRKT
          INTEGER :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
          INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
          REAL(KIND=8) :: XTST, XVALA, XVALB
    !
    !
          NVAL = MIN (SIZE(XVALT), SIZE(IRNGT))
          NUNI = NVAL
    !
          SELECT CASE (NVAL)
          CASE (:0)
             RETURN
          CASE (1)
             IRNGT (1) = 1
             RETURN
          CASE DEFAULT
             CONTINUE
          END SELECT
    !
    !  FILL-IN THE INDEX ARRAY, CREATING ORDERED COUPLES
    !
          DO IIND = 2, NVAL, 2
             IF (XVALT(IIND-1) < XVALT(IIND)) THEN
                IRNGT (IIND-1) = IIND - 1
                IRNGT (IIND) = IIND
             ELSE
                IRNGT (IIND-1) = IIND
                IRNGT (IIND) = IIND - 1
             END IF
          END DO
          IF (MODULO(NVAL, 2) /= 0) THEN
             IRNGT (NVAL) = NVAL
          END IF
    !
    !  WE WILL NOW HAVE ORDERED SUBSETS A - B - A - B - ...
    !  AND MERGE A AND B COUPLES INTO     C   -   C   - ...
    !
          LMTNA = 2
          LMTNC = 4
    !
    !  FIRST ITERATION. THE LENGTH OF THE ORDERED SUBSETS GOES FROM 2 TO 4
    !
          DO
             IF (NVAL <= 4) EXIT
    !
    !   LOOP ON MERGES OF A AND B INTO C
    !
             DO IWRKD = 0, NVAL - 1, 4
                IF ((IWRKD+4) > NVAL) THEN
                   IF ((IWRKD+2) >= NVAL) EXIT
    !
    !   1 2 3
    !
                   IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) EXIT
    !
    !   1 3 2
    !
                   IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
                      IRNG2 = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNG2
    !
    !   3 1 2
    !
                   ELSE
                      IRNG1 = IRNGT (IWRKD+1)
                      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNG1
                   END IF
                   EXIT
                END IF
    !
    !   1 2 3 4
    !
                IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) CYCLE
    !
    !   1 3 X X
    !
                IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                   IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
    !   1 3 2 4
                      IRNGT (IWRKD+3) = IRNG2
                   ELSE
    !   1 3 4 2
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+4) = IRNG2
                   END IF
    !
    !   3 X X X
    !
                ELSE
                   IRNG1 = IRNGT (IWRKD+1)
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                   IF (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) THEN
                      IRNGT (IWRKD+2) = IRNG1
                      IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
    !   3 1 2 4
                         IRNGT (IWRKD+3) = IRNG2
                      ELSE
    !   3 1 4 2
                         IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                         IRNGT (IWRKD+4) = IRNG2
                      END IF
                   ELSE
    !   3 4 1 2
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+3) = IRNG1
                      IRNGT (IWRKD+4) = IRNG2
                   END IF
                END IF
             END DO
    !
    !  THE CS BECOME AS AND BS
    !
             LMTNA = 4
             EXIT
          END DO
    !
    !  ITERATION LOOP. EACH TIME, THE LENGTH OF THE ORDERED SUBSETS
    !  IS DOUBLED.
    !
          DO
             IF (2*LMTNA >= NVAL) EXIT
             IWRKF = 0
             LMTNC = 2 * LMTNC
    !
    !   LOOP ON MERGES OF A AND B INTO C
    !
             DO
                IWRK = IWRKF
                IWRKD = IWRKF + 1
                JINDA = IWRKF + LMTNA
                IWRKF = IWRKF + LMTNC
                IF (IWRKF >= NVAL) THEN
                   IF (JINDA >= NVAL) EXIT
                   IWRKF = NVAL
                END IF
                IINDA = 1
                IINDB = JINDA + 1
    !
    !  ONE STEPS IN THE C SUBSET, THAT WE CREATE IN THE FINAL RANK ARRAY
    !
    !  MAKE A COPY OF THE RANK ARRAY FOR THE ITERATION
    !
                JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
                XVALA = XVALT (JWRKT(IINDA))
                XVALB = XVALT (IRNGT(IINDB))
    !
                DO
                   IWRK = IWRK + 1
    !
    !  WE STILL HAVE UNPROCESSED VALUES IN BOTH A AND B
    !
                   IF (XVALA > XVALB) THEN
                      IRNGT (IWRK) = IRNGT (IINDB)
                      IINDB = IINDB + 1
                      IF (IINDB > IWRKF) THEN
    !  ONLY A STILL WITH UNPROCESSED VALUES
                         IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                         EXIT
                      END IF
                      XVALB = XVALT (IRNGT(IINDB))
                   ELSE
                      IRNGT (IWRK) = JWRKT (IINDA)
                      IINDA = IINDA + 1
                      IF (IINDA > LMTNA) EXIT! ONLY B STILL WITH UNPROCESSED VALUES
                      XVALA = XVALT (JWRKT(IINDA))
                   END IF
    !
                END DO
             END DO
    !
    !  THE CS BECOME AS AND BS
    !
             LMTNA = 2 * LMTNA
          END DO
    !
    !   LAST MERGE OF A AND B INTO C, WITH REMOVAL OF DUPLICATES.
    !
          IINDA = 1
          IINDB = LMTNA + 1
          NUNI = 0
    !
    !  ONE STEPS IN THE C SUBSET, THAT WE CREATE IN THE FINAL RANK ARRAY
    !
          JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
          IF (IINDB <= NVAL) THEN
            XTST = NEARLESS (MIN(XVALT(JWRKT(1)), XVALT(IRNGT(IINDB))))
          ELSE
            XTST = NEARLESS (XVALT(JWRKT(1)))
          ENDIF
          DO IWRK = 1, NVAL
    !
    !  WE STILL HAVE UNPROCESSED VALUES IN BOTH A AND B
    !
             IF (IINDA <= LMTNA) THEN
                IF (IINDB <= NVAL) THEN
                   IF (XVALT(JWRKT(IINDA)) > XVALT(IRNGT(IINDB))) THEN
                      IRNG = IRNGT (IINDB)
                      IINDB = IINDB + 1
                   ELSE
                      IRNG = JWRKT (IINDA)
                      IINDA = IINDA + 1
                   END IF
                ELSE
    !
    !  ONLY A STILL WITH UNPROCESSED VALUES
    !
                   IRNG = JWRKT (IINDA)
                   IINDA = IINDA + 1
                END IF
             ELSE
    !
    !  ONLY B STILL WITH UNPROCESSED VALUES
    !
                IRNG = IRNGT (IWRK)
             END IF
             IF (XVALT(IRNG) > XTST) THEN
                XTST = XVALT (IRNG)
                NUNI = NUNI + 1
                IRNGT (NUNI) = IRNG
             END IF
    !
          END DO
    !
          RETURN
    !
    END SUBROUTINE R_UNIRNK

    SUBROUTINE I_UNIRNK (XVALT, IRNGT, NUNI)
    ! __________________________________________________________
    !   UNIRNK = MERGE-SORT RANKING OF AN ARRAY, WITH REMOVAL OF
    !   DUPLICATE ENTRIES.
    !   THE ROUTINE IS SIMILAR TO PURE MERGE-SORT RANKING, BUT ON
    !   THE LAST PASS, IT DISCARDS INDICES THAT CORRESPOND TO
    !   DUPLICATE ENTRIES.
    !   FOR PERFORMANCE REASONS, THE FIRST 2 PASSES ARE TAKEN
    !   OUT OF THE STANDARD LOOP, AND USE DEDICATED CODING.
    ! __________________________________________________________
    ! __________________________________________________________
          INTEGER, DIMENSION (:), INTENT (IN) :: XVALT
          INTEGER, DIMENSION (:), INTENT (OUT) :: IRNGT
          INTEGER, INTENT (OUT) :: NUNI
    ! __________________________________________________________
          INTEGER, DIMENSION (SIZE(IRNGT)) :: JWRKT
          INTEGER :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
          INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
          INTEGER :: XTST, XVALA, XVALB
    !
    !
          NVAL = MIN (SIZE(XVALT), SIZE(IRNGT))
          NUNI = NVAL
    !
          SELECT CASE (NVAL)
          CASE (:0)
             RETURN
          CASE (1)
             IRNGT (1) = 1
             RETURN
          CASE DEFAULT
             CONTINUE
          END SELECT
    !
    !  FILL-IN THE INDEX ARRAY, CREATING ORDERED COUPLES
    !
          DO IIND = 2, NVAL, 2
             IF (XVALT(IIND-1) < XVALT(IIND)) THEN
                IRNGT (IIND-1) = IIND - 1
                IRNGT (IIND) = IIND
             ELSE
                IRNGT (IIND-1) = IIND
                IRNGT (IIND) = IIND - 1
             END IF
          END DO
          IF (MODULO(NVAL, 2) /= 0) THEN
             IRNGT (NVAL) = NVAL
          END IF
    !
    !  WE WILL NOW HAVE ORDERED SUBSETS A - B - A - B - ...
    !  AND MERGE A AND B COUPLES INTO     C   -   C   - ...
    !
          LMTNA = 2
          LMTNC = 4
    !
    !  FIRST ITERATION. THE LENGTH OF THE ORDERED SUBSETS GOES FROM 2 TO 4
    !
          DO
             IF (NVAL <= 4) EXIT
    !
    !   LOOP ON MERGES OF A AND B INTO C
    !
             DO IWRKD = 0, NVAL - 1, 4
                IF ((IWRKD+4) > NVAL) THEN
                   IF ((IWRKD+2) >= NVAL) EXIT
    !
    !   1 2 3
    !
                   IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) EXIT
    !
    !   1 3 2
    !
                   IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
                      IRNG2 = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNG2
    !
    !   3 1 2
    !
                   ELSE
                      IRNG1 = IRNGT (IWRKD+1)
                      IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                      IRNGT (IWRKD+2) = IRNG1
                   END IF
                   EXIT
                END IF
    !
    !   1 2 3 4
    !
                IF (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) CYCLE
    !
    !   1 3 X X
    !
                IF (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) THEN
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                   IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
    !   1 3 2 4
                      IRNGT (IWRKD+3) = IRNG2
                   ELSE
    !   1 3 4 2
                      IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+4) = IRNG2
                   END IF
    !
    !   3 X X X
    !
                ELSE
                   IRNG1 = IRNGT (IWRKD+1)
                   IRNG2 = IRNGT (IWRKD+2)
                   IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                   IF (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) THEN
                      IRNGT (IWRKD+2) = IRNG1
                      IF (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) THEN
    !   3 1 2 4
                         IRNGT (IWRKD+3) = IRNG2
                      ELSE
    !   3 1 4 2
                         IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                         IRNGT (IWRKD+4) = IRNG2
                      END IF
                   ELSE
    !   3 4 1 2
                      IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                      IRNGT (IWRKD+3) = IRNG1
                      IRNGT (IWRKD+4) = IRNG2
                   END IF
                END IF
             END DO
    !
    !  THE CS BECOME AS AND BS
    !
             LMTNA = 4
             EXIT
          END DO
    !
    !  ITERATION LOOP. EACH TIME, THE LENGTH OF THE ORDERED SUBSETS
    !  IS DOUBLED.
    !
          DO
             IF (2*LMTNA >= NVAL) EXIT
             IWRKF = 0
             LMTNC = 2 * LMTNC
    !
    !   LOOP ON MERGES OF A AND B INTO C
    !
             DO
                IWRK = IWRKF
                IWRKD = IWRKF + 1
                JINDA = IWRKF + LMTNA
                IWRKF = IWRKF + LMTNC
                IF (IWRKF >= NVAL) THEN
                   IF (JINDA >= NVAL) EXIT
                   IWRKF = NVAL
                END IF
                IINDA = 1
                IINDB = JINDA + 1
    !
    !  ONE STEPS IN THE C SUBSET, THAT WE CREATE IN THE FINAL RANK ARRAY
    !
    !  MAKE A COPY OF THE RANK ARRAY FOR THE ITERATION
    !
                JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
                XVALA = XVALT (JWRKT(IINDA))
                XVALB = XVALT (IRNGT(IINDB))
    !
                DO
                   IWRK = IWRK + 1
    !
    !  WE STILL HAVE UNPROCESSED VALUES IN BOTH A AND B
    !
                   IF (XVALA > XVALB) THEN
                      IRNGT (IWRK) = IRNGT (IINDB)
                      IINDB = IINDB + 1
                      IF (IINDB > IWRKF) THEN
    !  ONLY A STILL WITH UNPROCESSED VALUES
                         IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                         EXIT
                      END IF
                      XVALB = XVALT (IRNGT(IINDB))
                   ELSE
                      IRNGT (IWRK) = JWRKT (IINDA)
                      IINDA = IINDA + 1
                      IF (IINDA > LMTNA) EXIT! ONLY B STILL WITH UNPROCESSED VALUES
                      XVALA = XVALT (JWRKT(IINDA))
                   END IF
    !
                END DO
             END DO
    !
    !  THE CS BECOME AS AND BS
    !
             LMTNA = 2 * LMTNA
          END DO
    !
    !   LAST MERGE OF A AND B INTO C, WITH REMOVAL OF DUPLICATES.
    !
          IINDA = 1
          IINDB = LMTNA + 1
          NUNI = 0
    !
    !  ONE STEPS IN THE C SUBSET, THAT WE CREATE IN THE FINAL RANK ARRAY
    !
          JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
          IF (IINDB <= NVAL) THEN
            XTST = NEARLESS (MIN(XVALT(JWRKT(1)), XVALT(IRNGT(IINDB))))
          ELSE
            XTST = NEARLESS (XVALT(JWRKT(1)))
          ENDIF
          DO IWRK = 1, NVAL
    !
    !  WE STILL HAVE UNPROCESSED VALUES IN BOTH A AND B
    !
             IF (IINDA <= LMTNA) THEN
                IF (IINDB <= NVAL) THEN
                   IF (XVALT(JWRKT(IINDA)) > XVALT(IRNGT(IINDB))) THEN
                      IRNG = IRNGT (IINDB)
                      IINDB = IINDB + 1
                   ELSE
                      IRNG = JWRKT (IINDA)
                      IINDA = IINDA + 1
                   END IF
                ELSE
    !
    !  ONLY A STILL WITH UNPROCESSED VALUES
    !
                   IRNG = JWRKT (IINDA)
                   IINDA = IINDA + 1
                END IF
             ELSE
    !
    !  ONLY B STILL WITH UNPROCESSED VALUES
    !
                IRNG = IRNGT (IWRK)
             END IF
             IF (XVALT(IRNG) > XTST) THEN
                XTST = XVALT (IRNG)
                NUNI = NUNI + 1
                IRNGT (NUNI) = IRNG
             END IF
    !
          END DO
    !
          RETURN
    !
    END SUBROUTINE I_UNIRNK
    
    FUNCTION R_NEARLESS (XVAL) RESULT (R_NL)
    !  NEAREST VALUE LESS THAN GIVEN VALUE
    ! __________________________________________________________
          REAL(KIND=8), INTENT (IN) :: XVAL
          REAL(KIND=8) :: R_NL
    ! __________________________________________________________
          R_NL = NEAREST (XVAL, -1.0)
          RETURN
    !
    END FUNCTION R_NEARLESS
    FUNCTION I_NEARLESS (XVAL) RESULT (I_NL)
    !  NEAREST VALUE LESS THAN GIVEN VALUE
    ! __________________________________________________________
          INTEGER, INTENT (IN) :: XVAL
          INTEGER :: I_NL
    ! __________________________________________________________
          I_NL = XVAL - 1
          RETURN
    !
    END FUNCTION I_NEARLESS

    SUBROUTINE I_PROB2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)
  
    IMPLICIT NONE
    
    INTEGER, DIMENSION(0:N-1)                 :: D_X, D_Y
    INTEGER                                   :: N
    INTEGER                                   :: NBINS_X, NBINS_Y
    REAL(KIND=8), DIMENSION(0:NBINS_X-1,0:NBINS_Y-1)  :: PROB
    REAL(KIND=8), DIMENSION(0:NBINS_X)                :: BINS_X
    REAL(KIND=8), DIMENSION(0:NBINS_Y)                :: BINS_Y
    
    INTEGER               :: I, J, X, Y
    REAL(KIND=8)                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL(KIND=8)                  :: BIN_DELTA_X, BIN_DELTA_Y
    REAL(KIND=8)                  :: P
    
    P = 1.0 / N
        
    MIN_D_X = MINVAL(D_X)
    MAX_D_X = MAXVAL(D_X)
    MIN_D_Y = MINVAL(D_Y)
    MAX_D_Y = MAXVAL(D_Y)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    
    BINS_X    = (/ ( MIN_D_X + I * BIN_DELTA_X, I=0,NBINS_X) /)
    BINS_Y    = (/ ( MIN_D_Y + I * BIN_DELTA_Y, I=0,NBINS_Y) /)
    FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    
    DO I=0,N-1
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X),NBINS_X-1)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y),NBINS_Y-1)
        PROB(X,Y) = PROB(X,Y) + P
    END DO
        
    RETURN
    
    END SUBROUTINE I_PROB2D
    
    SUBROUTINE R_PROB2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    REAL(KIND=8)   , DIMENSION(0:N-1)                 :: D_X, D_Y
    INTEGER                                   :: N
    INTEGER                                   :: NBINS_X, NBINS_Y
    REAL(KIND=8), DIMENSION(0:NBINS_X-1,0:NBINS_Y-1)  :: PROB
    REAL(KIND=8), DIMENSION(0:NBINS_X)                :: BINS_X
    REAL(KIND=8), DIMENSION(0:NBINS_Y)                :: BINS_Y
    
    INTEGER               :: I, J, X, Y
    REAL(KIND=8)                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL(KIND=8)                  :: BIN_DELTA_X, BIN_DELTA_Y
    REAL(KIND=8)                  :: P
    
    P = 1.0 / N
    
    MIN_D_X = MINVAL(D_X)
    MAX_D_X = MAXVAL(D_X)
    MIN_D_Y = MINVAL(D_Y)
    MAX_D_Y = MAXVAL(D_Y)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    
    BINS_X    = (/ ( MIN_D_X + I * BIN_DELTA_X, I=0,NBINS_X) /)
    BINS_Y    = (/ ( MIN_D_Y + I * BIN_DELTA_Y, I=0,NBINS_Y) /)
    FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    
    DO I=0,N-1
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X),NBINS_X-1)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y),NBINS_Y-1)
        PROB(X,Y) = PROB(X,Y) + P
    END DO
    
    RETURN
        
    END SUBROUTINE R_PROB2D

    SUBROUTINE II_PROBDEF2D(D_X,D_Y,BINS_X,BINS_Y,N,NBINS_X,NBINS_Y,PROB)

    IMPLICIT NONE
    
    INTEGER, DIMENSION(0:N-1), INTENT(IN)                  :: D_X, D_Y
    INTEGER, INTENT(IN)                                    :: N
    INTEGER, INTENT(IN)                                    :: NBINS_X, NBINS_Y
    REAL(KIND=8), DIMENSION(0:NBINS_X-1,0:NBINS_Y-1) ,INTENT(OUT)  :: PROB
    REAL(KIND=8), DIMENSION(0:NBINS_X), INTENT(IN)                 :: BINS_X
    REAL(KIND=8), DIMENSION(0:NBINS_Y), INTENT(IN)                 :: BINS_Y
    
    INTEGER               :: I, J, K, X, Y
    REAL(KIND=8)                  :: P 
    
    P = 1.0 / N
    
    FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    DO I=0,N-1
        X = NBINS_X-1
        Y = NBINS_Y-1
        DO K=1,NBINS_X-1
            IF ( BINS_X(K) > D_X(I) )  THEN
                X = K - 1
                EXIT
            END IF
        END DO
        DO K=1,NBINS_Y-1
            IF ( BINS_Y(K) > D_Y(I) ) THEN
                Y = K - 1
                EXIT
            END IF
        END DO
        PROB(X,Y) = PROB(X,Y) + P
    END DO
    
    END SUBROUTINE II_PROBDEF2D
    
    SUBROUTINE RR_PROBDEF2D(D_X,D_Y,BINS_X,BINS_Y,N,NBINS_X,NBINS_Y,PROB)

    IMPLICIT NONE
    
    REAL(KIND=8), DIMENSION(0:N-1), INTENT(IN)                     :: D_X, D_Y
    INTEGER, INTENT(IN)                                    :: N
    INTEGER, INTENT(IN)                                    :: NBINS_X, NBINS_Y
    REAL(KIND=8), DIMENSION(0:NBINS_X-1,0:NBINS_Y-1) ,INTENT(OUT)  :: PROB
    REAL(KIND=8), DIMENSION(0:NBINS_X), INTENT(IN)                 :: BINS_X
    REAL(KIND=8), DIMENSION(0:NBINS_Y), INTENT(IN)                 :: BINS_Y
    
    INTEGER               :: I, J, K, X, Y
    REAL(KIND=8)                  :: P 
    
    P = 1.0 / N
    
    FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    DO I=0,N-1
        X = NBINS_X-1
        Y = NBINS_Y-1
        DO K=1,NBINS_X-1
            IF ( BINS_X(K) > D_X(I) )  THEN
                X = K - 1
                EXIT
            END IF
        END DO
        DO K=1,NBINS_Y-1
            IF ( BINS_Y(K) > D_Y(I) ) THEN
                Y = K - 1
                EXIT
            END IF
        END DO
        PROB(X,Y) = PROB(X,Y) + P
    END DO
    
    END SUBROUTINE RR_PROBDEF2D
    
    SUBROUTINE RI_PROBDEF2D(D_X,D_Y,BINS_X,BINS_Y,N,NBINS_X,NBINS_Y,PROB)

    IMPLICIT NONE
    
    REAL(KIND=8),    DIMENSION(0:N-1), INTENT(IN)                  :: D_X
    INTEGER, DIMENSION(0:N-1), INTENT(IN)                  :: D_Y
    INTEGER, INTENT(IN)                                    :: N
    INTEGER, INTENT(IN)                                    :: NBINS_X, NBINS_Y
    REAL(KIND=8), DIMENSION(0:NBINS_X-1,0:NBINS_Y-1) ,INTENT(OUT)  :: PROB
    REAL(KIND=8), DIMENSION(0:NBINS_X), INTENT(IN)                 :: BINS_X
    REAL(KIND=8), DIMENSION(0:NBINS_Y), INTENT(IN)                 :: BINS_Y
    
    INTEGER               :: I, J, K, X, Y
    REAL(KIND=8)                  :: P 
    
    P = 1.0 / N
    
    FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    DO I=0,N-1
        X = NBINS_X-1
        Y = NBINS_Y-1
        DO K=1,NBINS_X-1
            IF ( BINS_X(K) > D_X(I) )  THEN
                X = K - 1
                EXIT
            END IF
        END DO
        DO K=1,NBINS_Y-1
            IF ( BINS_Y(K) > D_Y(I) ) THEN
                Y = K - 1
                EXIT
            END IF
        END DO
        PROB(X,Y) = PROB(X,Y) + P
    END DO
    
    END SUBROUTINE RI_PROBDEF2D
    
    SUBROUTINE IR_PROBDEF2D(D_X,D_Y,BINS_X,BINS_Y,N,NBINS_X,NBINS_Y,PROB)

    IMPLICIT NONE
    
    INTEGER, DIMENSION(0:N-1), INTENT(IN)                  :: D_X
    REAL(KIND=8),    DIMENSION(0:N-1), INTENT(IN)                  :: D_Y
    INTEGER, INTENT(IN)                                    :: N
    INTEGER, INTENT(IN)                                    :: NBINS_X, NBINS_Y
    REAL(KIND=8), DIMENSION(0:NBINS_X-1,0:NBINS_Y-1) ,INTENT(OUT)  :: PROB
    REAL(KIND=8), DIMENSION(0:NBINS_X), INTENT(IN)                 :: BINS_X
    REAL(KIND=8), DIMENSION(0:NBINS_Y), INTENT(IN)                 :: BINS_Y
    
    INTEGER               :: I, J, K, X, Y
    REAL(KIND=8)                  :: P 
    
    P = 1.0 / N
    
    FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    DO I=0,N-1
        X = NBINS_X-1
        Y = NBINS_Y-1
        DO K=1,NBINS_X-1
            IF ( BINS_X(K) > D_X(I) )  THEN
                X = K - 1
                EXIT
            END IF
        END DO
        DO K=1,NBINS_Y-1
            IF ( BINS_Y(K) > D_Y(I) ) THEN
                Y = K - 1
                EXIT
            END IF
        END DO
        PROB(X,Y) = PROB(X,Y) + P
    END DO
    
    END SUBROUTINE IR_PROBDEF2D
    
    SUBROUTINE R_MUTUALINFO_PROB(D,E1,BINS,NFRAMES,NREP,NBINS,M,EJ,PJ)

    INTEGER, INTENT(IN)                                :: NBINS
    INTEGER, INTENT(IN)                                :: NFRAMES
    INTEGER, INTENT(IN)                                :: NREP
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1)  :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)              :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)               :: BINS
    
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)             :: P_TEMP
    INTEGER                                          :: I, J, K, L
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                     :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: EJ
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NBINS-1,0:NBINS-1,0:NREP-1,0:NREP-1) :: PJ
    
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE R_MUTUALINFO_PROB"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D,E1,EJ,PJ,M,NFRAMES,NREP,NBINS,BINS)
    DO I = 0,NREP-2
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO 
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(K,L,I,J) = P_TEMP(K,L)
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE R_MUTUALINFO_PROB

    SUBROUTINE I_MUTUALINFO_PROB(D,E1,BINS,NFRAMES,NREP,NBINS,M,EJ,PJ)

    INTEGER, INTENT(IN)                                  :: NBINS
    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)                 :: BINS
    
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)        :: P_TEMP
    INTEGER                                             :: I, J, K, L
    INTEGER                                             :: X, Y
    REAL(KIND=8)                                        :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: EJ
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NBINS-1,0:NBINS-1,0:NREP-1,0:NREP-1) :: PJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE I_MUTUALINFO_PROB"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D,E1,EJ,PJ,M,NFRAMES,NREP,NBINS,BINS)
    DO I = 0,NREP-2
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO 
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(K,L,I,J) = P_TEMP(K,L)
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE I_MUTUALINFO_PROB
    
    SUBROUTINE R_MUTUALINFO_OTHER_PROB(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ,PJ)

    INTEGER, INTENT(IN)                                   :: NBINS1,NBINS2
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)         :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)         :: BINS2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)        :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)        :: E2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
        
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)           :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                     :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: EJ
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NBINS2-1,0:NBINS1-1,0:NREP2-1,0:NREP1-1) :: PJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE R_MUTUALINFO_OTHER_PROB"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(K,J,I,J) = P_TEMP(K,J)
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I)
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
     PJ = PJ / NFRAMES
     
    END SUBROUTINE R_MUTUALINFO_OTHER_PROB

    SUBROUTINE I_MUTUALINFO_OTHER_PROB(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ,PJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1) :: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1) :: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)        :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)        :: E2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)         :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)         :: BINS2

    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)      :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                     :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: EJ
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NBINS2-1,0:NBINS1-1,0:NREP2-1,0:NREP1-1) :: PJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE I_MUTUALINFO_OTHER_PROB"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(K,J,I,J) = P_TEMP(K,J)
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I)
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
     PJ = PJ / NFRAMES
    
    END SUBROUTINE I_MUTUALINFO_OTHER_PROB

    SUBROUTINE RI_MUTUALINFO_OTHER_PROB(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ,PJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1) :: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1) :: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)        :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)        :: E2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)         :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)         :: BINS2

    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)      :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                     :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: EJ
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NBINS2-1,0:NBINS1-1,0:NREP2-1,0:NREP1-1) :: PJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE RI_MUTUALINFO_OTHER_PROB"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(K,J,I,J) = P_TEMP(K,J)
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I)
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
     PJ = PJ / NFRAMES
    
    END SUBROUTINE RI_MUTUALINFO_OTHER_PROB
    
    SUBROUTINE IR_MUTUALINFO_OTHER_PROB(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ,PJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1) :: D2
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1) :: D1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)        :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)        :: E2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)         :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)         :: BINS2

    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                :: D_TEMP2
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)      :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                     :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: EJ
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NBINS2-1,0:NBINS1-1,0:NREP2-1,0:NREP1-1) :: PJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE IR_MUTUALINFO_OTHER_PROB"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(K,J,I,J) = P_TEMP(K,J)
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I)
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
     PJ = PJ / NFRAMES
    
    END SUBROUTINE IR_MUTUALINFO_OTHER_PROB
    
    SUBROUTINE I_MUTUALINFO_SIMP(D,E1,BINS,NFRAMES,NREP,NBINS,NREP1,M,EJ)

    INTEGER, INTENT(IN)                                  :: NREP1
    INTEGER, INTENT(IN)                                  :: NBINS
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN)                                  :: NFRAMES
    
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)                 :: BINS
    
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)             :: P_TEMP
    INTEGER                             :: I, K, L
    INTEGER                             :: X, Y
    REAL(KIND=8)                        :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1)  :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1)  :: EJ

    
    P = 1.0/FLOAT(NFRAMES)
    D_TEMP1(0:NFRAMES-1) = (/ (D(K,NREP1), K=0,NFRAMES-1) /)
    EJ(NREP1) = E1(NREP1)
    M(NREP1)  = E1(NREP1)
    
    WRITE (*,'(A)') "SUBROUTINE I_MUTUALINFO_SIMP"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,L,K) &
    !$OMP SHARED(D,E1,EJ,NREP1,NFRAMES,NBINS,BINS) 
    DO I = 0,NREP-1
        IF (I /= NREP1) THEN
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
            EJ(I) = 0.0
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO 
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(I) = EJ(I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(I) = E1(I) + E1(NREP1) - EJ(I)
        END IF
    END DO
    !$OMP END PARALLEL DO

    END SUBROUTINE I_MUTUALINFO_SIMP
    
    SUBROUTINE R_MUTUALINFO_SIMP(D,E1,BINS,NREP1,NFRAMES,NREP,NBINS,M,EJ)
    
    INTEGER, INTENT(IN)                               :: NREP1 
    INTEGER, INTENT(IN)                               :: NBINS
    INTEGER, INTENT(IN)                               :: NREP
    INTEGER, INTENT(IN)                               :: NFRAMES
    
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)             :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)              :: BINS

    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)             :: P_TEMP
    INTEGER                             :: I, K, L
    INTEGER                             :: X, Y
    REAL(KIND=8)                        :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1)   :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1)   :: EJ

    
    P = 1.0/FLOAT(NFRAMES)
    D_TEMP1(0:NFRAMES-1) = (/ (D(K,NREP1), K=0,NFRAMES-1) /)
    EJ(NREP1) = E1(NREP1)
    M(NREP1)  = E1(NREP1)
    
    WRITE (*,'(A)') "SUBROUTINE R_MUTUALINFO_SIMP"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,L,K) &
    !$OMP SHARED(D,E1,EJ,NREP1,NFRAMES,NBINS,BINS) 
    DO I = 0,NREP-1
        IF (I /= NREP1) THEN
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
            EJ(I) = 0.0
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO 
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(I) = EJ(I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(I) = E1(I) + E1(NREP1) - EJ(I)
        END IF
    END DO
    !$OMP END PARALLEL DO
        
    END SUBROUTINE R_MUTUALINFO_SIMP

    SUBROUTINE I_MUTUALINFO_SIMP_PROB(D,E1,BINS,NFRAMES,NREP,NBINS,NREP1,M,EJ,PJ)

    INTEGER, INTENT(IN)                                  :: NREP1
    INTEGER, INTENT(IN)                                  :: NBINS
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN)                                  :: NFRAMES
    
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)                 :: BINS
    
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)             :: P_TEMP
    INTEGER                             :: I, K, L
    INTEGER                             :: X, Y
    REAL(KIND=8)                        :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1)                     :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1)                     :: EJ
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NBINS-1,0:NBINS-1,0:NREP-1) :: PJ
    
    P = 1.0/FLOAT(NFRAMES)
    D_TEMP1(0:NFRAMES-1) = (/ (D(K,NREP1), K=0,NFRAMES-1) /)
    EJ(NREP1) = E1(NREP1)
    M(NREP1)  = E1(NREP1)
    
    WRITE (*,'(A)') "SUBROUTINE I_MUTUALINFO_SIMP_PROB"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,L,K) &
    !$OMP SHARED(D,E1,EJ,NREP1,NFRAMES,NBINS,BINS) 
    DO I = 0,NREP-1
        IF (I /= NREP1) THEN
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
            EJ(I) = 0.0
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO 
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    PJ(L,K,I) = P_TEMP(L,K)
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(I) = EJ(I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(I) = E1(I) + E1(NREP1) - EJ(I)
        END IF
    END DO
    !$OMP END PARALLEL DO
    PJ = PJ/NFRAMES

    END SUBROUTINE I_MUTUALINFO_SIMP_PROB
    
    SUBROUTINE R_MUTUALINFO_SIMP_PROB(D,E1,BINS,NREP1,NFRAMES,NREP,NBINS,M,EJ,PJ)
    
    INTEGER, INTENT(IN)                               :: NREP1 
    INTEGER, INTENT(IN)                               :: NBINS
    INTEGER, INTENT(IN)                               :: NREP
    INTEGER, INTENT(IN)                               :: NFRAMES
    
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)             :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)              :: BINS

    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)             :: P_TEMP
    INTEGER                             :: I, K, L
    INTEGER                             :: X, Y
    REAL(KIND=8)                        :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1)                     :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1)                     :: EJ
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NBINS-1,0:NBINS-1,0:NREP-1) :: PJ
    
    P = 1.0/FLOAT(NFRAMES)
    D_TEMP1(0:NFRAMES-1) = (/ (D(K,NREP1), K=0,NFRAMES-1) /)
    EJ(NREP1) = E1(NREP1)
    M(NREP1)  = E1(NREP1)
    
    WRITE (*,'(A)') "SUBROUTINE R_MUTUALINFO_SIMP_PROB"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,L,K) &
    !$OMP SHARED(D,E1,EJ,NREP1,NFRAMES,NBINS,BINS) 
    DO I = 0,NREP-1
        IF (I /= NREP1) THEN
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
            EJ(I) = 0.0
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO 
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    PJ(L,K,I) = P_TEMP(L,K)
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(I) = EJ(I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(I) = E1(I) + E1(NREP1) - EJ(I)
        END IF
    END DO
    !$OMP END PARALLEL DO
    PJ = PJ/NFRAMES
        
    END SUBROUTINE R_MUTUALINFO_SIMP_PROB

    SUBROUTINE R_MUTUALINFO(D,E1,BINS,NFRAMES,NREP,NBINS,M,EJ)

    INTEGER, INTENT(IN)                                :: NBINS
    INTEGER, INTENT(IN)                                :: NFRAMES
    INTEGER, INTENT(IN)                                :: NREP
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1)  :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)              :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)               :: BINS
    
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                       :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                       :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)               :: P_TEMP
    INTEGER                                            :: I, J, K, L
    INTEGER                                            :: X, Y
    REAL(KIND=8)                                       :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: EJ
    
    P = 1.0 / FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE R_MUTUALINFO"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D,E1,EJ,M,NFRAMES,P,NBINS,BINS)
    DO I = 0,NREP-2
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(X,Y) = P_TEMP(X,Y) + P
            END DO 
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE R_MUTUALINFO

    SUBROUTINE I_MUTUALINFO(D,E1,BINS,NFRAMES,NREP,NBINS,M,EJ)

    INTEGER, INTENT(IN)                                  :: NBINS
    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)                 :: BINS
    
    INTEGER, DIMENSION(0:NFRAMES-1)                      :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                      :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)         :: P_TEMP
    INTEGER                                              :: I, J, K, L
    INTEGER                                              :: X, Y
    REAL(KIND=8)                                         :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: EJ
    
    P  = 1.0 / FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE I_MUTUALINFO"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D,E1,EJ,M,NFRAMES,NBINS,BINS)
    DO I = 0,NREP-2
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(X,Y) = P_TEMP(X,Y) + P
            END DO 
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 90) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE I_MUTUALINFO

    SUBROUTINE R_MUTUALINFO_WEIGHT(D,E1,W,BINS,NFRAMES,NREP,NBINS,M,EJ)

    INTEGER, INTENT(IN)                                :: NBINS
    INTEGER, INTENT(IN)                                :: NFRAMES
    INTEGER, INTENT(IN)                                :: NREP
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1)  :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)              :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1)           :: W
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)               :: BINS
    
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                       :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                       :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)               :: P_TEMP
    
    REAL(KIND=8)                                               :: P, SUM_W, L_SUM_W
    INTEGER                                            :: I, J, K
    INTEGER                                            :: L, X, Y
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)    :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)    :: EJ
    
    P = 1.0 / FLOAT(NFRAMES)
    SUM_W = SUM(W)
    L_SUM_W = LOG(SUM_W)
    
    WRITE (*,'(A)') "SUBROUTINE R_MUTUALINFO_WEIGHT"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D,W,SUM_W,E1,EJ,M,NFRAMES,NBINS,BINS)
    DO I = 0,NREP-2
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(X,Y) = P_TEMP(X,Y) + P * W(K)
            END DO 
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(J,I) = EJ(J,I)/SUM_W + L_SUM_W
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE R_MUTUALINFO_WEIGHT

    SUBROUTINE I_MUTUALINFO_WEIGHT(D,E1,W,BINS,NFRAMES,NREP,NBINS,M,EJ)

    INTEGER, INTENT(IN)                                 :: NBINS
    INTEGER, INTENT(IN)                                 :: NFRAMES
    INTEGER, INTENT(IN)                                 :: NREP
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1):: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)               :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1)            :: W
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS)                :: BINS
    
    INTEGER, DIMENSION(0:NFRAMES-1)                    :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                    :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS-1,0:NBINS-1)               :: P_TEMP
    
    REAL(KIND=8)                                               :: P, SUM_W, L_SUM_W
    INTEGER                                            :: I, J, K
    INTEGER                                            :: L, X, Y
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)    :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)    :: EJ
    
    P = 1.0 / FLOAT(NFRAMES)
    SUM_W = SUM(W)
    L_SUM_W = LOG(SUM_W)
    
    WRITE (*,'(A)') "SUBROUTINE I_MUTUALINFO_WEIGHT"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D,W,SUM_W,E1,EJ,M,NFRAMES,NBINS,BINS)
    DO I = 0,NREP-2
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS-1, L=0:NBINS-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS-1
                Y = NBINS-1
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS-1
                    IF ( BINS(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(X,Y) = P_TEMP(X,Y) + P * W(K)
            END DO
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(J,I) = EJ(J,I)/SUM_W + L_SUM_W
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE I_MUTUALINFO_WEIGHT
        
    SUBROUTINE R_MUTUALINFO_OTHER(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)

    INTEGER, INTENT(IN)                                   :: NBINS1,NBINS2
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
        
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)           :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                     :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: EJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE R_MUTUALINFO_OTHER"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
     
    END SUBROUTINE R_MUTUALINFO_OTHER

    SUBROUTINE I_MUTUALINFO_OTHER(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1) :: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1) :: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2
    
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    ! DIMENSION OF PROB ARRAY ARE INVERTED BECAUSE THIS ROUTINE IS DERIVATIVE
    ! OF A ROUTINE THAT PASS PROB AS AN OUTPUT TO PYTHON
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)              :: P_TEMP  
    INTEGER                                             :: I, J, K, L
    INTEGER                                             :: X, Y
    REAL(KIND=8)                                        :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: EJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE I_MUTUALINFO_OTHER"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE I_MUTUALINFO_OTHER
    
    SUBROUTINE IR_MUTUALINFO_OTHER(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1) :: D1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2

    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                        :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)              :: P_TEMP
    INTEGER                                             :: I, J, K, L
    INTEGER                                             :: X, Y
    REAL(KIND=8)                                        :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: EJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE IR_MUTUALINFO_OTHER"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE IR_MUTUALINFO_OTHER

    SUBROUTINE RI_MUTUALINFO_OTHER(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1) :: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2

    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)      :: P_TEMP
    INTEGER                                             :: I, J, K, L
    INTEGER                                             :: X, Y
    REAL(KIND=8)                                        :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)   :: EJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE RI_MUTUALINFO_OTHER"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE RI_MUTUALINFO_OTHER
    
    SUBROUTINE R_MUTUALINFO_OTHER_WEIGTH(D1,D2,E1,E2,W,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)

    INTEGER, INTENT(IN)                                   :: NBINS1,NBINS2
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1)              :: W
        
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)           :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                     :: P, SUM_W, L_SUM_W
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)  :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)  :: EJ
    
    P = 1.0 / FLOAT(NFRAMES)
    SUM_W = SUM(W)
    L_SUM_W = LOG(SUM_W)
    
    WRITE (*,'(A)') "SUBROUTINE R_MUTUALINFO_OTHER_WEIGTH"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif    
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2,P,SUM_W,W)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P * W(K)
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(J,I) = EJ(J,I)/SUM_W + L_SUM_W
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO  
     
    END SUBROUTINE R_MUTUALINFO_OTHER_WEIGTH
    
    SUBROUTINE I_MUTUALINFO_OTHER_WEIGTH(D1,D2,E1,E2,W,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)

    INTEGER, INTENT(IN)                                   :: NBINS1,NBINS2
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1)              :: W
        
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)      :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                     :: P, SUM_W, L_SUM_W
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)  :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)  :: EJ
    
    P = 1.0 / FLOAT(NFRAMES)
    SUM_W = SUM(W)
    L_SUM_W = LOG(SUM_W)
    
    WRITE (*,'(A)') "SUBROUTINE I_MUTUALINFO_OTHER_WEIGTH"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif    
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2,P,SUM_W,W)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P * W(K)
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(J,I) = EJ(J,I)/SUM_W + L_SUM_W
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO  
     
    END SUBROUTINE I_MUTUALINFO_OTHER_WEIGTH
    
    SUBROUTINE IR_MUTUALINFO_OTHER_WEIGTH(D1,D2,E1,E2,W,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)

    INTEGER, INTENT(IN)                                   :: NBINS1,NBINS2
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1)              :: W
        
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)           :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                     :: P, SUM_W, L_SUM_W
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)  :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)  :: EJ
    
    P = 1.0 / FLOAT(NFRAMES)
    SUM_W = SUM(W)
    L_SUM_W = LOG(SUM_W)
    
    WRITE (*,'(A)') "SUBROUTINE IR_MUTUALINFO_OTHER_WEIGTH"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif    
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2,P,SUM_W,W)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P * W(K)
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(J,I) = EJ(J,I)/SUM_W + L_SUM_W
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO  
     
    END SUBROUTINE IR_MUTUALINFO_OTHER_WEIGTH
    
    SUBROUTINE RI_MUTUALINFO_OTHER_WEIGTH(D1,D2,E1,E2,W,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)

    INTEGER, INTENT(IN)                                   :: NBINS1,NBINS2
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1)              :: W
        
    REAL(KIND=8), DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL(KIND=8), DIMENSION(0:NBINS2-1,0:NBINS1-1)           :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    INTEGER                                          :: X, Y
    REAL(KIND=8)                                             :: P, SUM_W, L_SUM_W
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)  :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)  :: EJ
    
    P = 1.0 / FLOAT(NFRAMES)
    SUM_W = SUM(W)
    L_SUM_W = LOG(SUM_W)
        
    WRITE (*,'(A)') "SUBROUTINE RI_MUTUALINFO_OTHER_WEIGTH"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif    
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K) &
    !$OMP SHARED(D1,D2,E1,E2,EJ,M,NFRAMES,NBINS1,NBINS2,BINS1,BINS2,P,SUM_W,W)
    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(K,L)=0.0
            DO K=0,NFRAMES-1
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(Y,X) = P_TEMP(Y,X) + P * W(K)
            END DO
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(J,I) = EJ(J,I)/SUM_W + L_SUM_W
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    !$OMP END PARALLEL DO  
     
    END SUBROUTINE RI_MUTUALINFO_OTHER_WEIGTH

    FUNCTION OPT_BIN_TRAJ(D,NFRAMES)
    
    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1)          :: D
    
    REAL(KIND=8), ALLOCATABLE                                    :: OPT_BIN_TRAJ(:)
    INTEGER, DIMENSION(0:NFRAMES-1)                      :: I_BINS

    INTEGER                                              :: NBINS
    INTEGER                                              :: I
    
    CALL UNIRNK(D,I_BINS,NBINS)
    
    IF (ALLOCATED(OPT_BIN_TRAJ)) DEALLOCATE(OPT_BIN_TRAJ)
    ALLOCATE(OPT_BIN_TRAJ(0:NBINS))
    
    OPT_BIN_TRAJ(0)         = D(I_BINS(0))
    OPT_BIN_TRAJ(1:NBINS-1) = (/ ( ( D(I_BINS(I)) + D(I_BINS(I-1)) ) / 2.0 , I = 1,NBINS-1) /)
    OPT_BIN_TRAJ(NBINS)     = D(I_BINS(NBINS-1))
    
    END FUNCTION OPT_BIN_TRAJ

    SUBROUTINE I_DIGITALIZE(D,BINS,NREP,NDIM,NFRAMES,NBINS,O)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)  :: NREP
    INTEGER, INTENT(IN)  :: NDIM
    INTEGER, INTENT(IN)  :: NFRAMES
    INTEGER, INTENT(IN)  :: NBINS
    
    INTEGER, INTENT(IN)  :: D(0:NFRAMES-1,0:NDIM-1,0:NREP-1)
    REAL(KIND=8), INTENT(IN)     :: BINS(0:NBINS,0:NDIM-1)
    
    INTEGER, INTENT(OUT) :: O(0:NFRAMES-1,0:NDIM-1,0:NREP-1)
    
    INTEGER              :: I, K, L
    INTEGER              :: N, M

    WRITE (*,'(A)') "SUBROUTINE I_DIGITALIZE"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(L,K,M,N) &
    !$OMP SHARED(D,O,NFRAMES,NBINS,BINS)    
    DO K = 0,NREP-1
        DO L = 0,NDIM-1
            DO I = 0,NFRAMES-1
                M = NBINS - 1 
                DO N = 0,NBINS
                    IF ( BINS(N,L) > D(I,L,K) ) THEN
                        M = N - 1
                        EXIT
                    END IF
                END DO
            O(I,L,K) = M
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE I_DIGITALIZE


    SUBROUTINE R_DIGITALIZE(D,BINS,NREP,NDIM,NFRAMES,NBINS,O)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)  :: NREP
    INTEGER, INTENT(IN)  :: NDIM
    INTEGER, INTENT(IN)  :: NFRAMES
    INTEGER, INTENT(IN)  :: NBINS
    
    REAL(KIND=8), INTENT(IN)     :: D(0:NFRAMES-1,0:NDIM-1,0:NREP-1)
    REAL(KIND=8), INTENT(IN)     :: BINS(0:NBINS,0:NDIM-1)
    
    INTEGER, INTENT(OUT) :: O(0:NFRAMES-1,0:NDIM-1,0:NREP-1)
    
    INTEGER              :: I, K, L
    INTEGER              :: N, M
    
    WRITE (*,'(A)') "SUBROUTINE R_DIGITALIZE"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(L,K,M,N) &
    !$OMP SHARED(D,O,NFRAMES,NBINS,BINS)    
    DO K = 0,NREP-1
        DO L = 0,NDIM-1
            DO I = 0,NFRAMES-1
                M = NBINS - 1 
                DO N = 0,NBINS
                    IF ( BINS(N,L) > D(I,L,K) ) THEN
                        M = N - 1
                        EXIT
                    END IF
                END DO
            O(I,L,K) = M
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE R_DIGITALIZE

    SUBROUTINE I_TO1DIM(D,BINS,NFRAMES,NDIM,NREP,NBINS,O)
    
    INTEGER, INTENT(IN)  :: NREP
    INTEGER, INTENT(IN)  :: NDIM
    INTEGER, INTENT(IN)  :: NFRAMES
    !INTEGER, INTENT(IN)  :: NBINS(0:NDIM-1) ! TO IMPLEMENT A DIFFERENT NUMBER OF BINS PER DIM
    INTEGER, INTENT(IN)  :: NBINS
    INTEGER, INTENT(IN)  :: D(0:NFRAMES-1,0:NDIM-1,0:NREP-1)
    REAL(KIND=8), INTENT(IN)     :: BINS(0:NBINS,0:NDIM-1)
    
    INTEGER, INTENT(OUT) :: O(0:NFRAMES-1,0:NREP-1)
    
    INTEGER              :: PROD(0:NDIM)
    INTEGER              :: I, K, L, N ,M
    
    PROD(0)   = 1
    DO I = 1,NDIM-1
        PROD(I) = PROD(I-1) * NBINS
    END DO

    WRITE (*,'(A)') "SUBROUTINE I_TO1DIM"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(I,L,K,M,N) &
    !$OMP SHARED(D,O,NFRAMES,NBINS,BINS,NREP)    
    DO K = 0,NREP-1
        DO L = 0,NDIM-1
            DO I = 0,NFRAMES-1
                M = NBINS - 1 
                DO N = 0,NBINS
                    IF ( BINS(N,L) > D(I,L,K) ) THEN
                        M = N - 1
                        EXIT
                    END IF
                END DO
                O(I,K) = O(I,K) + INT( M * PROD(L) )
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE I_TO1DIM

    SUBROUTINE R_TO1DIM(D,BINS,NFRAMES,NDIM,NREP,NBINS,O)
    
    INTEGER, INTENT(IN)  :: NREP
    INTEGER, INTENT(IN)  :: NDIM
    INTEGER, INTENT(IN)  :: NFRAMES
    !INTEGER, INTENT(IN)  :: NBINS(0:NDIM-1) ! TO IMPLEMENT A DIFFERENT NUMBER OF BINS PER DIM
    INTEGER, INTENT(IN)  :: NBINS
    REAL(KIND=8), INTENT(IN)     :: D(0:NFRAMES-1,0:NDIM-1,0:NREP-1)
    REAL(KIND=8), INTENT(IN)     :: BINS(0:NBINS,0:NDIM-1)
    
    INTEGER, INTENT(OUT) :: O(0:NFRAMES-1,0:NREP-1)
    
    INTEGER              :: PROD(0:NDIM)
    INTEGER              :: I, K, L, N, M
    
    PROD(0)   = 1
    DO I = 1,NDIM-1
        PROD(I) = PROD(I-1) * NBINS
    END DO

    WRITE (*,'(A)') "SUBROUTINE R_TO1DIM"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(I,L,K,M,N) &
    !$OMP SHARED(D,O,NFRAMES,NBINS,BINS,NREP)     
    DO K = 0,NREP-1
        DO L = 0,NDIM-1
            DO I = 0,NFRAMES-1
                M = NBINS - 1 
                DO N = 0,NBINS
                    IF ( BINS(N,L) > D(I,L,K) ) THEN
                        M = N - 1
                        EXIT
                    END IF
                END DO
                O(I,K) = O(I,K) + INT( M * PROD(L) )
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE R_TO1DIM
    
    SUBROUTINE TRAJ(D,TIME,NBINS,NFRAMES,NDIM,NREP,O,O1)
    
    INTEGER, INTENT(IN)  :: NREP
    INTEGER, INTENT(IN)  :: NDIM
    INTEGER, INTENT(IN)  :: NFRAMES
    INTEGER, INTENT(IN)  :: TIME
    !INTEGER, INTENT(IN)  :: NBINS(0:NDIM-1) ! TO IMPLEMENT A DIFFERENT NUMBER OF BINS PER DIM
    INTEGER, INTENT(IN)  :: NBINS(0:NDIM-1)
    INTEGER, INTENT(IN)  :: D(0:NFRAMES-1,0:NDIM-1,0:NREP-1)
    
    INTEGER, INTENT(OUT) ::  O(0:NFRAMES-TIME-1,0:NREP-1)
    INTEGER, INTENT(OUT) :: O1(0:NFRAMES-TIME-1,0:NREP-1)
    
    INTEGER              :: PROD(0:NDIM)
    INTEGER              :: PROD_T(0:TIME)
    INTEGER              :: HASH, HASH_1
    INTEGER              :: I, K, L, N

    WRITE (*,'(A)') "SUBROUTINE TRAJ"
    PROD(0)   = 1
    PROD_T(0) = 1
    DO I = 1,NDIM
        PROD(I) = PROD(I-1) * NBINS(I-1)
    END DO
    DO I = 1,TIME
        PROD_T(I) = PROD_T(I-1) * PROD(NDIM)
    END DO
#ifdef OMP
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif      
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(HASH,HASH_1,K,I,L,N) &
    !$OMP SHARED(D,O,O1,PROD,PROD_T,TIME,NREP,NFRAMES,NDIM)    
    DO K = 0,NREP-1
        DO I = 0,NFRAMES-TIME-1
            HASH = 0
            HASH_1 = 0
            DO L = 0,NDIM-1
                DO N = 0,TIME-1
                    HASH   = HASH + ( D(I+N,L,K) * PROD(L) ) * PROD_T(N)
                END DO
                HASH_1 = HASH_1 + ( D(I+TIME,L,K) * PROD(L) ) * PROD_T(TIME)
            END DO
            O(I,K)  = HASH
            O1(I,K) = O(I,K) + HASH_1
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE TRAJ

    SUBROUTINE TRAJ_SIMP(D,TIME,NBINS,NFRAMES,NDIM,NREP,O)
    
    INTEGER, INTENT(IN)  :: NREP
    INTEGER, INTENT(IN)  :: NDIM
    INTEGER, INTENT(IN)  :: NFRAMES
    INTEGER, INTENT(IN)  :: TIME
    INTEGER, INTENT(IN)  :: NBINS(0:NDIM-1)
    INTEGER, INTENT(IN)  :: D(0:NFRAMES-1,0:NDIM-1,0:NREP-1)
    
    INTEGER, INTENT(OUT) ::  O(0:NFRAMES-TIME-1,0:NREP-1)
    
    INTEGER              :: PROD(0:NDIM)
    INTEGER              :: PROD_T(0:TIME-1)
    INTEGER              :: HASH
    INTEGER              :: I, K, L, N

    WRITE (*,'(A)') "SUBROUTINE TRAJ_SIMP"
    PROD(0)   = 1
    PROD_T(0) = 1
    DO I = 1,NDIM
        PROD(I) = PROD(I-1) * NBINS(I-1)
    END DO
    DO I = 1,TIME-1
        PROD_T(I) = PROD_T(I-1) * PROD(NDIM)
    END DO


#ifdef OMP
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif      
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(HASH,K,I,L,N) &
    !$OMP SHARED(D,O,PROD,PROD_T,TIME,NREP,NFRAMES,NDIM)    
    DO K = 0,NREP-1
        DO I = 0,NFRAMES-TIME-1
            HASH = 0
            DO L = 0,NDIM-1 ! MAYBE WE CAN INTERCHANGE THE "L" AND THE "I" LOOP FOR EFFICIENCY.
                            ! IN THIS CASE WE SHOULD GET RID OF THE "HASH" VARIABLE
                            ! AND DIRECTLY ACCUMULATE RESULTS ON THE "O" ARRAY
                            ! BUT THIS MUST ALSO INVOLVE A PREVIOUS INITIALIZATION OF THE "O" ARRAY
                DO N = 0,TIME-1
                    HASH   = HASH + ( D(I+N,L,K) * PROD(L) ) * PROD_T(N)
                END DO
            END DO
            O(I,K)  = HASH
        END DO
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE TRAJ_SIMP
    
    SUBROUTINE PROBDEF_TRAJ(D,BINS,N,NBINS,PROB)
    
    INTEGER, DIMENSION(N),INTENT(IN)                :: D
    INTEGER, INTENT(IN)                             :: N
    INTEGER, INTENT(IN)                             :: NBINS
    REAL(KIND=8), DIMENSION(0:NBINS-1),INTENT(INOUT)        :: PROB
    REAL(KIND=8), DIMENSION(0:NBINS),INTENT(INOUT)          :: BINS
    
    INTEGER               :: I, K, X
    REAL(KIND=8)                  :: P 
    
    P = 1.0 / N
    
    FORALL(I=0:NBINS-1) PROB(I)=0.0
    DO I=1,N
        X = NBINS-1
        DO K=0,NBINS-1
            IF ( BINS(K) > D(I) ) THEN
                X = K - 1
                EXIT
            END IF
        END DO
        PROB(X) = PROB(X) + P
    END DO
        
    RETURN
    
    END SUBROUTINE PROBDEF_TRAJ
    
    SUBROUTINE PROBDEF2D_TRAJ(D_X,D_Y,BINS_X,BINS_Y,N,NBINS_X,NBINS_Y,PROB)
    
    INTEGER, DIMENSION(N),INTENT(IN)                       :: D_X, D_Y
    INTEGER, INTENT(IN)                                    :: N
    INTEGER, INTENT(IN)                                    :: NBINS_X, NBINS_Y
    REAL(KIND=8), DIMENSION(0:NBINS_X-1,0:NBINS_Y-1) ,INTENT(OUT)  :: PROB
    REAL(KIND=8), DIMENSION(0:NBINS_X), INTENT(IN)                 :: BINS_X
    REAL(KIND=8), DIMENSION(0:NBINS_Y), INTENT(IN)                 :: BINS_Y
    
    INTEGER               :: I, J, K, X, Y
    REAL(KIND=8)                  :: P 
    
    P = 1.0 / N
    
    FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    DO I=1,N
        X = NBINS_X-1
        Y = NBINS_Y-1
        DO K=1,NBINS_X-1
            IF ( BINS_X(K) > D_X(I) )  THEN
                X = K - 1
                EXIT
            END IF
        END DO
        DO K=1,NBINS_Y-1
            IF ( BINS_Y(K) > D_Y(I) ) THEN
                Y = K - 1
                EXIT
            END IF
        END DO
        PROB(X,Y) = PROB(X,Y) + P
    END DO
        
    RETURN
    
    END SUBROUTINE PROBDEF2D_TRAJ
    
    SUBROUTINE ENTROPY_TRAJ(D,NFRAMES,NREP,E)
    
    INTEGER, INTENT(IN)             :: NFRAMES
    INTEGER, INTENT(IN)             :: NREP
    INTEGER, INTENT(IN)             :: D(0:NFRAMES-1,0:NREP-1)
    REAL(KIND=8), INTENT(OUT)       :: E(0:NREP-1)
    
    INTEGER, DIMENSION(NFRAMES)     :: D_TEMP
    INTEGER, DIMENSION(NFRAMES)     :: I_BINS
    INTEGER                         :: NBINS
    INTEGER                         :: I, J
    INTEGER                         :: K, X
    REAL(KIND=8)                    :: P
    
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: BINS
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: P_TEMP
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "ENTROPY_TRAJ"
#ifdef OMP
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif      
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP,BINS,P_TEMP,I_BINS,NBINS,K) &
    !$OMP SHARED(D,E,NFRAMES,P)
    DO I = 0,NREP-1
        D_TEMP(1:NFRAMES) = (/ (D(K,I), K=0,NFRAMES-1) /)
        CALL UNIRNK(D_TEMP,I_BINS,NBINS)
        ALLOCATE(BINS(0:NBINS))
        ALLOCATE(P_TEMP(0:NBINS-1))
        BINS(0)         = D_TEMP(I_BINS(1))
        BINS(1:NBINS-1) = (/ ( FLOAT( D_TEMP(I_BINS(K)) + D_TEMP(I_BINS(K-1)) ) / 2.0 , K = 2,NBINS) /)
        BINS(NBINS)     = D_TEMP(I_BINS(NBINS))
        CALL PROBDEF_TRAJ(D_TEMP,BINS,NFRAMES,NBINS,P_TEMP)
        FORALL(K=0:NBINS-1) P_TEMP(K)=0.0
        DO J=1,NFRAMES
            X = NBINS-1
            DO K=0,NBINS-1
                IF ( BINS(K) > D_TEMP(J) ) THEN
                    X = K - 1
                    EXIT
                END IF
            END DO
            P_TEMP(X) = P_TEMP(X) + P
        END DO
        DO K = 0,NBINS-1
            IF (P_TEMP(K) > 1) THEN
                E(I) = E(I) - LOG(P_TEMP(K)) * P_TEMP(K)
            END IF
        END DO
        DEALLOCATE(P_TEMP)
        DEALLOCATE(BINS)
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE ENTROPY_TRAJ
    
    SUBROUTINE ENTROPY_TRAJ_WEIGHT(D,W,NFRAMES,NREP,E)
    
    INTEGER, INTENT(IN)             :: NFRAMES
    INTEGER, INTENT(IN)             :: NREP
    INTEGER, INTENT(IN)             :: D(0:NFRAMES-1,0:NREP-1)
    REAL(KIND=8), INTENT(IN)                :: W(0:NFRAMES-1)
    REAL(KIND=8), INTENT(OUT)               :: E(0:NREP-1)
    
    INTEGER, DIMENSION(NFRAMES)          :: D_TEMP
    INTEGER, DIMENSION(NFRAMES)              :: I_BINS
    INTEGER                                  :: NBINS
    INTEGER                                  :: I, J
    INTEGER                                  :: K, X
    REAL(KIND=8)                             :: P, SUM_W, L_SUM_W
    
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE          :: BINS
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE          :: P_TEMP
    
    P = 1.0 / FLOAT(NFRAMES)
    SUM_W = SUM(W)
    L_SUM_W = LOG(SUM_W)
    
    WRITE (*,'(A)') "ENTROPY_TRAJ_WEIGHT"
#ifdef OMP
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif      
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP,BINS,P_TEMP,I_BINS,NBINS,K) &
    !$OMP SHARED(D,E,NFRAMES,P,W,SUM_W)
    DO I = 0,NREP-1
        D_TEMP(1:NFRAMES) = (/ (D(J,I), J=0,NFRAMES-1) /)
        CALL UNIRNK(D_TEMP,I_BINS,NBINS)
        ALLOCATE(BINS(0:NBINS))
        ALLOCATE(P_TEMP(0:NBINS-1))
        BINS(0)         = D_TEMP(I_BINS(1))
        BINS(1:NBINS-1) = (/ ( FLOAT( D_TEMP(I_BINS(K)) + D_TEMP(I_BINS(K-1)) ) / 2.0 , K = 2,NBINS) /)
        BINS(NBINS)     = D_TEMP(I_BINS(NBINS))
        FORALL(J=0:NBINS-1) P_TEMP(J)=0.0
        DO J=1,NFRAMES
            X = NBINS-1
            DO K=0,NBINS-1
                IF ( BINS(K) > D_TEMP(J) ) THEN
                    X = K - 1
                    EXIT
                END IF
            END DO
            P_TEMP(X) = P_TEMP(X) + P * W(J)
        END DO
        DO K = 0,NBINS-1
            IF (P_TEMP(K) > 0) THEN
                E(I) = E(I) - LOG(P_TEMP(K)) * P_TEMP(K)
            END IF
        END DO
        E(I) = E(I)/SUM_W + L_SUM_W
        DEALLOCATE(P_TEMP)
        DEALLOCATE(BINS)
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE ENTROPY_TRAJ_WEIGHT
    
    SUBROUTINE MUTUALINFO_TRAJ(D,E1,NFRAMES,NREP,M,EJ)

    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE          :: BINS_X
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE          :: BINS_Y
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE        :: P_TEMP
    INTEGER, DIMENSION(NFRAMES)          :: D_TEMP1
    INTEGER, DIMENSION(NFRAMES)          :: D_TEMP2
    INTEGER, DIMENSION(NFRAMES)              :: I_BINS_X
    INTEGER, DIMENSION(NFRAMES)              :: I_BINS_Y

    INTEGER                                  :: NBINS_X, NBINS_Y
    INTEGER                                  :: I, J, K, L
    INTEGER                                  :: X, Y
    REAL(KIND=8)                             :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: EJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE MUTUALINFO_TRAJ"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K,NBINS_X,NBINS_Y,BINS_X,BINS_Y,I_BINS_X,I_BINS_Y) &
    !$OMP SHARED(D,E1,EJ,M,NFRAMES)
    DO I = 0,NREP-2
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        !
        ! SINCE WE ARE USING THE UNIRNK ROUTINE (FROM THE DATAPACK LIBRARY) 
        ! THE TEMPORARY TRAJECTORY HAS TO BE USED THE 1-BASED INDEXING
        ! THEN, TO CALCULATE JOINT DISTRIBUTION, WHE USE A MODIFIED 
        ! PROBDEF2D_TRAJ ROUTINE.
        !  
        D_TEMP1(1:NFRAMES) = (/ (D(K,I), K=0,NFRAMES-1) /) 
        CALL UNIRNK(D_TEMP1,I_BINS_X,NBINS_X)
        ALLOCATE(BINS_X(0:NBINS_X))
        BINS_X(0) = D_TEMP1(I_BINS_X(1))
        DO K = 2,NBINS_X
            BINS_X(K) = ( D_TEMP1(I_BINS_X(K)) + D_TEMP1(I_BINS_X(K-1)) ) / 2.0
        END DO
        BINS_X(NBINS_X) = D_TEMP1(I_BINS_X(NBINS_X)) 
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            D_TEMP2(1:NFRAMES) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            CALL UNIRNK(D_TEMP2,I_BINS_Y,NBINS_Y)
            ALLOCATE(BINS_Y(0:NBINS_Y))
            BINS_Y(0) = D_TEMP2(I_BINS_Y(1))
            DO K = 2,NBINS_Y
                BINS_Y(K) = ( D_TEMP2(I_BINS_Y(K)) + D_TEMP2(I_BINS_Y(K-1)) ) / 2.0
            END DO
            BINS_Y(NBINS_Y) = D_TEMP2(I_BINS_Y(NBINS_Y))
            ALLOCATE(P_TEMP(0:NBINS_X-1,0:NBINS_Y-1))
            FORALL(J=0:NBINS_Y-1, I=0:NBINS_X-1) P_TEMP(I,J)=0.0
            DO K=1,NFRAMES
                X = NBINS_X-1
                Y = NBINS_Y-1
                DO L=1,NBINS_X-1
                    IF ( BINS_X(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS_Y-1
                    IF ( BINS_Y(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(X,Y) = P_TEMP(X,Y) + P
            END DO
            DO K = 0,NBINS_Y-1
                DO L = 0,NBINS_X-1
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
            DEALLOCATE(BINS_Y)
            DEALLOCATE(P_TEMP)
        END DO
        DEALLOCATE(BINS_X)
    END DO
    !$OMP END PARALLEL DO
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE MUTUALINFO_TRAJ
    
    SUBROUTINE MUTUALINFO_TRAJ_WEIGHT(D,E1,W,NFRAMES,NREP,M,EJ)

    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1)             :: W
    
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE          :: BINS_X
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE          :: BINS_Y
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE        :: P_TEMP
    INTEGER, DIMENSION(NFRAMES)              :: D_TEMP1
    INTEGER, DIMENSION(NFRAMES)              :: D_TEMP2
    INTEGER, DIMENSION(NFRAMES)              :: I_BINS_X
    INTEGER, DIMENSION(NFRAMES)              :: I_BINS_Y

    INTEGER                                  :: NBINS_X, NBINS_Y
    INTEGER                                  :: I, J, K, L
    INTEGER                                  :: X, Y
    REAL(KIND=8)                             :: P, SUM_W, L_SUM_W
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: EJ
    
    P = 1.0 / FLOAT(NFRAMES)
    SUM_W = SUM(W)
    L_SUM_W = LOG(SUM_W)
    
    WRITE (*,'(A)') "SUBROUTINE MUTUALINFO_TRAJ_WEIGHT"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K,NBINS_X,NBINS_Y,BINS_X,BINS_Y,I_BINS_X,I_BINS_Y) &
    !$OMP SHARED(D,E1,EJ,M,NFRAMES)
    DO I = 0,NREP-2
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        !
        ! SINCE WE ARE USING THE UNIRNK ROUTINE (FROM THE DATAPACK LIBRARY) 
        ! THE TEMPORARY TRAJECTORY HAS TO BE USED THE 1-BASED INDEXING
        ! THEN, TO CALCULATE JOINT DISTRIBUTION, WHE USE A MODIFIED 
        ! PROBDEF2D_TRAJ ROUTINE.
        !  
        D_TEMP1(1:NFRAMES) = (/ (D(K,I), K=0,NFRAMES-1) /) 
        CALL UNIRNK(D_TEMP1,I_BINS_X,NBINS_X)
        ALLOCATE(BINS_X(0:NBINS_X))
        BINS_X(0) = D_TEMP1(I_BINS_X(1))
        DO K = 2,NBINS_X
            BINS_X(K) = ( D_TEMP1(I_BINS_X(K)) + D_TEMP1(I_BINS_X(K-1)) ) / 2.0
        END DO
        BINS_X(NBINS_X) = D_TEMP1(I_BINS_X(NBINS_X)) 
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            D_TEMP2(1:NFRAMES) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            CALL UNIRNK(D_TEMP2,I_BINS_Y,NBINS_Y)
            ALLOCATE(BINS_Y(0:NBINS_Y))
            BINS_Y(0) = D_TEMP2(I_BINS_Y(1))
            DO K = 2,NBINS_Y
                BINS_Y(K) = ( D_TEMP2(I_BINS_Y(K)) + D_TEMP2(I_BINS_Y(K-1)) ) / 2.0
            END DO
            BINS_Y(NBINS_Y) = D_TEMP2(I_BINS_Y(NBINS_Y))
            ALLOCATE(P_TEMP(0:NBINS_X-1,0:NBINS_Y-1))
            FORALL(L=0:NBINS_Y-1, K=0:NBINS_X-1) P_TEMP(K,L)=0.0
            DO K=1,NFRAMES
                X = NBINS_X-1
                Y = NBINS_Y-1
                DO L=1,NBINS_X-1
                    IF ( BINS_X(L) > D_TEMP1(K) )  THEN
                        X = K - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS_Y-1
                    IF ( BINS_Y(L) > D_TEMP2(K) ) THEN
                        Y = K - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(X,Y) = P_TEMP(X,Y) + P * W(K-1) ! 1-BASED INDEXING
            END DO
            DO K = 0,NBINS_Y-1
                DO L = 0,NBINS_X-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(J,I) = EJ(J,I)/SUM_W + L_SUM_W
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
            DEALLOCATE(BINS_Y)
            DEALLOCATE(P_TEMP)
        END DO
        DEALLOCATE(BINS_X)
    END DO
    !$OMP END PARALLEL DO
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE MUTUALINFO_TRAJ_WEIGHT
    
    SUBROUTINE MUTUALINFO_OTHER_TRAJ(D1,D2,E1,E2,NFRAMES,NREP1,NREP2,M,EJ)
    
    
    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN)                                  :: NREP1,NREP2
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1):: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1):: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)               :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)               :: E2
    
    INTEGER, DIMENSION(NFRAMES)                          :: D_TEMP1
    INTEGER, DIMENSION(NFRAMES)                          :: D_TEMP2
    INTEGER, DIMENSION(NFRAMES)                          :: I_BINS1
    INTEGER, DIMENSION(NFRAMES)                          :: I_BINS2
    INTEGER                                              :: NBINS1,NBINS2
    
    REAL(KIND=8), ALLOCATABLE                 :: BINS1(:)
    REAL(KIND=8), ALLOCATABLE                 :: BINS2(:)
    REAL(KIND=8), ALLOCATABLE                 :: P_TEMP(:,:)
    INTEGER                           :: I, J, K, L
    INTEGER                           :: X, Y
    REAL(KIND=8)                      :: P
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)    :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)    :: EJ
    
    P = 1.0/FLOAT(NFRAMES)
    
    WRITE (*,'(A)') "SUBROUTINE MUTUALINFO_OTHER_TRAJ"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO & 
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K,X,Y,NBINS1,NBINS2,BINS1,BINS2,I_BINS1,I_BINS2) &
    !$OMP SHARED(D1,D2,E1,EJ,M,NFRAMES,P)
    DO I = 0,NREP1-1
        D_TEMP1(1:NFRAMES) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        CALL UNIRNK(D_TEMP1,I_BINS1,NBINS1)
        ALLOCATE(BINS1(0:NBINS1))
        BINS1(0) = D_TEMP1(I_BINS1(1))
        DO K = 2,NBINS1
            BINS1(K) = ( D_TEMP1(I_BINS1(K)) + D_TEMP1(I_BINS1(K-1)) ) / 2.0
        END DO
        BINS1(NBINS1) = D_TEMP1(I_BINS1(NBINS1)) 
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(1:NFRAMES) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            CALL UNIRNK(D_TEMP2,I_BINS2,NBINS2)
            ALLOCATE(BINS2(0:NBINS2))
            BINS2(0) = D_TEMP2(I_BINS2(1))
            DO K = 2,NBINS2
                BINS2(K) = ( D_TEMP2(I_BINS2(K)) + D_TEMP2(I_BINS2(K-1)) ) / 2.0
            END DO
            BINS2(NBINS2) = D_TEMP2(I_BINS2(NBINS2))
            ALLOCATE(P_TEMP(0:NBINS1-1,0:NBINS2-1))
            FORALL(J=0:NBINS2-1, I=0:NBINS1-1) P_TEMP(I,J)=0.0
            DO K=1,NFRAMES
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = L - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = L - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(X,Y) = P_TEMP(X,Y) + P
            END DO
            DO K = 0,NBINS2-1
                DO L = 0,NBINS1-1
                    IF (P_TEMP(L,K) > 1) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
            DEALLOCATE(BINS2)
            DEALLOCATE(P_TEMP)
        END DO
        DEALLOCATE(BINS1)
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE MUTUALINFO_OTHER_TRAJ
    
    SUBROUTINE MUTUALINFO_OTHER_TRAJ_WEIGHT(D1,D2,E1,E2,W,NFRAMES,NREP1,NREP2,M,EJ)
    
    
    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN)                                  :: NREP1,NREP2
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1):: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1):: D2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP1-1)               :: E1
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NREP2-1)               :: E2
    REAL(KIND=8), INTENT(IN), DIMENSION(0:NFRAMES-1)             :: W
    
    INTEGER, DIMENSION(NFRAMES)                          :: D_TEMP1
    INTEGER, DIMENSION(NFRAMES)                          :: D_TEMP2
    INTEGER, DIMENSION(NFRAMES)                          :: I_BINS1
    INTEGER, DIMENSION(NFRAMES)                          :: I_BINS2
    INTEGER                                              :: NBINS1,NBINS2
    
    REAL(KIND=8), ALLOCATABLE                 :: BINS1(:)
    REAL(KIND=8), ALLOCATABLE                 :: BINS2(:)
    REAL(KIND=8), ALLOCATABLE                 :: P_TEMP(:,:)
    INTEGER                           :: I, J, K, L
    INTEGER                           :: X, Y
    REAL(KIND=8)                      :: P, SUM_W, L_SUM_W
    
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)    :: M
    REAL(KIND=8), INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)    :: EJ
    
    P = 1.0 / FLOAT(NFRAMES)
    SUM_W = SUM(W)
    L_SUM_W = LOG(SUM_W)
    
    WRITE (*,'(A)') "SUBROUTINE MUTUALINFO_OTHER_TRAJ"
#ifdef OMP
    WRITE (*,'(A,I5)') "LAUNCHING THREADS : ", NUM_THREADS
    CALL OMP_SET_NUM_THREADS(NUM_THREADS)
#endif
    !$OMP PARALLEL DO & 
    !$OMP PRIVATE(D_TEMP1,D_TEMP2,P_TEMP,I,J,L,K,X,Y,NBINS1,NBINS2,BINS1,BINS2,I_BINS1,I_BINS2) &
    !$OMP SHARED(D1,D2,E1,EJ,M,NFRAMES,SUM_W,L_SUM_W,W)
    DO I = 0,NREP1-1
        D_TEMP1(1:NFRAMES) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        CALL UNIRNK(D_TEMP1,I_BINS1,NBINS1)
        ALLOCATE(BINS1(0:NBINS1))
        BINS1(0) = D_TEMP1(I_BINS1(1))
        DO K = 2,NBINS1
            BINS1(K) = ( D_TEMP1(I_BINS1(K)) + D_TEMP1(I_BINS1(K-1)) ) / 2.0
        END DO
        BINS1(NBINS1) = D_TEMP1(I_BINS1(NBINS1)) 
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(1:NFRAMES) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            CALL UNIRNK(D_TEMP2,I_BINS2,NBINS2)
            ALLOCATE(BINS2(0:NBINS2))
            BINS2(0) = D_TEMP2(I_BINS2(1))
            DO K = 2,NBINS2
                BINS2(K) = ( D_TEMP2(I_BINS2(K)) + D_TEMP2(I_BINS2(K-1)) ) / 2.0
            END DO
            BINS2(NBINS2) = D_TEMP2(I_BINS2(NBINS2))
            ALLOCATE(P_TEMP(0:NBINS1-1,0:NBINS2-1))
            FORALL(K=0:NBINS2-1, L=0:NBINS1-1) P_TEMP(L,K)=0.0
            DO K=1,NFRAMES
                X = NBINS1-1
                Y = NBINS2-1
                DO L=1,NBINS1-1
                    IF ( BINS1(L) > D_TEMP1(K) )  THEN
                        X = K - 1
                        EXIT
                    END IF
                END DO
                DO L=1,NBINS2-1
                    IF ( BINS2(L) > D_TEMP2(K) ) THEN
                        Y = K - 1
                        EXIT
                    END IF
                END DO
                P_TEMP(X,Y) = P_TEMP(X,Y) + P * W(K-1) ! 1-BASED INDEXING
            END DO
            DO K = 0,NBINS2-1
                DO L = 0,NBINS1-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - LOG(P_TEMP(L,K)) * P_TEMP(L,K)
                    END IF
                END DO
            END DO
            EJ(J,I) = EJ(J,I)/SUM_W + L_SUM_W
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
            DEALLOCATE(BINS2)
            DEALLOCATE(P_TEMP)
        END DO
        DEALLOCATE(BINS1)
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE MUTUALINFO_OTHER_TRAJ_WEIGHT
      
END MODULE
