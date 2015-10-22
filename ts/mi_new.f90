MODULE MI
    INTEGER, PARAMETER :: KDP = SELECTED_REAL_KIND(15)
    INTERFACE UNIRNK
      MODULE PROCEDURE D_UNIRNK, R_UNIRNK, I_UNIRNK
    END INTERFACE UNIRNK
    INTERFACE NEARLESS
      MODULE PROCEDURE D_NEARLESS, R_NEARLESS, I_NEARLESS
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
    INTERFACE MUTUALINFO_PROB
        MODULE PROCEDURE R_MUTUALINFO_PROB, I_MUTUALINFO_PROB
    END INTERFACE MUTUALINFO_PROB
    INTERFACE MUTUALINFO_OTHER
        MODULE PROCEDURE R_MUTUALINFO_OTHER, I_MUTUALINFO_OTHER, RI_MUTUALINFO_OTHER, IR_MUTUALINFO_OTHER
    END INTERFACE MUTUALINFO_OTHER
    INTERFACE MUTUALINFO_SIMP
        MODULE PROCEDURE I_MUTUALINFO_SIMP, R_MUTUALINFO2_SIMP
    END INTERFACE MUTUALINFO_SIMP
    INTERFACE MUTUALINFO2_SIMP
        MODULE PROCEDURE I_MUTUALINFO2_SIMP, R_MUTUALINFO2_SIMP
    END INTERFACE MUTUALINFO2_SIMP
    CONTAINS

    SUBROUTINE D_UNIRNK (XVALT, IRNGT, NUNI)
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
          REAL (KIND=KDP), DIMENSION (:), INTENT (IN) :: XVALT
          INTEGER, DIMENSION (:), INTENT (OUT) :: IRNGT
          INTEGER, INTENT (OUT) :: NUNI
    ! __________________________________________________________
          INTEGER, DIMENSION (SIZE(IRNGT)) :: JWRKT
          INTEGER :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
          INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
          REAL (KIND=KDP) :: XTST, XVALA, XVALB
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
    END SUBROUTINE D_UNIRNK
    
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
          REAL, DIMENSION (:), INTENT (IN) :: XVALT
          INTEGER, DIMENSION (:), INTENT (OUT) :: IRNGT
          INTEGER, INTENT (OUT) :: NUNI
    ! __________________________________________________________
          INTEGER, DIMENSION (SIZE(IRNGT)) :: JWRKT
          INTEGER :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
          INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
          REAL :: XTST, XVALA, XVALB
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
    
    FUNCTION D_NEARLESS (XVAL) RESULT (D_NL)
    !  NEAREST VALUE LESS THAN GIVEN VALUE
    ! __________________________________________________________
          REAL (KIND=KDP), INTENT (IN) :: XVAL
          REAL (KIND=KDP) :: D_NL
    ! __________________________________________________________
          D_NL = NEAREST (XVAL, -1.0_KDP)
          RETURN
    !
    END FUNCTION D_NEARLESS
    FUNCTION R_NEARLESS (XVAL) RESULT (R_NL)
    !  NEAREST VALUE LESS THAN GIVEN VALUE
    ! __________________________________________________________
          REAL, INTENT (IN) :: XVAL
          REAL :: R_NL
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
    
    INTEGER, DIMENSION(N)                     :: D_X, D_Y
    INTEGER                                   :: N
    INTEGER                                   :: NBINS_X, NBINS_Y
    REAL, DIMENSION(0:NBINS_X-1,0:NBINS_Y-1)  :: PROB
    REAL, DIMENSION(0:NBINS_X)                :: BINS_X
    REAL, DIMENSION(0:NBINS_Y)                :: BINS_Y
    
    INTEGER               :: I, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = MINVAL(D_X)
    MAX_D_X = MAXVAL(D_X)
    MIN_D_Y = MINVAL(D_Y)
    MAX_D_Y = MAXVAL(D_Y)
    
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y    
    
    BINS_X    = (/ ( MIN_D_X + I * BIN_DELTA_X, I=0,NBINS_X) /)
    BINS_Y    = (/ ( MIN_D_Y + I * BIN_DELTA_Y, I=0,NBINS_Y) /)
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    
    DO I=0,N-1
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X),NBINS_X-1)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y),NBINS_Y-1)
        PROB(X,Y) = PROB(X,Y) + 1.0 / N
    END DO
    
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=PROB(I,J)/N
    
    RETURN
    
    END SUBROUTINE I_PROB2D
    
    SUBROUTINE R_PROB2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    REAL   , DIMENSION(N)                       :: D_X, D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(0:NBINS_X-1,0:NBINS_Y-1)          :: PROB
    REAL, DIMENSION(0:NBINS_X)                :: BINS_X
    REAL, DIMENSION(0:NBINS_Y)                :: BINS_Y
    
    INTEGER               :: I, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = MINVAL(D_X)
    MAX_D_X = MAXVAL(D_X)
    MIN_D_Y = MINVAL(D_Y)
    MAX_D_Y = MAXVAL(D_Y)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    
    BINS_X    = (/ ( MIN_D_X + I * BIN_DELTA_X, I=0,NBINS_X) /)
    BINS_Y    = (/ ( MIN_D_Y + I * BIN_DELTA_Y, I=0,NBINS_Y) /)
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    
    DO I=0,N-1
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X),NBINS_X-1)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y),NBINS_Y-1)
        PROB(X,Y) = PROB(X,Y) + 1.0 / N 
    END DO
    
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=PROB(I,J)/N
    
    RETURN
        
    END SUBROUTINE R_PROB2D

    SUBROUTINE II_PROBDEF2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    INTEGER, DIMENSION(N)                       :: D_X, D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(0:NBINS_X-1,0:NBINS_Y-1)    :: PROB
    REAL, DIMENSION(0:NBINS_X)                :: BINS_X
    REAL, DIMENSION(0:NBINS_Y)                :: BINS_Y
    
    INTEGER               :: I, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = BINS_X(0)
    MAX_D_X = BINS_X(NBINS_X)
    MIN_D_Y = BINS_Y(0)
    MAX_D_Y = BINS_Y(NBINS_Y)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    
    DO I=0,N-1
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X),NBINS_X-1)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y),NBINS_Y-1)
        PROB(X,Y) = PROB(X,Y) + 1.0 / N
    END DO
    
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=PROB(I,J)/N
    
    RETURN
    
    END SUBROUTINE II_PROBDEF2D
    
    SUBROUTINE RR_PROBDEF2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)
    IMPLICIT NONE
    
    REAL, DIMENSION(N)                          :: D_X, D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(0:NBINS_X-1,0:NBINS_Y-1)          :: PROB
    REAL, DIMENSION(0:NBINS_X)                :: BINS_X
    REAL, DIMENSION(0:NBINS_Y)                :: BINS_Y
    
    INTEGER               :: I, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = BINS_X(0)
    MAX_D_X = BINS_X(NBINS_X)
    MIN_D_Y = BINS_Y(0)
    MAX_D_Y = BINS_Y(NBINS_Y)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(J,I)=0.0
    
    DO I=1,N
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X),NBINS_X-1)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y),NBINS_Y-1)
        PROB(X,Y) = PROB(X,Y) + 1.0 / N
    END DO
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(J,I)=PROB(J,I)/N
    RETURN
    
    END SUBROUTINE RR_PROBDEF2D
    
    SUBROUTINE RI_PROBDEF2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    REAL,    DIMENSION(N)                       :: D_X
    INTEGER, DIMENSION(N)                       :: D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(0:NBINS_X-1,0:NBINS_Y-1)          :: PROB
    REAL, DIMENSION(0:NBINS_X)                :: BINS_X
    REAL, DIMENSION(0:NBINS_Y)                :: BINS_Y
    
    INTEGER               :: I, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = BINS_X(0)
    MAX_D_X = BINS_X(NBINS_X)
    MIN_D_Y = BINS_Y(0)
    MAX_D_Y = BINS_Y(NBINS_Y)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    
    DO I=0,N-1
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X),NBINS_X-1)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y),NBINS_Y-1)
        PROB(X,Y) = PROB(X,Y) + 1.0 / N
    END DO
    
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=PROB(I,J)/N
    
    RETURN
    
    END SUBROUTINE RI_PROBDEF2D
    
    SUBROUTINE IR_PROBDEF2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    INTEGER, DIMENSION(N)                       :: D_X
    REAL,    DIMENSION(N)                       :: D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(0:NBINS_X-1,0:NBINS_Y-1)          :: PROB
    REAL, DIMENSION(0:NBINS_X)                :: BINS_X
    REAL, DIMENSION(0:NBINS_Y)                :: BINS_Y
    
    INTEGER               :: I, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = BINS_X(0)
    MAX_D_X = BINS_X(NBINS_X)
    MIN_D_Y = BINS_Y(0)
    MAX_D_Y = BINS_Y(NBINS_Y)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0
    
    DO I=0,N-1
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X),NBINS_X-1)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y),NBINS_Y-1)
        PROB(X,Y) = PROB(X,Y) + 1.0 / N
    END DO
    
    !FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=PROB(I,J)/N
    
    RETURN
    
    END SUBROUTINE IR_PROBDEF2D

    SUBROUTINE PROBDEF2D_TRAJ(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    INTEGER, DIMENSION(N)                       :: D_X, D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(0:NBINS_X-1,0:NBINS_Y-1)    :: PROB
    REAL, DIMENSION(0:NBINS_X)                :: BINS_X
    REAL, DIMENSION(0:NBINS_Y)                :: BINS_Y
    
    INTEGER               :: I, J, K, X, Y
    
    FORALL(I=0:NBINS_X-1, J=0:NBINS_Y-1) PROB(I,J)=0.0

    DO I=0,N-1
        X = NBINS_X-1
        Y = NBINS_Y-1
        DO K=0,NBINS_X-1
            IF ( BINS_X(K) > D_X(I) )  THEN
                X = K - 1
                EXIT
            END IF
        END DO
        DO K=0,NBINS_Y-1
            IF ( BINS_Y(K) > D_Y(I) ) THEN
                Y = K - 1
                EXIT
            END IF
        END DO
        PROB(X,Y) = PROB(X,Y) + 1.0 / N
    END DO
        
    RETURN
    
    END SUBROUTINE PROBDEF2D_TRAJ    
    
    SUBROUTINE R_MUTUALINFO_PROB(D,E1,BINS,NFRAMES,NREP,NBINS,M,EJ,PJ)

    INTEGER, INTENT(IN)                                :: NBINS
    INTEGER, INTENT(IN)                                :: NFRAMES
    INTEGER, INTENT(IN)                                :: NREP
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1)  :: D
    REAL, INTENT(IN), DIMENSION(0:NREP-1)              :: E1
    REAL, INTENT(IN), DIMENSION(0:NBINS)               :: BINS
    
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS-1,0:NBINS-1)             :: P_TEMP
    INTEGER                                           :: I, J, K
    
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: EJ
    REAL, INTENT(OUT), DIMENSION(0:NBINS-1,0:NBINS-1,0:NREP-1,0:NREP-1) :: PJ

    DO I = 0,NREP-2
        !WRITE(*,'(A,I3)') "    REP :",I+1
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            !WRITE(*,'(A,I3)') "vs. REP :",J+1
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    !WRITE(*,'(A,2I3)') "INITIALIZING :",L,K
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS,NBINS,P_TEMP,BINS,BINS)
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(L,K,I,J) = P_TEMP(L,K)
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                        !WRITE(*,'(A,4I3,A,F8.3,A,F8.3)') "SCANNING :",L,K,I,J," : ",P_TEMP(L,K)," -> ", LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE R_MUTUALINFO_PROB

    SUBROUTINE I_MUTUALINFO_PROB(D,E1,BINS,NFRAMES,NREP,NBINS,M,EJ,PJ)

    INTEGER, INTENT(IN)                                  :: NBINS
    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL, INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    REAL, INTENT(IN), DIMENSION(0:NBINS)                 :: BINS
    
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS-1,0:NBINS-1)                :: P_TEMP
    INTEGER                                             :: I, J, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: EJ
    REAL, INTENT(OUT), DIMENSION(0:NBINS-1,0:NBINS-1,0:NREP-1,0:NREP-1) :: PJ

    DO I = 0,NREP-2
        !WRITE(*,'(A,I3)') "    REP :",I+1
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            !WRITE(*,'(A,I3)') "vs. REP :",J+1
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    !WRITE(*,'(A,2I3)') "INITIALIZING :",L,K
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS,NBINS,P_TEMP,BINS,BINS)
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(L,K,I,J) = P_TEMP(L,K)
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                        !WRITE(*,'(A,4I3,A,F8.3,A,F8.3)') "SCANNING :",L,K,I,J," : ",P_TEMP(L,K)," -> ", LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE I_MUTUALINFO_PROB
    
    SUBROUTINE R_MUTUALINFO_OTHER_PROB(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ,PJ)

    INTEGER, INTENT(IN)                                   :: NBINS1,NBINS2
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL, INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL, INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2
    REAL, INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL, INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
        
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS2-1,0:NBINS1-1)           :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: EJ
    REAL, INTENT(OUT), DIMENSION(0:NBINS2-1,0:NBINS1-1,0:NREP2-1,0:NREP1-1) :: PJ

    DO I = 0,NREP1-1
        !WRITE(*,'(A,I3)') "    REP :",I+1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            !WRITE(*,'(A,I3)') "vs. REP :",J+1
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    !WRITE(*,'(A,2I3)') "INITIALIZING :",L,K
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS2,NBINS1,P_TEMP,BINS2,BINS1)
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(L,K,I,J) = P_TEMP(L,K)
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                        !WRITE(*,'(A,4I3,A,F8.3,A,F8.3)') "SCANNING :",L,K,I,J," : ",P_TEMP(L,K)," -> ", LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    
    END SUBROUTINE R_MUTUALINFO_OTHER_PROB

    SUBROUTINE I_MUTUALINFO_OTHER_PROB(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ,PJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1) :: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1) :: D2
    REAL, INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL, INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    REAL, INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL, INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2

    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS2-1,0:NBINS1-1)              :: P_TEMP
    INTEGER                                             :: I, J, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: EJ
    REAL, INTENT(OUT), DIMENSION(0:NBINS2-1,0:NBINS1-1,0:NREP2-1,0:NREP1-1) :: PJ

    DO I = 0,NREP1-1
        !WRITE(*,'(A,I3)') "    REP :",I+1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            !WRITE(*,'(A,I3)') "vs. REP :",J+1
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    !WRITE(*,'(A,2I3)') "INITIALIZING :",L,K
                    PJ(L,K,J,I) = 0.0
                END DO
            END DO 
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS2,NBINS1,P_TEMP,BINS2,BINS1)
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    PJ(L,K,J,I) = P_TEMP(L,K)
                    PJ(L,K,I,J) = P_TEMP(L,K)
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                        !WRITE(*,'(A,4I3,A,F8.3,A,F8.3)') "SCANNING :",L,K,I,J," : ",P_TEMP(L,K)," -> ", LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    
    END SUBROUTINE I_MUTUALINFO_OTHER_PROB

    SUBROUTINE I_MUTUALINFO_SIMP(D,E1,BINS,NFRAMES,NREP,NBINS,NREP1,M,EJ)

    INTEGER, INTENT(IN)                                  :: NREP1
    INTEGER, INTENT(IN)                                  :: NBINS
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN)                                  :: NFRAMES
    
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL, INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    REAL, INTENT(IN), DIMENSION(0:NBINS)                 :: BINS
    
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS-1,0:NBINS-1)             :: P_TEMP
    INTEGER                                          :: I, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP-1)                     :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP-1)                     :: EJ

    D_TEMP1(0:NFRAMES-1) = (/ (D(K,NREP1), K=0,NFRAMES-1) /)
    EJ(NREP1) = E1(NREP1)
    M(NREP1)  = E1(NREP1)
    
    DO I = 0,NREP-1
        !WRITE(*,'(A,I3)') "    REP :",I+1
        IF (I /= NREP1) THEN
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
            EJ(I) = 0.0
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS,NBINS,P_TEMP,BINS,BINS)
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(I) = EJ(I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(I) = E1(I) + E1(NREP1) - EJ(I)
        END IF
    END DO
    
    END SUBROUTINE I_MUTUALINFO_SIMP
    
    SUBROUTINE R_MUTUALINFO_SIMP(D,E1,BINS,NREP1,NFRAMES,NREP,NBINS,M,EJ)
    
    INTEGER, INTENT(IN)                               :: NREP1 
    INTEGER, INTENT(IN)                               :: NBINS
    INTEGER, INTENT(IN)                               :: NREP
    INTEGER, INTENT(IN)                               :: NFRAMES
    
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL, INTENT(IN), DIMENSION(0:NREP-1)             :: E1
    REAL, INTENT(IN), DIMENSION(0:NBINS)              :: BINS

    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS-1,0:NBINS-1)             :: P_TEMP
    INTEGER                                          :: I, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP-1)                     :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP-1)                     :: EJ

    D_TEMP1(0:NFRAMES-1) = (/ (D(K,NREP1), K=0,NFRAMES-1) /)
    EJ(NREP1) = E1(NREP1)
    M(NREP1)  = E1(NREP1)
    
    DO I = 0,NREP-1
        !WRITE(*,'(A,I3)') "    REP :",I+1
        IF (I /= NREP1) THEN
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
            EJ(I) = 0.0
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS,NBINS,P_TEMP,BINS,BINS)
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(I) = EJ(I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(I) = E1(I) + E1(NREP1) - EJ(I)
        END IF
    END DO
    
    END SUBROUTINE R_MUTUALINFO_SIMP

    SUBROUTINE I_MUTUALINFO2_SIMP(D,E1,BINSX,BINSY,NFRAMES,NREP,NBINSX,NBINSY,NREP1,M,EJ,PJ)
    
    INTEGER, INTENT(IN)                               :: NREP1 
    INTEGER, INTENT(IN)                               :: NBINSX, NBINSY
    INTEGER, INTENT(IN)                               :: NREP
    INTEGER, INTENT(IN)                               :: NFRAMES
    
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL, INTENT(IN), DIMENSION(0:NREP-1)             :: E1
    REAL, INTENT(IN), DIMENSION(0:NBINSX)             :: BINSX
    REAL, INTENT(IN), DIMENSION(0:NBINSY)             :: BINSY
    
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINSY-1,0:NBINSX-1)             :: P_TEMP
    INTEGER                                          :: I, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP-1)                     :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP-1)                     :: EJ
    REAL, INTENT(OUT), DIMENSION(0:NBINSY-1,0:NBINSX-1,0:NREP-1) :: PJ

    D_TEMP1(0:NFRAMES-1) = (/ (D(K,NREP1), K=0,NFRAMES-1) /)
    EJ(NREP1) = E1(NREP1)
    M(NREP1)  = E1(NREP1)
    
    DO I = 0,NREP-1
        !WRITE(*,'(A,I3)') "    REP :",I+1
        IF (I /= NREP1) THEN
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
            EJ(I) = 0.0
            DO K = 0,NBINSX-1
                DO L = 0,NBINSY-1
                    !WRITE(*,'(A,2I3)') "INITIALIZING :",L,K
                    PJ(L,K,I) = 0.0
                END DO
            END DO 
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINSY,NBINSX,P_TEMP,BINSX,BINSY)
            DO K = 0,NBINSX-1
                DO L = 0,NBINSY-1
                    PJ(L,K,I) = P_TEMP(L,K)
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(I) = EJ(I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(I) = E1(I) + E1(NREP1) - EJ(I)
        END IF
    END DO
    
    END SUBROUTINE I_MUTUALINFO2_SIMP

    SUBROUTINE R_MUTUALINFO2_SIMP(D,E1,BINSX,BINSY,NFRAMES,NREP,NBINSX,NBINSY,NREP1,M,EJ,PJ)
    
    INTEGER, INTENT(IN)                               :: NREP1 
    INTEGER, INTENT(IN)                               :: NBINSX, NBINSY
    INTEGER, INTENT(IN)                               :: NREP
    INTEGER, INTENT(IN)                               :: NFRAMES
    
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL, INTENT(IN), DIMENSION(0:NREP-1)             :: E1
    REAL, INTENT(IN), DIMENSION(0:NBINSX)             :: BINSX
    REAL, INTENT(IN), DIMENSION(0:NBINSY)             :: BINSY
    
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINSY-1,0:NBINSX-1)             :: P_TEMP
    INTEGER                                          :: I, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP-1)                     :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP-1)                     :: EJ
    REAL, INTENT(OUT), DIMENSION(0:NBINSY-1,0:NBINSX-1,0:NREP-1) :: PJ

    D_TEMP1(0:NFRAMES-1) = (/ (D(K,NREP1), K=0,NFRAMES-1) /)
    EJ(NREP1) = E1(NREP1)
    M(NREP1)  = E1(NREP1)
    
    DO I = 0,NREP-1
        !WRITE(*,'(A,I3)') "    REP :",I+1
        IF (I /= NREP1) THEN
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
            EJ(I) = 0.0
            DO K = 0,NBINSX-1
                DO L = 0,NBINSY-1
                    !WRITE(*,'(A,2I3)') "INITIALIZING :",L,K
                    PJ(L,K,I) = 0.0
                END DO
            END DO 
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINSY,NBINSX,P_TEMP,BINSX,BINSY)
            DO K = 0,NBINSX-1
                DO L = 0,NBINSY-1
                    PJ(L,K,I) = P_TEMP(L,K)
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(I) = EJ(I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(I) = E1(I) + E1(NREP1) - EJ(I)
        END IF
    END DO
    
    END SUBROUTINE R_MUTUALINFO2_SIMP

    SUBROUTINE R_MUTUALINFO(D,E1,BINS,NFRAMES,NREP,NBINS,M,EJ)

    INTEGER, INTENT(IN)                                :: NBINS
    INTEGER, INTENT(IN)                                :: NFRAMES
    INTEGER, INTENT(IN)                                :: NREP
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1)  :: D
    REAL, INTENT(IN), DIMENSION(0:NREP-1)              :: E1
    REAL, INTENT(IN), DIMENSION(0:NBINS)               :: BINS
    
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS-1,0:NBINS-1)             :: P_TEMP
    INTEGER                                           :: I, J, K
    
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)                     :: EJ

    DO I = 0,NREP-2
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS,NBINS,P_TEMP,BINS,BINS)
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE R_MUTUALINFO

    SUBROUTINE I_MUTUALINFO(D,E1,BINS,NFRAMES,NREP,NBINS,M,EJ)

    INTEGER, INTENT(IN)                                  :: NBINS
    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL, INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    REAL, INTENT(IN), DIMENSION(0:NBINS)                 :: BINS
    
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS-1,0:NBINS-1)                :: P_TEMP
    INTEGER                                             :: I, J, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)     :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)     :: EJ

    WRITE(*,'(A)')    "SUBROUTINE I_MUTUALINFO"
    WRITE(*,'(A,I10)') " NREP  : " , NREP
    WRITE(*,'(A,I10)') " NBINS : " , NBINS
    WRITE(*,'(A,I10)') " FRAMES: " , NFRAMES
    WRITE(*,'(A,I4,A,I4)') " D: ", SIZE(D,1), "x", SIZE(D,2)
    WRITE(*,'(A,I4,A,I4)') " E: ", SIZE(E1,1)

    DO I = 0,NREP-2
        WRITE(*,'(A,I10)') " I : " , I
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO J = I+1,NREP-1
            WRITE(*,'(A,I10)') " J : " , J
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS,NBINS,P_TEMP,BINS,BINS)
            DO K = 0,NBINS-1
                DO L = 0,NBINS-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE I_MUTUALINFO
    
    SUBROUTINE MUTUALINFO_TRAJ(D,E1,NFRAMES,NREP,M,EJ)

    INTEGER, INTENT(IN)                                  :: NFRAMES
    INTEGER, INTENT(IN)                                  :: NREP
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP-1) :: D
    REAL, INTENT(IN), DIMENSION(0:NREP-1)                :: E1
    
    REAL, DIMENSION(:), ALLOCATABLE          :: BINS_X
    REAL, DIMENSION(:), ALLOCATABLE          :: BINS_Y
    REAL, DIMENSION(:,:), ALLOCATABLE        :: P_TEMP
    INTEGER, DIMENSION(0:NFRAMES-1)          :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)          :: D_TEMP2
    INTEGER, DIMENSION(0:NFRAMES-1)          :: I_BINS_X
    INTEGER, DIMENSION(0:NFRAMES-1)          :: I_BINS_Y

    INTEGER                                  :: NBINS_X, NBINS_Y
    INTEGER                                  :: I, J, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP-1,0:NREP-1)      :: EJ

    WRITE(*,'(A)')     "SUBROUTINE MUTUALINFO_TRAJ"
    WRITE(*,'(A,I10)') " NREP  : " , NREP
    WRITE(*,'(A,I10)') " FRAMES: " , NFRAMES

    DO I = 0,NREP-2
        WRITE(*,'(A,I10)') " I : " , I
        EJ(I,I) = E1(I)
        M(I,I)  = E1(I)
        D_TEMP1(0:NFRAMES-1) = (/ (D(K,I), K=0,NFRAMES-1) /)
        CALL UNIRNK(D_TEMP1,I_BINS_X,NBINS_X)
        ALLOCATE(BINS_X(0:NBINS_X-1))
        DO K = 0,NBINS_X-1
            BINS_X(K) = D_TEMP1(I_BINS_X(K))
        END DO
        DO J = I+1,NREP-1
            WRITE(*,'(A,I10)') " J : " , J
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            CALL UNIRNK(D_TEMP2,I_BINS_Y,NBINS_Y)
            ALLOCATE(BINS_Y(0:NBINS_Y-1))
            DO K = 0,NBINS_Y-1
                BINS_Y(K) = D_TEMP2(I_BINS_Y(K))
            END DO
            ALLOCATE(P_TEMP(0:NBINS_Y-1,0:NBINS_X-1))
            CALL PROBDEF2D_TRAJ(D_TEMP2,D_TEMP1,NFRAMES,NBINS_X,NBINS_Y,P_TEMP,BINS_X,BINS_Y)
            DO K = 0,NBINS_X-1
                DO L = 0,NBINS_Y-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            EJ(I,J) = EJ(J,I) 
            M(J,I)  = E1(I) + E1(J) - EJ(J,I)
            M(I,J)  = M(J,I)
        END DO
    END DO
    
    EJ(NREP-1,NREP-1) = E1(NREP-1)
    M(NREP-1,NREP-1)  = E1(NREP-1)
    
    END SUBROUTINE MUTUALINFO_TRAJ
    
    SUBROUTINE R_MUTUALINFO_OTHER(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)

    INTEGER, INTENT(IN)                                   :: NBINS1,NBINS2
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL, INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL, INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2
    REAL, INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL, INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
        
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS2-1,0:NBINS1-1)           :: P_TEMP
    INTEGER                                          :: I, J, K, L 
    
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: EJ

    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS2,NBINS1,P_TEMP,BINS2,BINS1)
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    
    END SUBROUTINE R_MUTUALINFO_OTHER

    SUBROUTINE I_MUTUALINFO_OTHER(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1) :: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1) :: D2
    REAL, INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL, INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    REAL, INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL, INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2

    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS2-1,0:NBINS1-1)              :: P_TEMP
    INTEGER                                             :: I, J, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: EJ

    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS2,NBINS1,P_TEMP,BINS2,BINS1)
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    
    END SUBROUTINE I_MUTUALINFO_OTHER
    
    SUBROUTINE IR_MUTUALINFO_OTHER(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1) :: D1
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1)    :: D2
    REAL, INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL, INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    REAL, INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL, INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2

    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP1
    REAL, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS2-1,0:NBINS1-1)              :: P_TEMP
    INTEGER                                             :: I, J, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: EJ

    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS2,NBINS1,P_TEMP,BINS2,BINS1)
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    
    END SUBROUTINE IR_MUTUALINFO_OTHER

    SUBROUTINE RI_MUTUALINFO_OTHER(D1,D2,E1,E2,BINS1,BINS2,NFRAMES,NREP1,NREP2,NBINS1,NBINS2,M,EJ)
    
    INTEGER, INTENT(IN)                                   :: NBINS1, NBINS2
    INTEGER, INTENT(IN)                                   :: NFRAMES
    INTEGER, INTENT(IN)                                   :: NREP1, NREP2
    REAL, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP1-1)    :: D1
    INTEGER, INTENT(IN), DIMENSION(0:NFRAMES-1,0:NREP2-1) :: D2
    REAL, INTENT(IN), DIMENSION(0:NREP1-1)                :: E1
    REAL, INTENT(IN), DIMENSION(0:NREP2-1)                :: E2
    REAL, INTENT(IN), DIMENSION(0:NBINS1)                 :: BINS1
    REAL, INTENT(IN), DIMENSION(0:NBINS2)                 :: BINS2

    REAL, DIMENSION(0:NFRAMES-1)                        :: D_TEMP1
    INTEGER, DIMENSION(0:NFRAMES-1)                     :: D_TEMP2
    REAL, DIMENSION(0:NBINS2-1,0:NBINS1-1)              :: P_TEMP
    INTEGER                                             :: I, J, K, L
    
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: M
    REAL, INTENT(OUT), DIMENSION(0:NREP2-1,0:NREP1-1)                       :: EJ

    DO I = 0,NREP1-1
        D_TEMP1(0:NFRAMES-1) = (/ (D1(K,I), K=0,NFRAMES-1) /)
        DO J = 0,NREP2-1
            EJ(J,I) = 0.0
            D_TEMP2(0:NFRAMES-1) = (/ (D2(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP2,D_TEMP1,NFRAMES,NBINS2,NBINS1,P_TEMP,BINS2,BINS1)
            DO K = 0,NBINS1-1
                DO L = 0,NBINS2-1
                    IF (P_TEMP(L,K) > 0) THEN
                        EJ(J,I) = EJ(J,I) - P_TEMP(L,K) * LOG(P_TEMP(L,K))
                    END IF
                END DO
            END DO
            M(J,I) = E1(I) + E2(J) - EJ(J,I)
        END DO
    END DO
    
    END SUBROUTINE RI_MUTUALINFO_OTHER
      
END MODULE
