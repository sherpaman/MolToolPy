MODULE MI

    INTERFACE PROB2D
        MODULE PROCEDURE I_PROB2D, R_PROB2D
    END INTERFACE PROB2D
    INTERFACE PROBDEF2D
        MODULE PROCEDURE II_PROBDEF2D, RR_PROBDEF2D, RI_PROBDEF2D, IR_PROBDEF2D
    END INTERFACE PROBDEF2D
    INTERFACE MUTUALINFO
        MODULE PROCEDURE R_MUTUALINFO
    END INTERFACE MUTUALINFO
    
    CONTAINS

    SUBROUTINE I_PROB2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)
  
    IMPLICIT NONE
    
    INTEGER, DIMENSION(N)                       :: D_X, D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(NBINS_X,NBINS_Y)          :: PROB
    REAL, DIMENSION(NBINS_X+1)                :: BINS_X
    REAL, DIMENSION(NBINS_Y+1)                :: BINS_Y
    
    INTEGER               :: I, J, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = MINVAL(D_X)
    MAX_D_X = MAXVAL(D_X)
    MIN_D_Y = MINVAL(D_Y)
    MAX_D_Y = MAXVAL(D_Y)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    
    BINS_X    = (/ ( MIN_D_X + (I-1) * BIN_DELTA_X, I=1,NBINS_X+1) /)
    BINS_Y    = (/ ( MIN_D_Y + (I-1) * BIN_DELTA_Y, I=1,NBINS_Y+1) /)
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=0.0
    
    DO I=1,N
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X)+1,NBINS_X)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y)+1,NBINS_Y)
        PROB(X,Y) = PROB(X,Y) + 1.0 
    END DO
    
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=PROB(I,J)/N
    
    RETURN
    
    END SUBROUTINE I_PROB2D
    
    SUBROUTINE R_PROB2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    REAL   , DIMENSION(N)                       :: D_X, D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(NBINS_X,NBINS_Y)          :: PROB
    REAL, DIMENSION(NBINS_X+1)                :: BINS_X
    REAL, DIMENSION(NBINS_Y+1)                :: BINS_Y
    
    INTEGER               :: I, J, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = MINVAL(D_X)
    MAX_D_X = MAXVAL(D_X)
    MIN_D_Y = MINVAL(D_Y)
    MAX_D_Y = MAXVAL(D_Y)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    
    BINS_X    = (/ ( MIN_D_X + (I+1) * BIN_DELTA_X, I=1,NBINS_X+1) /)
    BINS_Y    = (/ ( MIN_D_Y + (I+1) * BIN_DELTA_Y, I=1,NBINS_Y+1) /)
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=0.0
    
    DO I=1,N
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X)+1,NBINS_X)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y)+1,NBINS_Y)
        PROB(X,Y) = PROB(X,Y) + 1.0 
    END DO
    
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=PROB(I,J)/N
    
    RETURN
        
    END SUBROUTINE R_PROB2D

    SUBROUTINE II_PROBDEF2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    INTEGER, DIMENSION(N)                       :: D_X, D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(NBINS_X,NBINS_Y)          :: PROB
    REAL, DIMENSION(NBINS_X+1)                :: BINS_X
    REAL, DIMENSION(NBINS_Y+1)                :: BINS_Y
    
    INTEGER               :: I, J, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = BINS_X(1)
    MAX_D_X = BINS_X(NBINS_X+1)
    MIN_D_Y = BINS_Y(1)
    MAX_D_Y = BINS_Y(NBINS_Y+1)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=0.0
    
    DO I=1,N
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X)+1,NBINS_X)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y)+1,NBINS_Y)
        PROB(X,Y) = PROB(X,Y) + 1.0 
    END DO
    
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=PROB(I,J)/N
    
    RETURN
    
    END SUBROUTINE II_PROBDEF2D
    
    SUBROUTINE RR_PROBDEF2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)
    IMPLICIT NONE
    
    REAL, DIMENSION(N)                          :: D_X, D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(NBINS_X,NBINS_Y)          :: PROB
    REAL, DIMENSION(NBINS_X+1)                :: BINS_X
    REAL, DIMENSION(NBINS_Y+1)                :: BINS_Y
    
    INTEGER               :: I, J, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = BINS_X(1)
    MAX_D_X = BINS_X(NBINS_X+1)
    MIN_D_Y = BINS_Y(1)
    MAX_D_Y = BINS_Y(NBINS_Y+1)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=0.0
    !WRITE(*,'(A)') "PROB_XY INITIALIZED"
    DO I=1,N
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X)+1,NBINS_X)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y)+1,NBINS_Y)
        !WRITE(*,'(I8,I3,I3)'), I,X,Y
        PROB(X,Y) = PROB(X,Y) + 1.0 
    END DO
    !WRITE(*,'(A)') "PROB_XY DONE"
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=PROB(I,J)/N
    !WRITE(*,'(A)') "PROB_XY NORMALIZED"
    RETURN
    
    END SUBROUTINE RR_PROBDEF2D
    
    SUBROUTINE RI_PROBDEF2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    REAL,    DIMENSION(N)                       :: D_X
    INTEGER, DIMENSION(N)                       :: D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(NBINS_X,NBINS_Y)            :: PROB
    REAL,    DIMENSION(NBINS_X+1)               :: BINS_X
    REAL,    DIMENSION(NBINS_Y+1)               :: BINS_Y
    
    INTEGER               :: I, J, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = BINS_X(1)
    MAX_D_X = BINS_X(NBINS_X+1)
    MIN_D_Y = BINS_Y(1)
    MAX_D_Y = BINS_Y(NBINS_Y+1)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=0.0
    
    DO I=1,N
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X)+1,NBINS_X)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y)+1,NBINS_Y)
        PROB(X,Y) = PROB(X,Y) + 1.0 
    END DO
    
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=PROB(I,J)/N
    
    RETURN
    
    END SUBROUTINE RI_PROBDEF2D
    
    SUBROUTINE IR_PROBDEF2D(D_X,D_Y,N,NBINS_X,NBINS_Y,PROB,BINS_X,BINS_Y)

    IMPLICIT NONE
    
    INTEGER, DIMENSION(N)                       :: D_X
    REAL,    DIMENSION(N)                       :: D_Y
    INTEGER                                     :: N
    INTEGER                                     :: NBINS_X, NBINS_Y
    REAL, DIMENSION(NBINS_X,NBINS_Y)          :: PROB
    REAL, DIMENSION(NBINS_X+1)                :: BINS_X
    REAL, DIMENSION(NBINS_Y+1)                :: BINS_Y
    
    INTEGER               :: I, J, X, Y
    REAL                  :: MIN_D_X, MIN_D_Y, MAX_D_X, MAX_D_Y
    REAL                  :: BIN_DELTA_X, BIN_DELTA_Y
    
    MIN_D_X = BINS_X(1)
    MAX_D_X = BINS_X(NBINS_X+1)
    MIN_D_Y = BINS_Y(1)
    MAX_D_Y = BINS_Y(NBINS_Y+1)
        
    BIN_DELTA_X = ( MAX_D_X - MIN_D_X ) / NBINS_X
    BIN_DELTA_Y = ( MAX_D_Y - MIN_D_Y ) / NBINS_Y
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=0.0
    
    DO I=1,N
        X = MIN(INT((D_X(I) - MIN_D_X) / BIN_DELTA_X)+1,NBINS_X)
        Y = MIN(INT((D_Y(I) - MIN_D_Y) / BIN_DELTA_Y)+1,NBINS_Y)
        PROB(X,Y) = PROB(X,Y) + 1.0 
    END DO
    
    FORALL(I=1:NBINS_X, J=1:NBINS_Y) PROB(I,J)=PROB(I,J)/N
    
    RETURN
    
    END SUBROUTINE IR_PROBDEF2D
    
    SUBROUTINE R_MUTUALINFO(D,E1,BINS,NFRAMES,NREP,NBINS,M,E2,P2)

    INTEGER, INTENT(IN)                          :: NBINS
    INTEGER, INTENT(IN)                          :: NFRAMES
    INTEGER, INTENT(IN)                          :: NREP
    REAL, INTENT(IN), DIMENSION(NFRAMES,NREP)    :: D
    REAL, INTENT(IN), DIMENSION(NREP)            :: E1
    REAL, INTENT(IN), DIMENSION(NBINS)           :: BINS
    
    REAL, DIMENSION(NFRAMES)                     :: D_TEMP1
    REAL, DIMENSION(NFRAMES)                     :: D_TEMP2
    REAL, DIMENSION(NBINS,NBINS)                 :: P_TEMP
    INTEGER, DIMENSION(2)                        :: S
    INTEGER                                      :: I, J, K
    
    REAL, INTENT(OUT), DIMENSION(NREP,NREP)             :: M
    REAL, INTENT(OUT), DIMENSION(NREP,NREP)             :: E2
    REAL, INTENT(OUT), DIMENSION(NBINS,NBINS,NREP,NREP) :: P2

    DO I = 0,NREP-1
        E2(I,I) = 2 * E1(I)
        M(I,I) = 2 * E1(I)
        D_TEMP1(1:NFRAMES) = (/ (D(K,I), K=0,NFRAMES-1) /)
        DO L = 0,NBINS-1
            DO K = 0,NBINS-1
                P2(L,K,J,I) = 1.0
            END DO
        END DO 
        DO J = I+1,NREP-1
            D_TEMP2(1:NFRAMES) = (/ (D(K,J), K = 0,NFRAMES-1) /)
            CALL PROBDEF2D(D_TEMP1,D_TEMP2,S(1),NBINS,NBINS,P_TEMP,BINS,BINS)
            DO L = 0,NBINS-1
                DO K = 0,NBINS-1
                    P2(K,L,J,I) = P_TEMP(K+1,L+1)
                    P2(K,L,I,J) = P_TEMP(K+1,L+1)
                    E2(J,I) = E2(J,I) + P_TEMP(K+1,L+1) * LOG(P_TEMP(K+1,L+1)) / LOG(2.0)
                END DO
            END DO
            M(J,I) = E1(I) + E1(J) - E2(J,I)
        END DO
    END DO
    
    END SUBROUTINE R_MUTUALINFO

END MODULE
