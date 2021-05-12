
!***********************************************************************************************************************************
!**                                       S U B R O U T I N E   T O T A L  D I S S O L V E D  G A S                               **
!***********************************************************************************************************************************

SUBROUTINE TOTAL_DISSOLVED_GAS (P,NSG,N,T,C)
  USE TDGAS; USE STRUCTURES; USE GLOBAL
  INTEGER  :: N, NSG
  REAL     :: T,P,C,SAT,DB,DA,TDG

  SAT = EXP(7.7117-1.31403*(LOG(T+45.93)))*P
  IF (NSG == 0) THEN
    IF (EQSP(N) == 1) THEN
      TDG = AGASSP(N)*.035313*QSP(N)+BGASSP(N)
      IF (TDG > 145.0) TDG = 145.0
      C = SAT
      IF (TDG >= 100.0) C = TDG*SAT/100.0
    ELSE IF (EQSP(N) == 2) THEN
      TDG = AGASSP(N)+BGASSP(N)*EXP(0.035313*QSP(N)*CGASSP(N))
      IF (TDG > 145.0)  TDG = 145.0
      C = SAT
      IF (TDG >= 100.0) C = TDG*SAT/100.0
    ELSE
      DA = SAT-C                                                                               ! MM 5/21/2009 DA: Deficit upstream
      DB = DA/(1.0+0.38*AGASSP(N)*BGASSP(N)*CGASSP(N)*(1.0-0.11*CGASSP(N))*(1.0+0.046*T))      ! DB: deficit downstream
      C  = SAT-DB
    END IF
  ELSE IF (EQGT(N) == 1) THEN
    TDG = AGASGT(N)*0.035313*QGT(N)+BGASGT(N)
    IF (TDG > 145.0) TDG = 145.0
    C = SAT
    IF (TDG >= 100.0) C = TDG*SAT/100.0
  ELSE IF (EQGT(N) == 2) THEN
    TDG = AGASGT(N)+BGASGT(N)*EXP(.035313*QGT(N)*CGASGT(N))
    IF (TDG > 145.0) TDG = 145.0
    C = SAT
    IF (TDG >= 100.0) C = TDG*SAT/100.0
  ELSE
    DA = SAT-C                                                                               ! MM 5/21/2009 DA: Deficit upstream
    DB = DA/(1.0+0.38*AGASSP(N)*BGASSP(N)*CGASSP(N)*(1.0-0.11*CGASSP(N))*(1.0+0.046*T))      ! DB: deficit downstream
    C  = SAT-DB
  END IF
END SUBROUTINE TOTAL_DISSOLVED_GAS
