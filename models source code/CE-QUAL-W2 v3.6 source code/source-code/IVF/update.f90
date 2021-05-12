subroutine update

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  !v3.5
  EXTERNAL RESTART_OUTPUT

!***********************************************************************************************************************************
!*                                       Task 2.7: Variable updates for next timestep                                             **
!***********************************************************************************************************************************

    SZ     = Z
    SKTI   = KTI
    SBKT   = BKT
    VOLUH2 = QUH1  *DLT
    VOLDH2 = QDH1  *DLT
    TSSUH2 = TSSUH1*DLT
    TSSDH2 = TSSDH1*DLT
    DO JW=1,NWB
     KT = KTWB(JW)
      ELKT(JW) = ELWS(DS(BS(JW)))     !EL(KT,DS(BS(JW)))-Z(DS(BS(JW)))*COSA(BS(JW))
      DO JB=BS(JW),BE(JW)
 ! CODE MOVED to after wse computation       ELWS(CUS(JB):DS(JB)+1) = EL(KT,CUS(JB):DS(JB)+1)-Z(CUS(JB):DS(JB)+1)*COSA(JB)
        DO I=US(JB)-1,DS(JB)
          AVHR(KT,I) = H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                               !SW 07/29/04
        END DO
        AVHR(KT,DS(JB)+1)=H1(KT,DS(JB)+1)                                                                              !SW 03/08/05
        DO I=CUS(JB)-1,DS(JB)+1
          FORALL(K=KT:KB(I))           !DO K=KTWB(JW),KB(I)
            QSS(K,I)   = 0.0
            TSS(K,I)   = 0.0
            SU(K,I)    = U(K,I)
            SW(K,I)    = W(K,I)
            T2(K,I)    = T1(K,I)
            SAZ(K,I)   = AZ(K,I)
            H2(K,I)    = H1(K,I)
            BH2(K,I)   = BH1(K,I)
            BHR2(K,I)  = BHR1(K,I)
            AVH2(K,I)  = AVH1(K,I)
            SAVH2(K,I) = AVH2(K,I)
            SAVHR(K,I) = AVHR(K,I)
          END FORALL                                                    !END DO
        END DO
      END DO
    END DO
    IF (CONSTITUENTS) THEN
      CSSUH2 = CSSUH1*DLT
      CSSDH2 = CSSDH1*DLT
      DO JW=1,NWB
        KT = KTWB(JW)
        DO JB=BS(JW),BE(JW)
          DO JC=1,NAC
            DO I=US(JB)-1,DS(JB)+1
              DO K=KT,KB(I)
                DO JE=1,NEP
                  IF (EPIPHYTON_CALC(JW,JE)) EPD(K,I,JE) = MAX(EPD(K,I,JE),0.0)
                END DO
! v3.5 start
                IF (SEDIMENT_CALC(JW))then
                  SED(K,I) = MAX(SED(K,I),0.0)
                  SEDp(K,I) = MAX(SEDp(K,I),0.0)
                  SEDn(K,I) = MAX(SEDn(K,I),0.0)
                  SEDc(K,I) = MAX(SEDc(K,I),0.0)
                end if
! v3.5 end
                CSSB(K,I,CN(JC)) = 0.0
                C1S(K,I,CN(JC))  = C1(K,I,CN(JC))
                C2(K,I,CN(JC))   = MAX(C1(K,I,CN(JC)),0.0)
              END DO
            END DO
          END DO
        END DO
      END DO
    END IF
! v3.5 start
    do jw = 1,nwb
      KT = KTWB(JW)
      do m=1,nmc
        IF (macrophyte_CALC(jw,m)) THEN
          do jb=bs(jw),be(jw)
            DO I=US(JB),DS(JB)
              DO K=kt,kb(i)
                smac(k,i,m)=mac(k,i,m)
                do j=1,kmx
                  smacrc(j,K,I,m)=macrc(j,K,I,m)
                  smacrm(j,K,I,m)=macrm(j,K,I,m)
                end do
              END DO
            END DO
          end do
        end if
      end do
    end do
! v3.5 end
  DO JW = 1,nwb !mlm
    IF (ULTIMATE(JW)) THEN   ! SR 5/15/06
      IF(LAYERCHANGE(JW) == .TRUE.)THEN
      DO K=KTWB(JW),KMX    ! only need to update this for KT - if layer change then update for all variables to be safe                                               !DO K=2,KMX
      RATZ(K,JW)  =  AVH2(K-1,DS(BE(JW)))/AVH2(K,DS(BE(JW)))                                         ! SW 5/20/05
      CURZ1(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K-1,DS(BE(JW)))   ! SW 5/20/05
      CURZ2(K,JW) = -2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))*AVH2(K,DS(BE(JW))))                        ! SW 5/20/05
      CURZ3(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K,DS(BE(JW)))     ! SW 5/20/05
      END DO
      ELSE
      DO K=KTWB(JW),KTWB(JW)+1
      RATZ(K,JW)  =  AVH2(K-1,DS(BE(JW)))/AVH2(K,DS(BE(JW)))                                         ! SW 5/20/05
      CURZ1(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K-1,DS(BE(JW)))   ! SW 5/20/05
      CURZ2(K,JW) = -2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))*AVH2(K,DS(BE(JW))))                        ! SW 5/20/05
      CURZ3(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K,DS(BE(JW)))     ! SW 5/20/05
      END DO
      ENDIF
    END IF
    END DO
    NIT     =  NIT+1
    ELTM    =  ELTM+DLT
    ELTMS   =  ELTMS+DLT
    ELTMF   =  ELTMF+DLT
    JDAY    =  ELTM/DAY
    ELTMJD  =  JDAY-TMSTRT
    END_RUN =  JDAY >= TMEND
    DLT     =  MAX(DLTMIN,DLTF(DLTDP)*CURMAX)
    DLT     =  MIN(DLT,1.1*DLTS)
    DLTAV   = (ELTM-TMSTRT*DAY)/NIT
    IF (DLT <  MINDLT) THEN
      MINDLT = DLTS
      JDMIN  = JDAY
    END IF
    IF (JDAY >= DLTD(DLTDP+1)) DLTDP = DLTDP+1
    IF (DLT  >  DLTMAX(DLTDP)) DLT   = DLTMAX(DLTDP)
    CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)
    IF (INT(JDAY) == JDAYNX) THEN
      JDAYG  = JDAYG+1
      JDAYNX = JDAYNX+1
    END IF
    IF (JDAYG > 300) WINTER = .TRUE.
    IF (JDAYG < 40)  WINTER = .FALSE.
    WRITE (GDCH,'(I3)') GDAY
    CALL GREGORIAN_DATE
    UPDATE_KINETICS = .FALSE.
    IF (MOD(NIT,CUF) == 0) UPDATE_KINETICS = .TRUE.
    return
    end subroutine update
