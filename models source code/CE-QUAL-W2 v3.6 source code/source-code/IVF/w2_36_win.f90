! CE-QUAL-W2 computations
INTEGER(4) FUNCTION CE_QUAL_W2 (DLG)

! CVF specific code

  USE DFLOGM; USE MSCLIB; USE DFWIN, RENAMED => DLT;

 !DEC$ATTRIBUTES STDCALL   :: ce_qual_w2
 !DEC$ATTRIBUTES REFERENCE :: Dlg
 USE IFPORT                 ! to get current working directory
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  !v3.5
  include "omp_lib.h"      ! OPENMP directive to adjust the # of processors
EXTERNAL RESTART_OUTPUT
  TYPE (DIALOG) :: DLG
  INTEGER       :: RESULT
!***********************************************************************************************************************************
!**                                                       Task 1: Inputs                                                          **
!***********************************************************************************************************************************

INTEGER(4) length,istatus
character*255 dirc
!  call omp_set_num_threads(4)   ! set # of processors to NPROC  Moved to INPUT subroutine


call get_command_argument(1,dirc,length,istatus)
dirc=TRIM(dirc)

! if(istatus.ne.0)write(*,*)'get_command_argument failed: status=',istatus

if(length /= 0)then
    istatus=CHDIR(dirc)
    select case(istatus)
      case(2)  ! ENOENT
        write(*,*)'The directory does not exist:',dirc
      case(20)   ! ENOTDIR
        write(*,*)'This is not a directory:', dirc
      case(0)    ! no error
    end select
endif

moddir = FILE$CURDRIVE              !  Get current directory
length = GETDRIVEDIRQQ(moddir)


! Open control file

  OPEN (CON,FILE=CONFN,STATUS='OLD',IOSTAT=I)
  IF (I /= 0) THEN
    TEXT = 'Could not open w2_con.npt'
    GO TO 240
  END IF

CALL INPUT
  RESTART_IN   =  RSIC == '      ON'.OR. RESTART_PUSHED
! Restart data
  IF (RESTART_PUSHED) RSIFN = 'rso.opt'

  JDAY = TMSTRT
  IF (RESTART_IN) THEN
    VERT_PROFILE = .FALSE.
    LONG_PROFILE = .FALSE.
    OPEN  (RSI,FILE=RSIFN,FORM='UNFORMATTED',STATUS='OLD')
    READ  (RSI) NIT,    NV,     KMIN,   IMIN,   NSPRF,  CMBRT,  ZMIN,   IZMIN,  START,  CURRENT
    READ  (RSI) DLTDP,  SNPDP,  TSRDP,  VPLDP,  PRFDP,  CPLDP,  SPRDP,  RSODP,  SCRDP,  FLXDP,  WDODP
    READ  (RSI) JDAY,   YEAR,   ELTM,   ELTMF,  DLT,    DLTAV,  DLTS,   MINDLT, JDMIN,  CURMAX
    READ  (RSI) NXTMSN, NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS, NXTMSC, NXTMSP, NXTMFL, NXTMWD
    READ  (RSI) VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTRB, VOLDT,  VOLWD,  VOLEV,  VOLSBR, VOLTR, VOLSR
    READ  (RSI) TSSEV,  TSSPR,  TSSTR,  TSSDT,  TSSWD,  TSSIN,  TSSOUT, TSSS,   TSSB,   TSSICE
    READ  (RSI) TSSUH,  TSSDH,  TSSUH2, TSSDH2, CSSUH2, CSSDH2, VOLUH2, VOLDH2, QUH1
    READ  (RSI) ESBR,   ETBR,   EBRI
    READ  (RSI) Z,      SZ,     ELWS,   SAVH2,  SAVHR,  H2
    READ  (RSI) KTWB,   KTI,    SKTI,   SBKT
    READ  (RSI) ICE,    ICETH,  CUF,    QSUM
    READ  (RSI) U,      W,      SU,     SW,     AZ,     SAZ,    DLTLIM
    READ  (RSI) T1,     T2,     C1,     C2,     C1S,    EPD,    SED,    KFS,    CSSK
    READ  (RSI) SEDC, SEDN, SEDP, ZOO, CD  ! mlm 10/06
    READ  (RSI) sdkv                       ! cb 11/30/06
    read  (rsi) tke                        ! sw 10/4/07
    CLOSE (RSI)
  END IF
  CE_QUAL_W2 =  1
CALL INIT
  IF (.NOT. RESTART_IN) THEN
    LINE    = CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)
    RESULT  = DLGSET (DLG,STARTING_TIME,TRIM(LINE))                                                   !Display starting time
    RESULT  = DLGSET (DLG,STATUS,'Executing')                                                         !Display execution status
    CURRENT = 0.0
  ELSE
    CALL CPU_TIME (CURRENT)
  END IF

  CALL OUTPUTINIT
    IF (RESTART_IN) THEN
    DO JW=1,NWB
      IF (SCREEN_OUTPUT(JW))CALL SCREEN_UPDATE(DLG)
    ENDDO
    ENDIF
!  ngraph=0
!  DO JH=1,NHY
!    IF (HYDRO_PLOT(JH))ngraph=ngraph+1
!  END DO
!  DO JC=1,NCT
!    IF (CONSTITUENT_PLOT(JC)) ngraph=ngraph+1
!  END DO
!  DO JD=1,NDC
!    IF (DERIVED_PLOT(JD))ngraph=ngraph+1
!  END DO

!  CALL GRAPH(ngraph)
!  DO JH=1,NHY
!    IF (HYDRO_PLOT(JH))       CALL GRAPH_OPEN (JH,HYD(:,:,JH),     HNAME(JH), HYMIN(JH),HYMAX(JH),1.0,       LNAME(JH))
!  END DO
!  DO JC=1,NCT
!    IF (CONSTITUENT_PLOT(JC)) CALL GRAPH_OPEN (JH+JC,C2(:,:,JC),   CNAME(JC), CMIN(JC), CMAX(JC), CMULT(JC), LNAME(JH+JC))
!  END DO
!  DO JD=1,NDC
!    IF (DERIVED_PLOT(JD))     CALL GRAPH_OPEN (JH+JC+JD,CD(:,:,JD),CDNAME(JD),CDMIN(JD),CDMAX(JD),CDMULT(JD),LNAME(JH+JC+JD))
!  END DO
!  INITIALIZE_GRAPH = .FALSE.

  IF (.NOT. RESTART_IN) CALL CPU_TIME (START)

  if (macrophyte_on.and.constituents) call porosity  !v3.5
!!  IF(SELECTC == '      ON')CALL SELECTIVEINIT   ! new subroutine for selecting water temperature target

!***********************************************************************************************************************************
!**                                                   Task 2: Calculations                                                        **
!***********************************************************************************************************************************
  DO WHILE (.NOT. END_RUN.AND. .NOT. STOP_PUSHED)
    IF (JDAY >= NXTVD) CALL READ_INPUT_DATA (NXTVD)
    CALL INTERPOLATE_INPUTS
    DLTTVD = (NXTVD-JDAY)*DAY
    DLT    =  MIN(DLT,DLTTVD+1.0)
    DLTS1  =  DLT
    IF (DLT <= DLTTVD+0.999) THEN
      DLTS = DLT
    ELSE
      KLOC = 1
      ILOC = 1
    END IF

210 continue   ! timestep violation entry point
!! IF(SELECTC == '      ON')CALL SELECTIVE   ! new subroutine for selecting water temperature target
CALL HYDROINOUT

!***********************************************************************************************************************************
!**                                           Task 2.2: Hydrodynamic calculations                                                 **
!***********************************************************************************************************************************

    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)

!***********************************************************************************************************************************
!**                                Task 2.2.1: Boundary concentrations, temperatures, and densities                               **
!***********************************************************************************************************************************

        IUT = IU
        IDT = ID
        IF (UP_FLOW(JB)) THEN
          IF (.NOT. INTERNAL_FLOW(JB)) THEN
            DO K=KT,KB(IU)
                IF (QIND(JB)+QINSUM(JB).GT.0.0) THEN
                  TIN(JB)               = (TINSUM(JB)               *QINSUM(JB)+TIND(JB)          *QIND(JB))/(QIND(JB)+QINSUM(JB))
                  CIN(CN(1:NAC),JB)     =  MAX((CINSUM(CN(1:NAC),JB)*QINSUM(JB)+CIND(CN(1:NAC),JB)*QIND(JB))/(QIND(JB)+QINSUM(JB)),&
                                                0.0)
                  T1(K,IU-1)            =  TIN(JB)
                  T2(K,IU-1)            =  TIN(JB)
                  C1S(K,IU-1,CN(1:NAC)) =  CIN(CN(1:NAC),JB)
                  QIN(JB)               =  QIND(JB)+QINSUM(JB)
                ELSE
                  QIN(JB)               =  0.0
                  TIN(JB)               =  TIND(JB)
                  T1(K,IU-1)            =  TIND(JB)
                  T2(K,IU-1)            =  TIND(JB)
                  C1S(K,IU-1,CN(1:NAC)) =  CIND(CN(1:NAC),JB)
                END IF
            END DO
          ELSE IF (.NOT. DAM_INFLOW(JB)) THEN                                                                          !TC 08/03/04
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              TIN(JB)           = T1(KT,UHS(JB))
              CIN(CN(1:NAC),JB) = MAX(C1S(KT,UHS(JB),CN(1:NAC)),0.0)
              FORALL(K=KT:KB(IU))           !DO K=KT,KB(IU)
                T1(K,IU-1)            = T1(K,UHS(JB))
                T2(K,IU-1)            = T1(K,UHS(JB))
                C1S(K,IU-1,CN(1:NAC)) = C1S(K,UHS(JB),CN(1:NAC))
                C1(K,IU-1,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
                C2(K,IU-1,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
              END FORALL                    !DO
            ELSE
              CALL UPSTREAM_WATERBODY
              TIN(JB)           = T1(KT,IU-1)
              CIN(CN(1:NAC),JB) = MAX(C1(KT,IU-1,CN(1:NAC)),0.0)
            END IF
          ELSE
            TIN(JB)           = TINSUM(JB)
            QIN(JB)           = QINSUM(JB)
            CIN(CN(1:NAC),JB) = MAX(CINSUM(CN(1:NAC),JB),0.0)
            DO K=KT,KB(ID)
              T1(K,IU-1)            = TIN(JB)
              T2(K,IU-1)            = TIN(JB)
              C1S(K,IU-1,CN(1:NAC)) = CIN(CN(1:NAC),JB)
            END DO
           END IF
        END IF
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
              T1(K,ID+1)            = T2(K,ID)
              T2(K,ID+1)            = T2(K,ID)
              C1S(K,ID+1,CN(1:NAC)) = C1S(K,ID,CN(1:NAC))
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          IUT = IU-1
          IF (UH_INTERNAL(JB)) THEN
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              DO K=KT,KB(IUT)
                RHO(K,IUT)           = RHO(K,UHS(JB))
                T1(K,IUT)            = T2(K,UHS(JB))
                T2(K,IUT)            = T2(K,UHS(JB))
                C1S(K,IUT,CN(1:NAC)) = C1S(K,UHS(JB),CN(1:NAC))
                C1(K,IUT,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
                C2(K,IUT,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
              END DO
            ELSE
              CALL UPSTREAM_WATERBODY
            END IF
            DO K=KT,KB(IUT)
              RHO(K,IUT) = DENSITY(T2(K,IUT),MAX(TDS(K,IUT),0.0),MAX(TISS(K,IUT),0.0))
            END DO
          ELSE IF (UH_EXTERNAL(JB)) THEN
            DO K=KT,KB(IUT)
              RHO(K,IUT)           = DENSITY(TUH(K,JB),MAX(TDS(K,IUT),0.0),MAX(TISS(K,IUT),0.0))
              T1(K,IUT)            = TUH(K,JB)
              T2(K,IUT)            = TUH(K,JB)
              C1S(K,IUT,CN(1:NAC)) = CUH(K,CN(1:NAC),JB)
              C1(K,IUT,CN(1:NAC))  = CUH(K,CN(1:NAC),JB)
              C2(K,IUT,CN(1:NAC))  = CUH(K,CN(1:NAC),JB)
            END DO
          END IF
        END IF
        IF (DN_HEAD(JB)) THEN
          IDT = ID+1
          IF (DH_INTERNAL(JB)) THEN
            IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
              DO K=KT,KB(IDT)
                RHO(K,IDT)           = RHO(K,DHS(JB))
                T1(K,IDT)            = T2(K,DHS(JB))
                T2(K,IDT)            = T2(K,DHS(JB))
                C1S(K,IDT,CN(1:NAC)) = C1S(K,DHS(JB),CN(1:NAC))
                C1(K,IDT,CN(1:NAC))  = C1S(K,DHS(JB),CN(1:NAC))
                C2(K,IDT,CN(1:NAC))  = C1S(K,DHS(JB),CN(1:NAC))
              END DO
            ELSE
              CALL DOWNSTREAM_WATERBODY
            END IF
            DO K=KT,KB(ID)
              RHO(K,IDT) = DENSITY(T2(K,IDT),MAX(TDS(K,IDT),0.0),MAX(TISS(K,IDT),0.0))
            END DO
          ELSE IF (DH_EXTERNAL(JB)) THEN
            DO K=KT,KB(IDT)
              RHO(K,IDT)           = DENSITY(TDH(K,JB),MAX(TDS(K,IDT),0.0),MAX(TISS(K,IDT),0.0))
              T1(K,IDT)            = TDH(K,JB)
              T2(K,IDT)            = TDH(K,JB)
              C1S(K,IDT,CN(1:NAC)) = CDH(K,CN(1:NAC),JB)
              C1(K,IDT,CN(1:NAC))  = CDH(K,CN(1:NAC),JB)
              C2(K,IDT,CN(1:NAC))  = CDH(K,CN(1:NAC),JB)
            END DO
          END IF
        END IF

!***********************************************************************************************************************************
!**                                                 Task 2.2.2: Momentum terms                                                    **
!***********************************************************************************************************************************

!****** Density pressures

        DO I=IUT,IDT
          DO K=KT,KB(I)
            P(K,I) = P(K-1,I)+RHO(K,I)*G*H(K,JW)*COSA(JB)
          END DO
        END DO

!****** Horizontal density gradients

        DO I=IUT,IDT-1
          HDG(KT,I) = DLXRHO(I)*(BKT(I)+BKT(I+1))*0.5*H(KT,JW)*(P(KT,I+1)-P(KT,I))
          DO K=KT+1,KBMIN(I)
            HDG(K,I) = DLXRHO(I)*BHR2(K,I)*((P(K-1,I+1)-P(K-1,I))+(P(K,I+1)-P(K,I)))
          END DO
        END DO

!****** Adjusted wind speed and surface wind shear drag coefficient

        DO I=IU-1,ID+1
          WIND10(I) = WIND(JW)*WSC(I)*LOG(10.0/Z0(JW))/LOG(WINDH(JW)/Z0(JW))     ! older  version z0=0.01                      ! SW 11/28/07
          FETCH(I)  = FETCHD(I,JB)
          IF (COS(PHI(JW)-PHI0(I)) < 0.0) FETCH(I) = FETCHU(I,JB)
          IF (FETCH(I) <= 0.0) FETCH(I) = DLX(I)
          IF (FETCH_CALC(JW)) THEN
            ZB        = 0.8*LOG(FETCH(I)*0.5)-1.0718
            WIND10(I) = WIND10(I)*(5.0*ZB+4.6052)/(3.0*ZB+9.2103)
          END IF

          IF(WIND10(I) >= 15.0)THEN                     ! SW 1/19/2008
          CZ(I) = 0.0026
          ELSEIF(WIND10(I) >= 4.0)THEN
          CZ(I) = 0.0005*SQRT(WIND10(I))
          ELSEIF(WIND10(I) >= 0.5)THEN
          CZ(I)= 0.0044*WIND10(I)**(-1.15)
          ELSE
          CZ(I)= 0.01
          ENDIF

  !        CZ(I) = 0.0
  !        IF (WIND10(I) >= 1.0)  CZ(I) = 0.0005*SQRT(WIND10(I))
  !        IF (WIND10(I) >= 4.0) CZ(I) = 0.0005*SQRT(WIND10(I))
  !        IF (WIND10(I) >= 15.0) CZ(I) = 0.0026
        END DO

!****** Longitudinal and lateral surface wind shear and exponential decay

        DO I=IUT,IDT-1
          WSHX(I) = CZ(I)*WIND10(I)**2*RHOA/RHOW*    COS(PHI(JW)-PHI0(I))* ICESW(I)
          WSHY(I) = CZ(I)*WIND10(I)**2*RHOA/RHOW*ABS(SIN(PHI(JW)-PHI0(I)))*ICESW(I)
          WWT     = 0.0
          IF (WIND10(I) /= 0.0) WWT = 6.95E-2*(FETCH(I)**0.233)*WIND10(I)**0.534
          DFC = -8.0*PI*PI/(G*WWT*WWT+NONZERO)
          DO K=KT,KBMIN(I)
            DECAY(K,I) = EXP(MAX(DFC*DEPTHB(K,I),-30.0))
          END DO

!******** Branch inflow lateral shear and friction

          DO JJB=1,NBR
            IF (I == UHS(JJB) .AND. .NOT. INTERNAL_FLOW(JJB)) THEN
              BETABR = (PHI0(I)-PHI0(US(JJB)))
              IF (JJB >= BS(JW) .AND. JJB <= BE(JW)) THEN
                DO K=KT,KBMIN(I)
                  IF (U(K,US(JJB)) < 0.0) THEN
                    UXBR(K,I) = UXBR(K,I)+ABS(U(K,US(JJB)))*COS(BETABR)     *VOLUH2(K,JJB)/(DLT*DLX(I))
                    UYBR(K,I) = UYBR(K,I)              +ABS(SIN(BETABR))*ABS(VOLUH2(K,JJB))/DLT
                  END IF
                END DO
              ELSE
                CALL UPSTREAM_BRANCH
              END IF
            END IF
            IF (I == DHS(JJB)) THEN
              BETABR = (PHI0(I)-PHI0(DS(JJB)))
              IF (I == US(JB) .AND. UHS(JB) /= DS(JJB)) THEN
                IF (JJB >= BS(JW) .AND. JJB <= BE(JW)) THEN
                  DO K=KT,KBMIN(I)
                    IF (U(K,DS(JJB)) >= 0.0) THEN
                      UXBR(K,I) = UXBR(K,I)+U(K,DS(JJB))*   COS(BETABR) *VOLDH2(K,JJB)/(DLT*DLX(I))
                      UYBR(K,I) = UYBR(K,I)            +ABS(SIN(BETABR))*VOLDH2(K,JJB)/DLT
                    END IF
                  END DO
                ELSE
                  CALL DOWNSTREAM_BRANCH
                END IF
              ELSE IF (I /= US(JB)) THEN
                IF (JJB >= BS(JW) .AND. JJB <= BE(JW)) THEN
                  DO K=KT,KBMIN(I)
                    IF (U(K,DS(JJB)) >= 0.0) THEN
                      UXBR(K,I) = UXBR(K,I)+U(K,DS(JJB))*   COS(BETABR) *VOLDH2(K,JJB)/(DLT*DLX(I))
                      UYBR(K,I) = UYBR(K,I)            +ABS(SIN(BETABR))*VOLDH2(K,JJB)/DLT
                    END IF
                  END DO
                ELSE
                  CALL DOWNSTREAM_BRANCH
                END IF
              END IF
            END IF
          END DO
          DO K=KT,KBMIN(I)
            FRICBR(K,I) = (FI(JW)/8.0)*RHO(K,I)*(UYBR(K,I)/(DLX(I)*H2(K,I)))**2
          END DO
        END DO

!****** Vertical eddy viscosities/diffusivities
        FIRSTI(JW) = IUT
		LASTI(JW) = IDT
        DO I=IUT,IDT-1
          CALL CALCULATE_AZ
          IF (KBMIN(I) <= KT+1 .AND. KB(I) > KBMIN(I)) THEN
            AZ(KBMIN(I),I) = AZMIN
            DZ(KBMIN(I),I) = DZMIN
          END IF
        END DO
        IF (AZC(JW) == '     TKE'.OR.AZC(JW) == '    TKE1') THEN
          AZT(:,IDT-1)  = AZ(:,IDT-1)
          DO I=IUT,IDT-2
            DO K=KT,KBMIN(I)
              AZT(K,I)  = 0.5*(AZ(K,I)+AZ(K,I+1))
            END DO
          AZ(KBMIN(I),I) = AZMIN              !SG 10/4/07 SW 10/4/07
          END DO
          AZ(KT:KMX-1,IUT:IDT-1)=AZT(KT:KMX-1,IUT:IDT-1)
        END IF
        DO JWR=1,NIW
        IF (WEIR_CALC) AZ(KTWR(JWR)-1:KBWR(JWR),IWR(1:NIW)) = 0.0
        END DO

!****** Average eddy diffusivities

        IF(AZC(JW) /= '     TKE'.AND.AZC(JW) /= '    TKE1')THEN
        DZ(KT:KB(IDT)-1,IDT) = DZT(KT:KB(IDT)-1,IDT-1)    ! DZT is only used for non-TKE algorithms
        ELSE
        DZ(KT:KB(IDT)-1,IDT) = DZ(KT:KB(IDT)-1,IDT-1)
        ENDIF
        DO I=IUT,IDT-1
          DO K=KT,KB(I)-1
            IF (K >= KBMIN(I)) THEN
              IF (KB(I-1) >= KB(I) .AND. I /= IUT) THEN
                DZ(K,I) = DZ(K,I-1)
              ELSE
                DZ(K,I) = DZMIN
              END IF
            ELSE
              IF(AZC(JW) /= '     TKE'.AND.AZC(JW) /= '    TKE1')THEN
                 if(i == iut)then                             ! SW 10/20/07
                    dz(k,i)=dzt(k,i)
                 else
                    DZ(K,I) = (DZT(K,I)+DZT(K,I-1))*0.5        ! SW 10/20/07  (DZT(K,I)+DZT(K+1,I))*0.5 ! For non-TKE algorithms, average DZ from edges to cell centers
                 endif
              ENDIF
            END IF
          END DO
        END DO

!****** Density inversions

        DO I=IUT,IDT
          DO K=KT,KB(I)-1
            DZQ(K,I) = MIN(1.0E-2,DZ(K,I))                                    !MIN(1.0E-4,DZ(K,I)) No reason to limit DZ in rivers/estuaries-used in ULTIMATE scheme
            IF (RHO(K,I) > RHO(K+1,I)) DZ(K,I) = DZMAX
          END DO
        END DO

!****** Wind, velocity, and bottom shear stresses @ top and bottom of cell

        SB(:,IUT:IDT-1) = 0.0
        DO I=IUT,IDT-1
          ST(KT,I) = WSHX(I)*BR(KTI(I),I)
          DO K=KT+1,KBMIN(I)
            ST(K,I) = WSHX(I)*DECAY(K-1,I)*BR(K,I)
            IF (.NOT. IMPLICIT_VISC(JW)) ST(K,I) = ST(K,I)+AZ(K-1,I)*(BR(K-1,I)+BR(K,I))*0.5*(U(K-1,I)-U(K,I))/((AVH2(K-1,I)       &
                                                   +AVH2(K-1,I+1))*0.5)
          END DO
          GC2 = 0.0
          IF (FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))
! v3.5 start
          hrad=BHR2(KT,I)/(BR(KTI(I),I)-BR(KT+1,I)+2.0*AVHR(KT,I))
          if(macrophyte_on.and.mannings_n(jw))then
            call macrophyte_friction(hrad,fric(i),effric,kt,i)
            gc2=g*effric*effric/hrad**0.33333333
          else if(.not.macrophyte_on.and.mannings_n(jw))then
            gc2=g*fric(i)*fric(i)/hrad**0.33333333
          end if
          IF (ONE_LAYER(I)) THEN
            SB(KT,I) = ST(KT+1,I)+GC2*(BR(KTI(I),I)+2.0*AVHR(KT,I))*U(KT,I)*ABS(U(KT,I))
          ELSE
            SB(KT,I) = GC2*(BR(KTI(I),I)-BR(KT+1,I)+2.0*AVHR(KT,I))*U(KT,I)*ABS(U(KT,I))
            DO K=KT+1,KBMIN(I)-1
              hrad=(BHR2(K,I)/(BR(K,I)-BR(K+1,I)+2.0*H(K,JW)))
              if(macrophyte_on.and.mannings_n(jw))then
                call macrophyte_friction(hrad,fric(i),effric,k,i)
                gc2=g*effric*effric/hrad**0.33333333
              else if(.not.macrophyte_on.and.mannings_n(jw))then
                gc2=g*fric(i)*fric(i)/hrad**0.33333333
              end if
              SB(K,I) = GC2*(BR(K,I)-BR(K+1,I)+2.0*H(K,JW))*U(K,I)*ABS(U(K,I))
            END DO
            IF (KT /= KBMIN(I)) THEN
              hrad=(BHR2(KBMIN(I),I)/(BR(KBMIN(I),I)+2.0*H(KBMIN(I),JW)))
              if(macrophyte_on.and.mannings_n(jw))then
                call macrophyte_friction(hrad,fric(i),effric,kbmin(i),i)
                gc2=g*effric*effric/hrad**0.33333333
              else if(.not.macrophyte_on.and.mannings_n(jw))then
                gc2=g*fric(i)*fric(i)/hrad**0.33333333
              end if
! v3.5 end
              IF (KBMIN(I) /= KB(I)) THEN
                SB(KBMIN(I),I) = GC2*(BR(KBMIN(I),I)-BR(KBMIN(I)+1,I)+2.0*H2(K,I))*U(KBMIN(I),I)*ABS(U(KBMIN(I),I))
              ELSE
                SB(KBMIN(I),I) = GC2*(BR(KBMIN(I),I)+2.0*H2(K,I))*U(KBMIN(I),I)*ABS(U(KBMIN(I),I))
              END IF
            END IF
          END IF
          DO K=KT,KBMIN(I)-1
            SB(K,I) = SB(K,I)+ST(K+1,I)
          END DO
          SB(KBMIN(I),I) = SB(KBMIN(I),I)+WSHX(I)*DECAY(KBMIN(I),I)*(BR(KBMIN(I)-1,I)+BR(KBMIN(I),I))*0.5
        END DO

!****** Horizontal advection of momentum

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            UDR       = (1.0+SIGN(1.0,(U(K,I)+U(K,I+1))*0.5))*0.5
            UDL       = (1.0+SIGN(1.0,(U(K,I)+U(K,I-1))*0.5))*0.5
            ADMX(K,I) = (BH2(K,I+1)*(U(K,I+1)+U(K,I))*0.5*(UDR*U(K,I)+(1.0-UDR)*U(K,I+1))-BH2(K,I)*(U(K,I)+U(K,I-1))               &
                        *0.5*(UDL*U(K,I-1)+(1.0-UDL)*U(K,I)))/DLXR(I)
          END DO
        END DO

!****** Horizontal dispersion of momentum

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            DM(K,I) = AX(JW)*(BH2(K,I+1)*(U(K,I+1)-U(K,I))/DLX(I+1)-BH2(K,I)*(U(K,I)-U(K,I-1))/DLX(I))/DLXR(I)
          END DO
        END DO

!****** Vertical advection of momentum

        DO I=IU,ID-1
          DO K=KT,KB(I)-1
            AB        = (1.0+SIGN(1.0,(W(K,I+1)+W(K,I))*0.5))*0.5
            ADMZ(K,I) = (BR(K,I)+BR(K+1,I))*0.5*(W(K,I+1)+W(K,I))*0.5*(AB*U(K,I)+(1.0-AB)*U(K+1,I))
          END DO
        END DO

!****** Gravity force due to channel slope

        DO I=IU-1,ID
          GRAV(KT,I) = AVHR(KT,I)*(BKT(I)+BKT(I+1))*0.5*G*SINA(JB)                                                     !SW 09/09/04
          DO K=KT+1,KB(I)                                                                                              !SW 09/09/04
            GRAV(K,I) = BHR2(K,I)*G*SINA(JB)
          END DO
        END DO

!***********************************************************************************************************************************
!**                                            Task 2.2.3: Water surface elevation                                                **
!***********************************************************************************************************************************

!****** Tridiagonal coefficients

        BHRHO(IU-1:ID+1) = 0.0; D(IU-1:ID+1) = 0.0; F(IU-1:ID+1) = 0.0
        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            BHRHO(I) = BHRHO(I)+(BH2(K,I+1)/RHO(K,I+1)+BH2(K,I)/RHO(K,I))
          END DO
          DO K=KT,KB(I)
            D(I) = D(I)+(U(K,I)*BHR2(K,I)-U(K,I-1)*BHR2(K,I-1)-QSS(K,I)+(UXBR(K,I)-UXBR(K,I-1))*DLT)
            F(I) = F(I)+(-SB(K,I)+ST(K,I)-ADMX(K,I)+DM(K,I)-HDG(K,I)+GRAV(K,I))
          END DO
        END DO

!****** Boundary tridiagonal coefficients

        D(IU) = 0.0
        DO K=KT,KB(IU)
          D(IU) = D(IU)+(U(K,IU)*BHR2(K,IU)-QSS(K,IU))+UXBR(K,IU)*DLT
        END DO
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
            D(ID) = D(ID)-U(K,ID-1)*BHR2(K,ID-1)-QSS(K,ID)+(UXBR(K,ID)-UXBR(K,ID-1))*DLT+QOUT(K,JB)
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          DO K=KT,KBMIN(IU-1)
            BHRHO(IU-1) = BHRHO(IU-1)+(BH2(K,IU)/RHO(K,IU)+BH2(K,IU-1)/RHO(K,IU-1))
          END DO
          DO K=KT,KB(IU)
            D(IU)   = D(IU)-U(K,IU-1)*BHR2(K,IU-1)
            F(IU-1) = F(IU-1)-(SB(K,IU-1)-ST(K,IU-1)+HDG(K,IU-1)-GRAV(K,IU-1))
          END DO
        END IF
        IF (DN_HEAD(JB)) THEN
          DO K=KT,KBMIN(ID)
            BHRHO(ID) = BHRHO(ID)+(BH2(K,ID+1)/RHO(K,ID+1)+BH2(K,ID)/RHO(K,ID))
          END DO
          DO K=KT,KB(ID)
            D(ID) = D(ID)+(U(K,ID)*BHR2(K,ID)-U(K,ID-1)*BHR2(K,ID-1)-QSS(K,ID))+(UXBR(K,ID)-UXBR(K,ID-1))*DLT
            F(ID) = F(ID)+(-SB(K,ID)+ST(K,ID)-HDG(K,ID)+GRAV(K,ID))
          END DO
        END IF
      END DO
    END DO
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)
        IF (INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                         !TC 08/03/04
          QIN(JB) = 0.0
          DO K=KTWB(JWUH(JB)),KB(UHS(JB))
            QIN(JB) = QIN(JB)+U(K,UHS(JB))*BHR2(K,UHS(JB))
          END DO
        END IF
        IF (UP_FLOW(JB)) D(IU) = D(IU)-QIN(JB)

!****** Boundary surface elevations

        IF (UH_INTERNAL(JB)) THEN
          Z(IU-1)    = ((-EL(KTWB(JWUH(JB)),UHS(JB))+Z(UHS(JB))*COSA(JBUH(JB)))+EL(KT,IU-1)+SINA(JB)*DLXR(IU-1))/COSA(JB)
          ELWS(IU-1) = EL(KT,IU-1)-Z(IU-1)*COSA(JB)
          KTI(IU-1)  = 2
          DO WHILE (EL(KTI(IU-1),IU-1) > ELWS(IU-1))
            KTI(IU-1) = KTI(IU-1)+1
          END DO
          KTI(IU-1) = MAX(KTI(IU-1)-1,2)
        END IF
        IF (UH_EXTERNAL(JB)) Z(IU-1) = (EL(KT,IU-1)-(ELUH(JB)+SINA(JB)*DLX(IU)*0.5))/COSA(JB)
        IF (DH_INTERNAL(JB)) THEN
          Z(ID+1)    = ((-EL(KTWB(JWDH(JB)),DHS(JB))+Z(DHS(JB))*COSA(JBDH(JB)))+EL(KT,ID+1))/COSA(JB)
          ELWS(ID+1) = EL(KT,ID+1)-Z(ID+1)*COSA(JB)
          KTI(ID+1)  = 2
          DO WHILE (EL(KTI(ID+1),ID+1) > ELWS(ID+1))
            KTI(ID+1) = KTI(ID+1)+1
          END DO
          KTI(ID+1) = MAX(KTI(ID+1)-1,2)
          IF (KTI(ID+1) >= KB(ID)) THEN
            Z(ID+1)    = Z(ID)-SLOPE(JB)*DLX(ID)/2.0
            ELWS(ID+1) = EL(KT,ID+1)-Z(ID+1)*COSA(JB)
            KTI(ID+1)  = 2
            DO WHILE (EL(KTI(ID+1),ID+1) > ELWS(ID+1))
              KTI(ID+1) = KTI(ID+1)+1
            END DO
            KTI(ID+1) = MAX(KTI(ID+1)-1,2)
          END IF
        END IF
        IF (DH_EXTERNAL(JB)) Z(ID+1) = (EL(KT,ID+1)-(ELDH(JB)-SINA(JB)*DLX(ID)*0.5))/COSA(JB)

!****** Implicit water surface elevation solution

        FORALL(I=IU:ID)                       !DO I=IU,ID
          A(I) = -RHO(KT,I-1)*G*COSA(JB)*DLT*DLT* BHRHO(I-1)*0.5/DLXR(I-1)
          C(I) = -RHO(KT,I+1)*G*COSA(JB)*DLT*DLT* BHRHO(I)  *0.5/DLXR(I)
          V(I) =  RHO(KT,I)  *G*COSA(JB)*DLT*DLT*(BHRHO(I)  *0.5/DLXR(I)+BHRHO(I-1)*0.5/DLXR(I-1))+DLX(I)*BI(KT,I)
          D(I) =  DLT*(D(I)+DLT*(F(I)-F(I-1)))+DLX(I)*BI(KT,I)*Z(I)
        END FORALL                             !DO
        IF (UP_HEAD(JB)) D(IU) = D(IU)-A(IU)*Z(IU-1)
        IF (DN_HEAD(JB)) D(ID) = D(ID)-C(ID)*Z(ID+1)
        BTA(IU) = V(IU)
        GMA(IU) = D(IU)
        DO I=IU+1,ID
          BTA(I) = V(I)-A(I)/BTA(I-1)*C(I-1)
          GMA(I) = D(I)-A(I)/BTA(I-1)*GMA(I-1)
        END DO
        Z(ID) = GMA(ID)/BTA(ID)
        DO I=ID-1,IU,-1
          Z(I) = (GMA(I)-C(I)*Z(I+1))/BTA(I)
        END DO

!****** Boundary water surface elevations

        IF (UP_FLOW(JB) .AND. .NOT. HEAD_FLOW(JB)) Z(IU-1) = Z(IU)
        IF (UP_FLOW(JB) .AND.       HEAD_FLOW(JB)) Z(IU-1) = (-EL(KTWB(JWUH(JB)),UHS(JB))+Z(UHS(JB))*COSA(JBUH(JB))+EL(KT,IU-1)  &
                                                             +SINA(JBUH(JB))*DLXR(IU-1))/COSA(JBUH(JB))
        IF (DN_FLOW(JB))                           Z(ID+1) = Z(ID)

!****** Updated surface layer and geometry

        IF (.NOT. TRAPEZOIDAL(JW)) THEN                                                                                !SW 07/16/04
          DO I=IU-1,ID+1
            IF (EL(KT,I)-Z(I)*COSA(JB) > EL(KTI(I),I)) THEN
              DO WHILE ( EL(KT,I)-Z(I)*COSA(JB) > EL(KTI(I),I) .AND. KTI(I) /= 2)
                Z(I)   = (EL(KT,I)-EL(KTI(I),I)-(EL(KT,I)-EL(KTI(I),I)-Z(I)*COSA(JB))*(B(KTI(I),I)/B(KTI(I)-1,I)))/COSA(JB)
! v3.5 start
                if(macrophyte_on)then
                  ktip=kti(i)
!c  keeping track if column kti has macrophytes
                  if(ktip.gt.2)kticol(i)=.false.
                end if
! v3.5 end
                KTI(I) =  MAX(KTI(I)-1,2)
              END DO
            ELSE IF (EL(KT,I)-Z(I)*COSA(JB) < EL(KTI(I)+1,I)) THEN
              DO WHILE (EL(KT,I)-Z(I)*COSA(JB) < EL(KTI(I)+1,I))
                Z(I)   = (EL(KT,I)-EL(KTI(I)+1,I)-(EL(KT,I)-EL(KTI(I)+1,I)-Z(I)*COSA(JB))*(B(KTI(I),I)/B(KTI(I)+1,I)))/COSA(JB)
                KTI(I) =  KTI(I)+1
                if(macrophyte_on)kticol(i)=.true.   !v3.5
                IF (KTI(I) >= KB(I)) EXIT
              END DO
            END IF
            BI(KT:KB(I),I) =  B(KT:KB(I),I)
            BI(KT,I)       =  B(KTI(I),I)
            H1(KT,I)       =  H(KT,JW)-Z(I)
            AVH1(KT,I)     = (H1(KT,I)+H1(KT+1,I))*0.5
            IF (KT == KTI(I) .OR. KTI(I) >= KB(I)) THEN
              BH1(KT,I) = B(KT,I)*H1(KT,I)
            ELSE
              BH1(KT,I) = BI(KT,I)*(EL(KT,I)-Z(I)*COSA(JB)-EL(KTI(I)+1,I))/COSA(JB)
            END IF
            DO K=KTI(I)+1,KT
              BH1(KT,I) = BH1(KT,I)+Bnew(K,I)*H(K,JW) !Bnew(K,I)*H(K,JW)   ! sw 1/23/06
            END DO
            BKT(I)    = BH1(KT,I)/H1(KT,I)
            if(kbi(i) < kb(i))bkt(i)=bh1(kt,i)/(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
            VOL(KT,I) = BH1(KT,I)*DLX(I)
          END DO
          DO I=IU-1,ID
            AVHR(KT,I) = H1(KT,I)  +(H1(KT,I+1) -H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                          !SW 07/29/04
            if(kbi(i) < kb(i))avhr(kt,i)=(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))  &
               +(H1(KT,I+1)-(el(kbi(i)+1,i+1)-el(kb(i)+1,i+1)) -H1(KT,I)+(el(kbi(i)+1,i)&
               -el(kb(i)+1,i)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)        ! SW 1/23/06
            BHR1(KT,I) = BH1(KT,I)+(BH1(KT,I+1)-BH1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                          !SW 07/29/04
          END DO
          AVHR(KT,ID+1) = H1(KT,ID+1)
          BHR1(KT,ID+1) = BH1(KT,ID+1)
          DLVOL(JB)        = 0.0
        ELSE                                                                                                           !SW 07/16/04
          DO I=IU-1,ID+1
            BI(KT:KB(I),I) =  B(KT:KB(I),I)
            CALL GRID_AREA2
            H1(KT,I)   =  H(KT,JW)-Z(I)
            AVH1(KT,I) = (H1(KT,I)+H1(KT+1,I))*0.5
            CALL GRID_AREA1 (EL(KT,I)-Z(I),EL(KT+1,I),BH1(KT,I),BI(KT,I))
            BKT(I)    = BH1(KT,I)/H1(KT,I)
            if(kbi(i) < kb(i))bkt(i)=bh1(kt,i)/(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
            VOL(KT,I) = BH1(KT,I)*DLX(I)
          END DO
          DO I=IU-1,ID
            AVHR(KT,I) = H1(KT,I)  +(H1(KT,I+1) -H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                          !SW 07/29/04
            if(kbi(i) < kb(i))avhr(kt,i)=(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i))) &
               +(H1(KT,I+1)-(el(kbi(i)+1,i+1)-el(kb(i)+1,i+1)) -H1(KT,I)+(el(kbi(i)+1,i)&
               -el(kb(i)+1,i)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                                     ! SW 1/23/06
            BHR1(KT,I) = BH1(KT,I)+(BH1(KT,I+1)-BH1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                          !SW 07/29/04
          END DO
          AVHR(KT,ID+1) = H1(KT,ID+1)
          BHR1(KT,ID+1) = BH1(KT,ID+1)
          DLVOL(JB)     = 0.0
        END IF
        ELWS(CUS(JB):DS(JB)+1) = EL(KT,CUS(JB):DS(JB)+1)-Z(CUS(JB):DS(JB)+1)*COSA(JB)
        DO I=IU,ID
          DLVOL(JB) = DLVOL(JB)+(BH1(KT,I)-BH2(KT,I))*DLX(I)
          IF (KT == 2 .AND. H1(KT,I) > H(2,JW) .AND. .NOT. SURFACE_WARNING) THEN
            WRITE (WRN,'(A,I0,A,F0.3)') 'Water surface is above the top of layer 2 in segment ',I,' at day ',JDAY
            WARNING_OPEN    = .TRUE.
            SURFACE_WARNING = .TRUE.
          END IF
        END DO

! v3.5
        if(macrophyte_on)then
!c  if depth in kti layer becomes greater than threshold, setting
!c      macrophyte conc. in kti column to initial conc.
          do i=iu,id
            depkti=ELWS(i)-el(kti(i)+1,i)

!******* macrophytes, setting conc. of macrophytes in new columns to
!********* initial concentration if column depth is greater than 'thrkti'
            if(.not.kticol(i).and.depkti.ge.thrkti)then
              kticol(i)=.true.
              jt=kti(i)
              mact(jt,kt,i)=0.0
              do m=1,nmc
                macrc(jt,kt,i,m)=macwbci(jw,m)
                colb=el(kti(i)+1,i)
                coldep=ELWS(i)-colb
                macrm(jt,kt,i,m)=macwbci(jw,m)*coldep*cw(jt,i)*dlx(i)
                mact(jt,kt,i)=mact(jt,kt,i)+macwbci(jw,m)
                maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(jt,kt,i,m)
              end do
            end if

!****** macrophytes, when column depth is less than 'thrkti', zeroing out conc.
            if(kticol(i).and.depkti.lt.thrkti)then
              kticol(i)=.false.
              jt=kti(i)
              mact(jt,kt,i)=0.0
              do m=1,nmc
                maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(jt,kt,i,m)
                macrc(jt,kt,i,m)=0.0
                macrm(jt,kt,i,m)=0.0
              end do
            end if
          end do
        end if
! v3.5 end

!***********************************************************************************************************************************
!**                                             Task 2.2.4: Longitudinal velocities                                               **
!***********************************************************************************************************************************

        IUT = IU
        IDT = ID
        IF (UP_HEAD(JB)) IUT = IU-1
        IF (DN_HEAD(JB)) IDT = ID+1

!****** Pressures

        DO I=IUT,IDT
          DO K=KT,KB(I)
            P(K,I) = P(K-1,I)+RHO(K,I)*G*H1(K,I)*COSA(JB)
          END DO
        END DO

!****** Horizontal pressure gradients

        DO I=IUT,IDT-1
          HPG(KT,I) = DLXRHO(I)*(BKT(I)+BKT(I+1))*0.5*(H1(KT,I+1)*P(KT,I+1)-H1(KT,I)*P(KT,I))
          DO K=KT+1,KBMIN(I)
            HPG(K,I) = DLXRHO(I)*BHR2(K,I)*((P(K-1,I+1)-P(K-1,I))+(P(K,I+1)-P(K,I)))
          END DO
        END DO

!****** Boundary horizontal velocities

        IF (UP_FLOW(JB)) THEN
          IF (.NOT. HEAD_FLOW(JB)) THEN
            QINF(:,JB) = 0.0
            IF (PLACE_QIN(JW)) THEN

!************ Inflow layer

              K     = KT
              SSTOT = 0.0
              DO JC=NSSS,NSSE
                SSTOT = SSTOT+CIN(JC,JB)
              END DO
              RHOIN = DENSITY(TIN(JB),MAX(CIN(1,JB),0.0),MAX(SSTOT,0.0))
              DO WHILE (RHOIN > RHO(K,IU) .AND. K < KB(IU))
                K = K+1
              END DO
              KTQIN(JB) = K
              KBQIN(JB) = K

!************ Layer inflows

              VQIN  =  QIN(JB)*DLT
              VQINI =  VQIN
              QINFR =  1.0
              INCR  = -1
              DO WHILE (QINFR > 0.0)
                V1 = VOL(K,IU)
                IF (K <= KB(IU)) THEN
                  IF (VQIN > 0.5*V1) THEN
                    QINF(K,JB) = 0.5*V1/VQINI
                    QINFR      = QINFR-QINF(K,JB)
                    VQIN       = VQIN-QINF(K,JB)*VQINI
                    IF (K == KT) THEN
                      K    = KBQIN(JB)
                      INCR = 1
                    END IF
                  ELSE
                    QINF(K,JB) = QINFR
                    QINFR      = 0.0
                  END IF
                  IF (INCR < 0) KTQIN(JB) = K
                  IF (INCR > 0) KBQIN(JB) = MIN(KB(IU),K)
                  K = K+INCR
                ELSE
                  QINF(KT,JB) = QINF(KT,JB)+QINFR
                  QINFR       = 0.0
                END IF
              END DO
            ELSE
              KTQIN(JB) = KT
              KBQIN(JB) = KB(IU)
              BHSUM     = 0.0
              DO K=KT,KB(IU)
                BHSUM = BHSUM+BH1(K,IU)
              END DO
              DO K=KT,KB(IU)
                QINF(K,JB) = BH1(K,IU)/BHSUM
              END DO
            END IF
            DO K=KT,KB(IU)
              U(K,IU-1) = QINF(K,JB)*QIN(JB)/BHR1(K,IU-1)
            END DO
          ELSE
            KTQIN(JB) = KT
            KBQIN(JB) = KB(IU)
            IF (JBUH(JB) <= BE(JW) .AND. JBUH(JB) >= BS(JW)) THEN
              DO K=KT,KB(IU)
                U(K,IU-1) = U(K,UHS(JB))*BHR1(K,UHS(JB))/BHR1(K,IU-1)
              END DO
            ELSE
              CALL UPSTREAM_VELOCITY
            END IF
          END IF
        END IF
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
            U(K,ID) = QOUT(K,JB)/BHR1(K,ID)
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          DO K=KT,KB(IU-1)
            U(K,IU-1) = (BHR2(K,IU-1)*U(K,IU-1)+DLT*(-SB(K,IU-1)+ST(K,IU-1)-HPG(K,IU-1)+GRAV(K,IU-1)))/BHR1(K,IU-1)
          END DO
        END IF
        IF (DN_HEAD(JB)) THEN
          DO K=KT,KB(ID+1)
            U(K,ID) = (BHR2(K,ID)*U(K,ID)+DLT*(-SB(K,ID)+ST(K,ID)-HPG(K,ID)+GRAV(K,ID)))/BHR1(K,ID)
          END DO
        END IF

!****** Horizontal velocities

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            U(K,I) = (BHR2(K,I)*U(K,I))/BHR1(K,I)+(DLT*(-SB(K,I)+ST(K,I)-ADMZ(K,I)+ADMZ(K-1,I)-ADMX(K,I)+DM(K,I)-HPG(K,I)+GRAV(K,I)&
                     +UXBR(K,I)/H2(K,I)))/BHR1(K,I)
            IF (INTERNAL_WEIR(K,I)) U(K,I) = 0.0
          END DO
        END DO

!****** Implicit vertical eddy viscosity

        IF (IMPLICIT_VISC(JW)) THEN
          AT = 0.0; CT = 0.0; VT = 0.0; DT = 0.0
          DO I=IUT,IDT-1
            DO K=KT,KBMIN(I)            !KB(I)  SW 10/7/07
              AT(K,I) = -DLT/BHR1(K,I)*AZ(K-1,I)*(BHR1(K-1,I)/AVHR(K-1,I)+BR(K,I))  /(AVH1(K-1,I)+AVH1(K-1,I+1))
              CT(K,I) = -DLT/BHR1(K,I)*AZ(K,I)  *(BHR1(K,I)  /AVHR(K,I)  +BR(K+1,I))/(AVH1(K,I)  +AVH1(K,I+1))
              VT(K,I) =  1.0-AT(K,I)-CT(K,I)
              DT(K,I) =  U(K,I)
            END DO
            CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KBMIN(I),KMX,U(:,I))
          END DO
        END IF

!****** Corrected horizontal velocities

        IF (UP_HEAD(JB)) THEN
          IS    =  ID
          IE    =  IU-1
          INCR  = -1
          Q(IS) =  0.0
          DO K=KT,KB(ID)
            Q(IS) = Q(IS)+U(K,IS)*BHR1(K,IS)
          END DO
          QSSUM(IS) = 0.0
          DO K=KT,KB(IS)
            QSSUM(IS) = QSSUM(IS)+QSS(K,IS)
          END DO
        ELSE
          IS   = IU-1
          IE   = ID
          INCR = 1
          IF (DN_FLOW(JB)) IE = ID-1
          Q(IS) = 0.0
          DO K=KT,KB(IU)
            Q(IS) = Q(IS)+U(K,IS)*BHR1(K,IS)
          END DO
        END IF
        QC(IS) = Q(IS)
        DO I=IS+INCR,IE,INCR
          QSSUM(I) = 0.0
          DO K=KT,KB(I)
            QSSUM(I) = QSSUM(I)+QSS(K,I)
          END DO
          BHRSUM = 0.0
          Q(I)   = 0.0
          DO K=KT,KBMIN(I)
            IF (.NOT. INTERNAL_WEIR(K,I)) THEN
              BHRSUM = BHRSUM+BHR1(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR1(K,I)
            END IF
          END DO
          IF (UP_HEAD(JB)) THEN
            QC(I) = QC(I+1)+(BH1(KT,I+1)-BH2(KT,I+1))*DLX(I+1)/DLT-QSSUM(I+1)
          ELSE
            QC(I) = QC(I-1)-(BH1(KT,I)  -BH2(KT,I))  *DLX(I)  /DLT+QSSUM(I)
          END IF
          DO K=KT,KBMIN(I)
            IF (INTERNAL_WEIR(K,I)) THEN
              U(K,I) = 0.0
            ELSE
              U(K,I) =  U(K,I)+(QC(I)-Q(I))/BHRSUM
              IF (Q(I) /= 0.0) QERR(I) = (Q(I)-QC(I))/Q(I)*100.0
            END IF
          END DO
        END DO

!****** Head boundary flows

        IF (UP_HEAD(JB)) QUH1(KT:KB(IU-1),JB) = U(KT:KB(IU-1),IU-1)*BHR1(KT:KB(IU-1),IU-1)
        IF (DN_HEAD(JB)) QDH1(KT:KB(ID+1),JB) = U(KT:KB(ID+1),ID)  *BHR1(KT:KB(ID+1),ID)

!***********************************************************************************************************************************
!**                                              Task 2.2.5: Vertical velocities                                                  **
!***********************************************************************************************************************************

        DO I=IU,ID
          DO K=KB(I)-1,KT,-1
            WT1    =  W(K+1,I)*BB(K+1,I)
            WT2    = (BHR(K+1,I)*U(K+1,I)-BHR(K+1,I-1)*U(K+1,I-1)-QSS(K+1,I))/DLX(I)
            W(K,I) = (WT1+WT2)/BB(K,I)
          END DO
        END DO
      END DO
    END DO

!***********************************************************************************************************************************
!**                                                  Task 2.2.6: Autostepping                                                     **
!***********************************************************************************************************************************

    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB),DS(JB)
          IF (H1(KT,I) < 0.0) THEN
            WRITE (WRN,'(A,F0.3,A,I0/4(A,F0.3))') 'Computational warning at Julian day = ',JDAY,' at segment ',I,'timestep = ',DLT,&
                                                  ' water surface deviation [Z] = ',Z(I),' m  layer thickness = ',H1(KT,I),' m'
            WARNING_OPEN = .TRUE.
            IF (DLT > DLTMIN) THEN
              WRITE (WRN,'(A,I0/2(A,F0.3),A,I0)') 'Negative surface layer thickness in segment ',I,'  time step reduced to ',  &
                                                   DLTMIN,' s on day ',JDAY,' at iteration ',NIT
              WARNING_OPEN = .TRUE.
              CURMAX       =  DLTMIN
              GO TO 220
            ELSE
              WRITE (W2ERR,'(A,F0.3/A,I0)') 'Unstable water surface elevation on day ',JDAY,'negative surface layer thickness '//  &
                                            'using minimum timestep at iteration ',NIT
              WRITE (W2ERR,'(A)') 'Segment, Surface layer thickness, m'
              DO II=MAX(CUS(JB),I-2),MIN(DS(JB),I+2)
                WRITE (W2ERR,'(T6,I3,T21,F10.2)') II,H1(KT,I)
              END DO
              TEXT = 'Runtime error - see w2.err'
              ERROR_OPEN = .TRUE.
              GO TO 230
            END IF
          END IF
        END DO
        DO I=CUS(JB),DS(JB)
          IF (VISCOSITY_LIMIT(JW)) TAU1   = 2.0*MAX(AX(JW),DXI(JW))/(DLX(I)*DLX(I))
          IF (CELERITY_LIMIT(JW))  CELRTY = SQRT((ABS(RHO(KB(I),I)-RHO(KT,I)))/1000.0*G*DEPTHB(KBI(I),I)*0.5)               ! SW 1/23/06
          DO K=KT,KB(I)
            IF (VISCOSITY_LIMIT(JW) .AND. .NOT. IMPLICIT_VISC(JW)) TAU2 = 2.0*AZ(K,I)/(H1(K,JW)*H1(K,JW))
            QTOT(K,I) = (ABS(U(K,I))*BHR1(K,I)+ABS(U(K,I-1))*BHR1(K,I-1)+(ABS(W(K,I))*BB(K,I)+ABS(W(K-1,I))*BB(K-1,I))*DLX(I)      &
                        +DLX(I)*ABS(BH2(K,I)-BH1(K,I))/DLT+ABS(QSS(K,I)))*0.5
            DLTCAL    = 1.0/((QTOT(K,I)/BH1(K,I)+CELRTY)/DLX(I)+TAU1+TAU2+NONZERO)
            IF (DLTCAL < CURMAX) THEN
              KLOC   = K
              ILOC   = I
              CURMAX = DLTCAL
              IF (DLTF(DLTDP)*CURMAX < MINDLT) THEN
                KMIN = K
                IMIN = I
              END IF
            END IF
          END DO
        END DO
      END DO
    END DO

!** Restore timestep dependent variables and restart calculations

220 CONTINUE
    IF (CURMAX < DLT .AND. DLT > DLTMIN) THEN
      DLT = DLTF(DLTDP)*CURMAX
      IF (DLT <= DLTMIN) THEN
        WRITE (WRN,'(A,F0.3/A,F0.3,A)') 'Computational warning at Julian day = ',JDAY,' timestep = ',DLT,' sec'
        WARNING_OPEN = .TRUE.
        DLT          =  DLTMIN
      END IF
      NV        = NV+1
      Z         = SZ
      U         = SU
      W         = SW
      AZ        = SAZ
      AVH2      = SAVH2
      AVHR      = SAVHR
      KTI       = SKTI
      BKT       = SBKT
      QSS       = 0.0
      SB        = 0.0
      DLTS      = DLT

        do jw=1,nwb                                                                            ! SW 8/25/05
        do jb=bs(jw),be(jw)
        do i=us(jb)-1,ds(jb)+1
            VOL(KTWB(JW),I) = BH2(KTWB(JW),I)*DLX(I)
            BI(KTWB(JW),I) = B(KTI(I),I)
        end do
        end do
        end do


      CURMAX    = DLTMAX(DLTDP)/DLTF(DLTDP)
      IF (PIPES) THEN
        YS   = YSS
        VS   = VSS
        VST  = VSTS
        YST  = YSTS
        DTP  = DTPS
        QOLD = QOLDS
      END IF
! v3.5 start
!********** Macrophytes
      do jw=1,nwb
        do m=1,nmc
          IF (macrophyte_CALC(jw,m)) THEN
            KT = KTWB(JW)
              DO JB=BS(Jw),BE(Jw)
                DO I=CUS(JB),DS(JB)
                  DO K=kt,kb(i)
                    mac(k,i,m)=smac(k,i,m)
                    if(kticol(i))then
                      jt=kti(i)
                    else
                      jt=kti(i)+1
                    end if
                    je=kb(i)
                    do j=jt,je
                      macrc(j,K,I,m)=smacrc(j,K,I,m)
                      macrm(j,K,I,m)=smacrm(j,K,I,m)
                    end do
                  END DO
                END DO
             end do
          end if
        end do
      end do
! v3.5 end
      GO TO 210
    END IF
    DLTLIM(KMIN,IMIN) = DLTLIM(KMIN,IMIN)+1.0

!** Layer bottom and middle depths

    DO JW=1,NWB
      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB)-1,DS(JB)
          DEPTHB(KTWB(JW),I) = H1(KTWB(JW),I)
             if(kbi(i) < kb(i))depthb(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
          DEPTHM(KTWB(JW),I) = H1(KTWB(JW),I)*0.5
             if(kbi(i) < kb(i))depthm(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))*0.5    ! SW 1/23/06
          DO K=KTWB(JW)+1,KMX
            DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
            DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
          END DO
        END DO
      END DO
    END DO

CALL temperature

IF (CONSTITUENTS) CALL wqconstituents

CALL LAYERADDSUB
if(error_open)go to 230

CALL BALANCES

CALL UPDATE

CALL OUTPUTA
!**** Screen output
DO JW=1,NWB
      IF (SCREEN_OUTPUT(JW)) THEN
        IF (JDAY >= NXTMSC(JW) .OR. JDAY >= SCRD(SCRDP(JW)+1,JW)) THEN
          IF (JDAY >= SCRD(SCRDP(JW)+1,JW)) THEN
            SCRDP(JW)  = SCRDP(JW)+1
            NXTMSC(JW) = SCRD(SCRDP(JW),JW)
          END IF
          KT         = KTWB(JW)
          NXTMSC(JW) = NXTMSC(JW)+SCRF(SCRDP(JW),JW)
          CALL SCREEN_UPDATE (DLG)
          CALL DATE_AND_TIME (CDATE,CCTIME)
 !         DO JH=1,NHY
 !           IF (HYDRO_PLOT(JH))       CALL GRAPH_UPDATE (JH,HYD(:,:,JH),     HNAME(JH), HYMIN(JH),1.0,       LNAME(JH))
 !         END DO
 !         DO JC=1,NCT
 !           IF (CONSTITUENT_PLOT(JC)) CALL GRAPH_UPDATE (JH+JC,C2(:,:,JC),   CNAME(JC), CMIN(JC), CMULT(JC), LNAME(JH+JC))
 !         END DO
 !         DO JD=1,NDC
 !           IF (DERIVED_PLOT(JD))     CALL GRAPH_UPDATE (JH+JC+JD,CD(:,:,JD),CDNAME(JD),CDMIN(JD),CDMULT(JD),LNAME(JH+JC+JD))
 !         END DO
        END IF
      END IF
 END DO
END DO    ! END OF DO WHILE LOOP

230 CONTINUE
  IF (STOP_PUSHED) THEN
    TEXT  = 'Execution stopped at '//CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)//' on '//CDATE(5:6)//'/'//CDATE(7:8)//'/'        &
                                   //CDATE(3:4)
    CALL RESTART_OUTPUT ('rso.opt')
  END IF

CALL ENDSIMULATION

240 CONTINUE
!  CALL DEALLOCATE_GRAPH

  if(closec=='      ON' .and. END_RUN)then
  CALL EXITDIALOG(DLG,TEXT)
  else
  CALL STOP_W2 (DLG,TEXT)
  return
END FUNCTION CE_QUAL_W2
