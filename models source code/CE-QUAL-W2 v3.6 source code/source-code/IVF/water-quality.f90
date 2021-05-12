
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   K I N E T I C S                                              **
!***********************************************************************************************************************************

SUBROUTINE KINETICS
  USE SCREENC; USE GLOBAL; USE KINETIC; USE GEOMC; USE TVDC; USE LOGICC; USE SURFHE
  use macrophytec; use zooplanktonc; use MAIN, only:EPIPHYTON_CALC, BOD_CALC, ALG_CALC

! Type declarations

  REAL                                :: LAM1,   LAM2,   NH4PR,  NO3PR,  LIMIT,  LIGHT,  L, L0, L1
  REAL                                :: KW,     INCR,   OH,     K1,     K2, bicart
  real                                :: cart,alkt,t1k,s2,sqrs2,dh1,dh2,h2co3t,co3t,pht,f,hion,hco3t
  real                                :: ltcoefm, lavg,  macext, tmac
  real                                :: fetch, u2, coef1,coef2,coef3,coef4,hs,ts,coef,uorb,tau
  real                                :: epsilon, cbodset, dosat,o2ex,co2ex,sedsi,sedem, sedso,sedsip
  real                                :: sedsop,sedson,sedsoc,sedsic,sedsidk,sedsum,SEDSUMK,XDUM
  real                                :: blim, sedsin, colb,coldep,bmass,bmasstest,cvol
  real                                :: algex, ssext, totmac, zooext, totss0, fdpo4, zminfac, ssr
  real                                :: zgztot,cbodc,cbodn,cbodp,bodtot
  real                                :: algp,algn,zoop,zoon,tpss,xx   ! SW 4/5/09
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: OMTRM,  SODTRM, NH4TRM, NO3TRM, BIBH2
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: DOM,    POM,    PO4BOD, NH4BOD, TICBOD
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: lam2m  ! v3.5
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ATRM,   ATRMR,  ATRMF
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ETRM,   ETRMR,  ETRMF
  integer                             :: K, JA, JE, M, js, jt, jj, jjz, jg, jcb, jbod, llm,j,jd
  integer                             :: mi,jaf,n,iter,ibod
  SAVE

! Allocation declarations

  ALLOCATE (OMTRM(KMX,IMX),    SODTRM(KMX,IMX),    NH4TRM(KMX,IMX),    NO3TRM(KMX,IMX), DOM(KMX,IMX), POM(KMX,IMX))
  ALLOCATE (PO4BOD(KMX,IMX),   NH4BOD(KMX,IMX),    TICBOD(KMX,IMX))
  ALLOCATE (ATRM(KMX,IMX,NAL), ATRMR(KMX,IMX,NAL), ATRMF(KMX,IMX,NAL))
  ALLOCATE (ETRM(KMX,IMX,NEP), ETRMR(KMX,IMX,NEP), ETRMF(KMX,IMX,NEP))
  ALLOCATE (lam2m(KMX,kmx),    BIBH2(KMX,IMX))
RETURN

!***********************************************************************************************************************************
!**                                      T E M P E R A T U R E  R A T E  M U L T I P L I E R S                                    **
!***********************************************************************************************************************************

ENTRY TEMPERATURE_RATES
  DO I=IU,ID
    DO K=KT,KB(I)
      LAM1        = FR(T1(K,I),NH4T1(JW),NH4T2(JW),NH4K1(JW),NH4K2(JW))
      NH4TRM(K,I) = LAM1/(1.0+LAM1-NH4K1(JW))
      LAM1        = FR(T1(K,I),NO3T1(JW),NO3T2(JW),NO3K1(JW),NO3K2(JW))
      NO3TRM(K,I) = LAM1/(1.0+LAM1-NO3K1(JW))
      LAM1        = FR(T1(K,I),OMT1(JW),OMT2(JW),OMK1(JW),OMK2(JW))
      OMTRM(K,I)  = LAM1/(1.0+LAM1-OMK1(JW))
      LAM1        = FR(T1(K,I),SODT1(JW),SODT2(JW),SODK1(JW),SODK2(JW))
      SODTRM(K,I) = LAM1/(1.0+LAM1-SODK1(JW))
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        LAM1          = FR(T1(K,I),AT1(JA),AT2(JA),AK1(JA),AK2(JA))
        LAM2          = FF(T1(K,I),AT3(JA),AT4(JA),AK3(JA),AK4(JA))
        ATRMR(K,I,JA) = LAM1/(1.0+LAM1-AK1(JA))
        ATRMF(K,I,JA) = LAM2/(1.0+LAM2-AK4(JA))
        ATRM(K,I,JA)  = ATRMR(K,I,JA)*ATRMF(K,I,JA)
        ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        LAM1          = FR(T1(K,I),ET1(JE),ET2(JE),EK1(JE),EK2(JE))
        LAM2          = FF(T1(K,I),ET3(JE),ET4(JE),EK3(JE),EK4(JE))
        ETRMR(K,I,JE) = LAM1/(1.0+LAM1-EK1(JE))
        ETRMF(K,I,JE) = LAM2/(1.0+LAM2-EK4(JE))
        ETRM(K,I,JE)  = ETRMR(K,I,JE)*ETRMF(K,I,JE)
        endif
      END DO
      do m=1,nmc
      if(macrophyte_calc(jw,m))then
        LAM1    = FR(T1(K,I),mT1(m),mT2(m),mK1(m),mk2(m))
        LAM2    = FF(T1(K,I),mT3(m),mT4(m),mK3(m),mk4(m))
        MACTRMR(K,I,m) = LAM1/(1.0+LAM1-mK1(m))
        MACTRMF(K,I,m) = LAM2/(1.0+LAM2-mK4(m))
        MACTRM(K,I,m)  = macTRMR(K,I,m)*macTRMF(K,I,m)
      endif
      end do
      if(zooplankton_calc)then
	    DO JZ = 1, NZP
          LAM1       = FR(T1(K,I),zt1(jz),zt2(jz),zk1(jz),zk2(jz))
          LAM2       = FF(T1(K,I),zt3(jz),zt4(jz),zk3(jz),zk4(jz))
          zoormr(k,i,jz)= lam1/(1.+lam1-zk1(jz))
          zoormf(k,i,jz)= lam2/(1.+lam2-zk4(jz))
          zoorm(k,i,jz) = zoormr(k,i,jz)*zoormf(k,i,jz)
        END DO
	  end if
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                 K I N E T I C   R A T E S                                                     **
!***********************************************************************************************************************************

ENTRY KINETIC_RATES

! Decay rates
!$OMP PARALLEL DO
  DO I=IU,ID
    DO K=KT,KB(I)
      DO1(K,I)          = O2(K,I)/(O2(K,I)+KDO)                  !V3.5
      DO2(K,I)          = 1.0 - DO1(K,I)                         !O2(K,I)/(O2(K,I)+KDO)
      DO3(K,I)          = (1.0+SIGN(1.0,O2(K,I)-1.E-10)) *0.5
      SEDD(K,I)         =   SODTRM(K,I) *SDKV(K,I)   *SED(K,I) *DO3(K,I)   !CB 10/22/06
      SEDDP(K,I)         =  SODTRM(K,I) *SDKV(K,I)   *SEDP(K,I) *DO3(K,I)
      SEDDN(K,I)         =  SODTRM(K,I) *SDKV(K,I)   *SEDN(K,I) *DO3(K,I)
      SEDDC(K,I)         =  SODTRM(K,I) *SDKV(K,I)   *SEDC(K,I) *DO3(K,I)
      SEDBR(K,I)         =  SEDB(JW)    *SED(K,I)                           !CB 11/30/06
      SEDBRP(K,I)        =  SEDB(JW)    *SEDP(K,I)                          !CB 11/30/06
      SEDBRN(K,I)        =  SEDB(JW)    *SEDN(K,I)                          !CB 11/30/06
      SEDBRC(K,I)        =  SEDB(JW)    *SEDC(K,I)                          !CB 11/30/06
      NH4D(K,I)         =  NH4TRM(K,I) *NH4DK(JW) *NH4(K,I) *DO1(K,I)
      NO3D(K,I)         =  NO3TRM(K,I) *NO3DK(JW) *NO3(K,I) *DO2(K,I)
      LDOMD(K,I)        =  OMTRM(K,I)  *LDOMDK(JW)*LDOM(K,I)*DO3(K,I)
      RDOMD(K,I)        =  OMTRM(K,I)  *RDOMDK(JW)*RDOM(K,I)*DO3(K,I)
      LPOMD(K,I)        =  OMTRM(K,I)  *LPOMDK(JW)*LPOM(K,I)*DO3(K,I)
      RPOMD(K,I)        =  OMTRM(K,I)  *RPOMDK(JW)*RPOM(K,I)*DO3(K,I)
      LRDOMD(K,I)       =  OMTRM(K,I)  *LRDDK(JW) *LDOM(K,I)*DO3(K,I)
      LRPOMD(K,I)       =  OMTRM(K,I)  *LRPDK(JW) *LPOM(K,I)*DO3(K,I)
      CBODD(K,I,1:NBOD) =  KBOD(1:NBOD)*TBOD(1:NBOD)**(T1(K,I)-20.0)*DO3(K,I)
        IF(K == KB(I))THEN     ! SW 4/18/07
	  SODD(K,I)         =  SOD(I)/BH2(K,I)*SODTRM(K,I)*BI(K,I)
	    ELSE
      SODD(K,I)         =  SOD(I)/BH2(K,I)*SODTRM(K,I)*(BI(K,I)-BI(K+1,I))
	    ENDIF

! Inorganic suspended solids settling rates - P adsorption onto SS and Fe
      FPSS(K,I) = PARTP(JW)         /(PARTP(JW)*TISS(K,I)+PARTP(JW)*FE(K,I)*DO1(K,I)+1.0)
      FPFE(K,I) = PARTP(JW)*FE(K,I) /(PARTP(JW)*TISS(K,I)+PARTP(JW)*FE(K,I)*DO1(K,I)+1.0)
      SSSI(K,I) = SSSO(K-1,I)
      TOTSS0    = 0.0
      DO JS=1,NSS
        TOTSS0 = TOTSS0+SSS(JS)*FPSS(K,I)*SS(K,I,JS)
      END DO
      SSSO(K,I) = (TOTSS0+FES(JW)*FPFE(K,I))*BI(K,I)/BH2(K,I)*DO1(K,I)                ! SW 11/7/07
      FPSS(K,I) =  FPSS(K,I)*TISS(K,I)

! OM stoichiometry
        ORGPLD(K,I)=0.0
        ORGPRD(K,I)=0.0
        ORGPLP(K,I)=0.0
        ORGPRP(K,I)=0.0
        ORGNLD(K,I)=0.0
        ORGNRD(K,I)=0.0
        ORGNLP(K,I)=0.0
        ORGNRP(K,I)=0.0
        IF(CAC(NLDOMP) == '      ON')THEN
          IF(LDOM(K,I).GT.0.0)THEN
          ORGPLD(K,I)=LDOMP(K,I)/LDOM(K,I)
          ELSE
          ORGPLD(K,I)=ORGP(JW)
          ENDIF
        ELSE
          ORGPLD(K,I)=ORGP(JW)
        END IF
        IF(CAC(NRDOMP) == '      ON')THEN
          IF(RDOM(K,I).GT.0.0)THEN
          ORGPRD(K,I)=RDOMP(K,I)/RDOM(K,I)
          ELSE
          ORGPRD(K,I)=ORGP(JW)
          ENDIF
        ELSE
          ORGPRD(K,I)=ORGP(JW)
        END IF
        IF(CAC(NLPOMP) == '      ON')THEN
          IF(LPOM(K,I).GT.0.0)THEN
          ORGPLP(K,I)=LPOMP(K,I)/LPOM(K,I)
          ELSE
          ORGPLP(K,I)=ORGP(JW)
          ENDIF
        ELSE
          ORGPLP(K,I)=ORGP(JW)
        END IF
        IF(CAC(NRPOMP) == '      ON')THEN
          IF(RPOM(K,I).GT.0.0)THEN
          ORGPRP(K,I)=RPOMP(K,I)/RPOM(K,I)
          ELSE
          ORGPRP(K,I)=ORGP(JW)
          ENDIF
        ELSE
          ORGPRP(K,I)=ORGP(JW)
        END IF
        IF(CAC(NLDOMN) == '      ON')THEN
          IF(LDOM(K,I).GT.0.0)THEN
          ORGNLD(K,I)=LDOMN(K,I)/LDOM(K,I)
          ELSE
          ORGNLD(K,I)=ORGN(JW)
          ENDIF
        ELSE
          ORGNLD(K,I)=ORGN(JW)
        END IF
        IF(CAC(NRDOMN) == '      ON')THEN
          IF(RDOM(K,I).GT.0.0)THEN
          ORGNRD(K,I)=RDOMN(K,I)/RDOM(K,I)
          ELSE
          ORGNRD(K,I)=ORGN(JW)
          ENDIF
        ELSE
          ORGNRD(K,I)=ORGN(JW)
        END IF
        IF(CAC(NLPOMN) == '      ON')THEN
          IF(LPOM(K,I).GT.0.0)THEN
          ORGNLP(K,I)=LPOMN(K,I)/LPOM(K,I)
          ELSE
          ORGNLP(K,I)=ORGN(JW)
          ENDIF
        ELSE
          ORGNLP(K,I)=ORGN(JW)
        END IF
        IF(CAC(NRPOMP) == '      ON')THEN
          IF(RPOM(K,I).GT.0.0)THEN
          ORGNRP(K,I)=RPOMN(K,I)/RPOM(K,I)
          ELSE
          ORGNRP(K,I)=ORGN(JW)
          ENDIF
        ELSE
          ORGNRP(K,I)=ORGN(JW)
        END IF

! Light Extinction Coefficient
      IF (.NOT. READ_EXTINCTION(JW)) THEN
      ALGEX = 0.0; SSEXT = 0.0; ZOOEXT = 0.0                                                     ! SW 11/8/07
        DO JA=1,NAL
          IF(ALG_CALC(JA))ALGEX = ALGEX+EXA(JA)*ALG(K,I,JA)
        END DO
        DO JS=1,NSS
          SSEXT = SSEXT+EXSS(JW)*SS(K,I,JS)
        END DO
        TOTMAC=0.0
        DO M=1,NMC
          IF(MACROPHYTE_CALC(JW,M))THEN
            JT=KTI(I)
            JE=KB(I)
            DO JJ=JT,JE
              TOTMAC = EXM(M)*MACRM(JJ,K,I,M)+TOTMAC
            END DO
          END IF
        END DO
        MACEXT=TOTMAC/(BH2(K,I)*DLX(I))

	    IF(ZOOPLANKTON_CALC)THEN
	        DO JZ = 1,NZP
	        ZOOEXT = ZOOEXT + ZOO(K,I,JZ)*EXZ(JZ)
	        END DO
	    ENDIF
        IF(KTICOL(I))THEN
          JT=KTI(I)
        ELSE
          JT=KTI(I)+1
        END IF
        JE=KB(I)
        DO JJ=JT,JE
          TOTMAC=0.0
          DO M=1,NMC
            IF(MACROPHYTE_CALC(JW,M))THEN
              TOTMAC = EXM(M)*MACRM(JJ,K,I,M)+TOTMAC
            END IF
          END DO
          IF(CW(JJ,I).GT.0.0)THEN
            MACEXT=TOTMAC/(CW(JJ,I)*DLX(I)*H2(K,I))
          ELSE
            MACEXT=0.0
          END IF
          GAMMAJ(JJ,K,I) = EXH2O(JW)+SSEXT+EXOM(JW)*(LPOM(K,I)+RPOM(K,I))+ALGEX+MACEXT+ZOOEXT
        END DO
      GAMMA(K,I) = EXH2O(JW)+SSEXT+EXOM(JW)*(LPOM(K,I)+RPOM(K,I))+ALGEX+MACEXT+ZOOEXT        ! SW 11/8/07 MOVED
      ELSE
      GAMMA(K,I) = EXH2O(JW)
      END IF

! Zooplankton Rates
   IF(ZOOPLANKTON_CALC)THEN
      DO JZ=1,NZP
        TGRAZE(K,I,JZ)=PREFP(JZ)*LPOM(K,I)
        DO JJZ = 1, NZP
          TGRAZE(K,I,JZ) = TGRAZE(K,I,JZ) + PREFZ(JJZ,JZ)*ZOO(K,I,JJZ)          !CB 5/17/2007
      END DO
        DO JA=1,NAL
          IF(ALG_CALC(JA))TGRAZE(K,I,JZ)=PREFA(JA,JZ)*ALG(K,I,JA)+TGRAZE(K,I,JZ)
        END DO
        ZMINFAC  = (1.0+SIGN(1.0,ZOO(K,I,JZ)-ZOOMIN(JZ)))*0.5
        ZRT(K,I,JZ) =  ZOORMR(K,I,JZ)*ZR(JZ)*ZMINFAC*DO3(K,I)
        IF (TGRAZE(K,I,JZ) <= 0.0 .OR. O2(K,I) < 2.0) THEN
          ZMU(K,I,JZ)       = 0.0
          AGZ(K,I,1:NAL,JZ) = 0.0
		  ZGZ(K,I,JZ,:) = 0.0
          IF (O2(K,I) < 2.0) ZMINFAC = 2*ZMINFAC
        ELSE
          ZMU(K,I,JZ) = MAX(ZOORM(K,I,JZ)*ZG(JZ)*(TGRAZE(K,I,JZ)-ZOOMIN(JZ))/(TGRAZE(K,I,JZ)+ZS2P(JZ)), 0.0)
          DO JA=1,NAL
          IF(ALG_CALC(JA))AGZ(K,I,JA,JZ) = ZMU(K,I,JZ)*ZOO(K,I,JZ)*(ALG(K,I,JA)*PREFA(JA,JZ)/TGRAZE(K,I,JZ))                      !  KV 5/26/2007
          END DO
          DO JJZ = 1,NZP ! OMNIVOROUS ZOOPLANKTON
          ZGZ(K,I,JJZ,JZ)  = ZMU(K,I,JZ)*ZOO(K,I,JZ)*(ZOO(K,I,JJZ)*PREFZ(JJZ,JZ)/TGRAZE(K,I,JZ))         !KV 5/26/2007
          END DO
        END IF
        ZMT(K,I,JZ) = MAX(1.0-ZOORMF(K,I,JZ),0.02)*ZM(JZ)*ZMINFAC
      END DO   ! ZOOP LOOP
   ENDIF

    END DO ! K LOOP
  END DO   ! I LOOP
!$OMP END PARALLEL DO

! Algal rates
   DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
      do i=iu,id
!**** Limiting factor
      LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)/ASAT(JA)
      LAM1  =  LIGHT
      LAM2  =  LIGHT
      DO K=KT,KB(I)

!****** Limiting factor
        LAM1           = LAM2
        LAM2           = LAM1*EXP(-GAMMA(K,I)*H2(K,I))
        FDPO4          = 1.0-FPSS(K,I)-FPFE(K,I)
        ALLIM(K,I,JA)  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA(K,I)*H2(K,I))
        IF (AHSP(JA)  /= 0.0) APLIM(K,I,JA) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+AHSP(JA)+NONZERO)
        IF (AHSN(JA)  /= 0.0) ANLIM(K,I,JA) = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+AHSN(JA)+NONZERO)
        IF (AHSSI(JA) /= 0.0) ASLIM(K,I,JA) =  DSI(K,I)/(DSI(K,I)+AHSSI(JA)+NONZERO)
        LIMIT          = MIN(APLIM(K,I,JA),ANLIM(K,I,JA),ASLIM(K,I,JA),ALLIM(K,I,JA))

!****** Algal rates
        AGR(K,I,JA) =  ATRM(K,I,JA)*AG(JA)*LIMIT
        ARR(K,I,JA) =  ATRM(K,I,JA)*AR(JA)*DO3(K,I)
        AMR(K,I,JA) = (ATRMR(K,I,JA)+1.0-ATRMF(K,I,JA))*AM(JA)
        AER(K,I,JA) =  MIN((1.0-ALLIM(K,I,JA))*AE(JA)*ATRM(K,I,JA),AGR(K,I,JA))
        IF (AS(JA) >= 0.0) THEN
          IF(K == KT)THEN
          ASR(K,I,JA) =  AS(JA)*(-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ELSE
          ASR(K,I,JA) =  AS(JA)*(ALG(K-1,I,JA)-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ENDIF
        ELSE
          IF(K == KB(I))THEN
            ASR(K,I,JA) = -AS(JA)*(-ALG(K,I,JA)  *BI(K,I)/BH2(K,I))                                           !SW 11/8/07
          ELSEIF(K == KT)THEN
            ASR(K,I,JA) = -AS(JA)* ALG(K+1,I,JA)*BI(K+1,I)*DLX(I)/VOL(K,I)                                   !SW 11/8/07
          ELSE
            ASR(K,I,JA) = -AS(JA)*(ALG(K+1,I,JA)*BI(K+1,I)/BH2(K,I)-ALG(K,I,JA)*BI(K,I)/BH2(K,I))             !SP 8/27/07
          END IF
        END IF
      end do
    end do
    ENDIF
  END DO    ! ALGAE LOOP

! Macrophyte Light/Nutrient Limitation and kinetic rates
  do m=1,nmc
  if(macrophyte_calc(jw,m))then
    DO I=IU,ID
      LTCOEFm = (1.0-BETA(jw))*SRON(jw)*SHADE(I)
      if(kticol(i))then
        jt=kti(i)
      else
        jt=kti(i)+1
      end if
      je=kb(i)
      do jj=jt,je
        lam1=ltcoefm
        lam2m(jj,kt)=lam1*exp(-gammaj(jj,kt,i)*h2(kt,i))
        lavg=(lam1-lam2m(jj,kt))/(GAMMAj(jj,kt,i)*H2(kt,i))
        mLLIM(jj,kt,I,m) = lavg/(lavg+msat(m))
        IF (mHSP(m)  /= 0.0.and.psed(m) < 1.0)then
          mPLIM(kt,I,m) =  FDPO4*PO4(kt,I)/(FDPO4*PO4(kt,I)+mHSP(m)+nonzero)
        else
          mPLIM(kt,I,m)=1.0
        end if
        IF (mHSN(m)  /= 0.0.and.nsed(m) < 1.0)then
          mNLIM(kt,I,m) = NH4(kt,I)/(NH4(kt,I)+mHSN(m)+nonzero)
        else
          mNLIM(kt,I,m)=1.0
        end if
        IF (mHSc(m) /= 0.0)then
          mcLIM(kt,i,m) = co2(kt,I)/(co2(kt,I)+mHSc(m)+NONZERO)
        end if
        LIMIT          = MIN(mPLIM(kt,I,m),mNLIM(kt,I,m),mcLIM(kt,I,m),mLLIM(jj,kt,I,m))

!************* sources/sinks

        mGR(jj,Kt,I,m) = macTRM(Kt,I,m)*mG(m)*LIMIT

      end do

      mRR(Kt,I,m) = macTRM(Kt,I,m)*mR(m)*DO3(Kt,I)
      mMR(Kt,I,m) = (macTRMR(Kt,I,m)+1.0-mAcTRMF(Kt,I,m))*mM(m)

      DO K=KT+1,KB(I)
        jt=k
        je=kb(i)
        do jj=jt,je
          lam1=lam2m(jj,k-1)
          lam2m(jj,k)=lam1*exp(-gammaj(jj,k,i)*h2(k,i))
          lavg=(lam1-lam2m(jj,k))/(GAMMAj(jj,k,i)*H2(k,i))
          mLLIM(jj,K,I,m) = lavg/(lavg+msat(m))
          IF (mHSP(m)  /= 0.0.and.psed(m) < 1.0)then
            mPLIM(K,I,m) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+mHSP(m)+nonzero)
          else
            mPLIM(K,I,m)=1.0
          end if
          IF (mHSN(m)  /= 0.0.and.nsed(m) < 1.0)then
            mNLIM(K,I,m) = NH4(K,I)/(NH4(K,I)+mHSN(m)+nonzero)
          else
             mNLIM(K,I,m)=1.0
          end if
          IF (mHSc(m) /= 0.0)then
            mcLIM(k,i,m) = co2(K,I)/(co2(K,I)+mHSc(m)+NONZERO)
          end if
          LIMIT          = MIN(mPLIM(K,I,m),mNLIM(K,I,m),mcLIM(K,I,m),mLLIM(jj,K,I,m))

!************* sources/sinks

          mGR(jj,K,I,m) = macTRM(K,I,m)*mG(m)*LIMIT

        end do

        mRR(K,I,m) = macTRM(K,I,m)*mR(m)*DO3(K,I)
        mMR(K,I,m) = (macTRMR(K,I,m)+1.0-mAcTRMF(K,I,m))*mM(m)
      end do
    END DO
    ENDIF
  END DO

RETURN

!***********************************************************************************************************************************
!**                                             G E N E R I C   C O N S T I T U E N T                                             **
!***********************************************************************************************************************************

ENTRY GENERIC_CONST (JG)
xx=0.0
DO I=IU,ID
      DO K=KT,KB(I)

         IF (CGS(JG) > 0.0) THEN
          IF(K == KT)THEN
          xx =  CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)    ! AS(JA)*(-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ELSE
          xx =  CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)     !AS(JA)*(ALG(K-1,I,JA)-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ENDIF
         ELSEif(cgs(jg)<0.0)then
          IF(K == KB(I))THEN
            xx = -CGS(JG)*(-CG(K,I,JG))*BI(K,I)/BH2(K,I)    !-AS(JA)*(-ALG(K,I,JA)  *BI(K,I)/BH2(K,I))                                           !SW 11/8/07
          ELSEIF(K == KT)THEN
            xx = -CGS(JG)*CG(K+1,I,JG)*BI(K+1,I)*DLX(I)/VOL(K,I)    !-AS(JA)* ALG(K+1,I,JA)*BI(K+1,I)*DLX(I)/VOL(K,I)                                   !SW 11/8/07
          ELSE
            xx = -CGS(JG)*(CG(K+1,I,JG)*BI(K+1,I)/BH2(K,I)-CG(K,I,JG)*BI(K,I)/BH2(K,I))    !-AS(JA)*(ALG(K+1,I,JA)*BI(K+1,I)/BH2(K,I)-ALG(K,I,JA)*BI(K,I)/BH2(K,I))             !SP 8/27/07
          END IF
         ENDIF

         IF (CGQ10(JG) /= 0.0) THEN
             CGSS(K,I,JG) = -CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K,I,JG)+xx            ! SW 4/5/09 CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)
         ELSE
             CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CG(K,I,JG)+xx                                                                ! SW 4/5/09 CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)              
         ENDIF
     END DO
   END DO
 RETURN  



!  IF (CGQ10(JG) /= 0.0) THEN
!    DO I=IU,ID
!      DO K=KT,KB(I)
!        CGSS(K,I,JG) = -CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K,I,JG)+xx            ! SW 4/5/09 CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)
!      END DO
!    END DO
!  ELSE
!    DO I=IU,ID
!      DO K=KT,KB(I)
!        CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CG(K,I,JG)+                                                                   ! SW 4/5/09 CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)
!      END DO
!    END DO
!  END IF
!RETURN

!***********************************************************************************************************************************
!**                                               S U S P E N D E D   S O L I D S                                                 **
!***********************************************************************************************************************************

ENTRY SUSPENDED_SOLIDS (J)
  DO I=IU,ID
    SSR = 0.0
    IF (SEDIMENT_RESUSPENSION(J)) THEN
      FETCH = FETCHD(I,JB)
      IF (COS(PHI(JW)-PHI0(I)) < 0.0) FETCH = FETCHU(I,JB)
      FETCH = MAX(FETCH,BI(KT,I),DLX(I))
      U2    = WIND(JW)*WSC(I)*WIND(JW)*WSC(I)+NONZERO
      COEF1 = 0.53  *(G*DEPTHB(KT,I)/U2)**0.75
      COEF2 = 0.0125*(G*FETCH/U2)**0.42
      COEF3 = 0.833* (G*DEPTHB(KT,I)/U2)**0.375
      COEF4 = 0.077* (G*FETCH/U2)**0.25
      HS    = 0.283 *U2/G*0.283*TANH(COEF1)*TANH(COEF2/TANH(COEF1))
      TS    = 2.0*PI*U2/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))
      L0    = G*TS*TS/(2.0*PI)
    END IF
    SSSS(KT,I,J) = -SSS(J)*SS(KT,I,J)*BI(KT,I)/BH2(KT,I)+SSR
   ! DO K=KT-1,KB(I)-1                                             
    DO K=KT,KB(I)-1                                             ! JP 2/3/12
      IF (SEDIMENT_RESUSPENSION(J)) THEN
        L1 = L0
        L  = L0*TANH(2.0*PI*DEPTHB(K,I)/L1)
        DO WHILE (ABS(L-L1) > 0.001)
          L1 = L
          L  = L0*TANH(2.0*PI*DEPTHB(K,I)/L1)
        END DO
        COEF = MIN(710.0,2.0*PI*DEPTHB(K,I)/L)
        UORB = PI*HS/TS*100.0/SINH(COEF)
        TAU  = 0.003*UORB*UORB
        IF (TAU-TAUCR(J) > 0.0) EPSILON = MAX(0.0,0.008/49.0*(TAU-TAUCR(J))**3*10000.0/DLT)
		if(k == kb(i))then   ! SW 4/18/07
		SSR = EPSILON*DLX(I)*BI(K,I)/VOL(K,I)
		else
        SSR = EPSILON*DLX(I)*(BI(K,I)-BI(K+1,I))/VOL(K,I)
		endif
      END IF
      SSSS(K,I,J) = SSS(J)*(SS(K-1,I,J)-SS(K,I,J))*BI(K,I)/BH2(K,I)+SSR
    END DO
    IF (SEDIMENT_RESUSPENSION(J)) SSR = EPSILON*DLX(I)*BI(KB(I),I)/VOL(KB(I),I)
    SSSS(KB(I),I,J) = SSS(J)*(SS(KB(I)-1,I,J)-SS(KB(I),I,J))/H(KB(I),JW)+SSR
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      P H O S P H O R U S                                                      **
!***********************************************************************************************************************************

ENTRY PHOSPHORUS
  PO4AR(:,IU:ID) = 0.0; PO4AG(:,IU:ID) = 0.0; PO4ER(:,IU:ID) = 0.0; PO4EG(:,IU:ID) = 0.0; PO4BOD(:,IU:ID) = 0.0
  po4mr(:,iu:id) = 0.0; po4mg(:,iu:id) = 0.0; po4zr(:,iu:id)=0.0   !v3.5

  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        IF(BOD_CALC(JCB))PO4BOD(K,I) = PO4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODP(JCB)
      END DO
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        PO4AG(K,I) = PO4AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*AP(JA)
        PO4AR(K,I) = PO4AR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*AP(JA)
        ENDIF
      END DO
      DO JE=1,NEP
      IF (EPIPHYTON_CALC(JW,JE))then
        PO4EG(K,I) = PO4EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*EP(JE)
        PO4ER(K,I) = PO4ER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*EP(JE)
      endif
      END DO
      PO4EP(K,I)  = PO4ER(K,I)-PO4EG(K,I)
      PO4AP(K,I)  = PO4AR(K,I)-PO4AG(K,I)
! v3.5 start
      PO4POM(K,I) = ORGPLP(k,i)*LPOMD(K,I)+orgprp(k,i)*RPOMD(K,I)
      PO4DOM(K,I) = ORGPLD(k,i)*LDOMD(K,I)+orgprd(k,i)*RDOMD(K,I)
! v3.5 end
      PO4OM(K,I)  = PO4POM(K,I)+PO4DOM(K,I)
! v3.5 start!
      PO4SD(K,I)  = SEDDp(K,I)
! v3.5 end
      PO4SR(K,I)  = PO4R(JW)*SODD(K,I)*DO2(K,I)
      PO4NS(K,I)  = SSSI(K,I)*PO4(K-1,I)-SSSO(K,I)*PO4(K,I)
! v3.5 start
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            po4mg(k,i)= po4mg(k,i)+mgr(jj,k,i,m)*macrm(jj,k,i,m)*mp(m)*(1.0-psed(m))
            po4mr(k,i)= po4mr(k,i)+mrr(k,i,m)*macrm(jj,k,i,m)*mp(m)
          end do
        end if
      end do
      po4mr(k,i)=po4mr(k,i)/(dlx(i)*bh(k,i))
      po4mg(k,i)=po4mg(k,i)/(dlx(i)*bh(k,i))
      IF(ZOOPLANKTON_CALC)THEN
      do jz = 1,nzp
        PO4zR(K,I) = po4zr(k,i) + zrt(K,I,jz)*zoo(K,I,jz)*zp(jz)
	  end do
	  ENDIF
! v3.5 end

      PO4SS(K,I)  = PO4AP(K,I)+PO4EP(K,I)+PO4OM(K,I)+PO4SD(K,I)+PO4SR(K,I)+PO4NS(K,I)+PO4BOD(K,I)  &
                    +po4mr(k,i)-po4mg(k,i) +po4zr(k,i)     !v3.5

    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                        A M M O N I U M                                                        **
!***********************************************************************************************************************************

ENTRY AMMONIUM
  NH4AG(:,IU:ID) = 0.0; NH4AR(:,IU:ID) = 0.0; NH4ER(:,IU:ID) = 0.0; NH4EG(:,IU:ID) = 0.0; NH4BOD(:,IU:ID) = 0.0
  nh4mg(:,iu:id) = 0.0; nh4mr(:,iu:id) = 0.0; nh4zr(:,iu:id)=0.0   
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        IF(BOD_CALC(JCB))NH4BOD(K,I) =  NH4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODN(JCB)
      END DO
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        IF (ANEQN(JA).EQ.2) THEN
        NH4PR      = NH4(K,I)*NO3(K,I)/((ANPR(JA)+NH4(K,I))*(ANPR(JA)+NO3(K,I)))+NH4(K,I)*ANPR(JA)/((NO3(K,I)  &
                                        +NH4(K,I)+NONZERO)*(ANPR(JA)+NO3(K,I)))
        ELSE
        NH4PR = NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        ENDIF
        IF (AHSN(JA) > 0.0) NH4AG(K,I) = NH4AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*AN(JA)*NH4PR
        NH4AR(K,I) = NH4AR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*AN(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        IF (ENEQN(JE) == 2)THEN
        NH4PR = NH4(K,I)*NO3(K,I)/((ENPR(JE)+NH4(K,I))*(ENPR(JE)+NO3(K,I)))+NH4(K,I)*ENPR(JE)/((NO3(K,I)  &
                                        +NH4(K,I)+NONZERO)*(ENPR(JE)+NO3(K,I)))
        ELSE
        NH4PR = NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        ENDIF
        IF (EHSN(JE) > 0.0) NH4EG(K,I) = NH4EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*EN(JE)*NH4PR
        NH4ER(K,I) = NH4ER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*EN(JE)
        endif
      END DO
      NH4EP(K,I)  =  NH4ER(K,I) -NH4EG(K,I)
      NH4AP(K,I)  =  NH4AR(K,I) -NH4AG(K,I)

      NH4DOM(K,I) = LDOMD(K,I)*orgnld(k,i) +RDOMD(K,I)*ORGNrd(k,i)
      NH4POM(K,I) = LPOMD(K,I)*orgnlp(k,i) +RPOMD(K,I)*ORGNrp(k,i)

      NH4OM(K,I)  =  NH4DOM(K,I)+NH4POM(K,I)

      NH4SD(K,I)  =  SEDDn(K,I)
      NH4SR(K,I)  =  NH4R(JW) *SODD(K,I)*DO2(K,I)

      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            nh4mr(k,i)= nh4mr(k,i)+mrr(k,i,m)*macrm(jj,k,i,m)*mn(m)
            nh4mg(k,i)= nh4mg(k,i)+mgr(jj,k,i,m)*macrm(jj,k,i,m)*mn(m)*(1.0-nsed(m))
          end do
        end if
      end do
      nh4mr(k,i)=nh4mr(k,i)/(dlx(i)*bh(k,i))
      nh4mg(k,i)=nh4mg(k,i)/(dlx(i)*bh(k,i))
	  IF(ZOOPLANKTON_CALC)THEN
	  do jz = 1,nzp
	    nh4zr(k,i) = nh4zr(k,i) + zrt(k,i,jz)*zoo(k,i,jz)*zn(jz) ! Baulk
	  end do
	  ENDIF
      NH4SS(K,I)  =  NH4AP(K,I)+NH4EP(K,I)+NH4OM(K,I)+NH4SD(K,I)+NH4SR(K,I)+NH4BOD(K,I)-NH4D(K,I)  &
         +nh4mr(k,i)-nh4mg(k,i) +nh4zr(k,i)	     !v3.5
! v3.5 end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                          N I T R A T E                                                        **
!***********************************************************************************************************************************

ENTRY NITRATE
  NO3AG(:,IU:ID) = 0.0; NO3EG(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        NO3PR = 1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ANEQN(JA).EQ.2)  NO3PR      = 1.0-(NH4(K,I)*NO3(K,I)/((ANPR(JA)+NH4(K,I))*(ANPR(JA)+NO3(K,I)))+NH4(K,I)*ANPR(JA)       &
                                          /((NO3(K,I)+NH4(K,I)+NONZERO)*(ANPR(JA)+NO3(K,I))))
        IF (AHSN(JA).GT.0.0) NO3AG(K,I) = NO3AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*NO3PR*AN(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        NO3PR = 1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ENEQN(JE).EQ.2)  NO3PR      = 1.0-(NH4(K,I)*NO3(K,I)/((ENPR(JE)+NH4(K,I))*(ENPR(JE)+NO3(K,I)))+NH4(K,I)*ENPR(JE)       &
                                          /((NO3(K,I)+NH4(K,I)+NONZERO)*(ENPR(JE)+NO3(K,I))))
        IF (EHSN(JE).GT.0.0) NO3EG(K,I) = NO3EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*NO3PR*EN(JE)
        endif
      END DO
      if(k == kb(i)) then      ! SW 4/18/07
      NO3SED(K,I) = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I))/BH2(K,I)
	  else
      NO3SED(K,I) = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I)-BI(K+1,I))/BH2(K,I)
	  endif
      NO3SS(K,I)  = NH4D(K,I)-NO3D(K,I)-NO3AG(K,I)-NO3EG(K,I)-NO3SED(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  D I S S O L V E D   S I L I C A                                              **
!***********************************************************************************************************************************

ENTRY DISSOLVED_SILICA
  DSIAG(:,IU:ID) = 0.0; DSIEG(:,IU:ID) = 0.0                          !; DSIBOD = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        DSIAG(K,I) = DSIAG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*ASI(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))DSIEG(K,I) = DSIEG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*ESI(JE)
      END DO
      DSID(K,I)  =  PSIDK(JW)*PSI(K,I)
      DSISD(K,I) =  SEDD(K,I)*ORGSI(JW)
      DSISR(K,I) =  DSIR(JW)*SODD(K,I)*DO2(K,I)
      DSIS(K,I)  = (SSSI(K,I)*DSI(K-1,I)-SSSO(K,I)*DSI(K,I))*PARTSI(JW)
      DSISS(K,I) =  DSID(K,I)+DSISD(K,I)+DSISR(K,I)+DSIS(K,I)-DSIAG(K,I)-DSIEG(K,I)    !+DSIBOD
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                P A R T I C U L A T E   S I L I C A                                            **
!***********************************************************************************************************************************

ENTRY PARTICULATE_SILICA
  PSIAM(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        PSIAM(K,I) = PSIAM(K,I)+AMR(K,I,JA)*PSI(K,I)*ASI(JA)
      ENDIF
      END DO
      PSID(K,I)  = PSIDK(JW)*PSI(K,I)
      PSINS(K,I) = PSIS(JW)*(PSI(K-1,I)*DO1(K-1,I)-PSI(K,I)*DO1(K,I))*BI(K,I)/BH2(K,I)
      PSISS(K,I) = PSIAM(K,I)-PSID(K,I)+PSINS(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                            I R O N                                                            **
!***********************************************************************************************************************************

ENTRY IRON
  DO I=IU,ID
    DO K=KT,KB(I)
      FENS(K,I) = FES(JW)*(FE(K-1,I)*DO1(K-1,I)-FE(K,I)*DO1(K,I))*BI(K,I)/BH2(K,I)
      FESR(K,I) = FER(JW)*SODD(K,I)*DO2(K,I)
      FESS(K,I) = FESR(K,I)+FENS(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M                                                     **
!***********************************************************************************************************************************

ENTRY LABILE_DOM
  LDOMAP(:,IU:ID) = 0.0; LDOMEP(:,IU:ID) = 0.0; ldommac(:,iu:id)= 0.0  !v3.5
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LDOMAP(K,I) = LDOMAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LDOMEP(K,I) = LDOMEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)
      END DO
! v3.5 start
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            ldommac(k,i)=ldommac(k,i)+(1.0-mpom(m))*mmr(k,i,m)*macrm(jj,k,i,m)
          end do
        end if
      end do
      ldommac(k,i)=ldommac(k,i)/(dlx(i)*bh(k,i))
      LDOMSS(K,I) = LDOMAP(K,I)+LDOMEP(K,I)-LDOMD(K,I)-LRDOMD(K,I)+ldommac(k,i)
! v3.5 end

    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M                                                  **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMSS(K,I) = LRDOMD(K,I)-RDOMD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M                                                      **
!***********************************************************************************************************************************

ENTRY LABILE_POM
  LPOMAP(:,IU:ID) = 0.0; lpomep(:,iu:id) = 0.0;   lpommac(:,iu:id) = 0.0; lpzooin(:,iu:id)=0.0;lpzooout(:,iu:id)=0.0  !v3.5  ! cb 5/19/06
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LPOMAP(K,I) = LPOMAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))
      END DO
      DO JE=1,NEP                                                          ! cb 5/19/06
        IF (EPIPHYTON_CALC(JW,JE))LPOMEP(K,I) = LPOMEP(K,I)+EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))       ! cb 5/19/06
      END DO                                                               ! cb 5/19/06
      LPOMNS(K,I) = POMS(JW)*(LPOM(K-1,I)-LPOM(K,I))*BI(K,I)/BH2(K,I)
! v3.5 start
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          do jj=jt,je
            lpommac(k,i)=lpommac(k,i)+mpom(m)*lrpmac(m)*mmr(k,i,m)*macrm(jj,k,i,m)
          end do
        end if
      end do
      lpommac(k,i)=lpommac(k,i)/(dlx(i)*bh(k,i))
      IF(ZOOPLANKTON_CALC)THEN
      do jz = 1,nzp
        if(tgraze(k,i,jz) > 0.0)then
          lpzooout(k,i)=lpzooout(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))
          lpzooin(k,i)=lpzooin(k,i)+zoo(k,i,jz)*prefp(jz)*zmu(k,i,jz)*LPOM(k,i)/tgraze(k,i,jz)
        else
          lpzooout(k,i)=lpzooout(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))
          lpzooin(k,i)=0.0
        end if
      end do
      ENDIF
      LPOMSS(K,I) = LPOMAP(K,I)+lpomep(k,i)-LPOMD(K,I)+LPOMNS(K,I)-LRPOMD(K,I)+lpommac(k,i)+lpzooout(k,i)-lpzooin(k,i)       ! cb 5/19/06
! v3.5 end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M                                                  **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM
  rpommac(:,iu:id) = 0.0 !v3.5
  DO I=IU,ID
    DO K=KT,KB(I)
      RPOMNS(K,I) = POMS(JW)*(RPOM(K-1,I)-RPOM(K,I))*BI(K,I)/BH2(K,I)
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          do jj=jt,je
            rpommac(k,i)=rpommac(k,i)+mpom(m)*(1.0-lrpmac(m))*mmr(k,i,m)*macrm(jj,k,i,m)
          end do
        end if
      end do
      rpommac(k,i)=rpommac(k,i)/(dlx(i)*bh(k,i))
      RPOMSS(K,I) = LRPOMD(K,I)+RPOMNS(K,I)-RPOMD(K,I)+rpommac(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                         A L G A E                                                             **
!***********************************************************************************************************************************

ENTRY ALGAE (J)
  agzt(:,iu:id,j) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      IF(ZOOPLANKTON_CALC)THEN
      do jz = 1,nzp
	  agzt(k,i,j) = agzt(k,i,j) + agz(k,i,j,jz)                       ! cb 5/26/07
	  end do
	  ENDIF
      ASS(K,I,J) = ASR(K,I,J)+(AGR(K,I,J)-AER(K,I,J)-AMR(K,I,J)-ARR(K,I,J))*ALG(K,I,J)-agzt(k,i,j)	
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                            B I O C H E M I C A L   O 2   D E M A N D                                          **
!***********************************************************************************************************************************

ENTRY BIOCHEMICAL_O2_DEMAND(JBOD)
  if(jbod == 1)cbodns(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      cbodset = cbodS(jbod)*(cbod(K-1,I,jbod)-cbod(K,I,jbod))*BI(K,I)/BH2(K,I)
      cbodns(k,i)=cbodns(k,i)+cbodset
      CBODSS(K,I,JBOD) = -CBODD(K,I,JBOD)*CBOD(K,I,JBOD)+cbodset
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                D I S S O L V E D   O X Y G E N                                                **
!***********************************************************************************************************************************

ENTRY DISSOLVED_OXYGEN
  DOAP(:,IU:ID) = 0.0; DOAR(:,IU:ID) = 0.0; DOEP(:,IU:ID) = 0.0; DOER(:,IU:ID) = 0.0; DOBOD(:,IU:ID) = 0.0
  domp(:,iu:id) = 0.0; domr(:,iu:id) = 0.0; dozr(:,iu:id)=0.0   !v3.5

  DO I=IU,ID
    DOSS(KT,I) = 0.0
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        IF(BOD_CALC(JCB))DOBOD(K,I) = DOBOD(K,I)+RBOD(JCB)*CBODD(K,I,JCB)*CBOD(K,I,JCB)
      END DO
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        DOAP(K,I) = DOAP(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*O2AG(JA)
        DOAR(K,I) = DOAR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*O2AR(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        DOEP(K,I) = DOEP(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*O2EG(JE)
        DOER(K,I) = DOER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*O2ER(JE)
        endif
      END DO
! v3.5 start
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            domp(k,i)=domp(k,i)+mgr(jj,k,i,m)*macrm(jj,k,i,m)*o2mg(m)
            domr(k,i)=domr(k,i)+mrr(k,i,m)*macrm(jj,k,i,m)*o2mr(m)
          end do
        end if
      end do
      domp(k,i)=domp(k,i)/(dlx(i)*bh(k,i))
      domr(k,i)=domr(k,i)/(dlx(i)*bh(k,i))
      DOPOM(K,I) = (LPOMD(K,I)+RPOMD(K,I))*O2OM(JW)
      DODOM(K,I) = (LDOMD(K,I)+RDOMD(K,I))*O2OM(JW)
      DOOM(K,I)  =  DOPOM(K,I)+DODOM(K,I)+DOBOD(K,I)      
      DONIT(K,I) =  NH4D(K,I)*O2NH4(JW)
      DOSED(K,I) =  SEDD(K,I)*O2OM(JW)
      DOSOD(K,I) =  SODD(K,I)*DO3(K,I)
     IF(ZOOPLANKTON_CALC)THEN
     do jz = 1, nzp
      dozr(k,i)  = dozr(k,i)+zrt(k,i,jz)*zoo(k,i,jz)*o2zr(jz)
	 end do
	 ENDIF
    DOSS(K,I)  =  DOAP(K,I)+DOEP(K,I)-DOAR(K,I)-DOER(K,I)-DOOM(K,I)-DONIT(K,I)-DOSOD(K,I)-DOSED(K,I)  &
                    +domp(k,i)-domr(k,i)-dozr(k,i)
    END DO
    DOSAT = SATO(T1(KT,I),TDS(KT,I),PALT(I),SALT_WATER(JW))
    IF (.NOT. ICE(I)) THEN
      CALL GAS_TRANSFER
      O2EX       =  REAER(I)
      DOAE(KT,I) = (DOSAT-O2(KT,I))*O2EX*BI(KT,I)/BH2(KT,I)
      DOSS(KT,I) =  DOSS(KT,I)+DOAE(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                              I N O R G A N I C   C A R B O N                                                  **
!***********************************************************************************************************************************

ENTRY INORGANIC_CARBON
  TICAP(:,IU:ID) = 0.0; TICEP(:,IU:ID) = 0.0; TICBOD(:,IU:ID) = 0.0
  ticmc(:,iu:id) = 0.0; ticzr(:,iu:id)=0.0  !v3.5
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        IF(BOD_CALC(JCB))TICBOD(K,I) = TICBOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODC(JCB)
      END DO
      DO JA=1,NAL
        IF(ALG_CALC(JA))TICAP(K,I) = TICAP(K,I)+AC(JA)*(ARR(K,I,JA)-AGR(K,I,JA))*ALG(K,I,JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))TICEP(K,I) = TICEP(K,I)+EC(JE)*(ERR(K,I,JE)-EGR(K,I,JE))*EPC(K,I,JE)
      END DO
! v3.5 start
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            ticmc(k,i)=ticmc(k,i)+(mrr(k,i,m)-mgr(jj,k,i,m))*macrm(jj,k,i,m)*mc(m)
          end do
        end if
      end do
      ticmc(k,i)=ticmc(k,i)/(dlx(i)*bh(k,i))
      IF(ZOOPLANKTON_CALC)THEN
      do jz = 1,nzp
        ticzr(k,i)=ticzr(k,i)+zrt(k,i,jz)*zoo(k,i,jz)*zc(jz) !mlm
	  end do
	  ENDIF
      TICSS(K,I) = TICAP(K,I)+TICEP(K,I)+seddc(k,i)+ORGC(JW)*(LPOMD(K,I)+RPOMD(K,I)+LDOMD(K,I)+RDOMD(K,I))                          &
                   +CO2R(JW)*SODD(K,I)*DO3(K,I)+TICBOD(K,I)+ticmc(k,i)+ticzr(k,i)      
! v3.5 end
    END DO
    IF (.NOT. ICE(I)) THEN
      IF (REAER(I) == 0.0) CALL GAS_TRANSFER
      CO2EX       = REAER(I)*0.923
      TICSS(KT,I) = TICSS(KT,I)+CO2EX*(0.286*EXP(-0.0314*(T2(KT,I))*PALT(I))-CO2(KT,I))*BI(KT,I)/BH2(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T                                                          **
!***********************************************************************************************************************************

ENTRY SEDIMENT
  SEDAS(:,IU:ID) = 0.0; LPOMEP(:,IU:ID) = 0.0; SEDCB(:,IU:ID) = 0.0
  DO I=IU,ID
    sedsi=0.0
    DO K=KT,KB(I)
    IF(K == KB(I))THEN
    BIBH2(K,I)=BI(K,I)/BH2(K,I)
    ELSE
    BIBH2(K,I)=BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
    ENDIF
      DO JA=1,NAL
        IF(ALG_CALC(JA))SEDAS(K,I) = SEDAS(K,I)+MAX(AS(JA),0.0)*ALG(K,I,JA)*BIBH2(K,I)                !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      sedem = 0.0   ! cb 5/19/06
      DO JE=1,NEP
!        LPOMEP(K,I) = LPOMEP(K,I)+EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
        IF (EPIPHYTON_CALC(JW,JE))sedem = sedem+ebr(k,i,je)/h1(k,i)*EPC(K,I,JE)    ! cb 5/19/06
      END DO
      do jd=1,nbod
        IF(BOD_CALC(JD))SEDcb(K,I) = SEDcb(K,I)+MAX(cbods(jd),0.0)*cbod(K,I,Jd)*BIBH2(K,I)/O2OM(JW)           !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      end do
      SEDOMS(K,I) = pomS(JW)*(LPOM(K,I)+RPOM(K,I))*BIBH2(K,I)                        !cb 10/22/06
      IF(K==KB(I))THEN
      SEDSO       = 0.0
      ELSE
      SEDSO       = sedS(JW)*SED(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      ENDIF
      SEDNS(K,I)  = SEDSI-SEDSO
      SEDSI       = SEDSO
      SED(K,I)    = MAX(SED(K,I)+(sedem+SEDAS(K,I)+sedcb(k,i)+SEDOMS(K,I)+SEDNS(K,I)-SEDD(K,I)-sedbr(k,i))*DLT,0.0)   ! cb 11/30/06
    END DO
  END DO
RETURN


!***********************************************************************************************************************************
!**                                                      S E D I M E N T   P H O S P H O R U S                                    **
!***********************************************************************************************************************************

ENTRY SEDIMENTP
  SEDASp(:,IU:ID) = 0.0; LPOMEPp(:,IU:ID) = 0.0; sedcbp(:,IU:ID) = 0.0
  DO I=IU,ID
    sedsip=0.0
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))SEDASp(K,I) = SEDASp(K,I)+MAX(AS(JA),0.0)*ap(ja)*ALG(K,I,JA)*BIBH2(K,I)          !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LPOMEPp(K,I) = LPOMEPp(K,I)+EPOM(JE)*ep(je)*(EMR(K,I,JE)*EPC(K,I,JE))
      END DO
      do jd=1,nbod
        IF(BOD_CALC(JD))sedcbp(k,i)=sedcbp(k,i)+MAX(cbods(jd),0.0)*bodp(jd)*cbod(K,I,jd)*BIBH2(K,I)      !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      end do
      SEDOMSp(K,I) = poms(JW)*(LPOMp(K,I)+RPOMp(K,I))*BIBH2(K,I)                         !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !cb 10/22/06
      IF(K == KB(I))THEN
      SEDSOp       = 0.0
      ELSE
      SEDSOp       = seds(JW)*SEDp(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      ENDIF
      SEDNSp(K,I)  = SEDSIp-SEDSOp
      SEDSIp       = SEDSOp
      SEDp(K,I)    = MAX(SEDp(K,I)+(LPOMEPp(K,I)+SEDASp(K,I)+SEDOMSp(K,I)+sedcbp(k,i)+SEDNSp(K,I)-SEDDp(K,I)   &
                     -sedbrp(k,i))*DLT,0.0)                                                                 !cb 11/30/06
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   N I T R O G E N                                        **
!***********************************************************************************************************************************

ENTRY SEDIMENTn
  SEDASn(:,IU:ID) = 0.0; LPOMEPn(:,IU:ID) = 0.0; sedcbn(:,IU:ID) = 0.0
  DO I=IU,ID
    sedsin=0.0
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))SEDASn(K,I) = SEDASn(K,I)+MAX(AS(JA),0.0)*an(ja)*ALG(K,I,JA)*BIBH2(K,I)            !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LPOMEPn(K,I) = LPOMEPn(K,I)+EPOM(JE)*en(je)*(EMR(K,I,JE)*EPC(K,I,JE))
      END DO
      do jd=1,nbod
        IF(BOD_CALC(JD))sedcbn(k,i)=sedcbn(k,i)+MAX(cbods(jd),0.0)*bodn(jd)*cbod(K,I,jd)*BIBH2(K,I)        !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      end do
      SEDOMSn(K,I) = poms(JW)*(LPOMn(K,I)+RPOMn(K,I))*BIBH2(K,I)                           !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !cb 10/22/06
      
      if(k == kb(i)) then      ! SW 12/16/07
      sedNO3(K,I)  = FNO3SED(JW)*NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I))/BH2(K,I)
      SEDSOn       = 0.0
	  else
      sedNO3(K,I)  = FNO3SED(JW)*NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I)-BI(K+1,I))/BH2(K,I)
      SEDSOn       = seds(JW)*SEDn(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
	  endif
      SEDNSn(K,I)  = SEDSIn-SEDSOn
      SEDSIn       = SEDSOn
      SEDn(K,I)    = MAX(SEDn(K,I)+(LPOMEPn(K,I)+SEDASn(K,I)+SEDOMSn(K,I)+sedcbn(k,i)+SEDNSn(K,I)+sedno3(k,i)   &
                     -SEDDn(K,I)-sedbrn(k,i))*DLT,0.0)  !cb 11/30/06                    
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   C A R B O N                                            **
!***********************************************************************************************************************************

ENTRY SEDIMENTC
  SEDASc(:,IU:ID) = 0.0; LPOMEPc(:,IU:ID) = 0.0; sedcbc(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      sedsip=0.0
      DO JA=1,NAL
        IF(ALG_CALC(JA))SEDASc(K,I) = SEDASc(K,I)+MAX(AS(JA),0.0)*ac(ja)*ALG(K,I,JA)*BIBH2(K,I)             !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LPOMEPc(K,I) = LPOMEPc(K,I)+EPOM(JE)*ec(je)*(EMR(K,I,JE)*EPC(K,I,JE))
      END DO
      do jd=1,nbod
        IF(BOD_CALC(JD))sedcbc(k,i)=sedcbc(k,i)+MAX(cbods(jd),0.0)*bodc(jd)*cbod(K,I,jd)*BIBH2(K,I)         !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      end do
      SEDOMSc(K,I) = poms(JW)*orgc(jw)*(LPOM(K,I)+RPOM(K,I))*BIBH2(K,I)                     !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))   !cb 10/22/06
      IF(K == KB(I))THEN
      SEDSOc       = 0.0
      ELSE
      SEDSOc       = seds(JW)*SEDc(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      ENDIF
      SEDNSc(K,I)  = SEDSIc-SEDSOc
      SEDSIc       = SEDSOc
      SEDc(K,I)    = MAX(SEDc(K,I)+(LPOMEPc(K,I)+SEDASc(K,I)+SEDOMSc(K,I)+sedcbc(k,i)+SEDNSc(K,I)-SEDDc(K,I)    &
                     -sedbrc(k,i))*DLT,0.0)                                                                   !cb 11/30/06
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   D E C A Y    R A T E                                   **
!***********************************************************************************************************************************

ENTRY SEDIMENT_decay_rate
  DO I=IU,ID
    sedsidk=0.0
    DO K=KT,KB(I)
      sedsum=0.0
      sedsumk=0.0
      
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        xdum=MAX(AS(JA),0.0)*ALG(K,I,JA)*BIBH2(K,I)
        SEDsumk = sedsumk + xdum * lpomdk(jw)    
        sedsum  = sedsum  + xdum
        ENDIF
      END DO
      
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        xdum=EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
        SEDsumk = sedsumk + xdum * lpomdk(jw)
        sedsum  = sedsum  + xdum
        endif
      END DO
      
      do jd=1,nbod
        IF(BOD_CALC(JD))THEN
        xdum=MAX(cbods(jd),0.0)*cbod(K,I,Jd)*BIBH2(K,I)*RBOD(JD)/o2om(jw)
        sedsumk = sedsumk+xdum*CBODD(K,I,JD)               
        sedsum  = sedsum + xdum
        ENDIF
      end do
      
      sedsumk = sedsumk + poms(JW)*(LPOM(K,I)*lpomdk(jw)+RPOM(K,I)*rpomdk(jw))*BIBH2(K,I)        !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  ! cb 10/22/06
      sedsum  = sedsum  + poms(JW)*(LPOM(K,I)+RPOM(K,I))*BIBH2(K,I)
      
      sedsumk = sedsumk*DLT
      sedsum  = sedsum*DLT  
    
      if((sedsum+sed(k,i)) > 0.0)then
      sdkv(k,i)    = (sedsumk+sed(k,i) * sdkv(k,i))/(sedsum+ sed(k,i))
      else
      sdkv(k,i)=0.0
      endif
            
    END DO
  END DO
RETURN

! v3.5 end

!***********************************************************************************************************************************
!*                                                         E P I P H Y T O N                                                      **
!***********************************************************************************************************************************

ENTRY EPIPHYTON (J)
  DO I=IU,ID

!** Limiting factor

    LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)/ESAT(J)
    LAM2  =  LIGHT
    LAM1  =  LIGHT
    DO K=KT,KB(I)

!**** Limiting factor

      LAM1          = LAM2
      LAM2          = LAM1*EXP(-GAMMA(K,I)*H1(K,I))
      FDPO4         = 1.0-FPSS(K,I)-FPFE(K,I)
      ELLIM(K,I,J)  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA(K,I)*H1(K,I))
      IF (EHSP(J)  /= 0.0) EPLIM(K,I,J) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+EHSP(J)+NONZERO)
      IF (EHSN(J)  /= 0.0) ENLIM(K,I,J) = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+EHSN(J)+NONZERO)
      IF (EHSSI(J) /= 0.0) ESLIM(K,I,J) =  DSI(K,I)/(DSI(K,I)+EHSSI(J)+NONZERO)
      LIMIT         =  MIN(EPLIM(K,I,J),ENLIM(K,I,J),ESLIM(K,I,J),ELLIM(K,I,J))
      BLIM          =  1.0-EPD(K,I,J)/(EPD(K,I,J)+EHS(J))

!**** Sources/sinks

      EGR(K,I,J) =  MIN(ETRM(K,I,J)*EG(J)*LIMIT*BLIM,PO4(K,I)/(EP(J)*DLT*EPD(K,I,J)/H1(KT,I)+NONZERO),(NH4(K,I)+NO3(K,I))/(EN(J)   &
                    *DLT*EPD(K,I,J)/H1(K,I)+NONZERO))
      ERR(K,I,J) =  ETRM(K,I,J)*ER(J)*DO3(K,I)
      EMR(K,I,J) = (ETRMR(K,I,J)+1.0-ETRMF(K,I,J))*EM(J)
      EER(K,I,J) =  MIN((1.0-ELLIM(K,I,J))*EE(J)*ETRM(K,I,J),EGR(K,I,J))
!      EPD(K,I,J) =  MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,I,J)/(H1(K,I)*0.0025))*DLT,0.0)
      EPD(K,I,J) =  MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,I,J)/H1(K,I))*DLT,0.00)   ! cb 5/18/06
      if(k == kb(i)) then      ! SW 12/16/07
      EPM(K,I,J) =  EPD(K,I,J)*(BI(K,I)+2.0*H1(K,I))*DLX(I)
	  else
      EPM(K,I,J) =  EPD(K,I,J)*(BI(K,I)-BI(K+1,I)+2.0*H1(K,I))*DLX(I)
	  endif
      EPC(K,I,J) =  EPM(K,I,J)/VOL(K,I)
    END DO
  END DO
RETURN

! v3.5 start
!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M   P H O S P H O R U S                               **
!***********************************************************************************************************************************

ENTRY LABILE_DOM_P
  LDOMpAP(:,IU:ID) = 0.0; LDOMpEP(:,IU:ID) = 0.0; ldompmp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LDOMpAP(K,I) = LDOMpAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)*ap(ja)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LDOMpEP(K,I) = LDOMpEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)*ep(je)
      END DO
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            ldompmp(k,i)=ldompmp(k,i)+(1.0-mpom(m))*mmr(k,i,m)*macrm(jj,k,i,m)*mp(m)
          end do
        end if
      end do
      ldompmp(k,i)=ldompmp(k,i)/(dlx(i)*bh(k,i))
      LDOMpSS(K,I) = LDOMpAP(K,I)+LDOMpEP(K,I)+ldompmp(k,i)-(LDOMD(K,I)+LRDOMD(K,I))*orgpld(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M   P H O S P H O R U S                            **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM_P
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMpSS(K,I) = LRDOMD(K,I)*orgpld(k,i)-RDOMD(K,I)*orgprd(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M   P H O S P H O R U S                                **
!***********************************************************************************************************************************

ENTRY LABILE_POM_p
  LPOMpAP(:,IU:ID) = 0.0;lpompmp(:,iu:id)=0.0;lpzooinp(:,iu:id)=0.0; lpzoooutp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LPOMpAP(K,I) = LPOMpAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))*ap(ja)
      END DO
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          do jj=jt,je
            lpompmp(k,i)=lpompmp(k,i)+mpom(m)*lrpmac(m)*mmr(k,i,m)*macrm(jj,k,i,m)*mp(m)
          end do
        end if
      end do
      lpompmp(k,i)=lpompmp(k,i)/(dlx(i)*bh(k,i))
	IF(ZOOPLANKTON_CALC)THEN
	do jz = 1,nzp
      if(tgraze(k,i,jz) > 0.0)then
        lpzoooutp(k,i)=lpzoooutp(k,i) + zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))*zp(jz)
        lpzooinp(k,i)=lpzooinp(k,i) + zoo(k,i,jz)*zmu(k,i,jz)*prefp(jz)*LPOM(k,i)/tgraze(k,i,jz)*zp(jz)
      else
        lpzoooutp(k,i)=lpzoooutp(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))*zp(jz)
        lpzooinp(k,i)=0.0
      end if
    end do
    ENDIF
      LPOMpNS(K,I) = POMS(JW)*(LPOM(K-1,I)*orgplp(k-1,i)-LPOM(K,I)*orgplp(k,i))*BI(K,I)/BH2(K,I)
      LPOMpSS(K,I) = LPOMpAP(K,I)+lpompmp(k,i)-LPOMD(K,I)*orgplp(k,i)+LPOMpNS(K,I)-LRPOMD(K,I)*orgplp(k,i)
	  IF(ZOOPLANKTON_CALC)THEN
	  do jz = 1,nzp
	   LPOMpSS(K,I) =LPOMpSS(K,I) + lpzoooutp(k,i)-lpzooinp(k,i)
	  end do
	  ENDIF

	END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M   P H O S P H O R U S                            **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM_P
  rpompmp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          do jj=jt,je
            rpompmp(k,i)=rpompmp(k,i)+mpom(m)*(1.0-lrpmac(m))*mmr(k,i,m)*macrm(jj,k,i,m)*mp(m)
          end do
        end if
      end do
      rpompmp(k,i)=rpompmp(k,i)/(dlx(i)*bh(k,i))
      RPOMpNS(K,I) = POMS(JW)*(RPOM(K-1,I)*orgprp(k-1,i)-RPOM(K,I)*orgprp(k,i))*BI(K,I)/BH2(K,I)
      RPOMpSS(K,I) = LRPOMD(K,I)*orgplp(k,i)+RPOMpNS(K,I)-RPOMD(K,I)*orgprp(k,i)+rpompmp(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M   N I T R O G E N                                   **
!***********************************************************************************************************************************

ENTRY LABILE_DOM_n
  LDOMnAP(:,IU:ID) = 0.0; LDOMnEP(:,IU:ID) = 0.0; ldomnmp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LDOMnAP(K,I) = LDOMnAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)*an(ja)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LDOMnEP(K,I) = LDOMnEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)*en(je)
      END DO
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          if(k.eq.kt)then
            jt=kti(i)
          else
            jt=k
          end if
          je=kb(i)
          do jj=jt,je
            ldomnmp(k,i)=ldomnmp(k,i)+(1.0-mpom(m))*mmr(k,i,m)*macrm(jj,k,i,m)*mn(m)
          end do
        end if
      end do
      ldomnmp(k,i)=ldomnmp(k,i)/(dlx(i)*bh(k,i))
      LDOMnSS(K,I) = LDOMnAP(K,I)+LDOMnEP(K,I)+ldomnmp(k,i)-(LDOMD(K,I)+LRDOMD(K,I))*orgnld(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M   N I T R O G E N                                **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM_n
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMnSS(K,I) = LRDOMD(K,I)*orgnld(k,i)-RDOMD(K,I)*orgnrd(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M   N I T R O G E N                                    **
!***********************************************************************************************************************************

ENTRY LABILE_POM_n
  LPOMnAP(:,IU:ID) = 0.0;lpomnmp(:,iu:id)=0.0;lpzooinn(:,iu:id)=0.0; lpzoooutn(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LPOMnAP(K,I) = LPOMnAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))*an(ja)
      END DO
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          do jj=jt,je
            lpomnmp(k,i)=lpomnmp(k,i)+mpom(m)*lrpmac(m)*mmr(k,i,m)*macrm(jj,k,i,m)*mn(m)
          end do
        end if
      end do
      lpomnmp(k,i)=lpomnmp(k,i)/(dlx(i)*bh(k,i))
	IF(ZOOPLANKTON_CALC)THEN
	do jz = 1,nzp
      if(tgraze(k,i,jz) > 0.0)then
        lpzoooutn(k,i)=lpzoooutn(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))*zn(jz)
        lpzooinn(k,i)=lpzooinn(k,i)+zoo(k,i,jz)*prefp(jz)*zmu(k,i,jz)*LPOM(k,i)/tgraze(k,i,jz)*zn(jz)
      else
        lpzoooutn(k,i)=lpzoooutn(k,i)+zoo(k,i,jz)*(zmt(k,i,jz)+(zmu(k,i,jz)-(zmu(k,i,jz)*zeff(jz))))*zn(jz)
        lpzooinn(k,i)=0.0
      end if
	end do
	ENDIF
      LPOMnNS(K,I) = POMS(JW)*(LPOM(K-1,I)*orgnlp(k-1,i)-LPOM(K,I)*orgnlp(k,i))*BI(K,I)/BH2(K,I)
      LPOMnSS(K,I) = LPOMnAP(K,I)+lpomnmp(k,i)-LPOMD(K,I)*orgnlp(k,i)+LPOMnNS(K,I)-LRPOMD(K,I)*orgnlp(k,i) &
            + lpzoooutn(k,i)-lpzooinn(k,i)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M   N I T R O G E N                                **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM_n
  rpomnmp(:,iu:id)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      do m=1,nmc
        if(macrophyte_calc(jw,m))then
          jt=k
          je=kb(i)
          do jj=jt,je
            rpomnmp(k,i)=rpomnmp(k,i)+mpom(m)*(1.0-lrpmac(m))*mmr(k,i,m)*macrm(jj,k,i,m)*mn(m)
          end do
        end if
      end do
      rpomnmp(k,i)=rpomnmp(k,i)/(dlx(i)*bh(k,i))
      RPOMnNS(K,I) = POMS(JW)*(RPOM(K-1,I)*orgnrp(k-1,i)-RPOM(K,I)*orgnrp(k,i))*BI(K,I)/BH2(K,I)
      RPOMnSS(K,I) = LRPOMD(K,I)*orgnlp(k,i)+RPOMnNS(K,I)-RPOMD(K,I)*orgnrp(k,i)+rpomnmp(k,i)
    END DO
  END DO
RETURN


!************************************************************************
!**                          M A C R O P H Y T E                       **
!************************************************************************

ENTRY macrophyte(LLM)
  m=LLM
  DO I=IU,ID
    if(kticol(i))then
      jt=kti(i)
    else
      jt=kti(i)+1
    end if
    je=kb(i)
    do jj=jt,je
      if(jj.lt.kt)then
        colb=el(jj+1,i)
      else
        colb=el(kt+1,i)
      end if
      coldep=ELWS(i)-colb
      if(macrc(jj,Kt,I,m).gt.mmax(m))then
        mgr(jj,kt,i,m)=0.0
      end if
      macss(jj,Kt,I,m) = (mGR(jj,Kt,I,m)-mMR(Kt,I,m)-mRR(Kt,I,m))*macrc(jj,Kt,I,m)
      macrm(jj,kt,i,m)   = macrm(jj,kt,i,m)+macss(jj,kt,i,m)*dlt*coldep*cw(jj,i)*DLX(I)
    end do

    DO K=KT+1,KB(I)
      jt=k
      je=kb(i)
      do jj=jt,je
        if(macrc(jj,K,I,m).gt.mmax(m))then
          mgr(jj,k,i,m)=0.0
        end if
        macss(jj,K,I,m) = (mGR(jj,K,I,m)-mMR(K,I,m)-mRR(K,I,m))*macrc(jj,K,I,m)
        if(mact(jj,k,i).gt.mbmp(m).and.mact(jj,k-1,i).lt.mbmp(m).and.macss(jj,k,i,m).gt.0.0)then
          if(k-1.eq.kt)then
            bmass=macss(jj,k,i,m)*dlt*h2(k,i)*cw(jj,i)*DLX(I)
            macrm(jj,k-1,i,m)=macrm(jj,k-1,i,m)+bmass
            colb=el(kt+1,i)
            coldep=ELWS(i)-colb
            macss(jj,k-1,i,m)=bmass/dlt/(coldep*cw(jj,i)*DLX(I)) + macss(jj,k-1,i,m)
          else
            bmass=macss(jj,k,i,m)*dlt*h2(k,i)*cw(jj,i)*DLX(I)
            macrm(jj,k-1,i,m)=macrm(jj,k-1,i,m)+bmass
            macss(jj,k-1,i,m)=bmass/dlt/(h2(k-1,i)*cw(jj,i)*DLX(I))+ macss(jj,k-1,i,m)
          end if
          macss(jj,k,i,m)=0.0
        else
          bmasstest=macrm(jj,k,i,m)+macss(jj,k,i,m)*dlt*H2(k,i)*cw(jj,i)*DLX(I)
          if(bmasstest.ge.0.0)then
            macrm(jj,k,i,m)   = bmasstest
          else
            macss(jj,k,i,m)=-macrm(jj,k,i,m)/dlt/(H2(k,i)*cw(jj,i)*DLX(I))
            macrm(jj,k,i,m)=0.0
          end if
        end if
      end do
    END DO
  END DO
  DO I=IU,ID
    tmac=0.0
    cvol=0.0
    if(kticol(i))then
      jt=kti(i)
    else
      jt=kti(i)+1
    end if
    je=kb(i)

    do jj=jt,je
      if(jj.lt.kt)then
        colb=el(jj+1,i)
      else
        colb=el(kt+1,i)
      end if
      coldep=ELWS(i)-colb
      if(cw(jj,i).gt.0.0)then
        macrc(jj,kt,i,m)=macrm(jj,kt,i,m)/(cw(jj,i)*coldep*dlx(i))
      else
        macrc(jj,kt,i,m)=0.0
      end if
      tmac=tmac+macrm(jj,kt,i,m)
      cvol=cvol+cw(jj,i)*coldep*dlx(i)
    end do

    mac(kt,i,m)=tmac/cvol

    DO K=KT+1,KB(I)
      jt=k
      je=kb(i)
      tmac=0.0
      cvol=0.0
      do jj=jt,je
        if(cw(jj,i).gt.0.0)then
          macrc(jj,k,i,m)=macrm(jj,k,i,m)/(cw(jj,i)*h2(k,i)*dlx(i))
        else
          macrc(jj,k,i,m)=0.0
        end if
        tmac=tmac+macrm(jj,k,i,m)
        cvol=cvol+cw(jj,i)*h2(k,i)*dlx(i)
      end do
      mac(k,i,m)=tmac/cvol
    end do
  end do

  DO I=IU,ID
    tmac=0.0
    cvol=0.0
    do k=kt,kb(i)
      if(k.eq.kt)then
        jt=kti(i)
      else
        jt=k
      end if
      je=kb(i)
      do jj=jt,je
        mact(jj,k,i)=0.0
        do mi=1,nmc
          if(macrophyte_calc(jw,mi))then
            mact(jj,k,i)=macrc(jj,k,i,mi)+mact(jj,k,i)
          end if
        end do
      end do
    end do
  end do


  RETURN

! v3.5 end

!***********************************************************************************************************************************
!*                                                  K I N E T I C   F L U X E S                                                   **
!***********************************************************************************************************************************

ENTRY KINETIC_FLUXES
  DO JAF=1,NAF(JW)
    DO I=CUS(BS(JW)),DS(BE(JW))
      DO K=KT,KB(I)
        KFS(K,I,KFCN(JAF,JW)) = KFS(K,I,KFCN(JAF,JW))+KF(K,I,KFCN(JAF,JW))*VOL(K,I)*DLT
      END DO
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       p H   C O 2                                                             **
!***********************************************************************************************************************************

ENTRY PH_CO2

! pH and carbonate species

  DO I=IU,ID
    DO K=KT,KB(I)
      CART = TIC(K,I)/12000.0
      ALKT = ALK(K,I)/5.0E+04
      T1K  = T1(K,I)+273.15

!**** Ionic strength

      IF (FRESH_WATER(JW)) S2 = 2.5E-05*TDS(K,I)
      IF (SALT_WATER(JW))  S2 = 1.47E-3+1.9885E-2*TDS(K,I)+3.8E-5*TDS(K,I)*TDS(K,I)

!**** Debye-Huckel terms and activity coefficients

      SQRS2  =  SQRT(S2)
      DH1    = -0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03+4.160762E-02*S2-9.284843E-03*S2*S2
      DH2    = -2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02+9.715745E-02*S2-2.067746E-02*S2*S2
      H2CO3T =  10.0**(0.0755*S2)
      HCO3T  =  10.0**DH1
      CO3T   =  10.0**DH2
      OH     =  HCO3T

!**** Temperature adjustment

      KW = 10.0**(-283.971-0.05069842*T1K+13323.0/T1K+102.24447*LOG10(T1K)-1119669.0/(T1K*T1K))/OH
      K1 = 10.0**(-3404.71/T1K+14.8435-0.032786*T1K)*H2CO3T/HCO3T
      K2 = 10.0**(-2902.39/T1K+ 6.4980-0.023790*T1K)*HCO3T/CO3T

!**** pH evaluation

      PHT = -PH(K,I)-2.1
      IF (PH(K,I) <= 0.0) PHT = -14.0
      INCR = 10.0
      DO N=1,3
        F    = 1.0
        INCR = INCR/10.0
        ITER = 0
        DO WHILE (F > 0.0 .AND. ITER < 12)
          PHT    = PHT+INCR
          HION   = 10.0**PHT
          BICART = CART*K1*HION/(K1*HION+K1*K2+HION*HION)
          F      = BICART*(HION+2.0*K2)/HION+KW/HION-ALKT-HION/OH
          ITER   = ITER+1
        END DO
        PHT = PHT-INCR
      END DO

!**** pH, carbon dioxide, bicarbonate, and carbonate concentrations

      HION      =  10.0**PHT
      PH(K,I)   = -PHT
      CO2(K,I)  =  TIC(K,I)/(1.0+K1/HION+K1*K2/(HION*HION))
      HCO3(K,I) =  TIC(K,I)/(1.0+HION/K1+K2/HION)
      CO3(K,I)  =  TIC(K,I)/((HION*HION)/(K1*K2)+HION/K2+1.0)
    END DO
  END DO
RETURN

!**********************************************************
!**           SUBROUTINE ZOOPLANKTON                     **
!**********************************************************

ENTRY zooplankton
  DO I=IU,ID
    DO K=KT,KB(I)
	  do jz = 1, nzp
            zgztot=0.0                                                                                                   ! kv 5/9/2007
	        do jjz = 1,nzp
!             zgztot=zgztot+zgz(k,i,jz,jjz)*zoo(k,i,jz)                                                                   ! kv 5/9/2007
            zgztot=zgztot+zgz(k,i,jz,jjz)                                                                             ! cb 5/26/07
            end do
        zooss(k,i,jz)= (zmu(k,i,jz)*zeff(jz)-zrt(k,i,jz)-zmt(k,i,jz))*zoo(k,i,jz) - zgztot   ! omnivorous zooplankton    ! kv 5/9/2007
	  end do
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                              D E R I V E D   C O N S T I T U E N T S                                          **
!***********************************************************************************************************************************

ENTRY DERIVED_CONSTITUENTS
  APR = 0.0; ATOT = 0.0; TOTSS = 0.0; CHLA = 0.0; CBODU=0.0
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      DO I=CUS(JB),DS(JB)
        DO K=KT,KB(I)
          DO JA=1,NAL
            IF(ALG_CALC(JA))APR(K,I) = APR(K,I)+(AGR(K,I,JA)-ARR(K,I,JA))*ALG(K,I,JA)*H2(K,I)*DAY
          END DO
        END DO
        DO K=KT,KB(I)
          CBODC = 0.0; CBODN = 0.0; CBODP = 0.0; BODTOT = 0.0; ALGP = 0.0; ALGN = 0.0
          DO JA=1,NAL
            IF(ALG_CALC(JA))ATOT(K,I) = ATOT(K,I)+ALG(K,I,JA)
          END DO
          DO IBOD=1,NBOD
          IF(BOD_CALC(IBOD))THEN
            CBODC  = CBODC+CBOD(K,I,IBOD)*BODC(IBOD)
            CBODN  = CBODN+CBOD(K,I,IBOD)*BODN(IBOD)
            CBODP  = CBODP+CBOD(K,I,IBOD)*BODP(IBOD)
            BODTOT = BODTOT+CBOD(K,I,IBOD)
          ENDIF
          END DO
          DOM(K,I) = LDOM(K,I)+RDOM(K,I)
          POM(K,I) = LPOM(K,I)+RPOM(K,I)
          DOC(K,I) = DOM(K,I)*ORGC(JW)+CBODC
          POC(K,I) = POM(K,I)*ORGC(JW)
          DO JA=1,NAL
          IF(ALG_CALC(JA))THEN
            POC(K,I) = POC(K,I)+ALG(K,I,JA)*AC(JA)
            ALGP     = ALGP+ALG(K,I,JA)*AP(JA)
            ALGN     = ALGN+ALG(K,I,JA)*AN(JA)
          ENDIF
          END DO
! v3.5 start
          IF(ZOOPLANKTON_CALC)THEN
            do jz=1,nzp
                poc(k,i)=poc(k,i)+zc(jz)*zoo(k,i,jz) !mlm baulk
                zoop=zoo(k,i,jz)*zp(jz) !mlm baulk
                zoon=zoo(k,i,jz)*zn(jz) !mlm baulk
                cbodu(k,i) = cbodu(k,i) + o2om(jw)*zoo(k,i,jz)
	        end do
	      ENDIF
          TOC(K,I)   = DOC(K,I)+POC(K,I)
          DOP(K,I)   = LDOM(K,I)*ORGPLD(k,i)+RDOM(k,i)*orgprd(k,i)+CBODP
          DON(K,I)   = LDOM(K,I)*ORGNLD(k,i)+RDOM(k,i)*orgnrd(k,i)+CBODN
          POP(K,I)   = LPOM(K,I)*ORGPLP(k,i)+RPOM(k,i)*orgprp(k,i)+ALGP+zoop
          PON(K,I)   = LPOM(K,I)*ORGNLP(k,i)+RPOM(k,i)*orgnrp(k,i)+ALGN+zoop
! v3.5 end
          TOP(K,I)   = DOP(K,I)+POP(K,I)
          TON(K,I)   = DON(K,I)+PON(K,I)
          TKN(K,I)   = TON(K,I)+NH4(K,I)
! v3.5 start
          CBODU(K,I) = CBODU(K,I)+O2OM(JW)*(DOM(K,I)+POM(K,I)+ATOT(K,I))+BODTOT
! v3.5 end
          TPSS       = 0.0
          DO JS=1,NSS
            TPSS = TPSS+SS(K,I,JS)*PARTP(JW)
          END DO
          TP(K,I)   =  TOP(K,I)+PO4(K,I)+TPSS
          TN(K,I)   =  TON(K,I)+NH4(K,I)+NO3(K,I)
          O2DG(K,I) = (O2(K,I)/SATO(T1(K,I),TDS(K,I),PALT(I),SALT_WATER(JW)))*100.0
          DO JA=1,NAL
          IF(ALG_CALC(JA))THEN
            CHLA(K,I)  = CHLA(K,I) +ALG(K,I,JA)/ACHLA(JA)
            TOTSS(K,I) = TOTSS(K,I)+ALG(K,I,JA)
          ENDIF
          END DO
          TOTSS(K,I) = TOTSS(K,I)+TISS(K,I)+POM(K,I)
        END DO
      END DO
    END DO
  END DO
RETURN
ENTRY DEALLOCATE_KINETICS
  DEALLOCATE (OMTRM,  SODTRM, NH4TRM, NO3TRM, DOM, POM, PO4BOD, NH4BOD, TICBOD, ATRM,   ATRMR,  ATRMF, ETRM,   ETRMR,  ETRMF, BIBH2)
! v3.5 start
  deALLOCATE (lam2m)
! v3.5 end
RETURN
END SUBROUTINE KINETICS
