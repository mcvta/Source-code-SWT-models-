SUBROUTINE ENDSIMULATION

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  !v3.5
  EXTERNAL RESTART_OUTPUT

!***********************************************************************************************************************************
!*                                                    Task 3: End Simulation                                                      **
!***********************************************************************************************************************************

  CALL DATE_AND_TIME (CDATE,CCTIME)
  IF (.NOT. ERROR_OPEN) TEXT = 'Normal termination at '//CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)//' on '//CDATE(5:6)//'/'     &
                                                       //CDATE(7:8)//'/'//CDATE(3:4)
  TEXT = ADJUSTL (TRIM(TEXT))
  CALL CPU_TIME (CURRENT)
  DO JW=1,NWB
    IF (SNAPSHOT(JW)) THEN
      WRITE (SNP(JW),'(/A/)')            ADJUSTL(TRIM(TEXT))
      WRITE (SNP(JW),'(A)')             'Runtime statistics'
      WRITE (SNP(JW),'(2(A,I0))')       '  Grid                 = ', IMX,' x ',KMX
      WRITE (SNP(JW),'(A,I0)')          '  Maximum active cells = ', NTACMX,'  Minimum active cells = ',NTACMN
      WRITE (SNP(JW),'(3(A,F0.1))')     '  Segment lengths, m   = ', DLXMIN,'-',DLXMAX
      WRITE (SNP(JW),'(3(A,F0.1))')     '  Layer heights, m     = ', HMIN,  '-',HMAX
      WRITE (SNP(JW),'(A)')             '  Timestep'
      WRITE (SNP(JW),'(A,I0)')          '    Total iterations   = ', NIT
      WRITE (SNP(JW),'(A,I0)')          '    # of violations    = ', NV
      WRITE (SNP(JW),'(A,F0.2)')        '    % violations       = ', FLOAT(NV)/FLOAT(NIT)*100.0
      WRITE (SNP(JW),'(A,I0,A)')        '    Average timestep   = ', INT(DLTAV),' sec'
      WRITE (SNP(JW),'(A,I0,A,F0.2,A)') '  Simulation time      = ', INT(ELTMJD),' days ',(ELTMJD-INT(ELTMJD))*24.0,' hours'
      WRITE (SNP(JW),'(A,F0.2,A)')      '  Total CPU runtime    = ',(CURRENT-START)/60.0,' min'
      CLOSE (SNP(JW))
    END IF
    IF (VECTOR(JW))      CLOSE (VPL(JW))
    IF (PROFILE(JW))     CLOSE (PRF(JW))
    IF (SPREADSHEET(JW)) CLOSE (SPR(JW))
    IF (CONTOUR(JW))     CLOSE (CPL(JW))
  END DO
  IF (TIME_SERIES) THEN
    DO J=1,NIKTSR
      CLOSE (TSR(J))
    END DO
  END IF
  IF (WARNING_OPEN) THEN
    CLOSE (WRN)
  ELSE
    CLOSE (WRN,STATUS='DELETE')
  END IF
  IF (ERROR_OPEN) THEN
    CLOSE (W2ERR)
  ELSE
    CLOSE (W2ERR,STATUS='DELETE')
  END IF
  DO J=40,NOPEN
   CLOSE (J)
  END DO

  IF(ERROR_OPEN)THEN
  OPEN(W2ERR,FILE='W2Errordump.opt',status='unknown')
  WRITE(w2err,*)'JDAY',jday,'SZ',sz,'Z',z,'H2KT',h2(kt,1:imx),'H1KT',h1(kt,1:imx),'WSE',elws,'Q',q,'QC',qc,'T1',t1(kt,1:imx),'T2',t2(kt,1:imx),'SUKT',su(kt,1:imx),&
                          'UKT',u(kt,1:imx),'QIN',qin,'QTR',qtr,'QWD',qwd
  close(w2err)
  ENDIF

  DEALLOCATE (LAYERCHANGE,TECPLOT,X1)
  DEALLOCATE (TSR,    WDO,    ETSR,   IWDO,   ITSR,   TITLE,  CDAC,   WSC,    ESTR,   WSTR,   QSTR,   KTSW,   KBSW,   SINKC)
  DEALLOCATE (EBC,    MBC,    PQC,    EVC,    PRC,    WINDC,  QINC,   QOUTC,  HEATC,  SLHTC,  QINIC,  DTRIC,  TRIC,   WDIC)
  DEALLOCATE (EXC,    EXIC,   VBC,    METIC,  SLTRC,  THETA,  FRICC,  NAF,    ELTMF,  ZMIN,   IZMIN,  C2CH,   CDCH,   EPCH, KFCH)
  DEALLOCATE (CPLTC,  HPLTC,  CMIN,   CMAX,   HYMIN,  HYMAX,  CDMIN,  CDMAX,  JBDAM,  ILAT,   CDPLTC, QINSUM, TINSUM, TIND)
  DEALLOCATE (QOLD,   DTP,    DTPS,   QOLDS,  QIND,   JSS,    HDIC,   QNEW,   YSS,    VSS,    YS,     VS,     VSTS,   NSPRF)
  DEALLOCATE (LATGTC, LATSPC, LATPIC, DYNPIPE,LATPUC, DYNGTC, OPT,    CIND,   CINSUM, CDWBC,  KFWBC,  CPRWBC, CINBRC, CTRTRC, CDTBRC)   ! SW 5/10/10
  DEALLOCATE (YSTS,   YST,    VST,    ALLIM,  APLIM,  ANLIM,  ASLIM,  ELLIM,  EPLIM,  ENLIM,  ESLIM,  CSSK,   C1,     C2, Z0)
  DEALLOCATE (KFS,    AF,     EF,     HYD,    KF,     AZSLC,  STRIC,  CPRBRC, CD,     KBMAX,  ELKT,   WIND2,  VISC,   CELC)
  DEALLOCATE (QOAVR,  QIMXR,  QOMXR,  REAERC, LAT,    LONGIT, ELBOT,  BTH,    VPR,    LPR,    NISNP,  NIPRF,  NISPR,  DECL)
  DEALLOCATE (A00,    HH,     T2I,    KTWB,   KBR,    IBPR,   DLVR,   ESR,    ETR,    NBL,    LPRFN,  EXTFN,  BTHFN,  METFN)
  DEALLOCATE (SNPFN,  PRFFN,  SPRFN,  CPLFN,  VPLFN,  FLXFN,  FLXFN2, VPRFN,  AFW,    BFW,    CFW,    WINDH,  RHEVC,  FETCHC, JBDN)
  DEALLOCATE (KBI, macch, gridc, gma, bta, qtot, sedcip, sedcin, sedcic,sedcis)            ! SW 9/27/2007
! v3.5 start
  DEALLOCATE (SDK,    FSOD,   FSED,   SEDCI,  SEDCc,   SEDPRC, ICEC,   SLICEC, ICETHI, ALBEDO, HWI,    BETAI,  GAMMAI,ICEMIN)
  DEALLOCATE (seds,   sedb)    !cb 11/28/06
! v3.5 end
  DEALLOCATE (EXH2O,  BETA,   EXOM,   EXSS,   DXI,    CBHE,   TSED,   TSEDF,  FI,     ICET2,  AZC,    AZMAX,  QINT,   QOUTT)
  DEALLOCATE (AX,     WTYPEC, TAIR,   TDEW,   WIND,   PHI,    CLOUD,  CSHE,   SRON,   RANLW,    RB,     RC,     RE,     SHADE)
  DEALLOCATE (ET,     RS,     RN,     SNPC,   SCRC,   PRFC,   SPRC,   CPLC,   VPLC,   FLXC,   NXTMCP, NXTMVP, NXTMFL, GAMMA)
  DEALLOCATE (NXTMSN, NXTMSC, NXTMPR, NXTMSP, SNPDP,  SCRDP,  PRFDP,  SPRDP,  CPLDP,  VPLDP,  FLXDP,  NCPL,   NVPL,   NFLX)
  DEALLOCATE (NSNP,   NSCR,   NPRF,   NSPR,   NEQN,   PO4R,   PARTP,  NH4DK,  NH4R,   NO3DK,  NO3S,   FER,    FES,    CDSUM)
  DEALLOCATE (CO2R,   SROC,   O2ER,   O2EG,   CAQ10,  CADK,   CAS,    BODP,   BODN,   BODC,   KBOD,   TBOD,   RBOD,   DTRC)
  DEALLOCATE (LDOMDK, RDOMDK, LRDDK,  OMT1,   OMT2,   OMK1,   OMK2,   LPOMDK, RPOMDK, LRPDK,  POMS,   ORGP,   ORGN,   ORGC)
  DEALLOCATE (RCOEF1, RCOEF2, RCOEF3, RCOEF4, ORGSI,  NH4T1,  NH4T2,  NH4K1,  NH4K2,  NO3T1,  NO3T2,  NO3K1,  NO3K2,  NSTR)
  DEALLOCATE (DSIR,   PSIS,   PSIDK,  PARTSI, SODT1,  SODT2,  SODK1,  SODK2,  O2NH4,  O2OM,   O2AR,   O2AG,   CG1DK,  CGS)
  DEALLOCATE (CGQ10,  CG0DK,  CUNIT,  CUNIT1, CUNIT2, CAC,    INCAC,  TRCAC,  DTCAC,  PRCAC,  CNAME,  CNAME1, CNAME2, CMULT)
  DEALLOCATE (CN,     INCN,   DTCN,   PRCN,   CSUM,   DLTMAX, QWDO,   TWDO,   SSS,    SEDRC,  TAUCR,  XBR, FNO3SED)
  DEALLOCATE (QTAVB,  QTMXB,  BS,     BE,     JBUH,   JBDH,   TSSS,   TSSB,   TSSICE, ESBR,   ETBR,   EBRI,   QDTR,   EVBR)
  DEALLOCATE (QIN,    PR,     QPRBR,  TIN,    TOUT,   TPR,    TDTR,   TPB,    NACPR,  NACIN,  NACDT,  NACTR,  NACD,   ELDH)
  DEALLOCATE (QSUM,   NOUT,   KTQIN,  KBQIN,  ELUH,   NL,     NPOINT, SLOPE,  ALPHA,  COSA,   SINA,   TDHFN,  QOTFN,  PREFN)
  DEALLOCATE (CPRFN,  EUHFN,  TUHFN,  CUHFN,  EDHFN,  QINFN,  TINFN,  CINFN,  CDHFN,  QDTFN,  TDTFN,  CDTFN,  TPRFN,  VOLEV)
  DEALLOCATE (VOLWD,  VOLSBR, VOLTBR, DLVOL,  VOLG,   VOLSR,  VOLTR,  VOLB,   VOLPR,  VOLTRB, VOLDT,  VOLUH,  VOLDH,  VOLIN)
  DEALLOCATE (US,     DS,     CUS,    UHS,    DHS,    UQB,    DQB,    CDHS,   VOLOUT, TSSWD,  TSSUH,  TSSDH,  TSSIN,  TSSOUT)
  DEALLOCATE (TSSEV,  TSSPR,  TSSTR,  TSSDT,  SOD,    ELWS,   BKT,    REAER,  ICETH,  ICE,    ICESW,  Q,      QC,     QERR)
  DEALLOCATE (KTI,    SROSH,  SEG,    DLXRHO, QSSUM,  DLX,    DLXR,   QUH1,   QDH1,   BI,     JWUH,   JWDH)
  DEALLOCATE (A,      C,      D,      F,      V,      SKTI,   KBMIN,  EV,     QDT,    QPR,    SBKT,   BHRHO)
  DEALLOCATE (SZ,     WSHX,   WSHY,   WIND10, CZ,     FETCH,  PHI0,   FRIC,   ADZ,    HMULT,  FMTC,   FMTCD,  CNAME3, CDNAME3)
  DEALLOCATE (Z,      KB,     VNORM,  ANPR,   ANEQN,  APOM,   ACHLA,  AHSP,   AHSN,   AHSSI)
  DEALLOCATE (AC,     ASI,    AT1,    AT2,    AT3,    AT4,    AK1,    AK2,    AK3,    AK4,    EXA,    ASAT,   AP,     AN)
  DEALLOCATE (AG,     AR,     AE,     AM,     AS,     ENPR,   ENEQN,  EG,     ER,     EE,     EM,     EB,     ESAT,   EP)
  DEALLOCATE (EC,     ESI,    ECHLA,  EHSP,   EHSN,   EHSSI,  EPOM,   EHS,    EN,     ET4,    EK1,    EK2,    EK3,    EK4)
  DEALLOCATE (ET1,    ET2,    ET3,    HNAME,  FMTH,    KFAC,  KFNAME, KFNAME2,KFCN,   C2I,    TRCN,   CDN,    CDNAME, CDNAME2,CDMULT)
  DEALLOCATE (CMBRS,  CMBRT,  FETCHU, FETCHD, IPRF,   ISNP,   ISPR,   BL,     LFPR,   DO3,    SED,    TKE,    PALT)
  DEALLOCATE (ADX,    DO1,    DO2,    B,      CONV,   CONV1,  EL,     DZ,     DZQ,    DX,     SAZ,    T1,TSS,QSS,BNEW, ILAYER)   ! SW 1/23/06
  DEALLOCATE (P,      SU,     SW,     BB,     BR,     BH,     BHR,    VOL,    HSEG,   DECAY,  FPFE,   FRICBR, UXBR,   UYBR)
  DEALLOCATE (DEPTHB, DEPTHM, FPSS,   TUH,    TDH,    TSSUH1, TSSUH2, TSSDH1, TSSDH2, SEDVP,  H,      EPC)
  DEALLOCATE (TVP,    QINF,   QOUT,   KOUT,   VOLUH2, VOLDH2, CWDO,   CDWDO,  CWDOC,  CDWDOC, CDTOT,  CPR,    CPB,    COUT)
  DEALLOCATE (CIN,    CDTR,   RSOD,   RSOF,   DLTD,   DLTF,   TSRD,   TSRF,   WDOD,   WDOF,   SNPD,   SNPF,   SPRD,   SPRF)
  DEALLOCATE (SCRD,   SCRF,   PRFD,   PRFF,   CPLD,   CPLF,   VPLD,   VPLF,   FLXD,   FLXF,   EPIC,   EPICI,  EPIPRC, EPIVP)
  DEALLOCATE (CUH,    CDH,    EPM,    EPD,    C1S,    CSSB,   CVP,    CSSUH1, CSSUH2, CSSDH2, CSSDH1, LNAME,  IWR,    KTWR)
  DEALLOCATE (JWUSP,  JWDSP,  QSP,    KBWR,   KTWD,   KBWD,   JBWD,   GTA1,   GTB1,   GTA2,   GTB2,   BGT,    IUGT,   IDGT)
  DEALLOCATE (QTR,    TTR,    KTTR,   KBTR,   EGT,    AGASGT, BGASGT, CGASGT, GASGTC, PUGTC,  ETUGT,  EBUGT,  KTUGT,  KBUGT)
  DEALLOCATE (PDGTC,  ETDGT,  EBDGT,  KTDGT,  KBDGT,  A1GT,   B1GT,   G1GT,   A2GT,   B2GT,   G2GT,   JWUGT,  JWDGT,  QGT)
  DEALLOCATE (EQGT,   JBUGT,  JBDGT,  JBUPI,  JBDPI,  JWUPI,  JWDPI,  QPI,    IUPI,   IDPI,   EUPI,   EDPI,   WPI,    DLXPI, BP)   ! SW 5/5/10
  DEALLOCATE (ETUPI,  EBUPI,  KTUPI,  KBUPI,  PDPIC,  ETDPI,  EBDPI,  KTDPI,  KBDPI,  FPI,    FMINPI, PUPIC,  ETDSP,  EBDSP)
  DEALLOCATE (PUSPC,  ETUSP,  EBUSP,  KTUSP,  KBUSP,  PDSPC,  KTDSP,  KBDSP,  IUSP,   IDSP,   ESP,    A1SP,   B1SP,   A2SP)
  DEALLOCATE (B2SP,   AGASSP, BGASSP, CGASSP, EQSP,   GASSPC, JBUSP,  JBDSP,  STRTPU, ENDPU,  EONPU,  EOFFPU, QPU,    PPUC)
  DEALLOCATE (IUPU,   IDPU,   EPU,    ETPU,   EBPU,   KTPU,   KBPU,   JWUPU,  JWDPU,  JBUPU,  JBDPU,  PUMPON, KTW,    KBW)
  DEALLOCATE (IWD,    KWD,    QWD,    EWD,    ITR,    QTRFN,  TTRFN,  CTRFN,  ELTRT,  ELTRB,  TRC,    JBTR,   QTRF,   CLRB)
  DEALLOCATE (TTLB,   TTRB,   CLLB,   SRLB1,  SRRB1,  SRLB2,  SRRB2,  SRFJD1, SHADEI, SRFJD2, TOPO,   QSW,    CTR)    ! SW 10/17/05
  DEALLOCATE (H1,     H2,     BH1,    BH2,    BHR1,   BHR2,   AVH1,   AVH2,   SAVH2,  AVHR,   SAVHR,  CBODD)
  DEALLOCATE (OPEN_VPR,       OPEN_LPR,             POINT_SINK,         HPRWBC,           READ_EXTINCTION, READ_RADIATION)
  DEALLOCATE (DIST_TRIBS,     UPWIND,               ULTIMATE,           FRESH_WATER,      SALT_WATER,      LIMITING_FACTOR)
  DEALLOCATE (UH_EXTERNAL,    DH_EXTERNAL,          UH_INTERNAL,        DH_INTERNAL,      UQ_INTERNAL,     DQ_INTERNAL)
  DEALLOCATE (UQ_EXTERNAL,    DQ_EXTERNAL,          UP_FLOW,            DN_FLOW,          UP_HEAD,         DN_HEAD)
  DEALLOCATE (INTERNAL_FLOW,  DAM_INFLOW,           DAM_OUTFLOW,        HEAD_FLOW,        HEAD_BOUNDARY)      !TC 08/03/04
  DEALLOCATE (ISO_CONC,             VERT_CONC,          LONG_CONC,        VERT_SEDIMENT,   LONG_SEDIMENT)
  DEALLOCATE (ISO_SEDIMENT,   VISCOSITY_LIMIT,      CELERITY_LIMIT,     IMPLICIT_AZ,      ONE_LAYER,       IMPLICIT_VISC)
  DEALLOCATE (FETCH_CALC,     LIMITING_DLT,         TERM_BY_TERM,       MANNINGS_N,       PLACE_QTR,       SPECIFY_QTR)
  DEALLOCATE (PLACE_QIN,      PRINT_CONST,          PRINT_HYDRO,        PRINT_SEDIMENT,   ENERGY_BALANCE,  MASS_BALANCE)
  DEALLOCATE (VOLUME_BALANCE, DETAILED_ICE,         ICE_CALC,           ICE_IN,           ALLOW_ICE,       PH_CALC)
  DEALLOCATE (EVAPORATION,    PRECIPITATION,        RH_EVAP,            NO_INFLOW,        NO_OUTFLOW,      NO_HEAT)
  DEALLOCATE (ISO_TEMP,       VERT_TEMP,            LONG_TEMP,          VERT_PROFILE,     LONG_PROFILE,    NO_WIND)
  DEALLOCATE (SNAPSHOT,       PROFILE,              VECTOR,             CONTOUR,          SPREADSHEET,     INTERNAL_WEIR)
  DEALLOCATE (SCREEN_OUTPUT,  FLUX,                 DYNAMIC_SHADE,      TRAPEZOIDAL, BOD_CALC, ALG_CALC)
  DEALLOCATE (SEDIMENT_CALC,  EPIPHYTON_CALC,       PRINT_DERIVED,      PRINT_EPIPHYTON,  TDG_SPILLWAY,    TDG_GATE, DYNSEDK)
  DEALLOCATE (ISO_EPIPHYTON,  VERT_EPIPHYTON,       LONG_EPIPHYTON,     LATERAL_SPILLWAY, LATERAL_GATE,    LATERAL_PUMP)
  DEALLOCATE (INTERP_HEAD,    INTERP_WITHDRAWAL,    INTERP_EXTINCTION,  INTERP_DTRIBS,    LATERAL_PIPE,    INTERP_TRIBS)
  DEALLOCATE (INTERP_OUTFLOW, INTERP_INFLOW,        INTERP_METEOROLOGY, CONSTITUENT_PLOT, DERIVED_PLOT,    ZERO_SLOPE)
  DEALLOCATE (HYDRO_PLOT,     SEDIMENT_RESUSPENSION)
! v3.5 start
  deallocate (ORGPLD, ORGPRD, ORGPLP, ORGPRP, ORGNLD, ORGNRD, ORGNLP)
  deALLOCATE (icpl,tavg,tavgw)
  deallocate (ORGNRP)
  deallocate  (print_macrophyte, macrophyte_calc,macwbc,conv2)
  deallocate  (mac, macrc,mact, macrm, macss)
  deallocate  (mgr,mmr, mrr)
  deallocate  (smacrc, smacrm)
  deallocate  (smact, smac)
  deallocate  (mt1,mt2,mt3,mt4,mk1,mk2,mk3,mk4,mg,mr,mm)
  deallocate  (mp, mn, mc,psed,nsed,mhsp,mhsn,mhsc,msat)
  deallocate  (cddrag,kticol,armac,macwbci, anorm, dwv, dwsa)
  deallocate  (mbmp,mmax,mpom,lrpmac,o2mr,o2mg)
  deallocate  (macmbrs, macmbrt,ssmacmb)
  deallocate  (cw, bic)
  deallocate  (mactrmr, mactrmf,mactrm)
  deallocate  (mlfpr)
  deallocate  (mllim, mplim,mclim,mnlim)
  deALLOCATE  (GAMMAj)
  deallocate (por,VOLKTi,VOLi,vstem,vstemkt,sarea)
  DEALLOCATE (IWIND) ! MLM 08/12/05
  deallocate (zg,zm,zeff,prefp,zr,zoomin,zs2p,exz,zt1,zt2,zt3,zt4,zk1,zk2)
  deallocate (ldompmp,ldomnmp,lpompmp,lpomnmp,rpompmp,rpomnmp,o2zr) ! mlm 06/10/06
  deallocate (mprwbc)                                               ! mlm 06/10/06
  deallocate (exm)                                                  ! mlm 06/10/06
  DEALLOCATE (USTARBTKE,E,EROUGH, ARODI, STRICK, TKELATPRDCONST,AZT,DZT)
  DEALLOCATE(FIRSTI, LASTI, TKELATPRD, STRICKON, WALLPNT, IMPTKE, TKEBC)
! deallocate (sedbr,sedbrp,sedbrn,sedbrc)   ! SW 6/4/07 No need to deallocate pointers
  deallocate (zk3,zk4,zp,zn,zc,prefa,zmu,tgraze,zrt,zmt,zoorm,zoormr,zoormf) ! POINTERS ,zoo,zooss,
  deallocate (lpzooout,lpzooin,po4zr,nh4zr,dozr,ticzr,agz,agzt)
  deallocate (zgz,PREFZ) !omnivorous zooplankton
!  DEALLOCATE (LDOMPSS,  RDOMPSS,  LPOMPSS,  RPOMPSS,  LDOMNSS,  RDOMNSS)
!  DEALLOCATE (LPOMNSS,  RPOMNSS)
!  DEALLOCATE (LDOMPAP, LDOMPEP, LPOMPAP, LPOMPNS, RPOMPNS)
!  DEALLOCATE (LDOMnAP, LDOMnEP, LPOMnAP, LPOMnNS, RPOMnNS)
!  DEALLOCATE (ldompmp,ldomnmp,lpompmp,lpomnmp,rpompmp,rpomnmp)
  DEALLOCATE (lpzooinp,lpzooinn,lpzoooutp,lpzoooutn)
!  DEALLOCATE (SEDDp,  SEDASp,  SEDOMSp, SEDNSp, lpomepp)
!  DEALLOCATE (SEDDn,  SEDASn,  SEDOMSn, SEDNSn, lpomepn, sedno3)
!  DEALLOCATE (SEDDc,  SEDASc,  SEDOMSc, SEDNSc, lpomepc)
  DEALLOCATE (sedc, sedn, sedp)
  DEALLOCATE (sedvpc, sedvpp, sedvpn)
  DEALLOCATE (sdkv,seddktot)
  DEALLOCATE (cbods,KFJW)
! v3.5 end

  CALL DEALLOCATE_TIME_VARYING_DATA
  CALL DEALLOCATE_TRANSPORT
  CALL DEALLOCATE_KINETICS
  CALL DEALLOCATE_WATERBODY
  CALL DEALLOCATE_PIPE_FLOW
  CALL DEALLOCATE_OPEN_CHANNEL

  RETURN
  END SUBROUTINE ENDSIMULATION
