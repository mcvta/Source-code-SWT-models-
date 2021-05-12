subroutine input

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  !v3.5
  EXTERNAL RESTART_OUTPUT

  integer nproc                                                             ! SW 7/13/09
  character*1 char1
  character*8 AID

! Title and array dimensions

  ALLOCATE (TITLE(11))
  READ (CON,'(///(8X,A72))') (TITLE(J),J=1,10)
  READ (CON,'(//8X,5I8,2A8)')     NWB, NBR, IMX, KMX, NPROC, CLOSEC,SELECTC              ! SW 7/31/09
  READ (CON,'(//8X,8I8)')     NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU
  READ (CON,'(//8X,7I8,a8)')  NGC, NSS, NAL, NEP, NBOD, nmc, nzp
  READ (CON,'(//8X,I8)')      NOD

  if(NPROC == 0)NPROC=2                                                                 ! SW 7/31/09
  call omp_set_num_threads(NPROC)   ! set # of processors to NPROC  Moved to INPUT subroutine  TOGGLE FOR DEBUG
  if(SELECTC=='        ')then
     SELECTC='     OFF'
  endif


! Constituent numbers

  NTDS  = 1
  NGCS  = 2
  NGCE  = NGCS+NGC-1
  NSSS  = NGCE+1
  NSSE  = NSSS+NSS-1
  NPO4  = NSSE+1
  NNH4  = NPO4+1
  NNO3  = NNH4+1
  NDSI  = NNO3+1
  NPSI  = NDSI+1
  NFE   = NPSI+1
  NLDOM = NFE+1
  NRDOM = NLDOM+1
  NLPOM = NRDOM+1
  NRPOM = NLPOM+1
  NBODS = NRPOM+1
  NBODE = NBODS+NBOD-1
  NAS   = NBODE+1
  NAE   = NAS+NAL-1
  NDO   = NAE+1
  NTIC  = NDO+1
  NALK  = NTIC+1

! v3.5 start
  NZOOS = NALK + 1
  NZOOE = NZOOS + NZP - 1
  NLDOMP=NZOOE+1
  NRDOMP=nldomp+1
  NLPOMP=nrdomp+1
  NRPOMP=nlpomp+1
  NLDOMn=nrpomp+1
  NRDOMn=nldomn+1
  NLPOMn=nrdomn+1
  NRPOMn=nlpomn+1
  nct=nrpomn
! v3.5 end

! Constituent, tributary, and widthdrawal totals

  NTRT = NTR+NGT+NSP+NPI+NPU
  NWDT = NWD+NGT+NSP+NPI+NPU
  NEPT = MAX(NEP,1)
  Nmct = MAX(nmc,1)
  nzpt=max(nzp,1)
  ALLOCATE (CDAC(NDC), X1(IMX), TECPLOT(NWB))
  ALLOCATE (WSC(IMX),    KBI(IMX))
  ALLOCATE (VBC(NWB),    EBC(NWB),    MBC(NWB),    PQC(NWB),    EVC(NWB),    PRC(NWB))
  ALLOCATE (WINDC(NWB),  QINC(NWB),   QOUTC(NWB),  HEATC(NWB),  SLHTC(NWB))
  ALLOCATE (QINIC(NBR),  DTRIC(NBR),  TRIC(NTR),   WDIC(NWD),   HDIC(NBR),   METIC(NWB))
  ALLOCATE (EXC(NWB),    EXIC(NWB))
  ALLOCATE (SLTRC(NWB),  THETA(NWB),  FRICC(NWB),  NAF(NWB),    ELTMF(NWB), Z0(NWB))
  ALLOCATE (ZMIN(NWB),   IZMIN(NWB))
  ALLOCATE (C2CH(NCT),   CDCH(NDC),   EPCH(NEPT),  macch(nmct), KFCH(NFL))  !v3.5
  ALLOCATE (CPLTC(NCT),  HPLTC(NHY),  CDPLTC(NDC))
  ALLOCATE (CMIN(NCT),   CMAX(NCT),   HYMIN(NHY),  HYMAX(NHY),  CDMIN(NDC),  CDMAX(NDC))
  ALLOCATE (JBDAM(NBR),  ILAT(NWDT))
  ALLOCATE (QINSUM(NBR), TINSUM(NBR), TIND(NBR),   JSS(NBR),    QIND(NBR))
  ALLOCATE (QOLD(NPI),   DTP(NPI),    DTPS(NPI),   QOLDS(NPI))
  ALLOCATE (LATGTC(NGT), LATSPC(NSP), LATPIC(NPI), DYNPIPE(NPI),LATPUC(NPU), DYNGTC(NGT))
  ALLOCATE (OPT(NWB,7),         CIND(NCT,NBR),         CINSUM(NCT,NBR))
  ALLOCATE (CDWBC(NDC,NWB),     KFWBC(NFL,NWB),        CPRWBC(NCT,NWB),    CINBRC(NCT,NBR),     CTRTRC(NCT,NTR))
  ALLOCATE (CDTBRC(NCT,NBR),    CPRBRC(NCT,NBR))
  ALLOCATE (YSS(NNPIPE,NPI),    VSS(NNPIPE,NPI),       YS(NNPIPE,NPI),     VS(NNPIPE,NPI),      VSTS(NNPIPE,NPI))
  ALLOCATE (YSTS(NNPIPE,NPI),   YST(NNPIPE,NPI),       VST(NNPIPE,NPI))
  ALLOCATE (CBODD(KMX,IMX,NBOD))
  ALLOCATE (ALLIM(KMX,IMX,NAL), APLIM(KMX,IMX,NAL),    ANLIM(KMX,IMX,NAL), ASLIM(KMX,IMX,NAL))
  ALLOCATE (ELLIM(KMX,IMX,NEP), EPLIM(KMX,IMX,NEP),    ENLIM(KMX,IMX,NEP), ESLIM(KMX,IMX,NEP))
  ALLOCATE (CSSK(KMX,IMX,NCT),  C1(KMX,IMX,NCT),       C2(KMX,IMX,NCT),    CD(KMX,IMX,NDC),     KF(KMX,IMX,NFL))
  ALLOCATE (KFS(KMX,IMX,NFL),   AF(KMX,IMX,NAL,5),     EF(KMX,IMX,NEP,5),  HYD(KMX,IMX,NHY), KFJW(NWB,NFL))
  ALLOCATE (TKE(KMX,IMX,3), AZT(KMX,IMX),DZT(KMX,IMX))
  ALLOCATE (USTARBTKE(IMX),E(IMX),EROUGH(NWB), ARODI(NWB), STRICK(NWB), TKELATPRDCONST(NWB))
  ALLOCATE (FIRSTI(NWB), LASTI(NWB), TKELATPRD(NWB), STRICKON(NWB), WALLPNT(NWB), IMPTKE(NWB), TKEBC(NWB))
  ALLOCATE (HYDRO_PLOT(NHY),    CONSTITUENT_PLOT(NCT), DERIVED_PLOT(NDC))
  ALLOCATE (ZERO_SLOPE(NWB),    DYNAMIC_SHADE(IMX))
  ALLOCATE (AZSLC(NWB))
  ALLOCATE (NSPRF(NWB))
  ALLOCATE (KBMAX(NWB),  ELKT(NWB),   WIND2(IMX))
  ALLOCATE (VISC(NWB),   CELC(NWB),   REAERC(NWB))
  ALLOCATE (QOAVR(NWB),  QIMXR(NWB),  QOMXR(NWB))
  ALLOCATE (LAT(NWB),    LONGIT(NWB), ELBOT(NWB))
  ALLOCATE (BTH(NWB),    VPR(NWB),    LPR(NWB))
  ALLOCATE (NISNP(NWB),  NIPRF(NWB),  NISPR(NWB))
  ALLOCATE (icpl(NWB))                                 ! cb 1/26/09
  ALLOCATE (A00(NWB),    HH(NWB),     DECL(NWB))
  ALLOCATE (T2I(NWB),    KTWB(NWB),   KBR(NWB),    IBPR(NWB))
  ALLOCATE (DLVR(NWB),   ESR(NWB),    ETR(NWB),    NBL(NWB))
  ALLOCATE (LPRFN(NWB),  EXTFN(NWB),  BTHFN(NWB),  METFN(NWB),  VPRFN(NWB))
  ALLOCATE (SNPFN(NWB),  PRFFN(NWB),  SPRFN(NWB),  CPLFN(NWB),  VPLFN(NWB),  FLXFN(NWB),FLXFN2(NWB))
  ALLOCATE (AFW(NWB),    BFW(NWB),    CFW(NWB),    WINDH(NWB),  RHEVC(NWB),  FETCHC(NWB))
  ALLOCATE (SDK(NWB),    FSOD(NWB),   FSED(NWB),   SEDCI(NWB),  SEDCc(NWB),   SEDPRC(NWB), seds(nwb), sedb(nwb), DYNSEDK(NWB))  !cb 11/28/06
  ALLOCATE (ICEC(NWB),   SLICEC(NWB), ICETHI(NWB), ALBEDO(NWB), HWI(NWB),    BETAI(NWB),  GAMMAI(NWB), ICEMIN(NWB), ICET2(NWB))
  ALLOCATE (EXH2O(NWB),  BETA(NWB),   EXOM(NWB),   EXSS(NWB),   DXI(NWB),    CBHE(NWB),   TSED(NWB),   TSEDF(NWB),  FI(NWB))
  ALLOCATE (AX(NWB),     WTYPEC(NWB), JBDN(NWB),   AZC(NWB),    AZMAX(NWB),  QINT(NWB),   QOUTT(NWB),  GRIDC(NWB))     !SW 07/14/04
  ALLOCATE (TAIR(NWB),   TDEW(NWB),   WIND(NWB),   PHI(NWB),    CLOUD(NWB),  CSHE(IMX),   SRON(NWB),   RANLW(NWB))
  ALLOCATE (SNPC(NWB),   SCRC(NWB),   PRFC(NWB),   SPRC(NWB),   CPLC(NWB),   VPLC(NWB),   FLXC(NWB))
  ALLOCATE (NXTMSN(NWB), NXTMSC(NWB), NXTMPR(NWB), NXTMSP(NWB), NXTMCP(NWB), NXTMVP(NWB), NXTMFL(NWB))
  ALLOCATE (SNPDP(NWB),  SCRDP(NWB),  PRFDP(NWB),  SPRDP(NWB),  CPLDP(NWB),  VPLDP(NWB),  FLXDP(NWB))
  ALLOCATE (NSNP(NWB),   NSCR(NWB),   NPRF(NWB),   NSPR(NWB),   NCPL(NWB),   NVPL(NWB),   NFLX(NWB))
  ALLOCATE (NEQN(NWB),   PO4R(NWB),   PARTP(NWB))
  ALLOCATE (NH4DK(NWB),  NH4R(NWB))
  ALLOCATE (NO3DK(NWB),  NO3S(NWB), FNO3SED(NWB))
  ALLOCATE (FER(NWB),    FES(NWB))
  ALLOCATE (CO2R(NWB),   SROC(NWB))
  ALLOCATE (O2ER(NEPT),  O2EG(NEPT))
  ALLOCATE (CAQ10(NWB),  CADK(NWB),   CAS(NWB))
  ALLOCATE (BODP(NBOD),  BODN(NBOD),  BODC(NBOD))
  ALLOCATE (KBOD(NBOD),  TBOD(NBOD),  RBOD(NBOD))
  ALLOCATE (LDOMDK(NWB), RDOMDK(NWB), LRDDK(NWB))
  ALLOCATE (OMT1(NWB),   OMT2(NWB),   OMK1(NWB),   OMK2(NWB))
  ALLOCATE (LPOMDK(NWB), RPOMDK(NWB), LRPDK(NWB),  POMS(NWB))
  ALLOCATE (ORGP(NWB),   ORGN(NWB),   ORGC(NWB),   ORGSI(NWB))
  ALLOCATE (RCOEF1(NWB), RCOEF2(NWB), RCOEF3(NWB), RCOEF4(NWB))
  ALLOCATE (NH4T1(NWB),  NH4T2(NWB),  NH4K1(NWB),  NH4K2(NWB))
  ALLOCATE (NO3T1(NWB),  NO3T2(NWB),  NO3K1(NWB),  NO3K2(NWB))
  ALLOCATE (DSIR(NWB),   PSIS(NWB),   PSIDK(NWB),  PARTSI(NWB))
  ALLOCATE (SODT1(NWB),  SODT2(NWB),  SODK1(NWB),  SODK2(NWB))
  ALLOCATE (O2NH4(NWB),  O2OM(NWB))
  ALLOCATE (O2AR(NAL),   O2AG(NAL))
  ALLOCATE (CGQ10(NGC),  CG0DK(NGC),  CG1DK(NGC),  CGS(NGC))
  ALLOCATE (CUNIT(NCT),  CUNIT1(NCT), CUNIT2(NCT))
  ALLOCATE (CAC(NCT),    INCAC(NCT),  TRCAC(NCT),  DTCAC(NCT),  PRCAC(NCT))
  ALLOCATE (CNAME(NCT),  CNAME1(NCT), CNAME2(NCT), CNAME3(NCT), CMULT(NCT),  CSUM(NCT))
  ALLOCATE (CN(NCT))
  ALLOCATE (SSS(NSS),    TAUCR(NSS),  SEDRC(NSS))
  ALLOCATE (CDSUM(NDC))
  ALLOCATE (DTRC(NBR))
  ALLOCATE (NSTR(NBR),   XBR(NBR))
  ALLOCATE (QTAVB(NBR),  QTMXB(NBR))
  ALLOCATE (BS(NWB),     BE(NWB),     JBUH(NBR),   JBDH(NBR),   JWUH(NBR),   JWDH(NBR))
  ALLOCATE (TSSS(NBR),   TSSB(NBR),   TSSICE(NBR))
  ALLOCATE (ESBR(NBR),   ETBR(NBR),   EBRI(NBR))
  ALLOCATE (QIN(NBR),    PR(NBR),     QPRBR(NBR),  QDTR(NBR),   EVBR(NBR))
  ALLOCATE (TIN(NBR),    TOUT(NBR),   TPR(NBR),    TDTR(NBR),   TPB(NBR))
  ALLOCATE (NACPR(NBR),  NACIN(NBR),  NACDT(NBR),  NACTR(NTR),  NACD(NWB))
  ALLOCATE (QSUM(NBR),   NOUT(NBR),   KTQIN(NBR),  KBQIN(NBR),  ELUH(NBR),   ELDH(NBR))
  ALLOCATE (NL(NBR),     NPOINT(NBR), SLOPE(NBR),  ALPHA(NBR),  COSA(NBR),   SINA(NBR), ilayer(imx))    ! SW 1/23/06
  ALLOCATE (CPRFN(NBR),  EUHFN(NBR),  TUHFN(NBR),  CUHFN(NBR),  EDHFN(NBR),  TDHFN(NBR),  QOTFN(NBR),  PREFN(NBR))
  ALLOCATE (QINFN(NBR),  TINFN(NBR),  CINFN(NBR),  CDHFN(NBR),  QDTFN(NBR),  TDTFN(NBR),  CDTFN(NBR),  TPRFN(NBR))
  ALLOCATE (VOLWD(NBR),  VOLSBR(NBR), VOLTBR(NBR), DLVOL(NBR),  VOLG(NWB),   VOLSR(NWB),  VOLTR(NWB),  VOLEV(NBR))
  ALLOCATE (VOLB(NBR),   VOLPR(NBR),  VOLTRB(NBR), VOLDT(NBR),  VOLUH(NBR),  VOLDH(NBR),  VOLIN(NBR),  VOLOUT(NBR))
  ALLOCATE (US(NBR),     DS(NBR),     CUS(NBR),    UHS(NBR),    DHS(NBR),    UQB(NBR),    DQB(NBR),    CDHS(NBR))
  ALLOCATE (TSSEV(NBR),  TSSPR(NBR),  TSSTR(NBR),  TSSDT(NBR),  TSSWD(NBR),  TSSUH(NBR),  TSSDH(NBR),  TSSIN(NBR),  TSSOUT(NBR))
  ALLOCATE (ET(IMX),     RS(IMX),     RN(IMX),     RB(IMX),     RC(IMX),     RE(IMX),     SHADE(IMX))
  ALLOCATE (DLTMAX(NOD), QWDO(IMX),   TWDO(IMX))                                                                        ! SW 1/24/05
  ALLOCATE (SOD(IMX),    ELWS(IMX),   BKT(IMX),    REAER(IMX))
  ALLOCATE (ICETH(IMX),  ICE(IMX),    ICESW(IMX))
  ALLOCATE (Q(IMX),      QC(IMX),     QERR(IMX),   QSSUM(IMX))
  ALLOCATE (KTI(IMX),    SKTI(IMX),   SROSH(IMX),  SEG(IMX),    DLXRHO(IMX))
  ALLOCATE (DLX(IMX),    DLXR(IMX))
  ALLOCATE (A(IMX),      C(IMX),      D(IMX),      F(IMX),      V(IMX),      BTA(IMX),    GMA(IMX))
  ALLOCATE (KBMIN(IMX),  EV(IMX),     QDT(IMX),    QPR(IMX),    SBKT(IMX),   BHRHO(IMX))
  ALLOCATE (SZ(IMX),     WSHX(IMX),   WSHY(IMX),   WIND10(IMX), CZ(IMX),     FETCH(IMX),  PHI0(IMX),   FRIC(IMX))
  ALLOCATE (Z(IMX),      KB(IMX),     PALT(IMX))
  ALLOCATE (VNORM(KMX))
  ALLOCATE (ANPR(NAL),   ANEQN(NAL),  APOM(NAL))
  ALLOCATE (AC(NAL),     ASI(NAL),    ACHLA(NAL),  AHSP(NAL),   AHSN(NAL),   AHSSI(NAL))
  ALLOCATE (AT1(NAL),    AT2(NAL),    AT3(NAL),    AT4(NAL),    AK1(NAL),    AK2(NAL),    AK3(NAL),    AK4(NAL))
  ALLOCATE (AG(NAL),     AR(NAL),     AE(NAL),     AM(NAL),     AS(NAL),     EXA(NAL),    ASAT(NAL),   AP(NAL),   AN(NAL))
  ALLOCATE (ENPR(NEPT),  ENEQN(NEPT))
  ALLOCATE (EG(NEPT),    ER(NEPT),    EE(NEPT),    EM(NEPT),    EB(NEPT),    ESAT(NEPT),  EP(NEPT),    EN(NEPT))
  ALLOCATE (EC(NEPT),    ESI(NEPT),   ECHLA(NEPT), EHSP(NEPT),  EHSN(NEPT),  EHSSI(NEPT), EPOM(NEPT),  EHS(NEPT))
  ALLOCATE (ET1(NEPT),   ET2(NEPT),   ET3(NEPT),   ET4(NEPT),   EK1(NEPT),   EK2(NEPT),   EK3(NEPT),   EK4(NEPT))
  ALLOCATE (HNAME(NHY),  FMTH(NHY),   HMULT(NHY),  FMTC(NCT),   FMTCD(NDC))
  ALLOCATE (KFAC(NFL),   KFNAME(NFL), KFNAME2(NFL),KFCN(NFL,NWB))
  ALLOCATE (C2I(NCT,NWB),    TRCN(NCT,NTR))
  ALLOCATE (CDN(NDC,NWB),    CDNAME(NDC),     CDNAME2(NDC),    CDNAME3(NDC),    CDMULT(NDC))
  ALLOCATE (CMBRS(NCT,NBR),  CMBRT(NCT,NBR),  INCN(NCT,NBR),   DTCN(NCT,NBR),   PRCN(NCT,NBR))
  ALLOCATE (FETCHU(IMX,NBR), FETCHD(IMX,NBR))
  ALLOCATE (IPRF(IMX,NWB),   ISNP(IMX,NWB),   ISPR(IMX,NWB),   BL(IMX,NWB))
  ALLOCATE (H1(KMX,IMX),     H2(KMX,IMX),     BH1(KMX,IMX),    BH2(KMX,IMX),    BHR1(KMX,IMX),   BHR2(KMX,IMX),   QTOT(KMX,IMX))
  ALLOCATE (SAVH2(KMX,IMX),  AVH1(KMX,IMX),   AVH2(KMX,IMX),   AVHR(KMX,IMX),   SAVHR(KMX,IMX))
  ALLOCATE (LFPR(KMX,IMX),   BI(KMX,IMX), bnew(kmx,imx))        ! SW 1/23/06
  ALLOCATE (ADX(KMX,IMX),    ADZ(KMX,IMX),    DO1(KMX,IMX),    DO2(KMX,IMX),    DO3(KMX,IMX),    SED(KMX,IMX))
  ALLOCATE (B(KMX,IMX),      CONV(KMX,IMX),   CONV1(KMX,IMX),  EL(KMX,IMX),     DZ(KMX,IMX),     DZQ(KMX,IMX),    DX(KMX,IMX))
  ALLOCATE (P(KMX,IMX),      SU(KMX,IMX),     SW(KMX,IMX),     SAZ(KMX,IMX),    T1(KMX,IMX),     TSS(KMX,IMX),    QSS(KMX,IMX))
  ALLOCATE (BB(KMX,IMX),     BR(KMX,IMX),     BH(KMX,IMX),     BHR(KMX,IMX),    VOL(KMX,IMX),    HSEG(KMX,IMX),   DECAY(KMX,IMX))
  ALLOCATE (DEPTHB(KMX,IMX), DEPTHM(KMX,IMX), FPSS(KMX,IMX),   FPFE(KMX,IMX),   FRICBR(KMX,IMX), UXBR(KMX,IMX),   UYBR(KMX,IMX))
  ALLOCATE (QUH1(KMX,NBR),   QDH1(KMX,NBR),   VOLUH2(KMX,NBR), VOLDH2(KMX,NBR), TUH(KMX,NBR),    TDH(KMX,NBR))
  ALLOCATE (TSSUH1(KMX,NBR), TSSUH2(KMX,NBR), TSSDH1(KMX,NBR), TSSDH2(KMX,NBR))
  ALLOCATE (TVP(KMX,NWB),    SEDVP(KMX,NWB),  H(KMX,NWB))
  ALLOCATE (QINF(KMX,NBR),   QOUT(KMX,NBR),   KOUT(KMX,NBR))
  ALLOCATE (CT(KMX,IMX),     AT(KMX,IMX),     VT(KMX,IMX),     DT(KMX,IMX),     GAMMA(KMX,IMX))
  ALLOCATE (CWDO(NCT,NOD),   CDWDO(NDC,NOD),  CWDOC(NCT),      CDWDOC(NDC),     CDTOT(NDC))
  ALLOCATE (CIN(NCT,NBR),    CDTR(NCT,NBR),   CPR(NCT,NBR),    CPB(NCT,NBR),    COUT(NCT,NBR))
  ALLOCATE (RSOD(NOD),       RSOF(NOD),       DLTD(NOD),       DLTF(NOD))
  ALLOCATE (TSRD(NOD),       TSRF(NOD),       WDOD(NOD),       WDOF(NOD))
  ALLOCATE (SNPD(NOD,NWB),   SNPF(NOD,NWB),   SPRD(NOD,NWB),   SPRF(NOD,NWB))
  ALLOCATE (SCRD(NOD,NWB),   SCRF(NOD,NWB),   PRFD(NOD,NWB),   PRFF(NOD,NWB))
  ALLOCATE (CPLD(NOD,NWB),   CPLF(NOD,NWB),   VPLD(NOD,NWB),   VPLF(NOD,NWB),   FLXD(NOD,NWB),   FLXF(NOD,NWB))
  ALLOCATE (EPIC(NWB,NEPT),  EPICI(NWB,NEPT), EPIPRC(NWB,NEPT))
  ALLOCATE (EPIVP(KMX,NWB,NEP))
  ALLOCATE (CUH(KMX,NCT,NBR),     CDH(KMX,NCT,NBR))
  ALLOCATE (EPM(KMX,IMX,NEPT),    EPD(KMX,IMX,NEPT),    EPC(KMX,IMX,NEPT))
  ALLOCATE (C1S(KMX,IMX,NCT),     CSSB(KMX,IMX,NCT),    CVP(KMX,NCT,NWB))
  ALLOCATE (CSSUH1(KMX,NCT,NBR),  CSSUH2(KMX,NCT,NBR),  CSSDH2(KMX,NCT,NBR), CSSDH1(KMX,NCT,NBR))
  ALLOCATE (OPEN_VPR(NWB),        OPEN_LPR(NWB))
  ALLOCATE (READ_EXTINCTION(NWB), READ_RADIATION(NWB))
  ALLOCATE (DIST_TRIBS(NBR),      LIMITING_FACTOR(NAL))
  ALLOCATE (UPWIND(NWB),          ULTIMATE(NWB))
  ALLOCATE (STRIC(NST,NBR), ESTRT(NST,NBR), WSTRT(NST,NBR), KTSWT(NST,NBR), KBSWT(NST,NBR),SINKCT(NST,NBR))
  ALLOCATE (FRESH_WATER(NWB),     SALT_WATER(NWB),      TRAPEZOIDAL(NWB))                                              !SW 07/16/04
  ALLOCATE (UH_EXTERNAL(NBR),     DH_EXTERNAL(NBR),     UH_INTERNAL(NBR),    DH_INTERNAL(NBR))
  ALLOCATE (UQ_EXTERNAL(NBR),     DQ_EXTERNAL(NBR),     UQ_INTERNAL(NBR),    DQ_INTERNAL(NBR))
  ALLOCATE (UP_FLOW(NBR),         DN_FLOW(NBR),         UP_HEAD(NBR),        DN_HEAD(NBR))
  ALLOCATE (INTERNAL_FLOW(NBR),   DAM_INFLOW(NBR),      DAM_OUTFLOW(NBR),    HEAD_FLOW(NBR),      HEAD_BOUNDARY(NWB))  !TC 08/03/04
  ALLOCATE (ISO_CONC(NCT,NWB),    VERT_CONC(NCT,NWB),   LONG_CONC(NCT,NWB))
  ALLOCATE (ISO_SEDIMENT(NWB),    VERT_SEDIMENT(NWB),   LONG_SEDIMENT(NWB))
  ALLOCATE (VISCOSITY_LIMIT(NWB), CELERITY_LIMIT(NWB),  IMPLICIT_AZ(NWB))
  ALLOCATE (FETCH_CALC(NWB),      ONE_LAYER(IMX),       IMPLICIT_VISC(NWB))
  ALLOCATE (LIMITING_DLT(NWB),    TERM_BY_TERM(NWB),    MANNINGS_N(NWB))
  ALLOCATE (PLACE_QIN(NWB),       PLACE_QTR(NTRT),      SPECIFY_QTR(NTRT))
  ALLOCATE (PRINT_CONST(NCT,NWB), PRINT_HYDRO(NHY,NWB), PRINT_SEDIMENT(NWB))
  ALLOCATE (VOLUME_BALANCE(NWB),  ENERGY_BALANCE(NWB),  MASS_BALANCE(NWB))
  ALLOCATE (DETAILED_ICE(NWB),    ICE_CALC(NWB),        ICE_IN(NBR),          ALLOW_ICE(IMX))
  ALLOCATE (EVAPORATION(NWB),     PRECIPITATION(NWB),   RH_EVAP(NWB),         PH_CALC(NWB))
  ALLOCATE (NO_INFLOW(NWB),       NO_OUTFLOW(NWB),      NO_HEAT(NWB),         NO_WIND(NWB))
  ALLOCATE (ISO_TEMP(NWB),        VERT_TEMP(NWB),       LONG_TEMP(NWB),       VERT_PROFILE(NWB),  LONG_PROFILE(NWB))
  ALLOCATE (SNAPSHOT(NWB),        PROFILE(NWB),         VECTOR(NWB),          CONTOUR(NWB),       SPREADSHEET(NWB))
  ALLOCATE (SCREEN_OUTPUT(NWB),   FLUX(NWB))
  ALLOCATE (PRINT_DERIVED(NDC,NWB),  PRINT_EPIPHYTON(NWB,NEPT))
  ALLOCATE (SEDIMENT_CALC(NWB),      EPIPHYTON_CALC(NWB,NEPT), SEDIMENT_RESUSPENSION(NSS),BOD_CALC(NBOD),ALG_CALC(NAL))
  ALLOCATE (TDG_SPILLWAY(NWDT,NSP),  TDG_GATE(NWDT,NGT),       INTERNAL_WEIR(KMX,IMX))
  ALLOCATE (ISO_EPIPHYTON(NWB,NEPT), VERT_EPIPHYTON(NWB,NEPT), LONG_EPIPHYTON(NWB,NEPT))
  ALLOCATE (LATERAL_SPILLWAY(NSP),   LATERAL_GATE(NGT),        LATERAL_PUMP(NPU),        LATERAL_PIPE(NPI))
  ALLOCATE (INTERP_HEAD(NBR),        INTERP_WITHDRAWAL(NWD),   INTERP_EXTINCTION(NWB),   INTERP_DTRIBS(NBR))
  ALLOCATE (INTERP_OUTFLOW(NST,NBR), INTERP_INFLOW(NBR),       INTERP_METEOROLOGY(NWB),  INTERP_TRIBS(NTR))
  ALLOCATE (LNAME(NCT+NHY+NDC))
  ALLOCATE (IWR(NIW),    KTWR(NIW),   KBWR(NIW))
  ALLOCATE (JWUSP(NSP),  JWDSP(NSP),  QSP(NSP))
  ALLOCATE (KTWD(NWDT),  KBWD(NWDT),  JBWD(NWDT))
  ALLOCATE (GTA1(NGT),   GTB1(NGT),   GTA2(NGT),   GTB2(NGT))
  ALLOCATE (BGT(NGT),    IUGT(NGT),   IDGT(NGT),   EGT(NGT))
  ALLOCATE (QTR(NTRT),   TTR(NTRT),   KTTR(NTRT),  KBTR(NTRT))
  ALLOCATE (AGASGT(NGT), BGASGT(NGT), CGASGT(NGT), GASGTC(NGT))
  ALLOCATE (PUGTC(NGT),  ETUGT(NGT),  EBUGT(NGT),  KTUGT(NGT),  KBUGT(NGT))
  ALLOCATE (PDGTC(NGT),  ETDGT(NGT),  EBDGT(NGT),  KTDGT(NGT),  KBDGT(NGT))
  ALLOCATE (A1GT(NGT),   B1GT(NGT),   G1GT(NGT),   A2GT(NGT),   B2GT(NGT),   G2GT(NGT))
  ALLOCATE (EQGT(NGT),   JBUGT(NGT),  JBDGT(NGT),  JWUGT(NGT),  JWDGT(NGT),  QGT(NGT))
  ALLOCATE (JBUPI(NPI),  JBDPI(NPI),  JWUPI(NPI),  JWDPI(NPI),  QPI(NPI), BP(NPI))                              ! SW 5/10/10
  ALLOCATE (IUPI(NPI),   IDPI(NPI),   EUPI(NPI),   EDPI(NPI),   WPI(NPI),    DLXPI(NPI),  FPI(NPI),    FMINPI(NPI), PUPIC(NPI))
  ALLOCATE (ETUPI(NPI),  EBUPI(NPI),  KTUPI(NPI),  KBUPI(NPI),  PDPIC(NPI),  ETDPI(NPI),  EBDPI(NPI),  KTDPI(NPI),  KBDPI(NPI))
  ALLOCATE (PUSPC(NSP),  ETUSP(NSP),  EBUSP(NSP),  KTUSP(NSP),  KBUSP(NSP),  PDSPC(NSP),  ETDSP(NSP),  EBDSP(NSP))
  ALLOCATE (KTDSP(NSP),  KBDSP(NSP),  IUSP(NSP),   IDSP(NSP),   ESP(NSP),    A1SP(NSP),   B1SP(NSP),   A2SP(NSP))
  ALLOCATE (B2SP(NSP),   AGASSP(NSP), BGASSP(NSP), CGASSP(NSP), EQSP(NSP),   GASSPC(NSP), JBUSP(NSP),  JBDSP(NSP))
  ALLOCATE (IUPU(NPU),   IDPU(NPU),   EPU(NPU),    STRTPU(NPU), ENDPU(NPU),  EONPU(NPU),  EOFFPU(NPU), QPU(NPU),   PPUC(NPU))
  ALLOCATE (ETPU(NPU),   EBPU(NPU),   KTPU(NPU),   KBPU(NPU),   JWUPU(NPU),  JWDPU(NPU),  JBUPU(NPU),  JBDPU(NPU), PUMPON(NPU))
  ALLOCATE (IWD(NWDT),   KWD(NWDT),   QWD(NWDT),   EWD(NWDT),   KTW(NWDT),   KBW(NWDT))
  ALLOCATE (ITR(NTRT),   QTRFN(NTR),  TTRFN(NTR),  CTRFN(NTR),  ELTRT(NTRT), ELTRB(NTRT), TRC(NTRT),   JBTR(NTRT), QTRF(KMX,NTRT))
  ALLOCATE (TTLB(IMX),   TTRB(IMX),   CLLB(IMX),   CLRB(IMX))
  ALLOCATE (SRLB1(IMX),  SRRB1(IMX),  SRLB2(IMX),  SRRB2(IMX),  SRFJD1(IMX), SHADEI(IMX), SRFJD2(IMX))
  ALLOCATE (TOPO(IMX,IANG))                                                                                        ! SW 10/17/05
  ALLOCATE (QSW(KMX,NWDT),  CTR(NCT,NTRT), HPRWBC(NHY,NWB))
  ALLOCATE (RATZ(KMX,NWB),   CURZ1(KMX,NWB),  CURZ2(KMX,NWB),   CURZ3(KMX,NWB))   ! SW 5/15/06
! v3.5 start
  ALLOCATE (zg(NZPt),zm(NZPt),zeff(NZPt),prefp(NZPt),zr(NZPt),zoomin(NZPt),zs2p(NZPt),exz(NZPt),PREFZ(NZPt,nzpt))
  ALLOCATE (zt1(NZPt),zt2(NZPt),zt3(NZPt),zt4(NZPt),zk1(NZPt),zk2(NZPt),zk3(NZPt),zk4(NZPt),o2zr(nzpt))
  ALLOCATE (zp(NZPt),zn(NZPt),zc(NZPt))
  allocate (prefa(nal,nzpt))
  allocate (po4zr(kmx,imx),nh4zr(kmx,imx))
  allocate (zmu(kmx,imx,nzp),tgraze(kmx,imx,nzp),zrt(kmx,imx,nzp),zmt(kmx,imx,nzp)) ! MLM POINTERS:,zoo(kmx,imx,NZP),zooss(kmx,imx,NZP))
  allocate (zoorm(kmx,imx,nzp),zoormr(kmx,imx,nzp),zoormf(kmx,imx,nzp))
  allocate (lpzooout(kmx,imx),lpzooin(kmx,imx),dozr(kmx,imx),ticzr(kmx,imx))
  allocate (agz(kmx,imx,nal,nzp), zgz(kmx,imx,nzp,nzp),agzt(kmx,imx,nal)) !omnivorous zooplankton
  allocate (ORGPLD(kmx,imx), ORGPRD(kmx,imx), ORGPLP(kmx,imx), ORGPRP(kmx,imx), ORGNLD(kmx,imx), ORGNRD(kmx,imx), ORGNLP(kmx,imx))
  allocate (ORGNRP(kmx,imx))
  allocate (ldompmp(kmx,imx),ldomnmp(kmx,imx),lpompmp(kmx,imx),lpomnmp(kmx,imx),rpompmp(kmx,imx),rpomnmp(kmx,imx))
  allocate (lpzooinp(kmx,imx),lpzooinn(kmx,imx),lpzoooutp(kmx,imx),lpzoooutn(kmx,imx))
  allocate (SEDVPp(KMX,NWB),SEDVPc(KMX,NWB),SEDVPn(KMX,NWB))
  allocate (sedp(kmx,imx),sedn(kmx,imx),sedc(kmx,imx))
  allocate (sdkv(kmx,imx),seddktot(kmx,imx))
  allocate (sedcip(nwb),sedcin(nwb),sedcic(nwb),sedcis(nwb))
  ALLOCATE (cbods(NBOD), cbodns(kmx,imx), sedcb(kmx,imx), sedcbp(kmx,imx), sedcbn(kmx,imx), sedcbc(kmx,imx))
  allocate  (print_macrophyte(nwb,nmct), macrophyte_calc(nwb,nmct),macwbc(nwb,nmct),conv2(kmx,kmx),mprwbc(nwb,nmct))
  allocate  (mac(kmx,imx,nmct), macrc(kmx,kmx,imx,nmct),mact(kmx,kmx,imx), macrm(kmx,kmx,imx,nmct), macss(kmx,kmx,imx,nmct))
  allocate  (mgr(kmx,kmx,imx,nmct),mmr(kmx,imx,nmct), mrr(kmx,imx,nmct))
  allocate  (smacrc(kmx,kmx,imx,nmct), smacrm(kmx,kmx,imx,nmct))
  allocate  (smact(kmx,kmx,imx), smac(kmx,imx,nmct))
  allocate  (mt1(nmct),mt2(nmct),mt3(nmct),mt4(nmct),mk1(nmct),mk2(nmct),mk3(nmct),mk4(nmct),mg(nmct),mr(nmct),mm(nmct))
  allocate  (mbmp(nmct), mmax(nmct), cddrag(nmct), dwv(nmct), dwsa(nmct), anorm(nmct))
  allocate  (mp(nmct), mn(nmct), mc(nmct),psed(nmct),nsed(nmct),mhsp(nmct),mhsn(nmct),mhsc(nmct),msat(nmct),exm(nmct))
  allocate  (O2MG(nmct), O2MR(nmct),  LRPMAC(nmct),  MPOM(nmct))
  allocate  (kticol(imx),armac(imx),macwbci(nwb,nmct))
  allocate  (macmbrs(nbr,nmct), macmbrt(nbr,nmct),ssmacmb(nbr,nmct))
  allocate  (cw(kmx,imx), bic(kmx,imx))
  allocate  (mactrmr(kmx,imx,nmct), mactrmf(kmx,imx,nmct),mactrm(kmx,imx,nmct))
  allocate  (mlfpr(kmx,kmx,imx,nmct))
  allocate  (mllim(kmx,kmx,imx,nmct), mplim(kmx,imx,nmct),mclim(kmx,imx,nmct),mnlim(kmx,imx,nmct))
  ALLOCATE  (GAMMAj(kmx,KMX,IMX))	
  allocate (por(kmx,imx),VOLKTi(imx),VOLi(Kmx,Imx),vstem(kmx,imx,nmct),vstemkt(imx,nmct),sarea(nmct))
  ALLOCATE (IWIND(NWB))
  ALLOCATE (LAYERCHANGE(NWB))
! v3.5 end

! Allocate subroutine variables

  CALL TRANSPORT
  CALL KINETICS
  CALL WATERBODY
  CALL OPEN_CHANNEL_INITIALIZE
  CALL PIPE_FLOW_INITIALIZE

! State variables

  TDS  => C2(:,:,1);         PO4  => C2(:,:,NPO4);      NH4  => C2(:,:,NNH4);        NO3  => C2(:,:,NNO3);   DSI  => C2(:,:,NDSI)
  PSI  => C2(:,:,NPSI);      FE   => C2(:,:,NFE);       LDOM => C2(:,:,NLDOM);       RDOM => C2(:,:,NRDOM);  LPOM => C2(:,:,NLPOM)
  RPOM => C2(:,:,NRPOM);     O2   => C2(:,:,NDO);       TIC  => C2(:,:,NTIC);        ALK  => C2(:,:,NALK)
  CG   => C2(:,:,NGCS:NGCE); SS   => C2(:,:,NSSS:NSSE); CBOD => C2(:,:,NBODS:NBODE); ALG  => C2(:,:,NAS:NAE)
! v3.5 start
  ZOO  => C2(:,:,NZOOS:NZOOE)
  LDOMP  => C2(:,:,nldomp); RDOMP  => C2(:,:,nrdomp); LPOMP  => C2(:,:,nlpomp); RPOMP  => C2(:,:,nrpomp)
  LDOMN  => C2(:,:,nldomN); RDOMN  => C2(:,:,nrdomn); LPOMN  => C2(:,:,nlpomn); RPOMN  => C2(:,:,nrpomn)
! v3.5 end

! State variable source/sinks

  CGSS   => CSSK(:,:,NGCS:NGCE);   SSSS   => CSSK(:,:,NSSS:NSSE); PO4SS  => CSSK(:,:,NPO4);  NH4SS  => CSSK(:,:,NNH4)
  NO3SS  => CSSK(:,:,NNO3);        DSISS  => CSSK(:,:,NDSI);      PSISS  => CSSK(:,:,NPSI);  FESS   => CSSK(:,:,NFE)
  LDOMSS => CSSK(:,:,NLDOM);       RDOMSS => CSSK(:,:,NRDOM);     LPOMSS => CSSK(:,:,NLPOM); RPOMSS => CSSK(:,:,NRPOM)
  CBODSS => CSSK(:,:,NBODS:NBODE); ASS    => CSSK(:,:,NAS:NAE);   DOSS   => CSSK(:,:,NDO);   TICSS  => CSSK(:,:,NTIC)
! v3.5 start
  zooss  => cssk(:,:,NZOOS:NZOOE)
  LDOMPSS  => cssk(:,:,nldomp); RDOMPSS  => cssk(:,:,nrdomp); LPOMPSS  => cssk(:,:,nlpomp); RPOMPSS  => cssk(:,:,nrpomp)
  LDOMNSS  => cssk(:,:,nldomN); RDOMNSS  => cssk(:,:,nrdomn); LPOMNSS  => cssk(:,:,nlpomn); RPOMNSS  => cssk(:,:,nrpomn)
! v3.5 end

! Derived variables

  DOC   => CD(:,:,1);  POC  => CD(:,:,2);  TOC  => CD(:,:,3);  DON  => CD(:,:,4);  PON   => CD(:,:,5);  TON  => CD(:,:,6)
  TKN   => CD(:,:,7);  TN   => CD(:,:,8);  DOP  => CD(:,:,9);  POP  => CD(:,:,10); TOP   => CD(:,:,11); TP   => CD(:,:,12)
  APR   => CD(:,:,13); CHLA => CD(:,:,14); ATOT => CD(:,:,15); O2DG => CD(:,:,16); TOTSS => CD(:,:,17); TISS => CD(:,:,18)
  CBODU => CD(:,:,19); PH   => CD(:,:,20); CO2  => CD(:,:,21); HCO3 => CD(:,:,22); CO3   => CD(:,:,23)

! Kinetic fluxes

  SSSI   => KF(:,:,1);  SSSO   => KF(:,:,2);  PO4AR  => KF(:,:,3);  PO4AG  => KF(:,:,4);  PO4AP  => KF(:,:,5)
  PO4ER  => KF(:,:,6);  PO4EG  => KF(:,:,7);  PO4EP  => KF(:,:,8);  PO4POM => KF(:,:,9);  PO4DOM => KF(:,:,10)
  PO4OM  => KF(:,:,11); PO4SD  => KF(:,:,12); PO4SR  => KF(:,:,13); PO4NS  => KF(:,:,14); NH4D   => KF(:,:,15)
  NH4AR  => KF(:,:,16); NH4AG  => KF(:,:,17); NH4AP  => KF(:,:,18); NH4ER  => KF(:,:,19); NH4EG  => KF(:,:,20)
  NH4EP  => KF(:,:,21); NH4POM => KF(:,:,22); NH4DOM => KF(:,:,23); NH4OM  => KF(:,:,24); NH4SD  => KF(:,:,25)
  NH4SR  => KF(:,:,26); NO3D   => KF(:,:,27); NO3AG  => KF(:,:,28); NO3EG  => KF(:,:,29); NO3SED => KF(:,:,30)
  DSIAG  => KF(:,:,31); DSIEG  => KF(:,:,32); DSID   => KF(:,:,33); DSISD  => KF(:,:,34); DSISR  => KF(:,:,35)
  DSIS   => KF(:,:,36); PSIAM  => KF(:,:,37); PSINS  => KF(:,:,38); PSID   => KF(:,:,39); FENS   => KF(:,:,40)
  FESR   => KF(:,:,41); LDOMD  => KF(:,:,42); LRDOMD => KF(:,:,43); RDOMD  => KF(:,:,44); LDOMAP => KF(:,:,45)
  LDOMEP => KF(:,:,46); LPOMD  => KF(:,:,47); LRPOMD => KF(:,:,48); RPOMD  => KF(:,:,49); LPOMAP => KF(:,:,50)
  LPOMEP => KF(:,:,51); LPOMNS => KF(:,:,52); RPOMNS => KF(:,:,53); CBODDK => KF(:,:,54); DOAP   => KF(:,:,55)
  !DOEP   => KF(:,:,56); DOAR   => KF(:,:,57); DOER   => KF(:,:,58); DOPOM  => KF(:,:,59); DODOM  => KF(:,:,60)
  DOAR   => KF(:,:,56); DOEP   => KF(:,:,57); DOER   => KF(:,:,58); DOPOM  => KF(:,:,59); DODOM  => KF(:,:,60)   ! cb 6/2/2009
  DOOM   => KF(:,:,61); DONIT  => KF(:,:,62); DOBOD  => KF(:,:,63); DOAE   => KF(:,:,64); DOSED  => KF(:,:,65)
  DOSOD  => KF(:,:,66); TICAP  => KF(:,:,67); TICEP  => KF(:,:,68); SEDD   => KF(:,:,69); SEDAS  => KF(:,:,70)
  SEDOMS => KF(:,:,71); SEDNS  => KF(:,:,72); SODD   => KF(:,:,73)
! v3.5 start
  LDOMPAP => KF(:,:,74); LDOMPeP => KF(:,:,75); LPOMpAP => KF(:,:,76); LPOMPNS => KF(:,:,77); RPOMPNS => KF(:,:,78)
  LDOMnAP => KF(:,:,79); LDOMneP => KF(:,:,80); LPOMnAP => KF(:,:,81); LPOMnNS => KF(:,:,82); RPOMnNS => KF(:,:,83)
  SEDDp   => KF(:,:,84); SEDASp  => KF(:,:,85); SEDOMSp => KF(:,:,86); SEDNSp  => KF(:,:,87); lpomepp => KF(:,:,88)
  SEDDn   => KF(:,:,89); SEDASn  => KF(:,:,90); SEDOMSn => KF(:,:,91); SEDNSn  => KF(:,:,92); lpomepn => KF(:,:,93)
  SEDDc   => KF(:,:,94); SEDASc  => KF(:,:,95); SEDOMSc => KF(:,:,96); SEDNSc  => KF(:,:,97); lpomepc => KF(:,:,98)
  sedno3  => KF(:,:,99)
  po4mr   => KF(:,:,100);po4mg   => KF(:,:,101); nh4mr   => KF(:,:,102); nh4mg => KF(:,:,103); ldommac => KF(:,:,104)
  rpommac => KF(:,:,105);lpommac => KF(:,:,106); domp    => KF(:,:,107); domr  => KF(:,:,108); ticmc   => KF(:,:,109)
  cbodns  => KF(:,:,110);sedcb   => KF(:,:,111); sedcbp  => KF(:,:,112); sedcbn => KF(:,:,113); sedcbc  => KF(:,:,114)
  sedbr   => KF(:,:,115);sedbrp  => KF(:,:,116); sedbrn  => KF(:,:,117); sedbrc  => KF(:,:,118)
! v3.5 end

! Algal rate variables

  AGR => AF(:,:,:,1); ARR => AF(:,:,:,2); AER => AF(:,:,:,3); AMR => AF(:,:,:,4); ASR => AF(:,:,:,5)
  EGR => EF(:,:,:,1); ERR => EF(:,:,:,2); EER => EF(:,:,:,3); EMR => EF(:,:,:,4); EBR => EF(:,:,:,5)

! Hydrodynamic variables

  DLTLIM => HYD(:,:,1);  U   => HYD(:,:,2);  W    => HYD(:,:,3); T2   => HYD(:,:,4);  RHO => HYD(:,:,5);  AZ  => HYD(:,:,6)
  VSH    => HYD(:,:,7);  ST  => HYD(:,:,8);  SB   => HYD(:,:,9); ADMX => HYD(:,:,10); DM  => HYD(:,:,11); HDG => HYD(:,:,12)
  ADMZ   => HYD(:,:,13); HPG => HYD(:,:,14); GRAV => HYD(:,:,15)

! I/O units

  SNP => OPT(:,1); PRF => OPT(:,2); VPL => OPT(:,3); CPL => OPT(:,4); SPR => OPT(:,5); FLX => OPT(:,6);FLX2 => OPT(:,7)

! Zero variables

  ITR  = 0;   JBTR = 0;   KTTR = 0;   KBTR = 0;   QTR  = 0.0; TTR  = 0.0; CTR  = 0.0; QTRF = 0.0; SNPD  = 0.0; TSRD  = 0.0
  PRFD = 0.0; SPRD = 0.0; CPLD = 0.0; VPLD = 0.0; SCRD = 0.0; FLXD = 0.0; WDOD = 0.0; RSOD = 0.0; ELTRB = 0.0; ELTRT = 0.0

! Input file unit numbers

  NUNIT = 40
  DO JW=1,NWB
    BTH(JW) = NUNIT
    VPR(JW) = NUNIT+1
    LPR(JW) = NUNIT+2
    NUNIT   = NUNIT+3
  END DO
  GRF = NUNIT; NUNIT = NUNIT+1

! Time control cards

  READ (CON,'(//8X,2F8.0,I8)')         TMSTRT,   TMEND,    YEAR
  READ (CON,'(//8X,I8,8F8.0)')         NDLT,     DLTMIN
  READ (CON,'(//(:8X,9F8.0))')        (DLTD(J),            J =1,NDLT)
  READ (CON,'(//(:8X,9F8.0))')        (DLTMAX(J),          J =1,NDLT)
  READ (CON,'(//(:8X,9F8.0))')        (DLTF(J),            J =1,NDLT)
  READ (CON,'(//(8X,2A8))')           (VISC(JW), CELC(JW), JW=1,NWB)

! Grid definition cards

  READ (CON,'(//(8X,7I8,F8.3))')      (US(JB),  DS(JB),     UHS(JB),   DHS(JB), UQB(JB), DQB(JB),  NL(JB), SLOPE(JB), JB=1,NBR)
  READ (CON,'(//(8X,3F8.0,3I8))')     (LAT(JW), LONGIT(JW), ELBOT(JW), BS(JW),  BE(JW),  JBDN(JW),                    JW=1,NWB)

! Initial condition cards

  READ (CON,'(//(8X,2F8.0,2A8))')     (T2I(JW),    ICETHI(JW),  WTYPEC(JW),  GRIDC(JW),                               JW=1,NWB)
  READ (CON,'(//(8X,6A8))')           (VBC(JW),    EBC(JW),     MBC(JW),     PQC(JW),   EVC(JW),   PRC(JW),           JW=1,NWB)
  READ (CON,'(//(8X,4A8))')           (WINDC(JW),  QINC(JW),    QOUTC(JW),   HEATC(JW),                               JW=1,NWB)
  READ (CON,'(//(8X,3A8))')           (QINIC(JB),  DTRIC(JB),   HDIC(JB),                                             JB=1,NBR)
  READ (CON,'(//(8X,5A8,4F8.0))')     (SLHTC(JW),  SROC(JW),    RHEVC(JW),   METIC(JW), FETCHC(JW), AFW(JW),                       &
                                       BFW(JW),    CFW(JW),     WINDH(JW),                                            JW=1,NWB)
  READ (CON,'(//(8X,2A8,6F8.0))')     (ICEC(JW),   SLICEC(JW),  ALBEDO(JW),  HWI(JW),   BETAI(JW),  GAMMAI(JW),                    &
                                       ICEMIN(JW), ICET2(JW),                                                         JW=1,NWB)
  READ (CON,'(//(8X,A8,F8.0))')       (SLTRC(JW),  THETA(JW),                                                         JW=1,NWB)
  READ (CON,'(//(8X,6F8.0,A8,F8.0))')      (AX(JW),     DXI(JW),     CBHE(JW),    TSED(JW),  FI(JW),     TSEDF(JW),                     &
                                       FRICC(JW), Z0(JW),                                                                    JW=1,NWB)
  READ (CON,'(//(8X,2A8,F8.0,I8,F8.0,F8.0,F8.0,F8.0,A8))')     (AZC(JW),    AZSLC(JW),   AZMAX(JW),   TKEBC(JW),EROUGH(JW),       &
                                       ARODI(JW),STRICK(JW),TKELATPRDCONST(JW),IMPTKE(JW),JW=1,NWB)          !,PHISET(JW

  DO JW=1,NWB
  IF(Z0(JW) <= 0.0)Z0(JW)=0.001      ! SW 11/28/07
   DO JB=BS(JW),BE(JW)
    DO I=US(JB),DS(JB)
    E(I)=EROUGH(JW)
    ENDDO
   ENDDO
  ENDDO

! Inflow-outflow cards

  READ (CON,'(//(8X,I8))')            (NSTR(JB),      JB=1,NBR)
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9A8)')            (STRIC(JS,JB),  JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9I8)')            (KTSWT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9I8)')            (KBSWT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9A8)')            (SINKCT(JS,JB),JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9F8.0)')          (ESTRT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9F8.0)')          (WSTRT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(//(:8X,2I8,6F8.0,A8,A8))') (IUPI(JP),   IDPI(JP),   EUPI(JP),   EDPI(JP),    WPI(JP),                                   &
                                       DLXPI(JP),  FPI(JP),    FMINPI(JP), LATPIC(JP),  DYNPIPE(JP),JP=1,NPI)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PUPIC(JP),  ETUPI(JP),  EBUPI(JP),  KTUPI(JP),   KBUPI(JP),  JP=1,NPI)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PDPIC(JP),  ETDPI(JP),  EBDPI(JP),  KTDPI(JP),   KBDPI(JP),  JP=1,NPI)
  READ (CON,'(//(:8X,2I8,5F8.0,A8))') (IUSP(JS),   IDSP(JS),   ESP(JS),    A1SP(JS),    B1SP(JS),                                  &
                                       A2SP(JS),   B2SP(JS),   LATSPC(JS),                          JS=1,NSP)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PUSPC(JS),  ETUSP(JS),  EBUSP(JS),  KTUSP(JS),   KBUSP(JS),  JS=1,NSP)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PDSPC(JS),  ETDSP(JS),  EBDSP(JS),  KTDSP(JS),   KBDSP(JS),  JS=1,NSP)
  READ (CON,'(//(:8X,A8,I8,3F8.0))')  (GASSPC(JS), EQSP(JS),   AGASSP(JS), BGASSP(JS),  CGASSP(JS), JS=1,NSP)
  READ (CON,'(//(:8X,2I8,7F8.0,A8))') (IUGT(JG),   IDGT(JG),   EGT(JG),    A1GT(JG),    B1GT(JG),                                  &
                                       G1GT(JG),   A2GT(JG),   B2GT(JG),   G2GT(JG),    LATGTC(JG), JG=1,NGT)
  READ (CON,'(//(:8X,4F8.0,A8))')     (GTA1(JG),   GTB1(JG),   GTA2(JG),   GTB2(JG),    DYNGTC(JG), JG=1,NGT)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PUGTC(JG),  ETUGT(JG),  EBUGT(JG),  KTUGT(JG),   KBUGT(JG),  JG=1,NGT)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PDGTC(JG),  ETDGT(JG),  EBDGT(JG),  KTDGT(JG),   KBDGT(JG),  JG=1,NGT)
  READ (CON,'(//(:8X,A8,I8,3F8.0))')  (GASGTC(JG), EQGT(JG),   AGASGT(JG), BGASGT(JG),  CGASGT(JG), JG=1,NGT)
  READ (CON,'(//(:8X,2I8,6F8.0,A8))') (IUPU(JP),   IDPU(JP),   EPU(JP),    STRTPU(JP),  ENDPU(JP),                                 &
                                       EONPU(JP),  EOFFPU(JP), QPU(JP),    LATPUC(JP),              JP=1,NPU)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PPUC(JP),   ETPU(JP),   EBPU(JP),   KTPU(JP),    KBPU(JP),   JP=1,NPU)
  READ (CON,'(//(:8X,9I8))')          (IWR(JW),    JW=1,NIW)
  READ (CON,'(//(:8X,9I8))')          (KTWR(JW),   JW=1,NIW)
  READ (CON,'(//(:8X,9I8))')          (KBWR(JW),   JW=1,NIW)
  READ (CON,'(//(:8X,9A8))')          (WDIC(JW),   JW=1,NWD)
  READ (CON,'(//(:8X,9I8))')          (IWD(JW),    JW=1,NWD)
  READ (CON,'(//(:8X,9F8.0))')        (EWD(JW),    JW=1,NWD)
  READ (CON,'(//(:8X,9I8))')          (KTWD(JW),   JW=1,NWD)
  READ (CON,'(//(:8X,9I8))')          (KBWD(JW),   JW=1,NWD)
  READ (CON,'(//(:8X,9A8))')          (TRC(JT),    JT=1,NTR)
  READ (CON,'(//(:8X,9A8))')          (TRIC(JT),   JT=1,NTR)
  READ (CON,'(//(:8X,9I8))')          (ITR(JT),    JT=1,NTR)
  READ (CON,'(//(:8X,9F8.0))')        (ELTRT(JT),  JT=1,NTR)
  READ (CON,'(//(:8X,9F8.0))')        (ELTRB(JT),  JT=1,NTR)
  READ (CON,'(//(8X,A8))')            (DTRC(JB),   JB=1,NBR)

! Output control cards (excluding constituents)

  READ (CON,'(/)')
  DO JH=1,NHY
    READ (CON,'(:8X,9A8)')            (HPRWBC(JH,JW),JW=1,NWB)
  END DO
  READ (CON,'(//(8X,A8,2I8))')        (SNPC(JW), NSNP(JW), NISNP(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SNPD(J,JW),J=1,NSNP(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SNPF(J,JW),J=1,NSNP(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9I8)')            (ISNP(I,JW),I=1,NISNP(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (SCRC(JW), NSCR(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SCRD(J,JW),J=1,NSCR(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SCRF(J,JW),J=1,NSCR(JW))
  END DO
  READ (CON,'(//(8X,A8,2I8))')        (PRFC(JW), NPRF(JW), NIPRF(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (PRFD(J,JW),J=1,NPRF(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (PRFF(J,JW),J=1,NPRF(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9I8)')            (IPRF(J,JW),J=1,NIPRF(JW))
  END DO
  READ (CON,'(//(8X,A8,2I8))')        (SPRC(JW), NSPR(JW), NISPR(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SPRD(J,JW),J=1,NSPR(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SPRF(J,JW),J=1,NSPR(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9I8)')            (ISPR(J,JW), J=1,NISPR(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (VPLC(JW),  NVPL(JW),  JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (VPLD(J,JW), J=1,NVPL(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (VPLF(J,JW), J=1,NVPL(JW))
  END DO
  READ (CON,'(//(8X,A8,I8,A8))')      (CPLC(JW),   NCPL(JW), TECPLOT(JW),JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (CPLD(J,JW), J=1,NCPL(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (CPLF(J,JW), J=1,NCPL(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (FLXC(JW),   NFLX(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (FLXD(J,JW), J=1,NFLX(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (FLXF(J,JW), J=1,NFLX(JW))
  END DO
  READ (CON,'(//8X,A8,2I8)')           TSRC,    NTSR,    NIKTSR; ALLOCATE (ITSR(MAX(1,NIKTSR)), ETSR(MAX(1,NIKTSR)))
  READ (CON,'(//(:8X,9F8.0))')        (TSRD(J), J=1,NTSR)
  READ (CON,'(//(:8X,9F8.0))')        (TSRF(J), J=1,NTSR)
  READ (CON,'(//(:8X,9I8))')          (ITSR(J), J=1,NIKTSR)
  READ (CON,'(//(:8X,9F8.0))')        (ETSR(J), J=1,NIKTSR)
  READ (CON,'(//8X,A8,2I8)')           WDOC,    NWDO,    NIWDO;  ALLOCATE (IWDO(MAX(1,NIWDO)))
  READ (CON,'(//(:8X,9F8.0))')        (WDOD(J), J=1,NWDO)
  READ (CON,'(//(:8X,9F8.0))')        (WDOF(J), J=1,NWDO)
  READ (CON,'(//(8X,9I8))')           (IWDO(J), J=1,NIWDO)
  READ (CON,'(//8X,A8,I8,A8)')         RSOC,    NRSO,    RSIC
  READ (CON,'(//(:8X,9F8.0))')        (RSOD(J), J=1,NRSO)
  READ (CON,'(//(:8X,9F8.0))')        (RSOF(J), J=1,NRSO)

! Constituent control cards

  READ (CON,'(//8X,2A8,I8)')           CCC, LIMC, CUF
  READ (CON,'(//(2A8))')              (CNAME2(JC),  CAC(JC),      JC=1,NCT)
  READ (CON,'(/)')
  DO JD=1,NDC
    READ (CON,'(A8,(:9A8))')           CDNAME2(JD),(CDWBC(JD,JW), JW=1,NWB)
  END DO
  READ (CON,'(/)')
!  DO JF=1,NFL
  do jf=1,73   ! Fix this later
    READ (CON,'(A8,(:9A8))')         KFNAME2(JF),(KFWBC(JF,JW),  JW=1,NWB)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9F8.0)')          (C2I(JC,JW),    JW=1,NWB)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CPRWBC(JC,JW), JW=1,NWB)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CINBRC(JC,JB), JB=1,NBR)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CTRTRC(JC,JT), JT=1,NTR)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CDTBRC(JC,JB), JB=1,NBR)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(:8X,9A8)')            (CPRBRC(JC,JB), JB=1,NBR)
  END DO

! Kinetics coefficients

  READ (CON,'(//(8X,4F8.0,2A8))')     (EXH2O(JW),  EXSS(JW),   EXOM(JW),   BETA(JW),   EXC(JW),   EXIC(JW),    JW=1,NWB)
  READ (CON,'(//(8X,9F8.0))')         (EXA(JA),                                                                JA=1,NAL)
  READ (CON,'(//(8X,9F8.0))')         (EXZ(Jz),                                                                Jz=1,Nzpt)  !v3.5
  READ (CON,'(//(8X,9F8.0))')         (EXM(Jm),                                                                Jm=1,nmct)  !v3.5
  READ (CON,'(//(8X,4F8.0))')         (CGQ10(JG),  CG0DK(JG),  CG1DK(JG),  CGS(JG),                            JG=1,NGC)
  READ (CON,'(//(8X,F8.0,A,F8.0))')   (SSS(JS),    SEDRC(JS),  TAUCR(JS),                                      JS=1,NSS)
  READ (CON,'(//(8X,9F8.0))')         (AG(JA),     AR(JA),     AE(JA),     AM(JA),     AS(JA),                                     &
                                       AHSP(JA),   AHSN(JA),   AHSSI(JA),  ASAT(JA),                           JA=1,NAL)
  READ (CON,'(//(8X,8F8.0))')         (AT1(JA),    AT2(JA),    AT3(JA),    AT4(JA),    AK1(JA),   AK2(JA),                         &
                                       AK3(JA),    AK4(JA),                                                    JA=1,NAL)
  READ (CON,'(//(8X,6F8.0,I8,F8.0))') (AP(JA),     AN(JA),     AC(JA),     ASI(JA),    ACHLA(JA), APOM(JA),                        &
                                       ANEQN(JA),  ANPR(JA),   JA=1,NAL)
  READ (CON,'(//(8X,9A8))')           (EPIC(JW,1),                                                             JW=1,NWB)
  DO JE=2,NEPT
    READ (CON,'(8X,9A8)')             (EPIC(JW,JE),                                                            JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9A8))')           (EPIPRC(JW,1),                                                           JW=1,NWB)
  DO JE=2,NEPT
    READ (CON,'(8X,9A8)')             (EPIPRC(JW,JE),                                                          JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9F8.0))')         (EPICI(JW,1),                                                            JW=1,NWB)
  DO JE=2,NEPT
    READ (CON,'(8X,9F8.0)')           (EPICI(JW,JE),                                                           JW=1,NWB)
  END DO
  READ (CON,'(//(8X,8F8.0))')         (EG(JE),     ER(JE),     EE(JE),     EM(JE),     EB(JE),    EHSP(JE),                        &
                                       EHSN(JE),   EHSSI(JE),                                                  JE=1,NEP)
  READ (CON,'(//(8X,2F8.0,I8,F8.0))') (ESAT(JE),   EHS(JE),    ENEQN(JE),  ENPR(JE),                           JE=1,NEP)
  READ (CON,'(//(8X,8F8.0))')         (ET1(JE),    ET2(JE),    ET3(JE),    ET4(JE),    EK1(JE),   EK2(JE),                         &
                                       EK3(JE),    EK4(JE),                                                    JE=1,NEP)
  READ (CON,'(//(8X,6F8.0))')         (EP(JE),     EN(JE),     EC(JE),     ESI(JE),    ECHLA(JE), EPOM(JE),    JE=1,NEP)
! v3.5 start
  READ (CON,'(//(8X,7F8.0))')         (zg(jz),zr(jz),zm(jz),zeff(jz),PREFP(jz),ZOOMIN(jz),ZS2P(jz),            Jz=1,Nzpt)

  READ (CON,'(//(8X,8F8.0))')         (PREFA(ja,1),                                                            Ja=1,nal)          ! MM 7/13/06
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))')       (PREFA(ja,jz),                                                           Ja=1,nal)
  end do
  READ (CON,'(//(8X,8F8.0))')       (PREFz(jjz,1),                                                          Jjz=1,nzpt)
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))')       (PREFz(jjz,jz),                                                          Jjz=1,nzpt)           ! MM 7/13/06
  end do
  READ (CON,'(//(8X,8F8.0))')         (zT1(Jz),    zT2(Jz),    zT3(Jz),    zT4(Jz),    zK1(Jz),   zK2(Jz),                         &
                                       zK3(Jz),    zK4(Jz),                                                    Jz=1,Nzpt)
  READ (CON,'(//(8X,3F8.0))')         (zP(Jz),     zN(Jz),     zC(Jz),                                         Jz=1,Nzpt)
  READ (CON,'(//(8X,9A8))')           (macwbC(JW,1),                                                           JW=1,NWB)
  DO Jm=2,nmct
    READ (CON,'(8X,9A8)')             (macwbC(JW,Jm),                                                          JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9A8))')           (mprwbC(JW,1),                                                           JW=1,NWB)
  DO Jm=2,nmct
    READ (CON,'(8X,9A8)')             (mprwbC(JW,Jm),                                                          JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9F8.0))')         (macwbCI(JW,1),                                                          JW=1,NWB)
  DO Jm=2,nmct
    READ (CON,'(8X,9F8.0)')           (macwbcI(JW,Jm),                                                         JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9F8.0))')         (mG(jm), mR(jm), mM(jm), msat(jm),mhsp(jm),mhsn(jm),mhsc(jm),                           &
                                          mpom(jm),lrpmac(jm),     jm=1,nmct)
  READ (CON,'(//(8X,2F8.0))')         (psed(jm), nsed(jm),                                                     jm=1,nmct)
  READ (CON,'(//(8X,2F8.0))')         (mbmp(jm), mmax(jm),                                                     jm=1,nmct)
  READ (CON,'(//(8X,4F8.0))')         (cddrag(jm),dwv(jm),dwsa(jm),anorm(jm),                                 jm=1,nmct)  !cb 6/29/06
  READ (CON,'(//(8X,8F8.0))')         (mT1(Jm),    mT2(Jm),    mT3(Jm),    mT4(Jm),    mK1(Jm),   mK2(Jm),                         &
                                       mK3(Jm),    mK4(Jm),                                                    Jm=1,nmct)
  READ (CON,'(//(8X,3F8.0))')         (mP(Jm),     mN(Jm),     mC(Jm),                                         Jm=1,nmct)
  READ (CON,'(//(8X,3F8.0))')         (LDOMDK(JW), RDOMDK(JW), LRDDK(JW),                                      JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (LPOMDK(JW), RPOMDK(JW), LRPDK(JW),  POMS(JW),                           JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (ORGP(JW),   ORGN(JW),   ORGC(JW),   ORGSI(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (OMT1(JW),   OMT2(JW),   OMK1(JW),   OMK2(JW),                           JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (KBOD(JB),   TBOD(JB),   RBOD(JB), cbods(jb),                           JB=1,NBOD)  !v3.5
  READ (CON,'(//(8X,3F8.0))')         (BODP(JB),   BODN(JB),   BODC(JB),                                       JB=1,NBOD)
  READ (CON,'(//(8X,2F8.0))')         (PO4R(JW),   PARTP(JW),                                                  JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (NH4R(JW),   NH4DK(JW),                                                  JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (NH4T1(JW),  NH4T2(JW),  NH4K1(JW),  NH4K2(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,3F8.0))')         (NO3DK(JW),  NO3S(JW),   FNO3SED(JW),                                    JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (NO3T1(JW),  NO3T2(JW),  NO3K1(JW),  NO3K2(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (DSIR(JW),   PSIS(JW),   PSIDK(JW),  PARTSI(JW),                         JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (FER(JW),    FES(JW),                                                    JW=1,NWB)
  READ (CON,'(//(8X,F8.0))')          (CO2R(JW),                                                               JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (O2NH4(JW),  O2OM(JW),                                                   JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (O2AR(JA),   O2AG(JA),                                                   JA=1,NAL)
  READ (CON,'(//(8X,2F8.0))')         (O2ER(JE),   O2EG(JE),                                                   JE=1,NEPT)
  READ (CON,'(//(8X,F8.0))')          (O2zR(Jz),                                                               Jz=1,Nzpt)
  READ (CON,'(//(8X,2F8.0))')         (O2mR(Jm),   O2mG(jm),                                                   Jm=1,nmct)
  READ (CON,'(//(8X,F8.0))')           KDO
  READ (CON,'(//(8X,2A8,6F8.0,A8))')     (SEDCc(JW),   SEDPRC(JW), SEDCI(JW),  SDK(JW), seds(jw),   FSOD(JW),   FSED(JW), sedb(jw),DYNSEDK(JW),   JW=1,NWB)  ! cb 11/28/06
  READ (CON,'(//(8X,4F8.0))')         (SODT1(JW),  SODT2(JW),  SODK1(JW),  SODK2(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,9F8.0))')         (SOD(I),                                                                  I=1,IMX)
  READ (CON,'(//(8X,A8,I8,4F8.2))')   (REAERC(JW), NEQN(JW),   RCOEF1(JW), RCOEF2(JW), RCOEF3(JW), RCOEF4(JW), JW=1,NWB)

! Input filenames

  READ (CON,'(//(8X,A72))')  RSIFN
  READ (CON,'(//(8X,A72))')  QWDFN
  READ (CON,'(//(8X,A72))')  QGTFN
  READ (CON,'(//(8X,A72))')  WSCFN
  READ (CON,'(//(8X,A72))')  SHDFN
  READ (CON,'(//(8X,A72))') (BTHFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (METFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (EXTFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (VPRFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (LPRFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (QINFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TINFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CINFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (QOTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (QTRFN(JT), JT=1,NTR)
  READ (CON,'(//(8X,A72))') (TTRFN(JT), JT=1,NTR)
  READ (CON,'(//(8X,A72))') (CTRFN(JT), JT=1,NTR)
  READ (CON,'(//(8X,A72))') (QDTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TDTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CDTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (PREFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TPRFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CPRFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (EUHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TUHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CUHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (EDHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TDHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CDHFN(JB), JB=1,NBR)

! Output filenames

  READ (CON,'(//(8X,A72))') (SNPFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (PRFFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (VPLFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (CPLFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (SPRFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (FLXFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))')  TSRFN
  READ (CON,'(//(8X,A72))')  WDOFN
  CLOSE (CON)

! Bathymetry file

  DO JW=1,NWB
    OPEN (BTH(JW),FILE=BTHFN(JW),STATUS='OLD')
  !  READ (BTH(JW),*)
	READ  (BTH(JW),'(a1)')char1                                     ! New Bathymetry format option SW 6/22/09
      if(char1=='$')then
      READ  (BTH(JW),*)
      READ  (BTH(JW),*) AID,(DLX(I),  I=US(BS(JW))-1,DS(BE(JW))+1)
      READ  (BTH(JW),*) AID,(ELWS(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      READ  (BTH(JW),*) AID,(PHI0(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      READ  (BTH(JW),*) AID,(FRIC(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      READ  (BTH(JW),*)
      DO K=1,KMX
      READ  (BTH(JW),*) H(K,JW),(B(K,I),I=US(BS(JW))-1,DS(BE(JW))+1)
      END DO
	  DO I=US(BS(JW))-1,DS(BE(JW))+1
      H2(:,I) = H(:,JW)
      END DO
      else	
    READ (BTH(JW),'(//(10F8.0))') (DLX(I),  I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (ELWS(I), I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (PHI0(I), I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (FRIC(I), I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (H(K,JW), K=1,KMX)
    DO I=US(BS(JW))-1,DS(BE(JW))+1
      READ (BTH(JW),'(//(10F8.0))') (B(K,I), K=1,KMX)
      H2(:,I) = H(:,JW)
    END DO
	endif
    CLOSE (BTH(JW))
  END DO
  H1 = H2
  BI = B

! Output file unit numbers

  ALLOCATE (TSR(NIKTSR))
  ALLOCATE (WDO(NIWDO,4))
  DO J=1,7
    DO JW=1,NWB
      OPT(JW,J) = NUNIT; NUNIT = NUNIT+1
    END DO
  END DO
  DO J=1,NIKTSR
    TSR(J) = NUNIT; NUNIT = NUNIT+1
  END DO
  DO JW=1,NIWDO
    WDO(JW,1) = NUNIT; NUNIT = NUNIT+1
    WDO(JW,2) = NUNIT; NUNIT = NUNIT+1
    WDO(JW,3) = NUNIT; NUNIT = NUNIT+1
    WDO(JW,4) = NUNIT; NUNIT = NUNIT+1
  END DO

! Variable names, formats, multipliers, and Compaq Visual FORTRAN array viewer controls

  OPEN (GRF,FILE='graph.npt',STATUS='OLD')
  READ (GRF,'(///(A43,1X,A9,3F8.0,A8))') (HNAME(J),  FMTH(J),  HMULT(J),  HYMIN(J), HYMAX(J), HPLTC(J), J=1,NHY)
  READ (GRF,'(// (A43,1X,A9,3F8.0,A8))') (CNAME(J),  FMTC(J),  CMULT(J),  CMIN(J),  CMAX(J),  CPLTC(J), J=1,NCT)
  READ (GRF,'(// (A43,1X,A9,3F8.0,A8))') (CDNAME(J), FMTCD(J), CDMULT(J), CDMIN(J), CDMAX(J), CDPLTC(J),J=1,NDC)
  CLOSE (GRF)
  DO JC=1,NCT
    L3         = 1
    L1         = SCAN (CNAME(JC),',')+2
    L2         = SCAN (CNAME(JC)(L1:43),'  ')+L1
    CUNIT(JC)  = CNAME(JC)(L1:L2)
    CNAME1(JC) = CNAME(JC)(1:L1-3)
    CNAME3(JC) = CNAME1(JC)
    DO WHILE (L3 < L1-3)
      IF (CNAME(JC)(L3:L3) == ' ') CNAME3(JC)(L3:L3) = '_'
      L3 = L3+1
    END DO
    CUNIT1(JC) = CUNIT(JC)(1:1)
    CUNIT2(JC) = CUNIT(JC)
    IF (CUNIT(JC)(1:2) == 'mg') THEN
      CUNIT1(JC) = 'g'
      CUNIT2(JC) = 'g/m^3'
    END IF
    IF (CUNIT(JC)(1:2) /= 'g/' .AND. CUNIT(JC)(1:2) /= 'mg') CUNIT1(JC) = '  '
  END DO
  DO JC=1,NDC
    L1          = 1
    L2          = MAX(4,SCAN (CDNAME(JC),',')-1)
    CDNAME3(JC) = CDNAME(JC)(1:L2)
    DO WHILE (L1 < L2)
      IF (CDNAME(JC)(L1:L1) == ' ') CDNAME3(JC)(L1:L1) = '_'
      L1 = L1+1
    END DO
  END DO
  FMTH(1:NHY) = ADJUSTL (FMTH(1:NHY))

! Initialize logical control variables

  VERT_PROFILE = .FALSE.
  LONG_PROFILE = .FALSE.
  CONSTITUENTS =  CCC  == '      ON'
  DO JW=1,NWB
    ISO_TEMP(JW)         = T2I(JW)     >=  0
    VERT_TEMP(JW)        = T2I(JW)     == -1
    LONG_TEMP(JW)        = T2I(JW)     <  -1
    if(constituents)then                     ! cb 12/04/08
! v3.5 start
    ISO_SEDIMENT(JW)     = SEDCI(JW)   >=  0   .AND. SEDCc(JW)   == '      ON'
    VERT_SEDIMENT(JW)    = SEDCI(JW)   == -1.0 .AND. SEDCc(JW)   == '      ON'
    LONG_SEDIMENT(JW)    = SEDCI(JW)   <  -1.0 .AND. SEDCc(JW)   == '      ON'
! v3.5 end
    ISO_EPIPHYTON(JW,:)  = EPICI(JW,:) >=  0   .AND. EPIC(JW,:) == '      ON'
    VERT_EPIPHYTON(JW,:) = EPICI(JW,:) == -1.0 .AND. EPIC(JW,:) == '      ON'
    LONG_EPIPHYTON(JW,:) = EPICI(JW,:) <  -1.0 .AND. EPIC(JW,:) == '      ON'
    DO JC=1,NCT
      ISO_CONC(JC,JW)  = C2I(JC,JW) >=  0.0
      VERT_CONC(JC,JW) = C2I(JC,JW) == -1.0 .AND. CAC(JC) == '      ON'
      LONG_CONC(JC,JW) = C2I(JC,JW) <  -1.0 .AND. CAC(JC) == '      ON'
      IF (VERT_CONC(JC,JW)) VERT_PROFILE(JW) = .TRUE.
      IF (LONG_CONC(JC,JW)) LONG_PROFILE(JW) = .TRUE.
    END DO
    IF (VERT_SEDIMENT(JW))         VERT_PROFILE(JW) = .TRUE.
    IF (LONG_SEDIMENT(JW))         LONG_PROFILE(JW) = .TRUE.
    IF (ANY(VERT_EPIPHYTON(JW,:))) VERT_PROFILE(JW) = .TRUE.
    IF (ANY(LONG_EPIPHYTON(JW,:))) LONG_PROFILE(JW) = .TRUE.
    end if                          ! cb 12/04/08
    IF (VERT_TEMP(JW))             VERT_PROFILE(JW) = .TRUE.
    IF (LONG_TEMP(JW))             LONG_PROFILE(JW) = .TRUE.
!    IF (VERT_SEDIMENT(JW))         VERT_PROFILE(JW) = .TRUE.
!    IF (LONG_SEDIMENT(JW))         LONG_PROFILE(JW) = .TRUE.
!    IF (ANY(VERT_EPIPHYTON(JW,:))) VERT_PROFILE(JW) = .TRUE.
!    IF (ANY(LONG_EPIPHYTON(JW,:))) LONG_PROFILE(JW) = .TRUE.
    do m=1,nmc
      macrophyte_calc(jw,m) = constituents.and.macwbc(jw,m).eq.' ON'
      print_macrophyte(jw,m) = macrophyte_calc(jw,m).AND.mprwbc(jw,m).EQ.' ON'
      if(macrophyte_calc(jw,m))macrophyte_on=.true.
    end do
  END DO

! Open error and warning files

  OPEN (W2ERR,FILE='w2.err',STATUS='UNKNOWN'); OPEN (WRN,FILE='w2.wrn',STATUS='UNKNOWN')

return
end subroutine input
