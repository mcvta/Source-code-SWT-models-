!***********************************************************************************************************************************
!**                                                                                                                               **
!**                                                         CE-QUAL-W2-PRE                                                        **
!**                                                          Version 3.6                                                          **
!**                                                                                                                               **
!**                                                     A preprocessor code for                                                   **
!**                                                           CE-QUAL-W2                                                          **
!**                                                                                                                               **
!**                                                   Developed by: Thomas M. Cole                                                **
!**                                                   Water Quality Modeling Group                                                **
!**                                                   U.S. Army Corps of Engineers                                                **
!**                                                   Waterways Experiment Station                                                **
!**                                                   Vicksburg, Mississippi 39180                                                **
!**                                                                                                                               **
!**                                                              and                                                              **
!**                                                                                                                               **
!**                                                          Scott Wells                                                          **
!**                                          Department of Civil and Environmental Engineering                                    **
!**                                                    Portland State University                                                  **
!**                                                         P.O. Box 751                                                          **
!**                                                  Portland, Oregon  97207-0751                                                 **
!**                                                  Phone number: (503) 725-4276                                                 **
!**                                                  FAX   number: (503) 725-5950                                                 **
!**                                                    e-mail scott@cecs.pdx.edu                                                  **
!**                                                                                                                               **
!***********************************************************************************************************************************

! Module definitions

MODULE PREC
  INTEGER, PARAMETER                               :: I2=SELECTED_INT_KIND(3)
  INTEGER, PARAMETER                               :: R8=SELECTED_REAL_KIND(15)
END MODULE PREC
MODULE FILE_ERROR
 CONTAINS
   SUBROUTINE DUPLICATE_FILENAME(ERR,FILE)
     INTEGER                                       :: ERR
     CHARACTER(2)                                  :: L
     CHARACTER(5)                                  :: FMTA
     CHARACTER(72)                                 :: FILE
     WRITE (L,'(I0)') LEN_TRIM(FILE)
     FMTA = '(A)'
     CALL ERRORS
     WRITE (ERR,FMTA) 'Duplicate file name in control file for file '//TRIM(FILE)
   END SUBROUTINE DUPLICATE_FILENAME
END MODULE FILE_ERROR
MODULE CONTROL
  REAL,    TARGET,   ALLOCATABLE, DIMENSION(:,:,:) :: C2,     CD
  INTEGER                                          :: IMX,    KMX,    UNIT
  INTEGER                                          :: NERR=0, NWRN=0
  INTEGER                                          :: NASS,   NAAL
  INTEGER,           ALLOCATABLE, DIMENSION(:)     :: SSAC,   ALAC
  LOGICAL                                          :: DELETE_ERR, DELETE_WRN
END MODULE CONTROL
MODULE GRIDC
  REAL,              ALLOCATABLE, DIMENSION(:)     :: Z
  REAL,              ALLOCATABLE, DIMENSION(:)     :: CMULT,  DMULT, HMULT, CDMULT
  REAL,              ALLOCATABLE, DIMENSION(:)     :: ICETH,  ALPHA
  REAL,              ALLOCATABLE, DIMENSION(:,:)   :: H,      EL
  REAL,              POINTER,     DIMENSION(:,:)   :: T2
  INTEGER                                          :: JW,     NAC,    INI=10
  INTEGER,           ALLOCATABLE, DIMENSION(:)     :: IEPR,   KEPR,   NADC,   NACTR,  NACIN,  NACPR,  NACDT
  INTEGER,           ALLOCATABLE, DIMENSION(:)     :: KTWB,   KB,     CUS,    US,     DS,     BS,     BE
  INTEGER,           ALLOCATABLE, DIMENSION(:)     :: CN
  INTEGER,           ALLOCATABLE, DIMENSION(:,:)   :: IPR,    CDN
  LOGICAL,           ALLOCATABLE, DIMENSION(:)     :: ICE_CALC
  CHARACTER(9),      ALLOCATABLE, DIMENSION(:)     :: HFMT
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)     :: TITLE
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:,:)   :: HPRWBC, CDWBC,  CPRWBC
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:)   :: CONV
  CHARACTER(43)                                    :: NAME
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: LNAME                      ! SW 1/16/04 added section
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)       :: CDNAME2
  CHARACTER(6),      ALLOCATABLE, DIMENSION(:)       :: CUNIT,  CUNIT2
  CHARACTER(9),      ALLOCATABLE, DIMENSION(:)       :: FMTH,   FMTC,   FMTCD
  CHARACTER(19),     ALLOCATABLE, DIMENSION(:)       :: CNAME1
  CHARACTER(43),     ALLOCATABLE, DIMENSION(:)       :: CNAME,  CDNAME, HNAME
  REAL,          ALLOCATABLE, DIMENSION(:)   :: CMIN,   CMAX,   HYMIN,  HYMAX,  CDMIN,  CDMAX
END MODULE GRIDC
MODULE KINETICS
  REAL,              ALLOCATABLE, DIMENSION(:)     :: ACHLA,  ORGP,   ORGN,   ORGC,   ORGSI,  O2OM,   ALGC,   PARTP
  REAL,    POINTER,               DIMENSION(:,:)   :: TN,     TP,     TKN
  REAL,    POINTER,               DIMENSION(:,:)   :: DON,    DOP,    DOC
  REAL,    POINTER,               DIMENSION(:,:)   :: PON,    POP,    POC
  REAL,    POINTER,               DIMENSION(:,:)   :: TON,    TOP,    TOC
  REAL,    POINTER,               DIMENSION(:,:)   :: CHLA,   ATOT
  REAL,    POINTER,               DIMENSION(:,:)   :: O2DG,   CBODU,  TISS,   TOTSS
  REAL,    POINTER,               DIMENSION(:,:)   :: TDS,    COL,    NH4,    NO3,    PO4,    FE,     DSI,    PSI
  REAL,    POINTER,               DIMENSION(:,:)   :: LDOM,   RDOM,   LPOM,   RPOM,   O2,     TIC,    ALK
  REAL,    POINTER,               DIMENSION(:,:,:) :: SS,     ALG,    CG,     CBOD
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)     :: SROC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CBODS,sedcip,sedcin, sedcic,sedcis,volkti
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: ORGPLD,  ORGPRD,   ORGPLP,    ORGPRP,  ORGNLD,  ORGNRD, ORGNLP, ORGNRP,por,voli
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: LDOMPMP, LDOMNMP,  LPOMPMP,   LPOMNMP, RPOMPMP, RPOMNMP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: LPZOOINP,LPZOOINN, LPZOOOUTP, LPZOOOUTN
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SEDCAR,    SEDN, SEDP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SEDVPC,  SEDVPP, SEDVPN
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SDKV,    SEDDKTOT
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMP,  RDOMP,  LPOMP,  RPOMP,  LDOMN,  RDOMN,  LPOMN,  RPOMN
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMPSS,  RDOMPSS, LPOMPSS, RPOMPSS, LDOMNSS, RDOMNSS
  REAL,    POINTER,               DIMENSION(:,:)     :: LPOMNSS,  RPOMNSS
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMPAP,  LDOMPEP, LPOMPAP, LPOMPNS, RPOMPNS
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMNAP,  LDOMNEP, LPOMNAP, LPOMNNS, RPOMNNS
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDP,    SEDASP,  SEDOMSP, SEDNSP,  LPOMEPP
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDN,    SEDASN,  SEDOMSN, SEDNSN,  LPOMEPN, SEDNO3
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDC,    SEDASC,  SEDOMSC, SEDNSC,  LPOMEPC
  REAL,    POINTER,               DIMENSION(:,:)     :: CBODNS,   SEDCB,   SEDCBP,  SEDCBN,  SEDCBC
END MODULE KINETICS

MODULE MACROPHYTEC
  REAL,    POINTER,               DIMENSION(:,:)     :: NH4MR,  NH4MG,  LDOMMAC, RPOMMAC, LPOMMAC, DOMP, DOMR, TICMC
  REAL,    POINTER,               DIMENSION(:,:)     :: PO4MR,  PO4MG
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MG,     MR,     MM, MMAX,   MBMP
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MT1,    MT2,    MT3,    MT4,    MK1,    MK2,    MK3,    MK4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MP,     MN,     MC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: PSED,   NSED,   MHSP,   MHSN,   MHSC,   msat,   exm
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CDSTEM, dmv, dwsa,anorm
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ARMAC , sarea, iwind
  REAL,              ALLOCATABLE, DIMENSION(:)       :: O2MG,   O2MR,   LRPMAC,  MPOM
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: MACMBRS,MACMBRT,SSMACMB
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: CW,     BIC, macwbci,vstemkt
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MACTRMR,MACTRMF,MACTRM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MMR,    MRR
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MAC,    MACT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MPLIM,  MNLIM, MCLIM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: SMAC,   SMACT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: GAMMAJ,vstem
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MGR
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MACRC,  MACRM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MLLIM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MACSS
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: SMACRC, SMACRM
  INTEGER                                            :: NMC, nmct
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: KTICOL
  LOGICAL,           ALLOCATABLE, DIMENSION(:,:)     :: PRINT_MACROPHYTE, MACROPHYTE_CALC
  LOGICAL                                            :: MACROPHYTE_ON
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:,:)     :: mprwbc, macwbc
  CHARACTER(10),      ALLOCATABLE, DIMENSION(:,:)    :: CONV2
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:,:,:) :: MLFPR
  DATA                                                  SAVOLRAT /9000.0/, DEN /6.0E4/
END MODULE MACROPHYTEC
MODULE ZOOPLANKTONC
  LOGICAL                                            :: ZOOPLANKTON_CALC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: zg,zm,zeff,PREFP,zr,ZOOMIN,ZS2P,EXZ
  REAL,              ALLOCATABLE, DIMENSION(:)       :: Zt1,Zt2,Zt3,Zt4,Zk1,Zk2,Zk3,Zk4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ZP,ZN,ZC,o2zr
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: PREFA, PREFZ ! OMNIVOROUS ZOOPLANKTON
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: po4zr,NH4ZR,DOZR,TICZR,LPZOOOUT,LPZOOIN
  REAL,    POINTER,               DIMENSION(:,:,:)   :: ZOO, ZOOSS
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ZMU,TGRAZE,ZRT,ZMT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ZOORM,ZOORMR,ZOORMF
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: agzt
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: AGZ, ZGZ ! OMNIVOROUS ZOOPLANKTON
END MODULE ZOOPLANKTONC


!***********************************************************************************************************************************
!**                                          C E - Q U A L - W 2  P R E P R O C E S S O R                                         **
!***********************************************************************************************************************************

PROGRAM W2_PRE
  USE CONTROL; USE GRIDC; USE KINETICS; USE PREC; USE FILE_ERROR; use macrophytec; use zooplanktonc


! Variable declarations

  REAL          :: JDAY,   JDAYO
  INTEGER       :: WIN=6
  INTEGER       :: CON,    WRN,    ERR,    NBR, NPROC
  INTEGER       :: YEAR,   CUF
  INTEGER       :: NWB,    NTR,    NSP,    NWD,    NPU,    NGT,    NPI
  INTEGER       :: NEP,    NEPT
  CHARACTER(1)  :: ESC, char1
  CHARACTER(2)  :: DEG
  CHARACTER(3)  :: BRA,    TRA
  CHARACTER(5)  :: FMTA
  CHARACTER(7)  :: BLANK
  CHARACTER(8)  :: AID,     CCC,    LIMC, CLOSEC
  CHARACTER(8)  :: RSIC,   RSOC,   TSRC,   WDOC
  CHARACTER(11) :: FMTI
  CHARACTER(13) :: FMTF
  CHARACTER(18) :: FMTFI,  FMTIF,  FMT2FI, FMTF2I
  CHARACTER(21) :: FMT2IF
  CHARACTER(72) :: BTHFN1, METFN1, QOTFN1
  CHARACTER(72) :: QINFN1, QTRFN1, QDTFN1, PREFN1, EUHFN1, EDHFN1
  CHARACTER(72) :: TINFN1, TTRFN1, TDTFN1, TPRFN1, TUHFN1, TDHFN1
  CHARACTER(72) :: CINFN1, CTRFN1, CDTFN1, CPRFN1, CUHFN1, CDHFN1
  CHARACTER(72) :: SNPFN1, PRFFN1, VPLFN1, CPLFN1, SPRFN1, FLXFN1
  CHARACTER(72) :: QWDFN,  RSIFN,  QGTFN,  TSRFN,  WDOFN
  CHARACTER(72) :: SHDFN,  WSCFN                                                                                       !SW 04/03/02
  CHARACTER(72) :: CONFN='w2_con.npt',INIFN='pre.opt',WRNFN='pre.wrn',ERRFN='pre.err'
  CHARACTER(80) :: HEADER(3)
  LOGICAL       :: TRIBUTARIES,      WITHDRAWALS
  LOGICAL       :: PIPES,            GATES,            SPILLWAYS,   PUMPS
  LOGICAL       :: DISSOLVED_SOLIDS, ALGAE,            CBO_DEMAND,  GENERIC_CONST
  LOGICAL       :: LABILE_DOM,       REFRACTORY_DOM,   LABILE_POM,  REFRACTORY_POM
  LOGICAL       :: PHOSPHORUS,       AMMONIUM,         NITRATE,     IRON
  LOGICAL       :: DISSOLVED_OXYGEN, SUSPENDED_SOLIDS, TOT_IC
  LOGICAL       :: ACTIVE_SEGMENT,   RESTART_IN


! Allocation declarations


  REAL,          ALLOCATABLE, DIMENSION(:)     :: ESP,    A1SP,   B1SP,   A2SP,   B2SP,   ETUSP,  EBUSP,  ETDSP,  EBDSP
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AGASSP, BGASSP, CGASSP
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EUPI,   EDPI,   WPI,    DLXPI,  FPI,    FMINPI, ETUPI,  EBUPI,  ETDPI,  EBDPI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EPU,    STRTPU, ENDPU,  EONPU,  EOFFPU, QPU,    ETPU,   EBPU
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EGT,    ETUGT,  EBUGT,  ETDGT,  EBDGT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: GTA1,   GTB1,   GTA2,   GTB2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: A1GT,   B1GT,   G1GT,   A2GT,   B2GT,   G2GT,    BGATE
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AGASGT, BGASGT, CGASGT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETTR,   EBTR,   EWD
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ICEI,   ALBEDO, HWI,    BETAI,  GAMMAI, ICEMIN, ICET2,  EXH2O,  BETA,   EXOM, Z0    ! SW 3.6
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EXSS,   SSS,    COLQ10, COLDK,  DX,     CBHE, TAUCR         !TC 08/20/03 SW 1/16/04
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AFW,    BFW,    CFW,    WINDH,  TSED,   FI,     TSEDF, WSC
  REAL,          ALLOCATABLE, DIMENSION(:)     :: SEDCI,  SEDDK,  SODT1,  SODT2,  SODK1,  SODK2,  FSOD,   FSED, seds,sedbr  !SW 6/1/07
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AX,     AZMAX,  T2I,    ELBOT,  SLOPE,  LAT,    LONG
  REAL,          ALLOCATABLE, DIMENSION(:)     :: DLX,    ELWS,   SHADE,  PHI0,   FRIC
  REAL,          ALLOCATABLE, DIMENSION(:)     :: XBR,    XGR,    XI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QTRAVB, QTRMXB, QINAV,  QPRAV,  QDTAV,  QSTRAV
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QSTR,   QSTRO
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINS,   QPRS,   QDTS,   QSTRS
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINMXB, QOTMXB, QINAVB, QOTAVB, QWDMXB, QWDAVB,  QDNAVB
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINMX,  QDTMX,  QSTRMX, QPRMX
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINAVW, QOTAVW, QINMXW, QOTMXW
  REAL,          ALLOCATABLE, DIMENSION(:)     :: RSOD,   DLTD,   TSRKTD
  REAL,          ALLOCATABLE, DIMENSION(:)     :: RSOF,   DLTF,   TSRKTF
  REAL,          ALLOCATABLE, DIMENSION(:)     :: DLTMAX, THETA
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CIN,    CIN2                                                         !SW 01/07/01
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CGQ10,  CG0DK,  CG1DK,  CGS
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AG,     AR,     AE,     AM,     AS,     EXA,    ASAT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AHSP,   AHSN,   AHSSI,  AT1,    AT2,    AT3,    AT4
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AK1,    AK2,    AK3,    AK4
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ALGP,   ALGN,   ALGSI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EG,     ER,     EE,     EM,     EB
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EHSN,   EHSP,   EHSSI,  ESAT,   EHS,    ENPR
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EP,     EN,     EC,     ESI,    ECHLA,  EPOM
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ET1,    ET2,    ET3,    ET4
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EK1,    EK2,    EK3,    EK4
  REAL,          ALLOCATABLE, DIMENSION(:)     :: BODP,   BODN,   BODC                                                 !TC 01/15/02
  REAL,          ALLOCATABLE, DIMENSION(:)     :: LDOMDK, RDOMDK, LRDDK
  REAL,          ALLOCATABLE, DIMENSION(:)     :: LPOMDK, RPOMDK, LRPDK,  POMS,   APOM,   ANPR
  REAL,          ALLOCATABLE, DIMENSION(:)     :: OMT1,   OMT2,   OMK1,   OMK2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: KBOD,   TBOD,   RBOD
  REAL,          ALLOCATABLE, DIMENSION(:)     :: PO4R
  REAL,          ALLOCATABLE, DIMENSION(:)     :: NH4R,   NH4DK
  REAL,          ALLOCATABLE, DIMENSION(:)     :: NH4T1,  NH4T2,  NH4K1,  NH4K2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: NO3DK,  NO3S,FNO3SED
  REAL,          ALLOCATABLE, DIMENSION(:)     :: NO3T1,  NO3T2,  NO3K1,  NO3K2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: DSIR,   PSIS,   PSIDK,  PARTSI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: FER,    FES
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CO2R
  REAL,          ALLOCATABLE, DIMENSION(:)     :: O2NH4,  O2AR,   O2AG
  REAL,          ALLOCATABLE, DIMENSION(:)     :: O2EG,   O2ER
  REAL,          ALLOCATABLE, DIMENSION(:)     :: RCOEF1, RCOEF2, RCOEF3, RCOEF4
  REAL,          ALLOCATABLE, DIMENSION(:)     :: SOD,ARODI,STRCKLR,BOUNDFR,AZE
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETSR
  REAL,          ALLOCATABLE, DIMENSION(:)     :: TSRD,   TSRF,   WDOD,   WDOF
  REAL,          ALLOCATABLE, DIMENSION(:)     :: TTLB,   TTRB,   CLLB,   CLRB,   SRLB1,  SRRB1,  SRLB2,  SRRB2,  ANG
  REAL,          ALLOCATABLE, DIMENSION(:)     :: SRFJDAY1, SRFJDAY2, SHADEIN,cinnload,cinpload,ctrnload,ctrpload
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: EPICI
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: TOPO                                                                 !SW 04/03/02
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: JDQIN,  JDQDT,  JDQTR,  JDQPR,  JDQWD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: QIN,    QDT,    QTR,    QPR,    QWD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: B,      C2IWB
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SAGR,   SABR,   CVGR,   CVBR,   HGR,    HBR,    BGR,    BBR
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: ESTR,   WSTR
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SNPD,   SCRD,   PRFD,   SPRD,   CPLD,   VPLD,   FLXD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SNPF,   SCRF,   PRFF,   SPRF,   CPLF,   VPLF,   FLXF
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: QTRS,   QTRAV,  QTRMX,  QTRMXS
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: QWDS,   QWDAV,  QWDMX
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: CINAVG, CINMIN, CINMAX, CTRAVG, CTRMIN, CTRMAX,cinload,ctrload ,cdtload                  !SW 01/07/01
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: TVP
  REAL,          ALLOCATABLE, DIMENSION(:,:,:) :: CVP,    EPIVP
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: VOLB,   VOLG
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: BTH,    VPR,    LPR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NSNP,   NSCR,   NSPR,   NPRF,   NVPL,   NFLX,   NCPL               !TC 06/11/02
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NQIN,   NQDT,   NQTR,   NQPR,   NQOT
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NISNP,  NIPRF,  NISPR,  KBMAX
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IUPI,   IDPI,   KTUPI,  KBUPI,  KTDPI,  KBDPI
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUSP,  KBUSP,  KTDSP,  KBDSP,  IUSP,   IDSP
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUGT,  KBUGT,  KTDGT,  KBDGT,  IUGT,   IDGT
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IUPU,   IDPU,   KTPU,   KBPU
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: EQSP,   EQGT,   NEQN,   ANEQN,  ENEQN
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTI
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: UHS,    DHS,    UQB,    DQB,    JBDN,   NPOINT, NSTR, NL
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NTBR,   NWBR,   NUPIBR, NDPIBR, NUGTBR, NDGTBR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: JBTR,   JBWD,   JBUPI,  JBDPI,  JBUGT,  JBDGT
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWD,    KTWD,   KBWD
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWR,    KTWR,   KBWR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: ITR,FBC
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWDO,   ITSR
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: INCN,   TRCN, incdt                                                         !SW 01/07/01
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: ISNP,   IPRF,   ISPR
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: NCCGR,  NCCBR
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: KTSTR,  KBSTR
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: EXC,    EXIC                                                         !SW 12/04/01
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: EPIC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: VISC,   CELC,   FETCHC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: AZC,    AZSLC,  WTYPEC, BK,     SLICEC, REAERC, DTRC, TKECAL
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: GASSPC, GASGTC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: CAC,    CNAME2
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: PTRC,   PPUC,   PUGTC,  PDGTC,  PDSPC,  PUSPC,  PDPIC,  PUPIC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: ICEC,   RHEVC,  SEDC,   SDC,    PRNSC, DYNSEDK
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SNPC,   SCRC,   SPRC,   VPLC,   CPLC,   PRFC,   FLXC, TCPL
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: QINIC,  TRIC,   DTRIC,  WDIC,   HDIC,   METIC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SLTRC,  SLHTC,  FRICC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: WINDC,  QINC,   QOUTC,  HEATC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: VBC,    MBC,    EBC,    PQC,    EVC,    PRC
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: BTHFN,  METFN,  VPRFN,  LPRFN,  QOTFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: QINFN,  QTRFN,  QDTFN,  PREFN,  EUHFN,  EDHFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: TINFN,  TTRFN,  TDTFN,  TPRFN,  TUHFN,  TDHFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: CINFN,  CTRFN,  CDTFN,  CPRFN,  CUHFN,  CDHFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: SNPFN,  PRFFN,  VPLFN,  CPLFN,  SPRFN,  FLXFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: EXTFN                                                                !TC 12/12/01
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: HPLTC,  CPLTC,  CDPLTC            ! SW 1/16/04
  CHARACTER(4),  ALLOCATABLE, DIMENSION(:)     :: CUNIT1                            ! SW 1/16/04
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SEDRC                             ! SW 1/16/04
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: LATGTC, LATSPC, LATPIC, LATPUC, DYNGTC

  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: SINKC,  STRIC,  CINBRC, CINTRC, CDTBRC, CPRBRC, CFWBC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: EPIPRC
  LOGICAL                                      :: CONSTITUENTS,     BRANCH_FOUND                                       !SW 10/17/02
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: DIST_TRIB, ZERO_SLOPE                                                !SW 10/17/02
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: SHIFT_DEMAND, SEDIMENT_CALC, PRECIPITATION
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: UH_EXTERNAL,  DH_EXTERNAL,   UQ_EXTERNAL, DQ_EXTERNAL, DAM_FLOW
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: UH_INTERNAL,  DH_INTERNAL,   UQ_INTERNAL, DQ_INTERNAL
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: OPEN_VPR,     OPEN_LPR, OPEN_LPRC
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: ISO_TEMP,     VERT_TEMP,     LONG_TEMP,   LONG_SED
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: ISO_CONC,     VERT_CONC,     LONG_CONC
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: ISO_EPIPHYTON, VERT_EPIPHYTON, LONG_EPIPHYTON, EPIPHYTON_CALC

! Data declarations

  DATA BLANK /'       '/
  DATA DMO2  /2.04E-9/,  DMCO2 /1.63E-9/
  DATA CON   /21/, NPT /22/, ERR /370/, WRN /371/, NEW /372/
  DATA NDC   /23/, NHY /15/, NFL /73/
  DATA IANGLES /18/                                                                                                    !SW 04/03/02


! Variable initializations

  FMTA       = '(A)'
  DEG        =  CHAR(248)//'C'
  ESC        =  CHAR(027)
  UNIT       =  CON
  DELETE_ERR = .TRUE.
  DELETE_WRN = .TRUE.

! Initialize input/output files

  OPEN (CON,FILE=CONFN,STATUS='OLD',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open control file [CONFN="w2_con.npt"]'
    STOP
  END IF
  OPEN (ERR,FILE=ERRFN,STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open preprocessor error output file [ERRFN="pre.err"]'
    STOP
  END IF
  OPEN (WRN,FILE=WRNFN,STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open preprocessor warning output file [WRNFN="pre.wrn"]'
    STOP
  END IF
  OPEN (INI,FILE=INIFN,STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open preprocessor output file [INIFN="pre.opt"]'
    STOP
  END IF

!***********************************************************************************************************************************
!*                                                     Task 1:  Read all inputs                                                   **
!***********************************************************************************************************************************

! Title cards

  WRITE (WIN,*) 'Control file'
  ALLOCATE (TITLE(11))
  WRITE (WIN,*) '  title cards'
  READ (CON,'(//A8/(8X,A72))',ERR=400) AID, (TITLE(J),J=1,10)
  IF (AID /= 'TITLE C ')               GO TO 400

! Array dimensions

  WRITE (WIN,*) '  array dimensions'
  READ (CON,'(/A8/(8X,5I8,A8))')       AID, NWB, NBR, IMX, KMX, NPROC, CLOSEC
  IF (AID /= 'GRID    ')               GO TO 400
  READ (CON,'(/A8/(8X,10I8))')         AID, NTR, NST, NIW, NWD,  NGT, NSP, NPI, NPU
  IF (AID /= 'IN/OUTFL')               GO TO 400
  READ (CON,'(/A8/(8X,10I8))')         AID, NGC, NSS, NAL, NEP,  NBOD, NMC, NZP
  IF (AID /= 'CONSTITU')               GO TO 400
  READ (CON,'(/A8/(8X,10I8))')         AID, NOD
  IF (AID /= 'MISCELL ')               GO TO 400

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
  NZOOS = NALK+1
  NZOOE = NZOOS+NZP-1
  NLDOMP= NZOOE+1
  NRDOMP= NLDOMP+1
  NLPOMP= NRDOMP+1
  NRPOMP= NLPOMP+1
  NLDOMN= NRPOMP+1
  NRDOMN= NLDOMN+1
  NLPOMN= NRDOMN+1
  NRPOMN= NLPOMN+1


! Constituent, tributary, and widthdrawal totals

  NCT  = NRPOMN
  NTRT = NTR+NGT+NSP+NPI+NPU
  NWDT = NWD+NGT+NSP+NPI+NPU
  NEPT = MAX(NEP,1)
  Nmct = MAX(nmc,1)
  NZPt  = MAX(nzp,1)                                                                                    ! SW 6/1/07

! Allocation declarations

  ALLOCATE (EXC(NWB),    EXIC(NWB))
  ALLOCATE (SSAC(NSS),   ALAC(NAL))
  ALLOCATE (CUNIT(NCT),  CNAME(NCT),  CNAME1(NCT), CNAME2(NCT), CDNAME(NDC), NADC(NWB), CDMULT(NDC), CMULT(NCT), CUNIT1(NCT), CUNIT2(NCT))    ! SW 1/16/04 8/25/05
  ALLOCATE (DMULT(NDC))
  ALLOCATE (HNAME(NHY),  HFMT(NHY))
  ALLOCATE (SLTRC(NWB),  THETA(NWB),  FRICC(NWB))
  ALLOCATE (WINDC(NWB),  QINC(NWB),   QOUTC(NWB),  HEATC(NWB),  SLHTC(NWB))
  ALLOCATE (QINIC(NBR),  DTRIC(NBR),  TRIC(NTR),   WDIC(NWD),   HDIC(NBR),   METIC(NWB))
  ALLOCATE (VBC(NWB),    EBC(NWB),    MBC(NWB),    PQC(NWB),    EVC(NWB),    PRC(NWB))
  ALLOCATE (ICEC(NWB),   SLICEC(NWB), ICEI(NWB),   ALBEDO(NWB), HWI(NWB),    BETAI(NWB),  GAMMAI(NWB), ICEMIN(NWB), ICET2(NWB))
  ALLOCATE (EXH2O(NWB),  BETA(NWB),   EXOM(NWB),   EXSS(NWB),   SROC(NWB), Z0(NWB))
  ALLOCATE (CBHE(NWB),   TSED(NWB),   FI(NWB),     TSEDF(NWB),  AFW(NWB),    BFW(NWB),    CFW(NWB),    WINDH(NWB),  RHEVC(NWB))
  ALLOCATE (SEDC(NWB),   PRNSC(NWB),  SEDCI(NWB),  SEDDK(NWB),  FSOD(NWB),   FSED(NWB), seds(nwb), sedbr(nwb),DYNSEDK(NWB))   ! SW 6/1/07
  ALLOCATE (SODT1(NWB),  SODT2(NWB),  SODK1(NWB),  SODK2(NWB))
  ALLOCATE (FETCHC(NWB), SDC(NWB),    DX(NWB))
  ALLOCATE (SNPC(NWB),   PRFC(NWB),   SPRC(NWB),   CPLC(NWB),   VPLC(NWB),   FLXC(NWB),   SCRC(NWB))
  ALLOCATE (NSNP(NWB),   NPRF(NWB),   NSPR(NWB),   NCPL(NWB),   NVPL(NWB),   NFLX(NWB),   NSCR(NWB), TCPL(NWB))
  ALLOCATE (NISNP(NWB),  NIPRF(NWB),  NISPR(NWB))
  ALLOCATE (SNPFN(NWB),  PRFFN(NWB),  SPRFN(NWB),  VPLFN(NWB),  CPLFN(NWB),  FLXFN(NWB))
  ALLOCATE (BTHFN(NWB),  METFN(NWB),  VPRFN(NWB),  LPRFN(NWB))
  ALLOCATE (QINFN(NBR),  TINFN(NBR),  CINFN(NBR),  CDHFN(NBR),  QDTFN(NBR),  TDTFN(NBR),  CDTFN(NBR),  PREFN(NBR),  TPRFN(NBR))
  ALLOCATE (CPRFN(NBR),  EUHFN(NBR),  TUHFN(NBR),  CUHFN(NBR),  EDHFN(NBR),  TDHFN(NBR),  QOTFN(NBR))
  ALLOCATE (EXTFN(NWB))                                                                                                !TC 12/12/01
  ALLOCATE (BODP(NBOD),  BODN(NBOD),  BODC(NBOD))                                                                      !TC 01/15/02
  ALLOCATE (DTRC(NBR),   ALPHA(NBR))
  ALLOCATE (AX(NWB),     WTYPEC(NWB), JBDN(NWB),   AZSLC(NWB),  AZMAX(NWB),  KBMAX(NWB),  AZC(NWB),FBC(NWB),AZE(NWB),ARODI(NWB),STRCKLR(NWB),BOUNDFR(NWB),TKECAL(NWB))
  ALLOCATE (VISC(NWB),   CELC(NWB))
  ALLOCATE (BTH(NWB),    VPR(NWB),    LPR(NWB))
  ALLOCATE (T2I(NWB),    KTWB(NWB))                                                                                    !TC 06/11/02
  ALLOCATE (QINAVW(NWB), QOTAVW(NWB), QINMXW(NWB), QOTMXW(NWB))
  ALLOCATE (LAT(NWB),    LONG(NWB),   ELBOT(NWB),  BS(NWB),     BE(NWB))
  ALLOCATE (IEPR(NWB),   KEPR(NWB))
  ALLOCATE (REAERC(NWB), NEQN(NWB))
  ALLOCATE (COLQ10(NWB), COLDK(NWB))
  ALLOCATE (LDOMDK(NWB), RDOMDK(NWB), LRDDK(NWB))
  ALLOCATE (LPOMDK(NWB), RPOMDK(NWB), LRPDK(NWB),  POMS(NWB),   APOM(NAL),   ANPR(NAL),   ANEQN(NAL))
  ALLOCATE (ORGP(NWB),   ORGN(NWB),   ORGC(NWB),   ORGSI(NWB))
  ALLOCATE (OMT1(NWB),   OMT2(NWB),   OMK1(NWB),   OMK2(NWB))
  ALLOCATE (KBOD(NBOD),  TBOD(NBOD),  RBOD(NBOD))                                                                      !TC 09/23/02
  ALLOCATE (PO4R(NWB),   PARTP(NWB))
  ALLOCATE (NH4R(NWB),   NH4DK(NWB))
  ALLOCATE (NH4T1(NWB),  NH4T2(NWB),  NH4K1(NWB),  NH4K2(NWB))
  ALLOCATE (NO3DK(NWB),  NO3S(NWB),FNO3SED(NWB))
  ALLOCATE (NO3T1(NWB),  NO3T2(NWB),  NO3K1(NWB),  NO3K2(NWB))
  ALLOCATE (DSIR(NWB),   PSIS(NWB),   PSIDK(NWB),  PARTSI(NWB))
  ALLOCATE (FER(NWB),    FES(NWB))
  ALLOCATE (CO2R(NWB))
  ALLOCATE (O2NH4(NWB),  O2OM(NWB),   O2AR(NAL),   O2AG(NAL))                                                          !TC 09/23/02
  ALLOCATE (CGQ10(NGC),  CG0DK(NGC),  CG1DK(NGC),  CGS(NGC))
  ALLOCATE (CAC(NCT))
  ALLOCATE (RCOEF1(NWB), RCOEF2(NWB), RCOEF3(NWB), RCOEF4(NWB))
  ALLOCATE (CIN(NCT),    CIN2(NCT))                                                                                    !SW 01/07/01
  ALLOCATE (TTLB(IMX),   TTRB(IMX),   CLLB(IMX),   CLRB(IMX),   SRLB1(IMX),  SRRB1(IMX),   SRLB2(IMX),  SRRB2(IMX))
  ALLOCATE (SRFJDAY1(IMX),SHADEIN(IMX),SRFJDAY2(IMX),TOPO(IMX,IANGLES),ANG(IANGLES))                                   !SW 04/03/02
  ALLOCATE (VOLB(NBR),   VOLG(NWB))
  ALLOCATE (XBR(NBR),    XGR(NWB),    XI(IMX))
  ALLOCATE (US(NBR),     DS(NBR),     CUS(NBR),    UHS(NBR),    DHS(NBR),    UQB(NBR),     DQB(NBR))
  ALLOCATE (NL(NBR),     NPOINT(NBR), SLOPE(NBR))
  ALLOCATE (NWBR(NBR),   NTBR(NBR),   NSTR(NBR),   NUPIBR(NBR), NDPIBR(NBR), NUGTBR(NBR),  NDGTBR(NBR))
  ALLOCATE (QTRAVB(NBR), QTRMXB(NBR))
  ALLOCATE (QINAV(NBR),  QPRAV(NBR),  QDTAV(NBR),  QSTRAV(NBR))
  ALLOCATE (QINS(NBR),   QPRS(NBR),   QDTS(NBR),   QSTRS(NBR))
  ALLOCATE (QINMX(NBR),  QDTMX(NBR),  QSTRMX(NBR), QPRMX(NBR))
  ALLOCATE (QINMXB(NBR), QOTMXB(NBR))
  ALLOCATE (QINAVB(NBR), QOTAVB(NBR))
  ALLOCATE (QWDAVB(NBR), QWDMXB(NBR))
  ALLOCATE (QDNAVB(NBR))
  ALLOCATE (NQIN(NBR),   NQDT(NBR),   NQPR(NBR),   NQOT(NBR))
  ALLOCATE (Z(IMX),      DLX(IMX),    KB(IMX),     BK(IMX))
  ALLOCATE (ICETH(IMX),  KTI(IMX),    ELWS(IMX),   SHADE(IMX))
  ALLOCATE (SOD(IMX),    PHI0(IMX),   FRIC(IMX))
  ALLOCATE (CDN(NDC,NWB))
  ALLOCATE (CN(NCT),     INCN(NCT,NBR),   TRCN(NCT,NTR), incdt(nct,nbr))                                                               !SW 01/07/01
  ALLOCATE (SSS(NSS),    SEDRC(NSS),   TAUCR(NSS))                                                                      !TC 08/20/03 SW 1/16/04
  ALLOCATE (AG(NAL),     AR(NAL),     AE(NAL),     AM(NAL),     AS(NAL),     EXA(NAL),    ASAT(NAL),   ALGP(NAL),   ALGN(NAL))
  ALLOCATE (ALGC(NAL),   ALGSI(NAL),  ACHLA(NAL),  AHSP(NAL),   AHSN(NAL),   AHSSI(NAL))
  ALLOCATE (AT1(NAL),    AT2(NAL),    AT3(NAL),    AT4(NAL),    AK1(NAL),    AK2(NAL),    AK3(NAL),    AK4(NAL))
  ALLOCATE (EG(nept),     ER(nept),     EE(nept),     EM(nept),     EB(nept),     ESAT(nept),   EP(nept),     EN(nept))
  ALLOCATE (EC(nept),     ESI(nept),    ECHLA(nept),  EHSP(nept),   EHSN(nept),   EHSSI(nept),  EPOM(nept),   EHS(nept))
  ALLOCATE (ET1(nept),    ET2(nept),    ET3(nept),    ET4(nept),    EK1(nept),    EK2(nept),    EK3(nept),    EK4(nept))
  ALLOCATE (ENPR(nept),   ENEQN(nept),  O2ER(nept),   O2EG(nept))
  ALLOCATE (TSRKTD(NOD), TSRKTF(NOD), IWDO(NOD),   RSOD(NOD),   RSOF(NOD),   DLTD(NOD),   DLTF(NOD),   DLTMAX(NOD))
  ALLOCATE (IUPI(NPI),   IDPI(NPI),   EUPI(NPI),   EDPI(NPI),   WPI(NPI),    DLXPI(NPI),  FPI(NPI),    FMINPI(NPI), PUPIC(NPI))
  ALLOCATE (ETUPI(NPI),  EBUPI(NPI),  KTUPI(NPI),  KBUPI(NPI),  PDPIC(NPI),  ETDPI(NPI),  EBDPI(NPI),  KTDPI(NPI))
  ALLOCATE (KBDPI(NPI),  JBUPI(NPI),  JBDPI(NPI))
  ALLOCATE (PUSPC(NSP),  ETUSP(NSP),  EBUSP(NSP),  KTUSP(NSP),  KBUSP(NSP),  PDSPC(NSP), ETDSP(NSP), EBDSP(NSP))
  ALLOCATE (KTDSP(NSP),  KBDSP(NSP),  IUSP(NSP),   IDSP(NSP),   ESP(NSP),    A1SP(NSP),  B1SP(NSP),  A2SP(NSP))
  ALLOCATE (B2SP(NSP),   AGASSP(NSP), BGASSP(NSP), CGASSP(NSP), EQSP(NSP),   GASSPC(NSP))
  ALLOCATE (GTA1(NGT),   GTB1(NGT),   GTA2(NGT),   GTB2(NGT),   IUGT(NGT),   IDGT(NGT),   EGT(NGT),    A1GT(NGT),   B1GT(NGT))
  ALLOCATE (G1GT(NGT),   A2GT(NGT),   B2GT(NGT),   G2GT(NGT),   BGATE(NGT),  PUGTC(NGT),  ETUGT(NGT),  EBUGT(NGT),  KTUGT(NGT))
  ALLOCATE (KBUGT(NGT),  PDGTC(NGT),  ETDGT(NGT),  EBDGT(NGT),  KTDGT(NGT),  KBDGT(NGT),  AGASGT(NGT), BGASGT(NGT), CGASGT(NGT))
  ALLOCATE (GASGTC(NGT), EQGT(NGT),   JBUGT(NGT),  JBDGT(NGT))
  ALLOCATE (LATGTC(NGT), LATSPC(NSP), LATPIC(NPI), LATPUC(NPU), DYNGTC(NGT))
  ALLOCATE (IUPU(NPU),   IDPU(NPU),   EPU(NPU),    STRTPU(NPU), ENDPU(NPU),  EONPU(NPU),  EOFFPU(NPU), QPU(NPU),    PPUC(NPU))
  ALLOCATE (ETPU(NPU),   EBPU(NPU),   KTPU(NPU),   KBPU(NPU))
  ALLOCATE (IWR(NIW),    KTWR(NIW),   KBWR(NIW))
  ALLOCATE (IWD(NWDT),   JBWD(NWD),   EWD(NWDT),   KTWD(NWDT),  KBWD(NWDT))
  ALLOCATE (ITR(NTRT),   QTRFN(NTR),  TTRFN(NTR),  CTRFN(NTR),  ETTR(NTRT),  EBTR(NTRT),  PTRC(NTRT))
  ALLOCATE (JBTR(NTR),   NQTR(NTR),   NACTR(NTR),  NACIN(NBR),  NACPR(NWB),  NACDT(NBR))
  ALLOCATE (TSRD(NOD),   TSRF(NOD),   WDOD(NOD),   WDOF(NOD))
  ALLOCATE (T2(KMX,IMX))
  ALLOCATE (ITSR(IMX*KMX),   ETSR(IMX*KMX))
  ALLOCATE (ESTR(NST,NBR),   WSTR(NST,NBR),   SINKC(NST,NBR),  QSTR(NST),       QSTRO(NST),      STRIC(NST,NBR))
  ALLOCATE (KTSTR(NST,NBR),  KBSTR(NST,NBR))
  ALLOCATE (TDS(KMX,IMX),    COL(KMX,IMX),    NH4(KMX,IMX),    NO3(KMX,IMX),    PO4(KMX,IMX),    FE(KMX,IMX),     DSI(KMX,IMX))
  ALLOCATE (PSI(KMX,IMX),    LDOM(KMX,IMX),   RDOM(KMX,IMX),   LPOM(KMX,IMX),   RPOM(KMX,IMX),   O2(KMX,IMX),     TIC(KMX,IMX))
  ALLOCATE (ALK(KMX,IMX))
  ALLOCATE (TN(KMX,IMX),     TP(KMX,IMX))
  ALLOCATE (DON(KMX,IMX),    DOP(KMX,IMX),    DOC(KMX,IMX))
  ALLOCATE (PON(KMX,IMX),    POP(KMX,IMX),    POC(KMX,IMX))
  ALLOCATE (TON(KMX,IMX),    TOP(KMX,IMX),    TOC(KMX,IMX))
  ALLOCATE (CHLA(KMX,IMX),   ATOT(KMX,IMX))
  ALLOCATE (O2DG(KMX,IMX),   TOTSS(KMX,IMX),  TISS(KMX,IMX))
  ALLOCATE (CBODU(KMX,IMX),  TKN(KMX,IMX))
  ALLOCATE (CINBRC(NCT,NBR), CINTRC(NCT,NTR), CDTBRC(NCT,NBR), CPRBRC(NCT,NBR), CPRWBC(NCT,NWB), CFWBC(NFL,NWB),  HPRWBC(NHY,NWB))
  ALLOCATE (QTRS(NTR,NBR),   QTRAV(NTR,NBR),  QTRMXS(NTR,NBR), QTRMX(NTR,NBR),  CTRAVG(NCT,NTR), CTRMIN(NCT,NTR), CTRMAX(NCT,NTR))
  ALLOCATE (WSC(IMX),    SNPD(NOD,NWB),   SNPF(NOD,NWB),   SCRD(NOD,NWB),   SCRF(NOD,NWB),   PRFD(NOD,NWB))
  ALLOCATE (QWDS(NWD,NBR),   QWDAV(NWD,NBR),  QWDMX(NWD,NBR))
  ALLOCATE (PRFF(NOD,NWB),   SPRD(NOD,NWB),   SPRF(NOD,NWB),   CPLD(NOD,NWB),   CPLF(NOD,NWB))
  ALLOCATE (VPLD(NOD,NWB),   VPLF(NOD,NWB),   FLXD(NOD,NWB),   FLXF(NOD,NWB))
  ALLOCATE (NCCBR(KMX,NBR),  CVBR(KMX,NBR),   SABR(KMX,NBR),   HBR(KMX,NBR),    BBR(KMX,NBR))
  ALLOCATE (SAGR(KMX,NWB),   NCCGR(KMX,NWB),  CVGR(KMX,NWB),   HGR(KMX,NWB),    BGR(KMX,NWB))
  ALLOCATE (TVP(KMX,NWB),    H(KMX,NWB))
  ALLOCATE (B(KMX,IMX),      CONV(KMX,IMX),   EL(KMX,IMX))
  ALLOCATE (CDWBC(NDC,NWB))
  ALLOCATE (C2IWB(NCT,NWB),  CINAVG(NCT,NBR), CINMIN(NCT,NBR), CINMAX(NCT,NBR),cinload(nct,nbr),ctrload(nct,ntr),cdtload(nct,nbr),cinnload(nbr),cinpload(nbr),ctrnload(ntr),ctrpload(ntr))
  ALLOCATE (EPIC(NWB,NEPt), EPICI(NWB,NEPt),EPIPRC(NWB,NEPt))                                                       ! cb 6/16/06
  ALLOCATE (IPR(IMX,NWB),    ISNP(IMX,NWB),   IPRF(IMX,NWB),   ISPR(IMX,NWB))
  ALLOCATE (CVP(KMX,NCT,NWB))
  ALLOCATE (C2(KMX,IMX,NCT), CD(KMX,IMX,NDC), CG(KMX,IMX,NGC), ALG(KMX,IMX,NAL), SS(KMX,IMX,NSS), CBOD(KMX,IMX,NBOD))
  ALLOCATE (SHIFT_DEMAND(NWB), SEDIMENT_CALC(NWB))
  ALLOCATE (OPEN_VPR(NWB),     OPEN_LPR(NWB),    OPEN_LPRC(NWB),  ISO_TEMP(NWB),    VERT_TEMP(NWB))
  ALLOCATE (UH_EXTERNAL(NBR),  DH_EXTERNAL(NBR),   UH_INTERNAL(NBR), DH_INTERNAL(NBR), DAM_FLOW(NBR))
  ALLOCATE (UQ_EXTERNAL(NBR),  DQ_EXTERNAL(NBR),   UQ_INTERNAL(NBR), DQ_INTERNAL(NBR))
  ALLOCATE (ISO_CONC(NCT,NWB), VERT_CONC(NCT,NWB), LONG_TEMP(NWB),   LONG_SED(NWB),    LONG_CONC(NCT,NWB))
  ALLOCATE (ICE_CALC(NWB),     DIST_TRIB(NBR),     PRECIPITATION(NWB))
  ALLOCATE (ZERO_SLOPE(NWB))                                                                                           !SW 10/17/02
  ALLOCATE (CMIN(NCT),   CMAX(NCT),   HYMIN(NHY),  HYMAX(NHY),  CDMIN(NDC),  CDMAX(NDC))                               ! SW 1/16/04
  ALLOCATE (FMTH(NHY),   HMULT(NHY),  FMTC(NCT),   FMTCD(NDC))                                                         ! SW 1/16/04
  ALLOCATE (CPLTC(NCT),  HPLTC(NHY),  CDPLTC(NDC))                                                                     ! SW 1/16/04
! v3.5 start
  ALLOCATE (zg(NZPt),zm(NZPt),zeff(NZPt),prefp(NZPt),zr(NZPt),zoomin(NZPt),zs2p(NZPt),exz(NZPt),PREFZ(NZPt,nzpt),o2zr(nzpt))
  ALLOCATE (zt1(NZPt),zt2(NZPt),zt3(NZPt),zt4(NZPt),zk1(NZPt),zk2(NZPt),zk3(NZPt),zk4(NZPt))
  ALLOCATE (zp(NZPt),zn(NZPt),zc(NZPt))
  allocate (prefa(nal,nzpt))
  allocate (po4zr(kmx,imx),nh4zr(kmx,imx))
  allocate (zmu(kmx,imx,nzpt),tgraze(kmx,imx,nzpt),zrt(kmx,imx,nzpt),zmt(kmx,imx,nzpt)) ! MLM POINTERS:,zoo(kmx,imx,nzpt),zooss(kmx,imx,nzpt))
  allocate (zoorm(kmx,imx,nzpt),zoormr(kmx,imx,nzpt),zoormf(kmx,imx,nzpt))
  allocate (lpzooout(kmx,imx),lpzooin(kmx,imx),dozr(kmx,imx),ticzr(kmx,imx))
  allocate (agz(kmx,imx,nal,nzpt), zgz(kmx,imx,nzpt,nzpt),agzt(kmx,imx,nal)) !omnivorous zooplankton
  allocate (ORGPLD(kmx,imx), ORGPRD(kmx,imx), ORGPLP(kmx,imx), ORGPRP(kmx,imx), ORGNLD(kmx,imx), ORGNRD(kmx,imx), ORGNLP(kmx,imx))
  allocate (ORGNRP(kmx,imx))
  allocate (ldompmp(kmx,imx),ldomnmp(kmx,imx),lpompmp(kmx,imx),lpomnmp(kmx,imx),rpompmp(kmx,imx),rpomnmp(kmx,imx))
  allocate (lpzooinp(kmx,imx),lpzooinn(kmx,imx),lpzoooutp(kmx,imx),lpzoooutn(kmx,imx))
  allocate (SEDVPp(KMX,NWB),SEDVPc(KMX,NWB),SEDVPn(KMX,NWB))
  allocate (sedp(kmx,imx),sedn(kmx,imx),sedcar(kmx,imx))
  allocate (sdkv(kmx,imx),seddktot(kmx,imx))
  allocate (sedcip(nwb),sedcin(nwb),sedcic(nwb),sedcis(nwb))
  ALLOCATE (cbods(NBOD), cbodns(kmx,imx), sedcb(kmx,imx), sedcbp(kmx,imx), sedcbn(kmx,imx), sedcbc(kmx,imx))
  allocate  (print_macrophyte(nwb,nmct), macrophyte_calc(nwb,nmct),macwbc(nwb,nmct),conv2(kmx,kmx),mprwbc(nwb,nmct))
  allocate  (mac(kmx,imx,nmct), macrc(kmx,kmx,imx,nmct),mact(kmx,kmx,imx), macrm(kmx,kmx,imx,nmct), macss(kmx,kmx,imx,nmct))
  allocate  (mgr(kmx,kmx,imx,nmct),mmr(kmx,imx,nmct), mrr(kmx,imx,nmct))
  allocate  (smacrc(kmx,kmx,imx,nmct), smacrm(kmx,kmx,imx,nmct))
  allocate  (smact(kmx,kmx,imx), smac(kmx,imx,nmct))
  allocate  (mt1(nmct),mt2(nmct),mt3(nmct),mt4(nmct),mk1(nmct),mk2(nmct),mk3(nmct),mk4(nmct),mg(nmct),mr(nmct),mm(nmct))
  allocate  (mbmp(nmct), mmax(nmct), cdstem(nmct),dmv(nmct),dwsa(nmct),anorm(nmct))
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
! v3.5 end
  ALLOCATE (ISO_EPIPHYTON(NWB,NEPt), VERT_EPIPHYTON(NWB,NEPt), LONG_EPIPHYTON(NWB,NEPt), EPIPHYTON_CALC(NWB,NEPt))         !TC 10/25/02

! State variable pointers

  TDS  => C2(:,:,1);     CG   => C2(:,:,NGCS:NGCE); SS   => C2(:,:,NSSS:NSSE); PO4  => C2(:,:,NPO4);        NH4  => C2(:,:,NNH4)
  NO3  => C2(:,:,NNO3);  DSI  => C2(:,:,NDSI);      PSI  => C2(:,:,NPSI);      FE   => C2(:,:,NFE);         LDOM => C2(:,:,NLDOM)
  RDOM => C2(:,:,NRDOM); LPOM => C2(:,:,NLPOM);     RPOM => C2(:,:,NRPOM);     CBOD => C2(:,:,NBODS:NBODE); ALG  => C2(:,:,NAS:NAE)
  O2   => C2(:,:,NDO);   TIC  => C2(:,:,NTIC);      ALK  => C2(:,:,NALK);      zoo => C2(:,:,NZOOS:NZOOE)  ; LDOMP=> C2(:,:,NLDOMP)
  RDOMP=> C2(:,:,NRDOMP); LPOMP=> C2(:,:,NLPOMP); RPOMP=> C2(:,:,NRPOMP); LDOMN=> C2(:,:,NLDOMN);RDOMN=> C2(:,:,NRDOMN); LPOMN=> C2(:,:,NLPOMN); RPOMN=> C2(:,:,NRPOMN)

! Derived variable pointers

  DOC  => CD(:,:,1);  POC  => CD(:,:,2); TOC => CD(:,:,3); DON => CD(:,:,4);  PON => CD(:,:,5);  TON  => CD(:,:,6)
  TN   => CD(:,:,7);  DOP  => CD(:,:,8); POP => CD(:,:,9); TOP => CD(:,:,10); TP  => CD(:,:,11); CHLA => CD(:,:,13)
  ATOT => CD(:,:,14); O2DG => CD(:,:,15)

! Variable initializations

  DO JW=1,NWB
    KTWB(JW) = 2
    BTH(JW)  = NPT
    VPR(JW)  = NPT+1
    LPR(JW)  = NPT+2
    NPT      = NPT+3
  END DO

! Time control cards

  WRITE (WIN,*) '  time control cards'
  READ (CON,'(/A8/8X,2F8.0,I8)',       ERR=400)  AID, TMSTRT, TMEND, YEAR
  IF (AID /= 'TIME CON')                         GO TO 400
  READ (CON,'(/A8/8X,I8,F8.0)',        ERR=400)  AID, NDLT, DLTMIN
  IF (AID /= 'DLT CON ')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',       ERR=400)  AID, (DLTD(J),J=1,NDLT)
  IF (AID /= 'DLT DATE')                         GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',       ERR=400)  AID, (DLTMAX(J),J=1,NDLT)
  IF (AID /= 'DLT MAX ')                         GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',       ERR=400)  AID, (DLTF(J),J=1,NDLT)
  IF (AID /= 'DLT FRN ')                         GO TO 400
  READ (CON,'(/A8/(8X,2A8))',          ERR=400)  AID, (VISC(JW), CELC(JW), JW=1,NWB)
  IF (AID /= 'DLT LIMI')                         GO TO 400

! Grid definition cards

  WRITE (WIN,*) '  grid definition cards'
  READ (CON,'(/A8/(8X,7I8,F8.3))',    ERR=400)  AID, (US(JB),DS(JB),UHS(JB),DHS(JB),UQB(JB),DQB(JB),NL(JB),SLOPE(JB),JB=1,NBR)
  IF (AID /= 'BRANCH G')                        GO TO 400
  READ (CON,'(/A8/(8X,3F8.0,3I8))',   ERR=400)  AID, (LAT(JW),LONG(JW),ELBOT(JW),BS(JW),BE(JW),JBDN(JW),JW=1,NWB)
  IF (AID /= 'LOCATION')                        GO TO 400

! Initial condition cards

  WRITE (WIN,*) '  initial conditions cards'
  READ (CON,'(/A8/(8X,2F8.0,A8))',    ERR=400)  AID, (T2I(JW),    ICEI(JW),   WTYPEC(JW), JW=1,NWB)
  IF (AID /= 'INIT CND')                        GO TO 400
  READ (CON,'(/A8/(8X,6A8))',         ERR=400)  AID, (VBC(JW),    EBC(JW),    MBC(JW),    PQC(JW),   EVC(JW),     PRC(JW),         &
                                                      JW=1,NWB)
  IF (AID /= 'CALCULAT')                        GO TO 400
  READ (CON,'(/A8/(8X,4A8))',         ERR=400)  AID, (WINDC(JW),  QINC(JW),   QOUTC(JW),  HEATC(JW), JW=1,NWB)
  IF (AID /= 'DEAD SEA')                        GO TO 400
  READ (CON,'(/A8/(8X,3A8))',         ERR=400)  AID, (QINIC(JB),  DTRIC(JB),  HDIC(JB),   JB=1,NBR)
  IF (AID /= 'INTERPOL')                        GO TO 400
  READ (CON,'(/A8/(8X,5A8,4F8.0))',   ERR=400)  AID, (SLHTC(JW),  SROC(JW),   RHEVC(JW),  METIC(JW), FETCHC(JW), AFW(JW),          &
                                                      BFW(JW),    CFW(JW),    WINDH(JW),  JW=1,NWB)
  IF (AID /= 'HEAT EXC')                        GO TO 400
  READ (CON,'(/A8/(8X,2A8,6F8.0))',   ERR=400)  AID, (ICEC(JW),   SLICEC(JW), ALBEDO(JW), HWI(JW),   BETAI(JW),  GAMMAI(JW),       &
                                                      ICEMIN(JW), ICET2(JW),  JW=1,NWB)
  IF (AID /= 'ICE COVE')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,F8.0))',     ERR=400)  AID, (SLTRC(JW),  THETA(JW),  JW=1,NWB)
  IF (AID /= 'TRANSPOR')                        GO TO 400
  READ (CON,'(/A8/(8X,6F8.0,A8,F8.0))',    ERR=400)  AID, (AX(JW), DX(JW), CBHE(JW), TSED(JW), FI(JW), TSEDF(JW), FRICC(JW),Z0(JW), JW=1,NWB)
  IF (AID /= 'HYD COEF')                        GO TO 400

! Viscosity cards

  WRITE (WIN,*) '  viscosity/friction cards'
  READ (CON,'(/A8/(8X,2A8,F8.0,I8,4F8.0,A8))',    ERR=400)  AID, (AZC(JW), AZSLC(JW), AZMAX(JW),FBC(JW),AZE(JW),ARODI(JW),STRCKLR(JW),BOUNDFR(JW),TKECAL(JW),JW=1,NWB)
  IF (AID /= 'EDDY VIS')                        GO TO 400

! Inflow-outflow cards

  WRITE (WIN,*) '  inflow/outflow cards'
  READ (CON,'(/A8/(8X,I8))',          ERR=400)  AID, (NSTR(JB),    JB=1,NBR)
  IF (AID /= 'N STRUC ')                        GO TO 400
  READ (CON,'(/A8/(8X,9A8))',         ERR=400)  AID, (STRIC(JS,1), JS=1,NSTR(1))
  IF (AID /= 'STR INT ')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9A8)',            ERR=400) (STRIC(JS,JB),      JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KTSTR(JS,1), JS=1,NSTR(1))
  IF (AID /= 'STR TOP ')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9I8)',            ERR=400) (KTSTR(JS,JB),      JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KBSTR(JS,1), JS=1,NSTR(1))
  IF (AID /= 'STR BOT ')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9I8)',            ERR=400) (KBSTR(JS,JB),      JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (SINKC(JS,1), JS=1,NSTR(1))
  IF (AID /= 'STR SINK')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9A8)',            ERR=400) (SINKC(JS,JB),      JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (ESTR(JS,1),  JS=1,NSTR(1))
  IF (AID /= 'STR ELEV')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9F8.0)',          ERR=400) (ESTR(JS,JB),       JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (WSTR(JS,1),  JS=1,NSTR(1))
  IF (AID /= 'STR WIDT')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9F8.0)',          ERR=400) (WSTR(JS,JB),       JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(8X,2I8,6F8.0,a8))',   ERR=400)  AID, (IUPI(J),   IDPI(J),  EUPI(J),   EDPI(J),   WPI(J),    DLXPI(J),  FPI(J),     &
                                                      FMINPI(J), latpic(j),J=1,NPI)
  IF (AID /= 'PIPES   ')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,2F8.0,2I8))',ERR=400)  AID, (PUPIC(J),  ETUPI(J), EBUPI(J),  KTUPI(J),  KBUPI(J),  J=1,NPI)
  IF (AID /= 'PIPE UP ')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,2F8.0,2I8))',ERR=400)  AID, (PDPIC(J),  ETDPI(J), EBDPI(J),  KTDPI(J),  KBDPI(J),  J=1,NPI)
  IF (AID /= 'PIPE DOW')                        GO TO 400
  READ (CON,'(/A8/(:8X,2I8,5F8.0,a8))',  ERR=400)  AID, (IUSP(J),   IDSP(J),  ESP(J),    A1SP(J),   B1SP(J),   A2SP(J),   B2SP(J),    &
                                                      latSPc(j), J=1,NSP)
  IF (AID /= 'SPILLWAY')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,2F8.0,2I8))',ERR=400)  AID, (PUSPC(J),  ETUSP(J), EBUSP(J),  KTUSP(J),  KBUSP(J),  J=1,NSP)
  IF (AID /= 'SPILL UP')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,2F8.0,2I8))',ERR=400)  AID, (PDSPC(J),  ETDSP(J), EBDSP(J),  KTDSP(J),  KBDSP(J),  J=1,NSP)
  IF (AID /= 'SPILL DO')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,I8,3F8.0))', ERR=400)  AID, (GASSPC(J), EQSP(J),  AGASSP(J), BGASSP(J), CGASSP(J), J=1,NSP)
  IF (AID /= 'SPILL GA')                        GO TO 400
  READ (CON,'(/A8/(8X,2I8,7F8.0,a8))',   ERR=400)  AID, (IUGT(J),   IDGT(J),  EGT(J),    A1GT(J),   B1GT(J),   G1GT(J),   A2GT(J),    &
                                                      B2GT(J),   G2GT(J), latgtc(j), J=1,NGT)   ! sw 7/14/06
  IF (AID /= 'GATES   ')                        GO TO 400
  READ (CON,'(/A8/(8X,4F8.0,a8))',       ERR=400)  AID, (GTA1(J),   GTB1(J),  GTA2(J),   GTB2(J),  dyngtc(j),        J=1,NGT)  ! sw 7/14/06
  IF (AID /= 'GATE WEI')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,2F8.0,2I8))',ERR=400)  AID, (PUGTC(J),  ETUGT(J), EBUGT(J),  KTUGT(J),  KBUGT(J),  J=1,NGT)
  IF (AID /= 'GATE UP ')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,2F8.0,2I8))',ERR=400)  AID, (PDGTC(J),  ETDGT(J), EBDGT(J),  KTDGT(J),  KBDGT(J),  J=1,NGT)
  IF (AID /= 'GATE DOW')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,I8,3F8.0))', ERR=400)  AID, (GASGTC(J), EQGT(J),  AGASGT(J), BGASGT(J), CGASGT(J), J=1,NGT)
  IF (AID /= 'GATE GAS')                        GO TO 400
  READ (CON,'(/A8/(8X,2I8,6F8.0,a8))',   ERR=400)  AID, (IUPU(J),   IDPU(J),  EPU(J),    STRTPU(J), ENDPU(J),  EONPU(J),  EOFFPU(J),  &
                                                      QPU(J),  latpuc(j),   J=1,NPU)                                                           ! RA 8/14/07
  IF (AID /= 'PUMPS 1')                         GO TO 400
  READ (CON,'(/A8/(8X,A8,2F8.0,2I8))',ERR=400)  AID, (PPUC(J),   ETPU(J),  EBPU(J),   KTPU(J),   KBPU(J),  J=1,NPU)
  IF (AID /= 'PUMPS 2')                         GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (IWR(JW),  JW=1,NIW)
  IF (AID /= 'WEIR SEG')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KTWR(JW), JW=1,NIW)
  IF (AID /= 'WEIR TOP')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KBWR(JW), JW=1,NIW)
  IF (AID /= 'WEIR BOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (WDIC(JW), JW=1,NWD)
  IF (AID /= 'WD INT  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (IWD(JW),  JW=1,NWD)
  IF (AID /= 'WD SEG  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (EWD(JW),  JW=1,NWD)
  IF (AID /= 'WD ELEV ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KTWD(JW), JW=1,NWD)
  IF (AID /= 'WD TOP  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KBWD(JW), JW=1,NWD)
  IF (AID /= 'WD BOT  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (PTRC(JT), JT=1,NTR)
  IF (AID /= 'TRIB PLA')                        GO TO 400
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (TRIC(JW), JW=1,NTR)
  IF (AID /= 'TRIB INT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (ITR(JT),  JT=1,NTR)
  IF (AID /= 'TRIB SEG')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (ETTR(JT), JT=1,NTR)
  IF (AID /= 'TRIB TOP')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (EBTR(JT), JT=1,NTR)
  IF (AID /= 'TRIB BOT')                        GO TO 400
  READ (CON,'(/A8/(8X,A8))',          ERR=400)  AID, (DTRC(JB), JB=1,NBR)
  IF (AID /= 'DST TRIB')                        GO TO 400

! Output control cards (excluding constituents)

  WRITE (WIN,*) '  output control cards'
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (HPRWBC(1,JW),JW=1,NWB)
  IF (AID /= 'HYD PRIN')                        GO TO 400
  DO JH=2,NHY
    READ (CON,'(:8X,9A8)',            ERR=400)  (HPRWBC(JH,JW), JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,A8,2I8))',     ERR=400)  AID, (SNPC(J), NSNP(J), NISNP(J), J=1,NWB)
  IF (AID /= 'SNP PRIN')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SNPD(J,1),J=1,NSNP(1))
  IF (AID /= 'SNP DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SNPD(J,JW),     J=1,NSNP(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SNPF(J,1),J=1,NSNP(1))
  IF (AID /= 'SNP FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SNPF(J,JW),     J=1,NSNP(JW))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (ISNP(J,1),J=1,NISNP(1))
  IF (AID /= 'SNP SEG ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9I8)',            ERR=400) (ISNP(J,JW),     J=1,NISNP(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,I8))',       ERR=400)  AID, (SCRC(JW), NSCR(JW), JW=1,NWB)
  IF (AID /= 'SCR PRIN')                        GO TO 400
  READ (CON,'(/A8/(8X,9F8.0))',       ERR=400)  AID, (SCRD(J,1),J=1,NSCR(1))
  IF (AID /= 'SCR DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SCRD(J,JW),     J=1,NSCR(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SCRF(J,1),J=1,NSCR(1))
  IF (AID /= 'SCR FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SCRF(J,JW),     J=1,NSCR(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,2I8))',      ERR=400)  AID, (PRFC(JW), NPRF(JW), NIPRF(JW), JW=1,NWB)
  IF (AID /= 'PRF PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (PRFD(J,1),J=1,NPRF(1))
  IF (AID /= 'PRF DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (PRFD(J,JW),     J=1,NPRF(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (PRFF(J,1),J=1,NPRF(1))
  IF (AID /= 'PRF FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (PRFF(J,JW),     J=1,NPRF(JW))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (IPRF(J,1),J=1,NIPRF(1))
  IF (AID /= 'PRF SEG ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9I8)',            ERR=400) (IPRF(J,JW),     J=1,NIPRF(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,2I8))',      ERR=400)  AID, (SPRC(JW), NSPR(JW), NISPR(JW), JW=1,NWB)
  IF (AID /= 'SPR PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SPRD(J,1),J=1,NSPR(1))
  IF (AID /= 'SPR DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SPRD(J,JW),     J=1,NSPR(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SPRF(J,1),J=1,NSPR(1))
  IF (AID /= 'SPR FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SPRF(J,JW),     J=1,NSPR(JW))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (ISPR(J,1),J=1,NISPR(1))
  IF (AID /= 'SPR SEG ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9I8)',            ERR=400) (ISPR(J,JW),     J=1,NISPR(JW))
  END DO
  READ (CON,'(/A8/(:8X,A8,I8))',      ERR=400)  AID, (VPLC(JW), NVPL(JW), JW=1,NWB)
  IF (AID /= 'VPL PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (VPLD(J,1),J=1,NVPL(1))
  IF (AID /= 'VPL DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (VPLD(J,JW),     J=1,NVPL(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (VPLF(J,1),J=1,NVPL(1))
  IF (AID /= 'VPL FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (VPLF(J,JW),     J=1,NVPL(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,I8,A8))',       ERR=400)  AID, (CPLC(JW), NCPL(JW),TCPL(JW), JW=1,NWB)
  IF (AID /= 'CPL PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (CPLD(J,1),J=1,NCPL(1))
  IF (AID /= 'CPL DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (CPLD(J,JW),     J=1,NCPL(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (CPLF(J,1),J=1,NCPL(1))
  IF (AID /= 'CPL FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (CPLF(J,JW),     J=1,NCPL(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,I8))',       ERR=400)  AID, (FLXC(JW), NFLX(JW), JW=1,NWB)
  IF (AID /= 'FLUXES  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (FLXD(J,1),J=1,NFLX(1))
  IF (AID /= 'FLX DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (FLXD(J,JW),     J=1,NFLX(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (FLXF(J,1),J=1,NFLX(1))
  IF (AID /= 'FLX FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (FLXF(J,JW),     J=1,NFLX(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,2I8))',      ERR=400)  AID, TSRC, NTSR, NIKTSR
  IF (AID /= 'TSR PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (TSRD(J),  J=1,NTSR)
  IF (AID /= 'TSR DATE')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (TSRF(J),  J=1,NTSR)
  IF (AID /= 'TSR FREQ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (ITSR(J),  J=1,NIKTSR)
  IF (AID /= 'TSR SEG ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (ETSR(J),  J=1,NIKTSR)
  IF (AID /= 'TSR LAYE')                        GO TO 400
  READ (CON,'(/A8/8X,A8,2I8)',        ERR=400)  AID, WDOC, NWDO, NIWDO
  IF (AID /= 'WITH OUT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (WDOD(J),  J=1,NWDO)
  IF (AID /= 'WITH DAT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (WDOF(J),  J=1,NWDO)
  IF (AID /= 'WITH FRE')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (IWDO(J),  J=1,NIWDO)
  IF (AID /= 'WITH SEG')                        GO TO 400
  READ (CON,'(/A8/8X,A8,I8,A8)',      ERR=400)  AID, RSOC, NRSO, RSIC
  IF (AID /= 'RESTART ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (RSOD(J),  J=1,NRSO)
  IF (AID /= 'RSO DATE')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (RSOF(J),  J=1,NRSO)
  IF (AID /= 'RSO FREQ')                        GO TO 400

! Constituent control cards

  WRITE (WIN,*) '  constituent control cards'
  READ (CON,'(/A8/8X,2A8,I8)',        ERR=400)  AID, CCC, LIMC, CUF
  IF (AID /= 'CST COMP')                        GO TO 400
  READ (CON,'(/A8/(:2A8))',           ERR=400)  AID, (CNAME2(JC),CAC(JC), JC=1,NCT)
  IF (AID /= 'CST ACTI')                        GO TO 400
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CDWBC(1,JW),  JW=1,NWB)
  IF (AID /= 'CST DERI')                        GO TO 400
  DO JD=2,NDC
    READ (CON,'(:8X,9A8)',            ERR=400) (CDWBC(JD,JW),       JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CFWBC(1,JW),  JW=1,NWB)
  IF (AID /= 'CST FLUX')                        GO TO 400
  DO JF=2,NFL
    READ (CON,'(:8X,9A8)',            ERR=400) (CFWBC(JF,JW),       JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (C2IWB(1,JW),  JW=1,NWB)
  IF (AID /= 'CST ICON')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9F8.0)',          ERR=400) (C2IWB(JC,JW),       JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CPRWBC(1,JW), JW=1,NWB)
  IF (AID /= 'CST PRIN')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CPRWBC(JC,JW),      JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CINBRC(1,JB), JB=1,NBR)
  IF (AID /= 'CIN CON ')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CINBRC(JC,JB),      JB=1,NBR)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CINTRC(1,JT), JT=1,NTR)
  IF (AID /= 'CTR CON ')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CINTRC(JC,JT),      JT=1,NTR)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CDTBRC(1,JB), JB=1,NBR)
  IF (AID /= 'CDT CON ')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CDTBRC(JC,JB),      JB=1,NBR)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CPRBRC(1,JB), JB=1,NBR)
  IF (AID /= 'CPR CON ')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CPRBRC(JC,JB),      JB=1,NBR)
  END DO

! Kinetics coefficients

  WRITE (WIN,*) '  kinetic coefficient cards'
  READ (CON,'(/A8/(8X,4F8.0,2A8))',     ERR=400)  AID, (EXH2O(JW),  EXSS(JW),   EXOM(JW),   BETA(JW),   EXC(JW),    EXIC(JW),      &
                                                        JW=1,NWB)
  IF (AID /= 'EX COEF ')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EXA(JA),    JA=1,NAL)
  IF (AID /= 'ALG EX  ')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EXZ(Jz),    Jz=1,nzpt)
  IF (AID /= 'ZOO EX  ')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EXM(Jm),    Jm=1,nmct)
  IF (AID /= 'MACRO EX')                          GO TO 400
    WRITE (WIN,*) '      extinction cards'
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (CGQ10(JG),  CG0DK(JG),  CG1DK(JG),  CGS(JG),    JG=1,NGC)
  IF (AID /= 'GENERIC ')                          GO TO 400
      WRITE (WIN,*) '      generic cards'
  READ (CON,'(/A8/(:8X,F8.0,A8,F8.0))',   ERR=400)  AID, (SSS(JS),    SEDRC(JS),   TAUCR(JS),  JS=1,NSS)                  ! SW 1/16/04
  IF (AID /= 'S SOLIDS')                          GO TO 400
      WRITE (WIN,*) '      iss cards'
  READ (CON,'(/A8/(8X,9F8.0))',         ERR=400)  AID, (AG(JA),     AR(JA),     AE(JA),     AM(JA),     AS(JA),     AHSP(JA),      &
                                                        AHSN(JA),   AHSSI(JA),  ASAT(JA),   JA=1,NAL)                  !TC 09/01/01
  IF (AID /= 'ALGAL RA')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID, (AT1(JA),    AT2(JA),    AT3(JA),    AT4(JA),    AK1(JA),    AK2(JA),       &
                                                        AK3(JA),    AK4(JA),    JA=1,NAL)
  IF (AID /= 'ALGAL TE')                          GO TO 400
  READ (CON,'(/A8/(8X,6F8.0,I8,F8.0))', ERR=400)  AID, (ALGP(JA),   ALGN(JA),   ALGC(JA),   ALGSI(JA),  ACHLA(JA),  APOM(JA),      &
                                                        ANEQN(JA),  ANPR(JA),   JA=1,NAL)
  IF (AID /= 'ALG STOI')                          GO TO 400
      WRITE (WIN,*) '      algae cards'
  READ (CON,'(/A8/(:8X,9A8))',          ERR=400)  AID, (EPIC(JW,1),                                              JW=1,NWB)
  IF (AID /= 'EPIPHYTE')                          GO TO 400
  DO JE=2,nept
    READ (CON,'(8X,9A8)',               ERR=400) (EPIC(JW,JE),                                                   JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',          ERR=400)  AID, (EPIPRC(JW,1),                                            JW=1,NWB)
  IF (AID /= 'EPI PRIN')                          GO TO 400
  DO JE=2,nept
    READ (CON,'(8X,9A8)',               ERR=400) (EPIPRC(JW,JE),                                                 JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EPICI(JW,1),                                             JW=1,NWB)
  IF (AID /= 'EPI INIT')                          GO TO 400
  DO JE=2,nept
    READ (CON,'(8X,9F8.0)',             ERR=400) (EPICI(JW,JE),                                                  JW=1,NWB)
  END DO
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID, (EG(JE),     ER(JE),     EE(JE),     EM(JE),     EB(JE),                    &
                                                        EHSP(JE),   EHSN(JE),   EHSSI(JE),                       JE=1,nept)
  IF (AID /= 'EPI RATE')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0,I8,F8.0))', ERR=400)  AID, (ESAT(JE),   EHS(JE),    ENEQN(JE),  ENPR(JE),            JE=1,nept)
  IF (AID /= 'EPI HALF')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID, (ET1(JE),    ET2(JE),    ET3(JE),    ET4(JE),    EK1(JE),                   &
                                                        EK2(JE),    EK3(JE),    EK4(JE),                         JE=1,nept)
  IF (AID /= 'EPI TEMP')                          GO TO 400
  READ (CON,'(/A8/(8X,6F8.0))',         ERR=400)  AID, (EP(JE),     EN(JE),     EC(JE),     ESI(JE),    ECHLA(JE),                 &
                                                        EPOM(JE),                                                JE=1,nept)
  IF (AID /= 'EPI STOI')                          GO TO 400
    WRITE (WIN,*) '      epiphyton cards'

! v3.5 start
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)   AID, (zg(jz),zr(jz),zm(jz),zeff(jz),PREFP(jz),ZOOMIN(jz),ZS2P(jz),            Jz=1,nzpt)
  IF (AID /= 'ZOOP RAT')                          GO TO 400
   IF(nzpt /= 0)then                                                         !RA 12/12/2006
    READ (CON,'(/A8/(:8X,8F8.0))',         ERR=400)  AID,     (PREFA(ja,1),                                                           Ja=1,nal)
   ELSE                                                                     !RA 12/12/2006
    READ (CON,'(/A8/(:8X))',         ERR=400)  AID                          !RA 12/12/2006
   END IF
    IF (AID /= 'ZOOP ALG')                          GO TO 400
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))',         ERR=400)       (PREFA(ja,jz),                                                           Ja=1,nal)
  end do
    READ (CON,'(/A8/(:8X,8F8.0))',         ERR=400)   AID,    (PREFz(jjz,1),                                                          Jjz=1,nzpt)
  IF (AID /= 'ZOOP ZOO')                          GO TO 400
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))',         ERR=400)       (PREFz(jjz,jz),                                                          Jjz=1,nzpt)
  end do
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (zT1(Jz),    zT2(Jz),    zT3(Jz),    zT4(Jz),    zK1(Jz),   zK2(Jz),                         &
                                       zK3(Jz),    zK4(Jz),                                                    Jz=1,nzpt)
  IF (AID /= 'ZOOP TEM')                          GO TO 400
  READ (CON,'(/A8/(8X,3F8.0))',         ERR=400)  AID,        (zP(Jz),     zN(Jz),     zC(Jz),                                         Jz=1,nzpt)      !MLM 7/11/06 3f8.0
  IF (AID /= 'ZOOP STO')                          GO TO 400
      WRITE (WIN,*) '      zooplankton cards'

  READ (CON,'(/A8/(:8X,9A8))',         ERR=400)  AID,      (macwbC(JW,1),                                                           JW=1,NWB)  ! mlm 7/11/06
    IF (AID /= 'MACROPHY')                          GO TO 400
  DO Jm=2,nmct
    READ (CON,'(8X,9A8)',         ERR=400)             (macwbC(JW,Jm),                                                          JW=1,NWB)
  END DO

  READ (CON,'(/A8/(:8X,9A8))',         ERR=400)  AID,           (mprwbC(JW,1),                                                           JW=1,NWB)
     IF (AID /= 'MAC PRIN')                          GO TO 400

  DO Jm=2,nmct
    READ (CON,'(8X,9A8)',         ERR=400)             (mprwbC(JW,Jm),                                                          JW=1,NWB)
  END DO


  READ (CON,'(/A8/(:8X,9F8.0))',         ERR=400)  AID,         (macwbCI(JW,1),                                                          JW=1,NWB)
  IF (AID /= 'MAC INI ')                          GO TO 400

  DO Jm=2,nmct
    READ (CON,'(8X,9F8.0)',         ERR=400)           (macwbcI(JW,Jm),                                                         JW=1,NWB)
  END DO

  READ (CON,'(/A8/(8X,9F8.0))',         ERR=400)  AID,         (mG(jm), mR(jm), mM(jm), msat(jm),mhsp(jm),mhsn(jm),mhsc(jm),                           &
                                          mpom(jm),lrpmac(jm),     jm=1,nmct)
    IF (AID /= 'MAC RATE')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (psed(jm), nsed(jm),                                                     jm=1,nmct)
    IF (AID /= 'MAC SED ')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (mbmp(jm), mmax(jm),                                                     jm=1,nmct)
    IF (AID /= 'MAC DIST')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (cdstem(jm), dmv(jm),dwsa(jm),anorm(jm),                                 jm=1,nmct)
    IF (AID /= 'MAC DRAG')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (mT1(Jm),    mT2(Jm),    mT3(Jm),    mT4(Jm),    mK1(Jm),   mK2(Jm),                         &
                                       mK3(Jm),    mK4(Jm),                                                    Jm=1,nmct)
    IF (AID /= 'MAC TEMP')                          GO TO 400
  READ (CON,'(/A8/(8X,6F8.0))',         ERR=400)  AID,         (mP(Jm),     mN(Jm),     mC(Jm),                                         Jm=1,nmct)
    IF (AID /= 'MAC STOI')                          GO TO 400
    WRITE (WIN,*) '      macrophyte cards'

  READ (CON,'(/A8/(8X,3F8.0))',         ERR=400)  AID, (LDOMDK(JW), RDOMDK(JW), LRDDK(JW),  JW=1,NWB)
  IF (AID /= 'DOM     ')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (LPOMDK(JW), RPOMDK(JW), LRPDK(JW),  POMS(JW),   JW=1,NWB)
  IF (AID /= 'POM     ')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (ORGP(JW),   ORGN(JW),   ORGC(JW),   ORGSI(JW),  JW=1,NWB)
  IF (AID /= 'OM STOIC')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (OMT1(JW),   OMT2(JW),   OMK1(JW),   OMK2(JW),   JW=1,NWB)
  IF (AID /= 'OM RATE ')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (KBOD(JB),   TBOD(JB),   RBOD(JB), cbods(jb),  JB=1,NBOD)
  IF (AID /= 'CBOD    ')                          GO TO 400
  READ (CON,'(/A8/(8X,3F8.0))',         ERR=400)  AID, (BODP(JB),   BODN(JB),   BODC(JB),   JB=1,NBOD)                   !TC 01/15/02
  IF (AID /= 'CBOD STO')                          GO TO 400
      WRITE (WIN,*) '      OM-BOD cards'

  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (PO4R(JW),   PARTP(JW),  JW=1,NWB)
  IF (AID /= 'PHOSPHOR')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (NH4R(JW),   NH4DK(JW),  JW=1,NWB)
  IF (AID /= 'AMMONIUM')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (NH4T1(JW),  NH4T2(JW),  NH4K1(JW),  NH4K2(JW),  JW=1,NWB)
  IF (AID /= 'NH4 RATE')                          GO TO 400
  READ (CON,'(/A8/(8X,3F8.0))',          ERR=400)  AID, (NO3DK(JW),  NO3S(JW),FNO3SED(JW), JW=1,NWB)
  IF (AID /= 'NITRATE ')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (NO3T1(JW),  NO3T2(JW),  NO3K1(JW),  NO3K2(JW),  JW=1,NWB)
  IF (AID /= 'NO3 RATE')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (DSIR(JW),   PSIS(JW),   PSIDK(JW),  PARTSI(JW), JW=1,NWB)
  IF (AID /= 'SILICA  ')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (FER(JW),    FES(JW),    JW=1,NWB)
  IF (AID /= 'IRON    ')                          GO TO 400
      WRITE (WIN,*) '      N-P-Si-Fe cards'

  READ (CON,'(/A8/(8X,F8.0))',          ERR=400)  AID, (CO2R(JW),   JW=1,NWB)
  IF (AID /= 'SED CO2 ')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (O2NH4(JW),  O2OM(JW),   JW=1,NWB)
  IF (AID /= 'STOICH 1')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (O2AR(JA),   O2AG(JA),   JA=1,NAL)
  IF (AID /= 'STOICH 2')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (O2ER(JE),   O2EG(JE),   JE=1,nept)
  IF (AID /= 'STOICH 3')                          GO TO 400
  READ (CON,'(/A8/(8X,1F8.0))',         ERR=400)  AID, (O2ZR(JZ),   JZ=1,nzpt)
  IF (AID /= 'STOICH 4')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (O2MR(Jm),   O2MG(Jm),   Jm=1,nmct)
  IF (AID /= 'STOICH 5')                          GO TO 400

  READ (CON,'(/A8/(8X,F8.0))',          ERR=400)  AID, O2LIM
  IF (AID /= 'O2 LIMIT')                          GO TO 400
      WRITE (WIN,*) '      O2 cards'

  READ (CON,'(/A8/(8X,2A8,6F8.0,A8))',     ERR=400)  AID, (SEDC(JW),   PRNSC(JW),  SEDCI(JW),  seds(jw), SEDDK(JW),  FSOD(JW),                  &
                                                        FSED(JW), sedbr(jw), DYNSEDK(JW), JW=1,NWB)                                         ! SW 6/1/07
  IF (AID /= 'SEDIMENT')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (SODT1(JW),  SODT2(JW),  SODK1(JW),  SODK2(JW),  JW=1,NWB)
  IF (AID /= 'SOD RATE')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (SOD(I),     I=1,IMX)
  IF (AID /= 'S DEMAND')                          GO TO 400
      WRITE (WIN,*) '      sediment cards'

  READ (CON,'(/A8/(8X,A8,I8,4F8.2))',   ERR=400)  AID, (REAERC(JW), NEQN(JW),   RCOEF1(JW), RCOEF2(JW), RCOEF3(JW), RCOEF4(JW),    &
                                                        JW=1,NWB)
  IF (AID /= 'REAERATI')                          GO TO 400

! Input filenames

  WRITE (WIN,*) '  input filename cards'
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, RSIFN
  IF (AID /= 'RSI FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, QWDFN
  IF (AID /= 'QWD FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, QGTFN
  IF (AID /= 'QGT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, WSCFN                                                                       !TC 03/12/02                                                                                    !elo
  IF (AID /= 'WSC FILE')              GO TO 400                                                                        !TC 03/12/02
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, SHDFN                                                                       !SW 04/03/02                                                                            !elo
  IF (AID /= 'SHD FILE')              GO TO 400                                                                        !SW 04/03/02
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (BTHFN(JW),JW=1,NWB)
  IF (AID /= 'BTH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (METFN(JW),JW=1,NWB)
  IF (AID /= 'MET FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (EXTFN(JW), JW=1,NWB)                                                       !SW 12/12/01
  IF (AID /= 'EXT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (VPRFN(JW),JW=1,NWB)
  IF (AID /= 'VPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (LPRFN(JW),JW=1,NWB)
  IF (AID /= 'LPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (QINFN(JB),JB=1,NBR)
  IF (AID /= 'QIN FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TINFN(JB),JB=1,NBR)
  IF (AID /= 'TIN FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CINFN(JB),JB=1,NBR)
  IF (AID /= 'CIN FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (QOTFN(JB),JB=1,NBR)
  IF (AID /= 'QOT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (QTRFN(JT),JT=1,NTR)
  IF (AID /= 'QTR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TTRFN(JT),JT=1,NTR)
  IF (AID /= 'TTR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CTRFN(JT),JT=1,NTR)
  IF (AID /= 'CTR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (QDTFN(JB),JB=1,NBR)
  IF (AID /= 'QDT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TDTFN(JB),JB=1,NBR)
  IF (AID /= 'TDT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CDTFN(JB),JB=1,NBR)
  IF (AID /= 'CDT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (PREFN(JB),JB=1,NBR)
  IF (AID /= 'PRE FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TPRFN(JB),JB=1,NBR)
  IF (AID /= 'TPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CPRFN(JB),JB=1,NBR)
  IF (AID /= 'CPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (EUHFN(JB),JB=1,NBR)
  IF (AID /= 'EUH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TUHFN(JB),JB=1,NBR)
  IF (AID /= 'TUH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CUHFN(JB),JB=1,NBR)
  IF (AID /= 'CUH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (EDHFN(JB),JB=1,NBR)
  IF (AID /= 'EDH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TDHFN(JB),JB=1,NBR)
  IF (AID /= 'TDH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CDHFN(JB),JB=1,NBR)
  IF (AID /= 'CDH FILE')              GO TO 400

! Output filenames

  WRITE (WIN,*) '  output filename cards'
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (SNPFN(JW),JW=1,NWB)
  IF (AID /= 'SNP FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (PRFFN(JW),JW=1,NWB)
  IF (AID /= 'PRF FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (VPLFN(JW),JW=1,NWB)
  IF (AID /= 'VPL FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CPLFN(JW),JW=1,NWB)
  IF (AID /= 'CPL FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (SPRFN(JW),JW=1,NWB)
  IF (AID /= 'SPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (FLXFN(JW),JW=1,NWB)
  IF (AID /= 'FLX FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, TSRFN
  IF (AID /= 'TSR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, WDOFN
  IF (AID /= 'WDO FILE')              GO TO 400
  CLOSE (CON)

! Bathymetry definition

! Check bathymetry
  isum=0
  DO JW=1,NWB
    isum=isum+be(jw)-bs(jw)+1
  END DO
  if(isum /= nbr)then
  CALL ERRORS
  WRITE (ERR,'(A)') 'There is an error in BS and BE in w2_con.npt: The sum of BS-BE for all waterbodies does not equal NBR.'
  endif


  DO JW=1,NWB
    WRITE (WIN,*) 'Bathymetry file:', adjustl(trim(bthfn(jw)))
    UNIT = BTH(JW)
    OPEN (BTH(JW),FILE=BTHFN(JW),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      READ  (BTH(JW),'(a1)')char1
      if(char1=='$')then
      WRITE( WIN,*) ' NEW BTH FORMAT SINCE FIRST CHARACTER==$ CSV FILE TYPE'
      WRITE (WIN,*) '  segment lengths'
      READ  (BTH(JW),*)
      READ  (BTH(JW),*,ERR=400) AID,(DLX(I),  I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) 'Finished reading:', AID
      WRITE (WIN,*) '  water surface elevations'
      READ  (BTH(JW),*,ERR=400) AID,(ELWS(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) 'Finished reading:', AID
      WRITE (WIN,*) '  segment orientation'
      READ  (BTH(JW),*,ERR=400) AID,(PHI0(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) 'Finished reading:', AID
      WRITE (WIN,*) '  segment bottom friction'
      READ  (BTH(JW),*,ERR=400) AID,(FRIC(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) 'Finished reading:', AID
      WRITE (WIN,*) '  layer thickness and segment widths'
      READ  (BTH(JW),*)
      DO K=1,KMX
      READ  (BTH(JW),*,ERR=400) H(K,JW),(B(K,I),I=US(BS(JW))-1,DS(BE(JW))+1)
      END DO
      else
      WRITE (WIN,*) '  segment lengths'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(DLX(I),  I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) '  water surface elevations'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(ELWS(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) '  segment orientation'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(PHI0(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) '  segment bottom friction'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(FRIC(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) '  layer thickness'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(H(K,JW), K=1,KMX)
      WRITE (WIN,*) '  segment widths'
      DO I=US(BS(JW))-1,DS(BE(JW))+1
        READ (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(B(K,I),  K=1,KMX)
      END DO
      endif
      CLOSE (BTH(JW))
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTA) 'Could not open '//BTHFN(JW)
      WRITE (WIN,FMTA) 'Could not open '//BTHFN(JW)
      STOP
    END IF
  END DO

! Names/units text files
  WRITE (WIN,*) 'Graph file: graph.npt'
  OPEN (NPT,FILE='graph.npt',STATUS='OLD',IOSTAT=IERR)
  IF (IERR == 0) THEN                                            ! SW 1/16/04 replaced section
  READ  (NPT,'(///(A43,1X,A9,3F8.0,A8))') (HNAME(J),  FMTH(J),  HMULT(J),  HYMIN(J), HYMAX(J), HPLTC(J), J=1,NHY)
  READ  (NPT,'(// (A43,1X,A9,3F8.0,A8))') (CNAME(J),  FMTC(J),  CMULT(J),  CMIN(J),  CMAX(J),  CPLTC(J), J=1,NCT)
  READ  (NPT,'(// (A43,1X,A9,3F8.0,A8))') (CDNAME(J), FMTCD(J), CDMULT(J), CDMIN(J), CDMAX(J), CDPLTC(J),J=1,NDC)
  CLOSE (NPT)
  DO JC=1,NCT
    L1         = SCAN (CNAME(JC),',')+2
    L2         = SCAN (CNAME(JC)(L1:43),'  ')+L1
    CUNIT(JC)  = CNAME(JC)(L1:L2)
    CNAME1(JC) = CNAME(JC)(1:L1-3)
    CUNIT1(JC) = CUNIT(JC)(1:1)
    CUNIT2(JC) = CUNIT(JC)
    IF (CUNIT(JC)(1:2) == 'mg') THEN
      CUNIT1(JC) = 'g'
      CUNIT2(JC) = 'g/m^3'
    END IF
    IF (CUNIT(JC)(1:2) /= 'g/' .AND. CUNIT(JC)(1:2) /= 'mg') CUNIT1(JC) = '  '
  END DO
  FMTH(1:NHY) = ADJUSTL (FMTH(1:NHY))
  ELSE
    CALL ERRORS
    WRITE (ERR,FMTA) 'Could not open graph.npt'
    WRITE (WIN,FMTA) 'Could not open graph.npt'
  END IF
  CLOSE (NPT)

  ! WSC FILE

  WRITE (WIN,*) 'Wind sheltering file:', adjustl(trim(wscfn))
  OPEN (NPT,FILE=WSCFN,STATUS='OLD',IOSTAT=IERR)
  IF (IERR == 0) THEN
  READ (NPT,'(///10F8.0:/(8X,9F8.0))') SDAY,(WSC(I),I=1,IMX)
        IF (SDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,'(3(A,F0.3))') 'Starting date ',SDAY,' > simulation start date [TMSTRT=',TMSTRT,'] in '//WSCFN
        ENDIF
        DO I=1,IMX
        IF(WSC(I) <= 0.0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,F0.3,A,I4,A)') 'Julian date ',SDAY,': WSC AT SEG(I)=',I,' <= 0.0 in '//WSCFN
        ENDIF
        IF (WSC(I) > 2.0) THEN
          CALL WARNINGS
          WRITE (WRN,'(A,F0.3,A,I4,A)') 'Julian day ',SDAY,': WSC(I) AT SEG(I)=',I,' > 2.0 in '//WSCFN
        END IF
        IF (WSC(I) > 0.0 .and. wsc(i) < 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,'(A,F0.3,A,I4,A)') 'Julian day ',SDAY,': WSC(I) AT SEG(I)=',I,' < 0.5 in '//WSCFN
        END IF
        ENDDO

  SDAYO=SDAY

!  DO J=1,100
  28995 continue  ! cb 3/26/10
  READ (NPT,'(10F8.0:/(8X,9F8.0))',END=29000) SDAY,(WSC(I),I=1,IMX)
          IF (SDAY <= SDAYO) THEN
          CALL ERRORS
          WRITE (ERR,'(3(A,F0.3))') 'Julian date ',SDAY,' <= previous date of ',SDAYO,' in '//WSCFN
          END IF
        DO I=1,IMX
        IF(WSC(I) <= 0.0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,F0.3,A,I4,A)') 'Julian date ',SDAY,': WSC AT SEG(I)=',I,' <= 0.0 in '//WSCFN
        ENDIF
        IF (WSC(I) > 2.0) THEN
          CALL WARNINGS
          WRITE (WRN,'(A,F0.3,A,I4,A)') 'Julian day ',SDAY,': WSC(I) AT SEG(I)=',I,' > 2.0 in '//WSCFN
        END IF
        IF (WSC(I) > 0.0 .and. wsc(i) < 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,'(A,F0.3,A,I4,A)') 'Julian day ',SDAY,': WSC(I) AT SEG(I)=',I,' < 0.5 in '//WSCFN
        END IF
        ENDDO
        SDAYO=SDAY
  ! ENDDO

  go to 28995   ! cb 3/26/10
  29000 IF (SDAY < TMEND) THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,F0.3))') 'Ending time ',SDAY,' < ending simulation time [TMEND=',TMEND,'] in '//WSCFN
        END IF
        IF(I < IMX) THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,F0.3))') 'Not enough records (IMX) for last day ',SDAY,' in '//WSCFN
        END IF
  CLOSE (NPT)

  ELSE
    CALL ERRORS
    WRITE (ERR,'(A)') 'Could not open:'//adjustl(trim(wscfn))
    WRITE (WIN,'(A)') 'Could not open:'//adjustl(trim(wscfn))
  END IF
  CLOSE (NPT)


!***********************************************************************************************************************************
!*                                                 Task 2:  Variable Initialization                                               **
!***********************************************************************************************************************************

  NWBR   = 0;   NTBR    = 0;   KB      = 0;   KBMAX   = 0;   JBU     = 0;   JBD     = 0;   NASS    = 0;   NAAL    = 0; KEPR    = 0
  NCCBR  = 0;   NCCGR   = 0;   CN      = 0;   INCN    = 0;   incdt   = 0;   TRCN    = 0;   NAC     = 0;   NADC    = 0;   NACIN   = 0; NACTR   = 0
  NACDT  = 0;   NACPR   = 0;   NQINMX  = 0;   NQTRMX  = 0;   NQDTMX  = 0;   NQPRMX  = 0;   NQWDMX  = 0;   NCCGR   = 0
  XI     = 0.0; XBR     = 0.0; XGR     = 0.0; T2      = 0.0; C2      = 0.0; CD      = 0.0; EL      = 0.0; VOLB    = 0.0
  VOLG   = 0.0; CVBR    = 0.0; CVGR    = 0.0; SABR    = 0.0; SAGR    = 0.0; QINMX   = 0.0; QDTMX   = 0.0; QTRMX   = 0.0
  QWDMX  = 0.0; QWDMXB  = 0.0; QPRMX   = 0.0; QSTRMX  = 0.0; QINS    = 0.0; QTRS    = 0.0; QDTS    = 0.0; QWDS    = 0.0
  QPRS   = 0.0; QSTRS   = 0.0; QINAVB  = 0.0; QINAV   = 0.0; QDTAV   = 0.0; QPRAV   = 0.0; QTRAVB  = 0.0; QOTAVB  = 0.0
  QSTRAV = 0.0; QWDAVB  = 0.0; QDNAVB  = 0.0; QTRAVB  = 0.0; QINAVW  = 0.0; QINAVB  = 0.0; QOTAVW  = 0.0; QINMXB  = 0.0
  QOTMXB = 0.0; QINMXW  = 0.0; QOTMXW  = 0.0; QTRMXB  = 0.0; CINMAX  = 0.0; CINAVG  = 0.0; SAGR    = 0.0; CVGR    = 0.0; cin=0.0;cin2=0.0
  HGR    = 0.0; BGR     = 0.0; CVGR    = 0.0; HBR     = 0.0; cinload=0.0; ctrload=0.0;cinpload=0.0;cinnload=0.0;ctrpload=0.0;ctrnload=0.0;cdtload=0.0
  CONV   = BLANK
  FMTA   = '(A)'
  FMTFI  = '(3(A,F0.3,A,I0))'
  FMTIF  = '(3(A,I0,A,F0.3))'
  FMTF   = '(3(A,F0.3))'
  FMTI   = '(4(A,I0))'
  FMT2FI = '(2(A,F0.3),A,I0)'
  FMTF2I = '(A,F0.2,2(A,I0))'
  FMT2IF = '(3(2(A,I0),A,F0.3))'

! Initialize logical control variables

  SEDIMENT_CALC    = .FALSE.; DISSOLVED_SOLIDS = .FALSE.; GENERIC_CONST    = .FALSE.; SUSPENDED_SOLIDS = .FALSE.
  PHOSPHORUS       = .FALSE.; AMMONIUM         = .FALSE.; NITRATE          = .FALSE.; IRON             = .FALSE.
  LABILE_DOM       = .FALSE.; REFRACTORY_DOM   = .FALSE.; LABILE_POM       = .FALSE.; REFRACTORY_POM   = .FALSE.
  CBO_DEMAND       = .FALSE.; ALGAE            = .FALSE.; DISSOLVED_OXYGEN = .FALSE.; TOT_IC           = .FALSE.
  OPEN_VPR         = .FALSE.; OPEN_LPR         = .FALSE.; OPEN_LPRC         = .FALSE.
  WITHDRAWALS      =  NWD   >   0;   TRIBUTARIES =  NTR   >   0;   PIPES     =  NPI   >   0
  GATES            =  NGT   >   0;   SPILLWAYS   =  NSP   >   0;   PUMPS     =  NPU   >   0
  ISO_TEMP         =  T2I   >=  0.0; VERT_TEMP   =  T2I   == -1.0; LONG_TEMP =  T2I   == -2.0
  LONG_SED         =  SEDCI == -2.0
  OPEN_VPR         =  VERT_TEMP
  OPEN_LPR         =  LONG_TEMP
  PRECIPITATION    =  PRC  == '      ON'
  CONSTITUENTS     =  CCC  == '      ON'
  ICE_CALC         =  ICEC == '      ON'
  RESTART_IN       =  RSIC == '      ON'
  DO JW=1,NWB
    IF (LONG_SED(JW)) OPEN_LPRC(JW) = .TRUE.
    DO JC=1,NCT
      ISO_CONC(JC,JW)  = C2IWB(JC,JW) >=  0.0
      VERT_CONC(JC,JW) = C2IWB(JC,JW) == -1.0  .AND.  CAC(JC) == '      ON'
      LONG_CONC(JC,JW) = C2IWB(JC,JW) <  -1.0  .AND.  CAC(JC) == '      ON'
      IF (VERT_CONC(JC,JW)) OPEN_VPR(JW) = .TRUE.
      IF (LONG_CONC(JC,JW)) OPEN_LPRC(JW) = .TRUE.
    END DO
  END DO
  DO JB=1,NBR
    DAM_FLOW(JB)    = UHS(JB) <  -1; UH_EXTERNAL(JB) = UHS(JB) == -1; DH_EXTERNAL(JB) = DHS(JB) == -1
    UH_INTERNAL(JB) = UHS(JB) >   0; DH_INTERNAL(JB) = DHS(JB) >   0; UQ_EXTERNAL(JB) = UHS(JB) ==  0
    DQ_EXTERNAL(JB) = DHS(JB) ==  0; DQ_INTERNAL(JB) = DQB(JB) >   0
    UQ_INTERNAL(JB) = UQB(JB) >   0 .AND. .NOT. DAM_FLOW(JB)
    DIST_TRIB(JB)   = DTRC(JB) == '      ON'
    UHS(JB)         = ABS(UHS(JB))
    DHS(JB)         = ABS(DHS(JB))
  END DO

! Inflow/outflow branch locations

  IF (WITHDRAWALS) THEN
    JBWD = 0
    NWBR = 0
    DO JB=1,NBR
      DO JW=1,NWD
        IF (IWD(JW) >= US(JB) .AND. IWD(JW) <= DS(JB)) THEN
          JBWD(JW) = JB
          NWBR(JB) = NWBR(JB)+1
        END IF
      END DO
    END DO
  END IF
  IF (TRIBUTARIES) THEN
    JBTR = 0
    NTBR = 0
    DO JB=1,NBR
      DO JT=1,NTR
        IF (ITR(JT) >= US(JB) .AND. ITR(JT) <= DS(JB)) THEN
          JBTR(JT) = JB
          NTBR(JB) = NTBR(JB)+1
        END IF
      END DO
    END DO
  END IF
  IF (PIPES) THEN
    JBUPI  = 0
    NUPIBR = 0
    JBDPI  = 0
    NDPIBR = 0
    DO JB=1,NBR
      DO JP=1,NPI
        IF (IUPI(JP) >= US(JB) .AND. IUPI(JP) <= DS(JB)) THEN
          JBUPI(JP)  = JB
          NUPIBR(JB) = NUPIBR(JB)+1
        END IF
        IF (IDPI(JP) >= US(JB) .AND. IDPI(JP) <= DS(JB)) THEN
          JBDPI(JP)  = JB
          NDPIBR(JB) = NDPIBR(JB)+1
        END IF
      END DO
    END DO
  END IF
  IF (GATES) THEN
    JBUGT  = 0
    NUGTBR = 0
    JBDGT  = 0
    NDGTBR = 0
    DO JB=1,NBR
      DO JG=1,NGT
        IF (IUGT(JG) >= US(JB) .AND. IUGT(JG) <= DS(JB)) THEN
          JBUGT(JG)  = JB
          NUGTBR(JB) = NUGTBR(JB)+1
        END IF
        IF (IDGT(JG) >= US(JB) .AND. IDGT(JG) <= DS(JB)) THEN
          JBDGT(JG)  = JB
          NDGTBR(JB) = NDGTBR(JB)+1
        END IF
      END DO
    END DO
  END IF

! Active constituents

  IF (CAC(NTDS)  == '      ON') DISSOLVED_SOLIDS = .TRUE.
  IF (CAC(NPO4)  == '      ON') PHOSPHORUS       = .TRUE.
  IF (CAC(NNH4)  == '      ON') AMMONIUM         = .TRUE.
  IF (CAC(NNO3)  == '      ON') NITRATE          = .TRUE.
  IF (CAC(NFE)   == '      ON') IRON             = .TRUE.
  IF (CAC(NLDOM) == '      ON') LABILE_DOM       = .TRUE.
  IF (CAC(NRDOM) == '      ON') REFRACTORY_DOM   = .TRUE.
  IF (CAC(NLPOM) == '      ON') LABILE_POM       = .TRUE.
  IF (CAC(NRPOM) == '      ON') REFRACTORY_POM   = .TRUE.
  IF (CAC(NDO)   == '      ON') DISSOLVED_OXYGEN = .TRUE.
  IF (CAC(NTIC)  == '      ON') TOT_IC           = .TRUE.
  DO JG=NGCS,NGCE
    IF (CAC(JG)  == '      ON') THEN
      GENERIC_CONST = .TRUE.
      EXIT
    END IF
  END DO
  DO JS=NSSS,NSSE
    IF (CAC(JS) == '      ON') THEN
      SUSPENDED_SOLIDS = .TRUE.
      EXIT
    END IF
  END DO
  DO JB=NBODS,NBODE
    IF (CAC(JB) == '      ON') THEN
      CBO_DEMAND = .TRUE.
      EXIT
    END IF
  END DO
  DO JA=NAS,NAE
    IF (CAC(JA) == '      ON') THEN
      ALGAE = .TRUE.
      EXIT
    END IF
  END DO
  DO JW=1,NWB
    IF (SEDC(JW) == '      ON') SEDIMENT_CALC(JW) = .TRUE.
  END DO
  IF (CONSTITUENTS) THEN
    DO JC=1,NCT
      IF (CAC(JC) == '      ON') THEN
        NAC     = NAC+1
        CN(NAC) = JC
      END IF
      DO JB=1,NBR
        IF (CINBRC(JC,JB) == '      ON') THEN                                                                          !SW 01/07/01
          NACIN(JB)          = NACIN(JB)+1                                                                             !SW 01/07/01
          INCN(NACIN(JB),JB) = JC                                                                                      !SW 01/07/01
        END IF
        IF (CDTBRC(JC,JB) == '      ON')then
            NACDT(JB) = NACDT(JB)+1
            incdt(nacdt(jb),jb)=jc
        ENDIF
      END DO
      DO JT=1,NTR                                                                                                      !SW 01/07/01
        IF (CINTRC(JC,JT) == '      ON') THEN                                                                          !SW 01/07/01
          NACTR(JT)          = NACTR(JT)+1                                                                             !SW 01/07/01
          TRCN(NACTR(JT),JT) = JC                                                                                      !SW 01/07/01
        END IF
      END DO                                                                                                           !SW 01/07/01
    END DO
    DO JW=1,NWB
      DO JC=1,NCT
        IF (CPRBRC(JC,JW) == '      ON') NACPR(JW) = NACPR(JW)+1
      END DO
      DO JD=1,NDC
        IF (CDWBC(JD,JW) == '      ON') THEN
          NADC(JW)         = NADC(JW)+1
          CDN(NADC(JW),JW) = JD                                                                                        !TC 09/23/02
        END IF
      END DO
    END DO
    DO JS=1,NSS
      IF (CAC(NSSS) == '      ON') THEN
        NASS       = NASS+1
        SSAC(NASS) = JS
      END IF
    END DO
    DO JA=1,NAL
      IF (CAC(NAS) == '      ON') THEN
        NAAL       = NAAL+1
        ALAC(NAAL) = JA
      END IF
    END DO
  END IF





! Grid linkage

  DO JB=1,NBR

  ! Check if UHS is negative

  IF(DAM_FLOW(JB))THEN
    ! Find branch that this is connected to
    DO JJB=1,NBR
        IF (UHS(JB) == DS(JJB))THEN
          IF(NSTR(JJB) == 0)THEN
           CALL ERRORS
           WRITE (ERR,FMTI) 'Error grid linkage:Since UHS(JB) is < 0, # STR must be > 0 for connected branch(JJB). NSTR(JJB)=0. JB:',JB,' JJB:',JJB
          ENDIF
        EXIT
        ENDIF
        IF(JJB == NBR)THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Error grid linkage:Since UHS is < 0, this MUST be the DS of another branch. JB:',JB
        ENDIF
    END DO

  END IF

    JCHECK = 0
    IF (UHS(JB) /= 0 .AND. .NOT. UH_EXTERNAL(JB) .AND. .NOT. DH_EXTERNAL(JB)) THEN                                     !TC 08/31/01
      DO JJB=1,NBR
        IF (ABS(UHS(JB)) >= US(JJB) .AND. ABS(UHS(JB)) <= DS(JJB)) EXIT
      END DO
      if(jjb.gt.nbr)then                                                                                               ! SW 6/1/07
            CALL ERRORS
            WRITE (ERR,FMTI) 'Error in grid linkage:check[US],[DS],[UHS],[DHS] and/or[JBDN]for waterbody ',JW,' JB=',JB
            WRITE (WIN,FMTI) 'Error in grid linkage:check[US],[DS],[UHS],[DHS] and/or[JBDN]for waterbody ',JW,' JB=',JB
            STOP
      endif
      IF (ABS(UHS(JB)) == DS(JJB)) THEN
        IF (DHS(JJB) == US(JB)) JCHECK = 1
        IF (UHS(JB) < 0) THEN
          IF (JCHECK == 1) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Error in grid linkage:check[US],[DS],[UHS],[DHS] and/or[JBDN]for waterbody ',JW,' JB=',JB
            WRITE (WIN,FMTI) 'Error in grid linkage:check[US],[DS],[UHS],[DHS] and/or[JBDN]for waterbody ',JW,' JB=',JB
            STOP
          END IF
        END IF
      END IF
    END IF
  END DO
  ZERO_SLOPE = .TRUE.
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)                                                                                                !SW 06/12/01
      IF (SLOPE(JB) /= 0.0) ZERO_SLOPE(JW) = .FALSE.                                                                   !SW 10/17/02
    END DO                                                                                                             !SW 06/12/01
  END DO                                                                                                               !SW 06/12/01

!***********************************************************************************************************************************
!*                                                Task 3: Input/Output File Checks                                                **
!***********************************************************************************************************************************

! Layer elevations

  DO JW=1,NWB                                                                                                          !SW 10/17/02
    IF (ZERO_SLOPE(JW)) THEN
      DO I=US(BS(JW))-1,DS(BE(JW))+1
        EL(KMX,I) = ELBOT(JW)
        DO K=KMX-1,1,-1
          EL(K,I) = EL(K+1,I)+H(K,JW)
        END DO
      END DO
    ELSE
      NPOINT = 0                                                                                                       !TC 02/11/01
      NUP              = 0
      NNBR             = 1
      NCBP             = 0
      NINTERNAL        = 0
      JB               = JBDN(JW)
      NPOINT(JB)       = 1
      EL(KMX,DS(JB)+1) = ELBOT(JW)
      DO WHILE (NNBR <= (BE(JW)-BS(JW)+1))
        NCBP = NCBP+1
        IF (NCBP > NBR) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Model linkage not logical for waterbody ',JW
          EXIT
        END IF
        IF (NINTERNAL == 0) THEN
          IF (NUP == 0) THEN
            DO I=DS(JB),US(JB)-1,-1
              IF (I /= DS(JB)) THEN
                EL(KMX,I) = EL(KMX,I+1)+SLOPE(JB)*(DLX(I)+DLX(I+1))*0.5
              ELSE
                EL(KMX,I) = EL(KMX,I+1)                                                                                !bug?
              END IF
              DO K=KMX-1,1,-1
                EL(K,I) = EL(K+1,I)+H(K,JW)*COS(ATAN(SLOPE(JB)))
              END DO
            END DO
          ELSE
            DO I=US(JB),DS(JB)+1
              IF (I /= US(JB)) THEN
                EL(KMX,I) = EL(KMX,I-1)-SLOPE(JB)*(DLX(I)+DLX(I+1))*0.5
              ELSE
                EL(KMX,I) = EL(KMX,I-1)
              END IF
              DO K=KMX-1,1,-1
                EL(K,I) = EL(K+1,I)+H(K,JW)*COS(ATAN(SLOPE(JB)))
              END DO
            END DO
            NUP = 0
          END IF
          DO K=KMX,1,-1
            IF (.NOT.DAM_FLOW(JB)) THEN
              IF (UHS(JB) /= 0) THEN
                EL(K,US(JB)-1) = EL(K,US(JB))+SLOPE(JB)*DLX(US(JB))
              ELSE
                EL(K,US(JB)-1) = EL(K,US(JB))
              END IF
            ELSE
              EL(K,US(JB)-1) = EL(K,US(JB))
            END IF
            IF (DHS(JB) /= 0) THEN
              EL(K,DS(JB)+1) = EL(K,DS(JB))-SLOPE(JB)*DLX(DS(JB))
            ELSE
              EL(K,DS(JB)+1) = EL(K,DS(JB))
            END IF
          END DO
        ELSE
          DO K=KMX-1,1,-1
            EL(K,UHS(JJB)) = EL(K+1,UHS(JJB))+H(K,JW)*COS(ATAN(SLOPE(JB)))
          END DO
          DO I = UHS(JJB)+1, DS(JB)
            EL(KMX,I) = EL(KMX,I-1)-SLOPE(JB)*(DLX(I)+DLX(I-1))*0.5
            DO K=KMX-1,1,-1
              EL(K,I) = EL(K+1,I)+H(K,JW)*COS(ATAN(SLOPE(JB)))
            END DO
          END DO
          DO I=UHS(JJB)-1,US(JB),-1
            EL(KMX,I) = EL(KMX,I+1)+SLOPE(JB)*(DLX(I)+DLX(I+1))*0.5
            DO K=KMX-1,1,-1
              EL(K,I) = EL(K+1,I)+H(K,JW)*COS(ATAN(SLOPE(JB)))
            END DO
          END DO
          NINTERNAL = 0
        END IF
        IF (NNBR == (BE(JW)-BS(JW)+1)) EXIT
        NNBROLD = NNBR
        DO JB=BS(JW),BE(JW)
          IF (NPOINT(JB) /= 1) THEN
            DO JJB = BS(JW), BE(JW)
              IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1
                NNBR             = NNBR+1
                EL(KMX,DS(JB)+1) = EL(KMX,DHS(JB))+SLOPE(JB)*(DLX(DS(JB))+DLX(DHS(JB)))*0.5
                EXIT
              END IF
              IF (UHS(JJB) == DS(JB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1
                NNBR             = NNBR+1
                EL(KMX,DS(JB)+1) = EL(KMX,US(JJB))+SLOPE(JB)*(DLX(US(JJB))+DLX(DS(JB)))*0.5
                EXIT
              END IF
              IF (UHS(JJB) <= DS(JB) .AND. UHS(JJB) >= US(JB) .AND. NPOINT(JJB) == 1) THEN
                NNBR             = NNBR+1
                NINTERNAL        = 1
                NPOINT(JB)       = 1
                EL(KMX,UHS(JJB)) = EL(KMX,US(JJB))+SLOPE(JJB)*DLX(US(JJB))*0.5
                EXIT
              END IF
              IF (UHS(JB) <= DS(JJB) .AND. UHS(JB) >= US(JJB) .AND. NPOINT(JJB) == 1) THEN
                NNBR             = NNBR+1
                NPOINT(JB)       = 1
                EL(KMX,US(JB)-1) = EL(KMX,UHS(JB))-SLOPE(JB)*DLX(US(JB))*0.5
                IF (SLOPE(JB) /= 0.0) THEN
                  CALL WARNINGS
                  WRITE (WRN,FMTI) 'Branch ',JB,' has non-zero slope and is only'
                  WRITE (WRN,FMTI) 'connected through a UHS BC. This may not be'
                  WRITE (WRN,FMTI) 'a problem:if it is,change branch to a waterbody'                                   !TC 02/09/01
                END IF
                NUP = 1
                EXIT
              END IF
            END DO
            IF (NPOINT(JB) == 1) EXIT
          END IF
        END DO
        IF (NCBP > NBR .OR. NNBROLD == NNBR) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Error in grid linkage: check[US],[DS],[UHS],[DHS],[JBDN];each branch must be connected to another branch for waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! NPROC and closec

IF(NPROC < 0)THEN
 CALL ERRORS
 WRITE(ERR,*)'ERROR: NPROC < 0, MUST BE > 0'
ELSEIF(NPROC==0)THEN
 CALL ERRORS
 WRITE(ERR,*)'NPROC = 0, MUST BE > 0 [W2 CODE WILL SET NPROC=2 WHEN NPROC=0]'
ELSEIF(NPROC > 4)THEN
 CALL WARNINGS
 WRITE(WRN,*)'NPROC > 4. MAKE SURE THIS DOES NOT LEAD TO LONG RUN TIMES. TRY 2.'
ENDIF
IF(CLOSEC /= '     OFF' .AND. CLOSEC /= '      ON')THEN
 CALL ERRORS
 WRITE(ERR,*)'CLOSEC MUST BE "OFF" OR "ON". CLOSEC=',CLOSEC
ENDIF

! Grid geometry

  WRITE (WIN,*) 'Bathymetry file checks'
  WRITE (WIN,*) '  water surface elevations'
  DO JW=1,NWB
    ZMIN = -1000.0
    DO JB=BS(JW),BE(JW)
      ALPHA(JB) = ATAN(SLOPE(JB))
      IFLAG     = 0
      DO I=US(JB)-1,DS(JB)+1
        KTI(I) = 2
        DO WHILE (EL(KTI(I),I) > ELWS(I) .AND. KTI(I) < KMX)
          KTI(I) = KTI(I)+1
          IF (KTI(I) > KMX-1) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Initial water surface elevation below grid at segment ',I
            STOP
          END IF
        END DO
        Z(I)     = (EL(KTI(I),I)-ELWS(I))/COS(ALPHA(JB))
        ZMIN     =  MAX(ZMIN,Z(I))
        KTMAX    =  MAX(2,KTI(I))
        KTWB(JW) =  MAX(KTMAX,KTWB(JW))
        KTI(I)   =  MAX(KTI(I)-1,2)
        IF (KTWB(JW) >= KMX .AND. IFLAG == 0) THEN
          IFLAG =  1
          CALL ERRORS
          WRITE (ERR,FMTI) 'Water surface elevation too low in branch ',JB,' at segment ',I
        END IF
        DO K=2,KMX-1
          IF (B(K+1,I) == 0.0 .OR. K == KMX-1) THEN
            KB(I) = K
            EXIT
          END IF
        END DO
        KBMAX(JW) = MAX(KBMAX(JW),KB(I))
      END DO
      KB(US(JB)-1) = KB(US(JB))
      KB(DS(JB)+1) = KB(DS(JB))
    END DO
  END DO

! Boundary cells

  WRITE (WIN,*) '  widths'
  DO JB=1,NBR
    DO K=1,KMX
      IF (B(K,US(JB)-1) /= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMT2IF) 'Upstream boundary segment width [B(',K,',',US(JB)-1,')=',B(K,US(JB)-1),'] /= 0 for branch ',JB
      END IF
      IF (B(K,DS(JB)+1) /= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMT2IF) 'Downstream boundary segment width [B(',K,',',DS(JB)+1,')=',B(K,DS(JB)+1),'] /= 0 for branch ',JB
      END IF
    END DO
    DO I=US(JB),DS(JB)
      DO K=2,KB(I)
        IF (B(K+1,I) > B(K,I) .AND. B(K+1,I) /= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMT2IF) 'Cell width [B(',K+1,',',I,')=',B(K+1,I),'] > [B(',K,',',I,')=',B(K,I)
        END IF
        IF (B(K,I) < 5.0 .AND. B(K,I) > 0.0) THEN                                                                      !SW 01/23/01
          CALL WARNINGS
          WRITE (WRN,FMT2IF) 'Cell width [B(',K,',',I,')=',B(K,I),'] < 5m which can cause stability problems'
        END IF
      END DO
      IF (B(1,I) /= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Surface boundary layer cell width [B(1,',I,')=',B(1,I),'] /= 0'
      ELSE IF (B(KMX,I) /= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMT2IF) 'Bottom boundary layer cell width [B(',KMX,',',I,')=',B(KMX,I),'] /= 0'
      END IF
    END DO
  END DO

! Boundary widths

  DO JB=1,NBR
    DO K=2,KMX
      B(K,US(JB)-1) = B(K,US(JB))
      B(K,DS(JB)+1) = B(K,DS(JB))
    END DO
  END DO

! Width change

  DO JB=1,NBR
    DO I=US(JB),DS(JB)-1
      DO K=2,KB(I)
        IF (B(K,I+1) /= 0.0 .AND. B(K,I) /= 0.0) THEN
          IF (B(K,I)/B(K,I+1) >= 10.0 .OR. B(K,I+1)/B(K,I) >= 10.0) THEN
            CALL WARNINGS
            WRITE (WRN,FMT2IF) 'Cell width [B(',K,',',I,')=',B(K,I),'] < or > 10X width [B(',K,',',I+1,')=',B(K,I+1),']'

          END IF
        END IF
      END DO
    END DO
  END DO

! Layer thickness

  WRITE (WIN,*) '  layer thickness'
  DO JW=1,NWB
    DO K=1,KMX
      IF (H(K,JW) <= 0) THEN
        CALL   ERRORS
        WRITE (ERR,FMTF2I) 'Layer thickness [H =',H(K,JW),'] <= 0.0 m for layer ',K,' in waterbody ',JW
      ELSE IF (H(K,JW) < 0.1) THEN
        CALL   WARNINGS
        WRITE (WRN,FMTF2I) 'Layer thickness [H =',H(K,JW),'] <  0.1 m for layer ',K,' in waterbody ',JW
      ELSE IF (H(K,JW) > 3.0) THEN
        CALL   WARNINGS
        WRITE (WRN,FMTF2I) 'Layer thickness [H =',H(K,JW),'] >  3.0 m for layer ',K,' in waterbody ',JW
      END IF
    END DO
  END DO

! Meteorologic data

  WRITE (WIN,*) 'Input files'
  WRITE (WIN,*) '  meteorology'
  DO JW=1,NWB
    write(win,*) '    ',metfn(jw)
    OPEN (UNIT=NPT,FILE=METFN(JW),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      READ (NPT,'(//)')
      J     = 1
      JDAYO = 0.0
      DO WHILE (.TRUE.)
        READ (NPT,'(10F8.0/8X,9F8.0)',END=110) JDAY, TAIR, TDEW, WIND, PHI, CLOUD
        IF (J == 1 .AND. JDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Starting date ',JDAY,' > simulation start date [TMSTRT=',TMSTRT,'] in '//METFN(JW)
        ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date of ',JDAYO,' in '//METFN(JW)
        END IF
        IF (TAIR > 50.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF) 'Air temperature [TAIR=',TAIR,'] > 50 deg on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (TDEW > 50.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF) 'Dew point temperature [TDEW=',TDEW,'] > 50 deg on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (WIND > 20.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF) 'Wind speed [WIND= ',WIND,'] > 20 m/s on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (WIND < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Wind speed [WIND=',WIND,'] < 0 m/s on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (CLOUD > 10.0 .OR. CLOUD < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Cloud cover [CLOUD=',CLOUD,'] < 0 or > 10 on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (TAIR == 0.0 .AND. TDEW == 0.0 .AND. ABS(TAIR-TAIRO) > 5.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF) 'Possible missing value on day ',JDAY,' in '//METFN(JW)
        END IF
        J     = J+1
        JDAYO = JDAY
        TAIRO = TAIR
      END DO
 110  CONTINUE
      IF (JDAY < TMEND) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//METFN(JW)
      END IF
      CLOSE (NPT)
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTI) 'Could not open waterbody ',JW,' meteorological file '//METFN(JW)
    END IF
  END DO

! Number of records for inflow files

  DO JB=1,NBR
    OPEN (UNIT=NPT,FILE=QINFN(JB),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      DO J=1,1000000
        READ (NPT,*,END=130)
      END DO
130   CONTINUE
      NQIN(JB) = J-4
      NQINMX   = MAX(NQINMX,NQIN(JB))+5
      CLOSE (NPT)
    END IF
    OPEN (UNIT=NPT,FILE=QDTFN(JB),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      DO J=1,1000000
        READ (NPT,*,END=131)
      END DO
131   CONTINUE
      NQDT(JB) = J-4
      NQDTMX   = MAX(NQDTMX,NQDT(JB))+5
      CLOSE (NPT)
    END IF
    OPEN (UNIT=NPT,FILE=PREFN(JB),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      DO J=1,1000000
        READ (NPT,*,END=132)
      END DO
132   CONTINUE
      NQPR(JB) = J-4
      NQPRMX   = MAX(NQPRMX,NQPR(JB))+5
      CLOSE (NPT)
    END IF
  END DO
  DO JT=1,NTR
    OPEN (UNIT=NPT,FILE=QTRFN(JT),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      DO J=1,1000000
        READ (NPT,*,END=133)
      END DO
133   CONTINUE
      NQTR(JT) = J-3
      NQTRMX   = MAX(NQTRMX,NQTR(JT))+5
      CLOSE (NPT)
    END IF
  END DO
  OPEN (UNIT=NPT,FILE=QWDFN,STATUS='OLD',IOSTAT=IERR)
  IF (IERR == 0) THEN
    DO J=1,10000000
      READ (NPT,*,END=134)
    END DO
134 CONTINUE
    NQWD    = J-4
    NQWDMX  = MAX(NQWDMX,NQWD)+5
    CLOSE (NPT)
  END IF
  ALLOCATE (JDQIN(NQINMX,NBR),QIN(NQINMX,NBR))
  ALLOCATE (JDQDT(NQDTMX,NBR),QDT(NQDTMX,NBR))
  ALLOCATE (JDQPR(NQPRMX,NBR),QPR(NQPRMX,NBR))
  ALLOCATE (JDQTR(NQTRMX,NTR),QTR(NQTRMX,NTR))
  ALLOCATE (JDQWD(NQWDMX,NWD),QWD(NQWDMX,NWD))
  QIN = 0.0; QDT = 0.0; QPR = 0.0; QTR = 0.0; JDQIN = 0.0; JDQDT = 0.0; JDQPR = 0.0; JDQTR = 0.0

! Time varying data

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      WRITE (BRA,'(I0)') JB
      IF (UQ_EXTERNAL(JB)) THEN

!****** Inflows

        WRITE (WIN,*) '  inflow'
        OPEN (UNIT=NPT,FILE=QINFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
          READ  (NPT,'(A)') HEADER
          J       = 1
          JDAYO   = 0.0
          QINO    = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0/8X,9F8.0)',END=135) JDQIN(J,JB), QIN(J,JB)

            IF (QIN(J,JB)<0.0) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'QIN cannot be negative:',QIN(J,JB),' at date [JDAY=',JDQIN(J,JB),'] in '//QINFN(JB)
            END IF

            IF (JDAYO >= TMSTRT .AND. JDAYO <= TMEND) THEN
              QINMX(JB) = MAX(QIN(J,JB),QINMX(JB))
              IF (JDQIN(J,JB) > TMEND) THEN
                QINS(JB) = QINS(JB)+QINO*(TMEND-JDAYO)
              ELSE
                QINS(JB) = QINS(JB)+QINO*(JDQIN(J,JB)-JDAYO)
              END IF
            END IF
            IF ((J == 1) .AND. (JDQIN(J,JB) > TMSTRT)) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDQIN(J,JB),' > simulation start date [TMSTRT=',TMSTRT,'] in '//QINFN(JB)
            ELSE IF (JDQIN(J,JB) <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDQIN(J,JB),' <= previous date ',JDAYO,' in '//QINFN(JB)
            END IF
            QINO  = QIN(J,JB)
            JDAYO = JDQIN(J,JB)
            J     = J+1
          END DO
 135      CONTINUE
          IF (JDQIN(J-1,JB) < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDQIN(J-1,JB),' < ending simulation time [TMEND=',TMEND,'] in '//QINFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' inflow file '//QINFN(JB)
        END IF

!****** Inflow temperatures

        WRITE (WIN,*) '  inflow temperature'
        OPEN (UNIT=NPT,FILE=TINFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
          READ  (NPT,'(A)') HEADER
          J       = 1
          JDAYO   = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0/8X,9F8.0)',END=140) JDAY, TIN
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in ' //TINFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TINFN(JB)
            END IF
            JDAYO = JDAY
            J     = J+1
          END DO
 140      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TINFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' inflow temperature file '//TINFN(JB)
        END IF

!****** Inflow concentrations

        IF (CONSTITUENTS .AND. NACIN(JB) > 0) THEN
          WRITE (WIN,*) '  inflow concentrations'
          OPEN (UNIT=NPT,FILE=CINFN(JB),IOSTAT=IERR,STATUS='OLD')
          IF (IERR == 0) THEN
            READ (NPT,'(//)')
            J                               = 1
            JDAYO                           = 0.0
            CINMIN(INCN(1:NACIN(JB),JB),JB) = 1.0E10
            DO WHILE (.TRUE.)
              READ (NPT,'(F8.0,1000(F8.0))',END=150) JDAY, (CIN(INCN(JC,JB)),JC=1,NACIN(JB))
              IF (J == 1) TSTART = JDAY
              CINMIN(INCN(1:NACIN(JB),JB),JB) = MIN(CIN(INCN(1:NACIN(JB),JB)),CINMIN(INCN(1:NACIN(JB),JB),JB))
              CINMAX(INCN(1:NACIN(JB),JB),JB) = MAX(CIN(INCN(1:NACIN(JB),JB)),CINMAX(INCN(1:NACIN(JB),JB),JB))
              IF (J > 1) CINAVG(INCN(1:NACIN(JB),JB),JB) =  CINAVG(INCN(1:NACIN(JB),JB),JB)+(CIN(INCN(1:NACIN(JB),JB))             &
                                                            +CIN2(INCN(1:NACIN(JB),JB)))*0.5*(JDAY-TDAY2)

! compute mass loading from inflows
             if(j==1)then
             open(npt+1,file=qinfn(jb),status='old')
             READ  (NPT+1,'(A)') HEADER
             QINO    = 0.0
             jj=1
             READ (NPT+1,'(10F8.0/8X,9F8.0)') JDQIN(JJ,JB), QIN(JJ,JB)
             endif
                do while(jdqin(jj,jb) < jday)
                jj=jj+1
                READ (NPT+1,'(10F8.0/8X,9F8.0)',end=151) JDQIN(JJ,JB), QIN(JJ,JB)
                enddo
                if(jday == jdqin(jj,jb).or.jj==1)then
                qq=qin(jj,jb)
                else
                qq=qin(jj,jb)+((qin(jj,jb)-qin(jj-1,jb))/(jdqin(jj,jb)-jdqin(jj-1,jb)))*(jdqin(jj,jb)-jday)
                endif
                if(j>1)cinload(INCN(1:NACIN(JB),JB),jb)=cinload(INCN(1:NACIN(JB),JB),jb)+qq*(CIN(INCN(1:NACIN(JB),JB))+CIN2(INCN(1:NACIN(JB),JB)))*0.5*(JDAY-TDAY2)*86400./1000.  ! units of kg
                ! end of inflow loading

151           CIN2(INCN(1:NACIN(JB),JB)) = CIN(INCN(1:NACIN(JB),JB))
              TDAY2                      = JDAY
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CINFN(JB)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CINFN(JB)
              END IF
              JDAYO = JDAY
              J     = J+1
           END DO
 150        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CINFN(JB)
            END IF
            CLOSE (NPT)
            close(npt+1)
            cinload(INCN(1:NACIN(JB),JB),jb)=cinload(INCN(1:NACIN(JB),JB),jb)/(jday-tstart)
            CINAVG(INCN(1:NACIN(JB),JB),JB) = CINAVG(INCN(1:NACIN(JB),JB),JB)/(JDAY-TSTART)
          ELSE
            CALL WARNINGS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' inflow concentration file '//CINFN(JB)
          END IF
        END IF
      END IF

!**** Distributed tributaries

      IF (DIST_TRIB(JB)) THEN

!****** Inflows

        WRITE (WIN,*) '  distributed tributary inflow'
        OPEN (UNIT=NPT,FILE=QDTFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
          READ (NPT,'(//)')
          J     = 1
          QDTO  = 0.0
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0/8X,9F8.0)',END=160) JDQDT(J,JB), QDT(J,JB)
            IF (JDAYO >= TMSTRT .AND. JDQDT(J,JB) <= TMEND) THEN
              QDTMX(JB) = MAX(QDT(J,JB),QDTMX(JB))
              IF (JDQDT(J,JB) > TMEND) THEN
                QDTS(JB) = QDTS(JB)+QDTO*(TMEND-JDAYO)
              ELSE
                QDTS(JB) = QDTS(JB)+QDTO*(JDQDT(J,JB)-JDAYO)
              END IF
            END IF
            IF (J == 1 .AND. JDQDT(J,JB) > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDQDT(J,JB),'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QDTFN(JB)
            ELSE IF (JDQDT(J,JB) <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDQDT(J,JB),' <= previous date ',JDAYO,' in '//QDTFN(JB)
            END IF
            JDAYO = JDQDT(J,JB)
            QDTO  = QDT(J,JB)
            J     = J+1
          END DO
 160      CONTINUE
          IF (JDQDT(J-1,JB) < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDQDT(J,JB),' < ending simulation time [TMEND=',TMEND,'] in '//QDTFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' distributed tributary inflow file '//QDTFN(JB)
        END IF

!****** Inflow temperatures

        WRITE (WIN,*) '  distributed tributary inflow temperature'
        OPEN (UNIT=NPT,FILE=TDTFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
          READ (NPT,'(//)')
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0/8X,9F8.0)',END=170) JDAY, TDT
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TDTFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TDTFN(JB)
            END IF
            J     = J+1
            JDAYO = JDAY
          END DO
 170      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TDTFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' distributed tributary inflow temperature file '//TDTFN(JB)    !TC 09/30/03
        END IF

!****** Inflow concentrations

        IF (CONSTITUENTS .AND. NACDT(JB) > 0) THEN
        cin=0.0
        cin2=0.0
          WRITE (WIN,*) '  distributed tributary concentrations'
          OPEN (UNIT=NPT,FILE=CDTFN(JB),IOSTAT=IERR,STATUS='OLD')
          IF (IERR == 0) THEN
            READ (NPT,'(//)')
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ (NPT,'(F8.0,100(F8.0))',END=180) JDAY, (CIN(incdt(jc,jb)),JC=1,NACDT(JB))
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CDTFN(JB)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CDTFN(JB)
              END IF

                 ! compute mass loading from dist trib inflows

             if(j==1)then
             open(npt+1,FILE=QDTFN(JB),STATUS='OLD')
             READ  (NPT+1,'(A)') HEADER
             QINO    = 0.0
             jj=1
             READ (NPT+1,'(10F8.0/8X,9F8.0)')  JDQDT(J,JB), QDT(J,JB)
             endif
                do while(JDQDT(jj,JB) < jday)
                jj=jj+1
                READ (NPT+1,'(10F8.0/8X,9F8.0)',end=180) JDQDT(JJ,JB), QDT(JJ,JB)
                enddo
                if(jday == jdqdt(jj,jb).or.jj==1)then
                qq=qdt(jj,jb)
                else
                qq=qdt(jj,jb)+((qdt(jj,jb)-qdt(jj-1,jb))/(jdqdt(jj,jb)-jdqdt(jj-1,jb)))*(jdqdt(jj,jb)-jday)
                endif
                if(qq<0.0)qq=0.0
                if(j>1)cdtload(incdt(1:NACdt(Jb),Jb),jb)=cdtload(incdt(1:NACdt(Jb),Jb),jb)+qq*(Cin(incdt(1:NACdt(Jb),Jb))+CIN2(incdt(1:NACdt(Jb),Jb)))*0.5*(JDAY-TDAY2)*86400./1000.  ! units of kg
                ! end of inflow loading

              JDAYO = JDAY
              CIN2(incdt(1:NACdt(Jb),Jb)) = CIN(incdt(1:NACdt(Jb),Jb))
              TDAY2=JDAY
              J     = J+1
            END DO
 180        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CDTFN(JB)
            END IF
            CLOSE (NPT)
            close(npt+1)
            cdtload(incdt(1:NACdt(Jb),Jb),jb)=cdtload(incdt(1:NACdt(Jb),Jb),jb)/(jday-tstart)                        ! CB 5/10/10 Change units to kg/day
          ELSE
            CALL WARNINGS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' distributed tributary inflow concentration file '//CDTFN(JB)
          END IF
        END IF
      END IF

!**** Precipitation

      IF (PRECIPITATION(JW)) THEN

!****** Inflows

        WRITE (WIN,*) '  precipitation'
        OPEN (UNIT=NPT,FILE=PREFN(JB),STATUS='OLD')
        IF (IERR == 0) THEN
          READ (NPT,'(//)')
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0/8X,9F8.0)',END=190) JDQPR(J,JB), QPR(J,JB)
            IF (JDAYO >= TMSTRT .AND. JDQPR(J,JB) <= TMEND) THEN
              QPRMX(JB) = MAX(QPR(J,JB),QPRMX(JB))
              IF (JDQPR(J,JB) > TMEND) THEN
                QPRS(JB) = QPRS(JB)+QPR(J,JB)*(TMEND-JDAYO)
              ELSE
              QPRS(JB) = QPRS(JB)+QPR(J,JB)*(JDQPR(J,JB)-JDAYO)
              END IF
            END IF
            IF (J == 1 .AND. JDQPR(J,JB) > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDQPR(J,JB),'] > simulation start date [TMSTRT=',TMSTRT,'] in '//PREFN(JB)
            ELSE IF (JDQPR(J,JB) <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDQPR(J,JB),' <= previous date ',JDAYO,' in '//PREFN(JB)
            END IF
            JDAYO = JDQPR(J,JB)
            J     = J+1
          END DO
 190      CONTINUE
          IF (JDQPR(J-1,JB) < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDQPR(J,JB),' < ending simulation time [TMEND=',TMEND,'] in '//PREFN(JB)
          END IF
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' precipitation file '//PREFN(JB)
        END IF
        CLOSE (NPT)

!****** Inflow temperatures

        WRITE (WIN,*) '  precipitation temperature'
        OPEN (UNIT=NPT,FILE=TPRFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
          READ (NPT,'(//)')
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0/8X,9F8.0)',END=200) JDAY, TPR
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TPRFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date',JDAY,' <= previous date ',JDAYO,' in '//TPRFN(JB)
            END IF
            J     = J+1
            JDAYO = JDAY
          END DO
 200      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TPRFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' precipitation temperature file '//TPRFN(JB)
        END IF

!****** Inflow concentrations

        IF (CONSTITUENTS .AND. NACPR(JW) > 0) THEN
          WRITE (WIN,*) '  precipitation concentrations'
          OPEN (UNIT=NPT,FILE=CPRFN(JB),IOSTAT=IERR,STATUS='OLD')
          IF (IERR == 0) THEN
            READ (NPT,'(//)')
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              READ (NPT,'(F8.0,1000(F8.0))',END=210) JDAY, (CPR,JC=1,NACPR(JW))
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CPRFN(JB)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CPRFN(JB)
              END IF
              JDAYO = JDAY
              J     = J+1
            END DO
 210        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CPRFN(JB)
            END IF
            CLOSE (NPT)
          ELSE
            CALL WARNINGS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' precipitation concentration file '//CPRFN(JB)
          END IF
        END IF
      END IF

!**** Upstream head

      IF (UH_EXTERNAL(JB)) THEN

!****** Elevations

        IFLAGEU = 0
        WRITE (WIN,*) '  upstream head elevations'
        OPEN (UNIT=NPT,FILE=EUHFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
          READ (NPT,'(//)')
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0/8X,9F8.0)',END=220) JDAY, EUH
            IF (JDAY >= TMSTRT .AND. IFLAGEU == 0) THEN
              ELDIFF = ABS(EUH-ELWS(US(JB)))
              IF (ELDIFF > 0.1) THEN
                CALL WARNINGS
                WRITE (WRN,FMTI) 'Boundary elevation and downstream water surface elevation difference > 0.1m for branch ',JB
              END IF
              IFLAGEU = 1
            END IF
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//EUHFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//EUHFN(JB)
            END IF
            JDAYO = JDAY
            J     = J+1
          END DO
 220      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//EUHFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' external upstream head file '//EUHFN(JB)
        END IF

!****** Temperatures

        WRITE (WIN,*) '  upstream head temperatures'
        OPEN (UNIT=NPT,FILE=TUHFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
          READ (NPT,'(//)')
          J     = J+1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0/8X,9F8.0)',END=230) JDAY, (TUH,K=2,KB(US(JB)))
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TUHFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TUHFN(JB)
            END IF
            J     = J+1
            JDAYO = JDAY
          END DO
 230      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TUHFN(JB)
         END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' external upstream head temperature file '//TUHFN(JB)
        END IF

!****** Constituent concentrations

        IF (CONSTITUENTS .AND. NAC > 0) THEN
          WRITE (WIN,*) '  upstream head concentrations'
          OPEN (UNIT=NPT,FILE=CUHFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR == 0) THEN
            READ (NPT,'(//)')
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              DO JC=1,NAC
                READ (NPT,*,END=240) JDAY, (CUH,K=2,KB(US(JB)))
                IF (J == 1 .AND. JDAY > TMSTRT) THEN
                  CALL ERRORS
                  WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CUHFN(JB)
                ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                  CALL ERRORS
                  WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CUHFN(JB)
                END IF
              END DO
              JDAYO = JDAY
              J     = J+1
            END DO
 240        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CUHFN(JB)
           END IF
            CLOSE (NPT)
          ELSE
            CALL ERRORS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' external upstream head concentration file '//CUHFN(JB)
          END IF
        END IF
      END IF

!**** Downstream head

      IF (DH_EXTERNAL(JB)) THEN

!****** Elevations

        IFLAGED = 0
        WRITE (WIN,*) '  downstream head elevations'
        OPEN (UNIT=NPT,FILE=EDHFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
          READ (NPT,'(//)')
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0/8X,9F8.0)',END=250) JDAY, EDH
            IF (JDAY >= TMSTRT .AND. IFLAGED == 0) THEN
              ELDIFF = ABS(EDH-ELWS(DS(JB)))
              IF (ELDIFF > 0.1) THEN
                CALL WARNINGS
                WRITE (WRN,FMTI) 'Boundary elevation and downstream water surface elevation difference > 0.1m for branch ',JB
              END IF
              IFLAGED = 1
            END IF
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//EDHFN(JB)
            ELSE IF (JDAY < JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//EDHFN(JB)
            END IF
            JDAYO = JDAY
            J     = J+1
          END DO
 250      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//EDHFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' external downstream head file '//EDHFN(JB)
        END IF

!****** Temperatures

        WRITE (WIN,*) '  downstream head temperatures'
        OPEN (UNIT=NPT,FILE=TDHFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
          READ (NPT,'(//)')
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(10F8.0:/(8X,9F8.0))',END=260) JDAY, (TDH,K=2,KB(DS(JB)))
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TDHFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TDHFN(JB)
            END IF
            J     = J+1
            JDAYO = JDAY
          END DO
 260      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TDHFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' external downstream head temperature file '//TDHFN(JB)
        END IF

!****** Constituent concentrations

        IF (CONSTITUENTS .AND. NAC > 0) THEN
          WRITE (WIN,*) '  downstream head concentrations'
          OPEN (UNIT=NPT,FILE=CDHFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR == 0) THEN
            READ (NPT,'(A)') (AID,I=1,3)
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
              DO JC=1,NAC
                READ (NPT,'(10F8.0/(:8X,9F8.0))',END=270) JDAY, (CDH,K=2,KB(DS(JB)))
                IF (J == 1 .AND. JDAY > TMSTRT) THEN
                  CALL ERRORS
                  WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CDHFN(JB)
                ELSE IF (JDAY < JDAYO .AND. J /= 1) THEN
                  CALL ERRORS
                  WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CDHFN(JB)
                END IF
              END DO
              JDAYO = JDAY
              J     = J+1
            END DO
 270        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CDHFN(JB)
            END IF
            CLOSE (NPT)
          ELSE
            CALL ERRORS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' external downstream head concentration file '//CDHFN(JB)
          END IF
        END IF
      END IF
    END DO
  END DO

! Tributaries

  IF (TRIBUTARIES) THEN
    DO JT=1,NTR
      WRITE (TRA,'(I0)') JT

!**** Inflows

      WRITE (WIN,*) '  tributary inflows'
      OPEN (UNIT=NPT,FILE=QTRFN(JT),STATUS='OLD',IOSTAT=IERR)
      IF (IERR == 0) THEN
        READ  (NPT,'(A)') HEADER
        J       =  1
        JDAYO   =  0.0
        DO WHILE (.TRUE.)
          READ (NPT,'(10F8.0/8X,9F8.0)',END=280) JDQTR(J,JT), QTR(J,JT)
          IF (JDAYO >= TMSTRT .AND. JDQTR(J,JT) <= TMEND) THEN
            IF (JBTR(JT) /= 0) QTRMX(JT,JBTR(JT)) = MAX(QTR(J,JT),QTRMX(JT,JBTR(JT)))
            IF (JDQTR(J,JT) > TMEND) THEN
              QTRS(JT,JBTR(JT)) = QTRS(JT,JBTR(JT))+QTR(J,JT)*(TMEND-JDAYO)
            ELSE
              QTRS(JT,JBTR(JT)) = QTRS(JT,JBTR(JT))+QTR(J,JT)*(JDQTR(J,JT)-JDAYO)
            END IF
          END IF
          IF (J == 1 .AND. JDQTR(J,JT) > TMSTRT) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Starting date ',JDQTR(J,JT),'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QTRFN(JT)
          ELSE IF (JDQTR(J,JT) <= JDAYO .AND. J /= 1) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Julian date ',JDQTR(J,JT),' <= previous date ',JDAYO,' in '//QTRFN(JT)
          END IF
          JDAYO = JDQTR(J,JT)
          J     = J+1
        END DO
 280    CONTINUE
        IF (JDQTR(J-1,JT) < TMEND) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Ending time ',JDQTR(J-1,JT),' < ending simulation time [TMEND=',TMEND,'] in '//QTRFN(JT)
        END IF
        CLOSE (NPT)
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTI) 'Could not open tributary ',JT,' inflow file '//QTRFN(JT)
      END IF

!**** Inflow temperatures

      WRITE (WIN,*) '  tributary inflow temperatures'
      OPEN (UNIT=NPT,FILE=TTRFN(JT),STATUS='OLD',IOSTAT=IERR)
      IF (IERR == 0) THEN
        READ  (NPT,'(A)') HEADER
        J       = 1
        JDAYO   = 0.0
        DO WHILE (.TRUE.)
          READ (NPT,'(10F8.0/8X,9F8.0)',END=290) JDAY, TTR
          IF (J == 1 .AND. JDAY > TMSTRT) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Starting date ',JDAY,' > simulation start date [TMSTRT=',TMSTRT,'] in '//TTRFN(JT)
          ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TTRFN(JT)
          END IF
          JDAYO = JDAY
          J     = J+1
        END DO
 290    CONTINUE
        IF (JDAY < TMEND) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TTRFN(JT)
        END IF
        CLOSE (NPT)
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTI) 'Could not open tributary ',JT,' inflow temperature file '//TTRFN(JT)
      END IF

!**** Inflow concentrations
      cin=0.0
      cin2=0.0
      CTRMIN(TRCN(1:NACTR(JT),JT),JT) = 1000.0
      CTRMAX(TRCN(1:NACTR(JT),JT),JT) = 0.0
      CTRAVG(TRCN(1:NACTR(JT),JT),JT) = 0.0
      IF (CONSTITUENTS .AND. NACTR(JT) > 0) THEN
        WRITE (WIN,*) '  tributary inflow concentrations'
        OPEN (UNIT=NPT,FILE=CTRFN(JT),IOSTAT=IERR,STATUS='OLD')
        IF (IERR == 0) THEN
          READ (NPT,'(//)')
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            READ (NPT,'(F8.0,1000(F8.0))',END=300) JDAY, (CIN(TRCN(JC,JT)),JC=1,NACTR(JT))
            IF (J == 1) TSTART = JDAY
            CTRMIN(TRCN(1:NACTR(JT),JT),JT) = MIN(CIN(TRCN(1:NACTR(JT),JT)),CTRMIN(TRCN(1:NACTR(JT),JT),JT))
            CTRMAX(TRCN(1:NACTR(JT),JT),JT) = MAX(CIN(TRCN(1:NACTR(JT),JT)),CTRMAX(TRCN(1:NACTR(JT),JT),JT))

    ! compute mass loading from trib inflows

             if(j==1)then
             open(npt+1,FILE=QTRFN(JT),STATUS='OLD')
             READ  (NPT+1,'(A)') HEADER
             QINO    = 0.0
             jj=1
             READ (NPT+1,'(10F8.0/8X,9F8.0)') JDQTR(J,JT), QTR(J,JT)
             endif
                do while(jdqtr(jj,jt) < jday)
                jj=jj+1
                READ (NPT+1,'(10F8.0/8X,9F8.0)',end=299) JDQtr(JJ,Jt), Qtr(JJ,Jt)
                enddo
                if(jday == jdqtr(jj,jt).or.jj==1)then
                qq=qtr(jj,jt)
                else
                qq=qtr(jj,jt)+((qtr(jj,jt)-qtr(jj-1,jt))/(jdqtr(jj,jt)-jdqtr(jj-1,jt)))*(jdqtr(jj,jt)-jday)
                endif
                if(qq<0.0)qq=0.0
                if(j>1)ctrload(trcn(1:NACtr(Jt),Jt),jt)=ctrload(trcn(1:NACtr(Jt),Jt),jt)+qq*(Cin(TRCN(1:NACTR(JT),JT))+CIN2(TRCN(1:NACTR(JT),JT)))*0.5*(JDAY-TDAY2)*86400./1000.  ! units of kg
                ! end of inflow loading

299         IF (J > 1) CTRAVG(TRCN(1:NACTR(JT),JT),JT) =  CTRAVG(TRCN(1:NACTR(JT),JT),JT)+(CIN(TRCN(1:NACTR(JT),JT))               &
                                                          +CIN2(TRCN(1:NACTR(JT),JT)))*0.5*(JDAY-TDAY2)
            CIN2(TRCN(1:NACTR(JT),JT)) = CIN(TRCN(1:NACTR(JT),JT))
            TDAY2                      = JDAY
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CTRFN(JT)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CTRFN(JT)
            END IF
            JDAYO = JDAY
            J     = J+1
          END DO
 300      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CTRFN(JT)
          END IF
          CLOSE (NPT)
          close(npt+1)
          ctrload(trcn(1:NACtr(Jt),Jt),jt)=ctrload(trcn(1:NACtr(Jt),Jt),jt)/(JDAY-TSTART)                              !CB 5/11/10 convert to units of kg/day
          CTRAVG(TRCN(1:NACTR(JT),JT),JT) = CTRAVG(TRCN(1:NACTR(JT),JT),JT)/(JDAY-TSTART)                              !SW 01/07/01
        ELSE
          CALL WARNINGS
          WRITE (ERR,FMTI) 'Could not open tributary ',JT, ' inflow concentration file '//CTRFN(JT)
        END IF
      END IF
    END DO
  END IF

! Outflows

  WRITE (WIN,*) '  outflows'
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      WRITE (BRA,'(I0)') JW
      IF (DQ_EXTERNAL(JB)) THEN
        IF (NSTR(JB) /= 0) THEN
          OPEN (UNIT=NPT,FILE=QOTFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR == 0) THEN
            READ  (NPT,'(A)') HEADER
            J       = 1
            JDAYO   = 0.0
            QSTRO   = 0.0
            DO WHILE (.TRUE.)
              READ (NPT,'(10F8.0/(9X,9F8.0))',END=310) JDAY, (QSTR(JS),JS=1,NSTR(JB))

            DO JS=1,NSTR(JB)
            IF (QSTR(JS)<0.0) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'QSTR cannot be negative:',QSTR(JS),' at date [JDAY=',JDAY,'] in '//QOTFN(JB)
            END IF
            ENDDO

              IF (JDAYO >= TMSTRT .AND. JDAYO <= TMEND) THEN
                DO JS=1,NSTR(JB)
                  QSTRMX(JB) = MAX(QSTR(JS),QSTRMX(JB))
                  IF (JDAY > TMEND) THEN
                    QSTRS(JB) = QSTRS(JB)+QSTRO(JS)*(TMEND-JDAYO)
                  ELSE
                    QSTRS(JB) = QSTRS(JB)+QSTRO(JS)*(JDAY-JDAYO)
                  END IF
                END DO
              END IF
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QOTFN(JB)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//QOTFN(JB)
              END IF
              JDAYO = JDAY
              QSTRO = QSTR
              J     = J+1
            END DO
 310        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//QOTFN(JB)
            END IF
            CLOSE (NPT)
          ELSE
            CALL ERRORS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' outflow file '//QOTFN(JB)
          END IF
        END IF
      END IF
    END DO
  END DO

! Withdrawals

  WRITE (WIN,*) '  withdrawals'
  IF (WITHDRAWALS) THEN
    OPEN (UNIT=NPT,FILE=QWDFN,STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      READ (NPT,'(//)')
      J     = 1
      JDAYO = 0.0
      DO WHILE (.TRUE.)
        READ (NPT,'(10F8.0/8X,9F8.0)',END=320) JDAY, (QWD(J,JWD),JWD=1,NWD)
        IF (JDAYO >= TMSTRT .AND. JDAY <= TMEND) THEN
          JWD = 1
          DO JB=1,NBR
            DO JWB=1,NWBR(JB)
              QWDMX(JWB,JB) = MAX(QWD(J,JWD),QWDMX(JWB,JB))
              IF (JDAY > TMEND) THEN
                QWDS(JWB,JB) = QWDS(JWB,JB)+QWD(J,JWD)*(TMEND-JDAYO)
              ELSE
                QWDS(JWB,JB) = QWDS(JWB,JB)+QWD(J,JWD)*(JDAY-JDAYO)
              END IF
              JWD = JWD+1
            END DO
          END DO
        END IF
        IF (J == 1 .AND. JDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Starting date ',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QWDFN
        ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//QWDFN
        END IF
        JDAYO = JDAY
        J     = J+1
      END DO
 320  CONTINUE
      IF (JDAY < TMEND) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//QWDFN
      END IF
      CLOSE (NPT)
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTA) 'Could not open withdrawal file '//QWDFN
    END IF
  END IF

! Gates

  IF (GATES) THEN
    WRITE (WIN,*) '  gates'
    OPEN (UNIT=NPT,FILE=QGTFN,IOSTAT=IERR,STATUS='OLD')
    IF (IERR == 0) THEN
      READ (NPT,'(//)')
      J     = 1
      JDAYO = 0.0
      DO WHILE (.TRUE.)
        READ (NPT,'(10F8.0,(:/8X,9F8.0))',END=120) JDAY, (BGATE(JG),JG=1,NGT)
        IF (J == 1 .AND. JDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Starting date ',JDAY,' > simulation start date [TMSTRT=',TMSTRT,'] in '//QGTFN
        ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//QGTFN
        END IF
        J     = J+1
        JDAYO = JDAY
      END DO
 120  CONTINUE
      IF (JDAY < TMEND) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//QGTFN
      END IF
      CLOSE (NPT)
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTA) 'Could not open gate filename '//QGTFN
    END IF
  END IF

! Shading

  WRITE (WIN,*) '  shading'
  OPEN  (NPT,FILE=SHDFN,IOSTAT=IERR,STATUS='OLD')                                                                      !TC 10/22/02
  IF (IERR == 0) THEN
    READ  (NPT,'(//)')
    DO I=1,IMX
      READ (NPT,'(8X,F8.0)') SHADEIN(I)
      IF (SHADEIN(I) == 0.0) THEN
        CALL   WARNINGS
        WRITE (WRN,FMTI) '100% shade for segment ',I
      ELSE IF (SHADEIN(I) > 1.0) THEN
        CALL   ERRORS
        WRITE (ERR,FMTI) 'shade > 1.0 for segment ',I
      END IF
      IF (SHADEIN(I) < 0.0) THEN
        BACKSPACE (NPT)
        READ (NPT,'(16X,28F8.0)') TTLB(I),TTRB(I),CLLB(I),CLRB(I),SRLB1(I),SRLB2(I),SRRB1(I),SRRB2(I),(TOPO(I,J),J=1,IANGLES),     &
                                  SRFJDAY1(I),SRFJDAY2(I)
        IF (TTLB(I) < EL(KB(I),I)) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Left bank vegetative shade cover elevation less than bottom channel elevation for segment ',I
        END IF
        IF (TTRB(I) < EL(KB(I),I)) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Right bank vegetative shade cover elevation less than bottom channel elevation for segment ',I
        END IF
        IF (CLLB(I) <= 0.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Left bank distance to vegetative shade cover <= 0 m for segment ',I
        END IF
        IF (CLRB(I) <= 0.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Right bank distance to vegetative shade cover <= 0 m for segment ',I
        END IF
        IF (SRLB1(I) < 0.0 .OR. SRLB1(I) > 1.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shading reduction factor LB1 < 0 or > 1 for segment ',I
        END IF
        IF (SRLB2(I) < 0.0 .OR. SRLB2(I) > 1.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shading reduction factor LB2 < 0 or > 1 for segment ',I
        END IF
        IF (SRRB1(I) < 0.0 .OR. SRRB1(I) > 1.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shading reduction factor RB1 < 0 or > 1 for segment ',I
        END IF
        IF (SRRB2(I) < 0.0 .OR. SRRB2(I) > 1.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shading reduction factor RB2 < 0 or > 1 for segment ',I
        END IF
        DO JJ=1,IANGLES
          IF(TOPO(I,JJ) < 0 .OR. TOPO(I,JJ) > 3.1419/2.0) THEN
            CALL   ERRORS
            WRITE (ERR,FMTI) 'Shading topographic angle < 0 or > 3.1419/2 for segment ',I
          END IF
        END DO
        IF (SRFJDAY1(I) > SRFJDAY2(I)) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shade reduction factor JDAY1 (leaf out) > JDAY2 (leaf off) for segment ',I
        END IF                                                                                                         !SW 04/03/02
      END IF                                                                                                           !SW 04/03/02
    END DO                                                                                                             !SW 04/03/02
  ELSE                                                                                                                 !TC 10/22/02
    CALL ERRORS                                                                                                        !TC 10/22/02
    WRITE (ERR,FMTA) 'Could not open shade filename '//SHDFN                                                           !TC 10/22/02
  END IF                                                                                                               !TC 10/22/02
  CLOSE (NPT)                                                                                                          !SW 04/03/02

! Restart file

  WRITE (WIN,*) '  restart'
  IF (RESTART_IN) THEN
    OPEN (NPT,FILE=RSIFN,STATUS='OLD',IOSTAT=IERR)
    IF (IERR /= 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTA) 'Could not open restart file '//RSIFN
    ELSE
      CLOSE (NPT)
    END IF
  END IF

! Water balance

  JTR = 0; JWD = 0                                                                                                     !TC 07/10/03
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      DO JT=1,NTBR(JB)
        JTR           = JTR+1                                                                                          !TC 07/10/03
        QTRMXB(JB)    = MAX(QTRMXB(JB),QTRMX(JTR,JB))
        QINMXB(JB)    = MAX(QINMXB(JB),QTRMXB(JB))
        QTRAV(JTR,JB) = QTRS(JTR,JB)/(TMEND-TMSTRT)
        QTRAVB(JB)    = QTRAVB(JB)+QTRAV(JTR,JB)
      END DO
      DO J=1,NWBR(JB)
        JWD           = JWD+1                                                                                          !TC 07/10/03
        QWDMXB(JB)    = MAX(QWDMXB(JB),QWDMX(JWD,JB))
        QOTMXB(JB)    = MAX(QOTMXB(JB),QWDMXB(JB))
        QWDAV(JWD,JB) = QWDS(JWD,JB)/(TMEND-TMSTRT)
        QWDAVB(JB)    = QWDAVB(JB)+QWDAV(JWD,JB)
      END DO
      QINAV(JB)  = QINS(JB) /(TMEND-TMSTRT)
      QPRAV(JB)  = QPRS(JB) /(TMEND-TMSTRT)
      QDTAV(JB)  = QDTS(JB) /(TMEND-TMSTRT)
      QSTRAV(JB) = QSTRS(JB)/(TMEND-TMSTRT)
      QINAVB(JB) = QINAV(JB) +QDTAV(JB)+QPRAV(JB)+QTRAVB(JB)
      QOTAVB(JB) = QSTRAV(JB)+QWDAVB(JB)
      QOTAVW(JW) = QOTAVW(JW)+QOTAVB(JB)
      QINAVW(JW) = QINAVW(JW)+QINAVB(JB)
      QINMXB(JB) = MAX(QINMX(JB), QDTMX(JB))
      QOTMXB(JB) = MAX(QSTRMX(JB),QWDMXB(JB))
      QOTMXW(JW) = MAX(QOTMXW(JW),QOTMXB(JB))
      QINMXW(JW) = MAX(QINMXW(JW),QINMXB(JB))
    END DO
  END DO
!  DO JW=1,NWB
!    DO JB=BS(JW),BE(JW)
!      IF (UQ_INTERNAL(JB)) THEN
!        QINAVB(JB) = QINAVB(UQB(JB))+QDTAV(JB)+QPRAV(JB)+QTRAVB(JB)
!        QINMXB(JB) = MAX(QINMXB(UQB(JB)),QINMXB(JB))
!        QINMXW(JW) = MAX(QINMXB(UQB(JB)),QINMXB(JB))
!        QINAVW(JW) = QINAVB(UQB(JB))
!      END IF
!      IF (DAM_FLOW(JB)) THEN
!        QINAVB(JB) = QOTAVB(UQB(JB))
!        QINMXB(JB) = MAX(QOTMXB(UQB(JB)),QINMX(JB))
!        QINMXW(JW) = MAX(QINMXW(JW),QINMXB(JB))
!        QINAVW(JW) = QSTRAV(UQB(JB))+QDTAV(JB)+QPRAV(JB)+QTRAVB(JB)
!      END IF
!    END DO
!  END DO

!***********************************************************************************************************************************
!*                                                      Task 4: Control File Inputs                                               **
!***********************************************************************************************************************************

  WRITE (WIN,*) 'Control file'
  WRITE (WIN,*) '  time controls'
  IF (TMEND < TMSTRT) THEN
    CALL ERRORS
    WRITE (ERR,FMTF) 'Starting time [TMSTRT=',TMSTRT,'] < ending time [TMEND=',TMEND,']'
  END IF

! Autostepping controls

  IF (NDLT <= 0) THEN
    CALL ERRORS
    WRITE (ERR,FMTI) 'Number of timestep intervals [NDLT=',NDLT,'] <= 0'
  END IF
  IF (DLTD(1) > TMSTRT) THEN
    CALL ERRORS
    WRITE (ERR,FMTF) 'Timestep date [DLTD(1)=',DLTD(1),'] > starting time [TMSRTRT=',TMSTRT,']'
  END IF
  DO J=1,NDLT-1
    IF (DLTD(J+1) < DLTD(J)) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Timestep date [DLTD(',J,')=',DLTD(J),'] < timestep date [DLTD(',J+1,')=',DLTD(J+1),']'
    END IF
  END DO
  DO J=1,NDLT
    IF (DLTD(J) < 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Timestep date [DLTD(',J,     ')=',DLTD(J),  '] <  1.0'
    END IF
    IF (DLTMAX(J) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Maximum timestep [DLTMAX(',J,')=',DLTMAX(J),'] <= 0.0'
    END IF
    IF (DLTF(J) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Timestep fraction [DLTF(',J,  ')=',DLTF(J), '] <= 0.0'
    ELSE IF (DLTF(J) >= 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Timestep fraction [DLTF(',J,  ')=',DLTF(J), '] >= 1.0'
    END IF
  END DO
  IF (DLTMIN <= 0.0) THEN
    CALL ERRORS
    WRITE (ERR,FMTF) 'Minimum timestep [DLTMIN=',DLTMIN,'] <= 0.0'
  END IF
  IF (DLTMIN > 0.0 .AND. DLTMIN < 0.1) THEN
    CALL WARNINGS
    WRITE (WRN,FMTF) 'Minimum timestep [DLTMIN=',DLTMIN,'] <  0.1'
  END IF
  DO JW=1,NWB
    IF (VISC(JW) /= '      ON' .AND. VISC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Vertical eddy viscosity limitation control [VISC='//VISC(JW)(6:8)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (CELC(JW) /= '      ON' .AND. CELC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Internal gravity wave limitation control   [CELC='//CELC(JW)(6:8)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
  END DO

! Channel slope and minimum number of layers for active segment

  WRITE (WIN,*) '  channel slope'
  DO JB=1,NBR
    IF (SLOPE(JB) > 0.1) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Channel slope [SLOPE=',SLOPE(JB),'] > 0.1 for branch ',JB
    ELSE IF (SLOPE(JB) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Channel slope [SLOPE=',SLOPE(JB),'] < 0.0 for branch ',JB
    END IF
    IF (NL(JB) > 3) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Minimum number of layers for segment to be active [NL=',NL(JB),'] > 3'
    END IF
    IF (NL(JB) <= 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Minimum number of layers for segment to be active [NL=',NL(JB),'] <= 0'
    END IF
  END DO

! Branch linkage

  DO JW=1,NWB
    BRANCH_FOUND = .FALSE.
    DO JB=BS(JW),BE(JW)
      IF (JBDN(JW) == JB) THEN
        BRANCH_FOUND = .TRUE.
        EXIT
      END IF
    END DO
    IF (.NOT. BRANCH_FOUND) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Downstream branch [JBDN=',JBDN(JW),'] is not located in waterbody ',JW
    END IF
  END DO

! Temperature intitial condition

  WRITE (WIN,*) '  initial conditions'
  DO JW=1,NWB
    IF (T2I(JW) < -2.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Initial temperature [T2I=',T2I(JW),'] <= 0.0 for waterbody ',JW
    END IF
  END DO

! Ice initial condition

  DO JW=1,NWB
    IF (ICE_CALC(JW)) THEN
      IF (ICEI(JW) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Ice initial thickness [ICEI=',ICEI(JW),'] <= 0.0 for waterbody ',JW
      ELSE IF (ICEI(JW) > 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Ice initial thickness [ICEI=',ICEI(JW),'] >  0.1 for waterbody ',JW
      END IF
      IF (ICEI(JW) > 0.0 .AND. ICEI(JW) < ICEMIN(JW)) THEN
        CALL WARNINGS
        WRITE (WRN,FMT2FI) 'Ice initial thickness [ICEI=',ICEI(JW),'] < minimum ice thickness [ICEMIN=',ICEMIN(JW), '] for '//     &
                           'waterbody ',JW
      END IF
    END IF
  END DO

! Water body type

  WRITE (WIN,*) '  waterbody type'
  DO JW=1,NWB
    IF (WTYPEC(JW) == '    SALT') THEN
      IF (.NOT.CONSTITUENTS) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Since water type is SALT, constituents [CCC] must be set to "ON" for waterbody ',JW
      ELSE IF (.NOT.DISSOLVED_SOLIDS) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Since water type is SALT, CAC(4) must be set to "ON" for waterbody ',JW
      END IF
    ELSE IF (WTYPEC(JW) /= '   FRESH' .AND. WTYPEC(JW) /= '    SALT') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Waterbody type [WTYPEC(',JW,')='//WTYPEC(JW)//'] /= "   FRESH" or "    SALT" for waterbody ',JW
    END IF
  END DO

! Calculations

  WRITE (WIN,*) '  calculations'
  DO JW=1,NWB
    IF (VBC(JW) /= '      ON' .AND. VBC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Volume balance control   [VBC='//VBC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (MBC(JW) /= '      ON' .AND. MBC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Mass balance control     [MBC='//MBC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (EBC(JW) /= '      ON' .AND. EBC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Energy balance control   [EBC='//EBC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (PQC(JW) /= '      ON' .AND. PQC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Inflow placement control [PQC='//PQC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (EVC(JW) /= '      ON' .AND. EVC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Evaporation control      [EVC='//EVC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (PRC(JW) /= '      ON' .AND. PRC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Precipitation control    [PRC='//PRC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
  END DO

! Interpolation controls

  WRITE (WIN,*) '  interpolation'
  DO JW=1,NWB
    IF (METIC(JW) /= '      ON' .AND. METIC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Meteorology interpolation control [METIC='//METIC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
  END DO
  DO JB=1,NBR
    IF (QINIC(JB) /= '      ON' .AND. QINIC(JB) /= '     OFF' .AND. UQ_EXTERNAL(JB)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Inflow interpolation control [QINIC='//QINIC(JB)//'] /= " ON" or "OFF" for branch ',JB
    END IF
    IF (HDIC(JB) /= '      ON' .AND. HDIC(JB) /= '     OFF' .AND. (UH_EXTERNAL(JB) .OR. DH_EXTERNAL(JB))) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Head boundary interpolation control [HDIC='//HDIC(JB)//'] /= " ON" or "OFF" for branch ',JB
    END IF
    IF (DTRIC(JB) /= '      ON' .AND. DTRIC(JB) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Distrubuted tributary inflow interpolation control [DTRIC='//DTRIC(JB)//'] /= " ON" or "OFF" for branch ',JB
    END IF
    DO JS=1,NSTR(JB)
      IF (DQ_EXTERNAL(JB)) THEN
        IF (STRIC(JS,JB) /= '      ON' .AND. STRIC(JS,JB) /= '     OFF') THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Outflow interpolation control [QOUTIC='//STRIC(JS,JB)//'] /= " ON" or "OFF" for structure ',JS,        &
                           ' in branch ',JB
        END IF
      END IF
    END DO
  END DO
  DO JT=1,NTR
    IF (TRIC(JT) /= '      ON' .AND. TRIC(JT) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Tributary inflow interpolation control [TRIC='//TRIC(JT)//'] /= " ON" or "OFF" for tributary ',JT
    END IF
  END DO
  DO JWD=1,NWD
    IF (WDIC(JWD) /= '      ON' .AND. WDIC(JWD) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Withdrawal interpolation control [WDIC='//WDIC(JWD)//'] /= " ON" or "OFF" for withdrawal ',JWD
    END IF
  END DO

! Dead sea

  WRITE (WIN,*) '  dead sea'
  DO JW=1,NWB
    IF (WINDC(JW) /= '      ON' .AND. WINDC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Wind control [WINDC='//WINDC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (QINC(JW) /= '      ON' .AND. QINC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Inflow control [QINC='//QINC(JW)//'] /= " ON" or "OFF" for waterbody ', JW
    END IF
    IF (QOUTC(JW) /= '      ON' .AND. QOUTC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Outflow control [QOUTC='//QOUTC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (HEATC(JW) /= '      ON' .AND. HEATC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Surface heat exchange control [HEATC='//HEATC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
  END DO

!** Surface heat exchange

  WRITE (WIN,*) '  surface heat exchange'
  DO JW=1,NWB
    IF (SLHTC(JW) /= '    TERM' .AND. SLHTC(JW) /= '      ET') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Heat exchange solution control [SLHTC='//SLHTC(JW)//'] /= "    TERM" or "      ET" for waterbody ',JW
    END IF

! Radiation and evaporation

    IF (SROC(JW) /= '      ON' .AND. SROC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Solar radiation input control [SROC='//SROC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (AFW(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed function coefficient [AFW=',AFW(JW),'] < 0 for waterbody ',JW
    END IF
    IF (BFW(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed function coefficient [BFW=',BFW(JW),'] < 0 for waterbody ',JW
    END IF
    IF (CFW(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed function coefficient [CFW=',CFW(JW),'] < 0 for waterbody ',JW
    END IF
    IF (WINDH(JW) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed measurement height [WINDH=',WINDH(JW),'] <= 0m for waterbody ',JW
    END IF
    IF (WINDH(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed measurement height [WINDH=',WINDH(JW),'] < 0m for waterbody ',JW
    END IF
    IF (RHEVC(JW) /= '      ON' .AND. RHEVC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Ryan-Harleman evaporation control [RHEVC='//RHEVC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (FETCHC(JW) /= '      ON' .AND. FETCHC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Fetch computation control         [FETCHC='//FETCHC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
  END DO

! Ice cover

  WRITE (WIN,*) '  ice cover'
  DO JW=1,NWB
    IF (ICEC(JW) /= '      ON' .AND. ICEC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Ice control [ICEC='//ICEC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (SLICEC(JW) /= '  DETAIL' .AND. SLICEC(JW) /= '  SIMPLE') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Ice solution control [SLICEC='//SLICEC(JW)//'] /= "  DETAIL" or "  SIMPLE" for waterbody ',JW
    END IF
    IF (ALBEDO(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice albedo [ALBEDO=',ALBEDO(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (HWI(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Water/ice heat exchange [HWI=',HWI(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (BETAI(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice solar radiation absorption [BETAI=',BETAI(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (BETAI(JW) > 0.9) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Ice solar radiation absorption [BETAI=',BETAI(JW),'] > 0.9 for waterbody ',JW
    END IF
    IF (GAMMAI(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice solar radiation extinction [GAMMAI=',GAMMAI(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (ICEMIN(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Minimum ice thickness for ice formation [ICEMIN=',ICEMIN(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (ICEMIN(JW) > 0.1) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Minimum ice thickness for ice formation [ICEMIN=',ICEMIN(JW),'] > 0.1 for waterbody ',JW
    END IF
    IF (ICET2(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice formation temperature [ICET2=',ICET2(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (ICET2(JW) > 4.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice formation temperature [ICET2=',ICET2(JW),'] > 4.0 for waterbody ',JW
    END IF
  END DO

! Transport solution

  WRITE (WIN,*) '  transport solution'
  DO JW=1,NWB
    IF (SLTRC(JW) /= 'QUICKEST' .AND. SLTRC(JW) /= '  UPWIND' .AND. SLTRC(JW) /= 'ULTIMATE') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Transport solution [SLTRC='//SLTRC(JW)//'] /= "QUICKEST" or "  UPWIND" or "ULTIMATE" for waterbody ',JW
    END IF
    IF (THETA(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Vertical advection time weighting [THETA=',THETA(JW),'] < 0.0 for waterbody ',JW
    ELSE IF (THETA(JW) > 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Vertical advection time weighting [THETA=',THETA(JW),'] > 1.0 for waterbody ',JW
    END IF
  END DO

! Hydraulic coefficients

  WRITE (WIN,*) '  hydraulic coefficients'
  DO JW=1,NWB
    IF (AX(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Longitudinal eddy viscosity [AX=',AX(JW),'] < 0.0 for waterbody ',JW
    ELSE IF (AX(JW) < 0.01 .OR. AX(JW) > 100.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Longitudinal eddy viscosity [AX=',AX(JW),'] < 0.01 or > 100 for waterbody ',JW
    END IF
    IF (DX(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Longitudinal eddy diffusivity [DX=',DX(JW),'] < 0.0 for waterbody ',JW
    ELSE IF (DX(JW) < 0.01 .OR. DX(JW) > 100.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Longitudinal eddy diffusivity [DX=',DX(JW),'] < 0.01 or > 100 for waterbody ',JW
    END IF
      DXMIN=1.E6
      DO JB=BS(JW),BE(JW)
      DO I=US(JB),DS(JB)
      DXMIN=MIN(DLX(I),DXMIN)
      ENDDO
      ENDDO
      DXTHEORY=5.84e-4*DXMIN**1.1
      IF (DX(JW) > DXTHEORY*10.) THEN
      CALL WARNINGS
      WRITE (WRN,"('Longitudinal eddy diffusivity [DX=',f8.3,'] 10x > DX-theory=',f8.4,' for waterbody ',i4,':Model may be unstable')")dx(jw),dxtheory,JW
      ENDIF
      IF (AX(JW) > DXTHEORY*10.) THEN
      CALL WARNINGS
      WRITE (WRN,"('Longitudinal eddy viscosity [AX=',f8.3,'] 10x > AX-theory=',f8.4,' for waterbody ',i4,':Model may be unstable')")ax(jw),dxtheory,JW
      ENDIF
    IF (CBHE(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,'(A,E8.1,A,I0)') 'Coefficient of bottom heat exchange [CBHE=',CBHE(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (CBHE(JW) > 3.0) THEN
      CALL WARNINGS
      WRITE (WRN,'(A,E8.1,A,I0)') 'Coefficient of bottom heat exchange [CBHE=',CBHE(JW),'] > 3.0 for waterbody ',JW    !TC 09/30/03
    END IF
    IF (TSED(JW) > 30 .OR. TSED(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment temperature [TSED=',TSED(JW),'] < 0 or > 30 for waterbody ',JW
    END IF
    IF (FI(JW) > 0.05) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Sediment temperature adjustment [FI=',FI(JW),'] > 0.05 for waterbody ',JW
    END IF
    IF (FI(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment temperature adjustment [FI=',FI(JW),'] < 0 for waterbody ',JW
    END IF
    IF (TSEDF(JW) < 0.00 .OR. TSEDF(JW) > 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment temperature factor [TSEDF=',TSEDF(JW),'] < 0 or > 1 for waterbody ',JW
    END IF
  END DO

! Vertical eddy viscosity formulation

  WRITE (WIN,*) '  vertical eddy viscosity'
  DO JW=1,NWB
    !NOTE ADD CHECKS FOR FBC,E,ARODI,STRCKLR,BOUNDFR,TKECAL !!

    IF (AZC(JW) /= '    NICK' .AND. AZC(JW) /= '      W2' .AND. AZC(JW) /= '     W2N' .AND. AZC(JW) /= '     RNG' .AND.            &
        AZC(JW) /= '   PARAB' .and. azc(jw) /= '     TKE' .and. azc(jw) /= '    TKE1') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Vertical eddy viscosity formulation [AZC='//AZC(JW)//'] /= "      W2", "     W2N", "    NICK",'//          &
                       ' "     RNG", or "   PARAB" or "     TKE" or "    TKE1" for waterbody ',JW
    END IF
    IF (AZC(JW) == '    TKE1')THEN
        IF(FBC(JW) < 0 .OR. FBC(JW) >3)THEN
         CALL ERRORS
         WRITE (ERR,FMTI) 'FBC for the TKE1 model [FBC=',FBC(JW),'] must be either 1, 2, or 3'//          &
                       ' for waterbody ',JW
        ENDIF
        IF(AZE(JW) /= 9.535)THEN
         CALL WARNINGS
         WRITE (WRN,FMTFI) 'E for the TKE1 model [E=',AZE(JW),'] is set at a value other than the default of 9.535'//          &
                       ' for waterbody ',JW
        ENDIF
        IF(AZE(JW) < 0.0 .or. aze(jw) > 30.0)THEN
         CALL ERRORS
         WRITE (ERR,FMTFI) 'E for the TKE1 model [E=',AZE(JW),'] is set at a value < 0 or > 30;default is 9.535'//          &
                       ' for waterbody ',JW
        ENDIF
        IF(ARODI(JW) < 0.01 .OR. ARODI(JW) > 0.7)THEN
         CALL ERRORS
         WRITE (ERR,FMTFI) 'ARODI for the TKE1 model [ARODI=',ARODI(JW),'] is set at a value < 0.01 or > 0.7; default is 0.43 for FBC=1'//          &
                       ' and 0.07 for FBC=2 (ARODI not used for FBC=3) for waterbody ',JW
        END IF
        IF(STRCKLR(JW) < 0.0 .OR. STRCKLR(JW) > 50.0)THEN
         CALL ERRORS
         WRITE (ERR,FMTFI) 'STRCKLR for the TKE1 model [STRCKLR=',STRCKLR(JW),'] is set at a value < 0.0 or > 50.; default is 24.0'//          &
                       ' for waterbody ',JW
        END IF
        IF(BOUNDFR(JW) < 0.0 .OR. BOUNDFR(JW) > 20.0)THEN
         CALL ERRORS
         WRITE (ERR,FMTFI) 'BOUNDFR for the TKE1 model [BOUNDFR=',BOUNDFR(JW),'] is set at a value < 0.0 or > 20.; default is 10.0'//          &
                       ' for waterbody ',JW
        END IF
        IF (tkecal(jw) /= '     IMP'.and. tkecal(jw) /= '     EXP') THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'TKE1 calculation technique [TKECAL='//TKECAL(JW)//'] /= "     IMP" or "     EXP" for waterbody ',JW
        END IF

    ENDIF


    IF (AZSLC(JW) /= '     IMP' .AND. AZSLC(JW) /= '     EXP') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Vertical eddy viscosity numerical solution scheme [AZSLC='//AZSLC(JW)//'] /= "IMP" or "EXP" for '//        &
                       'waterbody ',JW
    END IF
    IF (AZMAX(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Maximum vertical eddy viscosity [AZMAX]=',AZMAX(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (AZMAX(JW) > 1.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Maximum vertical eddy viscosity [AZMAX]=',AZMAX(JW),'] > 1.0 for waterbody ',JW
    END IF
    DO JB=BS(JW),BE(JW)
      IF (SLOPE(JB) == 0.0) THEN
        IF (AZC(JW) /= '     W2N' .AND. AZC(JW) /= '      W2' .AND. (.NOT.UH_EXTERNAL(JB) .AND. .NOT.DH_EXTERNAL(JB))) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Vertical eddy viscosity formulation [AZC='//AZC(JW)//'] should probably be "      W2" or '//           &
                           '"     W2N" for waterbody ',JW
        END IF
        IF (AZMAX(JW) > 0.01 .AND. AZC(JW) == 'EXP') THEN
          CALL WARNINGS
          WRITE (WRN,'(A,1E10.3,A,I0)') 'Explicit vertical eddy viscosity [AZC="     EXP"] and a maximum vertical eddy '//         &
                                        'viscosity of [AZMAX=',AZMAX(JW),'] will result in low timesteps for waterbody ',JW
        END IF
      END IF
    END DO
  END DO

! Bottom friction type

  WRITE (WIN,*) '  bottom friction'
  DO JW=1,NWB
    IF (FRICC(JW) /= '   CHEZY' .AND. FRICC(JW) /= '    MANN') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Bottom friction type [FRICC='//FRICC(JW)//'] /= "CHEZY" or "MANN" for waterbody ',JW
    END IF
    DO JB=BS(JW),BE(JW)
      DO I=US(JB),DS(JB)
        IF (FRICC(JW) == '   CHEZY' .AND. (FRIC(I) < 30.0 .OR. FRIC(I) > 100.0)) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Chezy coefficient [FRIC=',FRIC(I),'] < 20 or > 100 for segment ',I
        ELSE IF (FRICC(JW) == '    MANN' .AND. (FRIC(I) < 0.0 .OR. FRIC(I) > 0.20)) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Mannings N [FRIC=',FRIC(I),'] < 0 or > 0.2 for segment ',I
        END IF
      END DO
    END DO
  END DO

  ! WIND ROUGHNESS HEIGHT Z0
  WRITE (WIN,*) '  wind roughness z0'
  DO JW=1,NWB
    IF (Z0(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,'("Wind roughness height Z0 < 0.0 m",1X," Z0=",F10.5," FOR JW=",I4)') Z0(JW),JW
    ELSEIF (Z0(JW) > 1.0)THEN
      CALL WARNINGS
      WRITE (WRN,'("Wind roughness height Z0 > 1.0 m",1X," Z0=",F10.5," FOR JW=",I4)') Z0(JW),JW
    ELSEIF(Z0(JW) < 1.0E-7)THEN
      CALL ERRORS
      WRITE (ERR,'("Wind roughness height Z0 < 1.E-7 m",1X," Z0=",F10.5," FOR JW=",I4)') Z0(JW),JW
    END IF
  ENDDO



! Outlet structures

  WRITE (WIN,*) '  inflow/outflow'
  WRITE (WIN,*) '    outlet structures'
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      IF (DQ_EXTERNAL(JB)) THEN
        ID = DS(JB)
        IF (NSTR(JB) == 0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'No structures defined [NSTR=0] for downstream flow boundary in branch ',JB
        END IF
        DO JS=1,NSTR(JB)
          KSTR = 2
          DO WHILE (EL(KSTR,DS(JB)) >= ESTR(JS,JB) .AND. KSTR < KB(ID)+1)
            KSTR = KSTR+1
          END DO
          KSTR = KSTR-1
          IF (KTSTR(JS,JB) < 2) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' selective withdrawal top layer [KTSTR=',KTSTR(JS,JB),'] < 2 for branch ',JB
          END IF
          IF (KTSTR(JS,JB) > KBSTR(JS,JB)) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' selective withdrawal top layer [KTSTR=',KTSTR(JS,JB),'] > bottom selective '//      &
                             'withdrawal layer [KBSTR=',KBSTR(JS,JB),'] in branch ',JB
          END IF
          IF (KTSTR(JS,JB) > KSTR) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' selective withdrawal top layer [KTSTR=',KTSTR(JS,JB),'] > centerline elevation '//  &
                             '[KSTR=',KSTR,'] in branch ',JB
          END IF
          IF (KBSTR(JS,JB) > KB(ID)) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' selective withdrawal bottom layer [KBSTR=',KBSTR(JS,JB),'] > bottom active layer '//&
                             '[KB=',KB(ID),'] in branch ',JB
          END IF
          IF (KBSTR(JS,JB) < KSTR) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' Selective withdrawal bottom layer [KBSTR=',KBSTR(JS,JB),'] < centerline layer'//    &
                             ' elevation [KSTR=',KSTR,'] in branch ',JB
          END IF
          IF (SINKC(JS,JB) /= '    LINE' .AND. SINKC(JS,JB) /= '   POINT') THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' sink type [SINKC='//SINKC(JS,JB)//'] /= "    LINE" or "   POINT" in branch ',JB
          END IF
          IF (SINKC(JS,JB) == '    LINE') THEN
            IF (WSTR(JS,JB) <= 0.0) THEN
              CALL ERRORS
              WRITE (ERR,FMTIF) 'Structure ',JS,' width [WSTR=',WSTR(JS,JB),'] <= 0 in branch ',JB
            END IF
            IF (WSTR(JS,JB) > B(KSTR,ID)) THEN
              CALL ERRORS
              WRITE (ERR,'(A,I0,A,2(F0.3,2(A,I0,A,I0,A,F0.2,A,I0)))') 'Structure ',JS,' width [WSTR=',WSTR(JS,JB),'] > cell width [B(',KSTR, ',',ID,&
                                                     ')=',B(KSTR,ID),'] in branch ',JB
            END IF
          END IF
          IF (ESTR(JS,JB) < EL(KB(ID+1)+1,ID)) THEN
            CALL ERRORS
            WRITE (ERR,FMTFI) 'Selective withdrawal elevation [ESTR=',ESTR(JS,JB),'] < the bottom active cell elevation [EL(',     &
                               KB(ID+1),')=',EL(KB(ID)+1,ID),'] for structure ',JS
          END IF
          IF (ESTR(JS,JB) > EL(2,ID)) THEN
            CALL ERRORS
            WRITE (ERR,FMT2FI) 'Selective withdrawal elevation [ESTR=',ESTR(JS,JB),'] > grid top layer elevation [EL(2)=',EL(2,ID),&
                               '] for structure ',JS
          END IF
        END DO
      ELSE IF (.NOT. DQ_EXTERNAL(JB) .AND. NSTR(JB) > 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Number of outlet structures [NSTR=',NSTR(JB),'] > 0 for branch ',JB,', but this branch has an internal'//&
                         ' head or flow downstream boundary'
      END IF
    END DO
  END DO

! Pipes

  WRITE (WIN,*) '    pipes'
  DO JP=1,NPI
    IF (JBUPI(JP) == 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream segment [IUPI=',IUPI(JP),'] is located in an inactive boundary segment for pipe ',JP
    END IF
    IF (JBUPI(JP) == JBDPI(JP)) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Upstream [IUPI=',IUPI(JP),'] and downstream [IDPI=',IDPI(JP),'] segments are in same branch for pipe ',JP
    END IF
    IF (IUPI(JP) > 0 .AND. IUPI(JP) < IMX) THEN
      IF (EUPI(JP) > EL(2,IUPI(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Upstream elevation [EUPI=',EUPI(JP),'] > top elevation of the grid [EL=',EL(1,IUPI(JP)),'] for pipe ',JP
      END IF
      IF (EUPI(JP) < EL(KB(IUPI(JP)+1),IUPI(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Upstream elevation [EUPI=',EUPI(JP),'] < elevation of the bottom active cell for pipe ',JP
      END IF
      IF (PUPIC(JP) == ' SPECIFY') THEN
        IF (ETUPI(JP) > EL(2,IUPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUPI=',ETUPI(JP),'] > top elevation of the grid for pipe ',JP
        END IF
        IF (ETUPI(JP) < EL(KB(IUPI(JP)+1),IUPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUPI=',ETUPI(JP),'] < elevation of the bottom active cell'//&
                            ' for pipe ',JP
        END IF
        IF (KBUPI(JP) > KB(IUPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Upstream selective withdrawal bottom layer [KBUPI=',KBUPI(JP),'] > segment bottom active layer [KB=',  &
                            KB(IUPI(JP)),'] for pipe ',JP
        END IF
      END IF
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream segment [IUPI=',IUPI(JP),'] is outside of computational grid for pipe ',JP
    END IF
    IF (WPI(JP) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Pipe diameter [WPI=',WPI(JP),'] <= 0 for pipe ',JP
    END IF
    IF (WPI(JP) > 5.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Pipe diameter [WPI=',WPI(JP),'] > 5 for pipe ',JP
    END IF
    IF (DLXPI(JP) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Pipe length [DLXPI=',DLXPI(JP),'] <= 0 for pipe ',JP
    END IF
    IF (DLXPI(JP) > 1000.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Pipe length [DLXPI=',DLXPI(JP),'] > 1000 m for pipe ',JP
    END IF
    IF (FPI(JP) <= 0.0 .OR. FPI(JP) > 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Friction factor [FPI=',FPI(JP),'] < 0 or > 1 for pipe ',JP
    END IF
    IF (FMINPI(JP) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Minor friction [FMINPI=',FMINPI(JP),'] <  0 for pipe ',JP
    END IF
    IF (FMINPI(JP) > 10.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Minor friction [FMINPI=',FMINPI(JP),'] > 10 for pipe ',JP
    END IF
    IF (PUPIC(JP) /= ' SPECIFY' .AND. PUPIC(JP) /= '   DISTR' .AND. PUPIC(JP) /= ' DENSITY') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream inflow placement [PUPIC='//PUPIC(JP)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for pipe ',JP
    END IF
    IF (KTUPI(JP) <= 1) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTUPI=',KTUPI(JP),'] < 2 for pipe ',JP
    END IF
    IF (KTUPI(JP) > KBUPI(JP)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal top layer [KTUPI=',KTUPI(JP),'] > bottom selective withdrawal layer [KBUPI=',&
                        KBUPI(JP),'] for pipe ',JP
    END IF
    IF (IDPI(JP) /= 0) THEN
      IF (IDPI(JP) == US(JBDPI(JP)-1) .OR. IDPI(JP) == DS(JBDPI(JP))+1 .OR. IDPI(JP) < 2 .OR. IDPI(JP) >= IMX) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream segment [IDPI=',IDPI(JP),'] is located in an inactive boundary segment for pipe ',JP
      END IF
      IF (KTDPI(JP) < 2) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal top layer [KTDPI=',KTDPI(JP),'] < 2 for pipe ',JP
      END IF
      IF (KTDPI(JP) > KBDPI(JP)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal top layer [KTDPI=',KTDPI(JP),'] > bottom selective withdrawal layer '//  &
                         '[KBDPI=',KBDPI(JP),'] for pipe ',JP
      END IF
      IF (KBDPI(JP) > KB(IDPI(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal bottom layer [KBDPI=',KBDPI(JP),'] = bottom active layer [KBDPI=',       &
                                 KB(IDPI(JP)),'] for pipe ',JP
      END IF
      IF (PDPIC(JP) /= ' SPECIFY' .AND. PDPIC(JP) /= '   DISTR' .AND. PDPIC(JP) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream inflow placement [PDPIC='//PDPIC(JP)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for pipe ',JP
      END IF
      IF (PDPIC(JP) == ' SPECIFY') THEN
        IF (ETDPI(JP) > EL(2,IDPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDPI=',ETDPI(JP),'] > grid top elevation for pipe ',JP
        END IF
        IF (ETDPI(JP) < EL(KB(IDPI(JP)+1),IDPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDPI=',ETDPI(JP),'] < elevation of the bottom active '//  &
                            'cell for pipe ',JP
        END IF
      END IF
    END IF
    IF ((LATPIC(JP) /= '    DOWN') .AND. (LATPIC(JP) /= '     LAT')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Pipe lateral control [LATPIC='//LATPIC(JP)(4:8)//'] must be either "DOWN" or "LAT" for pipe ',JP
    END IF
    IF(LATPIC(JP) == '    DOWN') THEN
      IF(IUPI(JP) /= DS(JBUPI(JP)))THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Since PIPE is specified as LATPIC='//LATPIC(JP)(4:8)//', it must be at end of a branch or change this to LAT for pipe ',JP
      ENDIF
    END IF
  END DO

! Spillways/weirs

  WRITE (WIN,*) '    spillways'
  DO JS=1,NSP
    DO JB=1,NBR
      IF (IUSP(JS) == (US(JB)-1) .OR. IUSP(JS) == (DS(JB)+1) .OR. IUSP(JS) < 2 .OR. IUSP(JS) > IMX) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream segment [IUSP=',IUSP(JS),'] is located in a boundary segment for spillway ',JS
      END IF
      IF ((IDSP(JS) == US(JB)-1 .OR. IDSP(JS) == DS(JB)+1 .OR. IDSP(JS) < 2 .OR. IDSP(JS) >= IMX) .AND. IDSP(JS) /= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream segment [IDSP=',IDSP(JS),'] is located in a boundary segment for spillway ',JS
      END IF
      DO JBB=1,NBR                                                                                                     !SW 05/16/02
        IF (IUSP(JS) >= US(JBB) .AND. IUSP(JS) <= DS(JBB)) JBU = JBB                                                   !SW 05/16/02
        IF (IDSP(JS) >= US(JBB) .AND. IDSP(JS) <= DS(JBB)) JBD = JBB                                                   !SW 05/16/02
      END DO                                                                                                           !SW 05/16/02
      IF (JBD == JBU .AND. JBU /= 0 .AND. IDSP(JS) /= 0) THEN                                                          !SW 05/16/02
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream [IUSP=',IUSP(JS),'] and downstream [IDSP=',IDSP(JS),'] are in same branch for spillway ',JS
      END IF
    END DO
    IF (ESP(JS) > EL(2,IUSP(JS))) THEN
      CALL ERRORS
      WRITE (ERR,FMT2FI) 'Upstream centerline elevation [ESP=',ESP(JS),'] > the top elevation of the grid [EL=',EL(2,IUSP(JS)),    &
                         '] for spillway ',JS
    END IF
    IF (IUSP(JS) > 0 .AND. IUSP(JS) < IMX) THEN
      IF (ESP(JS) < EL(KB(IUSP(JS))+1,IUSP(JS))) THEN                                                                  !TC 09/30/03
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Upstream centerline elevation [ESP=',ESP(JS),'] < the bottom active cell elevation [EL=',              &
                            EL(KB(IUSP(JS)),IUSP(JS)),'] for spillway ',JS                                             !SW 10/16/02
      END IF
    END IF
    IF (KTUSP(JS) < 2) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTUSP=',KTUSP(JS),'] < 2 for spillway',JS
    END IF
    IF (KTUSP(JS) > KBUSP(JS)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal top layer [KTUSP=',KTUSP(JS),'] > bottom selective withdrawal layer [KBUSP=',&
                        KBUSP(JS),'] for spillway ',JS
    END IF
    IF (KBUSP(JS) > KB(IUSP(JS))) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal bottom layer [KBUSP=',KBUSP(JS),'] >= segment bottom active layer [KB=',     &
                        KB(IUSP(JS)),'] for spillway ',JS
    END IF
    IF (IDSP(JS) /= 0) THEN
      IF (KTDSP(JS) < 2) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTDSP=',KTDSP(JS),'] < 2 for spillway ',JS
      END IF
      IF (KTDSP(JS) > KBDSP(JS)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal top layer [KTDSP=',KTDSP(JS),'] > bottom selective withdrawal layer '//  &
                         '[KBUSP=',KBDSP(JS),'] for spillway ',JS
      END IF
      IF (KBDSP(JS) > KB(IDSP(JS))) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal bottom layer [KBDSP=',KBDSP(JS),'] >= segment bottom active layer [KB=', &
                          KB(IDSP(JS)),'] for spillway ',JS
      END IF
      IF (PUSPC(JS) /= ' SPECIFY' .AND. PUSPC(JS) /= '   DISTR' .AND. PUSPC(JS) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream inflow placement [PUSPC='//PUSPC(JS)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for spillway ',JS
      END IF
      IF (PDSPC(JS) /= ' SPECIFY' .AND. PDSPC(JS) /= '   DISTR' .AND. PDSPC(JS) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream inflow placement [PDSPC='//PDSPC(JS)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for spillway ',JS
      END IF
      IF (PUSPC(JS) == ' SPECIFY') THEN
        IF (ETUSP(JS) > EL(2,IUSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUSP=',ETUSP(JS),'] > grid top elevationfor spillway ',JS
        END IF
        IF (ETUSP(JS) < EL(KB(IUSP(JS)+1),IUSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUSP=',ETUSP(JS),'] < elevation of the bottom active '//    &
                            'cell for spillway ',JS
        END IF
        IF (ETUSP(JS) < ESP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUSP=',ETUSP(JS),'] < spillway ',JS,' elevation'
        END IF
        IF (EBUSP(JS) > EL(2,IUSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUSP=',EBUSP(JS),'] > spillway ',JS,' top elevation'
        END IF
        IF (EBUSP(JS) < EL(KB(IUSP(JS)+1),IUSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUSP=',EBUSP(JS),'] < elevation of the bottom active '// &
                            'cell for spillway ',JS
        END IF
        IF (EBUSP(JS) > ESP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUSP=',EBUSP(JS),'] > elevation of the spillway for '//  &
                            'spillway ',JS
        END IF
        IF (ETUSP(JS) < EBUSP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUSP=',ETUSP(JS),'] < elevation of the bottom inflow '//    &
                            'placement for spillway ',JS
        END IF
      END IF
      IF (PDSPC(JS) == ' SPECIFY') THEN
        IF (ETDSP(JS) > EL(2,IDSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDSP=',ETDSP(JS),'] > grid top elevation for spillway ',JS
        END IF
        IF (ETDSP(JS) < EL(KB(IDSP(JS)+1),IDSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDSP=',ETDSP(JS),'] < bottom active cell elevation for '//&
                            'spillway ',JS
        END IF
        IF (ETDSP(JS) < ESP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDSP=',ETDSP(JS),'] < spillway elevation for spillway ',JS
        END IF
        IF (EBDSP(JS) > EL(2,IDSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDSP=',EBDSP(JS),'] > spillway top elevation for '//   &
                            'spillway ',JS
        END IF
        IF (EBDSP(JS) < EL(KB(IDSP(JS)+1),IDSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDSP=',EBDSP(JS),'] < bottom active cell elevation '// &
                            'for spillway ',JS
        END IF
        IF (EBDSP(JS) > ESP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDSP=',EBDSP(JS),'] > spillway elevation for '//       &
                            'spillway ',JS
        END IF
        IF (ETDSP(JS) < EBDSP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDSP=',ETDSP(JS),'] < bottom inflow placement elevation'//&
                            ' for spillway ',JS
        END IF
      END IF
    END IF
    IF ((LATSPC(JS) /= '    DOWN') .AND. (LATSPC(JS) /= '     LAT')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway lateral control [LATSPC='//LATSPC(JS)(4:8)//'] must be either "DOWN" or "LAT" for spillway ',JS
    END IF
    IF (LATSPC(JS) == '    DOWN') THEN     ! SW 3/26/10
      iflag=0
      DO JB=1,NBR
        IF (IUSP(JS) == DS(JB)) THEN
          iflag=1
          exit
        END IF
      END DO
        IF(IFLAG == 0)THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Since SPILLWAY is specified as LATSPC='//LATSPC(JS)(4:8)//', it must be at end of a branch or change this to LAT for spillway ',JS
        ENDIF
    END IF

!** Coefficients

    IF (A1SP(JS) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway coefficient a1 [A1SP] <= 0 for spillway ',JS
    END IF
    IF (A2SP(JS) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway coefficient a2 [A2SP] < 0 for spillway ',JS
    END IF
    IF (B1SP(JS) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway coefficient b1 [B1SP] < 0 for spillway ',JS
    END IF
    IF (B2SP(JS) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway coefficient b2 [B2SP] < 0 for spillway ',JS
    END IF
    IF (B1SP(JS) == 0.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Spillway coefficient b1 [B1SP] = 0 for spillway ',JS
    END IF
    IF (B1SP(JS) > 2.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Spillway coefficient b1 [B1SP] > 2 for spillway ',JS
    END IF

!** Total dissolved gas equations

    IF ((GASSPC(JS) /= '      ON') .AND. (GASSPC(JS) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway gas control [GASSPC='//GASSPC(JS)(6:8)//'] must be either " ON" or "OFF" for spillway ',JS
    END IF
    IF (IDSP(JS) /= 0 .AND. GASSPC(JS) == '      ON') THEN
      IF (EQSP(JS) == 1) THEN
        IF (AGASSP(JS) < 0.0 .OR. AGASSP(JS) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient a [AGASSP] < 0 or > 10 for spillway ',JS
        END IF
        IF (BGASSP(JS) < 100.0 .OR. BGASSP(JS) > 150.) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient b [BGASSP] < 100 or > 150 for spillway ',JS
        END IF
      ELSE IF (EQSP(JS) == 2) THEN
        IF (AGASSP(JS) < 100.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient a [AGASSP] < 100 for spillway ',JS
        END IF
        IF (BGASSP(JS) > 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient b [BGASSP] > 0 for spillway ',JS
        END IF
        IF (CGASSP(JS) > 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient c [CGASSP] > 0 for spillway ',JS
        END IF
      ELSE IF (EQSP(JS) == 3) THEN
        IF (AGASSP(JS) < 0.65 .OR. AGASSP(JS) > 1.85) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient a [AGASSP] < 0.65 or > 1.85 for spillway ',JS
        END IF
        IF (BGASSP(JS) > 1.0 .OR. BGASSP(JS) < 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient b [BGASSP] > 1 or < 0.05 for spillway ',JS
        END IF
        IF (CGASSP(JS) < 0.0 .OR. CGASSP(JS) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient c [CGASSP] < 0 or > 10 for spillway ',JS
        END IF
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTI) 'Equation number [EQSP]=',EQSP(JS),'] must be between 1 and 3 for spillway ',JS
      END IF
    END IF
  END DO

! Gates

  WRITE (WIN,*) '    gates'
  DO JT=1,NGT
    DO JB=1,NBR
      IF (IUGT(JT) == (US(JB)-1) .OR. IUGT(JT) == (DS(JB)+1) .OR. IUGT(JT) < 2 .OR. IUGT(JT) > IMX) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream segment [IUGT=',IUGT(JT),'] is located in a boundary segment for gate ',JT
      END IF
      IF ((IDGT(JT) == US(JB)-1 .OR. IDGT(JT) == DS(JB)+1 .OR. IDGT(JT) < 2 .OR. IDGT(JT) >= IMX) .AND. IDGT(JT) /= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream segment [IDGT=',IDGT(JT),'] is located in a boundary segment for gate ',JT
      END IF
      DO JBB=1,NBR                                                                                                     !SW 05/16/02
        IF (IUGT(JT) >= US(JBB) .AND. IUGT(JT) <= DS(JBB)) JBU = JBB                                                   !SW 05/16/02
        IF (IDGT(JT) >= US(JBB) .AND. IDGT(JT) <= DS(JBB)) JBD = JBB                                                   !SW 05/16/02
      END DO                                                                                                           !SW 05/16/02
      IF (JBD == JBU .AND. JBU /= 0 .AND. IDGT(JT) /= 0) THEN                                                          !SW 05/16/02
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream [IUGT=',IUGT(JT),'] and downstream [IDGT=',IDGT(JT),'] are in same branch for gate ',JT
      END IF
    END DO
    IF (EGT(JT) > EL(2,IUGT(JT))) THEN
      CALL ERRORS
      WRITE (ERR,FMT2FI) 'Upstream centerline elevation [EGT=',EGT(JT),'] > the top elevation of the grid [EL=',EL(2,IUGT(JT)),    &
                         '] for gate ',JT
    END IF
    IF (IUGT(JT) > 0 .AND. IUGT(JT) < IMX) THEN
      IF (EGT(JT) < EL(KB(IUGT(JT))+1,IUGT(JT))) THEN                                                                  !TC 09/30/03
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Upstream centerline elevation [EGT=',EGT(JT),'] < the bottom active cell elevation [EL=',              &
                            EL(KB(IUGT(JT)),IUGT(JT)),'] for gate ',JT                                             !SW 10/16/02
      END IF
    END IF
    IF (KTUGT(JT) < 2) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTUGT=',KTUGT(JT),'] < 2 for gate',JT
    END IF
    IF (KTUGT(JT) > KBUGT(JT)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal top layer [KTUGT=',KTUGT(JT),'] > bottom selective withdrawal layer [KBUGT=',&
                        KBUGT(JT),'] for gate ',JT
    END IF
    IF (KBUGT(JT) > KB(IUGT(JT))) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal bottom layer [KBUGT=',KBUGT(JT),'] >= segment bottom active layer [KB=',     &
                        KB(IUGT(JT)),'] for gate ',JT
    END IF
    IF (IDGT(JT) /= 0) THEN
      IF (KTDGT(JT) < 2) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTDGT=',KTDGT(JT),'] < 2 for gate ',JT
      END IF
      IF (KTDGT(JT) > KBDGT(JT)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal top layer [KTDGT=',KTDGT(JT),'] > bottom selective withdrawal layer '//  &
                         '[KBUGT=',KBDGT(JT),'] for gate ',JT
      END IF
      IF (KBDGT(JT) > KB(IDGT(JT))) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal bottom layer [KBDGT=',KBDGT(JT),'] >= segment bottom active layer [KB=', &
                          KB(IDGT(JT)),'] for gate ',JT
      END IF
      IF (PUGTC(JT) /= ' SPECIFY' .AND. PUGTC(JT) /= '   DISTR' .AND. PUGTC(JT) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream inflow placement [PUGTC='//PUGTC(JT)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for gate ',JT
      END IF
      IF (PDGTC(JT) /= ' SPECIFY' .AND. PDGTC(JT) /= '   DISTR' .AND. PDGTC(JT) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream inflow placement [PDGTC='//PDGTC(JT)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for gate ',JT
      END IF
      IF (PUGTC(JT) == ' SPECIFY') THEN
        IF (ETUGT(JT) > EL(2,IUGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUGT=',ETUGT(JT),'] > grid top elevation for gate ',JT
        END IF
        IF (ETUGT(JT) < EL(KB(IUGT(JT)+1),IUGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUGT=',ETUGT(JT),'] < elevation of the bottom active '//    &
                            'cell for gate  ',JT
        END IF
        IF (ETUGT(JT) < EGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUGT=',ETUGT(JT),'] < gate ',JT,' elevation'
        END IF
        IF (EBUGT(JT) > EL(2,IUGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUGT=',EBUGT(JT),'] > gate ',JT,' top elevation'
        END IF
        IF (EBUGT(JT) < EL(KB(IUGT(JT)+1),IUGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUGT=',EBUGT(JT),'] < elevation of the bottom active '// &
                            'cell for gate ',JT
        END IF
        IF (EBUGT(JT) > EGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUGT=',EBUGT(JT),'] > elevation of the gate for '//  &
                            'gate ',JT
        END IF
        IF (ETUGT(JT) < EBUGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUGT=',ETUGT(JT),'] < elevation of the bottom inflow '//    &
                            'placement for gate ',JT
        END IF
      END IF
      IF (PDGTC(JT) == ' SPECIFY') THEN
        IF (ETDGT(JT) > EL(2,IDGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDGT=',ETDGT(JT),'] > grid top elevation for gate ',JT
        END IF
        IF (ETDGT(JT) < EL(KB(IDGT(JT)+1),IDGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDGT=',ETDGT(JT),'] < bottom active cell elevation for '//&
                            'gate ',JT
        END IF
        IF (ETDGT(JT) < EGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDGT=',ETDGT(JT),'] < gate elevation for gate ',JT
        END IF
        IF (EBDGT(JT) > EL(2,IDGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDGT=',EBDGT(JT),'] > gate top elevation for '//   &
                            'gate ',JT
        END IF
        IF (EBDGT(JT) < EL(KB(IDGT(JT)+1),IDGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDGT=',EBDGT(JT),'] < bottom active cell elevation '// &
                            'for gate ',JT
        END IF
        IF (EBDGT(JT) > EGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDGT=',EBDGT(JT),'] > gate elevation for '//       &
                            'gate ',JT
        END IF
        IF (ETDGT(JT) < EBDGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDGT=',ETDGT(JT),'] < bottom inflow placement elevation'//&
                            ' for gate ',JT
        END IF
      END IF
    END IF


    IF ((LATGTC(JT) /= '    DOWN') .AND. (LATGTC(JT) /= '     LAT')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate lateral control [LATGTC='//LATGTC(JT)(4:8)//'] must be either "DOWN" or "LAT" for gate ',JT
    END IF
    IF ((DYNGTC(JT) /= '    FLOW') .AND. (DYNGTC(JT) /= '     ZGT').AND. (DYNGTC(JT) /= '       B')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate dynamic variable control [DYNGTC='//DYNGTC(JT)(4:8)//'] must be either "FLOW" or "ZGT" or "B" for gate ',JT
    END IF
    IF (LATGTC(JT) == '    DOWN') THEN                   ! SW 3/26/10
      IF(IUGT(JT) /= DS(JBUGT(JT)))THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Since GATE is specified as LATGTC='//LATGTC(JT)(4:8)//', it must be at end of a branch or change this to LAT for gate ',JT
      ENDIF
    END IF


!** Coefficients

    IF (A1GT(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient a1 [A1GT] <= 0 for gate ',JT
    END IF
    IF (A2GT(JT) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient a2 [A2GT] < 0 for gate ',JT
    END IF
    IF (B1GT(JT) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient b1 [B1GT] < 0 for gate ',JT
    END IF
    IF (B2GT(JT) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient b2 [B2GT] < 0 for gate ',JT
    END IF
    IF (B1GT(JT) == 0.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Gate coefficient b1 [B1GT] = 0 for gate ',JT
    END IF
    IF (B1GT(JT) > 2.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Gate coefficient b1 [B1GT] > 2 for gate ',JT
    END IF
       IF (G1GT(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient g1 [G1GT] <= 0 for gate ',JT
    END IF
       IF (G2GT(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient g2 [G2GT] <= 0 for gate ',JT
    END IF
     IF (GTA1(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient gta1 [GTA1] <= 0 for gate ',JT
    END IF
     IF (GTA2(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient gta2 [GTA2] <= 0 for gate ',JT
    END IF
      IF (GTB1(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient gtb1 [GTB1] <= 0 for gate ',JT
    END IF
     IF (GTB2(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient gtb2 [GTB2] <= 0 for gate ',JT
    END IF

!** Total dissolved gas equations

    IF ((GASGTC(JT) /= '      ON') .AND. (GASGTC(JT) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate gas control [GASGTC='//GASGTC(JT)(6:8)//'] must be either " ON" or "OFF" for gate ',JT
    END IF
    IF (IDGT(JT) /= 0 .AND. GASGTC(JT) == '      ON') THEN
      IF (EQGT(JT) == 1) THEN
        IF (AGASGT(JT) < 0.0 .OR. AGASGT(JT) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient a [AGASGT] < 0 or > 10 for gate ',JT
        END IF
        IF (BGASGT(JT) < 100.0 .OR. BGASGT(JT) > 150.) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient b [BGASGT] < 100 or > 150 for gate ',JT
        END IF
      ELSE IF (EQGT(JT) == 2) THEN
        IF (AGASGT(JT) < 100.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient a [AGASGT] < 100 for gate ',JT
        END IF
        IF (BGASGT(JT) > 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient b [BGASGT] > 0 for gate ',JT
        END IF
        IF (CGASGT(JT) > 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient c [CGASGT] > 0 for gate ',JT
        END IF
      ELSE IF (EQGT(JT) == 3) THEN
        IF (AGASGT(JT) < 0.65 .OR. AGASGT(JT) > 1.85) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient a [AGASGT] < 0.65 or > 1.85 for gate ',JT
        END IF
        IF (BGASGT(JT) > 1.0 .OR. BGASGT(JT) < 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient b [BGASGT] > 1 or < 0.05 for gate ',JT
        END IF
        IF (CGASGT(JT) < 0.0 .OR. CGASGT(JT) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient c [CGASGT] < 0 or > 10 for gate ',JT
        END IF
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTI) 'Equation number [EQGT]=',EQGT(JT),'] must be between 1 and 3 for gate ',JT
      END IF
    END IF
  END DO




! Pumps

  WRITE (WIN,*) '    pumps'
  DO JP=1,NPU
    DO JB=1,NBR
      IF (IUPU(JP) == US(JB)-1 .OR. IUPU(JP) == DS(JB)+1 .OR. IUPU(JP) <= 0 .OR. IUPU(JP) >= IMX) THEN                 ! RA 8/14/07
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream pump segment [IUPU=',IUPU(JP),'] is located in a boundary segment for pump ',JP
      END IF
      IF (IDPU(JP) == US(JB)-1 .OR. IDPU(JP) == DS(JB)+1 .OR. IDPU(JP) < 0 .OR. IDPU(JP) >= IMX) THEN                  !TC 09/30/03
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream pump segment [IDPU=',IDPU(JP),'] is located in a boundary segment for pump ',JP
      END IF
    END DO
    IF (EPU(JP) > EL(2,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Upstream pump elevation [EPU=',EPU(JP),'] > layer 2 for pump ',JP
    END IF
    IF (EPU(JP) < EL(KB(IUPU(JP))+1,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'Upstream pump elevation [EPU=',EPU(JP),'] < bottom active layer [KB=',KB(IUPU(JP)+1),'] for pump ',JP
    END IF
    IF (EONPU(JP) > EL(2,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Upstream elevation for turning pump on [EONPU=',EONPU(JP),'] < layer 2 for pump ',JP
    END IF
    IF (EONPU(JP) < EL(KB(IUPU(JP))+1,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'Upstream elevation for turning pump on [EONPU=',EONPU(JP),'] < bottom active layer [KB=',KB(IUPU(JP)),   &
                         '] for pump ',JP
    END IF
    IF (EONPU(JP) <= EOFFPU(JP)) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Upstream elevation for turning pump on [EONPU=',EONPU(JP),'] < the elevation for turning the pump off '// &
                        'for pump ',JP
    END IF
    IF (EOFFPU(JP) > EL(2,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Upstream elevation for turning pump off [EOFFPU=',EOFFPU(JP),'] < grid top elevation for pump ',JP
    END IF
    IF (EOFFPU(JP) < EL(KB(IUPU(JP))+1,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'Upstream elevation for turning pump off [EOFFPU=',EOFFPU(JP),'] > bottom active layer [KB=',KB(IUPU(JP)),&
                         '] for pump ',JP
    END IF
    IF (KTPU(JP) < 2) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Top selective withdrawal layer [KTPU=',KTPU,'] < layer 2 for pump ',JP
    END IF
    IF (KBPU(JP) > KB(IUPU(JP))+1) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Bottom selective withdrawal layer [KBPU=',KBPU,'] > bottom active layer [KB=',KB(IUPU(JP)),'] for pump ',JP
    END IF
    IF (KTPU(JP) > KBPU(JP)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Top selective withdrawal layer [KTPU=',KTPU,'] > bottom selective withdrawal layer [KBPU=',KBPU(JP),       &
                       '] for pump ',JP
    END IF
    IF (EPU(JP) > EONPU(JP)) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Centerline elevation [EPU=',EPU(JP),'] > on elevation for pump ',JP
    END IF
    IF (EPU(JP) > EOFFPU(JP)) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Centerline elevation [EPU=',EPU(JP),'] > off elevation for pump ',JP
    END IF
    IF (PPUC(JP) == ' SPECIFY' .AND. IDPU(JP) /= 0) THEN                                                               !TC 09/30/03
      IF (PPUC(JP) /= ' SPECIFY' .AND. PPUC(JP) /= '   DISTR' .AND. PPUC(JP) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Inflow placement [PPPUC='//PPUC(JP)//'] /= " SPECIFY" or "   DISTR" or " DENSITY" for pump ',JP
      END IF
      IF (ETPU(JP) > EL(2,IDPU(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETPU=',ETPU(JP),'] > grid top elevation for pump ',JP
      END IF
      IF (ETPU(JP) < EL(KB(IDPU(JP))+1,IDPU(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETPU=',ETPU(JP),'] < bottom active cell elevation for pump ',&
                           JP
      END IF
      IF (EBPU(JP) > EL(2,IDPU(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBPU=',EBPU(JP),'] > top elevation of pump ',JP
      END IF
      IF (EBPU(JP) < EL(KB(IDPU(JP))+1,IDPU(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBPU=',EBPU(JP),'] < bottom active cell elevation for '// &
                          'pump ',JP
      END IF
      IF (EBPU(JP) > EPU(JP)) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBPU=',EBPU(JP),'] > elevation of pump ',JP
      END IF
      IF (ETPU(JP) < EBPU(JP)) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETPU=',ETPU(JP),'] < inflow placement bottom elevation for'//&
                          ' pump ',JP
      END IF
    END IF
    IF ((LATPUC(JP) /= '    DOWN') .AND. (LATPUC(JP) /= '     LAT')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Pump lateral control [LATPUC='//LATPUC(JP)(4:8)//'] must be either "DOWN" or "LAT" for pump ',JP
    END IF

     IF (LATPUC(JP) == '    DOWN') THEN     ! SW 3/26/10
      iflag=0
      DO JB=1,NBR
        IF (IUPU(JP) == DS(JB)) THEN
          iflag=1
          exit
        END IF
      END DO
        IF(IFLAG == 0)THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Since PUMP is specified as LATPUC='//LATPUC(JP)(4:8)//', it must be at end of a branch or change this to LAT for pump ',JP
        ENDIF
    END IF


  END DO

! Internal weir

  WRITE (WIN,*) '    internal weir'
  DO JWR=1,NIW
    IFLAG = 0
    DO JW=1,NWB
      DO JB=BS(JW),BE(JW)
        IF (IWR(JWR) >= US(JB) .AND. IWR(JWR) <= DS(JB)) THEN
          KT    = KTWB(JW)
          IFLAG = 1
          EXIT
        END IF
      END DO
      IF (IFLAG == 1) EXIT
    END DO
    IF (KTWR(JWR) <= KT) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Internal weir surface layer [KTWR=',KTWR(JWR),'] < water surface layer [KT=',KT,'] for weir ',JWR
    END IF
    IF (KTWR(JWR) > KBWR(JWR)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Internal weir top layer [KTWR=',KTWR(JWR),'] > weir bottom layer [KBWR=',KBWR(JWR), '] for weir ',JWR
    END IF
    IF (KBWR(JWR) > KB(IWR(JWR))+1) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Internal weir bottom layer [KBWR=',KBWR(JWR),'] > bottom active layer [KB=',KB(IWR(JWR)),'] for weir ',JWR
    END IF
    DO JB=1,NBR
      IF (IWR(JWR) == US(JB)-1) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Internal weir segment [IWR=',IWR(JWR),'] is an upstream boundary segment in branch ',JB
      ELSE IF (IWR(JWR) == DS(JB)+1) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Internal weir segment [IWR=',IWR(JWR),'] is a downstream boundary segment in branch ',JB
      ELSE IF (IWR(JWR) == DS(JB)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Internal weir segment [IWR=',IWR(JWR),'] is a downstream segment in branch ',JB
      END IF
    END DO
  END DO

! Lateral withdrawals

  WRITE (WIN,*) '    lateral withdrawals'
  DO JW=1,NWD
    IF (EWD(JW) > EL(2,IWD(JW))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Withdrawal elevation [EWD=',EWD(JW),'] > grid top elevation for withdrawal ',JW
    ELSE IF (EWD(JW) < EL(KB(IWD(JW)+1),IWD(JW))) THEN
      CALL ERRORS
      WRITE (ERR,FMT2FI) 'Withdrawal elevation [EWD=',EWD(JW),'] < bottom active layer elevation [EL=',EL(KB(IWD(JW)+1),IWD(JW)),  &
                          '] for withdrawal ',JW
    END IF
    DO JB=1,NBR
      IF (IWD(JW) == US(JB)-1 .OR. (JB == 1 .AND. IWD(JW) == 0)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Withdrawal segment [IWD=',IWD(JW),'] is an upstream boundary segment in branch ',JB
      ELSE IF (IWD(JW) == DS(JB)+1) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Withdrawal segment [IWD=',IWD(JW),'] is a downstream boundary segment in branch ',JB
      END IF
    END DO
    IF (KTWD(JW) < 2) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Top selective withdrawal layer [KTWD=',KTWD(JW),'] < layer 2 for withdrawal ',JW
    END IF
    IF (KBWD(JW) > KB(IWD(JW))+1) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Bottom selective withdrawal layer [KBWD=',KBWD(JW),'] > bottom active layer [KB=',KB(IWD(JW)),             &
                       '] for withdrawal ',JW
    END IF
    IF (KTWD(JW) > KBWD(JW)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Top selective withdrawal layer [KTWD=',KTWD(JW),'] > bottom selective withdrawal layer [KBWD=',KBWD(JW),   &
                       '] for withdrawal ', JW
    END IF
  END DO

! Tributaries

  WRITE (WIN,*) '    tributaries'
  IF (TRIBUTARIES) THEN
    DO JT=1,NTR
      IF (PTRC(JT) /= ' SPECIFY' .AND. PTRC(JT) /= '   DISTR' .AND. PTRC(JT) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Inflow placement [PTRC='//PTRC(JT)//'] /= " SPECIFY" or "   DISTR" or " DENSITY" for tributary ',JP
      END IF
      ACTIVE_SEGMENT = .FALSE.
      DO JB=1,NBR
        IF (ITR(JT) >= US(JB) .AND. ITR(JT) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
      END DO
      IF (.NOT.ACTIVE_SEGMENT) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Tributary segment [ITR(',JT,')=',ITR(JT),'] is a boundary segment for tributary ',JT
      END IF
      IF (PTRC(JT) == ' SPECIFY') THEN
        IF (ETTR(JT) > EL(2,ITR(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement top elevation [ETTR=',ETTR(JT),'] > grid top elevation for tributary ',JT
        END IF
        IF (ETTR(JT) < EL(KB(ITR(JT))+1,ITR(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement top elevation [ETTR=',ETTR(JT),'] < bottom active cell elevation for tributary ',JT
        END IF
        IF (EBTR(JT) > EL(2,ITR(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement bottom elevation [EBTR=',EBTR(JT),'] > grid top elevation for tributary ',JT
        END IF
        IF (EBTR(JT) < EL(KB(ITR(JT))+1,ITR(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement bottom elevation [EBTR=',EBTR(JT),'] < bottom active cell elevation for tributary ',JT
        END IF
        IF (ETTR(JT) < EBTR(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement top elevation [ETTR=',ETTR(JT),'] < bottom elevation inflow placement for '//        &
                            'tributary ',JT
        END IF
      END IF
    END DO
  END IF

! Distributed tributaries

  WRITE (WIN,*) '    distributed tributaries'
  DO JB=1,NBR
    IF ((DTRC(JB) /= '      ON') .AND. (DTRC(JB) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Distributed tributary control [DTRC='//DTRC(JB)(6:8)//'] must be either " ON" or "OFF" for branch ',JB
    END IF
  END DO

! Hydrodynamic output

  WRITE (WIN,*) 'Output controls'
  WRITE (WIN,*) '  hydrodynamic variables'
  DO JW=1,NWB
    DO JH=1,NHY
      IF ((HPRWBC(JH,JW) /= '      ON') .AND. (HPRWBC(JH,JW) /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Hydrodynamic output control [HPRWBC='//HPRWBC(JH,JW)(6:8)//'] must be either " ON" or "OFF" for '//      &
                         'variable ',JH,' in waterbody ',JW
      END IF
    END DO
  END DO

! Snapshot output

  WRITE (WIN,*) '  snapshot'
  DO JW=1,NWB
    IF ((SNPC(JW) /= '      ON') .AND. (SNPC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Snapshot control [SNPC='//SNPC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (SNPC(JW) == '      ON') THEN
      IF (NSNP(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Number of snapshot dates [NSNP=',NSNP(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NSNP(JW)-1
        IF (SNPD(J+1,JW) <= SNPD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'Snapshot date [SNPD=',SNPD(J+1,JW),'] <= previous date [SNPD=',SNPD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NSNP(JW)
        IF (SNPF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Snapshot frequency [SNPF=',SNPF(J,JW),'] <= 0 for date [SNPD=',SNPD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      IF (NISNP(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Number of snapshot segment [NISNP=',NISNP(JW),'] s <= 0 for waterbody ',JW
      END IF
      ACTIVE_SEGMENT = .FALSE.
      DO J=1,NISNP(JW)
        DO JB=1,NBR
          IF (ISNP(J,JW) >= US(JB) .AND. ISNP(J,JW) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
        END DO
        IF (.NOT.ACTIVE_SEGMENT) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Snapshot segment [ISNP=',ISNP(J,JW),'] is a boundary segment in waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Screen output

  WRITE (WIN,*) '  screen'
  DO JW=1,NWB
    IF ((SCRC(JW) /= '      ON') .AND. (SCRC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Screen control [SCRC='//SCRC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (SCRC(JW) == '      ON') THEN
      IF (NSCR(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Number of screen dates [NSCR=',NSCR(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NSCR(JW)-1
        IF (SCRD(J+1,JW) <= SCRD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'Screen date [SCRD=',SCRD(J+1,JW),'] <= previous date [SCRD=',SCRD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NSCR(JW)
        IF (SCRF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Screen frequency [SCRF=',SCRF(J,JW),'] <= 0 for waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Time series output

  WRITE (WIN,*) '  time series'
  IF ((TSRC /= '      ON') .AND. (TSRC /= '     OFF')) THEN
    CALL ERRORS
    WRITE (ERR,FMTA) 'Time series control [TSRC='//TSRC(6:8),'] must be either " ON" or "OFF"'
  END IF
  IF (TSRC == '      ON') THEN
    IF (NTSR <= 0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Number of time series dates [NTSR=',NTSR,' <= 0'
    END IF
    DO J=1,NTSR-1
      IF (TSRD(J+1) <= TSRD(J)) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Time series date [TSRD=',TSRD(J+1),'] <= previous date [TSRD=',TSRD(J),']'
      END IF
    END DO
    DO J=1,NTSR
      IF (TSRF(J) <= 0.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTF) 'Time series frequency [TSRF=',TSRF(J),'] <= 0 for date [TSRD=',TSRD(J),']'
      END IF
    END DO
    DO JB=1,NBR                 ! SW 8/12/05
    DO J=1,NTSR                 ! SW 7/5/05
      IF (ITSR(J) == (US(JB)-1) .OR. ITSR(J) == (DS(JB)+1) .OR. ITSR(J) < 2 .OR. ITSR(J) > IMX) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Time series segment[ITSR=',ITSR(J),'] is located at a boundary segment for NTSR ',J
      END IF
    END DO
    END DO
  END IF

! Spreadsheet output

  WRITE (WIN,*) '  spreadsheet'
  DO JW=1,NWB
    IF ((SPRC(JW) /= '      ON') .AND. (SPRC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spreadsheet control [SPRC='//SPRC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (SPRC(JW) == '      ON') THEN
      IF (NSPR(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Number of spreadsheet dates [NSPR=',NSPR(JW),']  <= 0 for waterbody ',JW
      END IF
      DO J=1,NSPR(JW)-1
        IF (SPRD(J+1,JW) <= SPRD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'Spreadsheet date [SPRD=',SPRD(J+1,JW),'] <= previous date [SPRD=', SPRD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NSPR(JW)
        IF (SPRF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Spreadsheet frequency [SPRF=',SPRF(J,JW),'] <= 0 for date [SPRD=',SPRD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      IF (NISPR(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Number of spreadsheet segments [NISPR=',NISPR(JW),'] <= 0 for waterbody ',JW
      END IF
      ACTIVE_SEGMENT = .FALSE.
      DO J=1,NISPR(JW)
        DO JB=1,NBR
          IF (ISPR(J,JW) >= US(JB) .AND. ISPR(J,JW) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
        END DO
        IF (.NOT.ACTIVE_SEGMENT) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spreadsheet segment [ISPR=',ISPR(J,JW),'] is a boundary segment in waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Profile output

  WRITE (WIN,*) '  profile plots'
  DO JW=1,NWB
    IF ((PRFC(JW) /= '      ON') .AND. (PRFC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Profile control [PRFC='//PRFC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (PRFC(JW) == '      ON') THEN
      IF (NPRF(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Number of profile dates [NPRF=',NPRF(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NPRF(JW)-1
        IF (PRFD(J+1,JW) <= PRFD(J,JW)) THEN
          CALL ERRORS
        WRITE (ERR,FMT2FI) 'Profile date [PRFD=',PRFD(J+1,JW),'] <= previous date [PRFD=',PRFD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NPRF(JW)
        IF (PRFF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Profile frequency [PRFF=',PRFF(J,JW),'] <= 0 for date [PRFD=',PRFD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      IF (NIPRF(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Number of profile segments [NIPRF=',NIPRF(JW),'] <= 0 for waterbody ',JW
      END IF
      ACTIVE_SEGMENT = .FALSE.
      DO J=1,NIPRF(JW)
        DO JB=1,NBR
          IF (IPRF(J,JW) >= US(JB) .AND. IPRF(J,JW) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
        END DO
        IF (.NOT.ACTIVE_SEGMENT) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Profile segment [IPRF=',IPRF(J,JW),'] is a boundary segment in waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Vector output

  WRITE (WIN,*) '  vector plots'
  DO JW=1,NWB
    IF ((VPLC(JW) /= '      ON') .AND. (VPLC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Vector plot control [VPLC='//VPLC(JW)(6:8),'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (VPLC(JW) == '      ON') THEN
      IF (NVPL(JW) <= 0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTI) 'Number of vector plot dates [NVPL=',NVPL(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NVPL(JW)-1
        IF (VPLD(J+1,JW) <= VPLD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'Vector plot date [VPLD=',VPLD(J,JW),'] <= previous date [VPLD=',VPLD(J-1,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NVPL(JW)
        IF (VPLF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Vector plot frequency [VPLF=',VPLF(J,JW),'] <= 0 for date [VPLD=',VPLD(J,JW),'] in waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Contour output

  WRITE (WIN,*) '  contour plots'
  DO JW=1,NWB
    IF ((CPLC(JW) /= '      ON') .AND. (CPLC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Contour plot control [VPLC='//VPLC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF ((TCPL(JW) /= '      ON') .AND. (TCPL(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Contour plot control [TECPLOT='//TCPL(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (CPLC(JW) == '      ON') THEN
      IF (NCPL(JW) <= 0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTI) 'Number of contour plot dates [NCPL=',NCPL(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NCPL(JW)-1
        IF (CPLD(J+1,JW) <= CPLD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'Contour plot date [VPLD=',VPLD(J,JW),'] <= previous date [VPLD=',VPLD(J-1,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NCPL(JW)
        IF (CPLF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Contour plot frequency [VPLF=',VPLF(J,JW),'] <= 0 for date [VPLD=',VPLD(J,JW),'] in waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Withdrawal output

  WRITE (WIN,*) '  withdrawal'
  IF ((WDOC /= '      ON') .AND. (WDOC /= '     OFF')) THEN
    CALL ERRORS
    WRITE (ERR,FMTA) 'Withdrawal outflow control [WDOC='//WDOC(6:8)//'] must be either " ON" or "OFF"'
  END IF
  IF (WDOC == '      ON') THEN
    DO J=1,NIWDO
      ACTIVE_SEGMENT = .FALSE.
      DO JB=1,NBR
        IF (IWDO(J) >= US(JB) .AND. IWDO(J) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
      END DO
      IF (.NOT.ACTIVE_SEGMENT) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Withdrawal out segment [IWDO=',IWDO(J),'] is not an active segment'
      END IF
    END DO
    DO J=1,NWDO
      IF (WDOF(J) <= 0.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTF) 'Withdrawal out frequency [WDOF=',WDOF(J),'] <= 0'
      END IF
    END DO
  END IF

! Restart

  WRITE (WIN,*) '  restart'
  IF ((RSOC /= '      ON') .AND. (RSOC /= '     OFF')) THEN
    CALL ERRORS
    WRITE (ERR,FMTA) 'Restart output control [RSOC='//RSOC(6:8)//'] must be either " ON" or "OFF"'
  END IF
  IF (RSOC == '      ON') THEN
    IF (NRSO <= 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Number of restart output dates [NRSO=',NRSO,'] <= 0'
    END IF
    DO J=1,NRSO-1
      IF (RSOD(J+1) <= RSOD(J)) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Restart output date [RSOD=',RSOD(J+1),'] <= previous date [RSOD=',RSOD(J),']'
      END IF
    END DO
    DO J=1,NRSO
      IF (RSOF(J) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Restart output frequency [RSOF=',RSOF(J),'] <= 0'
      END IF
    END DO
  END IF

! Constituent controls

  WRITE (WIN,*) 'Constituents'
  WRITE (WIN,*) '  active constituents'
  IF ((CCC /= '      ON') .AND. (CCC /= '     OFF')) THEN
    CALL ERRORS
    WRITE (ERR,FMTA) 'Constituent computation control [CCC='//CCC(6:8)//'] must be either " ON" or "OFF"'
  END IF
  IF (CONSTITUENTS) THEN
    DO JW=1,NWB
      IF ((SEDC(JW) /= '      ON') .AND. (SEDC(JW) /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Sediment computation control [SEDC='//SEDC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
      END IF
      IF ((PRNSC(JW) /= '      ON') .AND. (PRNSC(JW) /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Sediment print control [PRNSC='//PRNSC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
      END IF
    END DO
    IF (CUF > 12) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Kinetic update frequency [CUF=',CUF,'] > 12'
    END IF

!** Active constituents

    DO JAC=1,NAC
      JC = CN(JAC)
      IF ((CAC(JC) /= '      ON') .AND. (CAC(JC) /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Active constituent computation control [CAC='//CAC(JC)(6:8)//'] must be either " ON" or "OFF" '//        &
                         'for constituent ',JC
      END IF
    END DO

!** pH calculations

    DO JW=1,NWB
      IF (CDWBC(20,JW) == '      ON') THEN
        IF (CAC(NTIC) == '     OFF' .OR. CAC(NALK) == '     OFF') THEN
          CALL ERRORS
          WRITE (ERR,FMTA) 'pH computations are =" ON", so TIC and alkalanity must also be turned " ON"'
        END IF
      END IF
    END DO

!** Initial constituent concentrations

    WRITE (WIN,*) '  initial concentrations'
    DO JW=1,NWB
      DO JAC=1,NAC
        JC = CN(JAC)
        IF (C2IWB(JC,JW) < -2.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Initial constituent concentration < -2 for constituent ',JC,' in waterbody ',JW
        END IF
        IF (CAC(JC) == '      ON') THEN
          IF (C2IWB(JC,JW) == 0.0) THEN
            NAME = TRIM(CNAME1(JC))
            IF (NAME /= 'Residence time' .AND. NAME /= 'Water age') THEN
              CALL WARNINGS
              WRITE (WRN,FMTI) 'Initial constituent concentration = 0 for constituent ',JC,' in waterbody ',JW
            END IF
          END IF
        END IF
      END DO
    END DO

!** Constituent output

    WRITE (WIN,*) '  constituent output'
    DO JW=1,NWB
      DO JAC=1,NAC
        JC = CN(JAC)
        IF (CPRWBC(JC,JW) /= '      ON' .AND. CPRWBC(JC,JW) /= '     OFF') THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Constituent print control [CPRWBC='//CPRWBC(JC,JW)(6:8)//'] must be either " ON" or "OFF" for '//      &
                           'constituent ',JC
        END IF
        IF (CPRWBC(JC,JW) == '      ON') THEN
          IF (CAC(JC) /= '      ON') THEN
            CALL WARNINGS
            WRITE (WRN,FMTI) 'Constituent print control [CPRWBC='//CPRWBC(JC,JW)(6:8)//'] is " ON" but active constituent'//       &
                             ' control [CAC='//CAC(JC)(6:8)//'] is not " ON" for constituent ',JC
          END IF
        END IF
      END DO
    END DO

!** Inflow constituent concentrations

    WRITE (WIN,*) '  branch inflow concentrations'
    DO JB=1,NBR
      DO JAC=1,NAC
        JC = CN(JAC)
        IF ((CINBRC(JC,JB) /= '      ON') .AND. (CINBRC(JC,JB) /= '     OFF')) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Inflow active constituent control [CINBRC='//CINBRC(JC,JB)(6:8)//'] must be either " ON" or "OFF" '//  &
                           'for constituent ',JC,' for branch ',JB
        END IF
        IF (CAC(JC) == '      ON') THEN
          IF (CINBRC(JC,JB) /= '      ON') THEN
            NAME = TRIM(CNAME1(JC))
            IF (NAME /= 'Residence time' .AND. NAME /= 'Water age') THEN
              CALL WARNINGS
              WRITE (WRN,FMTI) 'Constituent ',JC,' is " ON", but inflow constituent concentration is "OFF" for branch ',JB
            END IF
          END IF
        END IF
      END DO
    END DO

!** Tributary constituent concentrations

    WRITE (WIN,*) '  tributary inflow concentrations'
    IF (TRIBUTARIES) THEN
      DO JT=1,NTR
        DO JAC=1,NAC
          JC = CN(JAC)
          IF ((CINTRC(JC,JT) /= '      ON') .AND. (CINTRC(JC,JT) /= '     OFF')) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Inflow active constituent control [CINTRC='//CINTRC(JC,JT)(6:8)//'] must be either " ON" or "OFF" '//&
                             'for constituent ',JC,' for tributary ',JT
          END IF
          IF (CAC(JC) == '      ON') THEN
            IF (CINTRC(JC,JT) /= '      ON') THEN
              NAME = TRIM(CNAME1(JC))
              IF (NAME /= 'Residence time' .AND. NAME /= 'Water age') THEN
                CALL WARNINGS
                WRITE (WRN,FMTI) 'Constituent ',JC,' is " ON", but inflow constituent concentration is "OFF" for tributary ',JT
              END IF
            END IF
          END IF
        END DO
      END DO
    END IF

!** Distributed tributary concentrations

    WRITE (WIN,*) '  distributed tributary inflow concentrations'
    DO JB=1,NBR
      IF (DIST_TRIB(JB)) THEN
        DO JAC=1,NAC
          JC = CN(JAC)
          IF ((CDTBRC(JC,JB) /= '      ON') .AND. (CDTBRC(JC,JB) /= '     OFF')) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Distributed tributary inflow active constituent control [CDTBRC='//CDTBRC(JC,JB)(6:8)//'] must be '//&
                             'either " ON" or "OFF" for constituent ',JC,' for branch ',JB
          END IF
          IF (CAC(JC) == '      ON') THEN
            IF (CDTBRC(JC,JB) /= '      ON') THEN
              CALL WARNINGS
              WRITE (WRN,FMTI) 'Constituent ',JC,' is " ON", but distrubuted tributary inflow constituent concentration '//        &
                               'is "OFF" for branch ',JB
            END IF
          END IF
        END DO
      END IF
    END DO

!** Precipitation concentrations

    WRITE (WIN,*) '  precipitation concentrations'
    DO JW=1,NWB
      IF (PRECIPITATION(JW)) THEN
        DO JB=BS(JW),BE(JW)
          DO JAC=1,NAC
            JC = CN(JAC)
            IF ((CPRBRC(JC,JB) /= '      ON') .AND. (CPRBRC(JC,JB) /= '     OFF')) THEN
              CALL ERRORS
              WRITE (ERR,FMTI) 'Precipitation active constituent control [CINBRC='//CINBRC(JC,JB)(6:8)//'] must be either " ON" '//&
                               'or "OFF" for constituent ',JC,' for branch ',JB
            END IF
            IF (CAC(JC) == '      ON') THEN
              IF (CPRBRC(JC,JB) /= '      ON') THEN
                CALL WARNINGS
                WRITE (WRN,FMTI) 'Constituent ',JC,' is " ON", but precipitation constituent concentration is "OFF" for  branch ',JB
              END IF
            END IF
          END DO
        END DO
      END IF
    END DO
  END IF

! Light extinction

  WRITE (WIN,*) '  kinetic coefficients'
  WRITE (WIN,*) '    light extinction'
  DO JW=1,NWB
    IF (EXH2O(JW) > 0.95) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Pure water light extinction > 0.95 for waterbody ',JW
    ELSE IF (EXH2O(JW) < 0.25) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Pure water light extinction < 0.25 for waterbody ',JW
    END IF
    IF (BETA(JW) > 0.75) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Short wave solar radiation absorption [BETA=',BETA(JW),'] > 0.75 for waterbody ',JW
    ELSE IF (BETA(JW) <= 0.3 .AND. BETA(JW) >= 0.1) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Short wave solar radiation absorption [BETA=',BETA(JW),'] < 0.3 for waterbody ', JW
    ELSE IF (BETA(JW) < 0.1) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Short wave solar radiation absorption [BETA=',BETA(JW),'] < 0.1 for waterbody ', JW
    END IF
    IF (CONSTITUENTS) THEN
      IF (SUSPENDED_SOLIDS) THEN
        IF (EXSS(JW) > 0.1) THEN
          CALL WARNINGS
          WRITE (ERR,FMTFI) 'Suspended solids light extinction [EXSS=',EXSS(JW),'] > 0.1 for waterbody ', JW
        ELSE IF (EXSS(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Suspended solids light extinction [EXSS=',EXSS(JW),'] < 0 for waterbody ',   JW
        END IF
      END IF
      IF (LABILE_DOM .OR. REFRACTORY_DOM .OR. LABILE_POM .OR. REFRACTORY_POM) THEN
        IF (EXOM(JW) > 0.85) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Organic matter light extinction [EXOM=',EXOM(JW),'] > 0.85 for waterbody ',  JW
        ELSE IF (EXOM(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter light extinction [EXOM=',EXOM(JW),'] < 0 for waterbody ',     JW
        END IF
      END IF
    END IF
  END DO
  IF (ALGAE) THEN
    DO JA=1,NAAL
      IF (EXA(ALAC(JA)) > 0.8) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal light extinction [EXA=',EXA(ALAC(JA)),'] > 0.8 for algal group ',JW
      ELSE IF (EXA(ALAC(JA)) < 0.1) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal light extinction [EXA=',EXA(ALAC(JA)),'] < 0.1 for algal group ',JW
      END IF
    END DO
  END IF

! Constituents

  IF (CONSTITUENTS) THEN

!** Generic constituent

    WRITE (WIN,*) '    generic'
    IF (GENERIC_CONST) THEN
      DO JG=1,NGC
        IF (CGQ10(JG) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Temperature multiplier [CGQ10]=',CGQ10(JG),'] < 0.0 for generic constituent ', JG
        ELSE IF (CGQ10(JG) > 1.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Temperature multiplier [CGQ10]=',CGQ10(JG),'] > 1.5 for generic constituent ', JG
        END IF
        IF (CG0DK(JG) < 0.0 .AND. CG0DK(JG) /= -1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) '0-order decay rate [CG0DK]=',    CG0DK(JG),'] < 0.0 for generic constituent ', JG
        ELSE IF (CG0DK(JG) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) '0-order decay rate [CG0DK]=',    CG0DK(JG),'] > 10.0 for generic constituent ',JG
        END IF
        IF (CGS(JG) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Settling rate [CGS]=',           CGS(JG),'] < 0.0 for generic constituent ',   JG
        ELSE IF (CGS(JG) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Settling rate [CGS]=',           CGS(JG),'] > 10.0 for generic constituent ',  JG
        END IF
      END DO
    END IF

!** Suspended solids

    WRITE (WIN,*) '    suspended solids'
    IF (SUSPENDED_SOLIDS) THEN
      DO JS=1,NASS
        IF (SSS(SSAC(JS)) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Suspended solids settling rate [SSS=',SSS(JS),'] > 10.0 for suspended solids group ',JS
        ELSE IF (SSS(SSAC(JS)) == 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Suspended solids settling rate [SSS=',SSS(JS),']  = 0.0 for suspended solids group ',JS
        ELSE IF (SSS(SSAC(JS)) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Suspended solids settling rate [SSS=',SSS(JS),']  < 0.0 for suspended solids group ',JS
        END IF
        IF (SEDRC(JS) /= '      ON' .AND. SEDRC(JS) /= '     OFF') THEN                        ! SW 1/16/04
          CALL ERRORS
          WRITE (ERR,FMTI) 'Wind sediment resuspension control   [SEDRC='//SEDRC(JS)//'] /= " ON" or "OFF" for SS group ',JS
        END IF

	  END DO
    END IF

!** Algae
    WRITE (WIN,*) '    algae'
    DO J=1,NAAL
      JA = ALAC(J)
      IF (AG(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal growth rate [AG=',        AG(JA),'] < 0.0 for algal group ',    JA
      ELSE IF (AG(JA) < 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal growth rate [AG=',        AG(JA),'] < 0.1 for algal group ',    JA
      ELSE IF (AG(JA) > 5.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal growth rate [AG=',        AG(JA),'] > 5.0 for algal group ',    JA
      END IF
      IF (AM(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal mortality rate [AM=',     AM(JA),'] < 0.0 for algal group ',    JA
      ELSE IF (AM(JA) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal mortality rate [AM=',     AM(JA),'] < 0.001 for algal group ',  JA
      ELSE IF (AM(JA) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal mortality rate [AM=',     AM(JA),'] > 1.0 for algal group ',    JA
      END IF
      IF (AE(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 0.0 for algal group ',    JA
      ELSE IF (AE(JA) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 0.001 for algal group ',  JA
      ELSE IF (AE(JA) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 1.0 for algal group ',    JA
      END IF
      IF (AR(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 0.0 for algal group ',    JA
      ELSE IF (AR(JA) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 0.001 for algal group ',  JA
      ELSE IF (AR(JA) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] > 1.0 for algal group ',    JA
      END IF
      IF (AS(JA) > 5.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal settling rate [AS=',      AS(JA),'] > 5.0 for algal group ',    JA
      ELSE IF (AS(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal settling rate [AS=',      AS(JA),'] < 0.0 for algal group ',    JA
      END IF
      IF (ASAT(JA) < 20.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal light saturation [ASAT=', ASAT(JA),'] < 20.0 for algal group ', JA
      ELSE IF (ASAT(JA) > 300.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal light saturation [ASAT=', ASAT(JA),'] > 300.0 for algal group ',JA
      END IF

!**** Algal temperature rate constants

      IF (AT1(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature [AT1=',AT1(JA),'] <  0.0 for algal group ',JA
      END IF
      IF (AT2(JA) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal maximum temperature [AT2=',AT2(JA),'] > 40.0 for algal group ',JA
      ELSE IF (AT2(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature [AT2=',AT2(JA),'] <  0.0 for algal group ',JA
      END IF
      IF (AT3(JA) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal maximum temperature [AT3=',AT3(JA),'] > 40.0 for algal group ',JA
      ELSE IF (AT3(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature [AT3=',AT3(JA),'] <  0.0 for algal group ',JA
      END IF
      IF (AT4(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature [AT4=',AT4(JA),'] <  0.0 for algal group ',JA
      END IF
      IF (AK1(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature rate multiplier [AK1=',AK1(JA),'] < 0.0 for algal group ',JA
      END IF
      IF (AK2(JA) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature rate multiplier [AK2=',AK2(JA),'] > 1.0 for algal group ',JA
      ELSE IF (AK2(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature rate multiplier [AK2=',AK2(JA),'] < 0.0 for algal group ',JA
      END IF
      IF (AK3(JA) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature rate multiplier [AK3=',AK3(JA),'] > 1.0 for algal group ',JA
      ELSE IF (AK3(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature rate multiplier [AK3=',AK3(JA),'] < 0.0 for algal group ',JA
      END IF
      IF (AK4(JA) > 1) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature rate multiplier [AK4=',AK4(JA),'] > 1.0 for algal group ',JA
      ELSE IF (AK4(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature rate multiplier [AK4=',AK4(JA),'] < 0.0 for algal group ',JA
      END IF
      IF (AT1(JA) >= AT2(JA)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Algal temperature rate multiplier [AT1=',AT1(JA),'] > [AT2=',AT2(JA),'] for algal group ',JA
      END IF
      IF (AT2(JA) >= AT3(JA)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Algal temperature rate multiplier [AT2=',AT2(JA),'] > [AT3=',AT3(JA),'] for algal group ',JA
      END IF
      IF (AT3(JA) >= AT4(JA)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Algal temperature rate multiplier [AT3=',AT3(JA),'] > [AT4=',AT4(JA),'] for algal group ',JA
      END IF

!**** Algal stoichiometry

      IF (ALGP(JA) > 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Phosphorus to algal stoichiometry [ALGP=',ALGP(JA),'] > 0.1 for algal group ',   JA
      ELSE IF (ALGP(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Phosphorus to algal stoichiometry [ALGP=',ALGP(JA),'] < 0.0 for algal group ',   JA
      END IF
      IF (ALGN(JA) > 0.3) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Nitrogen to algal stoichiometry [ALGN=',  ALGN(JA),'] > 0.3 for algal group ',   JA
      ELSE IF (ALGN(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Nitrogen to algal stoichiometry [ALGN=',  ALGN(JA),'] < 0.0 for algal group ',   JA
      END IF
      IF (ALGC(JA) > 0.6) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Carbon to algal stoichiometry [ALGC=',    ALGC(JA),'] > 0.6 for algal group ',   JA
      ELSE IF (ALGC(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Carbon to algal stoichiometry [ALGC=',    ALGC(JA),'] < 0.0 for algal group ',   JA
      END IF
      IF (ALGSI(JA) > 0.5) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Silica to algal stoichiometry [ALGSI=',   ALGSI(JA),'] > 0.5 for algal group ',  JA
      ELSE IF (ALGSI(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Silica to algal stoichiometry [ALGSI=',   ALGSI(JA),'] < 0.0 for algal group ',  JA
      END IF
      IF (ACHLA(JA) > 0.8) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal to chlorophyll a ratio mg/ug [ACHLA=',    ACHLA(JA),'] > 0.8 for algal group ',JA
      ELSE IF (ACHLA(JA) < 0.01) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal to chlorophyll a ratio mg/ug [ACHLA=',    ACHLA(JA),'] < 0.01 for algal group ', JA
      ELSE IF (ACHLA(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal to chlorophyll a ratio mg/ug [ACHLA=',    ACHLA(JA),'] < 0.0 for algal group ',  JA
      END IF
      IF (APOM(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal mortality to POM fraction [APOM=',  APOM(JA),'] < 0.0 for algal group ',   JA
      ELSE IF (APOM(JA) < 0.5) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal mortality to POM fraction [APOM='  ,APOM(JA),'] < 0.5 for algal group ',   JA
      END IF
    END DO

!** Epiphyton

    WRITE (WIN,*) '    epiphyton'
    DO JE=1,NEP
      IF (EG(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton growth rate [EG=',        EG(JE),'] < 0.0 for epiphyton group ',    JE
      ELSE IF (EG(JE) < 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton growth rate [EG=',        EG(JE),'] < 0.1 for epiphyton group ',    JE
      ELSE IF (EG(JE) > 5.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton growth rate [EG=',        EG(JE),'] > 5.0 for epiphyton group ',    JE
      END IF
      IF (EM(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton mortality rate [EM=',     EM(JE),'] < 0.0 for epiphyton group ',    JE
      ELSE IF (EM(JE) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton mortality rate [EM=',     EM(JE),'] < 0.001 for epiphyton group ',  JE
      ELSE IF (EM(JE) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton mortality rate [EM=',     EM(JE),'] > 1.0 for epiphyton group ',    JE
      END IF
      IF (EE(JE) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton excretion rate [EE=',     EE(JE),'] < 0.0 for epiphyton group ',    JE
      ELSE IF (EE(JE) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton excretion rate [EE=',     EE(JE),'] < 0.001 for epiphyton group ',  JE
      ELSE IF (EE(JE) > 1.0) THEN
         CALL WARNINGS
         WRITE (WRN,FMTFI) 'Epiphyton excretion rate [EE=',    EE(JE),'] < 1.0 for epiphyton group ',    JE
      END IF
      IF (ER(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton respiration rate [ER=',   ER(JE),'] < 0.0 for epiphyton group ',    JE
      ELSE IF (ER(JE) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton respiration rate [ER=',   ER(JE),'] < 0.001 for epiphyton group ',  JE
      ELSE IF (ER(JE) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton respiration rate [ER=',   ER(JE),'] > 1.0 for epiphyton group ',    JE
      END IF
      IF (EB(JE) > 0.01) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton burial rate [ES=',        EB(JE),'] > 0.01 for epiphyton group ',   JE
      ELSE IF (EB(JE) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton burial rate [ES=',        EB(JE),'] < 0.0 for epiphyton group ',    JE
      END IF
      IF (ESAT(JE) < 20.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton light saturation [ESAT=', ESAT(JE),'] < 20.0 for epiphyton group ', JE
      ELSE IF (ESAT(JE) > 300.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton light saturation [ESAT=', ESAT(JE),'] > 300.0 for epiphyton group ',JE
      END IF

!**** Epiphyton temperature rate constants

      IF (ET1(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature [ET1=',ET1(JE),'] <  0.0 for epiphyton group ',JE
      END IF
      IF (ET2(JE) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton maximum temperature [ET2=',ET2(JE),'] > 40.0 for epiphyton group ',JE
      ELSE IF (ET2(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature [ET2=',ET2(JE),'] <  0.0 for epiphyton group ',JE
      END IF
      IF (ET3(JE) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton maximum temperature [ET3=',ET3(JE),'] > 40.0 for epiphyton group ',JE
      ELSE IF (ET3(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature [ET3=',ET3(JE),'] <  0.0 for epiphyton group ',JE
      END IF
      IF (ET4(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature [ET4=',ET4(JE),'] <  0.0 for epiphyton group ',JE
      END IF
      IF (EK1(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature rate multiplier [EK1=',EK1(JE),'] < 0.0 for epiphyton group ',JE
      END IF
      IF (EK2(JE) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature rate multiplier [EK2=',EK2(JE),'] > 1.0 for epiphyton group ',JE
      ELSE IF (EK2(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature rate multiplier [EK2=',EK2(JE),'] < 0.0 for epiphyton group ',JE
      END IF
      IF (EK3(JE) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature rate multiplier [EK3=',EK3(JE),'] > 1.0 for epiphyton group ',JE
      ELSE IF (EK3(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature rate multiplier [EK3=',EK3(JE),'] < 0.0 for epiphyton group ',JE
      END IF
      IF (EK4(JE) > 1) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature rate multiplier [EK4=',EK4(JE),'] > 1.0 for epiphyton group ',JE
      ELSE IF (EK4(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature rate multiplier [EK4=',EK4(JE),'] < 0.0 for epiphyton group ',JE
      END IF
      IF (ET1(JE) >= ET2(JE)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Epiphyton temperature rate multiplier [ET1=',ET1(JE),'] > [ET2=',ET2(JE),'] for epiphyton group ',JE
      END IF
      IF (ET2(JE) >= ET3(JE)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Epiphyton temperature rate multiplier [ET2=',ET2(JE),'] > [ET3=',ET3(JE),'] for epiphyton group ',JE
      END IF
      IF (ET3(JE) >= ET4(JE)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Epiphyton temperature rate multiplier [ET3=',ET3(JE),'] > [ET4=',ET4(JE),'] for epiphyton group ',JE
      END IF

!**** Epiphyton stoichiometry

      IF (EP(JE) > 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Phosphorus to epiphyton stoichiometry [EP=', EP(JE),'] > 0.1 for epiphyton group ',     JE
      ELSE IF (EP(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Phosphorus to epiphyton stoichiometry [EP=', EP(JE),'] < 0.0 for epiphyton group ',     JE
      END IF
      IF (EN(JE) > 0.3) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Nitrogen to epiphyton stoichiometry [EN=',   EN(JE),'] > 0.3 for epiphyton group ',     JE
      ELSE IF (EN(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Nitrogen to epiphyton stoichiometry [EN=',   EN(JE),'] < 0.0 for epiphyton group ',     JE
      END IF
      IF (EC(JE) > 0.6) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Carbon to epiphyton stoichiometry [EC=',     EC(JE),'] > 0.6 for epiphyton group ',     JE
      ELSE IF (EC(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Carbon to epiphyton stoichiometry [EC=',     EC(JE),'] < 0.0 for epiphyton group ',     JE
      END IF
      IF (ESI(JE) > 0.5) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Silica to epiphyton stoichiometry [ESI=',    ESI(JE),'] > 0.5 for epiphyton group ',    JE
      ELSE IF (ESI(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Silica to epiphyton stoichiometry [ESI=',    ESI(JE),'] < 0.0 for epiphyton group ',    JE
      END IF
      IF (ECHLA(JE) > 0.8) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton to chlorophyll a ratio mg/ug [ECHLA=',   ECHLA(JE),'] > 0.8 for epiphyton group ',JE
      ELSE IF (ECHLA(JE) < 0.01) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton to chlorophyll a ratio mg/ug [ECHLA=',   ECHLA(JE),'] < 0.01 for epiphyton group ', JE
      ELSE IF (ECHLA(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton to chlorophyll a ratio mg/ug [ECHLA=',   ECHLA(JE),'] < 0.0 for epiphyton group ',  JE
      END IF
      IF (EPOM(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton mortality to POM fraction [EPOM=', EPOM(JE),'] < 0.0 for epiphyton group ',   JE
      ELSE IF (EPOM(JE) < 0.5) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton mortality to POM fraction [EPOM=' ,EPOM(JE),'] < 0.5 for epiphyton group ',   JE
      END IF

      do jw=1,nwb
        IF (EPIC(JW,JE) /= '     OFF' .and. epic(JW,JE) /= '      ON') THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,A,A,I0,a,i0))') 'Epiphyton ON/OFF control [EPIC=', EPIC(JW,JE),'] must be ON or OFF for epiphyton group ',   JE,' and JW=',JW
        endif
        IF (EPIPRC(JW,JE) /= '     OFF' .and. epiPRc(JW,JE) /= '      ON') THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,A,A,I0,a,i0))') 'Epiphyton ON/OFF print control [EPIPRC=', EPIPRC(JW,JE),'] must be ON or OFF for epiphyton group ',   JE,' and JW=',JW
        endif
      enddo

    END DO

!** Macrophytes ******************

    WRITE (WIN,*) '    macropytes'
    DO JM=1,NMC

      do jw=1,nwb
        IF (MACWBC(JW,JM) /= '     OFF' .and. macwbc(JW,Jm) /= '      ON') THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,A,A,I0,a,i0))') 'Macrophyte ON/OFF control [MACWBC=', MACWBC(JW,JM),'] must be ON or OFF for macrophyte group ',   JM,' and JW=',JW
        endif
        IF (MPRWBC(JW,JM) /= '     OFF' .and. mprwbc(JW,Jm) /= '      ON') THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,A,A,I0,a,i0))') 'Macrophyte ON/OFF print control [MPRWBC=', MPRWBC(JW,JM),'] must be ON or OFF for macrophyte group ',   JM,' and JW=',JW
        endif
      enddo

      IF (MG(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte growth rate [MG=',        MG(JM),'] < 0.0 for macrophyte group ',    JM
      ELSE IF (MG(JM) < 0.01) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte growth rate [MG=',        MG(JM),'] < 0.01 for macrophyte group ',    JM
      ELSE IF (MG(JM) > 2.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte growth rate [MG=',        MG(JM),'] > 2.0 for macrophyte group ',    JM
      END IF
      IF (MM(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte mortality rate [MM=',     MM(JM),'] < 0.0 for macrophyte group ',    JM
      ELSE IF (MM(JM) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte mortality rate [MM=',     MM(JM),'] < 0.001 for macrophyte group ',  JM
      ELSE IF (MM(JM) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte mortality rate [MM=',     MM(JM),'] > 1.0 for macrophyte group ',    JM
      END IF
      IF (MR(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte respiration rate [MR=',   MR(JM),'] < 0.0 for macrophyte group ',    JM
      ELSE IF (MR(JM) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte respiration rate [MR=',   MR(JM),'] < 0.001 for macrophyte group ',  JM
      ELSE IF (MR(JM) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte respiration rate [MR=',   MR(JM),'] > 1.0 for macrophyte group ',    JM
      END IF

      IF (MSAT(JM) < 15.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte light saturation [MSAT=', MSAT(JM),'] < 15.0 for macrophyte group ', JM
      ELSE IF (MSAT(JM) > 300.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte light saturation [MSAT=', MSAT(JM),'] > 300.0 for macrophyte group ',JM
      END IF


      IF (PSED(JM) < 0.0 .or. PSED(JM) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte sediment P uptake fraction [PSED=', PSED(JM),'] < 0 or > 1 for macrophyte group ', JM
      END IF
      IF (NSED(JM) < 0.0 .or. NSED(JM) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte sediment N uptake fraction [NSED=', NSED(JM),'] < 0 or > 1 for macrophyte group ', JM
      END IF
      IF (MBMP(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte threshold conc to move upward [MBMP=', MBMP(JM),'] < 0 for macrophyte group ', JM
       ELSE IF (MBMP(JM) > 150.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte threshold conc to move upward [MBMP=', MBMP(JM),'] > 150 for macrophyte group ', JM
      END IF
      IF (MMAX(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum concentration [MMAX=', MMAX(JM),'] < or = 0 for macrophyte group ', JM
       ELSE IF (MMAX(JM) > 1500.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte maximum concentration [MMAX=', MMAX(JM),'] > 1500 for macrophyte group ', JM
      END IF

      IF (dmv(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte dry weight to wet volume ratio [DMV=', DMV(JM),'] < or = 0 for macrophyte group ', JM
       ELSE IF (dmv(JM) > 1.0e7) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte dry weight to wet volume ratio [DMV=', DMV(JM),'] > 1E7 for macrophyte group ', JM
      END IF
      IF (cdstem(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte CD coefficient [DMV=', DMV(JM),'] < or = 0 for macrophyte group ', JM
       ELSE IF (cdstem(JM) < 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte CD coefficient [DMV=', DMV(JM),'] < 1 for macrophyte group ', JM
       ELSE IF (cdstem(JM) > 3.5) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte CD coefficient [DMV=', DMV(JM),'] > 3.5 for macrophyte group ', JM
      END IF

        IF (dwsa(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte dry weight to surface area [DWSA=', DWSA(JM),'] < or = 0 for macrophyte group ', JM
       ELSE IF (dwsa(JM) > 100.) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte dry weight to surface area [DWSA=', DWSA(JM),'] > 100. for macrophyte group ', JM
      END IF

      IF (anorm(JM) <= 0.0 .or. anorm(jm) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte area fraction normal to flow [ANORM=', ANORM(JM),'] < or = 0 or > 1.0 for macrophyte group ', JM
      END IF

!**** Macrophyte temperature rate constants

      IF (MT1(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature [MT1=',MT1(JM),'] <  0.0 for macrophyte group ',JM
      END IF
      IF (MT2(JM) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte maximum temperature [MT2=',MT2(JM),'] > 40.0 for macrophyte group ',JM
      ELSE IF (MT2(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature [MT2=',MT2(JM),'] <  0.0 for macrophyte group ',JM
      END IF
      IF (MT3(JM) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte maximum temperature [MT3=',MT3(JM),'] > 40.0 for macrophyte group ',JM
      ELSE IF (MT3(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature [MT3=',MT3(JM),'] <  0.0 for macrophyte group ',JM
      END IF
      IF (MT4(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature [MT4=',MT4(JM),'] <  0.0 for macrophyte group ',JM
      END IF
      IF (MK1(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature rate multiplier [MK1=',MK1(JM),'] < 0.0 for macrophyte group ',JM
      END IF
      IF (MK2(JM) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature rate multiplier [MK2=',MK2(JM),'] > 1.0 for macrophyte group ',JM
      ELSE IF (MK2(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature rate multiplier [MK2=',MK2(JM),'] < 0.0 for macrophyte group ',JM
      END IF
      IF (MK3(JM) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature rate multiplier [MK3=',MK3(JM),'] > 1.0 for macrophyte group ',JM
      ELSE IF (MK3(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature rate multiplier [MK3=',MK3(JM),'] < 0.0 for macrophyte group ',JM
      END IF
      IF (MK4(JM) > 1) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature rate multiplier [MK4=',MK4(JM),'] > 1.0 for macrophyte group ',JM
      ELSE IF (MK4(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature rate multiplier [MK4=',MK4(JM),'] < 0.0 for macrophyte group ',JM
      END IF
      IF (MT1(JM) >= MT2(JM)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Macrophyte temperature rate multiplier [MT1=',MT1(JM),'] > [MT2=',MT2(JM),'] for macrophyte group ',JM
      END IF
      IF (MT2(JM) >= MT3(JM)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Macrophyte temperature rate multiplier [MT2=',MT2(JM),'] > [MT3=',MT3(JM),'] for macrophyte group ',JM
      END IF
      IF (MT3(JM) >= MT4(JM)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Macrophyte temperature rate multiplier [MT3=',MT3(JM),'] > [MT4=',MT4(JM),'] for macrophyte group ',JM
      END IF

!**** Macrophyte stoichiometry

      IF (MP(JM) > 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Phosphorus to macrophyte stoichiometry [MP=', MP(JM),'] > 0.1 for macrophyte group ',     JM
      ELSE IF (MP(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Phosphorus to macrophyte stoichiometry [MP=', MP(JM),'] < 0.0 for macrophyte group ',     JM
      END IF
      IF (MN(JM) > 0.3) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Nitrogen to macrophyte stoichiometry [MN=',   MN(JM),'] > 0.3 for macrophyte group ',     JM
      ELSE IF (MN(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Nitrogen to macrophyte stoichiometry [MN=',   MN(JM),'] < 0.0 for macrophyte group ',     JM
      END IF
      IF (MC(JM) > 0.6) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Carbon to macrophyte stoichiometry [MC=',     MC(JM),'] > 0.6 for macrophyte group ',     JM
      ELSE IF (MC(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Carbon to macrophyte stoichiometry [MC=',     MC(JM),'] < 0.0 for macrophyte group ',     JM
      END IF

      IF (MPOM(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte mortality to POM fraction [MPOM=', MPOM(JM),'] < 0.0 for macrophyte group ',   JM
      ELSE IF (MPOM(JM) < 0.5) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte mortality to POM fraction [MPOM=' ,MPOM(JM),'] < 0.5 for macrophyte group ',   JM
      END IF
    END DO

!**Zooplankton - to be written**************

!** Dissolved organic matter

    WRITE (WIN,*) '    dissolved organic matter'
    IF (LABILE_DOM) THEN
      DO JW=1,NWB
        IF (LDOMDK(JW) > 0.3) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Labile dissolved organic matter decay rate [LDOMDK=',LDOMDK(JW),'] > 0.3 for waterbody ',JW
        ELSE IF (LDOMDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Labile dissolved organic matter decay rate [LDOMDK=',LDOMDK(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (LRDDK(JW) > 0.01) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Labile to refractory DOM decay rate [LRDDK=',LRDDK(JW),'] > 0.01 for waterbody ',JW
        ELSE IF (LRDDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Labile to refractory DOM decay rate [LRDDK=',LRDDK(JW),'] <  0.0 for waterbody ',JW
        END IF
      END DO
    END IF
    IF (REFRACTORY_DOM) THEN
      DO JW=1,NWB
        IF (RDOMDK(JW) > 0.01) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Refractory dissolved organic matter decay rate [RDOMDK=',RDOMDK(JW),'] > 0.01 for waterbody ',JW
        ELSE IF (RDOMDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Refractory dissolved organic matter decay rate [RDOMDK=',RDOMDK(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
    END IF

!** Particulate organic matter

    WRITE (WIN,*) '    particulate organic matter'
    IF (LABILE_POM) THEN
      DO JW=1,NWB
        IF (LPOMDK(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Labile particulate organic matter decay rate [LPOMDK=',LPOMDK(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (LPOMDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Labile particulate organic matter decay rate [LPOMDK=',LPOMDK(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (POMS(JW) > 1.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Particulate organic matter settling rate [POMS=',POMS(JW),'] > 1.0 for waterbody ',JW
        ELSE IF (POMS(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Particulate organic matter settling rate [POMS=',POMS(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
    END IF
    IF (REFRACTORY_POM) THEN
      DO JW=1,NWB
        IF (RPOMDK(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Refractory particulate organic matter decay rate [RPOMDK=',RPOMDK(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (RPOMDK(JW) < 0.0) THEN
          WRITE (WRN,FMTFI) 'Refractory particulate organic matter decay rate [RPOMDK=',RPOMDK(JW),'] < 0.0 for waterbody ',JW
          CALL ERRORS
        END IF
      END DO
    END IF

!** Organic matter temperature rate multipliers

    IF (LABILE_DOM .OR. REFRACTORY_DOM .OR. LABILE_POM .OR. REFRACTORY_POM) THEN
      DO JW=1,NWB
        IF (OMT1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay minimum temperature [OMT1=',OMT1(JW),'] < 0.0 for waterbody ', JW
        ELSE IF (OMT1(JW) > 15) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Organic matter decay minimum temperature [OMT1=',OMT1(JW),'] > 15.0 for waterbody ',JW
        END IF
        IF (OMT2(JW) > 40.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Organic matter decay maximum temperature [OMT2=',OMT2(JW),'] > 40.0 for waterbody ',JW
        ELSE IF (OMT2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay maximum temperature [OMT2=',OMT2(JW),'] < 0.0 for waterbody ', JW
        END IF
        IF (OMK1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay minimum temperature rate multiplier [OMK1=',OMK1(JW),'] < 0.0 for waterbody ',JW
        ELSE IF (OMK1(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Organic matter decay minimum temperature rate multiplier [OMK1=',OMK1(JW),'] > 0.5 for waterbody ',JW
        END IF
        IF (OMK2(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay maximum temperature rate multiplier [OMK2=',OMK2(JW),'] > 1.0 for waterbody ',JW
        ELSE IF (OMK2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay maximum temperature rate multiplier [OMK2=',OMK2(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
    END IF

!** Sediments

    WRITE (WIN,*) '    sediments'
    DO JW=1,NWB
      IF (SEDIMENT_CALC(JW)) THEN
        IF (SEDDK(JW) > 0.3) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment decay rate [SEDDK=',SEDDK(JW),'] > 0.3 for waterbody ',JW
        ELSE IF (SEDDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay rate [SEDDK=',SEDDK(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (SEDs(JW) > 1.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment focusing velocity or settling rate [SEDS=',SEDS(JW),'] > 1.0 m/d for waterbody ',JW
        ELSE IF (SEDs(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment focusing velocity or settling rate [SEDS=',SEDS(JW),'] < 0.0 for waterbody ',JW
        END IF

        IF (SEDbr(JW) > 1.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment burial rate [SEDBR=',SEDBR(JW),'] > 1.0 d^-1 for waterbody ',JW
        ELSE IF (SEDbr(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment burial rate [SEDBR=',SEDBR(JW),'] < 0.0 for waterbody ',JW
        END IF

        IF (DYNSEDK(JW) /= '     OFF' .AND. DYNSEDK(JW) /= '      ON')THEN
          CALL ERRORS
          WRITE (ERR,"('Dynamic calculation of sediment decay control [DYNSEDK=',A8,'] must be ON or OFF for waterbody ',I4)")DYNSEDK(JW),JW
        END IF


        IF (SODT1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay minimum temperature [SODT1=',SODT1(JW),'] < 0.0 for waterbody ', JW
        ELSE IF (SODT1(JW) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment decay minimum temperature [SODT1=',SODT1(JW),'] > 10.0 for waterbody ',JW
        END IF
        IF (SODT2(JW) > 40.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment decay maximum temperature [SODT2=',SODT2(JW),'] > 40.0 for waterbody ',JW
        ELSE IF (SODT2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay maximum temperature [SODT2=',SODT2(JW),'] < 0.0 for waterbody ', JW
        END IF
        IF (SODK1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay minimum temperature rate multiplier [SODK1=',SODK1(JW),'] < 0.0 for waterbody ',JW
        ELSE IF (SODK1(JW) > 0.5) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay minimum temperature rate multiplier [SODK1=',SODK1(JW),'] > 0.5 for waterbody ',JW
        END IF
        IF (SODK2(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay maximum temperature rate multiplier [SODK1=',SODK1(JW),'] > 1.0 for waterbody ',JW
        ELSE IF (SODK2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay maximum temperature rate multiplier [SODK1=',SODK1(JW),'] < 0.0 for waterbody ',JW
        END IF
      END IF
    END DO

!** Carbonaceous biochemical oxygen demand

    WRITE (WIN,*) '    CBOD'
    IF (CBO_DEMAND) THEN
      DO JB=1,NBOD
        IF (KBOD(JB) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'CBOD decay rate [KBOD=',             KBOD(JB),'] < 0.0 for CBOD group ',JB
        END IF
        IF (TBOD(JB) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'CBOD temperature coefficient [TBOD=',TBOD(JB),'] < 0.0 for CBOD group ',JB
        END IF
        IF (RBOD(JB) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'CBOD5/ultimate CBOD ratio [RBOD=',   RBOD(JB),'] < 0.0 for CBOD group ',JB
        END IF
        IF (CBODS(JB) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'CBOD settling rate [CBODS=',   CBODS(JB),'] < 0.0 for CBOD group ',JB
        END IF
      END DO
    END IF

!** Phosphorus

    WRITE (WIN,*) '    phosphorus'
    IF (PHOSPHORUS) THEN
      DO JW=1,NWB
        IF (PO4R(JW) > 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Phosphorus sediment release rate [PO4R=',PO4R(JW),'] > 0.05 for waterbody ',JW
        ELSE IF (PO4R(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Phosphorus sediment release rate [PO4R=',PO4R(JW),'] < 0.0 for waterbody ', JW
        END IF
        IF (PARTP(JW) > 2.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Phosphorus/suspended solids partitioning coefficient [PARTP=',PARTP(JW),'] > 2.0 for waterbody ',JW
        ELSE IF (PARTP(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Phosphorus/suspended solids partitioning coefficient [PARTP=',PARTP(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
      DO J=1,NAAL
        JA = ALAC(J)
        IF (AHSP(JA) > 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Algal phosphorus half-saturation concentration [AHSP=',AHSP(JA),'] > 0.05 for algal group ',JA
        ELSE IF (AHSP(JA) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Algal phosphorus half-saturation concentration [AHSP=',AHSP(JA),'] < 0.00 for algal group ',JA
        END IF
      END DO
    END IF

!** Ammonium

    WRITE (WIN,*) '    ammonium'
    IF (AMMONIUM) THEN
      DO JW=1,NWB
        IF (NH4R(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Ammonium sediment release rate [NH4R=',NH4R(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (NH4R(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Ammonium sediment release rate [NH4R=',NH4R(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (NH4DK(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Nitrification rate [NH4DK=',NH4DK(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (NH4DK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification rate [NH4DK=',NH4DK(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
      DO J=1,NAAL
        JA = ALAC(JA)
        IF (AHSN(JA) > 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Nitrogen half-saturation concentration for algal uptake [AHSN=',AHSN(JA),'] > 0.05 for waterbody ',JW
        ELSE IF (AHSN(JA) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrogen half-saturation concentration for algal uptake [AHSN=',AHSN(JA),'] < 0.00 for waterbody ',JW
        END IF
      END DO

!**** Nitrification temperature rate multipliers

      DO JW=1,NWB
        IF (NH4T1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification minimum temperature [NH4T1=',NH4T1(JW),'] < 0.00 for waterbody ',JW
        ELSE IF (NH4T1(JW) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Nitrification minimum temperature [NH4T1=',NH4T1(JW),'] > 10.0 for waterbody ',JW
        END IF
        IF (NH4T2(JW) > 40.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Nitrification maximum temperature [NH4T2=',NH4T2(JW),'] > 40.0 for waterbody ',JW
        ELSE IF (NH4T2(JW) < 10.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification maximum temperature [NH4T2=',NH4T2(JW),'] < 10.0 for waterbody ',JW
        END IF
        IF (NH4K1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification minimum temperature rate multiplier [NH4K2=',NH4K1(JW),'] < 0.0 for waterbody ',JW
        ELSE IF (NH4K1(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification minimum temperature rate multiplier [NH4K2=',NH4K1(JW),'] > 1.0 for waterbody ',JW
        END IF
        IF (NH4K2(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification maximum temperature rate multiplier [NH4K2=',NH4K2(JW),'] > 1.0 for waterbody ',JW
        ELSE IF (NH4K2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification maximum temperature rate multiplier [NH4K2=',NH4K2(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
    END IF

!** Nitrate

    WRITE (WIN,*) '    nitrate'
    DO JW=1,NWB
      IF (NITRATE) THEN
        IF (NO3DK(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Denitrification rate [NO3DK=',NO3DK(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (NO3DK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification rate [NO3DK=',NO3DK(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (NO3S(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,"('Rate of diffusion of NO3-N into bottom muds [NO3S=',f8.4,'] > 0.5 m/d for waterbody ',i4)")NO3S(JW),JW
        ELSE IF (NO3S(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,"('Rate of diffusion of NO3-N into bottom muds [NO3S=',f8.4,']  < 0.0 for waterbody ',i4)")NO3S(JW),JW
        END IF
        IF (FNO3SED(JW) > 0.0.and. NO3S(JW)>0.0) THEN
          CALL WARNINGS
          WRITE (WRN,"('NO3 is being diffused into sediments and incorporated into sediment organic matter for waterbody ',i4)")JW
          WRITE (WRN,"('  [This may have unintended consequences for OM:SED-N ratio and unexpected high release of NH4 from sediments as SED-N increases]')")
        ENDIF

        IF (FNO3SED(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,"('Fraction[0-1] of NO3 diffused into sediments incorporated into OM [FNO3SED=',f8.4,']  > 1.0 for waterbody ',i4)")FNO3SED(JW),JW
        ELSE IF (FNO3SED(JW)> 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,"('Fraction of NO3 diffused into sediments incorporated into OM [FNO3SED=',f8.4,'] > 0.5 for waterbody ',i4)")FNO3SED(JW),JW
        ELSE IF (FNO3SED(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,"('Fraction[0-1] of NO3 diffused into sediments incorporated into OM [FNO3SED=',f8.4,']  < 0.0 for waterbody ',i4)")FNO3SED(JW),JW
        END IF

        IF (NO3T1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification minimum temperature [NO3T1=',NO3T1(JW),'] < 0.00 for waterbody ',JW
        END IF
        IF (NO3T2(JW) > 40.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Denitrification maximum temperature [NO3T2=',NO3T2(JW),'] > 40.0 for waterbody ',JW
        ELSE IF (NO3T2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification maximum temperature [NO3T2=',NO3T2(JW),'] < 0.0 for waterbody ', JW
        END IF
        IF (NO3K1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification minimum temperature rate multiplier [NO3K1=',NO3K1(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (NO3K2(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification maximum temperature rate multiplier [NO3K1=',NO3K1(JW),'] < 0.0 for waterbody ',JW
        ELSE IF (NO3K2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification maximum temperature rate multiplier [NO3K1=',NO3K1(JW),'] < 0.0 for waterbody ',JW
        END IF
      END IF
    END DO

!** Carbon dioxide sediment release

    WRITE (WIN,*) '    carbon dioxide'
    DO JW=1,NWB
      IF (TOT_IC) THEN
        IF (CO2R(JW) > 1.4) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Carbon dioxide sediment release rate [CO2R=',CO2R(JW),'] > 1.4 for waterbody ',JW
        ELSE IF (CO2R(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Carbon dioxide sediment release rate [CO2R=',CO2R(JW),'] < 0.0 for waterbody ',JW
        END IF
      END IF
    END DO

!** Iron

    WRITE (WIN,*) '    iron'
    DO JW=1,NWB
      IF (IRON) THEN
        IF (FER(JW) > 2.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Iron sediment release rate [FER=',FER(JW),'] > 2.0 for waterbody ',JW
        ELSE IF (FER(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Iron sediment release rate [FER=',FER(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (FES(JW) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Iron settling rate [FES=',       FES(JW),'] > 10.0 for waterbody ',JW
        ELSE IF (FES(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Iron settling rate [FES=',       FES(JW),'] < 0.0 for waterbody ', JW
        END IF
      END IF
    END DO

!** Stoichiometry

    WRITE (WIN,*) '    stoichiometry'
    DO JW=1,NWB
      IF (O2NH4(JW) /= 4.57) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to ammonium stoichiometry [O2NH4=',       O2NH4(JW),'] /= 4.57 for waterbody ', JW
      END IF
      IF (O2OM(JW) /= 1.4) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to organic matter stoichiometry [O2OM=',  O2OM(JW),'] /= 1.4 for waterbody ',   JW
      END IF
      IF (ORGP(JW) /= 0.005) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Phosphorus/organic matter stoichiometry [ORGP=',  ORGP(JW),'] /= 0.005 for waterbody ',JW
      END IF
      IF (ORGN(JW) /= 0.08) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Nitrogen/organic matter stoichiometry [ORGN=',    ORGN(JW),'] /= 0.08 for waterbody ', JW
      END IF
      IF (ORGC(JW) /= 0.45) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Carbon/organic matter stoichiometry [ORGC=',      ORGC(JW),'] /= 0.45 for waterbody ', JW
      END IF
    END DO
    DO JA=1,NAL
      IF (O2AR(JA) /= 1.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to algal respiration stoichiometry [O2AR=',O2AR(JA),'] /= 1.1 for algal group', JA
      END IF
      IF (O2AG(JA) /= 1.4) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to algal production stoichiometry [O2AG=', O2AG(JA),'] /= 1.4 for algal group', JA
      END IF
    END DO
    DO JE=1,NEP
      IF (O2ER(JE) /= 1.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to epiphyton respiration stoichiometry [O2ER=',O2ER(JE),'] /= 1.1 for epiphyton group', JE
      END IF
      IF (O2EG(JE) /= 1.4) THEN                                                                                        !CB 04/08/03
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to epiphyton production stoichiometry [O2EG=', O2EG(JE),'] /= 1.4 for epiphyton group', JE
      END IF
    END DO

!** Oxygen limit

    WRITE (WIN,*) '    DO half saturation coeff KDO'
    IF (DISSOLVED_OXYGEN) THEN
      IF (O2LIM > 2.0) THEN
        CALL WARNINGS
        WRITE (WRN,"('DO half-saturation coeff for anoxia [KDO=',f8.3,'] > 2.0 mg/l')")o2lim
      ELSE IF (O2LIM <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,"('DO half-saturation coeff for anoxia [KDO=',f8.3,'] <= 0.0')")o2lim
      END IF
    END IF
  END IF

! Sediment oxygen demand

  WRITE (WIN,*) '    sediment oxygen demand'
  IF (DISSOLVED_OXYGEN) THEN
    DO I=1,IMX
      IF (SOD(I) > 5.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Sediment oxygen demand [SOD=',SOD(I),'] > 5.0 for waterbody ',JW
      ELSE IF (SOD(I) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Sediment oxygen demand [SOD=',SOD(I),'] < 0.0 for waterbody ',JW
      END IF
    END DO
  END IF

! Reaeration formulae

  WRITE (WIN,*) '    reaeration'
  DO JW=1,NWB
    IF (REAERC(JW) == '   RIVER') THEN
      DO JB=BS(JW),BE(JW)
        IF (SLOPE(JB) == 0) THEN
          IF (NEQN(JW) == 3 .OR. NEQN(JW) == 7 .OR. NEQN(JW) == 8 .OR. NEQN(JW) == 9) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
          END IF
        END IF
      END DO
      IF (NEQN(JW) < 0 .OR. NEQN(JW) > 9) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
      END IF
      IF (NEQN(JW) == 9) THEN
        IF (RCOEF1(JW) == 0.0 .AND. RCOEF2(JW) == 0.0 .AND. RCOEF3(JW) == 0.0 .AND. RCOEF4(JW) == 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
        END IF
      END IF
    ELSE IF (REAERC(JW) == '    LAKE') THEN
      IF (NEQN(JW) < 1 .OR. NEQN(JW) > 14) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
      END IF
      IF (NEQN(JW) == 14) THEN
        IF (RCOEF1(JW) == 0.0 .AND. RCOEF2(JW) == 0.0 .AND. RCOEF3(JW) == 0.0 .AND. RCOEF4(JW) == 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
        END IF
      END IF
    ELSE IF (REAERC(JW) == ' ESTUARY') THEN
      IF (NEQN(JW) < 0 .OR. NEQN(JW) > 2) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
      END IF
    ELSE
      CALL ERRORS
      WRITE (ERR,*) JW, CONFN
      WRITE (ERR,FMTI) 'Reaeration type [REAERC='//REAERC(JW)//'] /= "RIVER", "LAKE" or "ESTUARY" for waterbody ',JW
    END IF
  END DO

! Filenames

  WRITE (WIN,*) 'Filenames'
  DO JW=1,NWB
    BTHFN1 = BTHFN(JW); METFN1 = METFN(JW); SNPFN1 = SNPFN(JW); PRFFN1 = PRFFN(JW); SPRFN1 = SPRFN(JW); VPLFN1 = VPLFN(JW)
    CPLFN1 = CPLFN(JW); FLXFN1 = FLXFN(JW)
    DO JJW=1,NWB
      IF (JJW /= JW .AND. BTHFN1 == BTHFN(JJW)) CALL DUPLICATE_FILENAME (ERR,BTHFN1)
      IF (JJW /= JW .AND. METFN1 == METFN(JJW)) CALL DUPLICATE_FILENAME (ERR,METFN1)
      IF (JJW /= JW .AND. SNPFN1 == SNPFN(JJW)) CALL DUPLICATE_FILENAME (ERR,SNPFN1)
      IF (JJW /= JW .AND. PRFFN1 == PRFFN(JJW)) CALL DUPLICATE_FILENAME (ERR,PRFFN1)
      IF (JJW /= JW .AND. SPRFN1 == SPRFN(JJW)) CALL DUPLICATE_FILENAME (ERR,SPRFN1)
      IF (JJW /= JW .AND. VPLFN1 == VPLFN(JJW)) CALL DUPLICATE_FILENAME (ERR,VPLFN1)
      IF (JJW /= JW .AND. CPLFN1 == CPLFN(JJW)) CALL DUPLICATE_FILENAME (ERR,CPLFN1)
      IF (JJW /= JW .AND. FLXFN1 == FLXFN(JJW)) CALL DUPLICATE_FILENAME (ERR,FLXFN1)
      IF (PRECIPITATION(JW)) THEN
        IF (JJW /= JW .AND. PREFN1 == PREFN(JJW)) CALL DUPLICATE_FILENAME (ERR,PREFN1)
        IF (JJW /= JW .AND. TPRFN1 == TPRFN(JJW)) CALL DUPLICATE_FILENAME (ERR,TPRFN1)
        IF (CONSTITUENTS) THEN
          IF (JJW /= JW .AND. CPRFN1 == CPRFN(JJW)) CALL DUPLICATE_FILENAME (ERR,CPRFN1)
        END IF
      END IF
    END DO
  END DO
  DO JB=1,NBR
    QINFN1 = QINFN(JB); TINFN1 = TINFN(JB); CINFN1 = CINFN(JB); QOTFN1 = QOTFN(JB); QDTFN1 = QDTFN(JB); TDTFN1 = TDTFN(JB)
    CDTFN1 = CDTFN(JB); PREFN1 = PREFN(JB); TPRFN1 = TPRFN(JB); CPRFN1 = CPRFN(JB); EUHFN1 = EUHFN(JB); TUHFN1 = TUHFN(JB)
    CUHFN1 = CUHFN(JB); EDHFN1 = EDHFN(JB); TDHFN1 = TDHFN(JB); CDHFN1 = CDHFN(JB)
    DO JJB=1,NBR
      IF (JJB /= JB .AND. QINFN1 == QINFN(JJB) .AND. UQ_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,QINFN1)
      IF (JJB /= JB .AND. TINFN1 == TINFN(JJB) .AND. UQ_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,TINFN1)
      IF (JJB /= JB .AND. QOTFN1 == QOTFN(JJB) .AND. DQ_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,QOTFN1)
      IF (JJB /= JB .AND. QDTFN1 == QDTFN(JJB) .AND. DIST_TRIB(JB))       CALL DUPLICATE_FILENAME (ERR,QDTFN1)
      IF (JJB /= JB .AND. TDTFN1 == TDTFN(JJB) .AND. DIST_TRIB(JB))       CALL DUPLICATE_FILENAME (ERR,TDTFN1)
      IF (JJB /= JB .AND. EUHFN1 == EUHFN(JJB) .AND. UH_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,EUHFN1)
      IF (JJB /= JB .AND. TUHFN1 == TUHFN(JJB) .AND. UH_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,TUHFN1)
      IF (JJB /= JB .AND. EDHFN1 == EDHFN(JJB) .AND. DH_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,EDHFN1)
      IF (JJB /= JB .AND. TDHFN1 == TDHFN(JJB) .AND. DH_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,TDHFN1)
      IF (CONSTITUENTS) THEN
        IF (JJB /= JB .AND. CINFN1 == CINFN(JJB) .AND. UQ_EXTERNAL(JB))   CALL DUPLICATE_FILENAME (ERR,CINFN1)
        IF (JJB /= JB .AND. CDTFN1 == CDTFN(JJB) .AND. DIST_TRIB(JB))     CALL DUPLICATE_FILENAME (ERR,CDTFN1)
!        IF (JJB /= JB .AND. CPRFN1 == CPRFN(JJB) .AND. PRECIPITATION(JB)) CALL DUPLICATE_FILENAME (ERR,CPRFN1)   ! SW 1/16/04
        IF (JJB /= JB .AND. CUHFN1 == CUHFN(JJB) .AND. UH_EXTERNAL(JB))   CALL DUPLICATE_FILENAME (ERR,CUHFN1)
        IF (JJB /= JB .AND. CDHFN1 == CDHFN(JJB) .AND. DH_EXTERNAL(JB))   CALL DUPLICATE_FILENAME (ERR,CDHFN1)
      END IF
    END DO
  END DO
  DO JT=1,NTR
    QTRFN1 = QTRFN(JT); TTRFN1 = TTRFN(JT); CTRFN1 = CTRFN(JT)
    DO JJT=1,NTR
      IF (JJT /= JT .AND. QTRFN1 == QTRFN(JJT)) CALL DUPLICATE_FILENAME (ERR,QTRFN1)
      IF (JJT /= JT .AND. TTRFN1 == TTRFN(JJT)) CALL DUPLICATE_FILENAME (ERR,TTRFN1)
      IF (JJT /= JT .AND. CTRFN1 == CTRFN(JJT) .AND. CONSTITUENTS) CALL DUPLICATE_FILENAME (ERR,CTRFN1)
    END DO
  END DO

!***********************************************************************************************************************************
!*                                            Task 5:  Calculate Volume-Area-Elevation Table                                      **
!***********************************************************************************************************************************

  WRITE (WIN,*) 'Volume-area-elevation'
  DO JW=1,NWB
    VOLG(JW) = 0.0
    DO JB=BS(JW),BE(JW)
      IU  = US(JB)
      ID  = DS(JB)
      IF (.NOT.DAM_FLOW(JB)) THEN
        IF (UHS(JB) /= 0 .AND. UHS(JB) > 1) THEN
          DIFF = ABS(ELWS(IU)-ELWS(UHS(JB)))
          IF (DIFF > 1.5) THEN
            CALL WARNINGS
            WRITE (WRN,FMTF2I) 'Water surface elevation difference of',DIFF,' m between segment ',IU,'and segment ',UHS(JB)
          END IF
        END IF
      END IF
      IF (DHS(JB) /= 0 .AND. DHS(JB) > 1) THEN
        DIFF = ABS(ELWS(ID)-ELWS(DHS(JB)))
        IF (DIFF > 1.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF2I) 'Water surface elevation difference of ',DIFF,' m between segment ',ID,'and segment ', DHS(JB)
        END IF
      END IF
      DO I=ID,IU,-1
        XI(I)   = XI(I+1)+(DLX(I)+DLX(I+1))*0.5
        XBR(JB) = XBR(JB)+DLX(I)
      END DO
      XGR(JW) = XGR(JW)+XBR(JB)

!**** Water surface and bottom layers

      DO I=IU-1,ID+1
        IF (ELWS(I) <= EL(KB(I)+1,I)) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Water surface elevation is below bottom elevation at segment ',I
        ELSE IF (ELWS(I) <=  EL(KB(I)+1,I)+0.4*H(KB(I),JW)) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Water surface elevation is close to bottom elevation at segment ',I
        END IF
        IF (I < ID+1) THEN
          DIFF = ABS(ELWS(I)-ELWS(I+1))
          IF (DIFF > 1.5) THEN
            CALL WARNINGS
            WRITE (WRN,FMTF2I) 'Water surface elevation difference of ',DIFF,'m between segment ',I,'and segment ',I+1
          END IF
        END IF
      END DO

!**** Upstream active segment and single layer

      IUT = IU
      DO I=IU,ID
        IF (KB(I)-KTWB(JW) < NL(JB)-1) IUT = I+1
      END DO
      IF (IUT > DS(JB)-1) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'No active segments in branch ',JB
      END IF
      CUS(JB) = IUT
      IF (CUS(JB) /= US(JB)) THEN
        DO JBB=1,NBR
          IF (DHS(JBB) == US(JB)) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Water level too low in branch ',JB,'[CUS(JB)=',CUS(JB),'] connected to branch ',JBB,' by '//         &
                             'downstream head'
          END IF
          IF (UHS(JBB) == US(JB)) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Water level too low in branch ',JB,' [CUS(JB)=',CUS(JB),'[ connected to branch ',JBB,' by '//        &
                             'upstream head'
          END IF
        END DO
      END IF

!**** Branch and grid total volume

      DO I=IU,ID
        VOLB(JB) = VOLB(JB)+DLX(I)*B(KTWB(JW),I)*(H(KTWB(JW),JW)-Z(I))
        VOLG(JW) = VOLG(JW)+DLX(I)*B(KTWB(JW),I)*(H(KTWB(JW),JW)-Z(I))
        DO K=KTWB(JW)+1,KB(I)
          VOLB(JB) = VOLB(JB)+DLX(I)*B(K,I)*H(K,JW)
          VOLG(JW) = VOLG(JW)+DLX(I)*B(K,I)*H(K,JW)
        END DO
      END DO

!**** Branch and grid area and volume by layer

      DO K=KMX-1,2,-1
        NCCBR(K,JB) = NCCBR(K+1,JB)
        CVBR(K,JB)  = CVBR(K+1,JB)
        DO I=IU,ID
          IF (K <= KB(I)) THEN
            SABR(K,JB)  = SABR(K,JB)+B(K,I)        *DLX(I)
            CVBR(K,JB)  = CVBR(K,JB)+B(K,I)*H(K,JW)*DLX(I)
            NCCBR(K,JB) = NCCBR(K,JB)+1
          END IF
        END DO
        SAGR(K,JW)  = SAGR(K,JW) +SABR(K,JB)
        NCCGR(K,JW) = NCCGR(K,JW)+NCCBR(K,JB)
        CVGR(K,JW)  = CVGR(K,JW) +CVBR(K,JB)
      END DO
    END DO
  END DO

! Average depths and widths

  DO K=KMX-1,2,-1
    DO JW=1,NWB
      IF (SAGR(K,JW) /= 0.0) HGR(K,JW) = CVGR(K,JW)/SAGR(K,JW)
      BGR(K,JW) = SAGR(K,JW)/XGR(JW)
      DO JB=BS(JW),BE(JW)
        IF (SABR(K,JB) /= 0.0) HBR(K,JB) = CVBR(K,JB)/SABR(K,JB)
        BBR(K,JB) = SABR(K,JB)/XBR(JB)
      END DO
    END DO
  END DO

! Beginning and ending segment and bottom layer for snapshots

  DO JW=1,NWB
    IEPR(JW) = NISNP(JW)
    DO I=1,IEPR(JW)
      IPR(I,JW) = ISNP(I,JW)
    END DO
    DO I=1,IEPR(JW)
      KEPR(JW) = MAX(KB(IPR(I,JW)),KEPR(JW))
    END DO
  END DO

!***********************************************************************************************************************************
!*                                           Task 6:  Initial Conditions for Simulation                                           **
!***********************************************************************************************************************************

! Temperature and constituents

  WRITE (WIN,*) 'Initial conditions'
  DO JW=1,NWB
    IF (OPEN_VPR(JW)) THEN
      OPEN (VPR(JW),FILE=VPRFN(JW),STATUS='OLD',IOSTAT=IERR)
      IF (IERR == 0) THEN
        READ (VPR(JW),*)
        IF (VERT_TEMP(JW)) READ (VPR(JW),'(/A8/(8X,9F8.0))')        AID, (TVP(K,JW),   K=KTWB(JW),KBMAX(JW))
        IF (CONSTITUENTS) THEN
          DO JC=1,NCT
            IF (VERT_CONC(JC,JW)) READ (VPR(JW),'(/A8/(8X,9F8.0))') AID, (CVP(K,JC,JW),K=KTWB(JW),KBMAX(JW))
          END DO
        END IF
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTA) 'Vertical initial conditions for temperature or constituents, but could not open '//VPRFN(JW)
      END IF
    END IF

!** Longitudinal/vertical initial profiles

    IF (OPEN_LPR(JW)) THEN
      OPEN (LPR(JW),FILE=LPRFN(JW),STATUS='OLD',IOSTAT=IERR)
      IF (IERR == 0) THEN
        READ (LPR(JW),*)
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTA) 'Longitudinal/vertical initial conditions for temperature or constituents, but could not open '//LPRFN(JW)
        WRITE(BRA,'(I0)') JW
        NAME = 'lpr_wb'//TRIM(ADJUSTL(BRA))//'.npt'
        OPEN  (LPR(JW),FILE=NAME,STATUS='UNKNOWN')
        WRITE (LPR(JW),FMTA) 'Longitudinal/vertical initial constituents concentrations'
        IF (LONG_TEMP(JW)) THEN
          DO JB=BS(JW),BE(JW)
            DO I=CUS(JB),DS(JB)
              WRITE(LPR(JW),'(/10A,3X,A,I0/(:8X,9F8.3))') 'TEMP    ',('     T2I',J=1,9),'Segment ',I,(0.0,K=KTWB(JW),KB(I))
            END DO
          END DO
        END IF
     IF(CONSTITUENTS)THEN
        DO JC=1,NCT
          IF (LONG_CONC(JC,JW)) THEN
            DO JB=BS(JW),BE(JW)
              DO I=CUS(JB),DS(JB)
                WRITE(LPR(JW),'(/10A,3X,A,1X,I0/(:8X,9F8.3))') CNAME2(JC),('     C2I',J=1,9),'Segment',I,(0.0,K=KTWB(JW),KB(I))
              END DO
            END DO
          END IF
        END DO
        IF (LONG_SED(JW)) THEN
          DO JB=BS(JW),BE(JW)
            DO I=CUS(JB),DS(JB)
              WRITE(LPR(JW),'(/10A,3X,A,I0/(:8X,9F8.3))') 'SED     ',('   SEDCI',J=1,9),'Segment ',I,(0.0,K=KTWB(JW),KB(I))
            END DO
          END DO
        END IF
        END IF
      END IF
    END IF

    IF (.NOT.OPEN_LPR(JW).AND.CONSTITUENTS.AND.OPEN_LPRC(JW)) THEN
      OPEN (LPR(JW),FILE=LPRFN(JW),STATUS='OLD',IOSTAT=IERR)
      IF (IERR == 0) THEN
        READ (LPR(JW),*)
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTA) 'Longitudinal/vertical initial conditions for constituents, but could not open '//LPRFN(JW)
        WRITE(BRA,'(I0)') JW
        NAME = 'lpr_wb'//TRIM(ADJUSTL(BRA))//'.npt'
        OPEN  (LPR(JW),FILE=NAME,STATUS='UNKNOWN')
        WRITE (LPR(JW),FMTA) 'Longitudinal/vertical initial constituents concentrations'
        IF (LONG_TEMP(JW)) THEN
          DO JB=BS(JW),BE(JW)
            DO I=CUS(JB),DS(JB)
              WRITE(LPR(JW),'(/10A,3X,A,I0/(:8X,9F8.3))') 'TEMP    ',('     T2I',J=1,9),'Segment ',I,(0.0,K=KTWB(JW),KB(I))
            END DO
          END DO
        END IF
        DO JC=1,NCT
          IF (LONG_CONC(JC,JW)) THEN
            DO JB=BS(JW),BE(JW)
              DO I=CUS(JB),DS(JB)
                WRITE(LPR(JW),'(/10A,3X,A,1X,I0/(:8X,9F8.3))') CNAME2(JC),('     C2I',J=1,9),'Segment',I,(0.0,K=KTWB(JW),KB(I))
              END DO
            END DO
          END IF
        END DO
        IF (LONG_SED(JW)) THEN
          DO JB=BS(JW),BE(JW)
            DO I=CUS(JB),DS(JB)
              WRITE(LPR(JW),'(/10A,3X,A,I0/(:8X,9F8.3))') 'SED     ',('   SEDCI',J=1,9),'Segment ',I,(0.0,K=KTWB(JW),KB(I))
            END DO
          END DO
        END IF
      END IF
    END IF

    DO JB=BS(JW),BE(JW)

!**** Ice thickness

      ICETH(CUS(JB)-1:DS(JB)+1) = ICEI(JW)

!**** Temperature

      DO I=CUS(JB),DS(JB)
        IF (LONG_TEMP(JW)) READ (LPR(JW),'(/A8/(8X,9F8.0))',IOSTAT=IERR) AID, (T2(K,I),K=KTWB(JW),KB(I))
        IF (IERR == 0) THEN
          DO K=KTWB(JW),KB(I)
            IF (ISO_TEMP(JW))  T2(K,I) = T2I(JW)
            IF (VERT_TEMP(JW)) T2(K,I) = TVP(K,JW)
          END DO
        END IF
      END DO
    END DO

!**** Constituents

    DO JB=BS(JW),BE(JW)
      DO JC=1,NAC
        DO I=CUS(JB),DS(JB)
          IF (LONG_CONC(CN(JC),JW)) READ (LPR(JW),'(/A8/(8X,9F8.0))',IOSTAT=IERR) AID,(C2(K,I,CN(JC)),K=KTWB(JW),KB(I))
          IF (IERR == 0) THEN
            DO K=KTWB(JW),KB(I)
              IF (ISO_CONC(CN(JC),JW))  C2(K,I,CN(JC)) = C2IWB(CN(JC),JW)
              IF (VERT_CONC(CN(JC),JW)) C2(K,I,CN(JC)) = CVP(K,CN(JC),JW)
            END DO
          END IF
        END DO
      END DO
    END DO
  END DO

!***********************************************************************************************************************************
!*                                                     Task 7:  Output section                                                    **
!***********************************************************************************************************************************

! Initial input and conditions

  WRITE (WIN,*) 'Preprocessor output'
  WRITE (INI,'(1X,A/)') 'CE-QUAL-W2 preprocessor - V3.6'
  WRITE (INI,11130) (TITLE(J),J=1,10)                                                                                  !TC 07/23/02
  WRITE (INI,11140) 'Time Control', TMSTRT, TMEND, YEAR
  WRITE (INI,11150)  NDLT, DLTMIN
  WRITE (INI,11160) (DLTD(J),J=1,NDLT)
  WRITE (INI,11170) (DLTMAX(J),J=1,NDLT)
  WRITE (INI,11180) (DLTF(J),J=1,NDLT)
  WRITE (INI,11190) 'Timestep limitation', (JW,VISC(JW)(6:8), CELC(JW)(6:8), JW=1,NWB)
  WRITE (INI,'(/1X,A)') 'Initial Conditions'
  DO JW=1,NWB
    WRITE (INI,FMTI) '   Waterbody ',JW
    IF (ISO_TEMP(JW))  WRITE (INI,'(5X,A25,F5.1,A3)') 'Temperature      [T2I] = ',T2I(JW),' øC'
    IF (VERT_TEMP(JW)) WRITE (INI,'(5X,A52)')         'Temperature      [T2I] = Downstream vertical profile'
    WRITE (INI,FMTA) '     Water type    [WTYPEC] = '//WTYPEC(JW)(4:8)//' water'
    WRITE (INI,FMTF) '     Ice thickness   [ICEI] = ',ICEI(JW),' m'
  END DO
  WRITE (INI,11250) 'Calculations', (JW,EVC(JW)(6:8),PRC(JW)(6:8),VBC(JW)(6:8),EBC(JW)(6:8),MBC(JW)(6:8),PQC(JW)(6:8),             &
                                     WINDC(JW)(6:8), QINC(JW)(6:8),QOUTC(JW)(6:8),HEATC(JW)(6:8),SLHTC(JW)(5:8),JW=1,NWB)
  DO JW=1,NWB
    WRITE (INI,11262) JW,SROC(JW)(6:8)
    WRITE (INI,11263) AFW(JW)
    WRITE (INI,11264) BFW(JW)
    WRITE (INI,11265) CFW(JW)
    WRITE (INI,11266) WINDH(JW)
    WRITE (INI,11267) RHEVC(JW)(6:8)
  END DO
  WRITE (INI,11270)  (JB,QINIC(JB)(6:8), DTRIC(JB)(6:8), HDIC(JB)(6:8), JB=1,NBR)
  WRITE (INI,11271)  (JW,METIC(JW)(6:8), JW =1,NWB)
  IF (TRIBUTARIES) WRITE (INI,11272)  (JT,TRIC(JT)(6:8),  JT =1,NTR)
  DO JB=1,NBR
    IF (NSTR(JB) /= 0) WRITE (INI,11273)  JB,(JS,STRIC(JS,JB)(6:8),JS=1,NSTR(JB))
  END DO
  WRITE (INI,11274) (JWD,WDIC(JWD)(6:8),JWD=1,NWD)
  WRITE (INI,FMTA) ' Meteorological Parameters'
  DO JW=1,NWB
    WRITE (INI,11300)  JW, LAT(JW), LONG(JW)
  END DO
  IEGR = 1
  WRITE (INI,11330)
  DO WHILE (IEGR < IMX-1)
    IBGR = IEGR+1
    IEGR = IEGR+19
    IF (IEGR > IMX-1) IEGR = IMX-1
    WRITE (INI,11340) (I,I=IBGR,IEGR)
    WRITE (INI,11350) (PHI0(I),I=IBGR,IEGR)
  END DO
  WRITE (INI,11360) (JW,SLTRC(JW),THETA(JW),JW=1,NWB)
  WRITE (INI,11365)
  DO JW=1,NWB
    WRITE (INI,11400) JW
    WRITE (INI,11370) AX(JW),DX(JW),TSED(JW),DEG,CBHE(JW),DEG
  END DO
  WRITE (INI,11375)
  DO JW=1,NWB
    WRITE (INI,11380) JW,ICEC(JW)(6:8),SLICEC(JW),ALBEDO(JW),HWI(JW),BETAI(JW),GAMMAI(JW)
  END DO
  WRITE (INI,11385)
  DO JW=1,NWB
    WRITE (INI,11390) JW,(HNAME(J),HPRWBC(J,JW)(6:8),J=1,NHY)
  END DO
  DO JW=1,NWB
    WRITE (INI,11400)  JW
    WRITE (INI,11410)  SNPC(JW)(6:8)
    IF (SNPC(JW)(6:8) == ' ON') WRITE (INI,11420)  NSNP(JW)
    IF (SNPC(JW)(6:8) == ' ON') WRITE (INI,11430) (SNPD(J,JW),J=1,NSNP(JW))
    IF (SNPC(JW)(6:8) == ' ON') WRITE (INI,11440) (SNPF(J,JW),J=1,NSNP(JW))
    WRITE (INI,11450)  SCRC(JW)(6:8)
    IF (SCRC(JW)(6:8) == ' ON') WRITE (INI,11455)  NSCR(JW)
    IF (SCRC(JW)(6:8) == ' ON') WRITE (INI,11460) (SCRD(J,JW),J=1,NSCR(JW))
    IF (SCRC(JW)(6:8) == ' ON') WRITE (INI,11470) (SCRF(J,JW),J=1,NSCR(JW))
    WRITE (INI,11520)  FLXC(JW)(6:8)
    IF (FLXC(JW)(6:8) == ' ON') WRITE (INI,11530)  NFLX(JW)
    IF (FLXC(JW)(6:8) == ' ON') WRITE (INI,11540) (FLXD(J,JW),J=1,NFLX(JW))
    IF (FLXC(JW)(6:8) == ' ON') WRITE (INI,11550) (FLXF(J,JW),J=1,NFLX(JW))
    WRITE (INI,11560)  VPLC(JW)(6:8)
    IF (VPLC(JW)(6:8) == ' ON') WRITE (INI,11570)  NVPL(JW)
    IF (VPLC(JW)(6:8) == ' ON') WRITE (INI,11580) (VPLD(J,JW),J=1,NVPL(JW))
    IF (VPLC(JW)(6:8) == ' ON') WRITE (INI,11590) (VPLF(J,JW),J=1,NVPL(JW))
    WRITE (INI,11600)  PRFC(JW)(6:8)
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11610)  NPRF(JW)
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11620)  NIPRF(JW)
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11630) (IPRF(I,JW),I=1,NIPRF(JW))
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11640) (PRFD(J,JW),J=1,NPRF(JW))
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11650) (PRFF(J,JW),J=1,NPRF(JW))
    WRITE (INI,11660)  SPRC(JW)(6:8)
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11670)  NSPR(JW)
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11680)  NISPR(JW)
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11690) (ISPR(I,JW),I=1,NISPR(JW))
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11700) (SPRD(J,JW),J=1,NSPR(JW))
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11710) (SPRF(J,JW),J=1,NSPR(JW))
    WRITE (INI,11720)  CPLC(JW)(6:8)
    IF (CPLC(JW)(6:8) == ' ON') WRITE (INI,11730)  NCPL(JW)
    IF (CPLC(JW)(6:8) == ' ON') WRITE (INI,11740) (CPLD(J,JW),J=1,NCPL(JW))
    IF (CPLC(JW)(6:8) == ' ON') WRITE (INI,11750) (CPLF(J,JW),J=1,NCPL(JW))
  END DO
  WRITE (INI,11480)  TSRC(6:8)
  IF (TSRC(6:8) == ' ON') WRITE (INI,11490)  NTSR
  IF (TSRC(6:8) == ' ON') WRITE (INI,11500) (TSRD(J),J=1,NTSR)
  IF (TSRC(6:8) == ' ON') WRITE (INI,11510) (TSRF(J),J=1,NTSR)
  WRITE (INI,11760)  RSOC(6:8), RSIC(6:8)
  IF (RSOC(6:8) == ' ON') WRITE (INI,11770)  NRSO
  IF (RSOC(6:8) == ' ON') WRITE (INI,11780) (RSOD(J),J =1,NRSO)
  IF (RSOC(6:8) == ' ON') WRITE (INI,11790) (RSOF(J),J =1,NRSO)
  WRITE (INI,11800) (JB,NSTR(JB),JB=1,NBR)
  DO JB=1,NBR
    IF (NSTR(JB) > 0)  THEN
      WRITE (INI,11805)  JB
      WRITE (INI,11810) (JS,SINKC(JS,JB),WSTR(JS,JB),ESTR(JS,JB),KBSTR(JS,JB),JS=1,NSTR(JB))
    END IF
  END DO
  WRITE (INI,11840)  NWD, (IWD(JW),JW=1,NWD)
  WRITE (INI,11860)  NTR, (ITR(JT),JT=1,NTR)
  IF (TRIBUTARIES) THEN
    WRITE (INI,11870) (PTRC(JT),     JT=1,NTR)
    WRITE (INI,11880) (ETTR(JT),     JT=1,NTR)
    WRITE (INI,11890) (EBTR(JT),     JT=1,NTR)
  END IF
  IF (SPILLWAYS) THEN
    WRITE (INI,11821)  NSP,(ESP(JS),JS=1,NSP)
    WRITE (INI,11822) (IUSP(JS),    JS=1,NSP)
    WRITE (INI,11823) (IDSP(JS),    JS=1,NSP)
    WRITE (INI,11824) (PUSPC(JS),   JS=1,NSP)
    WRITE (INI,11825) (PDSPC(JS),   JS=1,NSP)
    WRITE (INI,11826) (ETUSP(JS),   JS=1,NSP)
    WRITE (INI,11827) (ETDSP(JS),   JS=1,NSP)
    WRITE (INI,11828) (EBUSP(JS),   JS=1,NSP)
    WRITE (INI,11829) (EBDSP(JS),   JS=1,NSP)
    WRITE (INI,11830) (A1SP(JS),B1SP(JS),A2SP(JS),B2SP(JS),JS=1,NSP)
  END IF
  IF (GATES) THEN
    WRITE (INI,11891)  NGT,(EGT(JG),JG=1,NGT)
    WRITE (INI,11892) (IUGT(JG),    JG=1,NGT)
    WRITE (INI,11893) (IDGT(JG),    JG=1,NGT)
    WRITE (INI,11894) (PUGTC(JG),   JG=1,NGT)
    WRITE (INI,11895) (PDGTC(JG),   JG=1,NGT)
    WRITE (INI,11896) (ETUGT(JG),   JG=1,NGT)
    WRITE (INI,11897) (ETDGT(JG),   JG=1,NGT)
    WRITE (INI,11898) (EBUGT(JG),   JG=1,NGT)
    WRITE (INI,11899) (EBDGT(JG),   JG=1,NGT)
    WRITE (INI,11820) (A1GT(JG),B1GT(JG),G1GT(JG),A2GT(JG),B2GT(JG),G2GT(JG),JG=1,NGT)
  END IF
  IF (PIPES) THEN
    WRITE (INI,11831)  NPI,(EUPI(JP),JP=1,NPI)
    WRITE (INI,11841) (EDPI(JP),     JP=1,NPI)
    WRITE (INI,11842) (WPI(JP),      JP=1,NPI)
    WRITE (INI,11843) (DLXPI(JP),    JP=1,NPI)
    WRITE (INI,11844) (FPI(JP),      JP=1,NPI)
    WRITE (INI,11845) (FMINPI(JP),   JP=1,NPI)
    WRITE (INI,11832) (IUPI(JP),     JP=1,NPI)
    WRITE (INI,11833) (IDPI(JP),     JP=1,NPI)
    WRITE (INI,11834) (PUPIC(JP),    JP=1,NPI)
    WRITE (INI,11835) (PDPIC(JP),    JP=1,NPI)
    WRITE (INI,11836) (ETUPI(JP),    JP=1,NPI)
    WRITE (INI,11837) (ETDPI(JP),    JP=1,NPI)
    WRITE (INI,11838) (EBUPI(JP),    JP=1,NPI)
    WRITE (INI,11839) (EBDPI(JP),    JP=1,NPI)
  END IF
  IF (PUMPS) THEN
    WRITE (INI,11851)  NPU,(EPU(JP),JP=1,NPU)
    WRITE (INI,11852) (IUPU(JP),    JP=1,NPU)
    WRITE (INI,11853) (IDPU(JP),    JP=1,NPU)
    WRITE (INI,11854) (STRTPU(JP),  JP=1,NPU)
    WRITE (INI,11855) (ENDPU(JP),   JP=1,NPU)
    WRITE (INI,11856) (EONPU(JP),   JP=1,NPU)
    WRITE (INI,11857) (EOFFPU(JP),  JP=1,NPU)
    WRITE (INI,11858) (QPU(JP),     JP=1,NPU)
    WRITE (INI,11859) (PPUC(JP),    JP=1,NPU)
    WRITE (INI,11861) (ETPU(JP),    JP=1,NPU)
    WRITE (INI,11862) (EBPU(JP),    JP=1,NPU)
    WRITE (INI,11863) (KTPU(JP),    JP=1,NPU)
    WRITE (INI,11864) (KBPU(JP),    JP=1,NPU)
  END IF
  WRITE (INI,11900)
  WRITE (INI,11910) (JB,DTRC(JB)(6:8),JB=1,NBR)
  WRITE (INI,11930) CONFN, RSIFN, QWDFN
  DO JW=1,NWB
    WRITE (INI,11940) JW, BTHFN(JW), METFN(JW), VPRFN(JW), LPRFN(JW)
  END DO
  DO JB=1,NBR
    WRITE (INI,11950) JB, QINFN(JB), TINFN(JB), CINFN(JB), QOTFN(JB), QDTFN(JB), TDTFN(JB), CDTFN(JB), PREFN(JB), TPRFN(JB),       &
                          CPRFN(JB), EUHFN(JB), TUHFN(JB), CUHFN(JB), EDHFN(JB), TDHFN(JB), CDHFN(JB)
  END DO
  WRITE (INI,11970) (JT,QTRFN(JT),TTRFN(JT),CTRFN(JT),JT=1,NTR)
  WRITE (INI,11980)  ERRFN, WRNFN, TSRFN, WDOFN
  DO JW=1,NWB
    WRITE (INI,11990) JW, SNPFN(JW), FLXFN(JW), PRFFN(JW), VPLFN(JW), CPLFN(JW)  !put in shade, wsc, graph, etc. "to do"
  END DO
  IF (CONSTITUENTS) THEN
    WRITE (INI,12000) CCC(6:8), LIMC(6:8), CUF
    DO JW=1,NWB
      WRITE (INI,12030) JW
      WRITE (INI,12020)
      WRITE (INI,12040) (CNAME1(JC),CAC(JC)(6:8),C2IWB(JC,JW),CFWBC(JC,JW)(6:8),CPRWBC(JC,JW)(6:8),JC=1,NCT),'Sediments          ',&
                         SEDC(JW)(6:8), SEDCI(JW),'OFF', PRNSC(JW)(6:8)
      DO JB=BS(JW),BE(JW)
        WRITE (INI,12045) JB,(CNAME1(JC),CINBRC(JC,JB)(6:8),CDTBRC(JC,JB)(6:8),CPRBRC(JC,JB)(6:8),JC=1,NCT)
      END DO
      WRITE (INI,12050)
      WRITE (INI,12060)
      DO JD=1,NDC
        WRITE (INI,12070) CDNAME(JD),CDWBC(JD,JW)(6:8)
      END DO
    END DO
    DO JT=1,NTR
      WRITE (INI,12046) JT,(CNAME1(JC),CINTRC(JC,JT)(6:8),JC=1,NCT)
    END DO
    DO JW=1,NWB
      WRITE (INI,12030)  JW
      WRITE (INI,12080) 'Constituent Rates'
      WRITE (INI,12100) (CNAME1(JG+1), CGQ10(JG), CG0DK(JG), CG1DK(JG), CGS(JG),JG=1,NGC)
      WRITE (INI,12090) (SSS(NS),NS=1,NSS)
      WRITE (INI,12110)  LDOMDK(JW), LRDDK(JW)
      WRITE (INI,12120)  RDOMDK(JW)
      WRITE (INI,12130)  LPOMDK(JW), LRPDK(JW), POMS(JW)
      WRITE (INI,12150)  RPOMDK(JW)
      DO JA=1,NAL
        WRITE (INI,12160) JA, AG(JA), AM(JA), AE(JA), AR(JA), AS(JA), ALGP(JA), ALGN(JA), ALGC(JA), ALGSI(JA), ACHLA(JA), APOM(JA)
      END DO
      WRITE (INI,12180)  PO4R(JW)
      WRITE (INI,12190)  NH4DK(JW), NH4R(JW)
      WRITE (INI,12200)  NO3DK(JW),NO3S(JW),FNO3SED(JW)
      WRITE (INI,12210)  DSIR(JW), PSIDK(JW), PSIS(JW)
      WRITE (INI,12220)  SEDDK(JW)
      WRITE (INI,12230)  FES(JW), FER(JW)
      WRITE (INI,12240) (SOD(I),I=US(BS(JW)),DS(BE(JW)))
      WRITE (INI,12250)  FSOD(JW), FSED(JW)
      WRITE (INI,12255)  SEDS(jw), SEDBR(jw)
    END DO
    WRITE (INI,12330) 'Upper Temperature Bounds'
    DO JW=1,NWB
      WRITE (INI,12290)  NH4T1(JW), NH4T2(JW)
      WRITE (INI,12300)  NO3T1(JW), NO3T2(JW)
      WRITE (INI,12310)  OMT1(JW),  OMT2(JW)
      WRITE (INI,12320)  SODT1(JW), SODT2(JW)
    END DO
    DO JA=1,NAL
      WRITE (INI,12280) JA, AT1(JA), AT2(JA)
    END DO
    WRITE (INI,12270) 'Lower Temperature Bounds'
    DO JA=1,NAL
      WRITE (INI,12340) JA, AT3(JA), AT4(JA)
    END DO
    DO JW=1,NWB
      WRITE (INI,12350) 'Stoichiometric Equivalence', O2NH4(JW), O2OM(JW)
      WRITE (INI,12355) (JA,O2AR(JA),O2AG(JA),JA=1,NAL)
      WRITE (INI,12360)     ORGC(JW),ORGP(JW), ORGN(JW), ORGSI(JW)
      WRITE (INI,12380) 'Half Saturation', (JA,AHSP(JA), AHSN(JA),  AHSSI(JA),JA=1,NAL)
      WRITE (INI,12390) 'Light',               BETA(JW), EXH2O(JW), EXSS(JW), EXOM(JW)
      WRITE (INI,12400) (JA,EXA(JA), JA=1,NAL)
      WRITE (INI,12410) (JA,ASAT(JA),JA=1,NAL)
      WRITE (INI,12420) 'Diffusion',                DMO2,      DMCO2
      WRITE (INI,12430) 'Partitioning Coefficients',PARTP(JW), PARTSI(JW)
      WRITE (INI,12440) 'Miscellaneous Constants',  O2LIM,     CO2R(JW)
    END DO
    IF (NBOD > 0) WRITE (INI,12370)  (JB, KBOD(JB), TBOD(JB), RBOD(JB), JB=1,NBOD)                                     !TC 06/18/02

!** Inflow constituent statistics                                                                                      !SW 01/07/01

    WRITE (INI,'(/A)') '  Inflow Constituent Statistics'
    DO JW=1,NWB
      DO JB=BS(JW),BE(JW)
        IF (UQ_EXTERNAL(JB)) THEN
          WRITE (INI,FMTI) '    Branch ',JB
          WRITE (INI,FMTA)    '       Constituent name        Average     Maximum     Minimum ApproxLoading(kg/d)'
          DO JC=1,NACIN(JB)
            WRITE (INI,'(7X,A,3(2X,F10.3),2x,e12.5)') CNAME1(INCN(JC,JB))(1:19), CINAVG(INCN(JC,JB),JB), CINMAX(INCN(JC,JB),JB),            &
                                             CINMIN(INCN(JC,JB),JB),cinload(INCN(JC,JB),JB)
          END DO
        END IF
      END DO
    END DO

!** Distributary Tributary constituent statistics           - need to include average and max/min

    DO JB=1,NBR
      if(dist_trib(jb))then
      WRITE (INI,FMTI) '    Dist Tributary ',JB
      WRITE (INI,FMTA)    '       Constituent name     ApproxLoading(kg/day)'
      DO JC=1,NACDT(JB)
        WRITE (INI,'(7X,A,2x,e12.5,3(2X,F10.3))') CNAME1(INCDT(JC,JB))(1:19),cdtload(incdt(jc,jb),jb)            !CTRAVG(TRCN(JC,JT),JT), CTRMAX(TRCN(JC,JT),JT),CTRMIN(TRCN(JC,JT),JT)
      END DO
      endif
    END DO

!** Tributary constituent statistics                                                                                   !SW 01/07/01

    DO JT=1,NTR
      WRITE (INI,FMTI) '     Tributary ',JT
      WRITE (INI,FMTA)    '       Constituent name        Average     Maximum     Minimum ApproxLoading(kg/d)'
      DO JC=1,NACTR(JT)
        WRITE (INI,'(7X,A,3(2X,F10.3), 2x,e12.5)') CNAME1(TRCN(JC,JT))(1:19), CTRAVG(TRCN(JC,JT),JT), CTRMAX(TRCN(JC,JT),JT),                &
                                         CTRMIN(TRCN(JC,JT),JT),ctrload(trCN(JC,Jt),Jt)
      END DO
    END DO


!** Overall Approx Loading for all inflows, tribs, dist tribs
      WRITE (INI,FMTA) '    Overall approx loading for all inflows, tributaries, and distributed tributaries'
      WRITE (INI,FMTA)    '       Constituent name      ApproxLoading(kg/d)'
DO J=1,NCT
  xx=0.0
    do jb=1,nbr
      do jc=1,nacin(jb)
        IF(incn(jc,jb) == j)then
        xx=xx+cinload(incn(jc,jb),jb)
        endif
      enddo
      do jc=1,nacdt(jb)
        IF(incdt(jc,jb) == j)then
        xx=xx+cdtload(incdt(jc,jb),jb)
        endif
      enddo
    enddo
    do jt=1,ntr
        do jc=1,nactr(jt)
        IF(trcn(jc,jt) == j)then
        xx=xx+ctrload(trcn(jc,jt),jt)
        endif
        enddo
    enddo
   if(xx>0.0)WRITE (INI,'(7X,A,2x,e12.5)') CNAME1(j)(1:19),xx
 enddo

END IF

  WRITE (INI,'(/1X,A)') 'Water Balance Summary'
  DO JW=1,NWB
    WRITE (INI,12460) JW, QINAVW(JW), QINMXW(JW), QOTAVW(JW), QOTMXW(JW)
    DO JB=BS(JW),BE(JW)
      WRITE (INI,12470) JB, QINAVB(JB),QINMXB(JB),QINAV(JB), QINMX(JB), QTRAVB(JB),QTRMXB(JB),QDTAV(JB),QDTMX(JB),QPRAV(JB),       &
                            QPRMX(JB), QSTRAV(JB),QSTRMX(JB),QWDAVB(JB),QWDMXB(JB)
    END DO
  END DO

  WRITE (INI,12480) 'Geometry', IMX, KMX, NBR
  DO JW=1,NWB
    KBGR = 0
    KEGR = 0
    WRITE (INI,12490) JW, US(BS(JW))-1, DS(BE(JW))+1, BS(JW), BE(JW), ELBOT(JW), KTWB(JW)
    DO WHILE (KEGR < KMX)
      KBGR = KEGR+1
      KEGR = MIN(KEGR+19,KMX)
      WRITE (INI,12500) (K,K=KBGR,KEGR)
      WRITE (INI,12510) (H(K,JW),K=KBGR,KEGR)
    END DO
    DO JB=BS(JW),BE(JW)
      WRITE (INI,12520) JB, US(JB), DS(JB), UHS(JB), DHS(JB)
    END DO
  END DO

! Volume-area-elevation table

  DO JW=1,NWB
    IF (NBR > 1) THEN
      DO JB=BS(JW),BE(JW)
        WRITE (INI,12600) VOLB(JB)
        WRITE (INI,12620) JB
        WRITE (INI,12640)
        DO K=2,KMX-1
          IF (K == KTWB(JW)) THEN
            WRITE (INI,12650) K, EL(K,DS(JBDN(JW))), SABR(K,JB)/1.E6, CVBR(K,JB)/1.E6, NCCBR(K,JB), HBR(K,JB), BBR(K,JB)
          ELSE
            WRITE (INI,12655) K, EL(K,DS(JBDN(JW))), SABR(K,JB)/1.E6, CVBR(K,JB)/1.E6, NCCBR(K,JB), HBR(K,JB), BBR(K,JB)
          END IF
        END DO
      END DO
    END IF
    WRITE (INI,12670) JW, VOLG(JW)
    WRITE (INI,12690) JW
    WRITE (INI,12630)
    DO K=2,KMX-1
      IF (QOTAVW(JW) > 0.0) THEN
        IF (K == KTWB(JW)) THEN
          WRITE (INI,12650) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW), CVGR(K,JW) &
                            /(QOTAVW(JW)*86400.0)
        ELSE
          WRITE (INI,12655) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW), CVGR(K,JW) &
                            /(QOTAVW(JW)*86400.0)
        END IF
      ELSE IF (QINAVW(JW) > 0.0) THEN
        IF (K == KTWB(JW)) THEN
          WRITE (INI,12650) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW), CVGR(K,JW) &
                            /(QINAVW(JW)*86400.0)
        ELSE
          WRITE (INI,12655) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW), CVGR(K,JW) &
                            /(QINAVW(JW)*86400.0)
        END IF
      ELSE
        IF (K == KTWB(JW)) THEN
          WRITE (INI,12650) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW)
        ELSE
          WRITE (INI,12655) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW)
        END IF
      END IF
    END DO

!** Bathymetry

    IEGR = US(BS(JW))-1
    WRITE (INI,12530) JW
    DO WHILE (IEGR < DS(BE(JW)))
      IBGR = IEGR+1
      IEGR = IEGR+18
      IF (IEGR > DS(BE(JW))) IEGR = DS(BE(JW))
      WRITE (INI,'(/7X,21I10)') (I,I=IBGR,IEGR)
      WRITE (INI,'(7X,21I10)')  (INT(XI(I)),I=IBGR,IEGR)
      DO K=1,KMX
        DO I=IBGR,IEGR
          WRITE (BK(I),'(F8.0)') B(K,I)
        END DO
        IF (K == KTWB(JW)) THEN
          WRITE (INI,'(1X,I3,1X,"KT",1X,25A10)') K,(BK(I),I=IBGR,IEGR)
        ELSE
          WRITE (INI,'(1X,I3,4X,25A10)')         K,(BK(I),I=IBGR,IEGR)
        END IF
      END DO
    END DO
    CALL PRINT_INITIAL_GRID (KTWB(JW))
  END DO
  CLOSE (INI,ERR=10)

! Close error/warning files

10 CONTINUE
  IF (DELETE_WRN) THEN
    CLOSE (WRN,STATUS='DELETE')
    WRITE (WIN,'(/A)') 'Number of warnings = 0'
  ELSE
    CLOSE (WRN)
    WRITE (WIN,'(/A,I0,A)') 'Number of warnings = ',NWRN,' - see the file pre.wrn for more information'
  END IF
  IF (DELETE_ERR) THEN
    CLOSE (ERR,STATUS='DELETE')
    WRITE (WIN,FMTA)  'Number of errors   = 0'
  ELSE
    CLOSE (ERR)
    WRITE (WIN,FMTI)  'Number of errors   = ',NERR,' - see the file pre.err for more information'
  END IF
  WRITE (WIN,FMTA)
  WRITE (WIN,FMTA) 'Normal termination'
  GO TO 99999

! Error message for input cards

400 CONTINUE
  CALL INPUT_ERROR (UNIT,ERR,AID,'Either illegal value or incorrect card somewhere in the following cards')


!**** Output formats

11130 FORMAT (1X,A)
11140 FORMAT (/1X,A/                                                                                                               &
              3X,'Starting time (Julian day) [TMSTRT] =',F8.2/                                                                     &
              3X,'Ending time (Julian day)    [TMEND] =',F8.2/                                                                     &
              3X,'Year                         [YEAR] =',I8)
11150 FORMAT (3X,'# Timestep intervals         [NDLT] =',I8/                                                                       &
              3X,'Minimum timestep (sec)     [DLTMIN] =',F8.1)
11160 FORMAT (3X,'Timestep day (Julian day)    [DLTD] =',F8.1)
11170 FORMAT (3X,'Maximum timestep (sec)     [DLTMAX] =',7F8.1,' sec')
11180 FORMAT (3X,'Fraction of timestep         [DLTF] =',7F8.2)
11190 FORMAT (3X,A/                                                                                                                &
             (5X,'Waterbody ',I0/                                                                                                  &
              5X,'  Vertical eddy viscosity    [VISC] = ',A/                                                                       &
              5X,'  Internal gravity wave      [CELC] = ',A))
11250 FORMAT (/1X,A/                                                                                                               &
             (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Evaporation     [EVC] = ',A/                                                                                     &
              5X,'Precipitation   [PRC] = ',A/                                                                                     &
              5X,'Volume balance  [VBC] = ',A/                                                                                     &
              5X,'Energy balance  [EBC] = ',A/                                                                                     &
              5X,'Mass balance    [MBC] = ',A/                                                                                     &
              5X,'Place inflows   [PQC] = ',A/                                                                                     &
              5X,'Wind          [WINDC] = ',A/                                                                                     &
              5X,'Inflow         [QINC] = ',A/                                                                                     &
              5X,'Outflow       [QOUTC] = ',A/                                                                                     &
              5X,'Heat exchange [HEATC] = ',A/                                                                                     &
              5X,'Heat exchange [SLHTC] = ',A))
11262 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'read radiation               [SROC] =',7(5X,A):/                                                                 &
             (42X,7(5X,A)))
11263 FORMAT (5X,'wind function coefficient a   [AFW] =',7F8.2:/                                                                   &
             (42X,7(F8.2)))
11264 FORMAT (5X,'wind function coefficient b   [BFW] =',7F8.2:/                                                                   &
             (42X,7(F8.2)))
11265 FORMAT (5X,'wind function coefficient c   [CFW] =',7F8.2:/                                                                   &
             (42X,7(F8.2)))
11266 FORMAT (5X,'wind height                 [WINDH] =',7F8.2:/                                                                   &
             (42X,7(F8.2)))
11267 FORMAT (5X,'Ryan-Harleman evaporation   [RHEVC] =',7(5X,A):/                                                                 &
             (42X,7(5X,A)))
11270 FORMAT (/1X,'Input Interpolations'/                                                                                          &
             (3X,'Branch ',I0/                                                                                                     &
              5X,'Inflow                [QINIC] = ',A/                                                                             &
              5X,'Distributed tributary [DTRIC] = ',A/                                                                             &
              5X,'Head boundary          [HDIC] = ',A))
11271 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Meteorology           [METIC] = ',A)
11272 FORMAT (:3X,'Tributary ',I0,'              [TRIC] = ',A)
11273 FORMAT (3X,'Branch ',I0/                                                                                                     &
             (5X,'Structure ',I0,'           [STRIC] = ',A))
11274 FORMAT (:5X,'Withdrawal ',I0,'            [WDIC] = ',A)
11300 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Latitude                  [LAT] =',F8.2/                                                                         &
              5X,'Longitude                [LONG] =',F8.2)
11330 FORMAT (3X,'Axis orientation')
11340 FORMAT (5X,'Segment #    ',(1X,21I5))
11350 FORMAT (5X,'[PHI0] (rads)',(1X,21F5.2))
11360 FORMAT (/1X,'Transport Solution'/                                                                                            &
             (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Transport [SLTRC] = ',A/                                                                                         &
              5X,'Theta     [THETA] =',F9.2))
11365 FORMAT (/1X,'Hydraulic coefficients')
11370 FORMAT (3X,'Longitudinal eddy viscosity           [AX] =',F9.2,' m^2/sec'/                                                   &
              3X,'Longitudinal eddy diffusivity         [DX] =',F9.2,' m^2/sec'/                                                   &
              3X,'Sediment temperature                [TSED] =',F9.2,1X,A/                                                         &
              3X,'Coefficient of bottom heat exchange [CBHE] =',1PE9.2,1X,A,' m/s')
11375 FORMAT (/1X,'Ice cover')
11380 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Ice calculations          [ICEC] = ',5X,A/                                                                       &
              5X,'Solution                [SLICEC] = ',A/                                                                          &
              5X,'Albedo                  [ALBEDO] = ',F8.2/                                                                       &
              5X,'Ice-water heat exchange    [HWI] = ',F8.2/                                                                       &
              5X,'Light absorption         [BETAI] = ',F8.2/                                                                       &
              5X,'Light decay             [GAMMAI] = ',F8.2)
11385 FORMAT (/1X,'Output Control')
11390 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
             (5X,A41,' = ',A))
11400 FORMAT (3X,'Waterbody ',I0)
11410 FORMAT (5X,'Snapshot         [SNPC] = ',A)
11420 FORMAT (7X,'Number of time intervals [NSNP] =',I7)
11430 FORMAT (7X,'Date  (Julian day)       [SNPD] =',11F7.2)
11440 FORMAT (7X,'Frequency  (days)        [SNPF] =',11F7.2)
11450 FORMAT (5X,'Screen           [SCRC] = ',A)
11455 FORMAT (7X,'Number of time intervals [NSCR] =',I7)
11460 FORMAT (7X,'Date  (Julian day)       [SCRD] =',11F7.2)
11470 FORMAT (7X,'Frequency  (days)        [SCRF] =',11F7.2)
11480 FORMAT (5X,'Time series      [TSRC] = ',A)
11490 FORMAT (7X,'Number of time intervals [NTSR] =',I7)
11500 FORMAT (7X,'Date  (Julian day)       [TSRD] =',11F7.2)
11510 FORMAT (7X,'Frequency  (days)        [TSRF] =',11F7.2)
11520 FORMAT (5X,'Fluxes           [FLXC] = ',A)
11530 FORMAT (7X,'Number of time intervals [NFLX] =',I7)
11540 FORMAT (7X,'Date  (Julian day)       [FLXD] =',11F7.2)
11550 FORMAT (7X,'Frequency  (days)        [FLXF] =',11F7.2)
11560 FORMAT (5X,'Vector plot      [VPLC] = ',A)
11570 FORMAT (7X,'Number of time intervals [NVPL] =',I7)
11580 FORMAT (7X,'Date  (Julian day)       [VPLD] =',11F7.2)
11590 FORMAT (7X,'Frequency  (days)        [VPLF] =',11F7.2)
11600 FORMAT (5X,'Profile plot     [PRFC] = ',A)
11610 FORMAT (7X,'Number of time intervals [NPRF] =',I7)
11620 FORMAT (7X,'Number of stations      [NIPRF] =',I7)
11630 FORMAT (7X,'Segment location         [IPRF] =',11I7)
11640 FORMAT (7X,'Date  (Julian day)       [PRFD] =',11F7.2)
11650 FORMAT (7X,'Frequency  (days)        [PRFF] =',11F7.2)
11660 FORMAT (5X,'Spreadsheet plot [SPRC] = ',A)
11670 FORMAT (7X,'Number of time intervals [NSPR] =',I7)
11680 FORMAT (7X,'Number of stations      [NISPR] =',I7)
11690 FORMAT (7X,'Segment location         [ISPR] =',11I7)
11700 FORMAT (7X,'Date  (Julian day)       [SPRD] =',11F7.2)
11710 FORMAT (7X,'Frequency  (days)        [SPRF] =',11F7.2)
11720 FORMAT (5X,'Contour plot     [CPLC] = ',A)
11730 FORMAT (7X,'Number of time intervals [NCPL] =',I7)
11740 FORMAT (7X,'Date (Julian day)        [CPLD] =',11F7.2)
11750 FORMAT (7X,'Frequency  (days)        [CPLF] =',11F7.2)
11760 FORMAT (5X,'Restart out      [RSOC] = ',A/                                                                                   &
              5X,'Restart in       [RSIC] = ',A)
11770 FORMAT (7X,'Number of time intervals [NRSO] =',I7)
11780 FORMAT (7X,'Date (Julian day)        [RSOD] =',11F7.2)
11790 FORMAT (7X,'Frequency  (days)        [RSOF] =',11F7.2)
11800 FORMAT (/'Inflow/Outflow'/                                                                                                   &
              3X,'Selective Withdrawal'/                                                                                           &
              5X,'Branch ',3X,'# of structures [NSTR]':/                                                                           &
             (7X,I2,14X,I3))
11805 FORMAT (5X,'Branch ',I0)
11810 FORMAT (:7X,'Structure',3X,'Type',3X,'Width (m)',3X,'Elevation (m)',3X,'Bottom Layer'/                                       &
             (10X,I2,4X,A8,1X,F8.1,5X,F8.1,9X,I5))
11820 FORMAT (5X,'gate coefficients'/                                                                                              &
              7X,'[A1GT]',4X,'[B1GT]',4X,'[G1GT]',4X,'[A2GT]',4X,'[B2GT]',4X,'[G2GT]'/                                             &
              :(2X,6(2X,F8.2)))
11821 FORMAT (3X,'Number of spillways/external weirs [NSP] = ',I0:/                                                                &
              5X,'centerline elevation             [ESP] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11822 FORMAT (5X,'upstream segment number         [IUSP] = ',7I8,                                                                  &
             (:/T26,40I3))
11823 FORMAT (5X,'downstream segment number       [IDSP] = ',7I8,                                                                  &
             (:/T26,40I3))
11824 FORMAT (5X,'upstream inflow placement      [PUSPC] = ',7A8,                                                                  &
             (:T33,7A8))
11825 FORMAT (5X,'downstream inflow placement    [PDSPC] = ',7A8,                                                                  &
             (:T33,7A8))
11826 FORMAT (5X,'upstream top elevation         [EUTSP] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11827 FORMAT (5X,'downstream top elevation       [EDTSP] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11828 FORMAT (5X,'upstream bottom elevation      [EUBSP] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11829 FORMAT (5X,'downstream bottom elevation    [EDBSP] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11830 FORMAT (5X,'spillway/external weir coefficients'/                                                                            &
              7X,'[A1SP]',4X,'[B1SP]',4X,'[A2SP]',4X,'[B2SP]'/                                                                     &
              :(2X,4(2X,F8.2)))
11831 FORMAT (3X,'Number of pipes [NPI] = ',I0:/                                                                                   &
              5X,'upstream centerline elevation   [EUPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11841 FORMAT (5X,'downstream centerline elevation [EDPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11842 FORMAT (5X,'diameter                         [WPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11843 FORMAT (5X,'length                         [DLXPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11844 FORMAT (5X,'Mannings n                       [FPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11845 FORMAT (5X,'minor friction losses         [FMINPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11832 FORMAT (5X,'upstream segment number         [IUPI] = ',7I8,                                                                  &
             (:/T26,40I3))
11833 FORMAT (5X,'downstream segment number       [IDPI] = ',7I8,                                                                  &
             (:/T26,40I3))
11834 FORMAT (5X,'upstream inflow placement      [PUPIC] = ',7A8,                                                                  &
             (:T33,7A8))
11835 FORMAT (5X,'downstream inflow placement    [PDPIC] = ',7A8,                                                                  &
             (:T33,7A8))
11836 FORMAT (5X,'upstream top elevation         [EUTPI] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11837 FORMAT (5X,'downstream top elevation       [EDTPI] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11838 FORMAT (5X,'upstream bottom elevation      [EUBPI] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11839 FORMAT (5X,'downstream bottom elevation    [EDBPI] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11840 FORMAT (3X,'Number of withdrawals [NWD] = ',I0:/                                                                             &
              5X,'segment number [IWD] = ',10I3,                                                                                   &
             (:/T26,10I3))
11851 FORMAT (3X,'Number of pumps [NPU] = ',I0:/                                                                                   &
              5X,'centerline elevation             [EPU] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11852 FORMAT (5X,'upstream segment number         [IUPU] = ',7I8,                                                                  &
             (:/T26,40I3))
11853 FORMAT (5X,'downstream segment number       [IDPU] = ',7I8,                                                                  &
             (:/T26,40I3))
11854 FORMAT (5X,'starting time                 [STRTPU] = ',7F8.2,                                                                &
             (:T33,7A8))
11855 FORMAT (5X,'ending time                    [ENDPU] = ',7F8.2,                                                                &
             (:T33,7A8))
11856 FORMAT (5X,'starting elevation             [EONPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11857 FORMAT (5X,'ending elevaton               [EOFFPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11858 FORMAT (5X,'pump flowrate                    [QPU] = ',F8.2,                                                                 &
             (:T33,7F8.2))
11859 FORMAT (5X,'inflow placement                [PPUC] = ',7A8,                                                                  &
             (:T33,7F8.2))
11861 FORMAT (5X,'top elevation                   [ETPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11862 FORMAT (5X,'bottom elevation                [EBPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11863 FORMAT (5X,'top layer                       [KTPU] = ',7I8,                                                                  &
             (:T33,7F8.2))
11864 FORMAT (5X,'bottom layer                    [KBPU] = ',7I8,                                                                  &
             (:T33,7F8.2))
11860 FORMAT (3X,'Number of tributaries [NTR] = ',I0:/                                                                             &
              5X,'segment number     [ITR] = ',7I8,                                                                                &
             (:/T33,7I8))
11870 FORMAT (:5X,'Inflow placement  [PTRC] = ',7A8,                                                                               &
             (/:T33,7A8))
11880 FORMAT (:5X,'Top elevation     [ETTR] = ',7F8.2,                                                                             &
             (:/T33,7F8.2))
11890 FORMAT (:5X,'Bottom elevation  [EBTR] = ',7F8.2,                                                                             &
             (:/T33,7F8.2))
11891 FORMAT (3X,'Number of gates [NGT] =',I3:/                                                                                    &
              5X,'centerline elevation             [EGT] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11892 FORMAT (:5X,'upstream segment number         [IUGT] = ',7I8,                                                                 &
             (:/T26,40I3))
11893 FORMAT (:5X,'downstream segment number       [IDGT] = ',7I8,                                                                 &
             (:/T26,40I3))
11894 FORMAT (:5X,'upstream inflow placement      [PUGTC] = ',7A8,                                                                 &
             (:T33,7A8))
11895 FORMAT (:5X,'downstream inflow placement    [PDGTC] = ',7A8,                                                                 &
             (:T33,7A8))
11896 FORMAT (:5X,'upstream top elevation         [EUTGT] = ',7F8.2,                                                               &
             (:T33,7F8.2))
11897 FORMAT (:5X,'downstream top elevation       [EDTGT] = ',7F8.2,                                                               &
             (:T33,7F8.2))
11898 FORMAT (:5X,'upstream bottom elevation      [EUBGT] = ',7F8.2,                                                               &
             (:T33,7F8.2))
11899 FORMAT (:5X,'downstream bottom elevation    [EDBGT] = ',7F8.2,                                                               &
             (:T33,7F8.2))
11900 FORMAT (3X,'Distributed tributaries [DTRC]')
11910 FORMAT (5X,'Branch ',I0,' = ',A3)

11930 FORMAT (/1X,'Input Filenames'/3X,'Control    = ',A72/3X,'Restart    = ',A72/3X,'Withdrawal = ',A72)
11940 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Bathymetry           = ',A72/                                                                                    &
              5X,'Meteorology          = ',A72/                                                                                    &
              5X,'Vertical profile     = ',A72/                                                                                    &
              5X,'Longitudinal profile = ',A72)
11950 FORMAT (3X,'Branch ',I0/                                                                                                     &
              5X,'Inflow                               = ',A72/                                                                    &
              5X,'Inflow temperature                   = ',A72/                                                                    &
              5X,'Inflow concentrations                = ',A72/                                                                    &
              5X,'Outflow                              = ',A72/                                                                    &
              5X,'Distributed tributary inflows        = ',A72/                                                                    &
              5X,'Distributed tributary temperatures   = ',A72/                                                                    &
              5X,'Distributed tributary concentrations = ',A72/                                                                    &
              5X,'Precipitation                        = ',A72/                                                                    &
              5X,'Precipitation temperatures           = ',A72/                                                                    &
              5X,'Precipitation concentrations         = ',A72/                                                                    &
              5X,'Upstream head                        = ',A72/                                                                    &
              5X,'Upstream head temperatures           = ',A72/                                                                    &
              5X,'Upstream head concentrations         = ',A72/                                                                    &
              5X,'Downstream head                      = ',A72/                                                                    &
              5X,'Downstream head temperatures         = ',A72/                                                                    &
              5X,'Downstream head concentrations       = ',A72)
11970 FORMAT (:(3X,'Tributary ',I0/                                                                                                &
              5X,'Inflow               = ',A72/                                                                                    &
              5X,'Inflow temperature   = ',A72/                                                                                    &
              5X,'Inflow concentration = ',A72)/)
11980 FORMAT (1X,'Output Filenames'/                                                                                               &
              3X,'Error       = ',A72/                                                                                             &
              3X,'Warning     = ',A72/                                                                                             &
              3X,'Time series = ',A72/                                                                                             &
              3X,'Withdrawal  = ',A72)
11990 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Snapshot     = ',A72/                                                                                            &
              5X,'Fluxes       = ',A72/                                                                                            &
              5X,'Profile      = ',A72/                                                                                            &
              5X,'Vector plot  = ',A72/                                                                                            &
              5X,'Contour plot = ',A72)
12000 FORMAT (/1X,'Constituents [CCC] = ',A3/                                                                                      &
              3X,'Algal limiting nutrient  [LIMC] = ',A3/                                                                          &
              3X,'Kinetics update frequency [CUF] = ',I3)
12020 FORMAT (5X,'State Variables'/                                                                                                &
              7X,'Constituent',T29,'Computation',T42,'Initial Conc',T57,'Fluxes',T67,'Printout'/                                   &
              9X,'[CNAME]',T32,'[CAC]',T42,'[C2IWB,g/m^3]',T57,'[CFWBC]',T67,'[CPRWBC]')
12030 FORMAT (/3X,'Waterbody ',I0)
12040 FORMAT (7X,A19,T33,A3,T42,F10.3,T59,A3,T69,A3)
12045 FORMAT (/3X,'Branch ',I0/                                                                                                    &
              5X,'State Variables'/                                                                                                &
              7X,'Constituent',T32,'Inflow',T42,'Distributed trib',T62,'Precipitation'/                                            &
              9X,'[CNAME]',T31,'[CINBRC]',T46,'[CDTBRC]',T64,'[CPRBRC]'/                                                           &
             (7X,A19,T33,A3,T48,A3,T67,A3))
12046 FORMAT (/3X,'Tributary ',I0/                                                                                                 &
              5X,'State Variables'/                                                                                                &
              7X,'Constituent',T32,'Inflow'/                                                                                       &
              9X,'[CNAME]',T31,'[CINTRC]'/                                                                                         &
             (7X,A19,T33,A3))
12050 FORMAT (/5X,'Derived Variables')
12060 FORMAT (19X,'Constituent',T46,'Computation'/                                                                                 &
              20X,'[CDNAME]',T48,'[CDWBC]')
12070 FORMAT (7X,A,T50,1000A3)
12080 FORMAT (5X,A/                                                                                                                &
              7X,'Constituent',13X,'Rate/Coefficient')
12090 FORMAT (7X,'Suspended solids',T27,'Settling',T50,'[SSS] =',F6.3,' m/day')
12100 FORMAT (7X,A21,T27,'Temperature mult',T48,'[CGQ10] =',F6.3/                                                                  &
              T27,'0-Order Decay        [CG0DK] =',F6.3,'g/m^3/day'/                                                               &
              T27,'1-Order Decay        [CG1DK] =',F6.3,'/day'/                                                                    &
              T27,'Settling               [CGS] =',F6.3,' m/day')
12110 FORMAT (7X,'Labile DOM',T27,'Decay',T47,'[LDOMDK] =',F6.3,' /day'/                                                           &
              T27,'to refractory        [LRDDK] =',F6.3,' /day')
12120 FORMAT (7X,'Refractory DOM',T27,'Decay',T47,'[RDOMDK] =',F6.3,' /day')
12130 FORMAT (7X,'Labile POM',T27,'Decay',T47,'[LPOMDK] =',F6.3,' /day'/                                                           &
              T27,'to refractory        [LRPDK] =',F6.3,' /day'/                                                                   &
              T27,'Settling',T49,'[POMS] =',F6.3,' m/day')
12150 FORMAT (7X,'Refractory POM',T27,'Decay',T47,'[RPOMDK] =',F6.3,' /day')
12160 FORMAT (7X,'Algal group',I3,T27,'Growth',T51,'[AG] =',F6.3,' /day'/                                                          &
              T27,'Mortality               [AM] =',F6.3,' /day'/                                                                   &
              T27,'Excretion               [AE] =',F6.3,' /day'/                                                                   &
              T27,'Respiration             [AR] =',F6.3,' /day'/                                                                   &
              T27,'Settling                [AS] =',F6.3,' m/day'/                                                                  &
              T27,'Org-P                 [ALGP] =',F6.3/                                                                           &
              T27,'Org-N                 [ALGN] =',F6.3/                                                                           &
              T27,'Org-C                 [ALGC] =',F6.3/                                                                           &
              T27,'Org-Si               [ALGSI] =',F6.3/                                                                           &
              T27,'Algae/chl a ratio    [ACHLA] =',F6.2,' mg/ug'/                                                                  &
              T27,'Fraction algae to POM [APOM] =',F6.2/)      ! MLM 7/11/06
12180 FORMAT (7X,'Phosphorous',T27,'Release',T46,'[PO4R] =',F6.3,' g/m^2/day')
12190 FORMAT (7X,'Ammonium',T27,'Decay',T45,'[NH4DK] =',F6.3,' /day'/                                                              &
              T27,'Release',T46,'[NH4R] =',F6.3,' g/m^2/day')
12200 FORMAT (7X,'Nitrate-Nitrite',T27,'Decay',T45,'[NO3DK] =',F6.3,' /day'/   &
              T27,'Diffusion to sed',T45,'[NO3S] =',F6.3,' m/day'/                           &
              T27,'Frac to OM',T45,'[FNO3SED] =',F6.3,' [-]')
12210 FORMAT (7X,'Silica',T27,'Decay',T45,'[PSIDK] =',F6.3,' /day'/                                                                &
              T27,'Release',T46,'[DSIR] =',F6.3,' g/m^2/day'/                                                                      &
              T27,'Settling',T46,'[PSIS] =',F6.3,' m/day')
12220 FORMAT (7X,'Sediment',T27,'Decay',T45,'[SEDDK] =',F6.3,' /day')
12230 FORMAT (7X,'Iron',T27,'Settling',T47,'[FES] =',F6.3,' m/day'/                                                                &
              T27,'Release',T47,'[FER] =',F6.3,' g/m^2/day')
12240 FORMAT (7X,'Oxygen',T27,'Sediment demand     [SOD] = ',13F5.1:/                                                              &
             (T55,13F5.1))
12250 FORMAT (T27,'SOD fraction',T46,'[FSOD] = ',F5.1/                                                                             &
              T27,'Sediment fraction',T46,'[FSED] = ',F5.1)
12255 FORMAT (T27,'Sediment focusing ',T46,'[SEDS] = ',F7.3,' m/day'/                                                                             &
              T27,'Sediment burial ',T46,'[SEDBR] = ',F7.3,' /day')
12270 FORMAT (5X,A/                                                                                                                &
              9X,'Constituent',T24,'Rate',T37,'Lower',T50,'Max Lower')
12280 FORMAT (7X,'Algal group ',I0,T23,'Growth',T34,'[AT1] =',F5.1,T49,'[AT2] =',F5.1)
12290 FORMAT (7X,'Ammonium',T23,'Decay',T32,'[NH4T1] =',F5.1,T47,'[NH4T2] =',F5.1)
12300 FORMAT (7X,'Nitrate',T23,'Decay',T32,'[NO3T1] =',F5.1,T47,'[NO3T2] =',F5.1)
12310 FORMAT (7X,'Organic',T23,'Decay',T33,'[OMT1] =',F5.1,T48,'[OMT2] =',F5.1)
12320 FORMAT (7X,'Sediment',T23,'Decay',T32,'[SEDT1] =',F5.1,T47,'[SEDT2] =',F5.1)
12330 FORMAT (5X,A/                                                                                                                &
              9X,'Constituent',T24,'Rate',T37,'Upper',T50,'Max Upper')
12340 FORMAT (7X,'Algal group ',I0,T23,'Growth',T34,'[AT3] =',F5.1,T49,'[AT4] =',F5.1)
12350 FORMAT (5X,A/                                                                                                                &
              7X,'Oxygen'/                                                                                                         &
              9X,'Ammonium       [O2NH4] =',F5.2/                                                                                  &
              9X,'Organic matter  [O2OM] =',F5.2/)
12355 FORMAT (9X,'Respiration     [O2AR] =',F5.2/                                                                                  &
              9X,'Algal growth    [O2AG] =',F5.2)
12360 FORMAT (7X,'Organic Matter'/                                                                                                 &
              9X,'Carbon      [BIOC] =',F6.3/                                                                                      &
              9X,'Phosphorous [BIOP] =',F6.3/                                                                                      &
              9X,'Nitrogen    [BION] =',F6.3/                                                                                      &
              9X,'Silica     [BIOSI] =',F6.3)
12370 FORMAT (:5X,'CBOD ',I0/                                                                                                      &
              7X,'Decay rate                   [KBOD] =',F6.3,' /day'/                                                             &
              7X,'Temperature adjustment       [TBOD] =',F6.3/                                                                     &
              7X,'Ultimate CBOD to CBOD5 ratio [RBOD] =',F6.3)
12380 FORMAT (5X,A15/                                                                                                              &
             (7X,'Algal group ',I0/                                                                                                &
              9X,'Phosphorous [AHSP] =',F6.3,' g/m^3'/                                                                             &
              9X,'Nitrogen    [AHSN] =',F6.3,' g/m^3'/                                                                             &
              9X,'Silica     [AHSSI] =',F6.3,' g/m^3'))
12390 FORMAT (5X,A/                                                                                                                &
              7X,'Attenuation'/                                                                                                    &
              9X,'Surface layer      [BETA] =',F5.2/                                                                               &
              9X,'Water             [EXH2O] =',F5.2,' /m'/                                                                         &
              9X,'Inorganic solids   [EXSS] =',F5.2,' /m'/                                                                         &
              9X,'Organic solids     [EXOM] =',F5.2,' /m')
12400 FORMAT (9X,'Algal group ',I0,'       [EXA] =',F5.2,' /m')
12410 FORMAT (7X,'Saturation Intensity'/                                                                                           &
             (9X,'Algal group ',I0,' [ASAT] =',F6.1,' W/m^2'))
12420 FORMAT (5X,A/                                                                                                                &
              7X,'Oxygen          [DMO2] =',1PE10.3,' m^2/g'/                                                                      &
              7X,'Carbon dioxide [DMCO2] =',1PE10.3,' m^2/g')
12430 FORMAT (5X,A/                                                                                                                &
              7X,'Phosphorous [PARTP] =',F6.3,' m^3/g'/                                                                            &
              7X,'Silica     [PARTSI] =',F6.3,' m^3/g')
12440 FORMAT (5X,A/                                                                                                                &
              7X,'Half-saturation coeff anoxia[KDO] =',F6.2,' g/m^3'/                                                              &
              7X,'CO2 sediment release       [CO2R] =',F6.2,' g/m^2/day')
12460 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              T7,'total inflows',T25,'total outflows'/                                                                             &
              T5,'average  maximum',T24,'average  maximum'/                                                                        &
              T4,F8.2,T13,F8.2,T23,F8.2,T32,F8.2/)
12470 FORMAT (5X,'Branch ',I0/                                                                                                     &
              7X,'Inflows'/                                                                                                        &
              14X,'total'/                                                                                                         &
              9X,'average  maximum'/                                                                                               &
              9X,F8.2,1X,F8.2/                                                                                                     &
              13X,'upstream',9X,'tributaries',5X,'distributed tributaries',5X,'precipitation'/                                     &
              9X,'average  maximum',3X,'average  maximum',6X,'average  maximum',6X,'average  maximum'/                             &
              8X,F8.2,1X,F8.2,2X,F8.2,1X,F8.2,5X,F8.2,1X,F8.2,5X,F8.2,1X,F8.2/                                                     &
              7X,'Outflows'/                                                                                                       &
              14X,'outlets',9X,'withdrawals'/                                                                                      &
              9X,'average  maximum',3X,'average  maximum'/                                                                         &
              8X,F8.2,1X,F8.2,2X,F8.2,1X,F8.2)
12480 FORMAT (/A8/                                                                                                                 &
              3X,'Overall Grid'/                                                                                                   &
              5X,'Total'/7X,'segments [IMX] = ',I0/                                                                                &
              7X,'layers   [KMX] = ',I0/                                                                                           &
              7X,'branches [NBR] = ',I0)
12490 FORMAT (/3X,'Waterbody ',I0/                                                                                                 &
              5X,'Segments                 = ',I0,'-',I0/                                                                          &
              5X,'Branches                 = ',I0,'-',I0/                                                                          &
              5X,'Bottom elevation [ELBOT] =',F9.2,' m'/                                                                           &
              5X,'Surface layer       [KT] =',I9/                                                                                  &
              5X,'Vertical spacing     [H]')
12500 FORMAT (7X,'Layer',6X,14I5)
12510 FORMAT (7X,'Height (m) ',14F5.1)
12520 FORMAT (5X,'Branch ',I0/                                                                                                     &
              7X,'Upstream segment',T31,'[US] = ',I0,T44,'Downstream ','segment       [DS] = ',I0/                                 &
              7X,'Upstream head segment [UHS] = ',I0,T44,'Downstream ','head segment [DHS] = ',I0)
12530 FORMAT (8X,'Waterbody ',I0,' Bathymetry [B], m')
12600 FORMAT (/1X,'Initial Branch Volume [VOLB] = ',F0.1,' m^3'/)
12620 FORMAT (T34,'Branch ',I0,' Volume-Area-Elevation Table'/T34,'Note: Elevation is at top of layer'/)
12630 FORMAT (1X,'Layer',T12,'Elevation',T27,'Area',T39,'Volume',T48,'Active Cells',T62,'Average depth',T77,'Average width',T92,   &
                 'Residence time'/                                                                                                 &
              T15,'(m)',T24,'(1.0E6 m^2)',T37,'(1.0E6 m^3)',T66,'(m)',T82,'(m)',T96,'(days)'/)
12640 FORMAT (1X,'Layer',T12,'Elevation',T27,'Area',T39,'Volume',T48,'Active Cells',T62,'Average depth',T77,'Average width',/      &
              T15,'(m)',T24,'(1.0E6 m^2)',T37,'(1.0E6 m^3)',T66,'(m)',T82,'(m)'/)
12650 FORMAT (1X,I3,' KT',T12,F7.2,T24,F8.3,T36,F12.3,T50,I5,T62,F8.1,T78,F9.2,T93,F12.2)
12655 FORMAT (1X,I3,      T12,F7.2,T24,F8.3,T36,F12.3,T50,I5,T62,F8.1,T78,F9.2,T93,F12.2)
12670 FORMAT (/1X,'Waterbody ',I0,' Initial Volume [VOLG] = ',F0.1,' m^3'/)
12690 FORMAT (T34,'Waterbody ',I0,' Volume-Area-Elevation Table'/T34,'Note: Elevation is at top of layer'/)
99999 CONTINUE
END PROGRAM W2_PRE

!***********************************************************************************************************************************
!*                                         S U B R O U T I N E   P R I N T  I N I T I A L  G R I D                                **
!***********************************************************************************************************************************

SUBROUTINE PRINT_INITIAL_GRID (KT)
  USE CONTROL; USE GRIDC

! Variable initialization

  CONV = '          '

! Water surface variables

  WRITE (INI,'(//5X,A18/)')   'Water Surface [Z], m'
  DO I=1,IEPR(JW)
    WRITE (CONV(1,I),'(F10.4)') Z(IPR(I,JW))
  END DO
  WRITE (INI,'(2X,1000I10)')  (IPR(I,JW),I=1,IEPR(JW))
  WRITE (INI,'(3X,1000A10/)') (CONV(1,I),I=1,IEPR(JW))
  IF (ICE_CALC(JW)) WRITE (INI,'(/2X,A24/3X,1000F10.4/)') 'Ice Thickness (m)', (ICETH(IPR(I,JW)),I=1,IEPR(JW))

! Temperatures

  IF (HPRWBC(4,JW) == '      ON') THEN
    DO I=1,IEPR(JW)
      DO K=KT,KB(IPR(I,JW))
        WRITE (CONV(K,I),'(F10.3)') T2(K,IPR(I,JW))
      END DO
    END DO
    WRITE (INI,'(/2X,A35/)')    HNAME(4)
    WRITE (INI,'(2X,1000I10)') (IPR(I,JW),I=1,IEPR(JW))
    DO K=KT,KEPR(JW)
      WRITE (INI,'(1X,I3,1000A10)') K,(CONV(K,I),I=1,IEPR(JW))
    END DO
  END IF

! Constituent concentrations

  DO J=1,NAC
    IF (CPRWBC(CN(J),JW) == '      ON') THEN
      DO I=1,IEPR(JW)
        DO K=KT,KB(IPR(I,JW))
          WRITE (CONV(K,I),'(F10.3)') C2(K,IPR(I,JW),CN(J))*CMULT(CN(J))
        END DO
      END DO
      WRITE (INI,'(/5X,A37/)')    CNAME(CN(J))
      WRITE (INI,'(2X,1000I10)') (IPR(I,JW),I=1,IEPR(JW))
      DO K=KT,KEPR(JW)
        WRITE (INI,'(1X,I2,1000A10)') K,(CONV(K,I),I=1,IEPR(JW))
      END DO
    END IF
  END DO

! Derived constituent concentrations

  CALL DERIVED_CONSTITUENTS
  DO J=1,NADC(JW)
    IF (CDWBC(CDN(J,JW),JW) == '      ON') THEN
      DO I=1,IEPR(JW)
        DO K=KT,KB(IPR(I,JW))
          WRITE (CONV(K,I),'(F10.3)') CD(K,IPR(I,JW),CDN(J,JW))*cDMULT(CDN(J,JW))     ! RA 7/11/08
        END DO
      END DO
      WRITE (INI,'(/5X,A43/)')    CDNAME(CDN(J,JW))
      WRITE (INI,'(2X,1000I10)') (IPR(I,JW),I=1,IEPR(JW))
      DO K=KT,KEPR(JW)
        WRITE (INI,'(1X,I2,1000A10)') K,(CONV(K,I),I=1,IEPR(JW))
      END DO
    END IF
  END DO
END SUBROUTINE PRINT_INITIAL_GRID

!***********************************************************************************************************************************
!*                                  S U B R O U T I N E    D E R I V E D  C O N S T I T U E N T S                                 **
!***********************************************************************************************************************************

SUBROUTINE DERIVED_CONSTITUENTS
  USE CONTROL; USE KINETICS; USE GRIDC
  REAL, ALLOCATABLE, DIMENSION(:,:) :: DOM, POM
  ALLOCATE (DOM(KMX,IMX), POM(KMX,IMX))

  IU    =  CUS(BS(JW))
  ID    =  DS(BE(JW))
  TISS  =  0.0
  ATOT  =  0.0
  TISS  =  0.0
  TOTSS =  0.0
  CHLA  =  0.0
  PALT  = (1.0-((EL(KTWB(JW),(IU+ID)/2)-Z((IU+ID)/2)*COS(ALPHA(BE(JW))))/1000.0)/44.3)**5.25
  DO I=IU,ID
    DO K=KTWB(JW),KB(I)
      DO JS=1,NASS
        TISS(K,I) = TISS(K,I)+SS(K,I,SSAC(JS))
      END DO
    END DO
    DO K=KTWB(JW),KB(I)
      DO JA=1,NAAL
        ATOT(K,I) = ATOT(K,I)+ALG(K,I,ALAC(JA))
      END DO
      DOM(K,I) = LDOM(K,I)+RDOM(K,I)
      POM(K,I) = LPOM(K,I)+RPOM(K,I)
      DOC(K,I) = DOM(K,I)*ORGC(JW)
      POC(K,I) = POM(K,I)*ORGC(JW)
      DO JA=1,NAAL
        POC(K,I) = POC(K,I)+ALG(K,I,ALAC(JA))*ALGC(ALAC(JA))
      END DO
      TOC(K,I)   = DOC(K,I)+POC(K,I)
      DOP(K,I)   = DOM(K,I)*ORGP(JW)
      DON(K,I)   = DOM(K,I)*ORGN(JW)
      POP(K,I)   = POM(K,I)*ORGP(JW)
      PON(K,I)   = POM(K,I)*ORGN(JW)
      TOP(K,I)   = DOP(K,I)+POP(K,I)
      TON(K,I)   = DON(K,I)+PON(K,I)
      TKN(K,I)   = TON(K,I)+NH4(K,I)
      CBODU(K,I) = O2OM(JW)*(DOM(K,I)+POM(K,I)+ATOT(K,I))
      TPSS       = 0.0
      DO JS=1,NASS
        TPSS = TPSS+SS(K,I,SSAC(JS))*PARTP(JW)
      END DO
      TP(K,I)   =  TOP(K,I)+PO4(K,I)+TPSS
      TN(K,I)   =  TON(K,I)+NH4(K,I)+NO3(K,I)
      SATO      =  EXP(7.7117-1.31403*(LOG(T2(K,I)+45.93)))*PALT
      O2DG(K,I) = (O2(K,I)/SATO)*100.0
      DO JA=1,NAAL
        TOTSS(K,I) = TOTSS(K,I)+ALG(K,I,JA)
        IF (ACHLA(ALAC(JA)) /= 0.0) CHLA(K,I)  = CHLA(K,I)+ALG(K,I,ALAC(JA))/ACHLA(ALAC(JA))
      END DO
      TOTSS(K,I) = TOTSS(K,I)+TISS(K,I)+POM(K,I)
    END DO
  END DO
END SUBROUTINE DERIVED_CONSTITUENTS

!***********************************************************************************************************************************
!*                                            S U B R O U T I N E    W A R N I N G S                                              **
!***********************************************************************************************************************************

SUBROUTINE WARNINGS
  USE CONTROL
  NWRN       =  NWRN+1
  DELETE_WRN = .FALSE.
END SUBROUTINE WARNINGS

!***********************************************************************************************************************************
!*                                              S U B R O U T I N E    E R R O R S                                                **
!***********************************************************************************************************************************

SUBROUTINE ERRORS
  USE CONTROL
  NERR       =  NERR+1
  DELETE_ERR = .FALSE.
END SUBROUTINE ERRORS

!***********************************************************************************************************************************
!*                                           S U B R O U T I N E    I N P U T  E R R O R                                          **
!***********************************************************************************************************************************

SUBROUTINE INPUT_ERROR (UNIT,ERR,AID,MESSAGE)
  INTEGER       :: UNIT,    ERR,      WIN=5
  CHARACTER(8)  :: AID
  CHARACTER(70) :: MESSAGE, ERROR(5)
  ERROR = ' '
  DO J=1,5
    BACKSPACE (UNIT)
  END DO
  DO J=1,5
    READ (UNIT,'(A70)',ERR=1) ERROR(J)
  END DO
1 CONTINUE
  WRITE (ERR,'(A/(A))') MESSAGE,ERROR
  WRITE (WIN,'(A/(A))') MESSAGE,ERROR
  WRITE (ERR,'(/A/)')  'Last value of card header was '//AID
  WRITE (WIN,'(/A/)')  'Last value of card header was '//AID
  CLOSE (ERR)
  STOP
END SUBROUTINE INPUT_ERROR

