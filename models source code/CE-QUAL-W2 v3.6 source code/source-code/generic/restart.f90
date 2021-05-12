
!***********************************************************************************************************************************
!*                                       S U B R O U T I N E    R E S T A R T   O U T P U T                                       **
!***********************************************************************************************************************************

SUBROUTINE RESTART_OUTPUT (RSOFN)
!  USE GLOBAL; USE SCREENC; USE RSTART; USE GDAYC; USE GEOMC; USE KINETIC, ONLY:EPM, EPD; USE TVDC, ONLY:QSUM; USE KINETIC, ONLY:SED
  USE GLOBAL; USE SCREENC; USE RSTART; USE GDAYC; USE GEOMC; USE KINETIC, ONLY:EPM,EPD,SEDC,SEDN,SEDP,PH, sdkv; USE TVDC, ONLY:QSUM !MLM 6/10/07;10/06
  USE KINETIC, ONLY:SED ; USE ZOOPLANKTONC, ONLY: ZOO; USE EDDY, ONLY: TKE
  CHARACTER(*) :: RSOFN
  OPEN  (RSO,FILE=RSOFN,FORM='UNFORMATTED',STATUS='UNKNOWN')
  WRITE (RSO) NIT,    NV,     KMIN,   IMIN,   NSPRF,  CMBRT,  ZMIN,   IZMIN,  START,  CURRENT
  WRITE (RSO) DLTDP,  SNPDP,  TSRDP,  VPLDP,  PRFDP,  CPLDP,  SPRDP,  RSODP,  SCRDP,  FLXDP,  WDODP
  WRITE (RSO) JDAY,   YEAR,   ELTM,   ELTMF,  DLT,    DLTAV,  DLTS,   MINDLT, JDMIN,  CURMAX
  WRITE (RSO) NXTMSN, NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS, NXTMSC, NXTMSP, NXTMFL, NXTMWD
  WRITE (RSO) VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTRB, VOLDT,  VOLWD,  VOLEV,  VOLSBR, VOLTR, VOLSR
  WRITE (RSO) TSSEV,  TSSPR,  TSSTR,  TSSDT,  TSSWD,  TSSIN,  TSSOUT, TSSS,   TSSB,   TSSICE
  WRITE (RSO) TSSUH,  TSSDH,  TSSUH2, TSSDH2, CSSUH2, CSSDH2, VOLUH2, VOLDH2, QUH1
  WRITE (RSO) ESBR,   ETBR,   EBRI
  WRITE (RSO) Z,      SZ,     ELWS,   SAVH2,  SAVHR,  H2
  WRITE (RSO) KTWB,   KTI,    SKTI,   SBKT
  WRITE (RSO) ICE,    ICETH,  CUF,    QSUM
  WRITE (RSO) U,      W,      SU,     SW,     AZ,     SAZ,    DLTLIM
  WRITE (RSO) T1,     T2,     C1,     C2,     C1S,    EPD,    SED,    KFS,    CSSK
  WRITE (RSO) SEDC, SEDN, SEDP, ZOO, CD  ! mlm 10/06
  write(rso) sdkv                        ! mlm 6/10/07
  write(rso) tke                         ! SW 10/4/07
  CLOSE (RSO)
END SUBROUTINE RESTART_OUTPUT
