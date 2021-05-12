subroutine outputinit

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc
  EXTERNAL RESTART_OUTPUT

REAL    DIST
INTEGER JN
!***********************************************************************************************************************************
!*                                                           Task 1.5: Outputs                                                    **
!***********************************************************************************************************************************

! Open output files

  IF (RESTART_IN) THEN
    DO JW=1,NWB
      IF (SNAPSHOT(JW))    OPEN (SNP(JW),FILE=SNPFN(JW),POSITION='APPEND')
      IF (VECTOR(JW))      OPEN (VPL(JW),FILE=VPLFN(JW),POSITION='APPEND')
      IF (SPREADSHEET(JW)) OPEN (SPR(JW),FILE=SPRFN(JW),POSITION='APPEND')
      IF (CONTOUR(JW))     OPEN (CPL(JW),FILE=CPLFN(JW),POSITION='APPEND')
      IF (PROFILE(JW))     OPEN (PRF(JW),FILE=PRFFN(JW),POSITION='APPEND')
      IF (FLUX(JW))THEN
              OPEN (FLX(JW),FILE=FLXFN(JW),POSITION='APPEND')
              WRITE (SEGNUM,'(I0)') JW
              SEGNUM = ADJUSTL(SEGNUM)
              L      = LEN_TRIM(SEGNUM)
              OPEN (FLX2(JW), FILE='kflux_jw'//SEGNUM(1:L)//'.opt',POSITION='APPEND')
      ENDIF
      IF (SNAPSHOT(JW)) THEN
        REWIND (SNP(JW))
        DO WHILE (.TRUE.)
          READ (SNP(JW),'(A72)',END=100) LINE
          IF (LINE(26:28) == 'NIT') THEN
            BACKSPACE SNP(JW)
            READ (SNP(JW),'(31X,I10)',END=100) NIT1
            IF (NIT1 > NIT) THEN
              DO J=1,24
                BACKSPACE (SNP(JW))
              END DO
              EXIT
            END IF
          END IF
        END DO
      END IF
100   CONTINUE
      IF (SPREADSHEET(JW)) THEN
        REWIND (SPR(JW))
        READ (SPR(JW),*)
        DO WHILE (JDAY1 < JDAY)
          READ (SPR(JW),'(A,F10.0)',END=101) LINE(1:38),JDAY1
        END DO
        BACKSPACE (SPR(JW))
101     CONTINUE
        JDAY1 = 0.0
      END IF
      IF (PROFILE(JW)) THEN
        REWIND (PRF(JW))
        READ   (PRF(JW),'(A)')        (LINE,J=1,11)
        READ   (PRF(JW),'(8I8)')       I
        READ   (PRF(JW),'(10I8)')     (I,J=1,NIPRF(JW))
        READ   (PRF(JW),'(20(1X,A))')  LINE (1:8), (LINE (1:3), JC=1,NCT),(LINE (1:3), JD=1,NDC)
        READ   (PRF(JW),'(2A)')        LINE (1:26),(LINE (1:26),JC=1,NCT),(LINE (1:43),JD=1,NDC)
        DO WHILE (JDAY1 < JDAY)
          READ (PRF(JW),'(A72)',END=102) LINE
          L1 = 0
          L1 = SCAN(LINE,',')
          IF (L1 /= 0) THEN
            BACKSPACE (PRF(JW))
            READ (PRF(JW),'(F8.0)',END=102) JDAY1
          END IF
        END DO
        JDAY1 = 0.0
      END IF
102   CONTINUE
      JDAY1 = 0.0
      IF (CONTOUR(JW)) THEN
        REWIND (CPL(JW))
        DO WHILE (JDAY1 < JDAY)
          READ (CPL(JW),'(A72)',END=103) LINE
          IF (LINE(1:8) == 'New date') THEN
            BACKSPACE (CPL(JW))
            READ (CPL(JW),'(A,F12.4)',END=103) LINE(1:9),JDAY1
          END IF
        END DO
        BACKSPACE (CPL(JW))
        JDAY1 = 0.0
      END IF
103   CONTINUE
      IF (VECTOR(JW)) THEN
        REWIND (VPL(JW))
        DO WHILE (JDAY1 < JDAY)
          READ (VPL(JW),'(A72)',END=104) LINE
          IF (LINE(1:8) == 'New date') THEN
            BACKSPACE (VPL(JW))
            READ (VPL(JW),*,END=104) LINE(1:9),JDAY1
          END IF
        END DO
        BACKSPACE (VPL(JW))
        JDAY1 = 0.0
      END IF
104   CONTINUE
      IF (FLUX(JW)) THEN
        REWIND (FLX(JW))
        DO WHILE (JDAY1 < JDAY)
          READ (FLX(JW),'(A72)',END=105) LINE
          IF (LINE(1:8) == 'New date') THEN
            BACKSPACE (FLX(JW))
            READ (FLX(JW),'(8X,F10.0)',END=105) JDAY1
          END IF
        END DO
        BACKSPACE (FLX(JW))
        REWIND (FLX2(JW))
        DO WHILE (JDAY1 < JDAY)
          READ (FLX2(JW),'(A72)',END=105) LINE
          IF (LINE(1:8) == 'New date') THEN
            BACKSPACE (FLX2(JW))
            READ (FLX2(JW),'(8X,F10.0)',END=105) JDAY1
          END IF
        END DO
        BACKSPACE (FLX2(JW))
      END IF
105   CONTINUE
    END DO
    IF (DOWNSTREAM_OUTFLOW) THEN
      DO JWD=1,NIWDO
        WRITE (SEGNUM,'(I0)') IWDO(JWD)
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
        OPEN   (WDO(JWD,1),FILE='qwo_'//SEGNUM(1:L)//'.opt',POSITION='APPEND')
        REWIND (WDO(JWD,1))
        READ   (WDO(JWD,1),'(//)')
        DO WHILE (JDAY1 < JDAY)
          READ (WDO(JWD,1),'(F8.0)',END=106) JDAY1
        END DO
        BACKSPACE (WDO(JWD,1))
106     CONTINUE
        OPEN   (WDO(JWD,2),FILE='two_'//SEGNUM(1:L)//'.opt',POSITION='APPEND')
        REWIND (WDO(JWD,2))
        READ   (WDO(JWD,2),'(//)')
        DO WHILE (JDAY1 < JDAY)
          READ (WDO(JWD,2),'(F8.0)',END=107) JDAY1
        END DO
        BACKSPACE (WDO(JWD,2))
107     CONTINUE
        IF (CONSTITUENTS) THEN
          OPEN   (WDO(JWD,3),FILE='cwo_'//SEGNUM(1:L)//'.opt',POSITION='APPEND')
          REWIND (WDO(JWD,3))
          READ   (WDO(JWD,3),'(//)')
          DO WHILE (JDAY1 < JDAY)
            READ (WDO(JWD,3),'(F8.0)',END=108) JDAY1
          END DO
          BACKSPACE (WDO(JWD,3))
108       CONTINUE
        END IF
        IF (DERIVED_CALC) THEN
          OPEN   (WDO(JWD,4),FILE='dwo_'//SEGNUM(1:L)//'.opt',POSITION='APPEND')
          REWIND (WDO(JWD,4))
          READ   (WDO(JWD,4),'(//)')
          DO WHILE (JDAY1 < JDAY)
            READ (WDO(JWD,4),'(F8.0)',END=109) JDAY1
          END DO
          BACKSPACE (WDO(JWD,4))
109       CONTINUE
        END IF
      END DO
    END IF
    IF (TIME_SERIES) THEN
      L1 = SCAN(TSRFN,'.')
      DO J=1,NIKTSR
        WRITE (SEGNUM,'(I0)') ITSR(J)
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
        WRITE (SEGNUM2,'(I0)')J
        SEGNUM2 = ADJUSTL(SEGNUM2)
        L2      = LEN_TRIM(SEGNUM2)
        TSRFN  = TSRFN(1:L1-1)//'_'//SEGNUM2(1:L2)//'_seg'//SEGNUM(1:L)//'.opt'
        OPEN   (TSR(J),FILE=TSRFN,POSITION='APPEND')
        REWIND (TSR(J))
        READ   (TSR(J),'(A72)')   (LINE,I=1,11)
        READ   (TSR(J),'(/F10.3)',END=110) JDAYTS
        DO WHILE (JDAYTS < JDAY)
          READ (TSR(J),'(F10.0)',END=110) JDAYTS
        END DO
        BACKSPACE (TSR(J))
110     CONTINUE
      END DO
    END IF
  ELSE
    DO JW=1,NWB
      IF (SNAPSHOT(JW))    OPEN (SNP(JW),FILE=SNPFN(JW),STATUS='UNKNOWN')
      IF (VECTOR(JW))      OPEN (VPL(JW),FILE=VPLFN(JW),STATUS='UNKNOWN')
      IF (PROFILE(JW))     OPEN (PRF(JW),FILE=PRFFN(JW),STATUS='UNKNOWN')
      IF (SPREADSHEET(JW)) OPEN (SPR(JW),FILE=SPRFN(JW),STATUS='UNKNOWN')
      IF (CONTOUR(JW))     OPEN (CPL(JW),FILE=CPLFN(JW),STATUS='UNKNOWN')
      IF (FLUX(JW))THEN
              OPEN (FLX(JW),FILE=FLXFN(JW),STATUS='UNKNOWN')
              WRITE (SEGNUM,'(I0)') JW
              SEGNUM = ADJUSTL(SEGNUM)
              L      = LEN_TRIM(SEGNUM)
              OPEN (FLX2(JW), FILE='kflux_jw'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
              WRITE(FLX2(JW), '("JDAY,  ELTM,",1000(A,","))')(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW))
      ENDIF

!**** Output files

      IF (PROFILE(JW)) THEN
        TTIME = TMSTRT
        DO WHILE (TTIME <= TMEND)
          NDSP = NDSP+1
          TTIME = TTIME+PRFF(PRFDP(JW),JW)
          IF (TTIME >= PRFD(PRFDP(JW)+1,JW)) PRFDP(JW) = PRFDP(JW)+1
        END DO
        PRFDP(JW) = 1
        WRITE (PRF(JW),'(A)')         TITLE
        WRITE (PRF(JW),'(8I8,L2)')    KMX,NIPRF(JW),NDSP,NCT,NDC,NAC+NACD(JW)+1,PRFDP(JW),KTWB(JW),CONSTITUENTS
        WRITE (PRF(JW),'(10I8)')      IPRF(1:NIPRF(JW),JW)
        WRITE (PRF(JW),'(20(1X,A))') ' ON',CPRWBC(:,JW)(6:8),CDWBC(:,JW)(6:8)
        WRITE (PRF(JW),'(2A)')       'Temperature, øC                            ',ADJUSTL(CNAME),ADJUSTL(CDNAME)
        WRITE (PRF(JW),'(20I4)')      1,CN(1:NAC)+1,CDN(1:NACD(JW),JW)+NCT+1
        WRITE (PRF(JW),'(10F8.0)')    1.0,CMULT,CDMULT
        WRITE (PRF(JW),'(20I4)')      KB(IPRF(1:NIPRF(JW),JW))
        WRITE (PRF(JW),'(10F8.2)')    H
        DO JP=1,NIPRF(JW)
          NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
          WRITE (PRF(JW),'(A8,I4/(8(F10.2)))') 'TEMP    ',NRS,(T2(K,IPRF(JP,JW)),K=KTWB(JW),KB(IPRF(JP,JW)))
        END DO
        DO JC=1,NAC
          IF (PRINT_CONST(CN(JC),JW)) THEN
            DO JP=1,NIPRF(JW)
              NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
              WRITE (PRF(JW),'(A,I4/(8(E13.6,1x)))') ADJUSTL(CNAME2(CN(JC))),NRS,(C2(K,IPRF(JP,JW),CN(JC))*CMULT(CN(JC)),K=KTWB(JW),  &                    ! CB 1/25/05
                                                          KB(IPRF(JP,JW)))
            END DO
          END IF
        END DO
        DO JD=1,NACD(JW)
          DO JP=1,NIPRF(JW)
            NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
            WRITE (PRF(JW),'(A,I4/(8(E13.6,1x)))') ADJUSTL(CDNAME2(CDN(JD,JW))),NRS,(CD(K,IPRF(JP,JW),CDN(JD,JW))*CDMULT(CDN(JD,JW)), &                     ! CB 1/25/05
                                                        K=KTWB(JW),KB(IPRF(JP,JW)))
          END DO
        END DO
      END IF
      IF (SPREADSHEET(JW)) THEN
        DO J=1,NISPR(JW)
          WRITE (SEGNUM,'(I0)') ISPR(J,JW)
          SEGNUM = ADJUSTL(SEGNUM)
          L      = LEN_TRIM(SEGNUM)
          SEG(J) = 'Seg_'//SEGNUM(1:L)
        END DO
        WRITE (SPR(JW),'(A,27X,A,5X,A,1000(1X,"Elevation",2X,A7))') 'Constituent','Julian_day','Depth',(SEG(J),J=1,NISPR(JW))
      END IF
      IF (CONTOUR(JW)) THEN
        IF(TECPLOT(JW) /= '      ON')THEN
        WRITE (CPL(JW),'(A)')           TITLE
        WRITE (CPL(JW),'(8(I8,2X))')    NBR
        WRITE (CPL(JW),'(8(I8,2X))')    IMX,KMX
        DO JB=BS(JW),BE(JW)
          WRITE (CPL(JW),'(9(I8,2X))')  US(JB),DS(JB)
          WRITE (CPL(JW),'(9(I8,2X))')  KB(US(JB):DS(JB))
        END DO
        WRITE (CPL(JW),'(8(E13.6,2X))') DLX
        WRITE (CPL(JW),'(8(E13.6,2X))') H
        WRITE (CPL(JW),'(8(I8,2X))')    NAC
        WRITE (CPL(JW),'(A)')           CNAME1(CN(1:NAC))
        ELSE
        !c calculating longitudinal distance of segments
            IF(JW == 1)DIST=0.0
            do jb=BS(JW),BE(JW)
                x1(us(jb))=dist+dlx(us(jb))/2.
                DO I=US(JB)+1,DS(JB)
                    DIST=DIST+(DLX(I)+dlx(i-1))/2.0
                    X1(I)=DIST
                END DO
                DIST=DIST+DLX(DS(JB))
                X1(DS(JB)+1)=DIST
            ENDDO
        WRITE (CPL(JW), *)'TITLE="CE-QUAL-W2"'
		WRITE (CPL(JW),19233)(CNAME2(CN(JN)),JN=1,NAC)
		19233 FORMAT('VARIABLES="Distance, m","Elevation, m","U","W","T","RHO"',<NAC>(',"',A8,'"'))
        ENDIF
      END IF
      IF (VECTOR(JW)) THEN
        WRITE (VPL(JW),*)  TITLE
        WRITE (VPL(JW),*)  H,KB,US,DS,DLX
      END IF
    END DO
    IF (TIME_SERIES) THEN
      L1 = SCAN(TSRFN,'.')
      DO J=1,NIKTSR
        WRITE (SEGNUM,'(I0)') ITSR(J)
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
        WRITE (SEGNUM2,'(I0)')J
        SEGNUM2 = ADJUSTL(SEGNUM2)
        L2      = LEN_TRIM(SEGNUM2)
        TSRFN  = TSRFN(1:L1-1)//'_'//SEGNUM2(1:L2)//'_seg'//SEGNUM(1:L)//'.opt'
        OPEN  (TSR(J),FILE=TSRFN,STATUS='UNKNOWN')
        WRITE (TSR(J),'(A)') (TITLE(I),I=1,11)
        I = ITSR(J)                                                                                            ! SR 5/10/05
        DO JW=1,NWB
          IF (I >= US(BS(JW)) .AND. I <= DS(BE(JW))) EXIT
        END DO
! v3.5 start
        IF (ICE_COMPUTATION) THEN
          if(sediment_calc(jw))then
          WRITE (TSR(J),'(1000(2X,A))') '    JDAY','     DLT','    ELWS','      T2','       U','       Q','    SRON','     EXT',   &
                                        '   DEPTH','   WIDTH','   SHADE','   ICETH',(CNAME2(CN(JC)),JC=1,NAC),                     &
                                        ('     EPI',JE=1,NEP),('     MAC',Jm=1,nmc),'     SED','    SEDP','    SEDN','    SEDC',   &
                                        (CDNAME2(CDN(JD,JW)),JD=1,NACD(JW)),(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW))
          else
            WRITE (TSR(J),'(1000(2X,A))') '    JDAY','     DLT','    ELWS','      T2','       U','       Q','    SRON','     EXT', &
                                        '   DEPTH','   WIDTH','   SHADE','   ICETH',(CNAME2(CN(JC)),JC=1,NAC),                     &
                                        ('     EPI',JE=1,NEP),('     MAC',Jm=1,nmc),(CDNAME2(CDN(JD,JW)),JD=1,NACD(JW)),(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW))
          end if
        ELSE
          if(sediment_calc(jw))then                                                                                                    !mlm 7/25/06
          WRITE (TSR(J),'(1000(2X,A))') '    JDAY','     DLT','    ELWS','      T2','       U','       Q','    SRON','     EXT',   &
                                        '   DEPTH','   WIDTH','   SHADE',(CNAME2(CN(JC)),JC=1,NAC),                     &
                                        ('     EPI',JE=1,NEP),('     MAC',Jm=1,nmc),'     SED','    SEDP','    SEDN','    SEDC',   &
                                        (CDNAME2(CDN(JD,JW)),JD=1,NACD(JW)),(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW))
          else
            WRITE (TSR(J),'(1000(2X,A))') '    JDAY','     DLT','    ELWS','      T2','       U','       Q','    SRON','     EXT', &
                                        '   DEPTH','   WIDTH','   SHADE',(CNAME2(CN(JC)),JC=1,NAC),                     &
                                        ('     EPI',JE=1,NEP),('     MAC',Jm=1,nmc),(CDNAME2(CDN(JD,JW)),JD=1,NACD(JW)),(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW))
          end if
        END IF
! v3.5 end

      END DO
    END IF
    IF (DOWNSTREAM_OUTFLOW) THEN
      DO JWD=1,NIWDO
        WRITE (SEGNUM,'(I0)') IWDO(JWD)
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
        OPEN  (WDO(JWD,1),FILE='qwo_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
        OPEN  (WDO(JWD,2),FILE='two_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
        WRITE (WDO(JWD,1),'(A,I0//A)') 'Flow file for segment ',       IWDO(JWD),'    JDAY     QWD'
        WRITE (WDO(JWD,2),'(A,I0//A)') 'Temperature file for segment ',IWDO(JWD),'    JDAY       T'
        DO JW=1,NWB
          IF (IWDO(JWD) >= US(BS(JW)) .AND. IWDO(JWD) <= DS(BE(JW))) EXIT
        END DO
        IF (CONSTITUENTS) THEN
          OPEN  (WDO(JWD,3),FILE='cwo_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
          WRITE (WDO(JWD,3),'(A,I0//,(100A))') 'Concentration file for segment ',     IWDO(JWD),'    JDAY',CNAME2(CN(1:NAC))
        END IF
        IF (DERIVED_CALC) THEN
          OPEN  (WDO(JWD,4),FILE='dwo_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
          WRITE (WDO(JWD,4),'(A,I0//(100A))') 'Derived constituent file for segment ',IWDO(JWD),'    JDAY',                        &
                                               CDNAME2(CDN(1:NACD(JW),JW))
        END IF
      END DO
    END IF
  END IF

  return
  end subroutine outputinit
