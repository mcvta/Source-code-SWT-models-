subroutine outputa

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  !v3.5
  EXTERNAL RESTART_OUTPUT

  INTEGER :: JAD,JAF,IFLAG,JWWD,NLINES,ITOT,JJ
  REAL    :: TSUM,QSUMM,CGAS,CGASD,TGATE,TSPILL,XDUM
  REAL(R8):: DLVBR,DLE

!***********************************************************************************************************************************
!*                                                    Task 2.8: Output Results                                                    **
!***********************************************************************************************************************************

! v3.5 start  - keep this code for error checking, but delete later
!   if(macrophyte_on)then
!     if(nit.eq.1)open(689,file='totmac.dat',status='unknown')
!     if((nit/200)*200.eq.nit)then
!       bmass=0.0
!       totar=0.0
!       DO Jw=1,Nwb
!         KT = KTwb(Jw)
!         DO JB=BS(Jw),BE(Jw)
!           IU = CUS(JB)
!           ID = DS(JB)
!           DO I=IU,ID
!             smass=0.0
!             do k=kt,kb(i)
!               do m=1,nmc
!                 bmass=bmass+mac(k,i,m)*bh2(k,i)*dlx(i)
!                 smass=smass+mac(k,i,m)*bh2(k,i)*dlx(i)
!               end do
!             end do
!             totar=totar+dlx(i)*b(kti(i),i)
!             armac(i)=smass/(dlx(i)*b(kti(i),i))
!           end do
!         end do
!       end do
!       bmassar=bmass/totar
!       write(689,'(f10.4,2f10.3)')jday,bmass/1000.0,bmassar
!     end if
!   end if
! v3.5 end


!** Time series

    IF (TIME_SERIES) THEN
      IF (JDAY.GE.NXTMTS.OR.JDAY.GE.TSRD(TSRDP+1)) THEN
        IF (JDAY.GE.TSRD(TSRDP+1)) THEN
          TSRDP  = TSRDP+1
          NXTMTS = TSRD(TSRDP)
        END IF
        NXTMTS = NXTMTS+TSRF(TSRDP)
        DO J=1,NIKTSR
          I = ITSR(J)
          DO JW=1,NWB
            IF (I >= US(BS(JW))-1 .AND. I <= DS(BE(JW))+1) EXIT
          END DO
          IF (ETSR(J) < 0) THEN
            K = INT(ABS(ETSR(J)))
          ELSE
            DO K=KTWB(JW),KB(I)
              IF (DEPTHB(K,I) > ETSR(J)) EXIT
            END DO
            IF (K > KB(I)) CYCLE
          END IF
          DO JAC=1,NAC
            L = LEN_TRIM(FMTC(CN(JAC)))
            WRITE (C2CH(JAC),FMTC(CN(JAC))(1:L)) C2(K,I,CN(JAC))*CMULT(CN(JAC))
          END DO
          DO JF=1,NAF(JW)
            WRITE(KFCH(JF),'(E10.3)')KF(K,I,(KFCN(JF,JW)))*VOL(K,I)/1000./DAY
          ENDDO
          DO JAD=1,NACD(JW)
            L = LEN_TRIM(FMTCD(CDN(JAD,JW)))
            WRITE (CDCH(JAD),FMTCD(CDN(JAD,JW))(1:L)) CD(K,I,CDN(JAD,JW))*CDMULT(CDN(JAD,JW))
          END DO
          DO JE=1,NEP
            WRITE (EPCH(JE),'(F10.3)') EPD(K,I,JE)                                                    ! SW 8/13/06
          END DO
          DO Jm=1,Nmc
            WRITE (macCH(Jm),'(F10.3)') mac(K,I,Jm)                                                   ! SW 8/13/06
          END DO
! v3.5 start
          if(sediment_calc(jw))then
            write (sedch,'(F10.3)') sed(K,I)                                                          ! SW 8/13/06
            write (sedpch,'(F10.3)') sedp(K,I)
            write (sednch,'(F10.3)') sedn(K,I)
            write (sedcch,'(F10.3)') sedc(K,I)
          end if
          IF (ICE_COMPUTATION) THEN
            if(sediment_calc(jw))then
              WRITE (TSR(J),'(f10.3,11F10.2,1000A)') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,GAMMA(K,I),DEPTHB(KB(I),I),    &     ! SW 8/13/06
                                             BI(KTWB(JW),I),SHADE(I),ICETH(I),(ADJUSTR(C2CH(JAC)),JAC=1,NAC),                      &     ! CB 7/26/07
                                            (ADJUSTR(EPCH(JE)),JE=1,NEP),(ADJUSTR(macCH(Jm)),Jm=1,Nmc),sedch,sedpch,sednch,sedcch, &
                                            (ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW)),(ADJUSTR(KFCH(JF)),JF=1,NAF(JW))
            else
            WRITE (TSR(J),'(f10.3,11F10.2,1000A)') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,GAMMA(K,I),DEPTHB(KB(I),I),      &     ! SW 8/13/06
                                             BI(KTWB(JW),I),SHADE(I),ICETH(I),(ADJUSTR(C2CH(JAC)),JAC=1,NAC),                      &    ! CB 7/26/07
                                            (ADJUSTR(EPCH(JE)),JE=1,NEP),(ADJUSTR(macCH(Jm)),Jm=1,Nmc),                            &
                                            (ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW)),(ADJUSTR(KFCH(JF)),JF=1,NAF(JW))
            end if
          ELSE
            if(sediment_calc(jw))then
              WRITE (TSR(J),'(f10.3,10F10.2,1000A)') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,GAMMA(K,I),DEPTHB(KB(I),I),    &     ! SW 8/13/06
                                             BI(KTWB(JW),I),SHADE(I),(ADJUSTR(C2CH(JAC)),JAC=1,NAC),(ADJUSTR(EPCH(JE)),            &    ! CB 7/26/07
                                             JE=1,NEP),(ADJUSTR(macCH(Jm)),Jm=1,Nmc),sedch,sedpch,sednch,sedcch,                   &
                                             (ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW)),(ADJUSTR(KFCH(JF)),JF=1,NAF(JW))
            else
            WRITE (TSR(J),'(f10.3,10F10.2,1000A)') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,GAMMA(K,I),DEPTHB(KB(I),I),      &      ! SW 8/13/06
                                             BI(KTWB(JW),I),SHADE(I),(ADJUSTR(C2CH(JAC)),JAC=1,NAC),(ADJUSTR(EPCH(JE)),            &    ! CB 7/26/07
                                             JE=1,NEP),(ADJUSTR(macCH(Jm)),Jm=1,Nmc),(ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW)),    &
                                             (ADJUSTR(KFCH(JF)),JF=1,NAF(JW))
            end if
          END IF
! v3.5 end
        END DO
      END IF
    END IF


    DO JW=1,NWB

!**** Inactive segments

      JB       = BS(JW)
      NBL(JW)  = 1
      IBPR(JW) = 1
      DO I=1,NISNP(JW)-1
        IF (CUS(JB) > ISNP(I,JW)) THEN
          BL(NBL(JW),JW) = I
          NBL(JW)        = NBL(JW)+1
          IBPR(JW)       = I+1
        END IF
        IF (ISNP(I+1,JW) > DS(JB)) JB = JB+1
      END DO
      NBL(JW) = NBL(JW)-1

!**** Snapshots

      IF (SNAPSHOT(JW)) THEN
        IF (JDAY >= NXTMSN(JW) .OR. JDAY >= SNPD(SNPDP(JW)+1,JW)) THEN
          IF (JDAY >= SNPD(SNPDP(JW)+1,JW)) THEN
            SNPDP(JW)  = SNPDP(JW)+1
            NXTMSN(JW) = SNPD(SNPDP(JW),JW)
          END IF
          NXTMSN(JW) = NXTMSN(JW)+SNPF(SNPDP(JW),JW)
          WRITE (SNP(JW),10490) (TITLE(J),J=1,10)
          WRITE (SNP(JW),10500) 'Time Parameters',MONTH,GDAY,YEAR,INT(JDAY),(JDAY-INT(JDAY))*24.0,INT(ELTMJD),                     &
                                (ELTMJD-INT(ELTMJD))*24.0,INT(DLTS1),KLOC,ILOC,INT(MINDLT),INT(JDMIN),(JDMIN-INT(JDMIN))*24.0,     &
                                 KMIN,IMIN
          IF (LIMITING_DLT(JW))  WRITE (SNP(JW),10510) KMIN,IMIN
          WRITE (SNP(JW),10520)  INT(DLTAV),NIT,NV
          WRITE (SNP(JW),10530) 'Meteorological Parameters'
          WRITE (SNP(JW),10540)  TAIR(JW),DEG,TDEW(JW),DEG,PHI(JW),CLOUD(JW),ET(DS(1)),DEG,CSHE(DS(1)),SRON(JW),DEG
          WRITE (SNP(JW),10550) 'Inflows','Upstream inflows'
          DO JB=BS(JW),BE(JW)
            IF (UP_FLOW(JB)) WRITE (SNP(JW),10560) JB,KTQIN(JB),KBQIN(JB),QIN(JB),TIN(JB),DEG
          END DO
          DO JB=BS(JW),BE(JW)
            IF (DIST_TRIBS(JB)) THEN
              WRITE (SNP(JW),10570)
              WRITE (SNP(JW),10580) JB,QDTR(JB),TDTR(JB),DEG
            END IF
          END DO
          IF (TRIBUTARIES) THEN
            WRITE (SNP(JW),10590) (ITR(JT),          JT=1,JTT)
            WRITE (SNP(JW),10600) (KTTR(JT),KBTR(JT),JT=1,JTT)
            WRITE (SNP(JW),10610) (QTR(JT),          JT=1,JTT)
            WRITE (SNP(JW),10620) (TTR(JT),          JT=1,JTT)
          END IF
          WRITE (SNP(JW),10630)
          DO JB=BS(JW),BE(JW)
            IF (DN_FLOW(JB)) THEN
              WRITE (SNP(JW),10640)  JB,(QSTR(JS,JB),JS=1,JSS(JB))
              WRITE (SNP(JW),10650)  QSUM(JB),(K,K=KTWB(JW),KB(DS(JB)))
              WRITE (SNP(JW),10660) (QOUT(K,JB), K=KTWB(JW),KB(DS(JB)))
            END IF
          END DO
          IF (WITHDRAWALS) THEN
            DO JWD=1,JWW
              WRITE (SNP(JW),10670) MAX(CUS(JBWD(JWD)),IWD(JWD)),QWD(JWD)
              IF (QWD(JWD) /= 0.0) THEN
                WRITE (SNP(JW),10680) (K,         K=KTW(JWD),KBW(JWD))
                WRITE (SNP(JW),10690) (QSW(K,JWD),K=KTW(JWD),KBW(JWD))
              ELSE
                WRITE (SNP(JW),10680)
                WRITE (SNP(JW),10690) QWD(JWD)
              END IF
            END DO
          END IF
          IF (CONSTITUENTS) THEN
            WRITE (SNP(JW),10700) 'Constituent Inflow Concentrations'
            DO JB=BS(JW),BE(JW)
              IF (UP_FLOW(JB) .AND. NACIN(JB) > 0)    WRITE (SNP(JW),10710) JB,(CNAME1(INCN(JC,JB))(1:18),CIN(INCN(JC,JB),JB),     &
                                                                            CUNIT2(INCN(JC,JB)),JC=1,NACIN(JB))
              IF (DIST_TRIBS(JB) .AND. NACDT(JB) > 0) WRITE (SNP(JW),10730) JB,(CNAME1(DTCN(JC,JB))(1:18),CDTR(DTCN(JC,JB),JB),    &
                                                                            CUNIT2(DTCN(JC,JB)),JC=1,NACDT(JB))
            END DO
            DO JT=1,NTR
              IF (NACTR(JT) > 0) WRITE (SNP(JW),10720) JT,(CNAME1(TRCN(JC,JT))(1:18),CTR(TRCN(JC,JT),JT),CUNIT2(TRCN(JC,JT)),      &
                                                       JC=1,NACTR(JT))
            END DO
          END IF
          IF (EVAPORATION(JW) .OR. PRECIPITATION(JW)) WRITE (SNP(JW),10740)
          IF (EVAPORATION(JW)) THEN
            WRITE (SNP(JW),10750) (JB,EVBR(JB),JB=BS(JW),BE(JW))                                                              ! SW 9/15/05
            WRITE (SNP(JW),10755) (JB,-VOLEV(JB),JB=BS(JW),BE(JW))
          END IF
          IF (PRECIPITATION(JW)) WRITE (SNP(JW),10760) (JB,PR(JB),JB=BS(JW),BE(JW))
          IF (HEAD_BOUNDARY(JW)) THEN
            WRITE (SNP(JW),10770)
            DO JB=BS(JW),BE(JW)
              IF (UH_EXTERNAL(JB)) WRITE (SNP(JW),10780) JB,ELUH(JB)
              IF (DH_EXTERNAL(JB)) WRITE (SNP(JW),10790) JB,ELDH(JB)
            END DO
          END IF
          IF (VOLUME_BALANCE(JW)) THEN
            WRITE (SNP(JW),10800)
            WRITE (SNP(JW),10810) JW,VOLSR(JW),VOLTR(JW),VOLTR(JW)-VOLSR(JW),DLVR(JW)
            DO JB=BS(JW),BE(JW)
              IF (VOLSBR(JB) /= 0.0) DLVBR = (VOLTBR(JB)-VOLSBR(JB))/VOLSBR(JB)
              WRITE (SNP(JW),10820) JB,VOLSBR(JB),VOLTBR(JB),VOLTBR(JB)-VOLSBR(JB),DLVBR*100.0
            END DO
          END IF
          IF (ENERGY_BALANCE(JW)) THEN
            WRITE (SNP(JW),10830)
            IF (ESR(JW) /= 0.0) DLE = (ESR(JW)-ETR(JW))/ESR(JW)
            WRITE (SNP(JW),10840) JW,ESR(JW)*4.184E3,ETR(JW)*4.184E3,(ESR(JW)-ETR(JW))*4.184E3,DLE*100.0
            DO JB=BS(JW),BE(JW)
              WRITE (SNP(JW),10870) JB
              IF (ESBR(JB) /= 0.0) DLE = (ESBR(JB)-ETBR(JB))/ESBR(JB)
              WRITE (SNP(JW),10850) ESBR(JB)*4.184E3,ETBR(JB)*4.1843E3,(ESBR(JB)-ETBR(JB))*4.1843E3,DLE*100.0
            END DO
          END IF
          IF (MASS_BALANCE(JW)) THEN
            WRITE (SNP(JW),10860)
            DO JB=BS(JW),BE(JW)
              WRITE (SNP(JW),10870) JB
              DO JC=1,NAC
                IF (CMBRS(CN(JC),JB) /= 0.0) DLMR = (CMBRT(CN(JC),JB)-CMBRS(CN(JC),JB))/(CMBRS(CN(JC),JB)+NONZERO)*100.0
                WRITE (SNP(JW),10880) CNAME1(CN(JC)),CMBRS(CN(JC),JB),CUNIT1(CN(JC)),CMBRT(CN(JC),JB),CUNIT1(CN(JC)),              &
                                     (CMBRT(CN(JC),JB)-CMBRS(CN(JC),JB)),CUNIT1(CN(JC)),DLMR
              END DO
! v3.5 start
              DO m=1,Nmc
                if(macrophyte_calc(jw,m))then
                  IF (maCMBRS(JB,m).NE.0.0) THEN
                    DLMR = (maCMBRT(JB,m)-maCMBRS(JB,m))/(maCMBRS(JB,m)+NONZERO)
                  END IF
                  WRITE (SNP(Jw),3312) m,maCMBRS(JB,m),maCMBRT(JB,m),(maCMBRT(JB,m)-maCMBRS(JB,m)),DLMR*100.0
3312  FORMAT(5X,'Macrophyte spec ',i2,/7X,'Spatially integrated mass [MACMBRS] = ',1PE15.8E2,1X,'g ',/7X,   &
                'Temporally integrated mass [MACMBRT] = ',1PE15.8E2,1X,'g ',/7X,'Mass error                         = ',  &
                 1PE15.8E2,1X,'g ',/7X,'Percent error                      = ',1PE15.8E2,' %')
                end if
              END DO
! v3.5 end
            END DO
          END IF
          WRITE (SNP(JW),10890) 'Geometry',KTWB(JW),ELKT(JW)
          WRITE (SNP(JW),10900) (JB,CUS(JB),JB=BS(JW),BE(JW))
          CALL OUTPUT (JDAY,IBPR(JW),NISNP(JW),KBR(JW),ISNP,BL(1,JW),NBL(JW))
        END IF
      END IF

!**** Vertical profiles

      IF (PROFILE(JW)) THEN
        IF (JDAY >= NXTMPR(JW) .OR. JDAY >= PRFD(PRFDP(JW)+1,JW)) THEN
          IF (JDAY >= PRFD(PRFDP(JW)+1,JW)) THEN
            PRFDP(JW)  = PRFDP(JW)+1
            NXTMPR(JW) = PRFD(PRFDP(JW),JW)
          END IF
          NXTMPR(JW) = NXTMPR(JW)+PRFF(PRFDP(JW),JW)
          NSPRF(JW)  = NSPRF(JW)+1
          WRITE (PRF(JW),'(F8.3,1X,A3,I3,A,2I4,F8.4,I8)')JDAY,ADJUSTL(MONTH),GDAY,', ',YEAR,KTWB(JW),SNGL(Z(DS(BS(JW)))),NSPRF(JW)
          DO JP=1,NIPRF(JW)
            NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
            WRITE (PRF(JW),'(A8,I4/(8F10.2))') 'TEMP    ',NRS,(T2(K,IPRF(JP,JW)),K=KTWB(JW),KB(IPRF(JP,JW)))
          END DO
          DO JC=1,NAC
            IF (PRINT_CONST(CN(JC),JW)) THEN
              DO JP=1,NIPRF(JW)
                NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
                WRITE (PRF(JW),'(A8,I4/(8(E13.6,2X)))') ADJUSTL(CNAME2(CN(JC))),NRS,(C2(K,IPRF(JP,JW),CN(JC))*CMULT(CN(JC)),       &     ! CB 1/24/05
                                                            K=KTWB(JW),KB(IPRF(JP,JW)))
              END DO
            END IF
          END DO
          IF (CONSTITUENTS) THEN
            DO JD=1,NACD(JW)
              DO JP=1,NIPRF(JW)
                NRS = KB(IPRF(JP,JW))-KTWB(JW)+1
                WRITE (PRF(JW),'(A8,I4/(8(E13.6,2X)))') ADJUSTL(CDNAME2(CDN(JD,JW))),NRS,(CD(K,IPRF(JP,JW),CDN(JD,JW))             &      ! CB 1/24/05
                                                           *CDMULT(CDN(JD,JW)),K=KTWB(JW), KB(IPRF(JP,JW)))
              END DO
            END DO
          END IF
        END IF
      END IF

!**** Spreadsheet

      IF (SPREADSHEET(JW)) THEN
        IF (JDAY >= NXTMSP(JW) .OR. JDAY >= SPRD(SPRDP(JW)+1,JW)) THEN
          IF (JDAY >= SPRD(SPRDP(JW)+1,JW)) THEN
            SPRDP(JW)  = SPRDP(JW)+1
            NXTMSP(JW) = SPRD(SPRDP(JW),JW)
          END IF
          CONV1      = BLANK1
          NXTMSP(JW) = NXTMSP(JW)+SPRF(SPRDP(JW),JW)
          DO J=1,NISPR(JW)
            KBMAX(JW) = MAX(KB(ISPR(J,JW)),KBMAX(JW))
            DO K=KTWB(JW),KB(ISPR(J,JW))
              WRITE (CONV1(K,J),'(F10.2)') T2(K,ISPR(J,JW))
            END DO
          END DO
          DO K=KTWB(JW),KBMAX(JW)
            WRITE (SPR(JW),'(A18,20X,2F10.3,1000(F10.3,A))') 'Temperature       ',JDAY,DEPTHM(K,DS(BS(JW))),                      &
                                                             (ELWS(ISPR(J,JW))-DEPTHM(K,ISPR(J,JW)),CONV1(K,J),J=1,NISPR(JW))
          END DO
          DO JC=1,NAC
            IF (PRINT_CONST(CN(JC),JW)) THEN
              DO J=1,NISPR(JW)
                DO K=KTWB(JW),KB(ISPR(J,JW))
                  WRITE (CONV1(K,J),'(F10.3)') C2(K,ISPR(J,JW),CN(JC))*CMULT(CN(JC))                                                ! SW 8/13/06
                END DO
              END DO
              DO K=KTWB(JW),KBMAX(JW)
                WRITE (SPR(JW),'(A38,2F10.3,1000(F10.3,A))') CNAME3(CN(JC)),JDAY,DEPTHM(K,DS(BS(JW))),                            &
                                                            (ELWS(ISPR(J,JW))-DEPTHM(K,ISPR(J,JW)),CONV1(K,J),J=1,NISPR(JW))
              END DO
            END IF
          END DO
          IF (CONSTITUENTS) THEN
            DO JD=1,NACD(JW)
              IF (PRINT_DERIVED(CDN(JD,JW),JW)) THEN
                DO J=1,NISPR(JW)
                  DO K=KTWB(JW),KB(ISPR(J,JW))
                    WRITE (CONV1(K,J),'(F10.3)') CD(K,ISPR(J,JW),CDN(JD,JW))*CDMULT(CDN(JD,JW))                                     ! SW 8/13/06
                  END DO
                END DO
                DO K=KTWB(JW),KBMAX(JW)
                  WRITE (SPR(JW),'(A38,2F10.3,1000(F10.3,A))') CDNAME3(CDN(JD,JW)),JDAY,DEPTHM(K,DS(BS(JW))),                     &
                                                              (ELWS(ISPR(J,JW))-DEPTHM(K,ISPR(J,JW)),CONV1(K,J),J=1,NISPR(JW))
                END DO
              END IF
            END DO
          END IF
        END IF
      END IF

!**** Velocity vectors

      IF (VECTOR(JW)) THEN
        IF (JDAY >= NXTMVP(JW) .OR. JDAY >= VPLD(VPLDP(JW)+1,JW)) THEN
          IF (JDAY >= VPLD(VPLDP(JW)+1,JW)) THEN
            VPLDP(JW)  = VPLDP(JW)+1
            NXTMVP(JW) = VPLD(VPLDP(JW),JW)
          END IF
          NXTMVP(JW) = NXTMVP(JW)+VPLF(VPLDP(JW),JW)
          WRITE (VPL(JW),*)  'New date ',JDAY,MONTH//GDCH//',',YEAR,KTWB(JW),(US(JB),JB=BS(JW),BE(JW))
          WRITE (VPL(JW),*) ((Z(I)*COSA(BS(JW))),     I=US(BS(JW)),DS(BE(JW)))
          WRITE (VPL(JW),*) ((EL(K,I),K=KTWB(JW),KMX),I=US(BS(JW)),DS(BE(JW)))
          WRITE (VPL(JW),*) ((U(K,I), K=KTWB(JW),KMX),I=US(BS(JW)),DS(BE(JW)))
          WRITE (VPL(JW),*) ((W(K,I), K=KTWB(JW),KMX),I=US(BS(JW)),DS(BE(JW)))
        END IF
      END IF

!**** Contours

      IF (CONTOUR(JW)) THEN
        IF (JDAY >= NXTMCP(JW) .OR. JDAY >= CPLD(CPLDP(JW)+1,JW)) THEN
          IF (JDAY >= CPLD(CPLDP(JW)+1,JW)) THEN
            CPLDP(JW)  = CPLDP(JW)+1
            NXTMCP(JW) = CPLD(CPLDP(JW),JW)
          END IF
          NXTMCP(JW) = NXTMCP(JW)+CPLF(CPLDP(JW),JW)

         IF(TECPLOT(JW) /= '      ON')THEN
          WRITE (CPL(JW),'(A,F12.4,5X,A9,5X,I2,5X,I4)') 'New date ',JDAY,MONTH,GDAY,YEAR
          WRITE (CPL(JW),'(9(I8,2X))')                   KTWB(JW)
          WRITE (CPL(JW),'(9(E13.6,2X))')               (QTR(JT),JT=1,NTR)
          WRITE (CPL(JW),'(9(E13.6,2X))')               (TTR(JT),JT=1,NTR)
          DO JT=1,NTR
            DO JAC=1,NACTR(JT)
              IF (PRINT_CONST(TRCN(JAC,JT),JW)) WRITE (CPL(JW),'(9(E13.6,2X))') CTR(TRCN(JAC,JT),JT)
            END DO
          END DO
          DO JB=BS(JW),BE(JW)
            WRITE (CPL(JW),'(9(I8,2X))')             CUS(JB)
            WRITE (CPL(JW),'(9(E13.6,2X))')          QIN(JB), QSUM(JB)
            DO I=CUS(JB),DS(JB)
              WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'BHR', (BHR1(K,I),K=KTWB(JW)+1,KB(I))
            END DO
            DO I=CUS(JB),DS(JB)
              WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'U',   (U(K,I),   K=KTWB(JW),KB(I))
            END DO
            WRITE (CPL(JW),'(A38/(9(E13.6,2X)))')   'QC',  (QC(I),    I=CUS(JB),DS(JB))
            WRITE (CPL(JW),'(A38/(9(E13.6,2X)))')   'Z',   (Z(I),     I=CUS(JB),DS(JB))
            WRITE (CPL(JW),'(A38/(9(I8,2X)))')   'KTI',   (kti(I),     I=CUS(JB),DS(JB))  ! v3.5
            DO I=CUS(JB),DS(JB)
              WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'Temperature',(T2(K,I),K=KTWB(JW),KB(I))
            END DO
            DO JC=1,NAC
              IF (PRINT_CONST(CN(JC),JW)) THEN
                DO I=CUS(JB),DS(JB)
                  WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') CNAME(CN(JC)),(C2(K,I,CN(JC))*CMULT(CN(JC)),K=KTWB(JW),KB(I))
                END DO
              END IF
            END DO
            DO JE=1,NEP
              DO I=CUS(JB),DS(JB)
                IF (PRINT_EPIPHYTON(JW,JE)) WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'Epiphyton',(EPD(K,I,JE),K=KTWB(JW),KB(I))
              END DO
            END DO
! v3.5 start
            if(print_sediment(jw))then
              DO I=CUS(JB),DS(JB)
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment',(seD(K,I),K=KTWB(JW),KB(I))
              end do
              DO I=CUS(JB),DS(JB)
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment P',(seDp(K,I),K=KTWB(JW),KB(I))
              end do
              DO I=CUS(JB),DS(JB)
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment N',(seDn(K,I),K=KTWB(JW),KB(I))
              end do
              DO I=CUS(JB),DS(JB)
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment C',(seDc(K,I),K=KTWB(JW),KB(I))
              end do
            end if
            do m=1,nmc
              IF (print_macrophyte(jw,m)) THEN
                DO I=CUS(JB),DS(JB)
                  WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Macrophytes',((macrc(j,K,I,m),j=kti(i),kb(i)), K=KTwb(Jw),KB(I))
                end do
              end if
            end do
! v3.5 end
          IF (CONSTITUENTS) THEN
            DO JD=1,NACD(JW)
              IF (PRINT_DERIVED(CDN(JD,JW),JW)) THEN
                  WRITE (CPL(JW),'(A38/(9(F10.3,2X)))') CDNAME(CDN(JD,JW)),((CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)),             &        ! SW 8/12/06
                                                        K=KTWB(JW),KB(I)),I=CUS(JB),DS(JB))                                              ! CB 1/03/05
              END IF
            END DO
          END IF
          END DO
        ELSE
 !         ICPL=ICPL+1
          ICPL(jw)=ICPL(jw)+1                 ! cb 1/26/09
          itot=0
 !         do jb=1,nbr
          do jb=bs(jw),be(jw)                 ! cb 1/26/09
          itot=itot+ds(jb)-cus(jb)+2
          enddo
           WRITE (CPL(JW),9864)JDAY,KMX-KTWB(JW)+2,itot
9864       FORMAT('ZONE T="',f9.3,'"',' I=',I3,' J=',I3,' F=POINT')
		  do jb=bs(jw),be(jw)
			do i=cus(jb), ds(jb)+1			
			K=KTWB(JW)    ! PRINT AN EXTRA LINE FOR THE SURFACE
			if(i /= ds(jb)+1)then
			write (CPL(JW),9999) x1(i),ELWS(I),U(k,i),-W(k,i),T1(k,i),RHO(k,i),(C2(K,I,CN(JC)),JC=1,NAC)
			else
			xdum=-99.0
			write (CPL(JW),9999) x1(i),ELWS(I),xdum,xdum,xdum,xdum,(XDUM, JJ=1,NAC)
			endif
				do k=ktwb(jw),kmx-1
				    if(i /= ds(jb)+1 .AND. k <= kb(i))then
					write (CPL(JW),9999) x1(i),ELWS(I)-DEPTHM(K,I),U(k,i),-W(k,i),T1(k,i),rho(k,i),(C2(K,I,CN(JC)),JC=1,NAC)
					    IF(K == KB(I))THEN
					    write (CPL(JW),9999) x1(i),ELWS(I)-DEPTHB(K,I),U(k,i), -W(k,i),T1(k,i),rho(k,i),(C2(K,I,CN(JC)),JC=1,NAC)
					    ENDIF
					else
					xdum=-99.0
					write (CPL(JW),9999) x1(i),ELWS(I-1)-DEPTHM(K,I-1),xdum,xdum,xdum,xdum,(XDUM, JJ=1,NAC)
					    IF(K == KB(I))THEN
					    write (CPL(JW),9999) x1(i),ELWS(I-1)-DEPTHB(K,I-1),xdum,xdum,xdum,xdum,(XDUM, JJ=1,NAC)
					    ENDIF
					endif
              end do
              end do
            end do
 !          WRITE(CPL(JW),9899)ICPL,m,gday,year
           WRITE(CPL(JW),9899)ICPL(jw),m,gday,year                                                                        ! cb 1/26/09
9899              FORMAT('TEXT X=0.75, y=0.85, H=2.8,ZN=',i4,',',' C=BLACK,','T= "',i2,'/',i2,'/',i4,'"')
 !          WRITE(CPL(JW),9863)ICPL,jday
           WRITE(CPL(JW),9863)ICPL(jw),jday                                                                                 ! cb 1/26/09
9863              FORMAT('TEXT X=0.75, y=0.90, H=2.8,ZN=',i4,',',' C=BLACK,','T= "Julian Day ',f9.3,'"')
		9999 FORMAT (100(e13.6,1x))
        ENDIF
        END IF
      END IF

!**** Fluxes   KF is the instantaneous flux in g/m3/s, KFS is the summed flux in g eventually divided by elapsed time between calls to FLUX output and converted to kg below

      IF (FLUX(JW)) THEN
        IF (JDAY >= NXTMFL(JW) .OR. JDAY >= FLXD(FLXDP(JW)+1,JW)) THEN
          IF (JDAY >= FLXD(FLXDP(JW)+1,JW)) THEN
            FLXDP(JW)  = FLXDP(JW)+1
            NXTMFL(JW) = FLXD(FLXDP(JW),JW)
          END IF
          NXTMFL(JW) = NXTMFL(JW)+FLXF(FLXDP(JW),JW)
          CONV       = BLANK
          DO JAF=1,NAF(JW)
            DO JB=BS(JW),BE(JW)
              DO I=CUS(JB),DS(JB)
                DO K=KTWB(JW),KB(I)
                  KFS(K,I,KFCN(JAF,JW)) = KFS(K,I,KFCN(JAF,JW))/ELTMF(JW)*DAY
                  KFJW(JW,KFCN(JAF,JW)) = KFJW(JW,KFCN(JAF,JW))+KFS(K,I,KFCN(JAF,JW))/1000.   ! SUM UP FOR ENTIRE WATERBODY
                END DO
              END DO
            END DO
            DO I=1,NISNP(JW)
              DO K=KTWB(JW),KB(ISNP(I,JW))
                WRITE (CONV(K,I),'(E10.3)') KFS(K,ISNP(I,JW),KFCN(JAF,JW))/1000.0
              END DO
            END DO
            IF (NEW_PAGE) THEN
              WRITE (FLX(JW),'(/(A72))') (TITLE(J),J=1,11)
              NLINES   =  KMX-KTWB(JW)+14
              NEW_PAGE = .FALSE.
            END IF
            NLINES   = NLINES+KMX-KTWB(JW)+11
            NEW_PAGE = NLINES > 72
            WRITE (FLX(JW),'(/A,F10.3,X,3(A,I0),A,F0.2,A/)') 'New date ',JDAY,MONTH//' ',GDAY,', ',YEAR,'   Julian Date = ',       &
                                                              INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,                            &
                                                            ' hours           '//KFNAME(KFCN(JAF,JW))
            WRITE (FLX(JW),'(3X,2000I10)')                  (ISNP(I,JW),I=1,NISNP(JW))
            DO K=KTWB(JW),KBR(JW)
              WRITE (FLX(JW),'(1X,I2,200A)') K,(CONV(K,I),I=1,NISNP(JW))
            END DO
          END DO
          WRITE(FLX2(JW),'(F10.3,",",f8.3,",",1000(E12.4,","))')JDAY,ELTMF(JW)/DAY,(KFJW(JW,KFCN(K,JW)),K=1,NAF(JW))
          ELTMF(JW)                    = 0.0
          KF(:,CUS(BS(JW)):DS(BE(JW)),KFCN(1:NAF(JW),JW))  = 0.0
          KFS(:,CUS(BS(JW)):DS(BE(JW)),KFCN(1:NAF(JW),JW)) = 0.0
          KFJW(JW,KFCN(1:NAF(JW),JW))  = 0.0
        END IF
      END IF

END DO

!** Downstream flow, temperature, and constituent files

    IF (DOWNSTREAM_OUTFLOW) THEN
      IF (JDAY >= NXTMWD .OR. JDAY >= WDOD(WDODP+1)) THEN
        IF (JDAY >= WDOD(WDODP+1)) THEN
          WDODP  = WDODP+1
          NXTMWD = WDOD(WDODP)
        END IF
        NXTMWD = NXTMWD+WDOF(WDODP)
        DO J=1,NIWDO
          QWDO(J)    = 0.0
          TWDO(J)    = 0.0
          CWDO(:,J)  = 0.0
          CDWDO(:,J) = 0.0
          CDTOT      = 0.0
          DO JWD=1,JWW
            IF (IWD(JWD) == IWDO(J) .AND. QWD(JWD) /= 0.0 .AND. ILAT(JWD) == 0) THEN
              TSUM  = 0.0
              QSUMM = 0.0
              CSUM  = 0.0
              CDSUM = 0.0
              DO JWWD=1,NWB
                IF (JBWD(JWD) >= BS(JWWD) .AND. JBWD(JWD) <= BE(JWWD)) EXIT
              END DO
              DO K=KTW(JWD),KBW(JWD)
                QSUMM = QSUMM+QSW(K,JWD)
                TSUM  = TSUM+T2(K,IWD(JWD))*QSW(K,JWD)
                DO JC=1,NAC
                  if (cn(jc) .ne. ndo) then                 !cb 09/28/04
!				  IF (CN(JC) == NDO) THEN
!                    IFLAG = 0
!                    DO JS=1,NSP
!                      IF (TDG_SPILLWAY(JWD,JS)) THEN
!                        CALL TOTAL_DISSOLVED_GAS (PALT(IWD(JWD)),0,JS,T2(K,IWD(JWD)),CGAS)
!                        CGASD        = (CGAS/EXP(7.7117-1.31403*(LOG(T2(K,IWD(JWD))+45.93)))*PALT(IWD(JWD)))*100.0
!                        CDSUM(NDG)   =  CDSUM(NDG)+CGASD*QSW(K,JWD)
!                        CSUM(CN(JC)) =  CSUM(CN(JC))+CGAS*QSW(K,JWD)
!                        IFLAG        =  1; EXIT
!                      END IF
!                    END DO
!                    IF (IFLAG == 0) THEN
!                      DO JG=1,NGT
!                        IF (TDG_GATE(JWD,JG)) THEN
!                          CALL TOTAL_DISSOLVED_GAS (PALT(IWD(JWD)),1,JG,T2(K,IWD(JWD)),CGAS)
!                          CGASD        = (CGAS/EXP(7.7117-1.31403*(LOG(T2(K,IWD(JWD))+45.93)))*PALT(IWD(JWD)))*100.0
!                          CDSUM(NDG)   =  CDSUM(NDG)+CGASD*QSW(K,JWD)
!                          CSUM(CN(JC)) =  CSUM(CN(JC))+CGAS*QSW(K,JWD)
!                          IFLAG        =  1; EXIT
!                        END IF
!                      END DO
!                    END IF
!                    IF (IFLAG == 0) CSUM(CN(JC)) = CSUM(CN(JC))+C2(K,IWD(JWD),CN(JC))*QSW(K,JWD)
!                  ELSE
                    CSUM(CN(JC)) = CSUM(CN(JC))+C2(K,IWD(JWD),CN(JC))*QSW(K,JWD)
                  END IF
                END DO
                DO JC=1,NACD(JWWD)
                  if (cdn(jc,jwwd) .ne. ndg) then                !cb 09/28/04
!                  IF (CDN(JC,JWWD) == NDG) THEN
!                    IFLAG = 0
!                    DO JG=1,NGT
!                      IF (TDG_GATE(JWD,JG)) THEN
!                        IFLAG = 1; EXIT
!                      END IF
!                    END DO
!                    DO JS=1,NSP
!                      IF (TDG_SPILLWAY(JWD,JS)) THEN
!                        IFLAG = 1; EXIT
!                      END IF
!                    END DO
!                    IF (IFLAG == 0) CDSUM(CDN(JC,JWWD)) = CDSUM(CDN(JC,JWWD))+CD(K,IWD(JWD),CDN(JC,JWWD))*QSW(K,JWD)
!                  ELSE
                    CDSUM(CDN(JC,JWWD)) = CDSUM(CDN(JC,JWWD))+CD(K,IWD(JWD),CDN(JC,JWWD))*QSW(K,JWD)
                  END IF
                END DO
              END DO

              iflag = 0                                                                                          !cb 09/28/04
              do js=1,nsp                                                                                        !cb 09/28/04
                if (tdg_spillway(jwd,js)) then                                                                   !cb 09/28/04
                  tspill=tsum/qsp(js)                                                                            !cb 09/28/04
                  call total_dissolved_gas (PALT(IWD(JWD)),0,js,tspill,cgas)                                     !cb 09/28/04
                  cgasd        = (CGAS/EXP(7.7117-1.31403*(LOG(tspill+45.93)))*PALT(IWD(JWD)))*100.0     !cb 09/28/04
                  cdsum(ndg)   =  cdsum(ndg)+cgasd*qsp(js)                                                       !cb 09/28/04
                  csum(ndo) =  csum(ndo)+cgas*qsp(js)                                                            !cb 09/28/04
                  iflag        =  1; exit                                                                        !cb 09/28/04
                end if                                                                                           !cb 09/28/04
              end do                                                                                             !cb 09/28/04
              if (iflag == 0) then                                                                               !cb 09/28/04
                do jg=1,ngt                                                                                      !cb 09/28/04
                  if (tdg_gate(jwd,jg)) then                                                                     !cb 09/28/04
                    tgate=tsum/qgt(jg)                                                                           !cb 09/28/04
                    call total_dissolved_gas (PALT(IWD(JWD)),1,jg,tgate,cgas)                                    !cb 09/28/04
                    cgasd        = (cgas/exp(7.7117-1.31403*(log(tgate+45.93)))*PALT(IWD(JWD)))*100.0                      !cb 09/28/04
                    cdsum(ndg)   =  cdsum(ndg)+cgasd*qgt(jt)                                                     !cb 09/28/04
                    csum(ndo) =  csum(ndo)+cgas*qgt(jg)                                                          !cb 09/28/04
                    iflag        =  1; exit                                                                      !cb 09/28/04
                  end if                                                                                         !cb 09/28/04
                end do                                                                                           !cb 09/28/04
              end if                                                                                             !cb 09/28/04
              if (iflag == 0 .and. cac(ndo)    == '      ON')then                                                !cb 09/28/04
                do k=ktw(jwd),kbw(jwd)                                                                           !cb 09/28/04
                  csum(ndo) = csum(ndo)+c2(k,iwd(jwd),ndo)*qsw(k,jwd)                                            !cb 09/28/04
                end do                                                                                           !cb 09/28/04
              end if                                                                                             !cb 09/28/04
			
              iflag = 0                                                                                          !cb 09/28/04
              do jg=1,ngt                                                                                        !cb 09/28/04
                if (tdg_gate(jwd,jg)) then                                                                       !cb 09/28/04
                  iflag = 1; exit                                                                                !cb 09/28/04
                end if                                                                                           !cb 09/28/04
              end do                                                                                             !cb 09/28/04
              do js=1,nsp                                                                                        !cb 09/28/04
                if (tdg_spillway(jwd,js)) then                                                                   !cb 09/28/04
                  iflag = 1; exit                                                                                !cb 09/28/04
                end if                                                                                           !cb 09/28/04
              end do                                                                                             !cb 09/28/04
              if (iflag == 0 .and. cdwbc(ndg,jwwd) == '      ON' ) then                                     !cb 09/28/04
                do k=ktw(jwd),kbw(jwd)                                                                    !cb 09/28/04
                  cdsum(ndg) = cdsum(ndg)+cd(k,iwd(jwd),ndg)*qsw(k,jwd)                                   !cb 09/28/04
                end do                                                                                    !cb 09/28/04
              end if                                                                                      !cb 09/28/04

              QWDO(J)                         = QWDO(J)                    +QSUMM
              TWDO(J)                         = TSUM                       +TWDO(J)
              CWDO(CN(1:NAC),J)               = CSUM(CN(1:NAC))            +CWDO(CN(1:NAC),J)
              CDWDO(CDN(1:NACD(JWWD),JWWD),J) = CDSUM(CDN(1:NACD(JWWD),JWWD))+CDWDO(CDN(1:NACD(JWWD),JWWD),J)
            END IF
          END DO
          DO JW=1,NWB
            DO JB=BS(JW),BE(JW)
              IF (DS(JB) == IWDO(J)) THEN
                QWDO(J)           = QWDO(J)            +QSUM(JB)
                TWDO(J)           = (TOUT(JB)          *QSUM(JB))+TWDO(J)
                CWDO(CN(1:NAC),J) = (COUT(CN(1:NAC),JB)*QSUM(JB))+CWDO(CN(1:NAC),J)
                DO K=KTWB(JW),KB(DS(JB))
                  CDTOT(CDN(1:NACD(JW),JW)) = CDTOT(CDN(1:NACD(JW),JW))+CD(K,DS(JB),CDN(1:NACD(JW),JW))*QOUT(K,JB)
                END DO
                CDWDO(CDN(1:NACD(JW),JW),J) = CDWDO(CDN(1:NACD(JW),JW),J)+CDTOT(CDN(1:NACD(JW),JW))
              END IF
            END DO
          END DO
          IF (QWDO(J) /= 0.0) TWDO(J) = TWDO(J)/QWDO(J)
          DO JC=1,NAC
            IF (QWDO(J) /= 0.0) CWDO(CN(JC),J) = CWDO(CN(JC),J)/QWDO(J)
            WRITE (CWDOC(CN(JC)),'(G8.3)') CWDO(CN(JC),J)
            CWDOC(CN(JC)) = ADJUSTR(CWDOC(CN(JC)))
          END DO
          DO JW=1,NWB
            IF (IWDO(J) >= US(BS(JW)) .AND. IWDO(J) <= DS(BE(JW))) EXIT
          END DO
          DO JD=1,NACD(JW)
            IF (QWDO(J) /= 0.0) CDWDO(CDN(JD,JW),J) = CDWDO(CDN(JD,JW),J)/QWDO(J)
            WRITE (CDWDOC(CDN(JD,JW)),'(G8.3)') CDWDO(CDN(JD,JW),J)
            CDWDOC(CDN(JD,JW)) = ADJUSTR(CDWDOC(CDN(JD,JW)))
          END DO
          if(jday<10000.)then
          WRITE (WDO(J,1),'(F8.3,F8.2)') JDAY, QWDO(J)
          WRITE (WDO(J,2),'(F8.3,F8.2)') JDAY, TWDO(J)
          IF (CONSTITUENTS) WRITE (WDO(J,3),'(F8.3,1000A8)') JDAY,(CWDOC(CN(JC)),     JC=1,NAC)
          IF (DERIVED_CALC) WRITE (WDO(J,4),'(F8.3,1000A8)') JDAY,(CDWDOC(CDN(JD,JW)),JD=1,NACD(JW))
          else
          WRITE (WDO(J,1),'(F8.2,F8.2)') JDAY, QWDO(J)
          WRITE (WDO(J,2),'(F8.2,F8.2)') JDAY, TWDO(J)
          IF (CONSTITUENTS) WRITE (WDO(J,3),'(F8.2,1000A8)') JDAY,(CWDOC(CN(JC)),     JC=1,NAC)
          IF (DERIVED_CALC) WRITE (WDO(J,4),'(F8.2,1000A8)') JDAY,(CDWDOC(CDN(JD,JW)),JD=1,NACD(JW))
          endif
        END DO
      END IF
    END IF

!** Restart

    IF (RESTART_OUT) THEN
      IF (JDAY >= NXTMRS .OR. JDAY >= RSOD(RSODP+1)) THEN
        IF (JDAY >= RSOD(RSODP+1)) THEN
          RSODP  = RSODP+1
          NXTMRS = RSOD(RSODP)
        END IF
        NXTMRS = NXTMRS+RSOF(RSODP)
        WRITE (EXT,'(I0)') INT(JDAY)
        EXT   = ADJUSTL(EXT)
        L     = LEN_TRIM(EXT)
        RSOFN = 'rso'//EXT(1:L)//'.opt'
        CALL RESTART_OUTPUT (RSOFN)
      END IF
    END IF

! Snapshot formats

10490 FORMAT ('CE-QUAL-W2 V3.6'/                                                                                                   &
             (1X,A72))
10500 FORMAT (/1X,A/                                                                                                               &
              3X,'Gregorian date      [GDAY] =',A19,1X,I0,', ',I0/                                                                 &
              3X,'Julian date         [JDAY] =',I10,' days',F6.2,' hours'/                                                         &
              3X,'Elapsed time      [ELTMJD] =',I10,' days',F6.2,' hours'/                                                         &
              3X,'Timestep             [DLT] =',I10,' sec'/                                                                        &
              3X,'  at location  [KLOC,ILOC] = (',I0,',',I0,')'/                                                                   &
              3X,'Minimum timestep  [MINDLT] =',I10,' sec '/                                                                       &
              3X,'  at Julian day    [JDMIN] =',I10,' days',F6.2,' hours'/                                                         &
              3X,'  at location  [KMIN,IMIN] = (',I0,',',I0,')')
10510 FORMAT (3X,'Limiting timestep'/                                                                                              &
              3X,'  at location  [KMIN,IMIN] = (',I0,',',I0,')')
10520 FORMAT (3X,'Average timestep   [DLTAV] =',I10,' sec'/                                                                        &
              3X,'Number of iterations [NIT] =',I10/                                                                               &
              3X,'Number of violations  [NV] =',I10/)
10530 FORMAT (1X,A)
10540 FORMAT (3X,'Input'/                                                                                                          &
              3X,'  Air temperature          [TAIR] =',F9.2,1X,A/                                                                  &
              3X,'  Dewpoint temperature     [TDEW] =',F9.2,1X,A/                                                                  &
              3X,'  Wind direction            [PHI] =',F9.2,' rad'/                                                                &
              3X,'  Cloud cover             [CLOUD] =',F9.2/                                                                       &
              3X,'  Calculated'/                                                                                                   &
              5X,'  Equilibrium temperature    [ET] =',F9.2,1X,A/                                                                  &
              5X,'  Surface heat exchange    [CSHE] =',E9.2,' m/sec'/                                                              &
              5X,'  Net short wave radiation [SRON] =',E9.2,1X,A,' W/m^2'/)
10550 FORMAT (1X,A/                                                                                                                &
              3X,A)
10560 FORMAT (5X,'Branch ',I0/                                                                                                     &
              5X,'  Layer       [KQIN] = ',I0,'-',I0/                                                                              &
              5X,'  Inflow       [QIN] =',F8.2,' m^3/sec'/                                                                         &
              5X,'  Temperature  [TIN] =',F8.2,1X,A)
10570 FORMAT (/3X,'Distributed Tributaries')
10580 FORMAT (5X,'Branch ',I0/                                                                                                     &
              5X,'  Inflow      [QDTR] =',F8.2,' m^3/sec'/                                                                         &
              5X,'  Temperature [TDTR] =',F8.2,1X,A)
10590 FORMAT (:/3X,'Tributaries'/                                                                                                  &
              5X,'Segment     [ITR] =',11I8:/                                                                                      &
             (T25,11I8))
10600 FORMAT (:5X,'Layer      [KTWB] = ',11(I0,'-',I0,2X):/                                                                        &
             (T25,11(I0,'-',I0)))
10610 FORMAT (:5X,'Inflow      [QTR] =',11F8.2:/                                                                                   &
             (T25,11F8.1))
10620 FORMAT (:5X,'Temperature [TTR] =',11F8.2:/                                                                                   &
             (T25,11F8.1))
10630 FORMAT (/1X,'Outflows')
10640 FORMAT (3X,'Structure outflows [QSTR]'/                                                                                      &
              3X,'  Branch ',I0,' = ',11F8.2:/                                                                                     &
             (T16,11F8.2))
10650 FORMAT (:/3X,'Total outflow [QOUT] =',F8.2,' m^3/s'/                                                                         &
              5X,'Outlets'/                                                                                                        &
              5X,'  Layer             [KOUT] =',12I7:/                                                                             &
             (33X,12I7))
10660 FORMAT (:7X,'Outflow (m^3/sec) [QOUT] =',12F7.2:/                                                                            &
             (33X,12F7.2))
10670 FORMAT (:5X,'Withdrawals'/                                                                                                   &
              5X,'  Segment            [IWD] =',I7/                                                                                &
              5X,'  Outflow (m^3/sec)  [QWD] =',F7.2)
10680 FORMAT (5X,'  Layer              [KWD] =',12I7/                                                                              &
             (33X,12I7))
10690 FORMAT (:5X,'  Outflow (m^3/sec)  [QSW] =',12F7.2/                                                                           &
             (33X,12F7.2))
10700 FORMAT (/'1',A)
10710 FORMAT (3X,'Branch ',I0,' [CIN]'/                                                                                            &
             (5X,A,T25,'=',F9.3,1X,A))
10720 FORMAT (3X,'Tributary ',I0,' [CTR]'/                                                                                         &
             (5X,A,T25,'=',F9.3,1X,A))
10730 FORMAT (3X,'Distributed tributary ',I0,' [CDT]'/                                                                             &
             (5X,A,T25,'=',F9.3,1X,A))
10740 FORMAT (/'Surface calculations')
10750 FORMAT (3X,'Evaporation rate [EV]'/                                                                                          &
             (:3X,'  Branch ',I0,' = ',E10.3,' m^3/s'))                                                                                      ! SW 9/15/05 4/21/10
10755 format (3x,'Cumulative evaporation [VOLEV]'/                                                                                 &
             (:3X,'  Branch ',I0,' = ',F0.1,' m^3'))
10760 FORMAT (3X,'Precipitation [PR]'/                                                                                             &
             (3X,'  Branch ',I0,' = ',F8.6),' m/s')
10770 FORMAT (/1X,'External head boundary elevations'/)
10780 FORMAT (3X,'Branch ',I0/5X,'Upstream elevation   [ELUH] =',F8.3,' m')
10790 FORMAT (3X,'Branch ',I0/5X,'Downstream elevation [ELDH] =',F8.3,' m')
10800 FORMAT (/'Water Balance')
10810 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              3X,'  Spatial change  [VOLSR]  = ',E15.8,' m^3'/                                                                     &
              3X,'  Temporal change [VOLTR]  = ',E15.8,' m^3'/                                                                     &
              3X,'  Volume error             = ',E15.8,' m^3'/                                                                     &
              3X,'  Percent error            = ',E15.8,' %')
10820 FORMAT (3X,'Branch ',I0/                                                                                                     &
              3X,'  Spatial change  [VOLSBR] = ',E15.8,' m^3'/                                                                     &
              3X,'  Temporal change [VOLTBR] = ',E15.8,' m^3'/                                                                     &
              3X,'  Volume error             = ',E15.8,' m^3'/                                                                     &
              3X,'  Percent error            = ',E15.8,' %')
10830 FORMAT (/1X,'Energy Balance')
10840 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              3X,'  Spatially integrated energy   [ESR] = ',E15.8,' kJ'/                                                           &
              3X,'  Temporally integrated energy  [ETR] = ',E15.8,' kJ'/                                                           &
              3X,'  Energy error                        = ',E15.8,' kJ'/                                                           &
              3X,'  Percent error                       = ',E15.8,' %')
10850 FORMAT (3X,'  Spatially integrated energy  [ESBR] = ',E15.8,' kJ'/                                                           &
              3X,'  Temporally integrated energy [ETBR] = ',E15.8,' kJ'/                                                           &
              3X,'  Energy error                        = ',E15.8,' kJ'/                                                           &
              3X,'  Percent error                       = ',E15.8,' %')
10860 FORMAT (/1X,'Mass Balance')
10870 FORMAT (3X,'Branch ',I0)
10880 FORMAT (5X,A/                                                                                                                &
              5X,'  Spatially integrated mass  [CMBRS] = ',E15.8,1X,A/                                                             &
              5X,'  Temporally integrated mass [CMBRT] = ',E15.8,1X,A/                                                             &
              5X,'  Mass error                         = ',E15.8,1X,A/                                                             &
              5X,'  Percent error                      = ',E15.8,' %')
10890 FORMAT (/1X,A/                                                                                                               &
              3X,'Surface layer [KT] = ',I0/                                                                                       &
              3X,'Elevation   [ELKT] =',F10.3,' m')
10900 FORMAT (/3X,'Current upstream segment [CUS]'/                                                                                &
             (3X,'  Branch ',I0,' = ',I0))

             return
             end subroutine outputa
