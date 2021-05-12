!***********************************************************************************************************************************
!**                                              Task 1.4.4: Initial conditions                                                   **
!***********************************************************************************************************************************
Subroutine initcond
USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  !v3.5
  EXTERNAL RESTART_OUTPUT

  real    :: tmac,xsar


   BIC=B

  DO JW=1,NWB
    KT = KTWB(JW)
    IF (VERT_PROFILE(JW)) THEN

!**** Temperature and water quality

      OPEN (VPR(JW),FILE=VPRFN(JW),STATUS='OLD')
      READ (VPR(JW),*)
      IF (VERT_TEMP(JW)) READ (VPR(JW),'(//(8X,9F8.0))') (TVP(K,JW),K=KT,KBMAX(JW))
      IF (CONSTITUENTS) THEN
        DO JC=1,NCT
          IF (VERT_CONC(JC,JW))      READ (VPR(JW),'(//(8X,9F8.0))') (CVP(K,JC,JW),  K=KT,KBMAX(JW))
        END DO
        DO JE=1,NEP
          IF (VERT_EPIPHYTON(JW,JE)) READ (VPR(JW),'(//(8X,9F8.0))') (EPIVP(K,JW,JE),K=KT,KBMAX(JW))
        END DO
        IF (VERT_SEDIMENT(JW))       READ (VPR(JW),'(//(8X,9F8.0))') (SEDVP(K,JW),   K=KT,KBMAX(JW))
      END IF
    END IF

!** Longitudinal/vertical initial profiles

    IF (LONG_PROFILE(JW)) THEN
      OPEN (LPR(JW),FILE=LPRFN(JW),STATUS='OLD')
      READ (LPR(JW),*)
    END IF

!** Branch related variables

    IF (.NOT. RESTART_IN) THEN
      DO JB=BS(JW),BE(JW)

!****** Temperature

        DO I=CUS(JB),DS(JB)
          IF (LONG_TEMP(JW)) READ (LPR(JW),'(//(8X,9F8.0))') (T1(K,I),K=KT,KB(I))
          DO K=KT,KB(I)
            IF (ISO_TEMP(JW))  T1(K,I) = T2I(JW)
            IF (VERT_TEMP(JW)) T1(K,I) = TVP(K,JW)
            T2(K,I) = T1(K,I)
          END DO
        END DO
      END DO

!**** Constituents

      DO JC=1,NAC
        DO JB=BS(JW),BE(JW)
          DO I=CUS(JB),DS(JB)
            JAC = CN(JC)
            IF (LONG_CONC(JAC,JW)) READ (LPR(JW),'(//(8X,9F8.0))') (C2(K,I,JAC),K=KT,KB(I))
            DO K=KT,KB(I)
              IF (ISO_CONC(JAC,JW))  C2(K,I,JAC) = C2I(JAC,JW)
              IF (VERT_CONC(JAC,JW)) C2(K,I,JAC) = CVP(K,JAC,JW)
              C1(K,I,JAC)  = C2(K,I,JAC)
              C1S(K,I,JAC) = C1(K,I,JAC)
            END DO
          END DO
        END DO
      END DO

!**** Epiphyton

      DO JB=BS(JW),BE(JW)
        DO JE=1,NEP
          IF (EPIPHYTON_CALC(JW,JE)) THEN
            DO I=CUS(JB),DS(JB)
              IF (LONG_EPIPHYTON(JW,JE)) READ (LPR(JW),'(//(8X,9F8.0))') (EPD(K,I,JE),K=KT,KB(I))
              IF (ISO_EPIPHYTON(JW,JE))  EPD(:,I,JE) = EPICI(JW,JE)
              IF (VERT_EPIPHYTON(JW,JE)) EPD(:,I,JE) = EPIVP(:,JW,JE)                              ! CB 5/16/2009
            END DO
          END IF
        END DO
      END DO

!**** Sediments

      DO JB=BS(JW),BE(JW)
      SDKV(:,US(JB):DS(JB))=SDK(JW)
        IF (SEDIMENT_CALC(JW)) THEN
          DO I=CUS(JB),DS(JB)
            IF (LONG_SEDIMENT(JW)) READ (LPR(JW),'(//(8X,9F8.0))') (SED(K,I),K=KT,KB(I))
            DO K=KT,KB(I)
              IF (ISO_SEDIMENT(JW))  SED(K,I) = SEDCI(JW)
              IF (VERT_SEDIMENT(JW)) SED(K,I) = SEDVP(K,JW)
            END DO
            SED(KT,I)         = SED(KT,I)/H2(KT,I)
            SED(KT+1:KB(I),I) = SED(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
          END DO
        END IF
      END DO
        DO JB=BS(JW),BE(JW)
          IF (SEDIMENT_CALC(JW)) THEN
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                IF (ISO_SEDIMENT(JW))SEDp(K,I) = orgp(JW)*sedci(jw)
                IF (VERT_SEDIMENT(JW))SEDp(K,I) = SEDVP(K,JW)*orgp(jw)
                IF (LONG_SEDIMENT(JW)) sedp(k,i)=orgp(jw)*sed(k,i)
              END DO
              SEDp(KT,I)         = SEDp(KT,I)/H2(KT,I)
              SEDp(KT+1:KB(I),I) = SEDp(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
            END DO
          END IF
        END DO
        DO JB=BS(JW),BE(JW)
          IF (SEDIMENT_CALC(JW)) THEN
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                IF (ISO_SEDIMENT(JW))SEDn(K,I) = orgn(JW)*sedci(jw)
                IF (VERT_SEDIMENT(JW))SEDn(K,I) = SEDVP(K,JW)*orgn(jw)
                IF (LONG_SEDIMENT(JW)) sedn(k,i)=orgn(jw)*sed(k,i)
              END DO
              SEDn(KT,I)         = SEDn(KT,I)/H2(KT,I)
              SEDn(KT+1:KB(I),I) = SEDn(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
          END DO
        END IF
      END DO
        DO JB=BS(JW),BE(JW)
          IF (SEDIMENT_CALC(JW)) THEN
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                IF (ISO_SEDIMENT(JW))SEDc(K,I) = SEDCI(JW)*orgc(jw)
                IF (VERT_SEDIMENT(JW))SEDc(K,I) = SEDVP(K,JW)*orgc(jw)
                IF (LONG_SEDIMENT(JW)) sedc(k,i)=orgc(jw)*sed(k,i)
              END DO
              SEDc(KT,I)         = SEDc(KT,I)/H2(KT,I)
              SEDc(KT+1:KB(I),I) = SEDc(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
            END DO
          END IF
        END DO
! v3.5 end

      SED(:,US(BS(JW)):DS(BE(JW))) = SED(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      SEDp(:,US(BS(JW)):DS(BE(JW))) = SEDp(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      SEDn(:,US(BS(JW)):DS(BE(JW))) = SEDn(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      SEDc(:,US(BS(JW)):DS(BE(JW))) = SEDc(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)

      DO JB=BS(JW),BE(JW)
        do m=1,nmc
          IF (macrophyte_CALC(jw,m)) THEN

!c distributing initial macrophyte conc to bottom column cells; macwbci = g/m^3
            DO I=cus(jb),ds(jb)

              depkti=ELWS(i)-el(kti(i)+1,i)

              if(depkti.ge.thrkti)then
                kticol(i)=.true.
                jt=kti(i)
              else
                kticol(i)=.false.
                jt=kti(i)+1
              end if

              je=kb(i)
              DO j=jt,je
                if(j.le.kt)then
                  k=kt
                else
                  k=j
                end if
                macrc(j,K,I,m) = macwbci(jw,m)
                smacrc(j,K,I,m) = macwbci(jw,m)
              END DO
            END DO

            DO I=cus(jb),ds(jb)
              tmac=0.0
              xsar=0.0
              do k=kti(i),kt
                jt=k
                je=kb(i)
                colb=el(k+1,i)
                coldep=ELWS(i)-colb
                do j=jt,je
                  tmac=tmac+macrc(j,kt,i,m)*cw(j,i)*coldep
                  xsar=xsar+cw(j,i)*coldep
                end do
              end do
              mac(kt,i,m)=tmac/xsar
              smac(kt,i,m)=mac(kt,i,m)

              DO K=KT+1,KB(I)
                jt=k
                je=kb(i)
                tmac=0.0
                do j=jt,je
                  tmac=tmac+macrc(j,k,i,m)*cw(j,i)
                end do
                mac(k,i,m)=tmac/b(k,i)
                smac(k,i,m)=mac(k,i,m)
              end do
            end do

            do i=cus(jb),ds(jb)
              jt=kti(i)
              je=kb(i)
              do j=jt,je
                if(j.lt.kt)then
                  colb=el(j+1,i)
                else
                  colb=el(kt+1,i)
                end if
                coldep=ELWS(i)-colb
                macrm(j,kt,i,m)=macrc(j,kt,i,m)*coldep*cw(j,i)*dlx(i)
                smacrm(j,kt,i,m)=macrm(j,kt,i,m)
              end do

              do K=KT+1,KB(I)

                jt=k
                je=kb(i)

                do j=jt,je

                  macrm(j,k,i,m)=macrc(j,k,i,m)*h2(k,i)*cw(j,i)*dlx(i)
                  smacrm(j,k,i,m)=macrm(j,k,i,m)
                end do

              END DO
            end do

          END IF
        end do
      end do
! v3.5 end

!**** Energy

      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB),DS(JB)
          IF (ENERGY_BALANCE(JW)) THEN
            DO K=KT,KB(I)
              EBRI(JB) = EBRI(JB)+T2(K,I)*DLX(I)*BH2(K,I)
            END DO
          END IF
          DO K=KT,KB(I)
            CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)+C2(K,I,CN(1:NAC))*DLX(I)*BH2(K,I)
          END DO
        END DO

! v3.5 start
!c   initializing macrophyte temporal mass balance term....
          do m=1,nmc
            if(macrophyte_calc(jw,m))then
              DO I=CUS(JB),DS(JB)
                if(kticol(i))then
                  jt=kti(i)
                else
                  jt=kti(i)+1
                end if
                je=kb(i)
                do j=jt,je
                  maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(j,Kt,I,m)
                end do
                DO K=KT+1,KB(I)
                  jt=k
                  je=kb(i)
                  do j=jt,je
                    maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(j,K,I,m)
                  end do
                END DO
              END DO
            end if
          end do
! v3.5 end

!****** Ice cover

        IF (ICE_CALC(JW)) THEN
          ICETH(US(JB):DS(JB)) = ICETHI(JW)
          ICE(US(JB):DS(JB))   = ICETH(US(JB):DS(JB)) > 0.0
        END IF

!****** Vertical eddy viscosity

        IUT = CUS(JB)
        IDT = DS(JB)-1
        IF (UP_HEAD(JB)) IUT = IU-1
        IF (DN_HEAD(JB)) IDT = ID
        DO I=IUT,IDT
          DO K=KT,KB(I)-1
            AZ(K,I)    = AZMIN
            TKE(K,I,1) = 1.25E-7
            TKE(K,I,2) = 1.0E-9
          END DO
        END DO
        DO JWR=1,NIW
          IF (WEIR_CALC) AZ(MAX(KT,KTWR(JWR)-1):KBWR(JWR),IWR(JWR)) = 0.0
        END DO
      END DO
    END IF

!** Horizontal diffusivities

    DO JB=BS(JW),BE(JW)
      DO I=CUS(JB),DS(JB)-1
        DO K=KT,KBMIN(I)
          DX(K,I) = DXI(JW)
          IF (INTERNAL_WEIR(K,I)) DX(K,I) = 0.0
        END DO
      END DO
    END DO
    IF (VERT_PROFILE(JW)) CLOSE (VPR(JW))
    IF (LONG_PROFILE(JW)) CLOSE (LPR(JW))
  END DO

! Atmospheric pressure

  IF (CONSTITUENTS)PALT(:) = (1.0-ELWS(:)/1000.0/44.3)**5.25       ! SW 2/3/08


  RETURN

  End subroutine initcond
