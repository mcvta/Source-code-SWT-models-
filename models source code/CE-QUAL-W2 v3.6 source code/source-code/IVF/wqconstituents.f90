subroutine wqconstituents

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc
  EXTERNAL RESTART_OUTPUT


      if(macrophyte_on.and.update_kinetics)call porosity !v3.5
      DO JW=1,NWB
        KT = KTWB(JW)
        DO JB=BS(JW),BE(JW)
          IU = CUS(JB)
          ID = DS(JB)

!******** Kinetic sources/sinks

! v3.5 start
          IF (SEDIMENT_CALC(JW))then
            CALL SEDIMENT
            CALL SEDIMENTp
            CALL SEDIMENTn
            CALL SEDIMENTc
            IF(DYNSEDK(JW) == '      ON')CALL SEDIMENT_decay_rate
          end if
          do m=1,nmc
            if (macrophyte_calc(jw,m))then
              call macrophyte(m)
            end if
          end do
! v3.5 end
          IF (UPDATE_KINETICS) THEN
            IF (UPDATE_RATES) THEN
              CALL TEMPERATURE_RATES
              CALL KINETIC_RATES
            END IF
            DO JAC=1,NAC
              JC = CN(JAC)
              IF (JC == NPO4)                    CALL PHOSPHORUS
              IF (JC == NNH4)                    CALL AMMONIUM
              IF (JC == NNO3)                    CALL NITRATE
              IF (JC == NDSI)                    CALL DISSOLVED_SILICA
              IF (JC == NPSI)                    CALL PARTICULATE_SILICA
              IF (JC == NFE)                     CALL IRON
              IF (JC == NLDOM)                   CALL LABILE_DOM
              IF (JC == NRDOM)                   CALL REFRACTORY_DOM
              IF (JC == NLPOM)                   CALL LABILE_POM
              IF (JC == NRPOM)                   CALL REFRACTORY_POM
              IF (JC == NDO)                     CALL DISSOLVED_OXYGEN
              IF (JC >= NGCS  .AND. JC <= NGCE)  CALL GENERIC_CONST(JC-NGCS+1)
              IF (JC >= NSSS  .AND. JC <= NSSE)  CALL SUSPENDED_SOLIDS(JC-NSSS+1)
              IF (JC >= NAS   .AND. JC <= NAE)THEN
                IF(ALG_CALC(JC-NAS+1))CALL ALGAE(JC-NAS+1)
              ENDIF
              IF (JC >= NBODS .AND. JC <= NBODE)THEN
                IF(BOD_CALC(JC-NBODS+1)) CALL BIOCHEMICAL_O2_DEMAND(JC-NBODS+1)
              ENDIF
              IF (JC >= NZOOS  .AND. JC <= NZOOE .AND.ZOOPLANKTON_CALC)CALL zooplankton  		
              IF (JC == NLDOMP)                CALL LABILE_DOM_p
              IF (JC == NRDOMP)                CALL refractory_DOM_p
              IF (JC == NLPOMP)                CALL labile_POM_p
              IF (JC == NRPOMP)                CALL refractory_POM_p
              IF (JC == NLDOMN)                CALL LABILE_DOM_n
              IF (JC == NRDOMN)                CALL refractory_DOM_n
              IF (JC == NLPOMN)                CALL labile_POM_n
              IF (JC == NRPOMN)                CALL refractory_POM_n
            END DO
! v3.5 end
            IF (PH_CALC(JW)) CALL INORGANIC_CARBON
            IF (PH_CALC(JW)) CALL PH_CO2
          END IF
          DO JE=1,NEP   ! sw 5/16/06
            IF (EPIPHYTON_CALC(JW,JE)) CALL EPIPHYTON(JE)
          END DO

!******** External sources/sinks

          DO JAC=1,NAC
            JC = CN(JAC)
            IF (TRIBUTARIES) THEN
              DO JT=1,JTT
                IF (JB == JBTR(JT)) THEN
                  I = ITR(JT)
                  IF (I < CUS(JB)) I = CUS(JB)
                  DO K=KTTR(JT),KBTR(JT)
                    IF (QTR(JT) < 0.0) THEN
                      CSSB(K,I,JC) = CSSB(K,I,JC)+C1(K,I,JC)*QTR(JT)*QTRF(K,JT)
                    ELSE
                      CSSB(K,I,JC) = CSSB(K,I,JC)+CTR(JC,JT)*QTR(JT)*QTRF(K,JT)
                    END IF
                  END DO
                END IF
              END DO
            END IF
            IF (DIST_TRIBS(JB)) THEN
              DO I=IU,ID
                IF (QDT(I) < 0.0) THEN
                  CSSB(KT,I,JC) = CSSB(KT,I,JC)+C1(KT,I,JC)*QDT(I)
                ELSE
                  CSSB(KT,I,JC) = CSSB(KT,I,JC)+CDTR(JC,JB)*QDT(I)
                END IF
              END DO
            END IF
            IF (WITHDRAWALS) THEN
              DO JWD=1,JWW
                IF (QWD(JWD) /= 0.0) THEN
                  IF (JB == JBWD(JWD)) THEN
                    I = MAX(CUS(JBWD(JWD)),IWD(JWD))
                    FORALL(K=KTW(JWD):KBW(JWD))
                      CSSB(K,I,JC) = CSSB(K,I,JC)-C1S(K,I,JC)*QSW(K,JWD)
                    END FORALL
                  END IF
                END IF
              END DO
            END IF
            IF (PRECIPITATION(JW)) THEN
              FORALL (I=IU:ID)
                CSSB(KT,I,JC) = CSSB(KT,I,JC)+CPR(JC,JB)*QPR(I)
              END FORALL
            END IF
            IF (UP_FLOW(JB)) THEN
              DO K=KT,KB(IU)
                IF (.NOT. HEAD_FLOW(JB)) THEN
                  CSSB(K,IU,JC) = CSSB(K,IU,JC)+QINF(K,JB)*QIN(JB)*CIN(JC,JB)
                ELSE
                  IF (U(K,IU-1) >= 0.0) THEN
                    CSSB(K,IU,JC) = CSSB(K,IU,JC)+U(K,IU-1)*BHR1(K,IU-1)*C1S(K,IU-1,JC)
                  ELSE
                    CSSB(K,IU,JC) = CSSB(K,IU,JC)+U(K,IU-1)*BHR1(K,IU-1)*C1S(K,IU,JC)
                  END IF
                END IF
              END DO
            END IF
            IF (DN_FLOW(JB)) CSSB(KT:KB(ID),ID,JC) = CSSB(KT:KB(ID),ID,JC)-QOUT(KT:KB(ID),JB)*C1S(KT:KB(ID),ID,JC)
            IF (UP_HEAD(JB)) THEN
              DO K=KT,KB(IU)
                IUT = IU
                IF (QUH1(K,JB) >= 0.0) IUT = IU-1
                CSSUH1(K,JC,JB) = C1S(K,IUT,JC)*QUH1(K,JB)
                CSSB(K,IU,JC)   = CSSB(K,IU,JC)+CSSUH1(K,JC,JB)
              END DO
              IF (UH_INTERNAL(JB)) THEN
                IF (UHS(JB) /= DS(JBUH(JB)) .OR. DHS(JBUH(JB)) /= US(JB)) THEN
                  IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
                    I = UHS(JB)
                    DO K=KT,KB(IU)
                      CSSB(K,I,JC) = CSSB(K,I,JC)-CSSUH2(K,JC,JB)/DLT
                    END DO
                  ELSE
                    CALL UPSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
                  END IF
                END IF
              END IF
            END IF
            IF (DN_HEAD(JB)) THEN
              DO K=KT,KB(ID+1)
                IDT = ID+1
                IF (QDH1(K,JB) >= 0.0) IDT = ID
                CSSDH1(K,JC,JB) = C1S(K,IDT,JC)*QDH1(K,JB)
                CSSB(K,ID,JC)   = CSSB(K,ID,JC)-CSSDH1(K,JC,JB)
              END DO
              IF (DH_INTERNAL(JB)) THEN
                IF (DHS(JB) /= US(JBDH(JB)) .OR. UHS(JBDH(JB)) /= DS(JB)) THEN
                  IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
                    I = DHS(JB)
                    DO K=KT,KB(ID+1)
                      CSSB(K,I,JC) = CSSB(K,I,JC)+CSSDH2(K,JC,JB)/DLT
                    END DO
                  ELSE
                    CALL DOWNSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
                  END IF
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO

!**** Kinetic fluxes

      DO JW=1,NWB
        IF (FLUX(JW)) CALL KINETIC_FLUXES
      END DO

!**** Constituent transport

      DO JW=1,NWB
        KT = KTWB(JW)
        DO JB=BS(JW),BE(JW)
          IU = CUS(JB)
          ID = DS(JB)
          DO JAC=1,NAC
            JC   =  CN(JAC)
            COLD => C1S(:,:,JC)
!!$OMP PARALLEL SECTIONS
!!$OMP SECTION
            CALL HORIZONTAL_MULTIPLIERS
!!$OMP SECTION
            CALL VERTICAL_MULTIPLIERS
!!$OMP END PARALLEL SECTIONS
            CNEW => C1(:,:,JC)
            SSB  => CSSB(:,:,JC)
            SSK  => CSSK(:,:,JC)
            CALL HORIZONTAL_TRANSPORT
            DT = CNEW
            DO I=IU,ID
              CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,CNEW(:,I))
            END DO
          END DO
        END DO
      END DO
      IF (DERIVED_CALC) CALL DERIVED_CONSTITUENTS

    end subroutine wqconstituents
