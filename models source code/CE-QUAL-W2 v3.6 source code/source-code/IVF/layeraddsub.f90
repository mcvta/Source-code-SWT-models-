subroutine layeraddsub
USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc


INTEGER :: KTMAX,JBIZ,KKB
REAL    :: DUMMY,TMAC
REAL(R8):: W1,W2,W3

!***********************************************************************************************************************************
!**                                       Task 2.5: Layer - Segment Additions and Subtractions                                    **
!***********************************************************************************************************************************

!** Water surface minimum thickness

    DO JW=1,NWB
      KT       =  KTWB(JW)
      ZMIN(JW) = -1000.0
      KTMAX    =  2                                                                                                 ! SR 10/17/05
      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB),DS(JB)
          IF(KB(I) > KTMAX) KTMAX = KB(I)                                                                           ! SR 10/17/05
          IF (Z(I) > ZMIN(JW)) THEN
            IZMIN(JW) = I
            JBIZ      = JB
          END IF
          ZMIN(JW) = MAX(ZMIN(JW),Z(I))
        END DO
      END DO
      ADD_LAYER = ZMIN(JW) < -0.85*H(KT-1,JW) .AND. KT /= 2
      SUB_LAYER = ZMIN(JW) >  0.60*H(KT,JW)   .AND. KT < KTMAX                                                       ! SR 10/17/05
      IF (KTWB(JW) == KMX-1 .AND. SLOPE(JBIZ) > 0.0 .AND. SUB_LAYER .AND. ONE_LAYER(IZMIN(JW))) THEN
        IF (ZMIN(JW) > 0.99*H(KT,JW)) WRITE (WRN,'(A,I0,2(A,F0.3))') 'Low water in segment ',IZMIN(JW),&
                                              ' water surface deviation'//' = ',ZMIN(JW),' at day ',JDAY
        WARNING_OPEN = .TRUE.
        SUB_LAYER    = .FALSE.
      END IF

      IF(ADD_LAYER == .TRUE. .OR. SUB_LAYER == .TRUE.)THEN
      LAYERCHANGE(JW)=.TRUE.
      ELSE
      LAYERCHANGE(JW)=.FALSE.
      ENDIF

!**** Add layers

      DO WHILE (ADD_LAYER)
        IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') '   Add layer ',KT-1,&
                                                        ' at Julian day = ',JDAY,'   NIT = ',NIT,' IZMIN =',IZMIN(JW)   ! SW 1/23/06

!****** Variable initialization

        KTWB(JW) = KTWB(JW)-1
        KT       = KTWB(JW)
        ilayer = 0
        DO JB=BS(JW),BE(JW)
          IU = CUS(JB)
          ID = DS(JB)
          DO I=IU-1,ID+1
            Z(I)          =  H(KT,JW)+Z(I)
            H1(KT,I)      =  H(KT,JW)-Z(I)
            H1(KT+1,I)    =  H(KT+1,JW)
            H2(KT+1,I)    =  H(KT+1,JW)
            AVH1(KT,I)    = (H1(KT,I)  +H1(KT+1,I))*0.5
            AVH1(KT+1,I)  = (H1(KT+1,I)+H1(KT+2,I))*0.5
            IF (.NOT. TRAPEZOIDAL(JW)) THEN
              BH1(KT,I)   = BH1(KT+1,I)-Bnew(KT+1,I)*H1(KT+1,I)                              ! SW 1/23/06
              BH1(KT+1,I) = Bnew(KT+1,I)*H1(KT+1,I)                                          ! SW 1/23/06
            ELSE
              CALL GRID_AREA1(EL(KT,I)-Z(I),EL(KT+1,I),BH1(KT,I),DUMMY)                                                          !SW 08/03/04
              BH1(KT+1,I) = 0.25*H1(KT+1,JW)*(BB(KT,I)+2.*B(KT+1,I)+BB(KT+1,I))
            ENDIF
            VOL(KT,I)     = BH1(KT,I)  *DLX(I)
            VOL(KT+1,I)   = BH1(KT+1,I)*DLX(I)
            BKT(I)        = BH1(KT,I)/H1(KT,I)
            DEPTHB(KT,I)  = H1(KT,I)
            DEPTHM(KT,I)  = H1(KT,I)*0.5
            BI(KT:KB(I),I) =  B(KT:KB(I),I)   ! SW 8/26/05
            BI(KT,I)       =  B(KTI(I),I)
            T1(KT,I)      = T1(KT+1,I)
            sdkv(kt,i)    = sdkv(kt+1,i)      ! SW 1/18/08
            sed(kt,i)     = sed(kt+1,i)       ! SW 1/18/08
            sedn(kt,i)    = sedn(kt+1,i)      ! SW 1/18/08
            sedp(kt,i)    = sedp(kt+1,i)      ! SW 1/18/08
            sedc(kt,i)    = sedc(kt+1,i)      ! SW 1/18/08
!            RHO(KT,I)     = DENSITY(T1(KT,I),MAX(TDS(KT,I),0.0),MAX(TISS(KT,I),0.0))    ! SR 5/15/06
            DO K=KT+1,KMX
              DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
              DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
            END DO
            C1(KT,I,CN(1:NAC))             = C1(KT+1,I,CN(1:NAC))
            CSSK(KT,I,CN(1:NAC))           = CSSK(KT+1,I,CN(1:NAC))
            KF(KT,I,KFCN(1:NAF(JW),JW))    = KF(KT+1,I,KFCN(1:NAF(JW),JW))
            KFS(KT,I,KFCN(1:NAF(JW),JW))   = KF(KT+1,I,KFCN(1:NAF(JW),JW))
            KF(KT+1,I,KFCN(1:NAF(JW),JW))  = 0.0
            KFS(KT+1,I,KFCN(1:NAF(JW),JW)) = 0.0
            if(kt >= kbi(i))then         ! CB 5/24/06
            adx(kt+1,i)=0.0             ! CB 5/15/06
            c1(kt+1,i,cn(1:nac))=0.0    ! CB 5/15/06
            cssk(kt+1,i,cn(1:nac))=0.0  ! CB 5/15/06
            endif                       ! CB 5/15/06
            DO JE=1,NEP
              if(kt < kbi(i))then    ! CB 4/28/06
              EPD(KT,I,JE)   = EPD(KT+1,I,JE)
              EPM(KT,I,JE)   = EPD(KT,I,JE)*(Bi(KT,I)-B(KT+1,I)+2.0*H1(KT,I))*DLX(I)             ! SR 5/15/06
              EPM(KT+1,I,JE) = EPM(KT+1,I,JE)-EPM(KT,I,JE)
              EPC(KT,I,JE)   = EPM(KT,I,JE)/VOL(KT,I)
              EPC(KT+1,I,JE) = EPM(KT+1,I,JE)/VOL(KT+1,I)
              else
              EPD(KT,I,JE)   = EPD(KT+1,I,JE) ! SW 5/15/06
              EPM(KT,I,JE)   = EPM(KT+1,I,JE)
              EPC(KT,I,JE) =   EPC(KT+1,I,JE)
              EPD(KT+1,I,JE) = 0.0
              EPM(KT+1,I,JE) = 0.0
              EPC(KT+1,I,JE) = 0.0
              end if                ! CB 4/28/06
            END DO
          END DO
! v3.5 start
!********macrophytes...
          do i=iu,id
            jt=kti(i)
            je=kb(i)
            do j=jt,je
              if(j.lt.kt)then
                colb=el(j+1,i)
              else
                colb=el(kt+1,i)
              end if
              coldep=ELWS(i)-colb
              mact(j,kt,i)=mact(j,kt+1,i)
              do m=1,nmc
                if(macrophyte_calc(jw,m))then
                  macrc(j,kt,i,m)=macrc(j,kt+1,i,m)
                  macrm(j,kt,i,m)=macrc(j,kt,i,m)*cw(j,i)*coldep*dlx(i)
                end if
              end do
            end do

            jt=kt+1
            je=kb(i)
            do j=jt,je
              do m=1,nmc
                if(macrophyte_calc(jw,m))then
                  macrm(j,kt+1,i,m)=macrc(j,kt+1,i,m)*cw(j,i)*h(kt+1,jw)*dlx(i)
                end if
              end do
            end do
          end do
! v3.5 end
          DO I=IU-1,ID
            AVHR(KT+1,I) = H1(KT+1,I) +(H1(KT+1,I+1) -H1(KT+1,I)) /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            AVHR(KT,I)   = H1(KT,I)   +(H1(KT,I+1)   -H1(KT,I))   /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT,I)   = BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT+1,I) = BH1(KT+1,I)+(BH1(KT+1,I+1)-BH1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            U(KT,I)      = U(KT+1,I)
          END DO
          DO I=IU,ID
            IF (ONE_LAYER(I)) THEN
              W(KT,I) = 0.0
            ELSE
              W1      =  W(KT+1,I)*BB(KT+1,I)
              W2      = (BHR1(KT+1,I)*U(KT+1,I)-BHR1(KT+1,I-1)*U(KT+1,I-1))/DLX(I)
              W3      = (-QSS(KT+1,I)*BH1(KT+1,I)/(BH1(KT+1,I)+BH1(KT,I)))/DLX(I)
              W(KT,I) = (W1+W2+W3)/BB(KT,I)
            END IF
          END DO
          IF (UP_HEAD(JB)) THEN
            BHSUM                     = BHR1(KT,IU-1)            +BHR1(KT+1,IU-1)
            QUH1(KT,JB)               = QUH1(KT+1,JB)            *BHR1(KT,IU-1)  /BHSUM
            QUH1(KT+1,JB)             = QUH1(KT+1,JB)            *BHR1(KT+1,IU-1)/BHSUM
            TSSUH1(KT,JB)             = TSSUH1(KT+1,JB)          *BHR1(KT,IU-1)  /BHSUM
            TSSUH1(KT+1,JB)           = TSSUH1(KT+1,JB)          *BHR1(KT+1,IU-1)/BHSUM
            CSSUH1(KT,CN(1:NAC),JB)   = CSSUH1(KT+1,CN(1:NAC),JB)*BHR1(KT,IU-1)  /BHSUM
            CSSUH1(KT+1,CN(1:NAC),JB) = CSSUH1(KT+1,CN(1:NAC),JB)*BHR1(KT+1,IU-1)/BHSUM
          END IF
          IF (DN_HEAD(JB)) THEN
            BHSUM                     = BHR1(KT,ID)              +BHR1(KT+1,ID)
            QDH1(KT,JB)               = QDH1(KT+1,JB)            *BHR1(KT,ID)    /BHSUM
            QDH1(KT+1,JB)             = QDH1(KT+1,JB)            *BHR1(KT+1,ID)  /BHSUM
            TSSDH1(KT,JB)             = TSSDH1(KT+1,JB)          *BHR1(KT,ID)    /BHSUM
            TSSDH1(KT+1,JB)           = TSSDH1(KT+1,JB)          *BHR1(KT+1,ID)  /BHSUM
            CSSDH1(KT,CN(1:NAC),JB)   = CSSDH1(KT+1,CN(1:NAC),JB)*BHR1(KT,ID)    /BHSUM
            CSSDH1(KT+1,CN(1:NAC),JB) = CSSDH1(KT+1,CN(1:NAC),JB)*BHR1(KT+1,ID)  /BHSUM
          END IF
          DO I=IU,ID-1
            DX(KT,I) = DXI(JW)
            IF (INTERNAL_WEIR(KT,I)) DX(KT,I) = 0.0
          END DO
          IUT = IU
          IDT = ID-1
          IF (UP_HEAD(JB)) IUT = IU-1
          IF (DN_HEAD(JB)) IDT = ID
          DO I=IUT,IDT
            AZ(KT,I)  = AZMIN
            TKE(KT,I,1) = 1.25E-7 !sg 10/4/07
            TKE(KT,I,2) = 1.0E-9  !sg 10/4/07
            SAZ(KT,I) = AZMIN
            IF (INTERNAL_WEIR(KT,I)) THEN
              AZ(KT,I)  = 0.0
              TKE(KT,I,1) = 0.0 !sg  10/4/07
              TKE(KT,I,2) = 0.0 !sg  10/4/07
              SAZ(KT,I) = 0.0
            END IF
          END DO
          IF (CONSTITUENTS) THEN
            CALL TEMPERATURE_RATES
            CALL KINETIC_RATES
          END IF

!******** Upstream active segment

          IUT = US(JB)
          IF (SLOPE(JB) == 0.0) THEN
            DO I=US(JB),DS(JB)
              IF (KB(I)-KT < NL(JB)-1) IUT = I+1
            END DO
          ELSE
            DO I=US(JB)-1,DS(JB)+1
              IF (KB(I) > KBI(I)) THEN
                Bnew(KB(I),I)  = b(kb(i),i)                                                    ! SW 1/23/06                                           ! SW 3/2/05
                DX(KB(I),I) = 0.0
                KB(I)       = KB(I)-1
                ilayer(i) = 1
                u(kb(i)+1,i)=0.0                                                               ! SW 1/23/06
                WRITE (WRN,'(2(A,I8),A,F0.3)') 'Raising bottom layer at segment ',I,' at iteration ',NIT,' at Julian day ',JDAY
              END IF
            END DO                   ! SW 1/23/06
            DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
!              IF (KB(I)-KT < NL(JB)-1) IUT = I+1    ! SW 1/23/06
!                IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))                    ! SW 1/23/06                                 ! SW 3/2/05
                IF (I /= US(JB)-1) KBMIN(I-1) = MIN(KB(I-1),KB(I))                     ! SW 1/23/06
                if(kbi(i) < kb(i))then
                bkt(i)=bh1(kt,i)/(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
                depthb(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
                depthm(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))*0.5    ! SW 1/23/06
                avhr(kt,i)=(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))   +(H1(KT,I+1)-(el(kbi(i)+1,i+1)-el(kb(i)+1,i+1))&
                 -H1(KT,I)+(el(kbi(i)+1,i)-el(kb(i)+1,i)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                                                           ! SW 1/23/06
                end if
            ENDDO
        DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
         do k=kbmin(i)+1,kb(i)
         u(k,i)=0.0
         end do
        END DO    ! SW 11/9/07

        DO I=US(JB),DS(JB)     ! SW 11/9/07
         if(ilayer(i).eq.1.and.ilayer(i+1).eq.0)then
         bhrsum=0.0
         q(i)=0.0
          DO K=KT,KBMIN(I)
            IF (.NOT. INTERNAL_WEIR(K,I)) THEN
              BHRSUM = BHRSUM+BHR1(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR1(K,I)
            END IF
          END DO
          DO K=KT,KBMIN(I)
            IF (INTERNAL_WEIR(K,I)) THEN
              U(K,I) = 0.0
            ELSE
              U(K,I) =  U(K,I)+(QC(I)-Q(I))/BHRSUM
            END IF
          END DO
          elseif(ilayer(i).eq.1.and.ilayer(i-1).eq.0)then
          bhrsum=0.0
          q(i-1)=0.0
          DO K=KT,KBMIN(I-1)
            IF (.NOT. INTERNAL_WEIR(K,I-1)) THEN
              BHRSUM = BHRSUM+BHR1(K,I-1)
              Q(I-1)   = Q(I-1)+U(K,I-1)*BHR1(K,I-1)
            END IF
          END DO
          DO K=KT,KBMIN(I-1)
            IF (INTERNAL_WEIR(K,I-1)) THEN
              U(K,I-1) = 0.0
            ELSE
              U(K,I-1) =  U(K,I-1)+(QC(I-1)-Q(I-1))/BHRSUM
            END IF
          END DO
          endif

            END DO
          END IF

!******** Segment addition

          IF (IUT /= IU) THEN
            IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/17X,2(A,I0))') ' Add segments ',IUT,' through ',IU-1
            DO I=IUT-1,IU-1
              Z(I)         =  Z(IU)
              KTI(I)       =  KTI(IU)
              H1(KT+1,I)   =  H(KT+1,JW)
              AVH1(KT+1,I) = (H1(KT+1,I)+H1(KT+2,I))*0.5
              AVHR(KT+1,I) =  H1(KT+1,I)+(H1(KT+1,I+1)-H1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
              IF (.NOT. TRAPEZOIDAL(JW)) THEN
                BH1(KT+1,I) =  B(KT+1,I)*H(KT+1,JW)
                H1(KT,I)    =  H(KT,JW)-Z(I)
                BI(KT:KB(I),I) =  B(KT:KB(I),I)                                                                        ! SW 4/18/07
                BI(KT,I)    =  B(KTI(I),I)
                AVH1(KT,I)  = (H1(KT,I)+H1(KT+1,I))*0.5
                AVHR(KT,I)  =  H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                       !SW 07/29/04
                BH1(KT,I)   =  BI(KT,I)*(EL(KT,I)-Z(I)*COSA(JB)-EL(KTI(I)+1,I))/COSA(JB)
                IF (KTI(I) >= KB(I)) BH1(KT,I) = B(KT,I)*H1(KT,I)
                DO K=KTI(I)+1,KT
                  BH1(KT,I) = BH1(KT,I)+BH(K,I)
                END DO
              ELSE
                CALL GRID_AREA1 (EL(KT,I)-Z(I),EL(KT+1,I),BH1(KT,I),BI(KT,I))                                                    !SW 08/03/04
                BH1(KT+1,I) =  0.25*H(KT+1,JW)*(BB(KT,I)+2.*B(KT+1,I)+BB(KT+1,I))
                H1(KT,I)    =  H(KT,JW)-Z(I)
                AVH1(KT,I)  = (H1(KT,I)+H1(KT+1,I))*0.5
                AVHR(KT,I)  =  H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                       !SW 07/29/04
              END IF
              BKT(I)      = BH1(KT,I)/H1(KT,I)
              DEPTHB(KT,I) = H1(KT,I)
              DEPTHM(KT,I) = H1(KT,I)*0.5
              DO K=KT+1,KB(I)
                DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
                DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
              END DO
            END DO
            DO I=IUT-1,IU-1
              BHR1(KT+1,I) = BH1(KT+1,I)+(BH1(KT+1,I+1)-BH1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                !SW 07/29/04
              BHR1(KT,I)   = BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                !SW 07/29/04
            END DO
            DO I=IUT,IU-1
              ICE(I)   = ICE(IU)
              ICETH(I) = ICETH(IU)
              WIND2(I) = WIND2(IU)
              IF (DYNAMIC_SHADE(I)) CALL SHADING
              DO K=KT,KB(I)
                DX(K,I) = DXI(JW)
                U(K,I)  = U(K,IU)
                SDKV(K,I)=SDK(JW)     ! SW 1/18/08
                IF (INTERNAL_WEIR(K,I)) THEN
                  DX(K,I) = 0.0
                  U(K,I)  = 0.0
                END IF
                T1(K,I)           = T1(K,IU)
                T2(K,I)           = T1(K,IU)
                SU(K,I)           = U(K,IU)
                C1(K,I,CN(1:NAC)) = C1(K,IU,CN(1:NAC))
                C2(K,I,CN(1:NAC)) = C1(K,IU,CN(1:NAC))
                DO JE=1,NEP
                  EPD(K,I,JE) = 0.01
                  EPC(K,I,JE) = 0.01/H1(K,I)
                END DO
                CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)+C1(K,IU,CN(1:NAC))*DLX(I)*BH1(K,I)
                EBRI(JB)            = EBRI(JB)           +T1(K,IU)          *DLX(I)*BH1(K,I)
              END DO
              DO K=KT,KB(I)-1
                AZ(K,I)  = AZ(K,IU)
		TKE(K,I,1) = TKE(K,IU,1) !sg 10/4/07
                TKE(K,I,2) = TKE(K,IU,2) !sg 10/4/07
                SAZ(K,I) = AZ(K,IU)
                IF (INTERNAL_WEIR(K,I)) THEN
                  AZ(K,I)  = 0.0
	          TKE(K,I,1) = 0.0 !sg 10/4/07
                  TKE(K,I,2) = 0.0 !sg 10/4/07
                  SAZ(K,I) = 0.0
                END IF
              END DO
            END DO
!*********macrophytes
              do m=1,nmc
                if(macrophyte_calc(jw,m))then
                  jt=kti(i)
                  je=kb(i)
                  do j=jt,je
                    if(j.lt.kt)then
                      colb=el(j+1,i)
                    else
                      colb=el(kt+1,i)
                    end if
                    coldep=ELWS(i)-colb
                    macrc(j,kt,i,m)=macwbci(jw,m)
                    macrm(j,kt,i,m)=macrc(j,kt,i,m)*cw(j,i)*coldep*dlx(i)
                    maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(j,kt,i,m)
                  end do
                  DO K=KT+1,KB(I)
                    jt=k
                    je=kb(i)
                    do j=jt,je
                      macrc(j,k,i,m)=macwbci(jw,m)
                      macrm(j,k,i,m)=macrc(j,k,i,m)*cw(j,i)*h2(k,i)*dlx(i)
                      maCMBRT(JB,m) = maCMBRT(JB,m)+macrm(j,k,i,m)
                    end do
                  END DO
                end if
              end do
            U(KB(IUT):KB(IU),IU-1)  = 0.0
            SU(KB(IUT):KB(IU),IU-1) = 0.0
            ADX(KB(IUT):KB(IU),IU)  = 0.0
            IU                      = IUT
            CUS(JB)                 = IU
            IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)
            IF (UH_INTERNAL(JB)) THEN
              IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
                KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
              ELSE
                DO KKB=KT,KMX
                  IF (EL(KKB,IU) <= EL(KB(UHS(JB)),UHS(JB))) EXIT
                END DO
                KB(IU-1) = MIN(KKB,KB(IU))
              END IF
            END IF
            IF (UP_HEAD(JB)) THEN
              AZ(KT:KB(IU-1)-1,IU-1)  = AZMIN
              TKE(KT:KB(IU-1)-1,IU-1,1) = 1.25E-7 !sg 10/4/07
              TKE(KT:KB(IU-1)-1,IU-1,2) = 1.0E-9  !sg 10/4/07
              SAZ(KT:KB(IU-1)-1,IU-1) = AZMIN
            END IF
          END IF
          IF (CONSTITUENTS) THEN
            CALL TEMPERATURE_RATES
            CALL KINETIC_RATES
          END IF

!******** Total active cells and single layers

          DO I=IU,ID
            NTAC         = NTAC+1
            ONE_LAYER(I) = KTWB(JW) == KB(I)
          END DO
          NTACMX = MAX(NTAC,NTACMX)
        END DO
        CALL INTERPOLATION_MULTIPLIERS

!****** Additional layers

        ZMIN(JW) = -1000.0
        DO JB=BS(JW),BE(JW)
          DO I=CUS(JB),DS(JB)
            ZMIN(JW) = MAX(ZMIN(JW),Z(I))
          END DO
        END DO
        ADD_LAYER = ZMIN(JW) < -0.80*H(KT-1,JW) .AND. KT /= 2
      END DO

!**** Subtract layers

      DO WHILE (SUB_LAYER)
        IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,1x,13("*"))') 'Subtract layer ',&
        KT,' at Julian day = ', JDAY,' NIT = ',NIT,' IZMIN =',IZMIN(JW)                                                                    ! SW 1/23/06

!****** Variable initialization

        KTWB(JW) = KTWB(JW)+1
        KT       = KTWB(JW)
        ilayer=0   ! SW 1/23/06  11/7/07
        DO JB=BS(JW),BE(JW)
          IU = CUS(JB)
          ID = DS(JB)
          IF (CONSTITUENTS) DO1(KT-1,IU-1:ID+1) = 0.0
          DO I=IU-1,ID+1
            Z(I)         =  Z(I)-H(KT-1,JW)
            H1(KT-1,I)   =  H(KT-1,JW)
            H1(KT,I)     =  H(KT,JW)-Z(I)
            BI(KT,I)     =  B(KTI(I),I)
            BI(KT-1,I)   =  B(KT-1,I)
            AVH1(KT-1,I) = (H(KT-1,JW)+H(KT,JW))*0.5
            AVH1(KT,I)   = (H1(KT,I)+H1(KT+1,I))*0.5
            IF(.NOT. TRAPEZOIDAL(JW))THEN
              if(kb(i) >= kt)then                            ! sw 1/23/06
              BH1(KT,I)   =  BH1(KT-1,I)+BH1(KT,I)           ! sw 1/23/06
              BH1(KT-1,I) =  B(KT-1,I)*H(KT-1,JW)
              else  ! sw 1/23/06
              BH1(KT,I)   =  BH1(KT-1,I)         ! sw 1/23/06
              end if  ! sw 1/23/06
            ELSE
              CALL GRID_AREA1 (EL(KT,I)-Z(I),EL(KT+1,I),BH1(KT,I),DUMMY)                                                         !SW 08/03/04
              BH1(KT-1,I) = 0.25*H1(KT-1,JW)*(BB(KT-2,I)+2.*B(KT-1,I)+BB(KT-1,I))
            ENDIF
            VOL(KT,I)   =  BH1(KT,I)*DLX(I)
            VOL(KT-1,I) =  BH1(KT-1,I)*DLX(I)
            BKT(I)      =  BH1(KT,I)/H1(KT,I)
            if(kb(i) >= kt)then                             ! SW 1/23/06
              U(KT,I)  = (U(KT-1,I)*BHR1(KT-1,I)+U(KT,I)*BHR(KT,I))/(BHR1(KT-1,I)+BHR(KT,I))
              T1(KT,I) = (T1(KT-1,I)*(BH1(KT,I)-BH(KT,I))+T1(KT,I)*BH(KT,I))/BH1(KT,I)
            ELSE
!              EBRI(JB) = EBRI(JB)-T1(KT,I)*VOL(KT,I)   1/23/06
              u(kt,i) = u(kt-1,i)    ! SW 1/23/06
              t1(kt,i) = t1(kt-1,i)  ! SW 1/23/06
            END IF
            if(kb(i) >= kt)then       ! sw 1/23/06
            C1(KT,I,CN(1:NAC))             = (C1(KT-1,I,CN(1:NAC))  *(BH1(KT,I)-BH(KT,I))+C1(KT,I,CN(1:NAC))  *BH(KT,I))/BH1(KT,I)
            CSSK(KT,I,CN(1:NAC))           = (CSSK(KT-1,I,CN(1:NAC))*(BH1(KT,I)-BH(KT,I))+CSSK(KT,I,CN(1:NAC))*BH(KT,I))/BH1(KT,I)
            else
            C1(KT,I,CN(1:NAC))             = C1(KT-1,I,CN(1:NAC))    ! SW 1/23/06
            CSSK(KT,I,CN(1:NAC))           = CSSK(KT-1,I,CN(1:NAC))   ! sw 1/23/06
            endif   ! sw 1/23/06
            CSSB(KT,I,CN(1:NAC))           =  CSSB(KT-1,I,CN(1:NAC))+CSSB(KT,I,CN(1:NAC))
            KF(KT,I,KFCN(1:NAF(JW),JW))    =  KF(KT-1,I,KFCN(1:NAF(JW),JW))
            KFS(KT,I,KFCN(1:NAF(JW),JW))   =  KFS(KT-1,I,KFCN(1:NAF(JW),JW))
            C1(KT-1,I,CN(1:NAC))           =  0.0
            C2(KT-1,I,CN(1:NAC))           =  0.0
            CSSB(KT-1,I,CN(1:NAC))         =  0.0
            CSSK(KT-1,I,CN(1:NAC))         =  0.0
            KF(KT-1,I,KFCN(1:NAF(JW),JW))  =  0.0
            KFS(KT-1,I,KFCN(1:NAF(JW),JW)) =  0.0
            DO JE=1,NEP
              if(kt <= kbi(i))then    ! CB 4/28/06
              EPM(KT,I,JE)   = EPM(KT-1,I,JE)+EPM(KT,I,JE)
              EPD(KT,I,JE)   = EPM(KT,I,JE)/((BI(KT,I)-BI(KT+1,I)+2.0*H1(KT,I))*DLX(I))
              EPC(KT,I,JE)   = EPM(KT,I,JE)/VOL(KT,I)
              EPM(KT-1,I,JE) = 0.0
              EPD(KT-1,I,JE) = 0.0
              EPC(KT-1,I,JE) = 0.0
              else                   ! SW 5/15/06
              EPM(KT,I,JE)   = EPM(KT-1,I,JE)
              EPD(KT,I,JE)   = EPD(KT-1,I,JE)
              EPC(KT,I,JE)   = EPC(KT-1,I,JE)
              EPM(KT-1,I,JE) = 0.0
              EPD(KT-1,I,JE) = 0.0
              EPC(KT-1,I,JE) = 0.0
              endif                   ! CB 4/28/06
            END DO
          END DO
          do i=iu,id
            do m=1,nmc
              if(macrophyte_calc(jw,m))then
                if(kticol(i))then
                  jt=kti(i)
                else
                  jt=kti(i)+1
                end if
                je=kb(i)
                do j=jt,je
                  if(j.lt.kt)then
                    colb=el(j+1,i)
                  else
                    colb=el(kt+1,i)
                  end if
                  coldep=ELWS(i)-colb
                  if(j.lt.kt)then
                    macrm(j,kt,i,m)=macrm(j,kt-1,i,m)
                  else
                    macrm(j,kt,i,m)=macrm(j,kt-1,i,m)+macrm(j,kt,i,m)
                  end if
                  if(cw(j,i).gt.0.0)then
                    macrc(j,kt,i,m)=macrm(j,kt,i,m)/(cw(j,i)*coldep*dlx(i))
                  else
                    macrc(j,kt,i,m)=0.0
                  end if
                  macrm(j,kt-1,i,m)=0.0
                  macrc(j,kt-1,i,m)=0.0

                end do
              end if
            end do
            jt=kti(i)
            je=kb(i)
            do j=jt,je
              mact(j,kt,i)=0.0
              mact(j,kt-1,i)=0.0
            end do
            do m=1,nmc
              if(macrophyte_CALC(jw,m))then
                 do j=jt,je
                   mact(j,kt,i)=macrc(j,Kt,I,m)+mact(j,kt,i)
                 end do
              end if
            end do
            do m=1,nmc
              tmac=0.0
              if(macrophyte_CALC(jw,m))then
                jt=kti(i)
                je=kb(i)
                do j=jt,je
                  tmac=tmac+macrm(j,kt,i,m)
                end do
              end if
              mac(kt,i,m)=tmac/(bh1(kt,i)*dlx(i))
            end do
            do m=1,nmc
              if(macrophyte_CALC(jw,m))then
                mac(kt-1,i,m)=0.0
              end if
            end do
          end do

          DO I=IU-1,ID
            AVHR(KT-1,I) =  H1(KT-1,I)+( H1(KT-1,I+1)- H1(KT-1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            AVHR(KT,I)   =  H1(KT,I)  +( H1(KT,I+1)  - H1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT,I)   = BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT-1,I) = BH1(KT-1,I)+(BH1(KT-1,I+1)-BH1(KT-1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
          END DO
          U(KT-1,IU-1:ID+1)     = 0.0
          W(KT-1,IU-1:ID+1)     = 0.0
          SW(KT-1,IU-1:ID+1)    = 0.0                                                                                  !TC 3/9/05
          P(KT-1,IU-1:ID+1)     = 0.0
          AZ(KT-1,IU-1:ID+1)    = 0.0
	  TKE(KT-1,IU-1:ID+1,1) = 0.0 !sg 10/4/07
	  TKE(KT-1,IU-1:ID+1,2) = 0.0 !sg 10/4/07
          DZ(KT-1,IU-1:ID+1)    = 0.0
          ADMZ(KT-1,IU-1:ID+1)  = 0.0
          ADZ(KT-1,IU-1:ID+1)   = 0.0
          DECAY(KT-1,IU-1:ID+1) = 0.0
          IF (UP_HEAD(JB)) THEN
            QUH1(KT,JB)             = QUH1(KT,JB)+QUH1(KT-1,JB)
            TSSUH1(KT,JB)           = TSSUH1(KT-1,JB)          +TSSUH1(KT,JB)
            CSSUH1(KT,CN(1:NAC),JB) = CSSUH1(KT-1,CN(1:NAC),JB)+CSSUH1(KT,CN(1:NAC),JB)
          END IF
          IF (DN_HEAD(JB)) THEN
            QDH1(KT,JB)             = QDH1(KT,JB)+QDH1(KT-1,JB)
            TSSDH1(KT,JB)           = TSSDH1(KT-1,JB)          +TSSDH1(KT,JB)
            CSSDH1(KT,CN(1:NAC),JB) = CSSDH1(KT-1,CN(1:NAC),JB)+CSSDH1(KT,CN(1:NAC),JB)
          END IF

!******** Upstream active segment

          IUT = US(JB)
          IF (SLOPE(JB) /= 0.0) THEN
            DO I=US(JB)-1,DS(JB)+1
              IF (KB(I) < KT ) THEN                                                                                  ! SR 10/17/05
                KB(I)                 = KT
                Bnew(KB(I),I)         = 0.000001   ! sw 1/23/06
                DX(KB(I),I)           = DXI(JW)
                ilayer(i)=1
                T1(KB(I),I)           = T1(KT,I)                   !    SW 5/15/06    T1(KB(I)-1,I)
                C1(KB(I),I,CN(1:NAC)) = C1(KT,I,CN(1:NAC))         !    SW 5/15/06    C1(KB(I)-1,I,CN(1:NAC))
                WRITE (WRN,'(2(A,I8),A,F0.3,A,F0.3)') 'Lowering bottom segment ',I,' at iteration ',NIT,' at Julian day ',&
                                                       JDAY,' Z(I)=',Z(I)
              END IF
            ENDDO                    ! SW 1/23/06
            DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
!               IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))  ! SW 1/23/06
                IF (I /= US(JB)-1) KBMIN(I-1) = MIN(KB(I-1),KB(I))  ! SW 1/23/06
                if(kbi(i) < kb(i))then
                bkt(i)=bh1(kt,i)/(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
                depthb(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
                depthm(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))*0.5    ! SW 1/23/06
                avhr(kt,i)=(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))   +(H1(KT,I+1)-(el(kbi(i)+1,i+1)-el(kb(i)+1,i+1))&
                                           -H1(KT,I)+(el(kbi(i)+1,i)-el(kb(i)+1,i)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                                                           ! SW 1/23/06
                end if
            ENDDO

        DO I=US(JB),DS(JB)   ! SW 1/23/06   11/13/07   US(JB)-1,DS(JB)+1

         if(ilayer(i).eq.1.and.ilayer(i+1).eq.0)then  ! SW 1/23/06
         bhrsum=0.0
         q(i)=0.0
          DO K=KT,KBMIN(I)
            IF (.NOT. INTERNAL_WEIR(K,I)) THEN
              BHRSUM = BHRSUM+BHR1(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR1(K,I)
            END IF
          END DO
          DO K=KT,KBMIN(I)
            IF (INTERNAL_WEIR(K,I)) THEN
              U(K,I) = 0.0
            ELSE
              U(K,I) =  U(K,I)+(QC(I)-Q(I))/BHRSUM
            END IF
          END DO
          elseif(ilayer(i).eq.1.and.ilayer(i-1).eq.0)then
          bhrsum=0.0
          q(i-1)=0.0
          DO K=KT,KBMIN(I-1)
            IF (.NOT. INTERNAL_WEIR(K,I-1)) THEN
              BHRSUM = BHRSUM+BHR1(K,I-1)
              Q(I-1)   = Q(I-1)+U(K,I-1)*BHR1(K,I-1)
            END IF
          END DO
          DO K=KT,KBMIN(I-1)
            IF (INTERNAL_WEIR(K,I-1)) THEN
              U(K,I-1) = 0.0
            ELSE
              U(K,I-1) =  U(K,I-1)+(QC(I-1)-Q(I-1))/BHRSUM
            END IF
          END DO
          endif    ! SW 1/23/06

            END DO
          END IF
          DO I=US(JB),DS(JB)
            IF (KB(I)-KT < NL(JB)-1) IUT = I+1
            ONE_LAYER(I) = KTWB(JW) == KB(I)
          END DO
          IF (IUT > DS(JB)) THEN
            WRITE (W2ERR,'(A,I0/A,F0.2,2(A,I0))') 'Fatal error - insufficient segments in branch ',JB,'Julian day = ',JDAY,      &
                                                  ' at iteration ',NIT,' with water surface layer = ',KT
            WRITE (W2ERR,'(2(A,I0))')             'Minimum water surface located at segment ',IZMIN(JW),' with bottom layer at ',&
                                                   KB(IZMIN(JW))
            TEXT = 'Runtime error - see w2.err'
            ERROR_OPEN = .TRUE.
            RETURN
       ! Go to 230
          END IF

!******** Segment subtraction

          IF (IUT /= IU) THEN
            IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/17X,A,I0,A,I0)') ' Subtract segments ',IU,' through ',IUT-1
            DO I=IU,IUT-1
              DO K=KT,KB(I)
                EBRI(JB)            = EBRI(JB)-T1(K,I)*VOL(K,I)
                CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)-C1(K,I,CN(1:NAC))*VOL(K,I)+(CSSB(K,I,CN(1:NAC))+CSSK(K,I,CN(1:NAC))   &
                                      *VOL(K,I))*DLT
              END DO
            END DO
! v3.5 start
            do i=iu,iut-1
              do m=1,nmc
                if(macrophyte_calc(jw,m))then
                  jt=kti(i)
                  je=kb(i)
                  do j=jt,je
                    if(j.lt.kt)then
                      colb=el(j+1,i)
                    else
                      colb=el(kt+1,i)
                    end if
                    coldep=ELWS(i)-colb
!                    maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(j,kt,i,m)+(macSS(j,Kt,I,m)*coldep*cw(j,i)*DLX(I))*DLT
                    maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(j,kt,i,m)
                  end do
                  DO K=KT+1,KB(I)
                    jt=k
                    je=kb(i)
                    do j=jt,je
!                      maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(j,k,i,m)+(macSS(j,K,I,m)*H2(k,i)*cw(j,i)*DmX(I))*DLT
                      maCMBRT(JB,m) = maCMBRT(JB,m)-macrm(j,k,i,m)
                    end do
                  END DO
                end if
              end do
            end do

            F(IU-1:IUT-1)     =  0.0
            Z(IU-1:IUT-1)     =  0.0
            ICETH(IU-1:IUT-1) =  0.0
            BHRHO(IU-1:IUT-1) =  0.0
            ICE(IU-1:IUT-1)   = .FALSE.
            DO K=KT,KB(IUT)
              ADX(K,IU-1:IUT-1)            = 0.0
              DX(K,IU-1:IUT-1)             = 0.0
              AZ(K,IU-1:IUT-1)             = 0.0
	      TKE(K,IU-1:IUT-1,1)          = 0.0 !SG  10/4/07
	      TKE(K,IU-1:IUT-1,2)          = 0.0 !SG  10/4/07
              SAZ(K,IU-1:IUT-1)            = 0.0
              U(K,IU-1:IUT-1)              = 0.0
              SU(K,IU-1:IUT-1)             = 0.0
              T1(K,IU-1:IUT-1)             = 0.0
              TSS(K,IU-1:IUT-1)            = 0.0
              QSS(K,IU-1:IUT-1)            = 0.0
              C1(K,IU-1:IUT-1,CN(1:NAC))   = 0.0
              C2(K,IU-1:IUT-1,CN(1:NAC))   = 0.0
              C1S(K,IU-1:IUT-1,CN(1:NAC))  = 0.0
              CSSB(K,IU-1:IUT-1,CN(1:NAC)) = 0.0
              CSSK(K,IU-1:IUT-1,CN(1:NAC)) = 0.0
            END DO

            do m=1,nmc
              IF (macrophyte_CALC(jw,m)) THEN
                mac(k,i,m)=0.0
                mact(j,k,i)=0.0
              end if
            end do
            jt=kti(i)
            je=kb(i)
            do j=jt,je
              do m=1,nmc
                if(macrophyte_calc(jw,m))then
                  macrc(j,K,I,m) = 0.0
                end if
              end do
            end do

            IU           =  IUT
            CUS(JB)      =  IU
            Z(IU-1)      = (EL(KT,IU-1)-(EL(KT,IU)-Z(IU)*COSA(JB)))/COSA(JB)
            SZ(IU-1)     =  Z(IU)
            KTI(IU-1)    =  KTI(IU)
            IF (.NOT. TRAPEZOIDAL(JW)) THEN
              BI(KT,IU-1)  = B(KTI(IU-1),I)
              H1(KT,IU-1)  = H(KT,JW)-Z(IU-1)
              BH1(KT,IU-1) = Bnew(KTI(IU-1),IU-1)*(EL(KT,IU-1)-EL(KTI(IU-1)+1,IU-1)-Z(IU-1)*COSA(JB))/COSA(JB)   ! sw 1/23/06  Bnew(KTI(IU-1),IU-1)*(EL(KT,IU-1)-EL(KTI(IU-1)+1,IU-1)-Z(IU-1)*COSA(JB))/COSA(JB)     ! SR 10/17/05
              IF (KT >= KB(IU-1)) BH1(KT,IU-1) = Bnew(KT,IU-1)*H1(KT,IU-1)   ! sw 1/23/06
              DO K=KTI(IU-1)+1,KT
                BH1(KT,IU-1) = BH1(KT,IU-1)+BH1(K,IU-1)
              END DO
            ELSE
              CALL GRID_AREA1 (EL(KT,I)-Z(I),EL(KT+1,IU-1),BH1(KT,IU-1),BI(KT,IU-1))                                             !SW 08/03/04
              BH1(KT,I) = 0.25*H(KT,JW)*(BB(KT-1,I)+2.*B(KT,I)+BB(KT,I))
              H1(KT,I)  = H(KT,JW)-Z(I)
            END IF
            BKT(IU-1)     =  BH1(KT,IU-1)/H1(KT,IU-1)
            BHR1(KT,IU-1) =  BH1(KT,IU-1)+(BH1(KT,IU)-BH1(KT,IU-1))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                 !SW 07/29/04
            IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)
            IF (UH_INTERNAL(JB)) THEN
              IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
                KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
              ELSE
                DO KKB = KT, KMX
                  IF (EL(KKB,IU) <= EL(KB(UHS(JB)),UHS(JB))) EXIT
                END DO
                KB(IU-1) = MIN(KKB,KB(IU))
              END IF
            END IF
          END IF
          IF (CONSTITUENTS) THEN      ! SW 5/15/06
            CALL TEMPERATURE_RATES
            CALL KINETIC_RATES
          END IF


!******** Total active cells

          DO I=IU,ID
            NTAC = NTAC-1
          END DO
          NTACMN = MIN(NTAC,NTACMN)
        END DO
        CALL INTERPOLATION_MULTIPLIERS

!****** Additional layer subtractions

        ZMIN(JW) = -1000.0
        DO JB=BS(JW),BE(JW)
          DO I=CUS(JB),DS(JB)
            ZMIN(JW) = MAX(ZMIN(JW),Z(I))
          END DO
        END DO
        SUB_LAYER = ZMIN(JW) > 0.60*H(KT,JW) .AND. KT < KTMAX                                                         ! SR 10/17/05
      END DO
    END DO

!** Temporary downstream head segment

    DO JB=1,NBR
      IF (DHS(JB) > 0) THEN
        DO JJB=1,NBR
          IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB)) EXIT
        END DO
        IF (CUS(JJB) > DHS(JB)) CDHS(JB) = CUS(JJB)
      END IF
    END DO
return
end subroutine layeraddsub
