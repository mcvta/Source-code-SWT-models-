!v3.5 start
!************************************************************************
!**               S U B R O U T I N E    porosity                      **
!************************************************************************

      SUBROUTINE porosity

      use geomc;use global;use macrophytec;use porosityc;use LOGICC
      use screenc

      INTEGER :: K,M,KUP,IEXIT,K1,KDN,IUT
      REAL    :: B11,VSTOT,ELR,ELL,ELL1,ELL2,ELR2,EL1,EL2,ERR1,ERR2

    if(nit.eq.0)then

 1040 FORMAT((8X,i8,3F8.0))

      DO Jw=1,Nwb
        KT = KTwb(Jw)
        DO JB=BS(Jw),BE(Jw)
          IU   = CUS(JB)
          ID   = DS(JB)
          DO I=IU,ID
            DO K=2,KB(I)
              VOLi(K,I) = BH(K,I)*DLX(I)
            END DO
            VOLi(KT,I)    = BH2(KT,I)*DLX(I)
          end do
        end do
      end do

    end if

    do jb=1,nbr
      cosa(jb)=cos(alpha(jb))
    end do

!c  calculating # of macrophyte stems in each cell

    DO Jw=1,Nwb
      KT = KTwb(Jw)
      DO JB=BS(Jw),BE(Jw)
        IU   = CUS(JB)
        ID   = DS(JB)
        DO I=IU,ID
!          HkTi  = H(KT,jw)-Z(I)    replaced by h1(kt,i)
          if(kt.eq.kti(i))then
            volkti(i)=h1(kt,i)*bic(kt,i)*dlx(i)
          else
            volkti(i) = bic(KTI(I),I)*(EL(KT,i)-EL(KTI(I)+1,i)-Z(I)*cosa(jb))/cosa(jb)*dlx(i)
          end if
          DO K=KTI(I)+1,KT
            volkti(I) = volkti(I)+voli(k,i)
          END DO

          do m=1,nmc
            vstemkt(i,m)=(mac(kt,i,m)*volkti(i))/dwv(m)    !cb 6/29/06
          end do

          do k=kt+1,kb(i)
            do m=1,nmc
              vstem(k,i,m)=(mac(k,i,m)*voli(k,i))/dwv(m)   !cb 6/29/06
            end do
          end do
        end do
      end do
    end do

    por=1.0
    DO Jw=1,Nwb
      KT = KTwb(Jw)
      DO JB=BS(Jw),BE(Jw)
        IU = cUS(JB)
        ID = DS(JB)

        do i=iu,id
          do k=kt,kb(i)
            if(k.eq.kt)then
              vstot=0.0
              do m=1,nmc
                  vstot=vstot+vstemkt(i,m)
              end do
              por(kt,i)=(VOLkti(I)-vstot)/volkti(i)
            else
              vstot=0.0
              do m=1,nmc
                vstot=vstot+vstem(k,i,m)
              end do
              por(k,i)=(VOLi(K,I)-vstot)/voli(k,i)
            end if
          end do
        end do

        do i=iu,id
          do k=kti(i),kb(i)
            if(k.le.kt)then
              b(k,i)=por(kt,i)*bic(k,i)
            else
              b(k,i)=por(k,i)*bic(k,i)
            end if

          end do
        end do

      end do
    end do



! Boundary widths

  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU = US(JB)
      ID = DS(JB)
      DO I=IU-1,ID+1
        B(1,I) = B(2,I)
        DO K=KB(I)+1,KMX
          B(K,I) = B(KB(I),I)
        END DO
      END DO
    END DO
  END DO
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU    = US(JB)
      ID    = DS(JB)
      IEXIT = 0
      DO K=1,KMX-1
        B(K,IU-1) = B(K,IU)
        IF (UH_INTERNAL(JB) .OR. HEAD_FLOW(JB)) THEN
          IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
            B(K,IU-1) = B(K,UHS(JB))
          ELSE
            ELR = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
            ELL = EL(2,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
            IF (ELR >= ELL) THEN
              B(K,IU-1) = B(2,UHS(JB))
            ELSE
              DO KUP=2,KMX-1
                ELL1 = EL(KUP,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                ELL2 = EL(KUP+1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                IF (ELL1 > ELR .AND. ELL2 <= ELR) THEN
                  IF (KUP > KB(UHS(JB))) THEN
                    KB(IU-1)    = K-1
                    KBMIN(IU-1) = MIN(KB(IU),KB(IU-1))
                    IEXIT       = 1
                    EXIT
                  END IF
                  ELR2 = EL(K+1,IU)+SINA(JB)*DLX(IU)*0.5
                  IF (ELR2 >= ELL2) THEN
                    B(K,IU-1) = B(KUP,UHS(JB)); EXIT
                  ELSE
                    K1 = KUP+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL1 = ELR
                    EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(IU)*0.5
                    DO WHILE (ELR2 <= EL2)
                      B11 = B11+(EL1-EL2)*B(K1-1,UHS(JB))
                      EL1 = EL2
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL2 == ELR2) EXIT
                      EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                      IF (EL2 <= ELR2) EL2 = ELR2
                    END DO
                    B(K,IU-1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,UHS(JB)) > EL(K,IU)) B(K,IU-1) = B(K-1,IU-1)
              IF (B(K,IU-1) == 0.0) B(K,IU-1) = B(K-1,IU-1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO
      IEXIT = 0
      DO K=1,KMX-1
        B(K,ID+1) = B(K,ID)
        IF (DH_INTERNAL(JB)) THEN
          IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
            B(K,ID+1) = B(K,DHS(JB))
          ELSE
            ELL = EL(K,ID)-SINA(JB)*DLX(ID)*0.5
            ELR = EL(2,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
            IF (ELL >= ELR) THEN
              B(K,ID+1) = B(2,DHS(JB))
            ELSE
              DO KDN=2,KMX-1
                ERR1 = EL(KDN,DHS(JB))  +SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                ERR2 = EL(KDN+1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                IF (ERR1 >= ELL .AND. ERR2 < ELL) THEN
                  IF (KDN > KB(DHS(JB))) THEN
                    KB(ID+1)  = K-1
                    KBMIN(ID) = MIN(KB(ID),KB(ID+1))
                    IEXIT     = 1
                    EXIT
                  END IF
                  ELL2 = EL(K+1,ID)-SINA(JB)*DLX(ID)*0.5
                  IF (ELL2 >= ERR2) THEN
                    B(K,ID+1) = B(KDN,DHS(JB)); EXIT
                  ELSE
                    K1  = KDN+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL2 = ELL
                    EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                    DO WHILE (ELL2 <= EL1)
                      B11 = B11+(EL2-EL1)*B(K1-1,DHS(JB))
                      EL2 = EL1
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL1 == ELL2) EXIT
                      EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                      IF (EL1 <= ELL2) EL1 = ELL2
                    END DO
                    B(K,ID+1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,DHS(JB)) > EL(K,ID)) B(K,ID+1) = B(K-1,ID+1)
              IF (B(K,ID+1) == 0.0) B(K,ID+1) = B(K-1,ID+1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO

!**** Areas and bottom widths

      IF (.NOT. TRAPEZOIDAL(JW)) THEN                                                                                  !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BH2(K,I) = B(K,I)*H(K,JW)
            BH(K,I)  = B(K,I)*H(K,JW)
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5                                  !SW 08/02/04
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO
!****** Derived geometry

        DO I=IU-1,ID+1
          BH2(KT,I) = B(KTI(I),I)*(EL(KT,I)-EL(KTI(I)+1,I)-Z(I)*COSA(JB))/COSA(JB)
          IF (KT == KTI(I)) BH2(KT,I) = H2(KT,I)*B(KT,I)
          DO K=KTI(I)+1,KT
            BH2(KT,I) = BH2(KT,I)+BH(K,I)
          END DO
          BKT(I)   = BH2(KT,I)/H2(KT,I)
          BI(KT,I) = B(KTI(I),I)
        END DO
      ELSE                                                                                                             !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5
          END DO
          BB(KB(I),I) = B(KB(I),I)*0.5
          BH2(1,I)    = B(1,I)*H(1,JW)
          BH(1,I)     = BH2(1,I)
          DO K=2,KMX-1
            BH2(K,I) = 0.25*H(K,JW)*(BB(K-1,I)+2.*B(K,I)+BB(K,I))
            BH(K,I)  = BH2(K,I)
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO
        DO I=IU-1,ID+1
          CALL GRID_AREA1(EL(KT,I)-Z(I),EL(KT+1,I),BH2(KT,I),BI(KT,I))
          BKT(I) = BH2(KT,I)/H2(KT,I)
        END DO
      END IF
      DO I=IU-1,ID
        DO K=1,KMX-1
          AVH2(K,I) = (H2(K,I) +H2(K+1,I)) *0.5
          AVHR(K,I) =  H2(K,I)+(H2(K,I+1)-H2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                  !SW 07/29/04
        END DO
        AVH2(KMX,I) = H2(KMX,I)
        DO K=1,KMX
          BR(K,I)   = B(K,I)  +(B(K,I+1)  -B(K,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR(K,I)  = BH(K,I) +(BH(K,I+1) -BH(K,I)) /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR2(K,I) = BH2(K,I)+(BH2(K,I+1)-BH2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
        END DO
      END DO
      DO K=1,KMX-1
        AVH2(K,ID+1) = (H2(K,ID+1)+H2(K+1,ID+1))*0.5
        BR(K,ID+1)   =   B(K,ID+1)
        BHR(K,ID+1)  =   BH(K,ID+1)
      END DO
      AVH2(KMX,ID+1) = H2(KMX,ID+1)
      AVHR(KT,ID+1)  = H2(KT,ID+1)
      BHR2(KT,ID+1)  = BH2(KT,ID+1)
      IUT = IU
      IF (UP_HEAD(JB)) IUT = IU-1
      DO I=IUT,ID
        DO K=1,KMX-1
          VOL(K,I) = B(K,I)*H2(K,I)*DLX(I)
        END DO
        VOL(KT,I)    = BH2(KT,I)*DLX(I)
        DEPTHB(KT,I) = H2(KT,I)
        DEPTHM(KT,I) = H2(KT,I)*0.5
        DO K=KT+1,KMX
          DEPTHB(K,I) = DEPTHB(K-1,I)+ H2(K,I)
          DEPTHM(K,I) = DEPTHM(K-1,I)+(H2(K-1,I)+H2(K,I))*0.5
        END DO
      END DO
    END DO
  END DO

10 continue

    return
    end

!************************************************************************
!**               S U B R O U T I N E    macrophyte_friction           **
!************************************************************************

      SUBROUTINE macrophyte_friction(hrad,bedfr,effric,k,ii)

    use geomc;use global;use macrophytec;use porosityc
      INTEGER :: K,M,II
      REAL    :: SAVOLRAT,XSAREA,TSAREA,ARTOT,SCTOT,CDAVG,FRIN,HRAD,BEDFR,EFFRIC

  do m=1,nmc
    savolrat=dwv(m)/dwsa(m)     !cb 6/29/2006
    if(k.eq.kt)then
!      sarea(m)=vstemkt(ii,m)*savolrat/pi
      sarea(m)=vstemkt(ii,m)*savolrat*anorm(m)     !cb 6/29/2006
    else
!      sarea(m)=vstem(k,ii,m)*savolrat/pi
      sarea(m)=vstem(k,ii,m)*savolrat*anorm(m)     !cb 6/29/2006
    end if
  end do
  xsarea=bh2(k,ii)

  tsarea=0.0
  artot=0.0
  sctot=0.0
  do m=1,nmc
    artot=artot+sarea(m)
    sctot=sctot+cddrag(m)*sarea(m)
    tsarea=tsarea+sarea(m)
  end do

  if(artot.gt.0.0)then
    cdavg=sctot/artot
    frin=cdavg*tsarea*hrad**(4./3.)/(2.0*g*xsarea*dlx(ii)*bedfr**2)
    effric=bedfr*sqrt(1.0+frin)
  else
    effric=bedfr
  end if

  return
  end

! v3.5
