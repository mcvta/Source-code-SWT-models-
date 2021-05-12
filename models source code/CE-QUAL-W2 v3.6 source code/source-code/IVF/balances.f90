subroutine balances

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  !v3.5
  EXTERNAL RESTART_OUTPUT
  real volinjw,volprjw,voloutjw,volwdjw,volevjw,voldtjw,voltrbjw

!***********************************************************************************************************************************
!*                                                    Task 2.6: Balances                                                          **
!***********************************************************************************************************************************

    QINT  = 0.0
    QOUTT = 0.0
    VOLSR = 0.0
    VOLTR = 0.0
 ! output flow loading and pollutant loading debugging information output frequency at output of cpl

    if(nit==0)then
    open(9525,file='flowbal.opt',status='unknown')
    write(9525,*)'WB,JDAY,VOLIN,volpr,volout,volwd,volev,voldt,voltrb'
    endif

    DO JW=1,NWB
    volinjw=0.0
    volprjw=0.0
    voloutjw=0.0
    volwdjw=0.0
    volevjw=0.0
    voldtjw=0.0
    voltrbjw=0.0

      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IF (VOLUME_BALANCE(JW)) THEN
          VOLSBR(JB) = VOLSBR(JB)+DLVOL(JB)
          VOLTBR(JB) = VOLEV(JB)+VOLPR(JB)+VOLTRB(JB)+VOLDT(JB)+VOLWD(JB)+VOLUH(JB)+VOLDH(JB)+VOLIN(JB)+VOLOUT(JB)
          VOLSR(JW)  = VOLSR(JW)+VOLSBR(JB)
          VOLTR(JW)  = VOLTR(JW)+VOLTBR(JB)
          QINT(JW)   = QINT(JW) +VOLIN(JB)+VOLTRB(JB)+VOLDT(JB)+VOLPR(JB)
          QOUTT(JW)  = QOUTT(JW)-VOLEV(JB)-VOLWD(JB) -VOLOUT(JB)
          IF (ABS(VOLSBR(JB)-VOLTBR(JB)) > VTOL .AND. VOLTBR(JB) > 100.0*VTOL) THEN
            IF (VOLUME_WARNING) THEN
              WRITE (WRN,'(A,F0.3,3(:/A,E15.8,A))') 'Computational warning at Julian day = ',JDAY,'spatial change  =', VOLSBR(JB), &
                                                    ' m^3','temporal change =',VOLTBR(JB),' m^3','volume error    =',              &
                                                     VOLSBR(JB)-VOLTBR(JB),' m^3'
              WRITE(WRN,*)'LAYER CHANGE:',LAYERCHANGE(JW)
              WRITE(WRN,*)'SZ',sz,'Z',z,'H2KT',h2(kt,1:imx),'H1KT',h1(kt,1:imx),'WSE',elws,'Q',q,'QC',qc,'T1',t1(kt,1:imx),'T2',&
                           t2(kt,1:imx),'SUKT',su(kt,1:imx),'UKT',u(kt,1:imx),'QIN',qinsum,'QTR',qtr,'QWD',qwd
              WARNING_OPEN   = .TRUE.
              VOLUME_WARNING = .FALSE.
            END IF
          END IF
        END IF
        IF (VOLSR(JW) /= 0.0) DLVR(JW) = (VOLTR(JW)-VOLSR(JW))/VOLSR(JW)*100.0
            volinjw=volinjw+volin(jb)
            volprjw=volprjw+volpr(jb)
            voloutjw=voloutjw+volout(jb)
            volwdjw=volwdjw+volwd(jb)
            volevjw=volevjw+volev(jb)
            voldtjw=voldtjw+voldt(jb)
            voltrbjw=voltrbjw+voltrb(jb)
      END DO
        IF (CONTOUR(JW)) THEN
        IF (JDAY+(DLT/DAY) >= NXTMCP(JW) .OR. JDAY+(DLT/DAY) >= CPLD(CPLDP(JW)+1,JW))write(9525,'(i3,",",1x,f10.3,",",10(e16.8,",",1x))')JW,JDAY,volinjw,volprjw,voloutjw,volwdjw,volevjw,voldtjw,voltrbjw
        ENDIF
      IF (ENERGY_BALANCE(JW)) THEN
        ESR(JW) = 0.0
        ETR(JW) = 0.0
        DO JB=BS(JW),BE(JW)
          ETBR(JB) = EBRI(JB)+TSSEV(JB)+TSSPR(JB)+TSSTR(JB)+TSSDT(JB)+TSSWD(JB)+TSSUH(JB)+TSSDH(JB)+TSSIN(JB)+TSSOUT(JB)+TSSS(JB)  &
                     +TSSB(JB)+TSSICE(JB)
          ESBR(JB) = 0.0
          DO I=CUS(JB),DS(JB)
            DO K=KT,KB(I)
              ESBR(JB) = ESBR(JB)+T1(K,I)*DLX(I)*BH1(K,I)
            END DO
          END DO
          ETR(JW) = ETR(JW)+ETBR(JB)
          ESR(JW) = ESR(JW)+ESBR(JB)
        END DO
      END IF
      IF (MASS_BALANCE(JW)) THEN
        DO JB=BS(JW),BE(JW)
          DO JC=1,NAC
            CMBRS(CN(JC),JB) = 0.0
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                CMBRS(CN(JC),JB) = CMBRS(CN(JC),JB)+C1(K,I,CN(JC))*DLX(I)*BH1(K,I)
                CMBRT(CN(JC),JB) = CMBRT(CN(JC),JB)+(CSSB(K,I,CN(JC))+CSSK(K,I,CN(JC))*BH1(K,I)*DLX(I))*DLT
              END DO
            END DO
          END DO
! v3.5 start
          do m=1,nmc
            if(macrophyte_calc(jw,m))then
              maCMBRS(JB,m) = 0.0
              DO I=CUS(JB),DS(JB)
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
                  maCMBRS(JB,m) = maCMBRS(jb,m)+macrm(j,Kt,I,m)
                  maCMBRT(JB,m) = maCMBRT(JB,m)+(macSS(j,Kt,I,m)*coldep*cw(j,i)*DLX(I))*DLT
                end do
                DO K=KT+1,KB(I)
                  jt=k
                  je=kb(i)
                  do j=jt,je
                    maCMBRS(JB,m) =maCMBRS(JB,m)+macrm(j,K,I,m)
!                    maCMBRT(JB,m) = maCMBRT(JB,m)+(macSS(j,K,I,m)*H2(k,i)*cw(j,i)*DLX(I))*DLT
                    maCMBRT(JB,m) = maCMBRT(JB,m)+(macSS(j,K,I,m)*(cw(j,i)/b(k,i))*bh1(k,i)*DLX(I))*DLT
                  end do
                END DO
              END DO
            end if
          end do
! v3.5 end
        END DO
      END IF
    END DO

    return
    end subroutine balances
