
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   W I T H D R A W A L                                          **
!***********************************************************************************************************************************

SUBROUTINE WITHDRAWAL
  USE GLOBAL; USE GEOMC; USE TVDC; USE SELWC; USE LOGICC
  real :: hswt,hswb,elr,wsel,elstr,coef,ratio,ht,rhoft,dlrhot,hb,rhofb,dlrhob,vsum,dlrhomax,hwdt,hwdb,elwd,tempest,estrtest
  integer :: k,js,kstr,ktop,kbot,jwd,kwd

RETURN

!***********************************************************************************************************************************
!**                                             D O W N S T R E A M   W I T H D R A W A L                                         **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_WITHDRAWAL (JS)

! Variable initialization

  HSWT = 0.0; HSWB = 0.0; VNORM = 0.0; QNEW = 0.0

! Water surface elevation

  ELR  = SINA(JB)*DLX(ID)*0.5
  WSEL = ELWS(ID)-ELR                   !EL(KT,ID)-Z(ID)*COSA(JB)

! Structure layer

  DO K=KT,KB(ID)
    IF (EL(K,ID)-ELR < ESTR(JS,JB)) EXIT
  END DO
  KSTR = MAX(K-1,KT)
  KSTR = MIN(KSTR,KB(ID))

! Initial withdrawal limits

  KTOP = MAX(KTSW(JS,JB),KT)
  IF (KSTR < KTOP) KTOP = KSTR
  KBOT = MIN(KBSW(JS,JB),KB(ID))
  IF (KBOT <= KT .AND. KBOT /= KB(ID)) KBOT = KT+1
  IF (KBOT > KB(ID)) KBOT = KB(ID)
  ELSTR = ESTR(JS,JB)
  IF (ESTR(JS,JB) <= EL(KB(ID)+1,ID+1)-ELR) THEN
    KSTR  = KB(ID)
    ELSTR = EL(KB(ID),ID)-ELR
  END IF
  IF (ESTR(JS,JB) > EL(KT,ID)-ELR) ELSTR = WSEL
  IF (KBSW(JS,JB) < KSTR) THEN
    KSTR  = KT
    ELSTR = WSEL
  END IF

! Boundary interference

  COEF = 1.0
  IF ((WSEL-EL(KBOT,ID)-ELR) /= 0.0) THEN
    RATIO = (ELSTR-(EL(KBOT,ID)-ELR))/(WSEL-(EL(KBOT,ID)-ELR))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KSTR-1,KTOP,-1

!** Density frequency

    HT    = (EL(K,ID)-ELR)-ELSTR
    RHOFT = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HT*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWT = (COEF*QSTR(JS,JB)/RHOFT)**0.333333
    ELSE
      HSWT = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFT))
    END IF
    IF (HT >= HSWT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR+HSWT) < WSEL) THEN
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
  ELSE IF (WSEL == ELSTR) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*HSWT/(WSEL-ELSTR)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KSTR+1,KBOT

!** Density frequency

    HB    = ELSTR-(EL(K,ID)-ELR)
    RHOFB = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HB*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWB = (COEF*QSTR(JS,JB)/RHOFB)**0.333333
    ELSE
      HSWB = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFB))
    END IF
    IF (HB >= HSWB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR-HSWB) > EL(KBOT+1,ID)) THEN
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
  ELSE
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))*HSWB/(ELSTR-(EL(KBOT+1,ID)-ELR))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                      ! GH 1/31/08
  DO K=KTOP,KBOT
!    VNORM(K) = ABS(1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2)*BHR2(K,ID)
 	   IF(K.GT.KSTR)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
     VNORM(K) = 1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2
 	 IF(VNORM(K).GT.1.0) VNORM(K)=1.0                         !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                         !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,ID)
     VSUM     = VSUM+VNORM(K)
  END DO

! Outflows
  qsumjs=0.0                                                  ! SW 7/30/09
  DO K=KTOP,KBOT
    QNEW(K)    = (VNORM(K)/VSUM)*QSTR(JS,JB)
    QOUT(K,JB) =  QOUT(K,JB)+QNEW(K)
    tavg(js,jb)=tavg(js,jb)+qnew(k)*t2(k,id)                  ! SW 7/30/09
    qsumjs=qsumjs+qnew(k)
  END DO
if(qsumjs.gt.0.0)tavg(js,jb)=tavg(js,jb)/qsumjs               ! SW 7/30/09

! Inactive layers and total outflow

  IF (JS == NST) THEN
    WHERE (QOUT(:,JB) == 0.0) U(:,ID) = 0.0
  END IF
RETURN
!***********************************************************************************************************************************
!**                                             D O W N S T R E A M   W I T H D R A W A L  ESTIMATE                               **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_WITHDRAWAL_estimate(JS,tempest,estrtest)

! Variable initialization

  HSWT = 0.0; HSWB = 0.0; VNORM = 0.0; QNEW = 0.0

! Water surface elevation

  ELR  = SINA(JB)*DLX(ID)*0.5
  WSEL = EL(KT,ID)-Z(ID)*COSA(JB)-ELR

! Structure layer

  DO K=KT,KB(ID)
    IF (EL(K,ID)-ELR < estrtest) EXIT
  END DO
  KSTR = MAX(K-1,KT)
  KSTR = MIN(KSTR,KB(ID))

! Initial withdrawal limits

  KTOP = MAX(KTSW(JS,JB),KT)
  IF (KSTR < KTOP) KTOP = KSTR
  KBOT = MIN(KBSW(JS,JB),KB(ID))
  IF (KBOT <= KT .AND. KBOT /= KB(ID)) KBOT = KT+1
  IF (KBOT > KB(ID)) KBOT = KB(ID)                                                                                     !SW 06/03/02
  ELSTR = estrtest
  IF (estrtest <= EL(KB(ID)+1,ID+1)-ELR) THEN                                                                       !SW 10/17/01
    KSTR  = KB(ID)
    ELSTR = EL(KB(ID),ID)-ELR                                                                                          !SW 10/17/01
  END IF
  IF (estrtest > EL(KT,ID)-ELR) ELSTR = WSEL
  IF (KBSW(JS,JB) < KSTR) THEN
    KSTR  = KT
    ELSTR = WSEL                                                                                                       !SW 10/05/00
  END IF

! Boundary interference

  COEF = 1.0
  IF ((WSEL-EL(KBOT,ID)-ELR) /= 0.0) THEN
    RATIO = (ELSTR-(EL(KBOT,ID)-ELR))/(WSEL-(EL(KBOT,ID)-ELR))                                                         !SW 10/17/01
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KSTR-1,KTOP,-1

!** Density frequency

    HT    = (EL(K,ID)-ELR)-ELSTR
    RHOFT = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HT*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWT = (COEF*QSTR(JS,JB)/RHOFT)**0.333333
    ELSE
      HSWT = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFT))
    END IF
    IF (HT >= HSWT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR+HSWT) < WSEL) THEN
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
  ELSE IF (WSEL == ELSTR) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*HSWT/(WSEL-ELSTR)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KSTR+1,KBOT

!** Density frequency

    HB    = ELSTR-(EL(K,ID)-ELR)                                                                                       !SW 10/17/01
    RHOFB = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HB*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWB = (COEF*QSTR(JS,JB)/RHOFB)**0.333333
    ELSE
      HSWB = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFB))
    END IF
    IF (HB >= HSWB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR-HSWB) > EL(KBOT+1,ID)) THEN
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
  ELSE
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))*HSWB/(ELSTR-(EL(KBOT+1,ID)-ELR))                                           !SW 10/17/01
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)
  DO K=KTOP,KBOT
  IF(K.GT.KSTR)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
       VNORM(K) = 1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2
     IF(VNORM(K).GT.1.0) VNORM(K)=1.0                          !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                          !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,ID)
     VSUM= VSUM+VNORM(K)
  END DO

! Outflows

  tempest=0.0
  DO K=KTOP,KBOT
    tempest=tempest+t2(k,id)*(VNORM(K)/VSUM)*QSTR(JS,JB)
  END DO

  if(qstr(js,jb).gt.0.0)tempest=tempest/qstr(js,jb)

RETURN
!***********************************************************************************************************************************
!**                                                L A T E R A L   W I T H D R A W A L                                            **
!***********************************************************************************************************************************

ENTRY LATERAL_WITHDRAWAL (JWD)

! Variable initialization

  VNORM = 0.0; QSW(:,JWD) = 0.0; HWDT = 0.0; HWDB = 0.0

! Structure layer

  K = KT
  DO K=KT,KB(I)
    IF (EL(K,I) < EWD(JWD)) EXIT
  END DO
  KWD = MAX(K-1,KT)
  KWD = MIN(KWD,KB(I))

! Initial withdrawal limits

  KTOP = MAX(KTWD(JWD),KT)
  IF (KWD < KTOP) KTOP = KWD
  KBOT = MIN(KBWD(JWD),KB(I))
  IF (KBOT <= KT .AND. KB(I) /= KBOT) KBOT = KT+1
  IF (KBOT > KB(I)) KBOT = KB(I)
  ELWD = EWD(JWD)
  IF (EWD(JWD) <= EL(KB(I)+1,I)) THEN
    KWD  = KB(I)
    ELWD = EL(KB(I),I)
  END IF
  IF (EWD(JWD) > EL(KT,I)) ELWD = EL(KT,I)
  IF (KBWD(JWD) < KWD) THEN
    KWD  = KT
    ELWD = EL(KT,I)
  END IF

! Boundary interference

  COEF = 1.0
  IF (KT /= KBOT) THEN
    RATIO = (ELWD-EL(KBOT,I))/(EL(KT,I)-EL(KBOT,I))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KWD-1,KTOP,-1

!** Density frequency

    HT    = EL(K,I)-ELWD
    RHOFT = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HT*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDT = (COEF*QWD(JWD)/RHOFT)**0.333333
    IF (HT >= HWDT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD+HWDT) < EL(KT,I)) THEN
    DLRHOT = ABS(RHO(KWD,I)-RHO(KTOP,I))
  ELSE IF (EL(KT,I) == ELWD) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KWD,I)-RHO(KT,I))*HWDT/(EL(KT,I)-ELWD)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KWD+1,KBOT

!** Density frequency

    HB    = ELWD-EL(K,I)
    RHOFB = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HB*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDB = (COEF*QWD(JWD)/RHOFB)**0.333333
    IF (HB >= HWDB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD-HWDB) > EL(KBOT+1,I)) THEN
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))
  ELSE
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))*HWDB/(ELWD-EL(KBOT+1,I))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                                                                             ! SW 1/24/05
  DO K=KTOP,KBOT
!    VNORM(K) = ABS(1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2)*BHR2(K,I)
 	   IF(K.GT.KWD)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
     VNORM(K) = 1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2
 	 IF(VNORM(K).GT.1.0) VNORM(K)=1.0                         !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                         !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,I)
     VSUM     = VSUM+VNORM(K)
  END DO

! Outflows
  qsumwd=0.0                                                  ! SW 7/30/09

  DO K=KTOP,KBOT
    QSW(K,JWD) = QSW(K,JWD)+(VNORM(K)/VSUM)*QWD(JWD)
    tavgw(jwd)=tavgw(jwd)+(VNORM(K)/VSUM)*QWD(JWD)*t2(k,i)                  ! SW 7/30/09
    qsumwd=qsumwd+(VNORM(K)/VSUM)*QWD(JWD)
  END DO
  if(qsumwd.gt.0.0)tavgw(jwd)=tavgw(jwd)/qsumwd               ! SW 7/30/09
  KTW(JWD) = KTOP
  KBW(JWD) = KBOT
  return
!***********************************************************************************************************************************
!**                                                L A T E R A L   W I T H D R A W A L ESTIMATE                                   **
!***********************************************************************************************************************************

  ENTRY LATERAL_WITHDRAWAL_ESTIMATE (JWD,tempest,estrtest)

! Variable initialization

  VNORM = 0.0; QSW(:,JWD) = 0.0; HWDT = 0.0; HWDB = 0.0

! Structure layer

  K = KT
  DO K=KT,KB(I)
    IF (EL(K,I) < estrtest) EXIT
  END DO
  KWD = MAX(K-1,KT)
  KWD = MIN(KWD,KB(I))

! Initial withdrawal limits

  KTOP = MAX(KTWD(JWD),KT)
  IF (KWD < KTOP) KTOP = KWD
  KBOT = MIN(KBWD(JWD),KB(I))
  IF (KBOT <= KT .AND. KB(I) /= KBOT) KBOT = KT+1
  IF (KBOT > KB(I)) KBOT = KB(I)
  ELWD = estrtest
  IF (estrtest <= EL(KB(I)+1,I)) THEN
    KWD  = KB(I)
    ELWD = EL(KB(I),I)
  END IF
  IF (estrtest > EL(KT,I)) ELWD = EL(KT,I)
  IF (KBWD(JWD) < KWD) THEN
    KWD  = KT
    ELWD = EL(KT,I)
  END IF

! Boundary interference

  COEF = 1.0
  IF (KT /= KBOT) THEN
    RATIO = (ELWD-EL(KBOT,I))/(EL(KT,I)-EL(KBOT,I))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KWD-1,KTOP,-1

!** Density frequency

    HT    = EL(K,I)-ELWD
    RHOFT = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HT*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDT = (COEF*QWD(JWD)/RHOFT)**0.333333
    IF (HT >= HWDT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD+HWDT) < EL(KT,I)) THEN
    DLRHOT = ABS(RHO(KWD,I)-RHO(KTOP,I))
  ELSE IF (EL(KT,I) == ELWD) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KWD,I)-RHO(KT,I))*HWDT/(EL(KT,I)-ELWD)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KWD+1,KBOT

!** Density frequency

    HB    = ELWD-EL(K,I)
    RHOFB = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HB*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDB = (COEF*QWD(JWD)/RHOFB)**0.333333
    IF (HB >= HWDB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD-HWDB) > EL(KBOT+1,I)) THEN
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))
  ELSE
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))*HWDB/(ELWD-EL(KBOT+1,I))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                                                                             ! SW 1/24/05
  DO K=KTOP,KBOT
!    VNORM(K) = ABS(1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2)*BHR2(K,I)
 	   IF(K.GT.KWD)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
     VNORM(K) = 1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2
 	 IF(VNORM(K).GT.1.0) VNORM(K)=1.0                          !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                          !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,I)
     VSUM = VSUM+VNORM(K)
  END DO

! Outflows

  DO K=KTOP,KBOT
    tempest=tempest+t2(k,i)*(VNORM(K)/VSUM)*QWD(JWD)
  END DO
  if(qwd(jwd).gt.0.0)tempest=tempest/qwd(jwd)
  KTW(JWD) = KTOP
  KBW(JWD) = KBOT
  return


END SUBROUTINE WITHDRAWAL
