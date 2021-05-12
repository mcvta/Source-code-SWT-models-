!***********************************************************************************************************************************
!**                                            F U N C T I O N   D E N S I T Y                                                    **
!***********************************************************************************************************************************

REAL(R8) FUNCTION DENSITY (T,TDS,SS)
  use PREC
  USE LOGICC, ONLY: SUSP_SOLIDS, FRESH_WATER, SALT_WATER; USE GLOBAL, ONLY:JW
  REAL(R8) :: T,TDS,SS
                       DENSITY = ((((6.536332E-9*T-1.120083E-6)*T+1.001685E-4)*T-9.09529E-3)*T+6.793952E-2)*T+0.842594
  IF (SUSP_SOLIDS)     DENSITY = DENSITY+6.2E-4*SS
  IF (FRESH_WATER(JW)) DENSITY = DENSITY+TDS*((4.99E-8*T-3.87E-6)*T+8.221E-4)
  IF (SALT_WATER(JW))  DENSITY = DENSITY+TDS*((((5.3875E-9*T-8.2467E-7)*T+7.6438E-5)*T-4.0899E-3)*T+0.824493)                      &
                                 +((-1.6546E-6*T+1.0227E-4)*T-5.72466E-3)*TDS**1.5+4.8314E-4*TDS*TDS
  DENSITY = DENSITY+999.0
END FUNCTION DENSITY
