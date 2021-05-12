# -*- coding: utf-8 -*-

#ENTRY SURFACE_TERMS (TSUR)

import numpy as np
from operator import add 
from cmath import pi
import math
from numpy import *
import pylab as plt
import pandas as pd
import os
import psutil


data = pd.read_excel('meteorology_1989_2008penideR.xlsx', sheet_name='Sheet1', index_col=0, header=0).iloc[:,[0,1,2,3,-1]]
data.columns = ['AirTemperature',' wind','RelativeHumidity', 'Cloud', 'WaterTemp']

Xval = data.iloc[:,[0,1,2,3,-1]].as_matrix()
Yval=Xval.transpose()

AirTemperature = Yval[0]
wind = Yval[1]
RelativeHumidity = Yval[2]
Cloud = Yval[3]
WaterTemp = Yval[4]


N = 35041
#7305
I = np.linspace(1,N,N)


RE = np.zeros(N)
RC = np.zeros(N)


def Latent_Heat_EVAP(AirTemperature, wind, RelativeHumidity, Cloud, WaterTemp):
    BOWEN_CONSTANT = 0.47
    VapourPressure_ea = 0.611*np.exp(17.3*AirTemperature/(AirTemperature+237.3))*RelativeHumidity/100
    TDEW = (np.log(VapourPressure_ea)+0.4926)/(0.0708-0.00421*np.log(VapourPressure_ea)) #ºC
    #TDEW = DewPt*1.8+32 # graus F
    #WaterTemp=WaterTemp*1.8+32
    #! Partial water vapor pressure of air (mm hg)
    #!  EA = EXP(2.3026*(9.5*TDEW(JW)/(TDEW(JW)+265.5)+0.6609))         ! SW 6/10/2011
    #!  IF (TDEW(JW) > 0.0) EA = EXP(2.3026*(7.5*TDEW(JW)/(TDEW(JW)+237.3)+0.6609))
 
    EA = np.exp(2.3026*(7.5*TDEW/(TDEW+237.3)+0.6609))  #

    #Partial water vapor pressure at the water surface

    if WaterTemp<0.0:
        ES = np.exp(2.3026*(9.5*WaterTemp/(WaterTemp+265.5)+0.6609))
    else:
        ES = np.exp(2.3026*(7.5*WaterTemp/(WaterTemp+237.3)+0.6609))


    #! Wind function
    
    TAIRV = (AirTemperature+273.0)/(1.0-0.378*EA/760.0) #ºK
    DTV   = (WaterTemp+273.0)/(1.0-0.378*ES/760.0)-TAIRV #ºK
    DTVL  =  0.0084*wind**3
    
    if DTV < DTVL:
        DTV=DTVL
        FW = (3.59*DTV**0.3333+4.26*wind)
    else:
        FW = 9.2+0.46*wind**2
    # Evaporative flux

    RE = FW*(ES-EA)

    # Conductive flux

    RC = FW*BOWEN_CONSTANT*(WaterTemp-AirTemperature)

    return RE
    


def Sensible_Heat_SENS(AirTemperature, wind, RelativeHumidity, Cloud, WaterTemp):
    BOWEN_CONSTANT = 0.47
    VapourPressure_ea = 0.611*np.exp(17.3*AirTemperature/(AirTemperature+237.3))*RelativeHumidity/100
    TDEW = (np.log(VapourPressure_ea)+0.4926)/(0.0708-0.00421*np.log(VapourPressure_ea)) #ºC
    #TDEW = DewPt*1.8+32 # graus F
    #WaterTemp=WaterTemp*1.8+32
    #! Partial water vapor pressure of air (mm hg)
    #!  EA = EXP(2.3026*(9.5*TDEW(JW)/(TDEW(JW)+265.5)+0.6609))         ! SW 6/10/2011
    #!  IF (TDEW(JW) > 0.0) EA = EXP(2.3026*(7.5*TDEW(JW)/(TDEW(JW)+237.3)+0.6609))
 
    EA = np.exp(2.3026*(7.5*TDEW/(TDEW+237.3)+0.6609))  #

    #Partial water vapor pressure at the water surface

    if WaterTemp<0.0:
        ES = np.exp(2.3026*(9.5*WaterTemp/(WaterTemp+265.5)+0.6609))
    else:
        ES = np.exp(2.3026*(7.5*WaterTemp/(WaterTemp+237.3)+0.6609))


    #! Wind function
    
    TAIRV = (AirTemperature+273.0)/(1.0-0.378*EA/760.0) #ºK
    DTV   = (WaterTemp+273.0)/(1.0-0.378*ES/760.0)-TAIRV #ºK
    DTVL  =  0.0084*wind**3
    
    if DTV < DTVL:
        DTV=DTVL
        FW = (3.59*DTV**0.3333+4.26*wind)
    else:
        FW = 9.2+0.46*wind**2
    # Evaporative flux

    RE = FW*(ES-EA)

    # Conductive flux

    RC = FW*BOWEN_CONSTANT*(WaterTemp-AirTemperature)

    return RC




for i in range(1,N+1):
    RE[i-1] = Latent_Heat_EVAP(AirTemperature[i-1], wind[i-1], RelativeHumidity[i-1], Cloud[i-1], WaterTemp[i-1])
    RC[i-1] = Sensible_Heat_SENS(AirTemperature[i-1], wind[i-1], RelativeHumidity[i-1], Cloud[i-1], WaterTemp[i-1])#RC[i-1] = SurfaceTerms(AirTemperature[i-1], wind[i-1], RelativeHumidity[i-1], Cloud[i-1], WaterTemp[i-1])
    



plot = True
if plot:
    plt.plot(I,RE) # 
    plt.show()
    plt.plot(I,RC) # 
    plt.show()



#Export results to Excel

np.savetxt('Latent_Heat_EVABpenideANN.out', RE, delimiter=',')

np.savetxt('Sensible_Heat_SENpenideANN.out', RC, delimiter=',')

