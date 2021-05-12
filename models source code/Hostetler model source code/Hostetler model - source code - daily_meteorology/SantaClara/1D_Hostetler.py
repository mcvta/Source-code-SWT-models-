'''
Created on 18 de Abr de 2012

@author: Hidraulica
'''

import numpy as np
from operator import add 
from cmath import pi
import math
from numpy import *
import pylab as plt
import pandas as pd
import os
import psutil
import time
start = time.time()

"Implicit Crank Nicolson scheme - Centered differences in space and time - Thomas algorithm" 
"Eddy diffusion following Henderson - During late summer convection during the cooling period was represented by the procedure described in Orlob and Selna (1970)"
"The model considers Subsurface heating by the absorption of penetrating solar radiation, in accordance with Beer's law"

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#INPUT - METEOROLOGY 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

data = pd.read_excel('meteorology_1989_2008.xlsx', index_col=0, header=0).iloc[:,[0,1,2,-1]]
data.columns = ['AirTemperature',' wind','RelativeHumidity','Cloud']

Xval = data.iloc[:,[0,1,2,-1]].as_matrix()
Yval=Xval.transpose()

AirTemperature = Yval[0]
wind = Yval[1]
RelativeHumidity = Yval[2]
Cloud = Yval[3]

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#INPUT - LAKE BATHYMETRY
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

data = pd.read_excel('bathymetrySantaClara.xlsx', index_col=0, header=0).iloc[:,[0,-1]]
data.columns = ['Depth','Areat0']

Xval = data.iloc[:,[0,-1]].as_matrix()
Yval=Xval.transpose()

Depth = Yval[0]
Areat0 = Yval[1]


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#INPUT - INITIAL CONDITIONS
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Simulation time period (days)
N=7305

#Latitude (degrees)
Lat=37.51

# Number of layers
Nlayers=37

# Lake initial water temperature C
Lake_initial_water_temperature = 12.0
InitialWaterTemperature = np.full((Nlayers),Lake_initial_water_temperature)

#Layer depth, m
DZ=2.0

#time interval, seconds, s
DT=86400
#Water density, kg/m3
DH2O=1000.0
#Water specific Heat j/kg C
Cp=4186

s=DT/(2*DZ**2.0)

#Ligth extinction coefficient. m-1
ExtCoef = 0.45

#Fraction of shortwave radiation incident at the water surface that is reflected at the water surface (albedo)
Albedo = 0.07

# Proportion of shortwave radiation that is absorved in the surface layer
Beta = 0.45

SHADING = 1.0 #0.85

#------------------------------------------------------------------------------------------------------------------------------------
# Tri Diagonal Matrix Algorithm (Thomas algorithm) solver (source:https://gist.github.com/cbellei/8ab3ab8551b8dfc8b081c518ccd9ada9)
#-------------------------------------------------------------------------------------------------------------------------------------

def TDMAsolver(a, b, c, d):
    nf = len(d)
    ac, bc, cc, dc = map(np.array, (a, b, c, d))
    for it in xrange(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
    xc = bc
    xc[-1] = dc[-1]/bc[-1]
    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]
    del bc, cc, dc

    return xc

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Clear Sky radiation (S.L.Dingman, 2012 - Physical Hydrology, 2nd Ed.)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def ClearSky(i, Lat, AirTemperature, RelativeHumidity, Cloud):
    Latitude_radians = (Lat/360)*2*pi    
    SolarConstant = 4.921 #MJ/m2hr
    DustAttenuation = 0.05
    OpticalAirMass = 3.6
    AirT = AirTemperature
    RH = RelativeHumidity
    SurfaceAlbedo = 0.4
    VapourPressure_ea = 0.611*np.exp(17.3*AirT/(AirT+237.3))*RH/100
    DewPt = (np.log(VapourPressure_ea)+0.4926)/(0.0708-0.00421*np.log(VapourPressure_ea))
    PrecipitableWater = 1.12*np.exp(0.0614*DewPt)
    Tsa = np.exp((-0.124-0.0207*PrecipitableWater)+(-0.0682-0.0248*PrecipitableWater)*OpticalAirMass)
    T = Tsa-DustAttenuation
    Ts = np.exp((-0.0363-0.0084*PrecipitableWater)+(-0.0572-0.0173*PrecipitableWater)*OpticalAirMass)
    gs = 1-Ts+DustAttenuation   
    DayAngle = 2*pi*(i-1)/365
    Declination_degrees = (180/pi)*(0.006918-0.399912*np.cos(DayAngle)+0.070257*np.sin(DayAngle)-0.006758*np.cos(2*DayAngle)+0.000907*np.sin(2*DayAngle)-0.002697*np.cos(3*DayAngle)+0.00148*np.sin(3*DayAngle))
    Declination_radians = (Declination_degrees/360)*2*pi
    Eccentricity = 1.00011+0.034221*np.cos(DayAngle)+0.00128*np.sin(DayAngle)+0.000719*np.cos(2*DayAngle)+0.000077*np.sin(2*DayAngle)
    Sunrise = -np.arccos(-np.tan(Latitude_radians)*np.tan(Declination_radians))/0.2618
    Sunset = np.arccos(-np.tan(Latitude_radians)*np.tan(Declination_radians))/0.2618
    DayLength = 2*Sunset
    Extraterrestrial = 2*SolarConstant*Eccentricity*(np.cos(Declination_radians)*np.cos(Latitude_radians)*np.sin(0.2618*Sunset)/0.2618+np.sin(Declination_radians)*np.sin(Latitude_radians)*Sunset)#MJ/m2hr
    Direct = Extraterrestrial*T
    Diffuse = 0.5*gs*Extraterrestrial
    Global = Direct+Diffuse
    Backscatter = 0.5*gs*SurfaceAlbedo*Global
    ClearSky = Global+Backscatter # MJ/m2day
    ClearSkyRadiaton = ClearSky*11.575*(1-0.0065*Cloud**2)
    
    return ClearSkyRadiaton

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EQUILILBRIUM TEMPERATURE
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def EquilibriumTemp(AirTemperature, RelativeHumidity, wind, ClearSkyRadiaton, WaterTemperature):
    AirT_C = AirTemperature
    AirT_F = AirT_C*1.8+32
    WaterTemp_C = WaterTemperature
    WindZ_ms = wind
    WindZ_mph = WindZ_ms*2.23714
    Wind2_mph = WindZ_mph*np.log(2/0.003)/np.log(10/0.003)
    SCR_Wm2=ClearSkyRadiaton
    SCR_BTU_FT2_DAY = SCR_Wm2*7.60796*SHADING
    RH = RelativeHumidity
    VapourPressure_ea = 0.611*np.exp(17.3*AirT_C/(AirT_C+237.3))*RH/100
    DewPt = (np.log(VapourPressure_ea)+0.4926)/(0.0708-0.00421*np.log(VapourPressure_ea))
    DewPt_F = DewPt*1.8+32
    Aconv = 7.60796
    Bconv = 1.520411 #IF(CFW =2,1.520411,3.401062)
    AFW = 9.2
    BFW = 0.46
    CFW = 2
    
    ET = DewPt_F
    TSTAR = (ET+DewPt_F)*0.5
    BETA  =  0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
    FW = Aconv*AFW+Bconv*BFW*Wind2_mph**CFW
    CSHE = 15.7+(0.26+BETA)*FW
    RA = 3.1872E-08*(AirT_F+459.67)**4
    ETP = (SCR_BTU_FT2_DAY+RA-1801.0)/CSHE+(CSHE-15.7)*(0.26*AirT_F+BETA*DewPt_F)/(CSHE*(0.26+BETA))

    j = 0
    while (abs(ETP-ET) > 0.05 or j<10):
        ET = ETP
        TSTAR = (ET+DewPt_F)*0.5
        BETA = 0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
        CSHE = 15.7+(0.26+BETA)*FW
        ETP = (SCR_BTU_FT2_DAY+RA-1801.0)/CSHE+(CSHE-15.7)*(0.26*AirT_F+BETA*DewPt_F)/(CSHE*(0.26+BETA))
        EquilibriumTemperature = (ET-32)*5/9
        SurfaceHeatExchange = CSHE*0.23659
        Qn = (EquilibriumTemperature-WaterTemp_C)*SurfaceHeatExchange
        j+=1
        
    return Qn


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEPTH OF THE THERMOCLINE (source:https://rdrr.io/cran/rLakeAnalyzer/src/R/thermo.depth.R)
# This function calculates the location of the thermocline from a temperature profile.  
# It uses a special technique to estimate where the thermocline lies even between two temperature measurement depths, giving a potentially finer-scale estimate than usual techniques.
# @param wtr a numeric vector of water temperature in degrees C
# @param depths a numeric vector corresponding to the depths (in m) of the wtr measurements
# @param Smin Optional paramter defining minimum density gradient for thermocline
# @param seasonal a logical value indicating whether the seasonal thermocline should be returned. This is fed to thermo.depth, which is used as the starting point.  The seasonal thermocline is defined as the deepest density gradient found in the profile. If \code{FALSE}, the depth of the maximum density gradient is used as the starting point.
# @param index Boolean value indicated if index of the thermocline depth, instead of the depth value, should be returned.
# @param mixed.cutoff A cutoff (deg C) where below this threshold, thermo.depth and meta.depths are not calculated (NaN is returned). Defaults to 1 deg C.
# @return Depth of thermocline. If no thermocline found, value is NaN.
#----------------------------------------------------------------------------------------------------------------------------------------

def find_peaks(data_in, thresh = 0.0): #Finds the local peaks in a vector. Checks the optionally supplied threshold for minimum height.
    peaks = [thresh - 1]
    for trip in zip(*[data_in[n:] for n in range(3)]):
        if max(trip) == trip[1]:
            peaks.append(trip[1])
        else:
            peaks.append(thresh - 1)
    return [index for index, peak in enumerate(peaks) if peak > thresh]


def ThermoTemp (wtr, depths, Smin=0.1, seasonal=True, index=False, mixed_cutoff=1):

    rhoVar = np.zeros(len(depths))
       
    for i in range(len(depths)):
        rhoVar[i]=(1000*(1-(wtr[i]+288.9414)*(wtr[i]-3.9863)**2/(508929.2*(wtr[i]+68.12963))))
        
    #Calculate the first derivate of density
    drho_dz = np.zeros(len(depths)-1)
    for i in range(len(depths)-1):
        drho_dz[i] = ( rhoVar[i+1]-rhoVar[i] )/( depths[i+1] - depths[i] )
        
    #look for two distinct maximum slopes, lower one assumed to be seasonal 
    thermoInd = np.argmax(drho_dz)
    mDrhoZ = drho_dz[thermoInd]
    thermoD = np.mean(depths[(thermoInd):(thermoInd+2)])
    #                            
    dRhoPerc = 0.15
    numDepths = len(depths)

    if thermoInd > 0 and thermoInd < (numDepths):

        Sdn = -(depths[thermoInd+1] - depths[thermoInd])/(drho_dz[thermoInd+1] - drho_dz[thermoInd])
        Sup = (depths[thermoInd]-depths[thermoInd-1])/(drho_dz[thermoInd]-drho_dz[thermoInd-1])

        upD  = depths[thermoInd]
        dnD  = depths[thermoInd+1]
 
        if np.isinf(Sup) or np.isinf(Sdn):
            thermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup))
      
    dRhoCut = np.max([dRhoPerc*mDrhoZ,Smin])
    locs = find_peaks(drho_dz, dRhoCut)
    pks = drho_dz[locs]
    
   
    if (len(pks)==0):
        SthermoD = thermoD
        SthermoInd = thermoInd
        
    else:
        mDrhoZ = pks[len(pks)-1]
        SthermoInd = locs[len(pks)-1] 
        if (SthermoInd > (thermoInd + 1)):
            SthermoD = np.mean(depths[(SthermoInd):(SthermoInd+2)])
            
            if (SthermoInd > 0 and SthermoInd < (numDepths)):
                Sdn = -(depths[SthermoInd+1] - depths[SthermoInd])/(drho_dz[SthermoInd+1] - drho_dz[SthermoInd])
                Sup = (depths[SthermoInd] - depths[SthermoInd-1])/(drho_dz[SthermoInd] - drho_dz[SthermoInd-1])
                upD  = depths[SthermoInd]
                dnD  = depths[SthermoInd+1]
                
                if np.isinf(Sup) or np.isinf(Sdn):
                    SthermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup))      
        else:
            SthermoD = thermoD
            SthermoInd = thermoInd
                    
    if SthermoD < thermoD:
        SthermoD = thermoD
        SthermoInd = thermoInd
         
#Ok, which output was requested. Index or value, seasonal or non-seasonal

    if(index):
        if(seasonal):
            return(SthermoInd)
        else:
             return(thermoInd)
    
    else:
        if(seasonal):
            return(SthermoD)
        else:
            return(thermoD)


#-------------------------------------------------------------------------------------------------------------------------------------
#EDDY DIFFUSIVITY
#-------------------------------------------------------------------------------------------------------------------------------------
def Keddy(TempLake,Depth,wind,Areat0):
    #Latitude radianos 0.66
    Lat=38.0
    Depth[0] = 0.5
    if wind == 0.0:
        w=(1.2*10**-3)*0.1
        klat=0.51*math.sin(Lat)/0.1
    else:
        w=(1.2*10**-3)*wind
        klat=0.51*math.sin(Lat)/wind
    
    #water density
    Dens=np.zeros(len(TempLake))
    for i in range(len(TempLake)-1):
        Dens[i]=(1000*(1-(TempLake[i]+288.9414)/(508929.2*(TempLake[i]+68.12963))))*(TempLake[i]-3.9863)**2
   
    #Brunt-Vaisala frequency 
    N=np.zeros(len(TempLake))
    for i in range(len(TempLake)-1):
        N[i]=-(9.81/Dens[i])*(Dens[i+1]-Dens[i])/(Depth[i+1]-Depth[i])
    
        #iniciando Richardson number
    Ri=np.zeros(len(TempLake))
    secondterm=np.zeros(len(TempLake))
    for i in range(len(TempLake)-1):
        secondterm [i] = ((w**2.0)*math.exp(-2*klat*Depth[i]))
        if secondterm [i] == 0.0:
            Ri[i]=0.0
        else:
            Ri[i]=(-1+(1+40*(N[i]**2.0)*(0.4**2.0)*Depth[i]**2/(w**2*math.exp(-2*klat*Depth[i])))**0.5)/20
    
    KH =np.zeros(len(TempLake))
    km = 1.0*10**-6 
    for i in range(len(TempLake)-1):
        KH[i]=(0.4*w*Depth[i])*math.exp(-klat*Depth[i])*((1+37*Ri[i]**2)**-1) + km

        KH[len(TempLake)-1]=KH[len(TempLake)-2]
    Kdiffusivity=[a*b for a,b in zip(Areat0,KH)]
    
    return Kdiffusivity

#-------------------------------------------------------------------------------------------------------------------------------------
#FIX PROFILE
#-------------------------------------------------------------------------------------------------------------------------------------

def fix_profile(Nlayers,TempLake, area):
    var1=range(Nlayers)

    for i in var1:
        for j in range(Nlayers-1):
            if TempLake[i]>TempLake[j]:
                change_profile(TempLake,i,j,area,Nlayers)
                break
    return TempLake         
 
def change_profile(TempLake,i,inst,area,Nlayers):
    jmin=inst-1
    if jmin<0:
        jmin=0

    jmax=i
    if jmax>Nlayers:
        jmax=Nlayers

    A=0.0
    B=0.0
    for j in range(jmin,jmax+1):
        A+=TempLake[j]*area[j]
        B+=area[j]
        Tconst=(A/B)
    
    for j in range(jmin,jmax+1):
        TempLake[j]=Tconst
#-------------------------------------------------------------------------------------------------------------------------------------
#Run Clear Sky Radiation
#-------------------------------------------------------------------------------------------------------------------------------------


ClearSkyRadiation = np.zeros(N)
I = np.linspace(1,N,N)

for i in xrange(N):
    ClearSkyRadiation[i] = ClearSky(i, Lat, AirTemperature[i], RelativeHumidity[i],Cloud[i])
    
#-------------------------------------------------------------------------------------------------------------------------------------
# Subsurface heating by the absorption of penetrating solar radiation, in accordance with Beer's law
#-------------------------------------------------------------------------------------------------------------------------------------


def SubRad (Albedo,Beta,Cp,Depth,DH2O,DT,ExtCoef,Radiation,Nlayers):
    SubRad=np.zeros(Nlayers)
    RadSN = (1-Albedo)*Radiation
    Variable1 = 1/(Cp*DH2O)     
    
    for i in range(Nlayers-2):
        SubRad[0] = (1-Beta)*RadSN*math.exp(-ExtCoef*0.0)
        SubRad[i+1] = (1-Beta)*SubRad[i]*math.exp(-ExtCoef*Depth[i+1])
    return SubRad*Variable1*DT


#-------------------------------------------------------------------------------------------------------------------------------------
#Crank Nicolson scheme
#-------------------------------------------------------------------------------------------------------------------------------------

#Initializing coefficients np.arrays
a=np.zeros((N,Nlayers-1)) #40
b=np.zeros((N,Nlayers)) #41
c=np.zeros((N,Nlayers-1)) 
d=np.zeros((N,Nlayers))
alfa = np.zeros(Nlayers)#41

# Radiation array
Radiation = np.zeros(N)

# SubRadiation array
SubRadiation=np.zeros((N+1,Nlayers))

# Water temperature array
TempLake=np.zeros((N+1,Nlayers))
TempLake[0]=InitialWaterTemperature

# K=np.zeros((N+1,Nlayers))
K=np.full((N+1,Nlayers), 0.0000)

# Thermocline_Depth array
Thermocline_Depth = np.zeros(N)


# MAIN LOOP

for i in xrange(Nlayers):
    alfa[i]=1.0/(Areat0[i]+Areat0[i])
    

for i in xrange(1,N+1):
    I=np.linspace(1,N,N)
    Radiation[i-1] = EquilibriumTemp(AirTemperature[i-1], RelativeHumidity[i-1], wind[i-1], ClearSkyRadiation[i-1], TempLake[i-1][0])*Beta
    SubRadiation[i-1] = SubRad(Albedo,Beta,Cp,Depth,DH2O,DT,ExtCoef,Radiation[i-1],Nlayers)
    #Thermocline_Depth[i-1] = ThermoTemp(TempLake[i-1], Depth, Smin=0.1, seasonal=True, index=False, mixed_cutoff=1)
    K[i-1]=Keddy(TempLake[i-1],Depth,wind[i-1],Areat0)
    K[i]=Keddy(TempLake[i-1],Depth,wind[i-1],Areat0)
    
    #b0, c0, d0
    b[i-1][0]=1+alfa[0]*s*(0.0+2.0*K[i][0]+K[i][1])
    c[i-1][0]=-alfa[0]*s*(0.0+2.0*K[i][0]+K[i][1])
    d[i-1][0]=((1.0-alfa[0]*s*(0.0+2.0*K[i-1][0]+K[i-1][1]))*TempLake[i-1][0])+(alfa[0]*s*(0.0+2.0*K[i-1][0]+K[i-1][1])*TempLake[i-1][1])+alfa[0]*s*(0.0+K[i-1][0])*((2.0*Radiation[i-1]*DZ)/(Cp*DH2O*(K[i-1][0]/Areat0[0])))+alfa[0]*s*(0.0+K[i][0])*((2.0*Radiation[i-1]*DZ)/(DH2O*Cp*(K[i][0]/Areat0[0])))
  
    #Bottom
    #an, bn, dn
    a[i-1][-1]=-alfa[-1]*s*(K[i][-2]+2.0*K[i][-1]+0.0)
    b[i-1][-1]=1+alfa[-1]*s*(K[i][-2]+2.0*K[i][-1]+0.0)
    d[i-1][-1]=alfa[-1]*s*TempLake[i-1][-2]*(K[i-1][-2]+2.0*K[i-1][-1]+0.0)+(1-alfa[-1]*s*(K[i-1][-2]+2*K[i-1][-1]+0.0))*TempLake[i-1][-1] # 0,1,2,3,4,5,6,7,8,9 
    
    
    #ai, bi, ci, di
    for j in xrange(1,Nlayers-1):# 40
        a[i-1][j-1]=-alfa[j]*s*(K[i][j-1]+K[i][j]) # 40
        b[i-1][j]=1.0+alfa[j]*s*(K[i][j-1]+2.0*K[i][j]+K[i][j+1])
        c[i-1][j]=-alfa[j]*s*(K[i][j]+K[i][j+1])
        d[i-1][j]=alfa[j]*s*(K[i-1][j]+K[i-1][j-1])*TempLake[i-1][j-1]+alfa[j]*s*(K[i-1][j+1]+K[i-1][j])*TempLake[i-1][j+1]+(1-alfa[j]*s*(K[i-1][j-1]+2.0*K[i-1][j]+K[i-1][j+1]))*TempLake[i-1][j]

        
    TempTDMA=list(map(add, TDMAsolver(a[i-1], b[i-1], c[i-1], d[i-1]), SubRadiation[i-1]))
    # Real surface layer water temperature

    
    # Fix water temperature profile
    TempTDMA_Fix=fix_profile(Nlayers,TempTDMA,Areat0)
    TempLake[i]=TempTDMA_Fix
    Tfict=(alfa[0]*s*(0.0+K[i][0])*(2.0*Radiation[i-1]*DZ)/(DH2O*Cp*(K[i][0]/Areat0[0])))+TempLake[i][1]
    TempLake[i][0]=(TempLake[i][0]+Tfict)/2.0
    
    if  TempLake[i][0]<TempLake[i][1]:
        TempLake[i][0]=TempLake[i][1]
    else:
        TempLake[i][0]=TempLake[i][0]

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#OUTPUT TO EXCEL
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

export = True
if export:
    Kout=np.full((N+1,Nlayers), 0.0000)

    for i in xrange(1,N+1):
        Kout[i] = [a/b for a,b in zip(K[i],Areat0)]

    np.savetxt('Kdifussivity.out', Kout, delimiter=',')
np.savetxt('TempLake.out', TempLake, delimiter=',')
np.savetxt('Radiation.out', Radiation, delimiter=',')

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#PLOT
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

plot = False
if plot:
    plt.plot(I,ClearSkyRadiation) # ClearSkyRadiation
    plt.show()
    plt.plot(I,Radiation) # EquilibriumTemp
    plt.show()  
    plt.plot(I,Thermocline_Depth) # Thermocline Depth
    plt.show()  


# MAIN PLOT
plot = True
if plot:
    Z=TempLake 
    X,Y=np.meshgrid(range(Z.shape[0]+1),range(Z.shape[1]+1)) 
    plt.figure(1)
    plt.subplot(211)
    im = plt.pcolormesh(X,Y,Z.transpose(), cmap='jet') 
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.colorbar(im, orientation='horizontal') 
    plt.ylabel('Profundidade, m')
    plt.subplot(212)
    Zout=Z.transpose()
    plt.plot(np.linspace(1,N+1,N+1), Zout[0], lw=1)
    plt.xlabel('Dias')
    plt.ylabel('Temperatura - Z[0], C')
    plt.show()

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MEMORY USED BY PROCESS
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

process = psutil.Process(os.getpid())
print(process.memory_info().rss)*9.53674*10**-7, "megabytes" # in Megabytes 


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RUNNING TIME
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

end = time.time()
print(end - start),"seconds" # in seconds
