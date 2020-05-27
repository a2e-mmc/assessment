# -*- coding: utf-8 -*-
"""
Mesoscale Tendencies: Sensitivity to spatial and temporal averaging

Created on Tue Apr 19 04:53:49 2016

@author: cener
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime 
import utm
import glob
import sys
import netCDF4 
import pandas as pd

# Constants
g = 9.81    # [m s-2]tar -zxvf
P0 = 100000 # Reference pressure [Pa]
T0 = 300        # Reference temperature for perturbation temperature [K]
K  = 0.41       # von Karman constant
kappa = 0.2854  # Poisson constant (R/Cp)
R_air = 287.058  # Specific gas constant for dry air [J kg-1 K-1]
Cp_air = 1005   # Specific heat of air [J kg-1 K-1]
omega = 7.2921159e-5    # angular speed of the Earth [rad/s]

# Station coordinates
lat_s = 33.61054   # degrees N
lon_s = -102.05054    # degrees E
siteID = 'SWIFT'# identifier
fc = 2*omega*np.sin(lat_s*np.pi/180)
utmX_s, utmY_s, utm_zonenumber_s, utm_zoneletter_s = utm.from_latlon(lat_s,lon_s)

# Evaluation period
datefrom = datetime.datetime(2013,11,8,00,0,0)
dateto = datetime.datetime(2013,11,10,00,0,0)

# Vertical profiles 
t0 = datetime.datetime(2013,11,8,22,0,0,7320) 
#t0 = datetime.datetime(2013,11,9,6,0,0,0) 
#t0 = datetime.datetime(2013,11,8,18,0,0,0)

# ----------------------------------------------------------------------
# Load simulation data
# ----------------------------------------------------------------------
dirdata = '/projects/mmc/SWIFTRegion/8Nov2013/WRF_Rerun_Caro'
filenc = 'L0/SWIFT_all_divbyfc_w0_L0_tend_L0'
L = 0.0  # Length [m] of microscale domain Lx = Ly for spatial avaraging
dxmeso = 3000.0  # Horizontal resolution [m] of the mesoscale simulation 
dtmeso = 600.0     # temporal resolution [s]
Nav = (int(np.floor(L/dxmeso)) + 1 ) * 3

savenpz = 0     # Save tendencies to npz file 
savenc = 0      # Save tendencies to netCDF file 

# Resample using running time averages of different window sizes
#windowsize = np.array([1,31,61,121,181])   # number of (dtmeso) timesteps +1
#windowsize = np.array([19,7,1,13,4])  # take the output times into account
windowsize = np.array([4,19,1,7,13])
Nw = windowsize.shape[0]
#tav = dtmeso*(windowsize - 1)/60.0  # averaging time [min]

tav = (windowsize - 1)*10  # averaging time [min]

wrfoutfiles = sorted(glob.glob(dirdata+'/'+filenc+'_*'))
Nfiles = len(wrfoutfiles)

ifile = 0
cnt = 0
firstgoodfile = 1

while ifile < Nfiles:        
    f1 = netCDF4.Dataset(wrfoutfiles[ifile])
    print(wrfoutfiles[ifile])
    sys.stdout.flush()
    
    times1 = f1.variables['time'][:]
    dates1 = mdates.num2date(times1)
    # print(dates1)
    heights1 = f1.variables['z'][:]
    Nz = heights1.shape[0]
    
    idates = np.where(np.logical_and(times1 >= mdates.date2num(datefrom), 
                                     times1 <= mdates.date2num(dateto)))[0]                                 
    
    if idates.shape[0] > 0:
        times1 = times1[idates]
        Nt = times1.shape[0]
        U1 = pd.DataFrame(f1.variables['U'][:], index = dates1)
        V1 = pd.DataFrame(f1.variables['V'][:], index = dates1)
        W1 = pd.DataFrame(f1.variables['W'][:], index = dates1)
        Th1 = pd.DataFrame(f1.variables['Th'][:], index = dates1)# + 10.  # Th was calculated with T0 = 290 K instead of 300 K
        Ug1 = pd.DataFrame(f1.variables['Ug'][:], index = dates1)
        Vg1 = pd.DataFrame(f1.variables['Vg'][:], index = dates1)
        Uadv1 = pd.DataFrame(f1.variables['Uadv'][:], index = dates1)
        Vadv1 = pd.DataFrame(f1.variables['Vadv'][:], index = dates1)
        Thadv1 = pd.DataFrame(f1.variables['Thadv'][:], index = dates1)
        Ucor1 = pd.DataFrame(f1.variables['Ucor'][:], index = dates1)
        Vcor1 = pd.DataFrame(f1.variables['Vcor'][:], index = dates1)
        Uphys1 = pd.DataFrame(f1.variables['Uphys'][:], index = dates1)
        Vphys1 = pd.DataFrame(f1.variables['Vphys'][:], index = dates1)
        Utend1 = pd.DataFrame(f1.variables['Utend'][:], index = dates1)
        Vtend1 = pd.DataFrame(f1.variables['Vtend'][:], index = dates1)
        ust1 = pd.DataFrame(f1.variables['ust'][:], index = dates1)
        T21 = pd.DataFrame(f1.variables['T2'][:], index = dates1)
        TSK1 = pd.DataFrame(f1.variables['TSK'][:], index = dates1)
        HFX1 = pd.DataFrame(f1.variables['HFX'][:], index = dates1)
        LH1 = pd.DataFrame(f1.variables['LH'][:], index = dates1)
        Psfc1 = pd.DataFrame(f1.variables['Psfc'][:], index = dates1)

        rho1 = Psfc1/(R_air*T21)
        wt1 = HFX1/(rho1*Cp_air)
        beta01 = 1/T21
        L01 = -ust1**3/(K*g*beta01*wt1)   # Surface-layer Obukhov length [m]
        
        U1 = U1.loc[datefrom:dateto]
        V1 = V1.loc[datefrom:dateto]
        W1 = W1.loc[datefrom:dateto]
        Th1 = Th1.loc[datefrom:dateto]
        Ug1 = Ug1.loc[datefrom:dateto]
        Vg1 = Vg1.loc[datefrom:dateto]
        Uadv1 = Uadv1.loc[datefrom:dateto]
        Vadv1 = Vadv1.loc[datefrom:dateto]
        Thadv1 = Thadv1.loc[datefrom:dateto]
        Ucor1 = Ucor1.loc[datefrom:dateto]
        Vcor1 = Vcor1.loc[datefrom:dateto]
        Uphys1 = Uphys1.loc[datefrom:dateto]
        Vphys1 = Vphys1.loc[datefrom:dateto]
        Utend1 = Utend1.loc[datefrom:dateto]
        Vtend1 = Vtend1.loc[datefrom:dateto]
        ust1 = ust1.loc[datefrom:dateto]
        T21 = T21.loc[datefrom:dateto]
        TSK1 = TSK1.loc[datefrom:dateto]
        wt1 = wt1.loc[datefrom:dateto]
        LH1 = LH1.loc[datefrom:dateto]
        Psfc1 = Psfc1.loc[datefrom:dateto]
        L01 = L01.loc[datefrom:dateto]

        if firstgoodfile == 1:
            U = U1; V = V1; W = W1; Th = Th1; 
            Ug = Ug1; Vg = Vg1; Uadv = Uadv1; Vadv = Vadv1; 
            Utend = Utend1; Vtend = Vtend1; Thadv = Thadv1;
            Ucor = Ucor1; Vcor = Vcor1; Uphys = Uphys1; Vphys = Vphys1
            ust = ust1; T2 = T21; TSK = TSK1 
            wt = wt1; LH = LH1; Psfc = Psfc1; L0 = L01
            times = times1; heights = heights1 
            firstgoodfile = 0
            cnt = cnt + 1
        else:
            times = np.hstack((times,times1))
            heights = heights + heights1
            U = pd.concat([U,U1])
            V = pd.concat([V,V1])
            W = pd.concat([W,W1])
            Th = pd.concat([Th,Th1])
            Ug = pd.concat([Ug,Ug1])
            Vg = pd.concat([Vg,Vg1])
            Utend = pd.concat([Utend,Utend1])
            Vtend = pd.concat([Vtend,Vtend1])
            Uadv = pd.concat([Uadv,Uadv1])
            Vadv = pd.concat([Vadv,Vadv1])
            Ucor = pd.concat([Ucor,Ucor1])
            Vcor = pd.concat([Vcor,Vcor1])
            Uphys = pd.concat([Uphys,Uphys1])
            Vphys = pd.concat([Vphys,Vphys1])
            Thadv = pd.concat([Thadv,Thadv1])
            ust = pd.concat([ust,ust1])
            T2 = pd.concat([T2,T21])
            TSK = pd.concat([TSK,TSK1])
            wt = pd.concat([wt,wt1])
            LH = pd.concat([LH,LH1])
            Psfc = pd.concat([Psfc,Psfc1])
            L0 = pd.concat([L0,L01])
            
            cnt = cnt +1
            
    ifile = ifile + 1 
         
isort = np.argsort(times)      
times = times[isort]             

# Sort timeseries
U.sort_index(inplace=True); V.sort_index(inplace=True); 
W.sort_index(inplace=True); Th.sort_index(inplace=True) 
Utend.sort_index(inplace=True); Vtend.sort_index(inplace=True) 
Ug.sort_index(inplace=True); Vg.sort_index(inplace=True)
Uadv.sort_index(inplace=True); Vadv.sort_index(inplace=True)
Ucor.sort_index(inplace=True); Vcor.sort_index(inplace=True); 
Uphys.sort_index(inplace=True); Vphys.sort_index(inplace=True); 
Thadv.sort_index(inplace=True); ust.sort_index(inplace=True); 
T2.sort_index(inplace=True); TSK.sort_index(inplace=True); 
wt.sort_index(inplace=True); LH.sort_index(inplace=True); 
Psfc.sort_index(inplace=True); L0.sort_index(inplace=True)
heights = heights/cnt

# Fill initial NaN values in tendencies with next non-NaN values
Uadv = Uadv.bfill(); Vadv = Vadv.bfill()
Utend = Utend.bfill(); Vtend = Vtend.bfill()
Ucor = Ucor.bfill(); Vcor = Vcor.bfill()
Uphys = Uphys.bfill(); Vphys = Vphys.bfill()
Ug = Ug.bfill(); Vg = Vg.bfill()
Thadv = Thadv.bfill()

# Resample
Ugw = []; Vgw = []; Uadvw = []; Vadvw = []; Ucorw = []; Vcorw = []
Uphysw = []; Vphysw = []; Utendw = []; Vtendw = []; Thadvw = []
Uw = []; Vw = []; Ww = []; Thw = []
ustw = []; T2w = []; TSKw = []; wtw = []; LHw = []; Psfcw = []; 
L0w = []
    
for w in range(0,Nw):
    Uw.append(U.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(U))
    Vw.append(V.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(V))
    Ww.append(W.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(W))
    Thw.append(Th.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Th))
    Ugw.append(Ug.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Ug))
    Vgw.append(Vg.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Vg))
    Uadvw.append(Uadv.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Uadv))
    Vadvw.append(Vadv.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Vadv))
    Thadvw.append(Thadv.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Thadv))
    Ucorw.append(Ucor.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Ucor))
    Vcorw.append(Vcor.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Vcor))
    Uphysw.append(Uphys.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Uphys))
    Vphysw.append(Vphys.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Vphys))
    Utendw.append(Utend.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Utend))
    Vtendw.append(Vtend.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Vtend))
    ustw.append(ust.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(ust))
    T2w.append(T2.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(T2))
    TSKw.append(TSK.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(TSK))
    wtw.append(wt.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(wt))
    LHw.append(LH.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(LH))
    Psfcw.append(Psfc.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(Psfc))
    L0w.append(L0.rolling(window = windowsize[w], center = True, min_periods = 1).mean().fillna(L0))

# Save to text files
if savenpz == 1:
    for w in range(0,Nw):
        fileout = siteID + '_' + dirdata + '_w' + str(int(tav[w])) +'_L' + str(int(L)) 
        np.savez(fileout, U = Uw[w], V = Vw[w], W = Ww[w], Th = Thw[w],
                 Ug = Ugw[w], Vg = Vgw[w], 
                 Uadv = Uadvw[w], Vadv = Vadvw[w], Thadv = Thadvw[w],
                 Ucor = Ucorw[w], Vcor = Vcorw[w], 
                 Uphys = Uphysw[w], Vphys = Vphysw[w], 
                 Utend = Utendw[w], Vtend = Vtendw[w],
                 ust = ustw[w], T2 = T2w[w], TSK = TSKw[w], 
                 wt = wtw[w], LH = LHw[w], Psfc = Psfcw[w], L0 = L0w[w],
                 times = times, heights = heights)    
        print('Saving ' + fileout)
        sys.stdout.flush()
                 
if savenc == 1:
    for w in range(0,Nw):
#        fileout = siteID + '_' + dirdata + '_w' + str(int(tav[w])) +'_L' + str(int(L)) + '.nc'
        fileout = siteID + '_w' + str(int(tav[w])) +'_L' + str(int(L)) + '.nc'
        f = netCDF4.Dataset(fileout, 'w')
        f.history = f1.history +', tav = ' +  str(int(tav[w])) + ' min'
        
        f.createDimension('time', np.shape(times)[0])
        f.createDimension('z', Nz)
        f.createDimension('site', 1)
        
        lats = f.createVariable('lat', 'float', ('site',))
        lats.long_name = 'Site latitude'
        lats.units = 'degrees North'
        lats[:] = lat_s
        
        lons = f.createVariable('lon', 'float', ('site',))
        lons.long_name = 'Site longitude'
        lons.units = 'degrees East'
        lons[:] = lon_s
        
        fcs = f.createVariable('fc', 'float', ('site',))
        fcs.long_name = 'Coriolis parameter'
        fcs.units = 's-1'
        fcs[:] = fc
        
        timess = f.createVariable('time', 'float', ('time',))
        timess.long_name = 'Time'
        timess.units = 'Days since 001-01-01 00:00:00 UTC, plus one'
        timess[:] = times
        
        heightss = f.createVariable('z', 'float', ('z',))
        heightss.long_name = 'Height above ground level'
        heightss.units = 'm'
        heightss[:] = heights
        
        Us = f.createVariable('U', 'float', ('time','z',))
        Us.long_name = 'U velocity component' 
        Us.units = 'm s-1'
        Us[:] = Uw[w].as_matrix()
        
        Vs = f.createVariable('V', 'float', ('time','z',))
        Vs.long_name = 'V velocity component' 
        Vs.units = 'm s-1'
        Vs[:] = Vw[w].as_matrix()
        
        Ws = f.createVariable('W', 'float', ('time','z',))
        Ws.long_name = 'W velocity component' 
        Ws.units = 'm s-1'
        Ws[:] = Ww[w].as_matrix()
                
        Ths = f.createVariable('Th', 'float', ('time','z',))
        Ths.long_name = 'Potential temperature' 
        Ths.units = 'K'
        Ths[:] = Thw[w].as_matrix()
        
        Utends = f.createVariable('Utend', 'float', ('time','z',))
        Utends.long_name = 'Tendency momentum LHS term (divided by fc) U-component' 
        Utends.units = 'm s-1'
        Utends[:] = Utendw[w].as_matrix()
        
        Vtends = f.createVariable('Vtend', 'float', ('time','z',))
        Vtends.long_name = 'Tendency momentum LHS term (divided by fc) V-component' 
        Vtends.units = 'm s-1'
        Vtends[:] = Vtendw[w].as_matrix()
        
        Ugs = f.createVariable('Ug', 'float', ('time','z',))
        Ugs.long_name = 'Geostrophic wind U-component = Pressure gradient RHS term (divided by fc) V-component' 
        Ugs.units = 'm s-1'
        Ugs[:] = Ugw[w].as_matrix()
        
        Vgs = f.createVariable('Vg', 'float', ('time','z',))
        Vgs.long_name = 'Geostrophic wind V-component = - Pressure gradient RHS term (divided by fc) U-component' 
        Vgs.units = 'm s-1'
        Vgs[:] = Vgw[w].as_matrix()
        
        Uadvs = f.createVariable('Uadv', 'float', ('time','z',))
        Uadvs.long_name = 'Advective momentum RHS term (divided by fc) U-component' 
        Uadvs.units = 'm s-1'
        Uadvs[:] = Uadvw[w].as_matrix()
        
        Vadvs = f.createVariable('Vadv', 'float', ('time','z',))
        Vadvs.long_name = 'Advective momentum RHS term (divided by fc) V-component' 
        Vadvs.units = 'm s-1'
        Vadvs[:] = Vadvw[w].as_matrix()
        
        Ucors = f.createVariable('Ucor', 'float', ('time','z',))
        Ucors.long_name = 'Coriolis momentum RHS term (divided by fc) U-component = V velocity component' 
        Ucors.units = 'm s-1'
        Ucors[:] = Ucorw[w].as_matrix()
        
        Vcors = f.createVariable('Vcor', 'float', ('time','z',))
        Vcors.long_name = 'Coriolis momentum RHS term (divided by fc) V-component = -U velocity component' 
        Vcors.units = 'm s-1'
        Vcors[:] = Vcorw[w].as_matrix()
        
        Uphyss = f.createVariable('Uphys', 'float', ('time','z',))
        Uphyss.long_name = 'Vertical diffusion momentum RHS term (divided by fc) U-component' 
        Uphyss.units = 'm s-1'
        Uphyss[:] = Uphysw[w].as_matrix()
        
        Vphyss = f.createVariable('Vphys', 'float', ('time','z',))
        Vphyss.long_name = 'Vertical diffusion momentum RHS term (divided by fc) V-component' 
        Vphyss.units = 'm s-1'
        Vphyss[:] = Vphysw[w].as_matrix()
        
        Thadvs = f.createVariable('Thadv', 'float', ('time','z',))
        Thadvs.long_name = 'Potential temperature advective RHS term' 
        Thadvs.units = 'K s-1'
        Thadvs[:] = Thadvw[w].as_matrix()
        
        usts = f.createVariable('ust', 'float', ('time',))
        usts.long_name = 'Friction velocity at the surface' 
        usts.units = 'm s-1'
        usts[:] = ustw[w].as_matrix()
        
        T2s = f.createVariable('T2', 'float', ('time',))
        T2s.long_name = '2-m temperature' 
        T2s.units = 'K'
        T2s[:] = T2w[w].as_matrix()
        
        TSKs = f.createVariable('TSK', 'float', ('time',))
        TSKs.long_name = 'skin temperature' 
        TSKs.units = 'K'
        TSKs[:] = TSKw[w].as_matrix()
        
        wts = f.createVariable('wt', 'float', ('time',))
        wts.long_name = 'Kinematic Upward sensible heat flux at surface wt = HFX/(rho*Cp)' 
        wts.units = 'K m s-1'
        wts[:] = wtw[w].as_matrix()
        
        LHs = f.createVariable('LH', 'float', ('time',))
        LHs.long_name = 'Upward latent heat flux at surface' 
        LHs.units = 'W m-2'
        LHs[:] = LHw[w].as_matrix()    

        Psfcs = f.createVariable('Psfc', 'float', ('time',))
        Psfcs.long_name = 'Surface pressure' 
        Psfcs.units = 'Pa'
        Psfcs[:] = Psfcw[w].as_matrix()
        
        L0s = f.createVariable('L', 'float', ('time',))
        L0s.long_name = 'M-O length at the surface, based on ust, T2 and HFX' 
        L0s.units = 'm'
        L0s[:] = L0w[w].as_matrix()

        f.close()
        print('Saving' + fileout)
        sys.stdout.flush()
        
# Plot momentum budget
w = 2         # index to windowsize  
Lav = L/1000.0  # km

#Z = (Uw[w].T, Utendw[w].T, -Vgw[w].T, Uadvw[w].T, Ucorw[w].T, Uphysw[w].T)
#Zvarname = ('$U$','$U_{tend}$','$U_{pg}$','$U_{adv}$','$U_{cor}$','$U_{pbl}$')
#figname = siteID+'TimeHeight_'+"%.0f"%(Lav)+'km_'+"%.0f"%(tav[w])+'min.png'

#Z = (Vw[w].T, Vtendw[w-1].T, Ugw[w-1].T, Vadvw[w-1].T, Vcorw[w-1].T, Vphysw[w-1].T)
#Zvarname = ('$V$','$V_{tend}$','$V_{pg}$','$V_{adv}$','$V_{cor}$','$V_{pbl}$')
#figname = siteID+'_tend_'+"%.0f"%(Lav)+'km_'+"%.0f"%(tav[w])+'min_V-MB_tz.png'

# Wind speed:
ff = (Uw[w].T**2 + Vw[w].T**2)**0.5
fftend = (Utendw[w].T**2 + Vtendw[w].T**2)**0.5
ffg = ( (-Vgw[w].T)**2 + Ugw[w].T**2)**0.5
ffadv = (Uadvw[w].T**2 + Vadvw[w].T**2)**0.5
ffcor = (Ucorw[w].T**2 + Vcorw[w].T**2)**0.5
ffphys = (Uphysw[w].T**2 + Vphysw[w].T**2)**0.5

Z = (ff, fftend, ffg, ffadv, ffcor, ffphys)
Zvarname = ('$S$','$S_{tend}$','$S_{pg}$','$S_{adv}$','$S_{cor}$','$S_{pbl}$')
figname = siteID+'TimeHeight_Speed_'+"%.0f"%(Lav)+'km_'+"%.0f"%(tav[w])+'min.png'


# Direction for all of these:
#dd = 180 + np.arctan2(Uw[w].T,Vw[w].T)*180/np.pi
#ddtend = 180 + np.arctan2(Utendw[w].T,Vtendw[w].T)*180/np.pi
#ddg = 180 + np.arctan2(-Vgw[w].T,Ugw[w].T)*180/np.pi
#ddadv = 180 + np.arctan2(Uadvw[w].T,Vadvw[w].T)*180/np.pi
#ddcor = 180 + np.arctan2(Ucorw[w].T,Vcorw[w].T)*180/np.pi
#ddphys = 180 + np.arctan2(Uphysw[w].T,Vphysw[w].T)*180/np.pi

#Z = (dd, ddtend, ddg, ddadv, ddcor, ddphys)
#Zvarname = ('$D$','$D_{tend}$','$D_{pg}$','$D_{adv}$','$D_{cor}$','$D_{pbl}$')
#figname = siteID+'TimeHeight_Dir_'+"%.0f"%(Lav)+'km_'+"%.0f"%(tav[w])+'min.png'


# for Temperature advection == doesn't work
#Z = (Thw, Thadv, Thadvw)
#Zvarname = ('$Thw$','$Thadv$','$Thadvs$',)
#figname = siteID+'TimeHeight_Thadv_'+"%.0f"%(Lav)+'km_'+"%.0f"%(tav[w])+'min.png'

taxis_label = 'UTC time in hours since ' + datefrom.strftime('%Y-%m-%d %H:%M') + ', $L_{avg}$ = ' + "%i"%(Nav) + ' km, $t_{avg}$ = ' + "%.0f"%(tav[w]) + ' min'  
  
taxis_label = 'UTC time' 
      
zaxis = heights 
taxis = times
# Zlevels = np.linspace(0,360,361, endpoint=True)
Zlevels = np.linspace(0,21,22, endpoint=True)
zaxis_label = '$z$ [m]'
hoursFmt = mdates.DateFormatter('%H')

[X,Y]=np.meshgrid(taxis,zaxis)

zlim = 2000.0
#ticks = np.asarray(pd.date_range(datefrom,dateto,freq='6H'))
ticks = np.linspace(mdates.date2num(datefrom),mdates.date2num(dateto),7)
fig, ax = plt.subplots(2, 3, sharex='col', sharey='row')
Zcmap = plt.get_cmap('YlGnBu') # for dir: twilight # for wind: YlGnBu
for iax in range (0,6):
    ix,iy = np.unravel_index(iax,(2,3))
    #CS = ax[ix,iy].contour(X,Y,Z[iax], Zlevels, linewidths=0.5, colors='k')
    CF = ax[ix,iy].contourf(X,Y,Z[iax], Zlevels, cmap=Zcmap)
    ax[ix,iy].set_ylim([10, zlim])
    ax[ix,iy].set_yscale('linear')
    #ax[ix,iy].set_yscale('log')
    #ax[ix,iy].clabel(CS, Zlevels[0::4], inline=1, fontsize=8, fmt='%1.1f')
    #plt.colorbar(ticks=Zlevels)
    ax[ix,iy].xaxis.set_major_formatter(hoursFmt)    
    ax[ix,iy].set_title(Zvarname[iax])
    ax[ix,iy].set_xticks(ticks)

ax[0,0].set_ylabel(zaxis_label); ax[1,0].set_ylabel(zaxis_label)
ax[1,1].set_xlabel(taxis_label)

plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)

fig.subplots_adjust(right=0.92)
cbar_ax = fig.add_axes([0.95, 0.12, 0.025, 0.75])
cbar = fig.colorbar(CF, cax=cbar_ax)    
cbar.ax.set_xlabel('[m s$^{-1}$]',labelpad = -237, x = 1.5)
cbar.ax.xaxis.set_label_coords(1.05, 1.08)
plt.savefig(figname, dpi=300, bbox_inches='tight')

# Vertical profiles 
# t0 = datetime.datetime(2013,11,9,18,0,0)  # stable LLjet
#t0 = datetime.datetime(2006,7,1,19,0,0) 
w = np.array([0,1,2,3,4])

ZS = []; ZWD = []
for i in range(0,len(w)):
    U0 = Uw[w[i]].loc[t0].values
    V0 = Vw[w[i]].loc[t0].values
    S0 = (U0**2 + V0**2)**0.5
    WD0 = 180 + np.arctan2(U0,V0)*180/np.pi
    Upgf0 = -Vgw[w[i]].loc[t0].values
    Vpgf0 = Ugw[w[i]].loc[t0].values
    Spgf0 = (Upgf0**2 + Vpgf0**2)**0.5
    WDpgf0 = 180 + np.arctan2(Upgf0,Vpgf0)*180/np.pi
    Uadv0  = Uadvw[w[i]].loc[t0].values
    Vadv0  = Vadvw[w[i]].loc[t0].values
    Sadv0 = (Uadv0**2 + Vadv0**2)**0.5
    WDadv0 = 180 + np.arctan2(Uadv0,Vadv0)*180/np.pi
    Spgfadv0 = ((Upgf0+Uadv0)**2 + (Vpgf0+Vadv0)**2)**0.5
    WDpgfadv0 = 180 + np.arctan2(Upgf0+Uadv0,Vpgf0+Vadv0)*180/np.pi
#    WD0[WD0>180] = WD0[WD0>180] - 360
#    WDpgf0[WDpgf0>180] = WDpgf0[WDpgf0>180] - 360
#    WDadv0[WDadv0>180] = WDadv0[WDadv0>180] - 360
    WDpgfadv0[WDpgfadv0>180] = WDpgfadv0[WDpgfadv0>180] - 360
##    ZS.append((S0, Spgf0, Sadv0, Spgfadv0))
    ZS.append((S0, Spgf0, Sadv0))
#    ZWD.append((WD0, WDpgf0, WDadv0, WDpgfadv0))
ZSlabel = (('$S$ [$m s^{-1}$]','$S_{pg}$ [$m s^{-1}$]',
                  '$S_{adv}$ [$m s^{-1}$]'))
#                  '$S_{pg+adv}$ [$m s^{-1}$]'))
#ZWDlabel = (('$WD$ ['+u'\N{DEGREE SIGN}'+']', 
#                    '$WD_{pg}$ ['+u'\N{DEGREE SIGN}'+']', 
#                    '$WD_{adv}$ ['+u'\N{DEGREE SIGN}'+']', 
#                    '$WD_{pg+adv}$ ['+u'\N{DEGREE SIGN}'+']'))
    
figname = siteID+'Profiles_'+"%.0f"%(Lav)+'km_'+t0.strftime('%Y-%m-%d_%H:%M')+'.png'

linespec = ['b-','r-','k-','m-','g-','c-.'] 
fig,ax = plt.subplots(1, 3, sharey='row', figsize=(8,6))
fig.subplots_adjust(bottom=0.20)
Nticks = 6
##for iax in range (0,4):
for iax in range (0,3):
    for w in range(0,Nw):
        ax[iax].plot(ZS[w][iax][1:], heights[1:], 
                        linespec[w], linewidth = 1,
                        label = "%.0f"%(tav[w])+' min')
    ax[iax].set_xlabel(ZSlabel[iax],size=12)
    ax[iax].set_xlim([0, 12])
    ax[iax].set_xticks(np.linspace(0,20,Nticks))
    ax[iax].grid(which='major',color='k',linestyle='--')
 #   bx = ax[iax].twiny()
   # bx.plot(ZWD[i][iax][1:],heights[1:],color='grey',linestyle='-')
 #   bx.set_frame_on(True)
 #   bx.patch.set_visible(False)
 #   bx.xaxis.set_ticks_position('bottom')
 #   bx.xaxis.set_label_position('bottom')
 #   bx.spines['bottom'].set_position(('outward', 40))
#    bx.spines['bottom'].set_color('grey')
#    bx.tick_params(axis='x', colors='grey')
   # bx.set_xlabel(ZWDlabel[iax], color='grey',size=12)                       
 #   xlim = bx.get_xlim()
 #   bx.set_xticks(np.linspace(xlim[0],xlim[0]+np.ceil((xlim[1]-xlim[0])/(Nticks-1))*(Nticks-1),Nticks))
 #   ax[iax].set_ylim([0, zlim])

ax[0].set_ylabel('$z$ [$m$]',size=12)
plt.suptitle(t0.strftime('%Y-%m-%d %H:%M')+', $L_{avg}$ = ' + "%i"%(Nav) + ' km', size = 12); 
ax[2].legend(prop={'size':9})
plt.setp([a.get_yticklabels() for a in ax[1:]], visible=False)
plt.tight_layout()
fig.subplots_adjust(top=0.94)
plt.savefig(figname, dpi=300, bbox_inches='tight')
