"""
Helper functions specifically for the
profile asimilation coupled with obs study
"""

import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

# --------------------------------------------------------------
# manually add a2e-mmc repos to PYTHONPATH if needed
module_path = os.path.join(os.environ['HOME'],'tools','a2e-mmc')
if module_path not in sys.path:
    sys.path.append(module_path)
# --------------------------------------------------------------


# ---------------------------------------
# Dictionaries to define style of figures
# ---------------------------------------

sim_style={
    'case R':{'linestyle':'-','color':'tab:gray','linewidth':2},
    'case A0':{'linestyle':'-','color':'tab:orange','linewidth':2},
    'case A1':{'linestyle':'-','color':'tab:red','linewidth':2},
    'case An':{'linestyle':'-','color':'tab:green','linewidth':2},
    'case B0':{'linestyle':'-','color':'tab:purple','linewidth':2},
    'case B1':{'linestyle':'-','color':'tab:brown','linewidth':2},
    'case Bn':{'linestyle':'-','color':'tab:pink','linewidth':2},
    'case Ca':{'linestyle':'-','color':'tab:olive','linewidth':2},
    'case Cb':{'linestyle':'-','color':'tab:cyan','linewidth':2},
    'forcing data':{'linestyle':'--','color':'black','linewidth':2},
    'WRF mesoscale':{'linestyle':':','color':'tab:blue','linewidth':2},
}

my_fieldlabels={
    'wspd':r'$S$ [m s$^{-1}$]',
    'wdir':r'$D$ [$\mathrm{^\circ}$]',
    'TKE':r'$k$ [m$^{2}$ s$^{-2}$]',
    'TKEf':r'$k$ [m$^{2}$ s$^{-2}$]',
    'thetaw':r'$\overline{w^\prime \theta^\prime}$ [Km s$^{-1}$]',
    'tauz':r'$\tau_z$ [m$^{2}$ s$^{-2}$]',
}

# ------------------------------------------------------
# Functions to calculate certain quantities of interests
# ------------------------------------------------------

def calc_QOIs(df):
    """
    Calculate derived quantities
    """
    from mmctools.mmctools.helper_functions import calc_wind

    # Calculate horizontal wind speed and wind direction
    if not 'wspd' in df:
        df['wspd'],df['wdir'] = calc_wind(df)

    # If not present in dataset, set thetav = theta
    # E.g., - microscale simulations do not account for moisture
    #       - thetav was not extracted from WRF, and ref simulation was driven based on theta
    if not 'thetav' in df:
        df['thetav'] = df['theta']

    # Calculate the temperature gradient
    temperature_gradient_scale = 1000 # Used to convert temperature gradients from SI units K/m to K/km
    df['thetav_gradient'] = temperature_gradient_scale * calc_grad(df['thetav'])

    # Calculate some turbulence statistics if dataset contains turbulent stresses
    try:
        calc_turb_stats(df)
    except KeyError:
        pass

def calc_turb_stats(df):
    """
    Calculate turbulent statistics like friction velocity,
    turbulent kinetic energy, and turbulent intensity.
    """
    if not 'u*' in df: df['u*'] = (df['uw']**2 + df['vw']**2)**0.25
    if not 'TKE' in df: df['TKE'] = 0.5*(df['uu'] + df['vv'] + df['ww'])
    if not 'TI' in df:
        ang = np.arctan2(df['v'],df['u'])
        df['TI'] = df['uu']*np.cos(ang)**2 + 2*df['uv']*np.sin(ang)*np.cos(ang) + df['vv']*np.sin(ang)**2
        df['TI'] = np.sqrt(df['TI']) / df['wspd']
    
    if not 'tauz' in df: df['tauz'] = np.sqrt(df['uw']**2+df['vw']**2)

    # Calculate filtered TKE (30min) to filter out large peaks
    df['TKEf'] = df['TKE'].unstack().rolling(3,center=True).median().stack()


def calc_grad(df):
    """
    Calculate vertical gradient of specified field
    """
    from scipy.interpolate import interp1d

    # calculate at midpoints
    tv = df.unstack(level='datetime')
    dz = pd.Series(tv.reset_index()['height'].diff().values, index=tv.index) # dz and tv should have matching indices
    tvgrad = tv.diff().divide(dz,axis=0)

    # interpolate from midpoints back to original heights, for compatibility with original dataframe
    zorig = dz.index.values
    zmid = dz.index.values - dz/2
    tvgrad = tvgrad.set_index(zmid) # first index and row are NaN
    interpfun = interp1d(tvgrad.index, tvgrad, axis=0, bounds_error=False, fill_value=np.nan)
    for zi in zorig[:-1]:
        tvgrad.loc[zi] = interpfun(zi)
    zmid = zmid.values
    tvgrad.loc[zorig[0]] = tvgrad.loc[zmid[1]] # nearest values: second row (first row is NaN)
    tvgrad.loc[zorig[-1]] = tvgrad.loc[zmid[-1]] # nearest values: last row
    tvgrad = tvgrad.loc[zorig]
    
    tvgrad.index.name = 'height'    
    tvgrad = tvgrad.stack().reorder_levels(order=['datetime','height']).sort_index()
    return tvgrad

def calc_ablh(df,blending_height=150):
    """
    Calculate convective/stable boundary-layer heights
    """
    from mmctools.mmctools.helper_functions import estimate_ABL_height

    ablh = estimate_ABL_height(Tw=df['thetaw'])
    sblh = estimate_ABL_height(uw=np.sqrt(df['uw']**2 + df['vw']**2))
    ablh.loc[ablh < blending_height] = sblh[ablh < blending_height]

    # Apply rolling average (10-min data ==> 6 samples corresponds to 1h)
    return ablh.rolling(6, center=True).median()

def calc_surface_QOIs(df,zref,zmax):
    """
    Calculate derived quantities in the surface layer
    """
    # Store thetav at zref
    df_sl = sample(df,[zref,]).xs(zref,level='height')[['thetav']]
    # Calculate friction velocity
    df_sl['u*'] = calc_surface_friction_velocity(df['u*'],zmax)
    # Calculate surface heat flux
    df_sl['thetaw'] = calc_surface_heat_flux(df['thetaw'],zmax)
    # Store auxiliary thetaw field for zoom plot
    df_sl['thetaw_zoom'] = df_sl['thetaw']

    #Additional filtering (30min) as otherwise the signal contains some large peaks
    df_sl = df_sl.rolling(3, center=True).median()
    return df_sl

def calc_surface_friction_velocity(ustar,zmax):
    """
    Calculate surface friction velocity as the mamimum value below zmax
    """
    ustar = ustar.loc[ustar.index.get_level_values(1) <= zmax].unstack()
    return ustar.max(axis=1)

def calc_surface_heat_flux(thetaw,zmax):
    """
    Calculate surface heat flux as the min or max value below zmax
    (depending on the sign in the first grid cell)
    """
    # Store value in the first grid cell
    thetaw_s = thetaw.xs(thetaw.index.get_level_values(1)[0],level='height')
    # Limit to data below zmax
    thetaw = thetaw.loc[thetaw.index.get_level_values(1) <= zmax].unstack()
    # Normalise by value of first grid cell
    thetaw_norm = thetaw.divide(thetaw_s, axis=0)
    # Take max and rescale (if first grid cell is negative, this turns into a minimum function)
    return thetaw_s*thetaw_norm.max(axis=1)

def sample(df,newheights):
    """
    Interpolate to a new height
    """
    # First create a combined index including the original indices and the desired sampling times
    raw_index = df.index.levels[1]
    combined_index = raw_index.union(pd.Index(newheights,name='height'))
    # Resample to combined index to allow interpolation
    df = df.unstack(level=0).reindex(combined_index).interpolate(method='index').stack()
    return df.reorder_levels(['datetime','height']).sort_index()

# --------------------------------
# Functions to load reference data
# --------------------------------

def load_wrf_mesoscale_data(filename):
    """
    Load and preprocess netcdf file with WRF data generated in wrf_to_sowfa.ipynb
    """
    from mmctools.mmctools.helper_functions import calc_wind

    # Read netcdf file as xarray and convert to pandas dataframe
    wrf = xr.open_dataset(filename).to_dataframe()
    # Convert to standard names
    wrf.rename({'U':'u','V':'v','W':'w','UST':'u*','wt':'thetaw'},axis='columns',inplace=True)
    # Rename Time index
    wrf.rename_axis(index={'Time':'datetime'},inplace=True)
    # Reorder indices
    wrf = wrf.reorder_levels(['datetime','height']).sort_index()
    # Keep columns we are interested in
    wrf = wrf[['u','v','w','theta','u*','thetaw']]

    # Calculate quantities of interest
    calc_QOIs(wrf)

    return wrf

def load_TTU_tower_10min_reference_data(filename):
    """
    Load and preprocess TTU tower data (10-min statistics)
    Datasets are generated in process_TTU_tower.ipynb
    """
    from mmctools.mmctools.helper_functions import T_to_Tv, theta

    # Read csv file into pandas dataframe
    df = pd.read_csv(filename,parse_dates=True,index_col=['datetime','height'])
    # Calculate virtual temperature and virtual potential temperature
    df['Tv'] = T_to_Tv(df['T'],RH=df['RH'],p=df['p'])
    df['thetav'] = theta(df['Tv'], p=df['p'])

    # Calculate quantities of interest
    calc_QOIs(df)

    return df

def load_TTU_tower_1Hz_reference_data(filename):
    """
    Load and preprocess TTU tower data (1 Hz data)
    """
    from mmctools.mmctools.helper_functions import theta, calc_wind

    # Read csv file into pandas dataframe
    df = pd.read_csv(filename,parse_dates=True,index_col=['datetime','height'])
    # Calculate potential temperature
    df['theta'] = theta(df['Ts'], p=df['p'])
    # Calculate horizontal wind speed and wind direction
    df['wspd'],df['wdir'] = calc_wind(df)

    return df

def load_TTU_radar_reference_data(filename):
    """
    Load and preprocess TTU radar data
    Dataset generated in process_TTU_radar.ipynb
    """
    combined_data = pd.read_csv(filename,parse_dates=True,index_col=['datetime','height'])
    radar_scan0 = combined_data.loc[combined_data['scan_type']==0].copy()
    radar_scan1 = combined_data.loc[combined_data['scan_type']==1].copy()
    return radar_scan0, radar_scan1

def load_TTU_RASS_reference_data(filename):
    """
    Load and preprocess TTU RASS data
    Dataset generated in temperature_profile_reconstruction.ipynb
    """
    combined_data = pd.read_csv(filename,parse_dates=True,index_col=['datetime','height'])
    return combined_data.loc[combined_data['type']==1].copy()

def load_driving_data(filename):
    """
    Load and preprocess SOWFA driving data
    """
    # Load SOWFA forcing file
    df = load_sowfa_forcing_file(filename)
    # Calculate quantities of interest
    calc_QOIs(df)
    return df 


def load_sowfa_forcing_file(fpath):
    """
    Load a SOWFA forcing file (specific fileformat for SOWFA)
    """
    from mmctools.mmctools.windtools.windtools.openfoam import InputFile

    f = InputFile(fpath)

    heights = np.array(f['sourceHeightsMomentum'])
    assert np.all(heights == np.array(f['sourceHeightsTemperature']))
    srcMomX = np.array(f['sourceTableMomentumX'])
    srcMomY = np.array(f['sourceTableMomentumY'])
    srcMomZ = np.array(f['sourceTableMomentumZ'])
    srcTemp = np.array(f['sourceTableTemperature'])
    
    # Cast data to pandas dataframe
    dflist = []
    for i in range(srcMomX[:,0].size):
        data = {}
        data['height'] = f['sourceHeightsMomentum']
        data['u'] = srcMomX[i,1:]
        data['v'] = srcMomY[i,1:]
        data['w'] = srcMomZ[i,1:]
        data['theta'] = srcTemp[i,1:]
        df = pd.DataFrame(data=data)
        df['t'] = srcMomX[i,0]
        dflist.append(df)
    df = pd.concat(dflist)
    
    # Convert time in seconds to datetime
    df['t'] = pd.to_timedelta(df['t'],unit='s') + pd.to_datetime('2013-11-08 00:00:00')
    
    df['t'] = df['t'].dt.round('10min')
    
    # Convert to standard names
    df.rename({'t':'datetime'},axis='columns',inplace=True)
    
    # Set multi-index with levels datetime and height
    df.set_index(['datetime','height'],inplace=True)
    return df


# -------------------------------------
# Calculating rotor-averaged quantities
# -------------------------------------

def calc_rotor_average(df,zhub,diameter):
    """
    Calculate rotor-averaged quantities
    and return a new dataframe
    """
    from scipy.interpolate import interp1d

    zlow  = zhub - diameter/2.
    zhigh = zhub + diameter/2.

    data = df.unstack()
    heights = df.index.get_level_values(1).unique()
    dfout = {}
    for field in df.columns:
        dfout[field] = np.mean(interp1d(heights,data[field].values,axis=1,fill_value='extrapolate')(np.linspace(zlow,zhigh,21)),axis=1)

    dfout['alpha'] = calc_alpha(heights,data['wspd'].values,zhub,diameter)
    dfout['psi']   = calc_psi(heights,data['wdir'].values,zhub,diameter)
    return pd.DataFrame(dfout,index=data.index)


def calc_alpha(z,S,zh,D,ax=-1):
    '''
    Compute wind shear exponent by fitting a power law to
    the velocity profile over the rotor disk region
    Parameters
    ----------
    z: numpy 1D array of Nz
        heights
    S: numpy nD array
        wind speed [m/s]
    zh,D: float
        hub height and radius of rotor disk
    ax: int
        axis corresponding to the vertical direction
        default: last axis
    Returns
    -------
    alpha: numpy nD-1 array
        wind shear exponent
    '''
    from scipy.interpolate import interp1d
    from scipy.optimize import curve_fit

    z1 = zh - D/2.
    z2 = zh + D/2.
    zcc = np.linspace(z1,z2,10)

    #Move specified ax to last position, then reshape to 2d array and iterate
    Nz = S.shape[ax]
    N  = int(S.size/Nz)
    new_shape = np.moveaxis(S,ax,-1).shape[:-1]
    alpha = np.zeros((N))
    for i in range(N):
        Sint = interp1d(z,np.moveaxis(S,ax,-1).reshape(N,Nz)[i,:],fill_value='extrapolate')(zcc)
        Shub = interp1d(z,np.moveaxis(S,ax,-1).reshape(N,Nz)[i,:],fill_value='extrapolate')(zh)

        f = lambda x, alpha: Shub*(x/zh)**alpha
        popt,_ = curve_fit(f,zcc,Sint,1.0)
        alpha[i] = popt[0]
    return alpha.reshape(new_shape)


def calc_psi(z,WD,zh,D,ax=-1):
    '''
    Compute wind veer by fitting a line to
    the wind direction profile over the rotor disk region
    Parameters
    ----------
    z: numpy 1D array of Nz
        heights
    WD: numpy nD array
        wind direction [degrees]
    zh,D: float
        hub height and radius of rotor disk
    ax: int
        average wind veer over the rotor disk
        default: last axis
    Returns
    -------
    psi: numpy nD-1 array
        average wind veer over the rotor disk
    '''
    from scipy.interpolate import interp1d
    from scipy.optimize import curve_fit

    z1 = zh - D/2.
    z2 = zh + D/2.
    zcc = np.linspace(z1,z2,10)

    #Move specified ax to last position, then reshape to 2d array and iterate
    Nz = WD.shape[ax]
    N  = int(WD.size/Nz)
    new_shape = np.moveaxis(WD,ax,-1).shape[:-1]
    psi = np.zeros((N))
    for i in range(N):
        WDint = interp1d(z,np.moveaxis(WD,ax,-1).reshape(N,Nz)[i,:],fill_value='extrapolate')(zcc)
        WDhub = interp1d(z,np.moveaxis(WD,ax,-1).reshape(N,Nz)[i,:],fill_value='extrapolate')(zh)

        f = lambda x, *p: p[0]*x + p[1]
        popt,_ = curve_fit(f,zcc-zh,WDint-WDhub,[1.0,0.0])
        psi[i] = popt[0]
    return psi.reshape(new_shape)


# ------------------------------
# Calculating turbulence spectra
# ------------------------------

def reindex_if_needed(df,dt=None):
    """
    Check whether timestamps are equidistant with step dt (in seconds). If dt is not
    specified,  dt is equal to the minimal timestep in the dataframe. If timestamps
    are not equidistant, interpolate to equidistant time grid with step dt.
    """
    dts = np.diff(df.index.get_level_values(0).unique())/pd.to_timedelta(1,'s')

    # If dt not specified, take dt as the minimal timestep
    if dt is None:
        dt = np.min(dts)

    if not np.allclose(dts,dt):
        # df is missing some timestamps, which will cause a problem when computing spectra.
        # therefore, we first reindex the dataframe
        start = df.index.levels[0][0]
        end   = df.index.levels[0][-1]
        new_index = pd.date_range(start,end,freq=pd.to_timedelta(dt,'s'),name='datetime')
        return df.unstack().reindex(new_index).interpolate(method='index').stack()
    else:
        return df
