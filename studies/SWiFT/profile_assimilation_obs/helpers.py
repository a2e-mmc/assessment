"""
Helper functions specifically for the
profile asimilation coupled with Obs study
"""
import os, sys
import numpy as np
import pandas as pd
import xarray
from scipy.interpolate import interp1d

# manually add a2e-mmc repos to PYTHONPATH if needed
module_path = os.path.join(os.environ['HOME'],'tools','a2e-mmc')
if module_path not in sys.path:
    sys.path.append(module_path)

from mmctools.helper_functions import calc_wind, covariance, power_spectral_density, theta

# manually add NWTC/datatools repo to PYTHONPATH
module_path = os.path.join(os.environ['HOME'],'tools')
if module_path not in sys.path:
    sys.path.append(module_path)
            
from datatools.SOWFA6.postProcessing.averaging import PlanarAverages
from datatools.SOWFA6.postProcessing.probes import Probe
from datatools.SOWFA6.postProcessing.sourceHistory import SourceHistory
from datatools import openfoam_util



# ----------------------
# Loading reference data
# ----------------------

def load_wrf_reference_data(fpath):
    """
    Load WRF reference data
    """
    # Load data with xarray
    xa = xarray.open_dataset(fpath)
    # Convert to pandas dataframe
    wrf = xa.to_dataframe()
    # Convert to standard names
    wrf.rename({'U':'u','V':'v','W':'w','UST':'u*'},
               axis='columns',inplace=True)
    # Compute wind speed and wind direction
    wrf['wspd'], wrf['wdir'] = calc_wind(wrf)

    return wrf


def load_radar_reference_data(fpath):
    """
    Load TTU radar reference data
    """
    radar = pd.read_csv(fpath,parse_dates=True,index_col=['datetime','height'])
    # Extract scan types 0 and 1
    radar_scan0 = radar.loc[radar['scan_type']==0].copy()
    radar_scan1 = radar.loc[radar['scan_type']==1].copy()

    return radar_scan0, radar_scan1


def load_tower_reference_data(fpath):
    """
    Load TTU tower 10-min reference data
    """
    return pd.read_csv(fpath,parse_dates=True,index_col=['datetime','height'])


def load_tower_reference_spectra(fpath,times,heights,interval,window_size):
    """
    Load TTU tower 1-Hz data and compute spectra
    """
    # Load data
    tower = pd.read_csv(fpath,parse_dates=True,index_col=['datetime','height'])
    # Calculate some QoI
    tower['wspd'], tower['wdir'] = calc_wind(tower)
    tower['theta'] = theta(tower['Ts'],tower['p'])
    # Interpolate data to specified heights
    tower_hgt = interpolate_to_heights(tower,heights)
    # Reindex if needed
    tower_hgt = reindex_if_needed(tower_hgt)
    # Compute spectra
    tower_spectra = calc_spectra(tower_hgt,times,heights,interval,window_size)
    return tower_spectra


# -------------------------------------------------
# Calculating statistics and quantities of interest
# -------------------------------------------------

def calc_QOIs(df):
    """
    Calculate derived quantities (IN PLACE)
    """
    df['wspd'],df['wdir'] = calc_wind(df)
    df['u*'] = (df['uw']**2 + df['vw']**2)**0.25
    df['TKE'] = 0.5*(df['uu'] + df['vv'] + df['ww'])
    ang = np.arctan2(df['v'],df['u'])
    df['TI'] = df['uu']*np.cos(ang)**2 + 2*df['uv']*np.sin(ang)*np.cos(ang) + df['vv']*np.sin(ang)**2
    df['TI'] = np.sqrt(df['TI']) / df['wspd']


# ------------------------------
# Calculating turbulence spectra
# ------------------------------

def interpolate_to_heights(df,heights):
    """
    Interpolate data in dataframe to specified heights
    and return a new dataframe
    """
    # Unstack to single height index (= most time-consuming operation)
    unstacked = df.unstack(level='datetime')
    # Interpolate to specified heights
    f = interp1d(unstacked.index,unstacked,axis=0,fill_value='extrapolate')
    for hgt in heights:
        unstacked.loc[hgt] = f(hgt)
    # Restack and set index
    df_out = unstacked.loc[heights].stack().reset_index().set_index(['datetime','height']).sort_index()
    return df_out


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


def calc_spectra(df,times,heights,interval,window_size):
    """
    Calculate spectra for a given number of times and heights
    and return a new dataframe
    """
    dflist = []
    for tstart in times:
        for height in heights:
            spectra = power_spectral_density(df.xs(height,level='height'),
                                             tstart=pd.to_datetime(tstart),
                                             interval=interval,
                                             window_size=window_size)
            spectra['datetime'] = pd.to_datetime(tstart)
            spectra['height'] = height
            spectra.reset_index(inplace=True)
            dflist.append(spectra)
    df_spectra = pd.concat(dflist)
    return df_spectra.set_index(['datetime','height','frequency']).sort_index()
