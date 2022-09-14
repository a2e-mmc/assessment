"""
Helper functions specifically for the
profile asimilation coupled with WRF study
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# manually add a2e-mmc repos to PYTHONPATH if needed
module_path = os.path.join(os.environ['HOME'],'a2e-mmc')
if module_path not in sys.path:
    sys.path.append(module_path)

from mmctools.helper_functions import calc_wind, covariance, power_spectral_density, theta

# manually add windtools repo to PYTHONPATH
module_path = os.path.join(os.environ['HOME'])
if module_path not in sys.path:
    sys.path.append(module_path)
from windtools.SOWFA6.postProcessing.averaging import PlanarAverages
from windtools.SOWFA6.postProcessing.probes import Probe
from windtools.SOWFA6.postProcessing.sourceHistory import SourceHistory
from windtools.openfoam import InputFile



# Coriolis parameter at the SWiFT site
latitude = 33.61054 # [degrees]
fc = 2 * ( 2*np.pi/(24.0 * 3600.0) ) * np.sin(np.deg2rad(latitude))


# Reference date for SOWFA time (seconds since)
tref = '2013-11-08 00:00:00'



# ------------------
# Loading sowfa data
# ------------------

def load_sowfa_data(dpath,times,heights,interval='1h',window_size='10min'):
    """
    Load and process data of a particular SOWFA simulation

    Data loaded includes:
    - Planar averages
    - Virutal tower data
    - Source history data
    """
    # Load planar average data
    # ------------------------
    fpath = os.path.join(dpath,'postProcessing/planarAverages')
    df_pavg = reader_planar_average(fpath)
    
    # Calculate 10-min averages and quantities of interest
    df_pavg_10min = df_pavg.unstack().resample('10min').mean().stack()
    calc_QOIs(df_pavg_10min)
    
    # Calculate hourly averages
    df_pavg_60min = df_pavg_10min.unstack().resample('60min').mean().stack()
    
    
    # Load virtual tower data
    # -----------------------
    fpath = os.path.join(dpath,'postProcessing/probe1')
    df_prob = reader_probe(fpath)
    df_prob['wspd'], df_prob['wdir'] = calc_wind(df_prob)
    
    # Calculate 10-min statistics and quantities of interest
    df_prob_10min = calc_stats(df_prob)
    calc_QOIs(df_prob_10min)
    
    # Calculate hourly averages
    df_prob_60min = df_prob_10min.unstack().resample('60min').mean().stack()
    
    # Calculate frequency spectrum
    df_prob_hgts = interpolate_to_heights(df_prob,heights)
    df_prob_hgts = reindex_if_needed(df_prob_hgts)
    df_prob_spectra = calc_spectra(df_prob_hgts,times,heights,interval,window_size)
    
    
    # Load source history data
    # ------------------------
    fpath = os.path.join(dpath,'postProcessing/sourceHistory')
    df_srch = reader_source_history(fpath)
    
    # Rescale budget components
    rescale_budget_components(df_srch)
    
    # Calculate 10-min and hourly averages
    df_srch_10min = df_srch.unstack().resample('10min').mean().stack()
    df_srch_60min = df_srch.unstack().resample('60min').mean().stack()
    
    return df_prob_10min, df_prob_60min, df_prob_spectra, df_pavg_10min, df_pavg_60min, df_srch, df_srch_10min, df_srch_60min


def reader_probe(fpath):
    """
    Read sowfa probe file and convert to standard pandas dataframe
    """
    # Read in virtual tower data and convert to pandas DataFrame
    df = Probe(fpath).to_pandas()
    
    # Convert time in seconds to datetime
    df.reset_index(inplace=True)
    df['t'] = pd.to_timedelta(df['t'],unit='s') + pd.to_datetime(tref)
    
    # Rename columns
    df.columns = ['datetime', 'height', 'u', 'v', 'w', 'theta']
    
    # Set multi-index with levels datetime and height
    df.set_index(['datetime','height'],inplace=True)
    return df


def reader_planar_average(fpath):
    """
    Read sowfa planar average file and convert to standard pandas dataframe
    """
    # Read in planar average data and convert to pandas DataFrame
    df = PlanarAverages(fpath,varList=['U','UU','T']).to_pandas()
    
    # Convert time in seconds to datetime
    df.reset_index(inplace=True)
    df['t'] = pd.to_timedelta(df['t'],unit='s') + pd.to_datetime(tref)
    
    # Rename columns
    df.columns = ['datetime', 'height', 'u', 'v', 'w', 'uu', 'uv', 'uw', 'vv', 'vw', 'ww', 'theta']
    
    # Set multi-index with levels datetime and height
    df.set_index(['datetime','height'],inplace=True)
    return df


def reader_source_history(fpath):
    """
    Read sowfa source history file and convert to standard pandas dataframe
    """
    # Read in source history data and convert to pandas DataFrame
    df = SourceHistory(fpath).to_pandas()
    
    # Convert time in seconds to datetime
    df.reset_index(inplace=True)
    df['t'] = pd.to_timedelta(df['t'],unit='s') + pd.to_datetime(tref)
    
    # Rename columns
    df.columns = ['datetime', 'height', 'Fu', 'Fv', 'Fw', 'Ft']
    
    # Set multi-index with levels datetime and height
    df.set_index(['datetime','height'],inplace=True)
    return df


# ----------------------
# Loading reference data
# ----------------------

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

    Return both 10-min data and hourly averaged data
    """
    df_10min = pd.read_csv(fpath,parse_dates=True,index_col=['datetime','height'])
    df_60min = df_10min.unstack().resample('60min').mean().stack()
    return df_10min, df_60min


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


def load_wrf_reference_data(dpath):
    """
    Load WRF reference data from SOWFA input files
    """
    # Load wrf field data
    # ------------------
    wrf_pavg_10min = load_sowfa_input_file(os.path.join(dpath,'fieldTable'))

    # Drop value at z=0
    wrf_pavg_10min = wrf_pavg_10min.loc[(slice(None),wrf_pavg_10min.index.get_level_values(1)>0),:]

    # Calculate wind speed and direction
    wrf_pavg_10min['wspd'], wrf_pavg_10min['wdir'] = calc_wind(wrf_pavg_10min)

    # Calculate hourly averages
    wrf_pavg_60min = wrf_pavg_10min.unstack().resample('60min').mean().stack()

    # Load wrf forcing data
    # ---------------------
    wrf_srch_10min = load_sowfa_input_file(os.path.join(dpath,'forcingTable'))

    # Drop value at z=0
    wrf_srch_10min = wrf_srch_10min.loc[(slice(None),wrf_srch_10min.index.get_level_values(1)>0),:]

    # Rename columns and rescale budget components
    wrf_srch_10min.rename({'u':'Fu','v':'Fv','w':'Fw','theta':'Ft'},axis='columns',inplace=True)
    rescale_budget_components(wrf_srch_10min)

    # Calculate hourly averages
    wrf_srch_60min = wrf_srch_10min.unstack().resample('60min').mean().stack()

    return wrf_pavg_10min, wrf_pavg_60min, wrf_srch_10min, wrf_srch_60min


def load_sowfa_input_file(fpath):
    """
    Load a specific SOWFA input file
    """
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
    df['t'] = pd.to_timedelta(df['t'],unit='s') + pd.to_datetime(tref)
    
    df['t'] = df['t'].dt.round('10min')
    
    # Convert to standard names
    df.rename({'t':'datetime'},axis='columns',inplace=True)
    
    # Set multi-index with levels datetime and height
    df.set_index(['datetime','height'],inplace=True)
    return df


# -------------------------------------------------
# Calculating statistics and quantities of interest
# -------------------------------------------------

def calc_stats(df,offset='10min'):
    """
    Calculate statistics for a given data frame
    and return a new dataframe
    """
    # calculate statistical quantities on unstacked 
    unstacked = df.unstack()
    stats = unstacked.resample(offset).mean().stack()
    # - calculate variances
    stats['uu'] = unstacked['u'].resample(offset).var().stack()
    stats['vv'] = unstacked['v'].resample(offset).var().stack()
    stats['ww'] = unstacked['w'].resample(offset).var().stack()
    # - calculate covariances
    stats['uv'] = covariance(unstacked['u'], unstacked['v'], offset, resample=True).stack()
    stats['vw'] = covariance(unstacked['v'], unstacked['w'], offset, resample=True).stack()
    stats['uw'] = covariance(unstacked['u'], unstacked['w'], offset, resample=True).stack()
    stats['thetaw'] = covariance(unstacked['theta'], unstacked['w'], offset, resample=True).stack()
    return stats


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


def rescale_budget_components(df):
    """
    Rescale budget components 'Fu', 'Fv' and 'Ft'
    from [m/s2] to [m/s] and from [K/s] to [K/h] 
    """
    df['Fu'] = df['Fu']/fc
    df['Fv'] = df['Fv']/fc
    df['Ft'] = df['Ft']*3600.0


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
