#!/usr/bin/env python

import argparse
import pandas as pd
import xarray as xr
import numpy as np
from scipy.interpolate import interp1d
from mmctools.helper_functions import calc_wind, power_spectral_density

def calc_raw2spectra(rawdata,spectra_heights,spectra_times,interval='1h',window_size='30min'):
    # Read raw data and convert to pandas DataFrame
    df = xr.open_dataset(rawdata).to_dataframe()

    # Drop p and p_rgh as these are not used for the analysis
    df.drop(labels=['p','p_rgh'],axis=1,inplace=True)

    # Compute wspd and wdir so that corresponding spectra are calculated as well
    df['wspd'],df['wdir'] = calc_wind(df)
    
    # Interpolate to specific heights
    df_hgt = interpolate_to_heights(df,spectra_heights)
    # Reindex if some timestamps are missing
    df_hgt = reindex_if_needed(df_hgt)
    # Calculate spectra at specific times
    df_spectra = calc_spectra(df_hgt,spectra_times,spectra_heights,interval,window_size)

    return df_spectra.to_xarray()

def interpolate_to_heights(df,heights):
    # Unstack to single height index (= most time-consuming operation)
    unstacked = df.unstack(level='datetime')
    # Interpolate to specified heights
    f = interp1d(unstacked.index,unstacked,axis=0,fill_value='extrapolate')
    for hgt in heights:
        unstacked.loc[hgt] = f(hgt)
    # Restack and set index
    df_out = unstacked.loc[heights].stack().reset_index().set_index(['datetime','height']).sort_index()
    return df_out

def reindex_if_needed(df):
    dts = np.diff(df.index.get_level_values(0).unique())/pd.to_timedelta(1,'s')
    if not np.allclose(dts,dts[0]):
        # df is missing some timestamps, which will cause a problem when computing spectra.
        # therefore, we first reindex the dataframe
        start = df.index.levels[0][0]
        end   = df.index.levels[0][-1]
        new_index = pd.date_range(start,end,freq=pd.to_timedelta(dts[0],'s'),name='datetime')
        return df.unstack().reindex(new_index).interpolate(method='index').stack()
    else:
        return df

def calc_spectra(df,times,heights,interval,window_size):
    dflist = []
    for tstart in times:
        for height in heights:
            spectra = power_spectral_density(df.xs(height,level='height'),
                                             tstart=pd.to_datetime(tstart),
                                             interval=interval,
                                             window_size=window_size,
                                             window_type="hann"
                                             )
            spectra['datetime'] = pd.to_datetime(tstart)
            spectra['height'] = height
            dflist.append(spectra)
    df_spectra = pd.concat(dflist)
    df_spectra.reset_index(inplace=True)
    return df_spectra.set_index(['datetime','height','frequency']).sort_index()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('rawdata', type=str)
    parser.add_argument('outputfile', type=str)
    args = parser.parse_args()

    # For calculation of spectra for hourly intervals
    #xa = calc_raw2spectra(args.rawdata,
    #    spectra_heights=[10.,32.,80.,120.],
    #    spectra_times = pd.date_range(start='2013-11-08 12:00:00',end='2013-11-09 11:00:00',freq='1h')
    #).to_netcdf(args.outputfile)

    # For calculation of one spectrum for the entire diurnal cycle
    xa = calc_raw2spectra(args.rawdata,
        spectra_heights=[75.001],
        spectra_times = ['2013-11-08 12:00:00',],
        interval='24h',
        window_size='24h',
    ).to_netcdf(args.outputfile)

