#!/usr/bin/env python

import argparse
import pandas as pd
import xarray as xr
from mmctools.helper_functions import covariance

def calc_raw2stat(rawdata,offset='10min'):
    # Read raw data and convert to pandas DataFrame
    df = xr.open_dataset(rawdata).to_dataframe()

    # Drop p and p_rgh as these are not used for the analysis
    df.drop(labels=['p','p_rgh'],axis=1,inplace=True)

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
    return stats.to_xarray()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('rawdata', type=str)
    parser.add_argument('outputfile', type=str)
    args = parser.parse_args()

    xa = calc_raw2stat(args.rawdata).to_netcdf(args.outputfile)
