#!/usr/bin/env python

import argparse
import pandas as pd
import xarray as xr
from mmctools.windtools.windtools.SOWFA6.postProcessing.probes import Probe

def reader(datadir):
    # Read in virtual tower data and convert to pandas DataFrame
    df = Probe(datadir,varList=['U','T','p','p_rgh'],verbose=False).to_pandas()
    
    # Convert time in seconds to datetime
    df.reset_index(inplace=True)
    df['t'] = pd.to_timedelta(df['t'],unit='s') + pd.to_datetime('2013-11-08 00:00')
    
    # Rename columns
    df.columns = ['datetime', 'height', 'u', 'v', 'w', 'theta','p','p_rgh']
    
    # Set multi-index with levels datetime and height
    df.set_index(['datetime','height'],inplace=True)

    # Convert to xarray
    xa = df.to_xarray()

    # Set metadata
    xa.height.attrs['units'] = "m"
    xa.height.attrs['long_name'] = "height above ground level"
    xa.u.attrs['units'] = "m/s"
    xa.u.attrs['long_name'] = "west to east velocity component"
    xa.v.attrs['units'] = "m/s"
    xa.v.attrs['long_name'] = "south-to-north velocity component"
    xa.w.attrs['units'] = "m/s"
    xa.w.attrs['long_name'] = "vertical velocity component"
    xa.theta.attrs['units'] = "K"
    xa.theta.attrs['long_name'] = "potential temperature"
    xa.p.attrs['units'] = "m**2/s**2"
    xa.p.attrs['long_name'] = "kinematic pressure"
    xa.p_rgh.attrs['units'] = "m**2/s**2"
    xa.p_rgh.attrs['long_name'] = "modified kinematic pressure"

    return xa

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('datadir', type=str)
    parser.add_argument('outputfile', type=str)
    args = parser.parse_args()

    xa = reader(args.datadir).to_netcdf(args.outputfile)
