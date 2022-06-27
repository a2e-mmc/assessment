import numpy as np
import xarray as xr
import pandas as pd

def calc_streamwise (df, dfwdir, height=None, extrapolateHeights=True, interpolateTime = True):
    '''
    Calculate streamwise and cross stream velocity components
    
    Parameters
    ==========
    df : Dataset
        Dataset containing the u, v, w with at least datetime coordinate
    dfdir : Dataset, DataArray, or dataframe
        Dataset containing planar-average wind direction
    height: scalar
        If the input dataset only exists in one height (e.g. from VTKs),
        specify the height in this variable
    extrapolateHeights: bool
        Whether or not to extrapolate the wdir from heights given in dfwdir
        onto the heights asked in df. Useful if vertical slices that contain 
        all heights. True by default.
      
    '''
    
    # Datasets, DataArrays, or dataframes
    if not isinstance(dfwdir,xr.Dataset):
        if isinstance(dfwdir,pd.DataFrame):
            dfwdir = dfwdir.to_xarray( name = 'wdir' )
        elif isinstance(dfwdir,xr.DataArray):
            dfwdir = dfwdir.to_dataset( name = 'wdir' )
        else:
            raise ValueError(f'unsupported type: {type(dfwdir)}')
            
    
    # Drop variables from dfwdir
    varsToDrop = set(dfwdir.variables.keys()) - set(['datetime','height','z','wdir','time'])
    dfwdir = dfwdir.drop(list(varsToDrop))

    print(dfwdir.dims)
    
    # Interpolate wdir from planar average into the coordinates of df
    if interpolateTime:
        if 'height' in list(df.variables.keys()):
            heightInterp = height=df.height
        elif 'z' in list(df.variables.keys()):
            heightInterp = df.z.values[0]
        elif height != None:
            heightInterp = height
        else:
            raise NameError("The input dataset does not appear to have a 'height' or 'z' coordinate. Use `height=<scalar>` to specify one.")

        print(heightInterp, df.datetime)
        if extrapolateHeights:
            wdir_at_same_coords = dfwdir.interp(datetime=df.datetime, SLICES_Z=heightInterp, kwargs={"fill_value": "extrapolate"})
        else:
            wdir_at_same_coords = dfwdir.interp(datetime=df.datetime)
    else:
        print(dfwdir.dims)
        #wdir_at_same_coords = xr.DataArray( dfwdir.values, dims = 'datetime' )
        wdir_at_same_coords = dfwdir
    # Add wdir information to main dataset
    rotdf = xr.combine_by_coords([df, wdir_at_same_coords], compat = 'override')
    wdir = rotdf['wdir']
    
    # Rotate flowfield
    ustream = rotdf['u']*np.cos(np.deg2rad(270-wdir)) + rotdf['v']*np.sin(np.deg2rad(270-wdir))
    vcross =  rotdf['u']*np.sin(np.deg2rad(270-wdir)) - rotdf['v']*np.cos(np.deg2rad(270-wdir))

    return ustream, vcross, wdir

########################
# calc_mean_pert
########################
# Calculate bar and prime (mean and perturbation) quantities
#

def calc_mean_pert( ds, variable_list = ['U', 'V', 'W'], mean_wind_dir = 'periodic',\
                   mean_dim_periodic = ('nx', 'ny'), mean_dim_nonperiodic = 'nx'  ):
    '''
    Purpose of this function is to compute the mean and perturbation quantities for computing fluxes and stresses.
    
        ds: xarray Dataset. Contains the coords, dims, and variables (U,V,W) 
            that have been computed by the postprocessing function above
        variable_list: array-like. Contains variable names (strings) for mean/perturbation quantities.
            Must be 4-D variables using x/y/z coords, error-catches are not implemented.
        mean_wind_dir: either 'periodic' (default) or 'zonal' (i.e. mean wind dir is from west to east).
            periodic: compute means on x/y planes to get mean quantities as a function of time and height
            zonal: mean quantities will be computed on lines of constant x, so mean will also be a function of x.
                this means less statistical power, and some temporal averaging may be required, but that is not
                accounted for in this function (yet)
    '''
    
    mean_str_suff = '_bar'
    pert_str_suff = '_p'
    
    for vv in variable_list:
        print(vv)
        mean_str = vv + mean_str_suff
        pert_str = vv + pert_str_suff
        
        if mean_wind_dir == 'periodic':
            print("Periodic simulation")
            ds[mean_str] = ds[vv].mean(dim = mean_dim_periodic)
            ds[pert_str] = ds[vv] - ds[mean_str]
        elif mean_wind_dir == 'zonal':
            print("Zonal simulation, may need some temporal averaging for power")
            ds[mean_str] = ds[vv].mean(dim = mean_dim_nonperiodic )
            ds[pert_str] = ds[vv] - ds[mean_str]
    return ds








########################
# calc_stresses
########################
# Calculate resolved stress terms
#

def calc_stresses( ds, mean_dims = ('nx', 'ny'), do_uw = True, do_vw = False, do_uv = False, \
                 U_vname = 'U_p',
                 V_vname = 'V_p',
                 W_vname = 'W_p'):
    '''
    Calculate components of the Stress-Energy tensor relevant to shear production of turbulence
        ds: xarray dataset.
        do_uw: Boolean (default True). If true, calculates tau13 (the u'w' component of the stress energy tensor)
        do_vw: Boolean (default False). If true, calculates tau23 (the u'w' component of the stress energy tensor)
        do_uv: Boolean (default False). If true, calculates tau12 (the u'w' component of the stress energy tensor)
    '''
    if do_uw:
        print('calculating tau13...')
        ds['tau13'] = ( ds[U_vname] * ds[W_vname] ).mean(dim = mean_dims )
    if do_vw:
        print('calculating tau23...')
        ds['tau23'] = ( ds[V_vname] * ds[W_vname] ).mean(dim = mean_dims )
    if do_uv:
        print('calculating tau12...')
        ds['tau12'] = ( ds[U_vname] * ds[V_vname] ).mean(dim = mean_dims )
        
    return ds
    

########################
# calc_tke
########################
# Calculate resolved TKE
#

def calc_tke( ds, mean_dims = ('nx', 'ny'), \
                 U_vname = 'U_p',
                 V_vname = 'V_p',
                 W_vname = 'W_p'):
    import xarray as xr
    '''
    Calculates RESOLVE LES TKE. Does not compute the subgrid component.
        ds: xarray dataset.
        mean_dims: dimensions to take the mean over, default 'nx' and 'ny'
    '''
    print("calculating TKE...")
    ds['TKE'] = 1./2. * (  ( xr.ufuncs.square( ds[U_vname] ) ).mean(dim = mean_dims ) \
                         + ( xr.ufuncs.square( ds[V_vname] ) ).mean(dim = mean_dims ) \
                         + ( xr.ufuncs.square( ds[W_vname] ) ).mean(dim = mean_dims ) )
    
    return ds