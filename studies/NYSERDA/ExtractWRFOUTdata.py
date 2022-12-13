#!/glade/u/home/hawbecke/local/envs/mmc/bin/python

import pandas as pd
import numpy as np
import xarray as xr
from os import path
import sys
import glob
from mmctools.helper_functions import calc_wind
from netCDF4 import Dataset
from wrf import getvar

main_directory = '/glade/scratch/hawbecke/WRF/MMC/NYSERDA/SENSITIVITY_SUITE/production/'
#main_directory = '/glade/work/hawbecke/MMC/NYSERDA/WRF/SENSITIVITY_SUITE/'

interp_data = True

icbc_type = 'MERRA2'  # ERAI, ERA5, FNL

from NYSERDA_case_dict import case_dict

cases = [case_dict[x]['case_str'] for x in list(case_dict.keys())]

# Shallow only:
cases = cases[-2:-1]

# TESTING:
#cases = [cases[0]] + cases[3:6]# + [cases[-1]]

# SST Only:
#cases = cases[:-2]

eta_level_testing = False
if eta_level_testing:
    main_directory += 'eta_level_test/'
    cases = cases[-3:]
    time_step = 5.0
else:
    cases = cases[:]
    failed_cases = []
    for case in failed_cases:
        cases.remove(case)
    time_step = 6.0

sim_start = '2020-04-04 06:00:00'
#case_start = '2020-04-05 16:00:00'
#case_end   = '2020-04-06 06:00:00'
case_start = '2020-04-06 01:00:00'
case_end   = '2020-04-06 06:00:00'

best_case  = 'WPS3_WRF1'
n_profiles = 8

dom_dict = {1: {'dt':15,
                'marker': r'$ 1 $',
                'ls':'-'},
            2: {'dt':5,
                'marker': r'$ 2 $',
                'ls':'--'},
            3: {'dt':1,
                'marker': r'$ 3 $',
                'ls':'-'},
            4: {'dt':0.2,
                'marker': r'$ 4 $',
                'ls':'-.'},
            5: {'dt':0.04,
                'marker': r'$ 5 $',
                'ls':':'},
            }
doms_of_interest = [1,2,3,4,5]
#doms_of_interest = [1,2,3,4]

#for cc,case in enumerate(cases):
#    plt.plot(np.arange(0,10),np.arange(0,10)*cc,c=case_dict[case]['color'])


# Get obs for level information
obs_dir = '/glade/work/hawbecke/MMC/NYSERDA/obs/'
obs_res = 'high' # 'low' or 'high'
atlshrs_full = xr.open_dataset('{}ATLSHORES_lidars.nc'.format(obs_dir))
atlshrs_obs = atlshrs_full.sel(case='obs').squeeze(dim='Tlevels')
atlshrs_obs.T.data += 273.15
atlshrs_obs.sst.data += 273.15

if obs_res == 'low':
    nyserda_full = xr.open_dataset('{}NYSERDA_lidars.nc'.format(obs_dir))
    nyserda_obs = nyserda_full.sel(case='obs')
    
    #best_case = 'WPS3_WRF1'
    #nyserda_opt = nyserda_full.sel(case=best_case)
    hr_nyserda_obs = xr.open_dataset('{}HR_NYSERDA_lidars.nc'.format(obs_dir))
else:
    lr_nyserda_obs = xr.open_dataset('{}NYSERDA_lidars.nc'.format(obs_dir))
    lr_nyserda_obs = lr_nyserda_obs.sel(case='obs',Tlevels=2).drop('Tlevels')
    lr_nyserda_obs = lr_nyserda_obs.rename({'Ulevels':'level'})

    nyserda_obs = xr.open_dataset('{}HR_NYSERDA_lidars.nc'.format(obs_dir))
    nyserda_obs = nyserda_obs.drop('T')
    nyserda_obs = nyserda_obs.squeeze().drop('Tlevels')
    nyserda_obs = nyserda_obs.rename({'Ulevels':'level'})

    sst = lr_nyserda_obs.sst.sel(datetime=slice(nyserda_obs.datetime[0],nyserda_obs.datetime[-1]))
    nyserda_obs = xr.merge([nyserda_obs,sst])
    t2 = lr_nyserda_obs.T.sel(datetime=slice(nyserda_obs.datetime[0],nyserda_obs.datetime[-1])).squeeze()
    nyserda_obs = xr.merge([nyserda_obs,t2])




def get_data(ds,vars_to_extract=None):
    ua = (ds.U[:,:,1:] + ds.U[:,:,:-1])*0.5
    ua = ua.rename({'west_east_stag':'west_east'})
    va = (ds.V[:,1:,:] + ds.V[:,:-1,:])*0.5
    va = va.rename({'south_north_stag':'south_north'})
    wa = (ds.W[1:,:,:] + ds.W[:-1,:,:])*0.5
    wa = wa.rename({'bottom_top_stag':'bottom_top'})
    ds['ua'] = ua
    ds['va'] = va
    ds['wa'] = wa
    wspd,wdir = calc_wind(ds,u='ua',v='va')
    ds['wspd'] = wspd
    ds['wdir'] = wdir
    z = (ds.PH + ds.PHB)/9.81 - ds.HGT
    zs = (z[1:,:,:] + z[:-1,:,:])*0.5
    zs = zs.rename({'bottom_top_stag':'bottom_top'})
    #ds['z'] = z
    ds['zs'] = zs

    if 'TKE_PBL' in vars_to_extract:
        temp = (ds['TKE_PBL'][1:,:,:] + ds['TKE_PBL'][:-1,:,:])*0.5
        temp = temp.rename({'bottom_top_stag':'bottom_top'})
        del(ds['TKE_PBL'])
        ds['TKE_PBL'] = temp
    
    les_vars = ['m11','m22','m33']
    for varn in les_vars:
        if varn not in list(ds.data_vars):
            print('creating {}'.format(varn))
            ds[varn] = ds.wspd*0.0
            ds[varn].name = varn
    if vars_to_extract is not None:
        ds = ds[vars_to_extract]
    return(ds)

def get_locs_for_domain(loc_dict,
                        tower_dir,
                        case,
                        dom,
                        twr,
                        window=None):
    ts_fname = '{}NYSERDA_{}_towers_d0{}.nc'.format(tower_dir,case,dom)
    if dom >= 4:
        ts_fname = ts_fname.replace('.nc','_chunk0.nc')
    ts = xr.open_dataset(ts_fname)
    ts = ts.sel(station=twr)
    ts_i = ts.i.data
    ts_j = ts.j.data
    if window is not None:
        if type(window) is int:
            half_window = int((window - 1)*0.5)
            ts_i = np.arange(ts_i-half_window,ts_i+half_window+1)
            ts_j = np.arange(ts_j-half_window,ts_j+half_window+1)
        elif type(window) is list:
            ts_i = np.asarray(window)
            ts_j = np.asarray(window)
            
    loc_dict[dom] = {'i':ts_i,'j':ts_j}
    return(loc_dict)

wrfout_file_dict = {}

wrfout_s = pd.to_datetime('2020-04-06 01:00')
wrfout_e = pd.to_datetime('2020-04-06 06:00')

for dd,dom in enumerate([1,2,3,4,5]):
    file_list = []
    wrfouts = sorted(glob.glob('{}{}/wrfout_d0{}*'.format(main_directory,cases[0],dom)))
    for ww,wrfout in enumerate(wrfouts):
        wrf_time = pd.to_datetime(' '.join((wrfout.split('/')[-1]).split('_')[2:]))
        if (wrf_time >= wrfout_s) and (wrf_time <= wrfout_e):
            file_list += [wrfout.split('/')[-1]]
    wrfout_file_dict[dom] = file_list
    
vars_to_extract = ['wspd','wdir','T','TSK','ZNT','zs','TKE','QKE',
                   'T2','TKE_PBL','UST','HFX','LH','QFX','SST',
                   'XLAT','XLONG','ua','va','wa','m11','m22','m33']

loc_dict = {}
tower_dir = '/glade/scratch/hawbecke/WRF/MMC/NYSERDA/SENSITIVITY_SUITE/production/tower_netCDFs/'

twr_of_interest = 'E06'

large_domain = False

if large_domain:
    domain_window_dict = {1:None,
                          2:None,
                          3:24,
                          4:[],
                          5:[]}
else:
    domain_window_dict = {1:None,
                          2:None,
                          3:5,
                          4:25,
                          5:125}

twr_lat = 39.546772
twr_lon = -73.428892

if large_domain:
    for dd,dom in enumerate([4,5]):
        wrf_dom = xr.open_dataset('{}{}/{}'.format(main_directory,cases[0],wrfout_file_dict[dom][0]),decode_times=False).isel(Time=0)
        dom_lat = wrf_dom.XLAT
        dom_lon = wrf_dom.XLONG
        dist = np.sqrt((dom_lat - twr_lat)**2 + (dom_lon - twr_lon)**2)
        twr_j,twr_i = np.where(dist==np.min(dist))
        twr_j,twr_i = twr_j[0],twr_i[0]
        dom_u10 = wrf_dom.U10
        window_xe = len(wrf_dom.west_east)
        window_xs = int(window_xe/2)
        window_ye = len(wrf_dom.south_north)
        window_ys = int(window_ye/2)

        window_xe -= 10
        window_ye -= 10

        wrf_dom.close()
        domain_window_dict[dom] = list(np.arange(window_xs,window_xe+1))


for cc,case in enumerate(cases):
    
    if interp_data:
        if large_domain:
            case_ds_fname = '{}extracted_data/{}_interpolated_extracted_data_largeSubsection.nc'.format(main_directory,case)
        else:
            case_ds_fname = '{}extracted_data/{}_interpolated_extracted_data.nc'.format(main_directory,case)
    else:
        case_ds_fname = '{}extracted_data/{}_extracted_data.nc'.format(main_directory,case)
    print(case)
    if not path.exists(case_ds_fname):
        for dd,dom in enumerate([1,2,3,4,5]):
            print('-----------')
            print(dom)
            wrfouts = sorted(glob.glob('{}{}/wrfout_d0{}*'.format(main_directory,case,dom)))
            wrfouts = wrfout_file_dict[dom]
            for ww,wrfout in enumerate(wrfouts[:]):
                wrfout = '{}{}/{}'.format(main_directory,case,wrfout)
                wrf_time = pd.to_datetime(' '.join((wrfout.split('/')[-1]).split('_')[2:]))
                print(wrf_time)
                ds = xr.open_dataset(wrfout,decode_times=False).isel(Time=0).squeeze()
                if dom not in list(loc_dict.keys()):
                    print('Getting location for d0{}'.format(dom))
                    loc_dict = get_locs_for_domain(loc_dict,
                                                   tower_dir,
                                                   case,
                                                   dom,
                                                   twr_of_interest,
                                                   window=domain_window_dict[dom]
                                                   )

                ds = get_data(ds,vars_to_extract=vars_to_extract)

                ts_i = loc_dict[dom]['i']
                ts_j = loc_dict[dom]['j']
                ds = ds.sel(west_east=ts_i,south_north=ts_j)
                if domain_window_dict[dom] is not None:
                    ds_pert = ds.copy()
                    mean_dims = ['west_east','south_north']
                    ds_lat = ds.XLAT.mean(dim=mean_dims)
                    ds_lon = ds.XLONG.mean(dim=mean_dims)
                    ds = ds.mean(dim=mean_dims)
                    ds_pert -= ds
                    up2_mean = (ds_pert.ua**2).mean(dim=mean_dims)
                    vp2_mean = (ds_pert.va**2).mean(dim=mean_dims)
                    wp2_mean = (ds_pert.wa**2).mean(dim=mean_dims)
                    tke = 0.5*(up2_mean + vp2_mean + wp2_mean)
                    tke.name='TKEres'
                    ds['TKEres'] = tke
                    ds['XLAT'] = ds_lat
                    ds['XLONG'] = ds_lon
                else:
                    ds['TKEres'] = ds.TKE*0.0

                ds = ds.rename({'TKE':'TKEsgs'})
                if dom >= 3:
                    ds['TKEsgs'] += ds.m11 + ds.m22 + ds.m33

                if interp_data:
                    ds = ds.assign_coords({'level':ds.zs}).rename({'bottom_top':'level'})
                    ds = ds.interp(coords={'level':nyserda_obs.level.data},method='linear')
                ds = ds.expand_dims({'datetime':[wrf_time]},axis=0)
                ds = ds.drop(['XTIME'])
                if ww == 0:
                    ds_f = ds.copy()
                else:
                    ds_f = xr.merge([ds_f,ds])

                del(ds)

            lat = ds_f.XLAT
            lon = ds_f.XLONG
            if domain_window_dict[dom] is None:
                lat = lat.drop(['XLAT','XLONG'])
                lon = lon.drop(['XLAT','XLONG'])
            ds_f = ds_f.drop(['XLAT','XLONG'])
            ds_f['lat'] = lat
            ds_f['lon'] = lon

            ds_f = ds_f.expand_dims({'dom':[dom]})
            if dd == 0:
                case_ds = ds_f.copy()
            else:
                case_ds = xr.merge((case_ds,ds_f))
            del(ds_f)

        case_ds.to_netcdf(case_ds_fname)
