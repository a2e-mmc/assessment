from os import path
import numpy as np
import xarray as xr
from mmctools.wrf.utils import tsout_seriesReader

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - - - - - - - -   USER SETTINGS   - - - - - - - - - - - - - - - - #
# Base location of simulations (2 levels below individual restarts)
wrf_dir   = '/glade/scratch/hawbecke/WRF/MMC/FINO1/20100512_to_20100526/'
# Name of the case (e.g. GFS_YSU_KF_CHRN_OST) (1 level below individual restarts)
cases     = ['ERAI_YSU_KF_CHRN_OST']

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# - - - - - - - - - - - - - - END USER SETTINGS - - - - - - - - - - - - - - - - #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


# These should not change: setup script was written to have this be standard
restarts  = ['FINO_2010051212', 'FINO_2010051412', 'FINO_2010051612', 'FINO_2010051812',
             'FINO_2010052012', 'FINO_2010052212', 'FINO_2010052412']
wrf_start = ['2010-05-12 12:00:00', '2010-05-14 12:00:00', '2010-05-16 12:00:00', '2010-05-18 12:00:00',
             '2010-05-20 12:00:00', '2010-05-22 12:00:00', '2010-05-24 12:00:00']
ncases = len(cases)


wrf_twrs = {}
for case in cases:
    case_dir = '{}{}/'.format(wrf_dir,case)
    wrf_twrs[case] = {'FN1': [],
                      'FN2': [],
                      'FN3': []}

    for ff, fino in enumerate([1,2,3]):
        twr_path = '{}/FINO{}.nc'.format(case_dir,fino)
        if path.exists(twr_path):
            print('loading in full dataset!')
            wrf_twrs[case]['FN{}'.format(fino)] = xr.open_dataset(twr_path)
        else:
            if read_in_obs:
                obs_ds = fino_obs['FINO{}'.format(fino)]
                obs_levels = np.unique(sorted(np.concatenate([obs_ds.spd_levels.data, obs_ds.dir_levels.data,obs_ds.tmp_levels.data])))
            else:
                if fino == 1:
                    obs_levels = np.asarray([30.,  33.,  40.,  50.,  60.,  70.,  80.,  90., 100.])
                elif fino == 2:
                    obs_levels = np.asarray([30.,  31.,  32.,  40.,  42.,  50.,  51.,  52.,  62.,  70.,  71.,  72.,  82.,  91.,  92.,  99., 102.])
                elif fino == 3:
                    obs_levels = np.asarray([28.,  29.,  30.,  40.,  50.,  55.,  60.,  70.,  80.,  90.,  95.,  100., 106.])
                else:
                    print('fino can only be 1-3')
            wrf_twrs[case]['FN{}'.format(fino)] = tsout_seriesReader(case_dir,restarts,wrf_start,'d03',structure='unordered',
                                                                     time_step=10.0,select_tower=['FN{}'.format(fino)],
                                                                     heights=obs_levels,height_var='ph')
            wrf_twrs[case]['FN{}'.format(fino)].to_netcdf(twr_path)
            print('New file created: {}'.format(twr_path))
