#!/glade/u/home/hawbecke/local/envs/mmc/bin/python
import xarray as xr
from os import path
from mmctools.wrf.utils import tsout_seriesReader
import glob

# - - - - USER SETTINGS - - - - #
# don't forget the trailing slash!
main_directory = '/glade/scratch/hawbecke/WRF/MMC/NYSERDA/SENSITIVITY_SUITE/production/' 
# - - - END USER SETTINGS - - - #
from NYSERDA_case_dict import case_dict
cases = [case_dict[x]['case_str'] for x in list(case_dict.keys())]

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

wrf_data = {}
avg_twr_dict = {}

for dd,dom in enumerate(doms_of_interest):
    
    time_step = dom_dict[dom]['dt']
    dom_str = 'd0{}'.format(dom)

    #for cc,case in enumerate(cases):
    for cc,case in enumerate(cases[:1]):
        if dd == 0:
            wrf_data[case] = {}
        print(case,dom)

        tower_f = '{0}tower_netCDFs/NYSERDA_{1}_towers_{2}.nc'.format(main_directory,case,dom_str)

        sim_start = '2020-04-04 06:00:00'
        
        if dom >= 3:
            restarts = sorted(glob.glob('{}{}/RESTART_?'.format(main_directory,cases[0])))[3:-1]
            restarts = [rst.split('/')[-1] for rst in restarts]
            f_dir = '{}{}/'.format(main_directory,case)
            get_avg_profile = True
        else:
            restarts = ['RESTART_A']
            f_dir = '{}{}/'.format(main_directory,case)
            get_avg_profile = False

        if path.exists(tower_f):
            wrf_data[case][dom_str] = xr.open_dataset(tower_f)
        else:
            print('Reading tslist output...')
            tower_dat = tsout_seriesReader(fdir=f_dir,
                                           restarts=restarts,
                                           simulation_start_time=[sim_start]*len(restarts),
                                           domain_of_interest=dom_str,
                                           time_step=time_step,
                                           structure='unordered',
                                          )
            print('Saving to file: {}'.format(tower_f))
            tower_dat.to_netcdf(tower_f)
            wrf_data[case][dom_str] = tower_dat
            
        if get_avg_profile:
            avg_tower_f = tower_f.replace('.nc','_avg.nc') 
            if path.exists(avg_tower_f):
                print('Average already exists {}'.format(avg_tower_f))
            else:
                print('Taking average of tslist output')
                avg_twr = wrf_data[case][dom_str].mean(dim='station')
                avg_twr['lon'] = wrf_data[case][dom_str].lon.mean(dim='station')
                avg_twr['lat'] = wrf_data[case][dom_str].lat.mean(dim='station')
                avg_twr['zsurface'] = wrf_data[case][dom_str].zsurface.mean(dim='station')

                avg_twr = avg_twr.assign_coords({'station':'E06',
                                                 'lat':avg_twr.lat,
                                                 'lon':avg_twr.lon,
                                                 'zsurface':avg_twr.zsurface}).expand_dims({'station':1})

                avg_twr_dict[dom_str] = avg_twr
                print('Saving to file: {}'.format(avg_tower_f))
                avg_twr.to_netcdf(avg_tower_f)