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

# # # FOR TESTING # # #
cases = [cases[0]] + cases[3:6] # + [cases[-1]]
for cc,case in enumerate(cases): print(cc,case)
# # # # # # # # # # # #

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

def averageTslistData(data):
    avg_twr = data.mean(dim='station')
    avg_twr['lon'] = data.lon.mean(dim='station')
    avg_twr['lat'] = data.lat.mean(dim='station')
    avg_twr['zsurface'] = data.zsurface.mean(dim='station')

    avg_twr = avg_twr.assign_coords({'station':'E06',
                                     'lat':avg_twr.lat,
                                     'lon':avg_twr.lon,
                                     'zsurface':avg_twr.zsurface}).expand_dims({'station':1})
    return(avg_twr)


avg_twr_dict = {}

for dd,dom in enumerate(doms_of_interest):
    
    time_step = dom_dict[dom]['dt']
    dom_str = 'd0{}'.format(dom)

    for cc,case in enumerate(cases):
        tower_f = '{0}tower_netCDFs/NYSERDA_{1}_towers_{2}.nc'.format(main_directory,case,dom_str)
        avg_tower_f = tower_f.replace('.nc','_avg.nc') 

        sim_start = '2020-04-04 06:00:00'
        
        if dom >= 3:
            restarts = sorted(glob.glob('{}{}/RESTART_*'.format(main_directory,case)))[3:]
            restarts = [rst.split('/')[-1] for rst in restarts]
            # Get rid of duplicate restart directories... the last one is always good:
            dup_restarts = []
            for rst in restarts:
                if len(rst) > len('RESTART_X'):
                    dup_restarts += [rst]
            head_restarts = []
            for rst in dup_restarts:
                head_rst = rst[:-2]
                if head_rst not in head_restarts:
                    head_restarts += [head_rst]
            for head_rst in head_restarts:
                dup_rsts = []
                for dup_rst in dup_restarts:
                    if head_rst in dup_rst: dup_rsts += [dup_rst]
                dup_rsts = [head_rst] + dup_rsts
                dup_rsts = sorted(dup_rsts)
                keep_rst = dup_rsts[-1]
                dup_rsts = dup_rsts[:-1]
                add_ind = restarts.index(dup_rsts[0])
                for rst in dup_rsts:
                    restarts.remove(rst)
                restarts = restarts[:add_ind] + [keep_rst] + restarts[add_ind+1:]
            f_dir = '{}{}/'.format(main_directory,case)
            get_avg_profile = True
        else:
            restarts = ['RESTART_A']
            f_dir = '{}{}/'.format(main_directory,case)
            get_avg_profile = False

        if path.exists(tower_f):
            tower_dat = xr.open_dataset(tower_f)

        else:
            print('Reading tslist output...')
            if dom != 5:
                tower_dat = tsout_seriesReader(fdir=f_dir,
                                               restarts=restarts,
                                               simulation_start_time=[sim_start]*len(restarts),
                                               domain_of_interest=dom_str,
                                               time_step=time_step,
                                               structure='unordered',
                                              )
                
                print('Saving to file: {}'.format(tower_f))
                tower_dat.to_netcdf(tower_f)
            else:
                get_avg_profile = False
                chunked_list = list()
                chunk_size = 2
                for i in range(0, len(restarts), chunk_size):
                    chunked_list.append(restarts[i:i+chunk_size])

                for cc,chunk in enumerate(chunked_list):
                    chunk_f = tower_f.replace('.nc','_chunk{}.nc'.format(cc))
                    if not path.exists(chunk_f):
                        tower_dat = tsout_seriesReader(fdir=f_dir,
                                                       restarts=chunk,
                                                       simulation_start_time=[sim_start]*len(chunk),
                                                       domain_of_interest=dom_str,
                                                       time_step=time_step,
                                                       structure='unordered',
                                                      )
                        print('Saving to file: {}'.format(chunk_f))
                        tower_dat.to_netcdf(chunk_f)
                    else:
                        tower_dat = xr.open_dataset(chunk_f)

                    grid_stns = []
                    for stn in tower_dat.station.data:
                        if stn[0] == 'T':
                            grid_stns += [str(stn)]

                    avg_twr = averageTslistData(tower_dat.sel(station=grid_stns))
                    del(tower_dat)
                    if cc == 0:
                        avg_twr_f = avg_twr
                    else:
                        avg_twr_f = xr.merge([avg_twr_f,avg_twr])
                if not path.exists(avg_tower_f):
                    print('Saving to file: {}'.format(avg_tower_f))
                    avg_twr_f.to_netcdf(avg_tower_f)
            

        if get_avg_profile:
            if path.exists(avg_tower_f):
                print('Average already exists {}'.format(avg_tower_f))
            else:
                print('Taking average of tslist output')
                grid_stns = []
                for stn in tower_dat.station.data:
                    if stn[0] == 'T':
                        grid_stns += [str(stn)]
                avg_twr = averageTslistData(tower_dat.sel(station=grid_stns))

                avg_twr_dict[dom_str] = avg_twr
                print('Saving to file: {}'.format(avg_tower_f))
                avg_twr.to_netcdf(avg_tower_f)
print('Finished.')
