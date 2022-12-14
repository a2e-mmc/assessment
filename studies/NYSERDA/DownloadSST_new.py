import os
import glob
import xarray as xr
import numpy as np

#file_names = ['MUR-JPL-L4-GLOB-v4.1','UKMO-L4HRfnd-GLOB-OSTIA','ABOM-L4LRfnd-GLOB-GAMSSA_28km','AVHRR_OI-NCEI-L4-GLOB-v2.0','CMC0.1deg-CMC-L4-GLOB-v3.0',
#              'Geo_Polar_Blended-OSPO-L4-GLOB-v1.0','JPL_OUROCEAN-L4UHfnd-GLOB-G1SST','K10_SST-NAVO-L4-GLOB-v01']

save_dir_base = '/glade/scratch/hawbecke/WRF/MMC/NYSERDA/SST/'

file_names = {'MUR-JPL-L4-GLOB-v4.1':'./MUR',
              'UKMO-L4HRfnd-GLOB-OSTIA':'./OSTIA',
              'ABOM-L4LRfnd-GLOB-GAMSSA_28km':'./GAMSSA',
              'AVHRR_OI-NCEI-L4-GLOB-v2.0':'./NCEI',
              'CMC0.1deg-CMC-L4-GLOB-v3.0':'./CMC',
              'Geo_Polar_Blended-OSPO-L4-GLOB-v1.0':'./OSPO',
              'JPL_OUROCEAN-L4UHfnd-GLOB-G1SST':'./G1SST',
              'K10_SST-NAVO-L4-GLOB-v01':'./NAVO'}
#file_names = ['UKMO-L4HRfnd-GLOB-OSTIA']
#file_names = ['AVHRR_OI-NCEI-L4-GLOB-v2.0']

#date_start_old = '20200405'
#date_end_old   = '20200408'
date_start = '2020-04-05T00:00:00Z'
date_end   = '2020-04-08T00:00:00Z'

wrf_dir = '/glade/work/hawbecke/MMC/NYSERDA/met_em/MERRA2/orig/'
met_files = glob.glob('{}met_em.d01*'.format(wrf_dir))

met = xr.open_dataset(met_files[0])
min_lon_orig = int(np.floor(float(np.min(met.XLONG_M))) - 2)
min_lat_orig = int(np.floor(float(np.min(met.XLAT_M))) - 2)
max_lon_orig = int(np.floor(float(np.max(met.XLONG_M))) + 2)
max_lat_orig = int(np.floor(float(np.max(met.XLAT_M))) + 2)

# Example call:
#./subset_dataset.py -s 20100101 -f 20100201 -b -140 -110 20 30 -x MUR-JPL-L4-GLOB-v4.1
for fn in file_names:
    save_dir = save_dir_base + file_names[fn]
    if 'NAVO' in fn:
        min_lat = max_lat_orig*-1
        max_lat = min_lat_orig*-1
    else:
        min_lat = min_lat_orig
        max_lat = max_lat_orig
    min_lon = min_lon_orig
    max_lon = max_lon_orig
    #old_cmd = 'subset_dataset.py -s {0} -f {1} -b {2} {3} {4} {5} -x {6}'.format(
    #                        date_start_old,date_end_old,min_lon,max_lon,min_lat,max_lat,fn)
    cmd = '~/local/envs/mmc/bin/podaac-data-downloader -b="{0},{1},{2},{3}" -c {4} -d {5} --start-date {6} --end-date {7}'.format(min_lon,min_lat,max_lon,max_lat,fn,save_dir,date_start,date_end)
    os.system(cmd)
