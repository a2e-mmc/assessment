#!/usr/bin/env python
# coding: utf-8

# # Read_WRF_slices
# Use mmctools wrfseries_slice_reader to read FINO offshore microscale simulations.<br>
# Save dataseries to pdata file

# ## Set file paths, import 'sys' module for setting filepath

# In[7]:


#wrf_dir = '/p/lscratchh/lassman1/a2e-mmc/WRF/FINO/ERA5_grid2.0/'
#cases = ['microscale_setup_v2_CONTROL', 'microscale_setup_v2_CPM']

#wrf_dir = '/p/lustre1/lassman1/a2e-mmc/FINO/WRF/resims/' 
#cases = ['microscale_setup_v2_CONTROL', 'microscale_setup_v2_CPM']

wrf_dir = '/p/lustre1/lassman1/a2e-mmc/FINO/WRF/resims/'
cases = ['CONTROL', 'CPM']



pdata_save_dir = '/p/lustre1/lassman1/a2e-mmc/FINO/pdata_resim/'



import sys
import os

#module_path = os.path.join(os.environ['HOME'],'code/Python/a2e-mmc/')
module_path = os.path.join(os.environ['HOME'],'mmc/mmc_github_clones/mmctools/')                                                          
if module_path not in sys.path:
    sys.path.append(module_path)


# In[8]:


# WRF File name template
# Done in 20-minute incriments because the model ran with 20 minute restarts, and some files had different numbers of slices
# Also this will avoid a memory error (hopefully)
list_file_filters = [#'mmc_d06_2010-05-16_01:[01][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_01:[23][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_01:[45][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_02:[01][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_02:[23][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_02:[45][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_03:[01][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_03:[23][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_03:[45][0-9]:[0-5][0-9]', \
                     'mmc_d06_2010-05-16_04:[01][0-9]:[0-5][0-9]', \
                     'mmc_d06_2010-05-16_04:[23][0-9]:[0-5][0-9]', \
                     'mmc_d06_2010-05-16_04:[45][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_05:[01][0-9]:[0-5][0-9]', \
                     #'mmc_d06_2010-05-16_05:[23][0-9]:[0-5][0-9]', \
             ]


hgt = 80


# ## Import the rest of the modules

# In[9]:


# import modules
import numpy as np
from netCDF4 import Dataset as nc
import xarray as xr


import pandas as pd
import datetime

import wrf as wrfpy
import glob


# In[10]:


#from mmctools.plotting import plot_timehistory_at_height, plot_profile, plot_spectrum
from mmctools.helper_functions import calc_wind, theta, model4D_calcQOIs
from mmctools.wrf.utils import wrfout_seriesReader,wrfout_slices_seriesReader


# In[11]:


def read_range(  current_ind, 
               list_file_filters = list_file_filters,
               case_id = 0,
               n_read = 3, 
               save_sfx  = 'all_times', 
               save_data = False ):
        print("Starting seriesReader for first slice")
        ds_slice =  wrfout_slices_seriesReader( wrf_dir + cases[case_id] + '/auxout/', list_file_filters[current_ind], specified_height = [hgt] )
    
        ct = current_ind
        for ff in list_file_filters[current_ind+1:current_ind + n_read]:
            print(ct, ' of ', len(list_file_filters))
            ds_temp = wrfout_slices_seriesReader( wrf_dir + cases[case_id] + '/auxout/', ff, specified_height = [hgt] )
            ds_slice = xr.concat([ ds_slice, ds_temp], dim = 'datetime')
            ct += 1    
        print(ds_slice.datetime.shape)
        if save_data:
            save_str = pdata_save_dir + cases[case_id] + 'slice_out_'+str(hgt) + save_sfx + '.nc'
            print("saving...", save_str)
            ds_slice.to_netcdf( save_str )
        print("Finished")
        return current_ind, ds_slice
        
    


# In[ ]:


print("HR1")
current_ind, slice_hr1 = read_range( 0, case_id = 1, save_sfx = '_HR4', save_data = True)
#print("HR2")
#current_ind = read_range( current_ind, save_sfx = '_HR2')
#print("HR3")
#current_ind = read_range( current_ind, save_sfx = '_HR3')
#print("HR4")
#current_ind = read_range( current_ind, save_sfx = '_HR4')
import sys
sys.exit()

# In[ ]:


ds_slice2_control = wrfout_slices_seriesReader( wrf_dir + cases[0] + '/auxout/', list_file_filters[0], specified_height = [hgt] )
ds_slice2_cpm = wrfout_slices_seriesReader( wrf_dir + cases[1] + '/auxout/', list_file_filters[0], specified_height = [hgt] )


# In[ ]:


ds_slice4_control = wrfout_slices_seriesReader( wrf_dir + cases[0] + '/auxout/', list_file_filters[1], specified_height = [hgt] )
ds_slice4_cpm = wrfout_slices_seriesReader( wrf_dir + cases[1] + '/auxout/', list_file_filters[1], specified_height = [hgt] )

#ds40 = xr.DataArray()
ct = 1
for ff in list_file_filters:
    print(ct, ' of ', len(list_file_filters))
    ds_temp = wrfout_slices_seriesReader( wrf_dir + cases[0] + '/auxout/', ff, specified_height = [hgt] )
    ds_slice = xr.concat([ ds_slice, ds_temp], dim = 'datetime')
    #ct += 1ds_slice.to_netcdf( pdata_save_dir + cases[0] + 'slice_out_'+str(hgt)+'.nc')
# In[ ]:


ds_slice2_cpm.coords


# In[ ]:


from matplotlib import pyplot as plt

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,11))
axs = axs.flatten()
 
vmin=10
vmax=18
t_hour_0 = '2010-05-16 02:00:00'
t_hour_1 = '2010-05-16 04:00:00'
 
# --------- SOWFA ------------
currds = ds_slice2_cpm.isel(datetime=1)
#print(currds.u.values)
#print( (currds.u*currds.u).values )
#print( currds.v.values  )
#print( ( currds.v*currds.v).values  )
#print( np.sqrt( currds.u*currds.u + currds.v * currds.v ).values )
currds['wspd'] = np.sqrt( currds.u * currds.u + currds.v * currds.v )[0,:,:]
print(currds.wspd.shape)
p = axs[0].pcolormesh(currds.x, currds.y, currds.wspd, vmin=vmin, vmax=vmax, cmap='viridis', shading='auto')
axs[0].set_title('WRF-CPM 02:00 UTC')
 
# --------- WRF--CPM ------------
currds = ds_slice2_control.isel(datetime=1)
currds['wspd'] = np.sqrt( currds.u*currds.u + currds.v * currds.v )[0,:,:]
p = axs[1].pcolormesh(currds.x, currds.y, currds.wspd, vmin=vmin, vmax=vmax, cmap='viridis', shading='auto')
axs[1].set_title('WRF-CONTROL 02:00 UTC')
 
# --------- WRF Synthetic Turb ------------
currds = ds_slice4_cpm.isel(datetime=1)
currds['wspd'] = np.sqrt( currds.u*currds.u + currds.v * currds.v )[0,:,:]
p = axs[2].pcolormesh(currds.x, currds.y, currds.wspd, vmin=vmin, vmax=vmax, cmap='viridis', shading='auto')
axs[2].set_title('WRF-CPM 04:00 UTC')
 
# --------- WRF No pertubations ------------
currds = ds_slice4_control.isel(datetime=1)
currds['wspd'] = np.sqrt( currds.u*currds.u + currds.v * currds.v )[0,:,:]
p = axs[3].pcolormesh(currds.x, currds.y, currds.wspd, vmin=vmin, vmax=vmax, cmap='viridis', shading='auto')
axs[3].set_title('WRF-CONTROL 04:00 UTC')
 
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.3, 0.02, 0.4])
cb = fig.colorbar(p, cax=cbar_ax)
cb.set_label('wind speed [m/s]')
 
for ax in axs:
    ax.set_xlim([0,6000]); ax.set_ylim([0,6000])
    ax.set_aspect('equal', 'box')
    ax.axis('off')
 
plt.savefig('WRF_slice_pcolormesh.png', dpi = 300 )


# In[ ]:


print(sys.path)


# In[ ]:




