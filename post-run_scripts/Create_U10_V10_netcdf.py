#Import necessary libraries                                            
# Make sure to run the following to load necessary libraries:                           
# source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh       

import requests
import xarray
import numpy as np
import scipy.sparse
import sys


input_file_directory='/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SCREAMv0.SCREAM-DY2.ne1024pg2.20201127/run/'
input_h1_prefix='SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h1.2020-'
input_h5_prefix='SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h5.2020-'
input_file_suffix='-00000.nc'

time_slice=np.arange(0,96,12)

#days=['01-20','01-21','01-22','01-23','01-24','01-25','01-26','01-27','01-28','01-29','01-30','01-31','02-01','02-02','02-03','02-04','02-05','02-06','02-07','02-08','02-09','02-10','02-11','02-12','02-13','02-14','02-15','02-16','02-17','02-18','02-19',]
days=['02-20','02-21','02-22','02-23','02-24','02-25','02-26','02-27','02-28','03-01',]

for k in days:
    print(k)
    h1_input_filename=''.join([input_file_directory,input_h1_prefix,k,input_file_suffix])
    h5_input_filename=''.join([input_file_directory,input_h5_prefix,k,input_file_suffix])

    ds_h1 = xarray.open_dataset( h1_input_filename )
    ds_h5 = xarray.open_dataset( h5_input_filename )

    # Load variables
    WINDSPD_10M_3hrly = ds_h1['WINDSPD_10M'].isel(time=time_slice)

    U_bot=ds_h5['U'].isel(lev=127)
    V_bot=ds_h5['V'].isel(lev=127)

    U_bot_array=np.array(U_bot)
    V_bot_array=np.array(V_bot)

    ratio_vOu=V_bot/U_bot
    ratio_uOv=U_bot/V_bot
    ratio_vOu=np.nan_to_num(ratio_vOu,nan=999)  # remove nans where u goes to 0
    ratio_uOv=np.nan_to_num(ratio_uOv,nan=999)  # same for v

    u10=np.sign(np.array(U_bot)) * WINDSPD_10M_3hrly / (1+ratio_vOu**2)**0.5 #Calculate u10
    v10=np.sign(np.array(V_bot)) * WINDSPD_10M_3hrly / (1+ratio_uOv**2)**0.5

    #Add attributes to data
    u10.name='U_10M'
    u10.attrs["units"]='m/s'
    u10.attrs["long_name"]='10-m zonal wind speed'
    u10.attrs["description"]='Derived from WINDSPD_10M and bottom model level U'

    v10.name='V_10M'
    v10.attrs["units"]='m/s'
    v10.attrs["long_name"]='10-m meridional wind speed'
    v10.attrs["description"]='Derived from WINDSPD_10M and bottom model level V'

    output_file_name="".join([input_file_directory,"Additional_output/","SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.U_10M.2020-",k,"-00000.nc"])
    u10.to_netcdf(output_file_name)
    output_file_name="".join([input_file_directory,"Additional_output/","SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.V_10M.2020-",k,"-00000.nc"])
    v10.to_netcdf(output_file_name)





    
