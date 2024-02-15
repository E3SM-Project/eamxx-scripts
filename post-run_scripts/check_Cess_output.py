#!/usr/bin/env latest_e3sm_unified 

#libraries to be used by routines
import xarray
import subprocess


# Load SCREAMv1 file outputs and check if open_mfdataset works.
# If it throws an error, either there's missing data, files with 0 timesteps, or overlapping data.

directory_prefix='/global/cfs/cdirs/e3smdata/simulations/Cess/cess-plus4k.ne1024pg2_ne1024pg2.F2010-SCREAMv1.cess-oct2/run/'
output_streams=['3hourlyAVG_ne120.AVERAGE','3hourlyINST_ne1024.INSTANT',
                '3hourlyINST_ne120.INSTANT','50hourly_QcQiNcNi.INSTANT',
                '50hourly_QrNrQmBm.INSTANT','6hourlyAVG_ne30.AVERAGE',
                '6hourlyINST_ne30.INSTANT','ACIregions_2D.INSTANT',
                'ARMSite_2D.INSTANT','ARMSite_3D.INSTANT',
                'hourly2DVars.INSTANT','monthly_COSP_ne1024.AVERAGE',
                'monthly_ne1024.AVERAGE']

time_periods=['2019-08','2019-09','2019-10','2019-11','2019-12','2020-01','2020-02','2020-03','2020-04','2020-05','2020-06','2020-07','2020-08']

for i in time_periods:
    #print(i)
    for j in output_streams:
        #print(j)
        try:
            f_check=xarray.open_mfdataset(directory_prefix+'/output.scream.Cess.'+j+'*'+i+'*nc')
            print('Timesteps good with '+i+' '+j)
        except:
            print('Error with '+i+' '+j)
