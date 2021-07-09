#!/usr/bin/env python
"""
Python version of 4d bilinear interpolation
for SPA. Meant as a test code.

NOTES:
1. horizontal_interpolate takes in cell-center lats and lons and 
constructs cell edges by assuming they're halfway between cell centers.
We actually know cell corners, but this gives us different lats and lons 
for each of the 4 corners, which greatly complicates things. 
"""
#GET SET UP:
#===================
from netCDF4 import Dataset
import pylab as pl
import numpy as np
import horiz_interp
p0=1e5 #surf pressure constant for hybrid coordinates
var='CCN3' #just choosing one of the included vars at random

model_fi=Dataset('/Users/caldwell19/junk/20210708.master.F2000-SCREAM-LR.ne4pg2_ne4pg2.cori-knl.2x64x1.eam.h0.0001-01-02-00000.nc') #just grabbing a random ne4 file to get its grid info
model_lats=model_fi.variables['lat'][:]
model_lons=model_fi.variables['lon'][:]
#model_var=model_fi.variables['TMQ'][:] #prob won't use this var, just want its size.

#This version gives all 4 edges for each cell, which doesn't work with horiz_interp 
#model_fi=Dataset('~/junk/ne30pg2.nc') #from /global/homes/z/zender/data/grids/ne30pg2.nc on NERSC.
#model_lat_edges=model_fi.variables['grid_corner_lat'] #ncols x 4 corners
#model_lon_edges=model_fi.variables['grid_corner_lon'] #ncols x 4 corners

model_fi.close()

#STUFF TO DO WHEN RUN *STARTS*:
#===================

#GET FIRST 2 MONTHS OF DATA
#--------------------------------
#SPA data is a 13G file called spa_file_unified_and_clipped.nc.
#A copy is on NERSC at /global/cscratch1/sd/beydoun/SPA_files/

fi=Dataset('spa_file_unified_and_clipped.nc')
file_lats=fi.variables['lat'][:]
file_lons=fi.variables['lon'][:]

times = fi.variables['time']
file_time_old = times[0]
file_time_new = times[1]

file_var_old = fi.variables[var][0,:,:] #assuming time,lon,lat dims
file_var_new = fi.variables[var][0,:,:]

file_ps_old = fi.variables['PS'][0,:,:]
file_ps_new = fi.variables['PS'][0,:,:]

#since fi contains data for all months, I won't bother closing it.

#GET WEIGHTS FOR HORIZONTAL INTERP:
#----------------------------------
#On MPI, model_lats and model_lons would be local to the MPI proc, but
#file_lats and file_lons would need to be global.
weight_lats,weight_lons=horiz_interp.xy_interp_init(file_lats,file_lons,file_lats,file_lons)
count_lats,index_lats,count_lons,index_lons=horiz_interp.weights2indices(weight_lats,weight_lons)

#ANOTHER OPTION WOULD BE TO GRAB THE WEIGHTS FROM AN EXISTING MAP FILE.
#map_file = Dataset('map_file.nc')
#horiz_wts = file.variables['map'] #can we read this in an MPI way which just loads the weights needed by columns on a particular node?
#map_file.close()

#DO HORIZONTAL INTERPOLATION:
#============================
dest_data=xy_interp(weight_lats,count_lats,index_lats,weight_lons,count_lons,index_lons,file_ps_old)

#PLOT DATA ON OLD AND NEW GRID:
#============================
pl.figure(1)
pl.subplot(2,1,1)
LONS,LATS=pl.meshgrid(file_lons,file_lats)
pl.pcolor(LONS,LATS,file_ps_old)
pl.colorbar()
pl.title('PS on src grid')

pl.subplot(2,1,2)
pl.scatter(model_lons,model_lats,dest_data)
pl.colorbar()
pl.title('Interpolated data')

pl.show()

