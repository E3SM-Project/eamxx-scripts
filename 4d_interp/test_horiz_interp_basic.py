#!/usr/bin/env python
"""
Test horizontal interpolation using ESMF/Tempest weights (i.e. horiz_interp2.py)
"""
#GET SET UP:
#===================
from netCDF4 import Dataset
import pylab as pl
from matplotlib.tri import Triangulation
import numpy as np

#LOAD REMAP FILE:
#===================
#downloaded from https://web.lcrc.anl.gov/public/e3sm/mapping/maps/
map_fi=Dataset('map_ne30np4_to_ne4np4_aave.20160601.nc')
#map_fi=Dataset('map_ne30np4_to_ne16np4_aave.20160601.nc')

wts=map_fi.variables['S'][:] #1D index of weights. Needs decoder indexing
col=map_fi.variables['col'][:] #Pointer to Source Grid Element
row=map_fi.variables['row'][:] #Pointer to Destination Grid Element
#note that row and col assume 1st element is 1, not 0!!!
#Use this info via: dest_data[row[i]] = src_data[col[i]]*S[i]

src_lats=map_fi.variables['yc_a'][:]
src_lons=map_fi.variables['xc_a'][:]
dest_lats=map_fi.variables['yc_b'][:]
dest_lons=map_fi.variables['xc_b'][:]

map_fi.close()
       
#GET A SINGLE TIME SAMPLE OF SPA DATA
#====================================
#All SPA data is in a 13G file called spa_file_unified_and_clipped.nc.
#A copy is on NERSC at /global/cscratch1/sd/beydoun/SPA_files/
#For testing horiz regridding, we instead use a file containing a single
#field (CCN) in raw ne30 format from this same directory.

src_fi=Dataset('SPA_ne30_only_CCN.nc')
src_var = src_fi.variables['PS'][0,:] #shape=[time,ncols]
src_fi.close()

#DO THE INTERPOLATION:
#======================
dest_var=np.zeros(dest_lats.shape)
for i in range(len(wts)):
    #subtract 1 from each index since python indices start at 0
    dest_var[row[i]-1] += wts[i]*src_var[col[i]-1]

#PLOT DATA ON OLD AND NEW GRID:
#============================
levels=np.arange(40000.,110000.,5000.)

pl.figure(1)
ax=pl.subplot(2,1,1)
tri_src = Triangulation(src_lons, src_lats)
ax.tricontourf(tri_src, src_var,levels=levels)
#pl.triplot(tri_src,linewidth=0.1) #plot grid
pl.title('PS on src grid')

ax=pl.subplot(2,1,2)
pl.title('Interpolated data')
tri_dest = Triangulation(dest_lons, dest_lats)
ax.tricontourf(tri_dest, dest_var,levels=levels)
#note that we're using delaunay triangulation to plot
#unstructured data and it doesn't know about wrapping
#from lon=360 to lon=0 so will have missing data at
#the edges of the plot. This reflects the plotting, not
#a problem with the interpolation method.
#pl.triplot(tri_dest,linewidth=0.2) #plot grid

pl.show()

