#!/usr/bin/env python
"""
Extend horiz_interp_basic.py under the hypothetical
case that we use multiple MPI ranks to break up the 
memory structure of our arrays. 

NOTES: 
1. I compute my_global_indices as a separate step here
because I will need to remember this info every month
when I read in a new data file.
2. there's an annoying semantic problem: "col" refers to
column of the remap matrix, but we also have atmospheric
columns to loop over. I handle this by calling the first "col"
and the 2nd atm_col.
3. I'm using lists to handle the fact that each atm col may 
involve a different number of remap indices. This works really
well but may not translate to SCREAM. Also, I use np.where 
instead of scanning/looping to look for appropriate indices. 
We will need to cook up an alternative impl for C++. 
4. I'm still a bit unclear on what masterproc should do vs 
individual procs. What I have works but passes global copies
of weights and indices to each proc in the init step; global 
copies don't need to persist outside init. We could instead 
pass the list of local atm_col indices to master and have it 
just scatter the bits of info each proc needs, but that would
be worse for load balancing and passing local index lists back
to masterproc would be annoying. 
"""
#GET SET UP:
#===================
from netCDF4 import Dataset
import pylab as pl
from matplotlib.tri import Triangulation
import numpy as np

#EMULATE MPI MEMORY MANAGEMENT:
#===================
masterproc=0
global_atm_ncols=866 #len(dest_lats); only needed b/c python doesn't know grid
nranks=3 #1, 2, 433 evenly divide ne4's global_ncols

dest_data=np.zeros(global_atm_ncols) #just needed for plotting in python

for mpi_rank in range(nranks):
    #below, each mpi_rank gets interspersed data for extra challenge
    local_indices=range(mpi_rank,global_atm_ncols,nranks) 
    my_atm_ncols=len(local_indices)
        
    if 0==0:  #if masterproc:
    #the following loop should only be executed on masterproc, but python is
    #saying it can't find row and col when they're conditionally defined.
           
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

        #just used for python plotting
        src_lats=map_fi.variables['yc_a'][:]
        src_lons=map_fi.variables['xc_a'][:]
        dest_lats=map_fi.variables['yc_b'][:]
        dest_lons=map_fi.variables['xc_b'][:]

        map_fi.close()

        #SHIFT ALL THE ROW + COL INDICES DOWN 1 FOR PYTHON INDEXING
        #====================================
        row=row-1
        col=col-1

        #GET A SINGLE TIME SAMPLE OF SPA DATA
        #====================================
        #in SCREAM, we will need to load a different array every month
        #so this step will need to be included in a "if new month"
        #conditional. Not sure if just masterproc should load it then
        #broadcast, or whether we can pass masterproc my_global_indices
        #defined below and it can scatter just the data each mpi rank needs.
        src_fi=Dataset('SPA_ne30_only_CCN.nc')
        src_data = src_fi.variables['PS'][0,:] #shape=[time,ncols]
        src_fi.close()
        
        #would mpi_broadcast wts, col, row, and src_data here.
        #lat and lon info are just for python plotting
        
    #GET LOCAL INDICES
    #========================
    #get "my_global_indices" for each atm_col = the indices of the global wts,
    #row, and col arrays which that atm_col needs to know about. Need to
    #keep my_global_indices around for the whole run so we can use it to
    #extract new data at the start of each new month.
    
    #using list datatype for my_rows because each atm_col on this mpi_rank 
    #could involve a different number of indices in wts
    my_global_indices=[] 
    for ind in local_indices:
        #np.where returns all the indices of wts array
        #that involve the index for the destination column
        #"ind" on this mpi rank
        my_global_indices.append(np.where(row==ind)) 

    #EXTRACT ENTRIES OF wts, col, and src_data NEEDED BY EACH ATM_COL
    #========================
    #this is the only info that needs to be stored in each mpi rank
    #persistently across timesteps
    my_wts=[]
    my_cols=[]
    my_data=[]
    for i in range(my_atm_ncols):
        my_wts.append(wts[my_global_indices[i]])
        my_cols.append(col[my_global_indices[i]])
        #note the double redirection in the line below.
        #A bit hacky but works. 
        my_data.append(src_data[col[my_global_indices[i]]])

    #DO INTERPOLATION
    #=======================
    my_dest_data=np.zeros(my_atm_ncols)
    for i in range(my_atm_ncols):
        my_dest_data[i]=np.sum( my_wts[i]*my_data[i] )

    #STITCH DEST_DATA BACK INTO A GLOBAL ARRAY
    #=======================
    #this isn't needed for C++, just for plotting results in python.
    for i in range(my_atm_ncols):
        dest_data[local_indices[i]]=my_dest_data[i]
   
#PLOT DATA ON OLD AND NEW GRID:
#============================
levels=np.arange(40000.,110000.,5000.)

pl.figure()
ax=pl.subplot(2,1,1)
tri_src = Triangulation(src_lons, src_lats)
ax.tricontourf(tri_src, src_data,levels=levels)
#pl.triplot(tri_src,linewidth=0.1) #plot grid
pl.title('PS on src grid')

ax=pl.subplot(2,1,2)
pl.title('Interpolated data')
tri_dest = Triangulation(dest_lons, dest_lats)
ax.tricontourf(tri_dest, dest_data,levels=levels)
#note that we're using delaunay triangulation to plot
#unstructured data and it doesn't know about wrapping
#from lon=360 to lon=0 so will have missing data at
#the edges of the plot. This reflects the plotting, not
#a problem with the interpolation method.
#pl.triplot(tri_dest,linewidth=0.2) #plot grid

pl.show()

