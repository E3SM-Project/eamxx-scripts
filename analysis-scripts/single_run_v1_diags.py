#!/usr/bin/env python
"""
Explore first v1 output
"""

from netCDF4 import MFDataset, Dataset
import numpy as np
import pylab as pl

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def water_path(q,p_int):
    """
    Compute the liquid, ice, or vapor water path given a 3d variable q (which can 
    be qc, qi, or qv, respectively). Output units are g/m2.
    """
    g=9.8 #gravity
    dp=p_int[:,:,0:-1] - p_int[:,:,1:] #Note 0 is TOA and -1 is surf so dp should be neg.

    wp=-1000.*np.sum(q*dp/g,axis=-1)

    return wp
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#SAVE ALL FIGS AS 1 PDF:
#=========================
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages("/global/cfs/cdirs/e3sm/www/terai/SCREAM/v1_analysis/ne30_v1_PM_quickcrash_first2timesteps.pdf")

#GET MAPPING TOOLS:
#=========================
map_fi='/global/homes/z/zender/data/maps/map_ne30np4_to_fv128x256_aave.20160301.nc'
#map_fi='/global/homes/z/zender/data/maps/map_ne4np4_to_fv25x48_aave.20170401.nc'
#map_fi='/usr/gdata/climdat/maps/map_ne4np4_to_fv25x48_aave.20170401.nc'
#map_fi='/usr/gdata/climdat/maps/map_ne4np4_to_fv25x48_bilin.20170401.nc'
map_f=Dataset(map_fi)

#area-weighted ave = np.sum(var[0:ncol]*areas)/np.sum(areas)
areas=map_f['area_a'][:]

#lat and lon at cell centers
lats=map_f['yc_a'][:]
lons=map_f['xc_a'][:]

lats_edges=map_f['yv_a'][:] #dims=[len(data),5] since corners can have 5 vertices?
lons_edges=map_f['xv_a'][:]

map_f.close()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def plot_2d(data,units):
  
  pl.figure(figsize=[8.5,11])

  #GLOBAL AVE TIMESERIES:
  #=========================
  ave_data=np.dot(data,np.reshape(areas,[len(areas),1])).flatten()/np.sum(areas)

  pl.subplot(3,1,1)
  pl.plot(ti,ave_data)
  pl.title(var+' (x axis=days into run, y axis=global-ave, units=%s)'%(units))
  #pl.xlabel('time (days into run)')
  #pl.ylabel(var)
  
  #MAPS:
  #=========================
  #use delaunay triangulation to plot on native grid

  #use same color axis for both init and final times:
  data_a_prctile=np.percentile(data[0,:],[1,10,90,99])
  data_b_prctile=np.percentile(data[-1,:],[1,10,90,99])

  mn=min(data_a_prctile[0],data_b_prctile[0])
  mx=max(data_a_prctile[3],data_b_prctile[3])
  colormap='viridis'
  if np.sign(mx)+np.sign(mn)<1:
      colormap='RdBu_r'
      max_mag=max(-1*mn,mx)
      mn=-1*max_mag
      mx=max_mag
  pl.subplot(3,1,2)
  pl.tricontourf(lons,lats,data[0,:],levels=
                 np.linspace(mn,mx,11),cmap=colormap,extend='both')
  pl.colorbar()
  #pl.xlabel('lon')
  pl.ylabel('lat')
  pl.title('init time')
  
  pl.subplot(3,1,3)
  pl.tricontourf(lons,lats,data[-1,:],levels=
                 np.linspace(mn,mx,11),cmap=colormap,extend='both')
  pl.colorbar()
  pl.xlabel('lon')
  pl.ylabel('lat')
  pl.title('final time')

  pl.show()
  
  pdf.savefig(pl.gcf().number)

  return
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def plot_3d(data,units):

    hts_len=data.shape[2]
    pl.figure(figsize=[8.5,11])

    #TIME VS HEIGHT GLOBAL AVE BY LEVEL
    #==============================
    ave_data=np.zeros([len(ti),hts_len])
    for i in range(hts_len):
      ave_data[:,i]=np.dot(data[:,:,i],np.reshape(areas,[len(areas),1])).flatten()\
                     /np.sum(areas)

    pl.subplot(3,2,1)
    HTS,TI=pl.meshgrid(np.arange(hts_len),ti)
    pl.pcolor(TI,HTS,ave_data[:,-1::-1]) #note flip so surface at bottom of plot
    pl.ylabel('model level (surf at bottom)')
    pl.title(var+', time vs ht of glob ave')
    pl.colorbar()

    #LAT VS HT AVE OVER ALL TIMES
    #=============================
    time_ave_data=np.average(data,axis=0)

    lat_bnds=np.arange(-90,91.,10.)
    lat_bin_centers=(lat_bnds[0:-1]+lat_bnds[1:])/2.

    zonal_aves=np.zeros([len(lat_bin_centers),hts_len])
    for i in range(len(lat_bin_centers)):
      msk=np.logical_and(np.greater_equal(lats,lat_bnds[i]),\
                         np.less(lats,lat_bnds[i+1])).astype(int)
      #dims are [1,len(lats)] x [len(lats),hts_len]
      zonal_aves[i,:]=(np.dot(np.reshape(msk*areas,[1,len(areas)]),time_ave_data)\
                       /np.sum(msk*areas)).flatten()

    pl.subplot(3,2,2)
    HTS,LATS=pl.meshgrid(np.arange(hts_len),lat_bin_centers)
    pl.pcolor(LATS,HTS,zonal_aves[:,-1::-1]) #flip so surf at bottom
    #pl.xlabel('lat')
    #pl.ylabel('ht')
    pl.title('time-ave lat vs ht, units=%s'%(units))
    pl.colorbar()

    #LOWEST-LAYER MAPS:
    #=========================
    #use delaunay triangulation to plot on native grid

    #use same color axis for both init and final times:
    data_a_prctile=np.percentile(data[0,:,-1],[1,10,90,99])
    data_b_prctile=np.percentile(data[-1,:,-1],[1,10,90,99])

    mn=min(data_a_prctile[0],data_b_prctile[0])
    mx=max(data_a_prctile[3],data_b_prctile[3])

    if mn==mx:
        mx=mn+0.01

    colormap='viridis'
    if np.sign(mx)+np.sign(mn)<1:
        colormap='RdBu_r'
        max_mag=max(-1*mn,mx)
        mn=-1*max_mag
        mx=max_mag

    if mn==mx:
        print('Min and max for var=%s were identical, equal to %f\n'%(var,mn))
    else:
        pl.subplot(3,2,3)
        pl.tricontourf(lons,lats,data[0,:,-1],levels=np.linspace(mn,mx,11),cmap=colormap)
        pl.colorbar()
        #pl.xlabel('lon')
        pl.ylabel('lat')
        pl.title('lowest layer, init time')

        pl.subplot(3,2,5)
        pl.tricontourf(lons,lats,data[-1,:,-1],levels=np.linspace(mn,mx,11),cmap=colormap)
        pl.colorbar()
        pl.xlabel('lon')
        pl.ylabel('lat')
        pl.title('lowest layer, final time')

    #MIDPT LAYER MAPS:
    #=========================
    #use delaunay triangulation to plot on native grid

    #find layer with most extreme value:
    biggest=np.zeros(hts_len)
    vert_ind=0
    for i in range(hts_len):
      biggest[i]=max(np.max(np.abs(data[0,:,i])),np.max(np.abs(data[-1,:,i])))
      if biggest[i]>biggest[vert_ind]:
        vert_ind=i
    vert_ind=data.shape[2]-2

    #use same color axis for both init and final times:
    data_a_prctile=np.percentile(data[0,:,vert_ind],[1,10,90,99])
    data_b_prctile=np.percentile(data[-1,:,vert_ind],[1,10,90,99])

    mn=min(data_a_prctile[0],data_b_prctile[0])
    mx=max(data_a_prctile[3],data_b_prctile[3])
    if mn==mx:
        mx=mn+0.01

    colormap='viridis'
    if np.sign(mx)+np.sign(mn)<1:
        colormap='RdBu_r'
        max_mag=max(-1*mn,mx)
        mn=-1*max_mag
        mx=max_mag
        
    if mn==mx:
        print('Min and max for var=%s, vert_ind=%i were identical, equal to %f\n'%(var,vert_ind,mn))
    else: 
        pl.subplot(3,2,4)
        pl.tricontourf(lons,lats,data[0,:,vert_ind],levels=np.linspace(mn,mx,11),cmap=colormap)
        pl.colorbar()
        #pl.xlabel('lon')
        #pl.ylabel('lat')
        pl.title('Layer %i, init time'%(vert_ind))

        pl.subplot(3,2,6)
        pl.tricontourf(lons,lats,data[-1,:,vert_ind],levels=np.linspace(mn,mx,11),cmap=colormap)
        pl.colorbar()
        pl.xlabel('lon')
        #pl.ylabel('lat')
        pl.title('Layer %i, final time'%(vert_ind))


    pdf.savefig(pl.gcf().number)

    return
#=-=-=-=-=-=-=-=-=-=-=-=

#GET ACTUAL DATA:
#==========================

fi_dir='/global/cfs/cdirs/e3sm/terai/SCREAM/v1_analysis/run4/'
f=MFDataset(fi_dir+'*INSTANT*.nc')

#time is fractional days since start of simulation
ti=f['time'][:] #days since start of simulation
#edge pressures
p_int=f.variables['p_int'][:]

vars=f.variables.keys()

for var in vars: #['surf_latent_flux','surf_sens_flux','omega','qc','T_mid','qv']: #vars:
  print('Handling '+var)  
  data=f[var][:]
  units=f[var].units

  if var=="surf_latent_flux":
    data=data*2.5e6; #convert kg w/(m2*s) to W/m2
    units='W/m2'

  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  #FOR 2D VARIABLES:
  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  if data.shape==(len(ti),len(lats)):
    plot_2d(data,units)


  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  #FOR 3D VARIABLES:
  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  #handles both lev and ilev cases unlike: (len(ti),len(lats),len(levs)):
  if len(data.shape)==3: 
    plot_3d(data,units)

  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  #FOR 3D VARIABLES:
  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  #handles both lev and ilev cases unlike: (len(ti),len(lats),len(levs)):
  if len(data.shape)==4:
    data_u=np.squeeze(data[:,:,0,:])
    plot_3d(data_u,'U '+units)
    data_v=np.squeeze(data[:,:,1,:])
    plot_3d(data_v,'V '+units)

#LOOP OVER DERIVED VARIABLES:
#=================

"""

var_dict={'LWP':'qc','IWP':'qi','WVP':'qv'}
for var in ['LWP','IWP','WVP']:
  print('Handling '+var)  
  data=f.variables[var_dict[var]][:]
  data_vint=water_path(data,p_int)
  plot_2d(data_vint,'g/m2')
  
var='TBOT'
print('Handling '+var)  
data=f.variables['T_mid'][:,:,-1]
units=f.variables['T_mid'].units
plot_2d(data,units)

var='QBOT'
print('Handling '+var)  
data=f.variables['qv'][:,:,-1]
units=f.variables['qv'].units
plot_2d(data,units)

var='Theta'
print('Handling '+var)
T=f.variables['T_mid'][:]
p=f.variables['p_mid'][:]
data=T*(100000./p)**(287.1/1004.)
units=f.variables['T_mid'].units
plot_3d(data,units)

"""

#CLEAN UP/FINALIZE
#=================
pdf.close()
#pl.show()
print('*** DONE ***')
