#!/usr/bin/env python
"""
Helper functions for diagnostic routines
"""

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

def plot_2d(ti,data,var,units,lats,lons,areas,outfi):

    pl.figure(figsize=[8.5,11])

    #GLOBAL AVE TIMESERIES:
    #=========================
    ave_data=np.dot(data,np.reshape(areas,[len(areas),1])).flatten()/np.sum(areas)
    max_data=np.max(data,axis=1)
    min_data=np.min(data,axis=1)

    pl.subplot(3,1,1)
    pl.plot(ti,ave_data)
    pl.plot(ti,max_data,'r--')
    pl.plot(ti,min_data,'g-.')
    pl.legend(['ave','max','min'],loc='best',fontsize='small')
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

    #if 1% and 99% vals are identical, try min to max.
    #note this means some plots may show min-to-max and others 1% to 99% spread, which is awkward.
    if mn==mx:
        mn=min(np.min(data[0,:]),np.min(data[-1,:]))
        mx=max(np.max(data[0,:]),np.max(data[-1,:]))
        print('  In plot_3d, 1% and 99% vals were identical for '+var+'. Trying min and max.')

    #if min and max are truly identical, just give up.

    if mn==mx:
        print('  Min and max for var=%s were identical, equal to %f. Skipping.'%(var,mn))
    else: 

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

    outfi.savefig(pl.gcf().number)
    pl.close()

    return

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def plot_3d(ti,data,var,units,lats,lons,areas,outfi):

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

    #if 1% and 99% vals are identical, try min to max.
    #note this means some plots may show min-to-max and others 1% to 99% spread, which is awkward.
    if mn==mx:
        mn=min(np.min(data[0,:,-1]),np.min(data[-1,:,-1]))
        mx=max(np.max(data[0,:,-1]),np.max(data[-1,:,-1]))
        print('  In plot_3d, 1% and 99% vals were identical for '+var+'. Trying min and max.')

    #if min and max are truly identical, just give up.
    if mn==mx:
        print('  In plot_3d, min and max are identical for surf-lev var='+var+'. Skipping')

    else:
        colormap='viridis'
        if np.sign(mx)+np.sign(mn)<1:
            colormap='RdBu_r'
            max_mag=max(-1*mn,mx)
            mn=-1*max_mag
            mx=max_mag

        pl.subplot(3,2,3)
        pl.tricontourf(lons,lats,data[0,:,-1],levels=np.linspace(mn,mx,11),
                       cmap=colormap,extend='both')
        pl.colorbar()
        #pl.xlabel('lon')
        pl.ylabel('lat')
        pl.title('lowest layer, init time')

        pl.subplot(3,2,5)
        pl.tricontourf(lons,lats,data[-1,:,-1],levels=np.linspace(mn,mx,11),
                       cmap=colormap,extend='both')
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
        print('  In plot_3d, min and max are identical for maxlev var='+var+'. Skipping')

    else:
        colormap='viridis'
        if np.sign(mx)+np.sign(mn)<1:
            colormap='RdBu_r'
            max_mag=max(-1*mn,mx)
            mn=-1*max_mag
            mx=max_mag

        pl.subplot(3,2,4)
        pl.tricontourf(lons,lats,data[0,:,vert_ind],levels=np.linspace(mn,mx,11),
                       cmap=colormap,extend='both')
        pl.colorbar()
        #pl.xlabel('lon')
        #pl.ylabel('lat')
        pl.title('Layer %i, init time'%(vert_ind))

        pl.subplot(3,2,6)
        pl.tricontourf(lons,lats,data[-1,:,vert_ind],levels=np.linspace(mn,mx,11),
                       cmap=colormap,extend='both')
        pl.colorbar()
        pl.xlabel('lon')
        #pl.ylabel('lat')
        pl.title('Layer %i, final time'%(vert_ind))

    outfi.savefig(pl.gcf().number)
    pl.close()

    return
    #=-=-=-=-=-=-=-=-=-=-=-=
