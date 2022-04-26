#!/usr/bin/env python
"""
Explore first v1 output
"""

from netCDF4 import MFDataset, Dataset
import numpy as np
import pylab as pl
import xarray

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def water_path(q,p_int):
    """
    Compute the liquid, ice, or vapor water path given a 3d variable q (which can 
    be qc, qi, or qv, respectively). Output units are g/m2.
    """
    g=9.8 #gravity
    if p_int.shape[2]==73:
      dp=p_int[:,:,0:-1] - p_int[:,:,1:] #Note 0 is TOA and -1 is surf so dp should be neg.
      wp=-1000.*np.sum(np.array(q)*np.array(dp)/g,axis=-1)
    else:
      dp=p_int[:,0:-1,:] - p_int[:,1:,:]
      wp=-1000.*np.sum(np.array(q)*np.array(dp)/g,axis=1)

    return wp

def vertintegrated(q,p_int):
    """
    Compute the vertical integral given a 3d variable q (which can 
    be nc, ni, or nr, respectively). Output units are 1/m2.
    """
    g=9.8 #gravity
    if p_int.shape[2]==73:
      dp=p_int[:,:,0:-1] - p_int[:,:,1:] #Note 0 is TOA and -1 is surf so dp should be neg.
      wp=-1.*np.sum(np.array(q)*np.array(dp)/g,axis=-1)
    else:
      dp=p_int[:,0:-1,:] - p_int[:,1:,:]
      wp=-1.*np.sum(np.array(q)*np.array(dp)/g,axis=1)

    return wp

def toplevel(var):
    """
    Compute the top level outputs from a variable
    """
    if (var.shape[2]==72) | (var.shape[2]==73):
      wp=np.array(var[:,:,0])
    else:
      wp=np.array(var[:,0,:])
    return wp

def bottomlevel(var):
    """
    Compute the bottom level outputs from a variable
    """
    if (var.shape[2]==72) | (var.shape[2]==73):
      wp=np.array(var[:,:,-1])
    else:
      wp=np.array(var[:,-1,:])
    return wp

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#SAVE ALL FIGS AS 1 PDF:
#=========================
import matplotlib.backends.backend_pdf
#pdf = matplotlib.backends.backend_pdf.PdfPages("/global/cfs/cdirs/e3sm/www/terai/SCREAM/v1_analysis/ne30_v1_v0_nstep1_v0nocloudIC_20220411.pdf")
pdf = matplotlib.backends.backend_pdf.PdfPages("/global/cfs/cdirs/e3sm/www/terai/SCREAM/v1_analysis/ne30_v1_v0_step2_comparison_220421.pdf")

#GET MAPPING TOOLS:
#=========================

# v1
map_fi='/global/cfs/cdirs/e3sm/terai/mapping/map_ne30np4_to_cmip6_180x360_nco.20210525.nc'
map_f=Dataset(map_fi)

#area-weighted ave = np.sum(var[0:ncol]*areas)/np.sum(areas)
areas_v1=map_f['area_a'][:]

#lat and lon at cell centers
lats_v1=map_f['yc_a'][:]
lons_v1=map_f['xc_a'][:]

lats_edges_v1=map_f['yv_a'][:] #dims=[len(data),5] since corners can have 5 vertices?
lons_edges_v1=map_f['xv_a'][:]

map_f.close()

# v0
map_fi='/global/cfs/cdirs/e3sm/terai/mapping/map_ne30np4_to_cmip6_180x360_nco.20210525.nc' #'/global/cfs/cdirs/e3sm/terai/mapping/map_ne30pg2_to_cmip6_180x360_nco.20200901.nc'
map_f=Dataset(map_fi)

#area-weighted ave = np.sum(var[0:ncol]*areas)/np.sum(areas)
areas_v0=map_f['area_a'][:]

#lat and lon at cell centers
lats_v0=map_f['yc_a'][:]
lons_v0=map_f['xc_a'][:]

lats_edges_v0=map_f['yv_a'][:] #dims=[len(data),5] since corners can have 5 vertices?
lons_edges_v0=map_f['xv_a'][:]

map_f.close()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def plot_2d(lats_a,lons_a,data_a,lats_b,lons_b,data_b,units):
  
  pl.figure(figsize=[8.5,11])

  #GLOBAL stats_pdf:
  #=========================
  # Use 1st, 10th, 90th, and 99th percentiles to determine PDF and map colorscale edges
  data_a_prctile=np.percentile(data_a,[1,10,90,99])
  data_b_prctile=np.percentile(data_b,[1,10,90,99])

  ratio_small_vals_a=len(data_a[(data_a<1e-7) & (data_a>-1e-7)])/len(data_a)
  ratio_small_vals_b=len(data_b[(data_b<1e-7) & (data_b>-1e-7)])/len(data_b)
  print(ratio_small_vals_a)
  print(ratio_small_vals_b)
  if (ratio_small_vals_a>0.1) & (ratio_small_vals_b>0.1):
      data_a_wozero=data_a.copy()
      data_b_wozero=data_b.copy()
      data_a_wozero[np.absolute(data_a_wozero)<1e-9]=np.nan
      data_b_wozero[np.absolute(data_b_wozero)<1e-9]=np.nan
      data_a_wozero_prctile=np.nanpercentile(data_a_wozero,[0.1,10,90,99.9])
      data_b_wozero_prctile=np.nanpercentile(data_b_wozero,[0.1,10,90,99.9])
      bin_edges=np.min([data_a_wozero_prctile[0],data_b_wozero_prctile[0]])*(np.max([data_a_wozero_prctile[3],data_b_wozero_prctile[3]])/np.min([data_a_wozero_prctile[0],data_b_wozero_prctile[0]]))**(np.arange(51)/50)
      hist_a,bin_edges=np.histogram(data_a,bins=bin_edges)
      hist_b,bin_edges=np.histogram(data_b,bins=bin_edges)
      bin_means=(bin_edges[1:]*bin_edges[:-1])**0.5
      xscale='log'
  else:
      hist_a,bin_edges=np.histogram(data_a,bins=50,range=(np.min([data_a_prctile[0],data_b_prctile[0]]),np.max([data_a_prctile[3],data_b_prctile[3]])))
      hist_b,bin_edges=np.histogram(data_b,bins=50,range=(np.min([data_a_prctile[0],data_b_prctile[0]]),np.max([data_a_prctile[3],data_b_prctile[3]])))
      bin_means=(bin_edges[1:]+bin_edges[:-1])/2.
      xscale='linear'
  

  pl.subplot(4,1,1)
  pl.plot(bin_means,hist_a/len(data_a),label='v1')
  pl.plot(bin_means,hist_b/len(data_b),label='v0')
  pl.xscale(xscale)
  pl.legend(loc='upper right')
  pl.title(var+' Global PDF, units=%s)'%(units))
  pl.xlabel(var)
  pl.ylabel('Probability')
  
  #MAPS:
  #=========================
  #use delaunay triangulation to plot on native grid

  #use same color axis for both init and final times:
  mn=np.min([data_a_prctile[1],data_b_prctile[1]])
  mx=np.max([data_a_prctile[2],data_b_prctile[2]])

  mn=min(data_a_prctile[0],data_b_prctile[0])
  mx=max(data_a_prctile[3],data_b_prctile[3])
  colormap='viridis'
  if np.sign(mx)+np.sign(mn)<1:
      colormap='RdBu_r'
      max_mag=max(-1*mn,mx)
      mn=-1*max_mag
      mx=max_mag
  
  pl.subplot(3,1,2)
  pl.tricontourf(lons_a,lats_a,data_a,levels=
                 np.linspace(mn,mx,11),cmap=colormap,extend='both')
  pl.colorbar()
  #pl.xlabel('lon')
  pl.ylabel('lat')
  pl.title('v1 output')
  
  pl.subplot(3,1,3)
  pl.tricontourf(lons_b,lats_b,data_b,levels=
                 np.linspace(mn,mx,11),cmap=colormap,extend='both')
  pl.colorbar()
  pl.xlabel('lon')
  pl.ylabel('lat')
  pl.title('v0 output')

  pl.show()
  
  pdf.savefig(pl.gcf().number)

  return
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def plot_3d(data_a,data_b,units):

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
    mn=min(np.min(data[0,:,-1]),np.min(data[-1,:,-1]))
    mx=max(np.max(data[0,:,-1]),np.max(data[-1,:,-1]))

    mn=min(data_a_prctile[0],data_b_prctile[0])
    mx=max(data_a_prctile[3],data_b_prctile[3])
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

    #use same color axis for both init and final times:
    mn=min(np.min(data[0,:,vert_ind]),np.min(data[-1,:,vert_ind]))
    mx=max(np.max(data[0,:,vert_ind]),np.max(data[-1,:,vert_ind]))

    mn=min(data_a_prctile[0],data_b_prctile[0])
    mx=max(data_a_prctile[3],data_b_prctile[3])
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
        pl.tricontourf(lons,lats,data[0,:,vert_ind],levels=np.linspace(mn,mx,11) )
        pl.colorbar()
        #pl.xlabel('lon')
        #pl.ylabel('lat')
        pl.title('Layer %i, init time'%(vert_ind))

        pl.subplot(3,2,6)
        pl.tricontourf(lons,lats,data[-1,:,vert_ind],levels=np.linspace(mn,mx,11))
        pl.colorbar()
        pl.xlabel('lon')
        #pl.ylabel('lat')
        pl.title('Layer %i, final time'%(vert_ind))


    pdf.savefig(pl.gcf().number)

    return
#=-=-=-=-=-=-=-=-=-=-=-=

#GET ACTUAL DATA:
#==========================

#fi_dir='/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/SMS_D_Ln2_P675x1.ne30_ne30.F2010SCREAMv1.cori-knl_intel.20220412_140616_4io5m1/run/'
#f_v1=xarray.open_mfdataset(fi_dir+'SCREAMv1_output_nonhydrostatic.INSTANT.Steps_x1.np675.0001-01-01.001000.nc', drop_variables=('P3_input_dim', 'P3_output_dim'))

fi_dir='/global/cfs/cdirs/e3sm/terai/SCREAM/v1_analysis/run3/'
f_v1=xarray.open_mfdataset(fi_dir+'SCREAMv1_output_nonhydrostatic.INSTANT.Steps_x1.np12.2010-01-01.001000.nc', drop_variables=('P3_input_dim', 'P3_output_dim'))


fi_dir='/global/cscratch1/sd/terai/e3sm_scratch/cori-knl/master.ne30_ne30.F2010-SCREAM-LR.cori-knl_intel.42x33x4.v1_nocloudnospa_novertdiff.20220415-1324/run/'
f_v0=xarray.open_mfdataset(fi_dir+'master.ne30_ne30.F2010-SCREAM-LR.cori-knl_intel.42x33x4.v1_nocloudnospa_novertdiff.20220415-1324.eam.h1.2010-01-01-00600.nc', drop_variables=('P3_input_dim', 'P3_output_dim'))


#time is fractional days since start of simulation
ti=f_v1['time'][:] #days since start of simulation
#edge pressures
fi_dir='/global/cfs/cdirs/e3sm/terai/SCREAM/v1_analysis/run4/'
f_v1_sp=xarray.open_mfdataset(fi_dir+'SCREAMv1_output_nonhydrostatic.INSTANT.Steps_x1.np12.2010-01-01.001000.nc', drop_variables=('P3_input_dim', 'P3_output_dim'))

p_int=f_v1_sp.variables['p_int'][:]

PS_v0=f_v0['PS'][:]
hyai_v0=f_v0['hyai']
hybi_v0=f_v0['hybi']
P0_v0=f_v0['P0']
aiP=np.array(hyai_v0)*np.array(P0_v0)
aiP_reshape=aiP.reshape(1,73,1)
PS_reshape=np.array(PS_v0).reshape(PS_v0.shape[0],1,PS_v0.shape[1])
bi_reshape=np.tile(np.array(hybi_v0).reshape(1,73,1),(PS_reshape).shape)
Pi_v0=aiP_reshape+bi_reshape*PS_reshape
print('Pi_v0[0,0,0]')
print(Pi_v0[0,0,0])
print('p_int[0,0,0]')
print(np.array(p_int)[0,0,0])
vars=f_v1.variables.keys()

var_dict={'surf_latent_flux':'LHFLX','surf_sens_flux':'SHFLX','ps':'PS','phis':'PHIS','precip_liq_surf':'PRECT','qc':'CLDLIQ','qr':'RAINQM','qi':'CLDICE','qv':'Q','nc':'NUMLIQ','nr':'NUMRAI','ni':'NUMICE', \
          'LW_flux_up':'FUL','LW_flux_dn':'FDL','SW_flux_up':'FUS','SW_flux_dn':'FDS','T_mid':'T','horiz_winds':'U','omega':'OMEGA'} #'qv':'TMQ',

for var in ['surf_latent_flux','surf_sens_flux','precip_liq_surf','ps','qc','qr','qi','qv','T_mid','omega','horiz_winds']: #vars:   ,'nc','nr','ni'
  if var in var_dict.keys():
    print('Handling '+var)  
    data_v1=f_v1[var][:]
    units=f_v1[var].units

    if var=="surf_latent_flux":
      data_v1=data_v1*2.5e6; #convert kg w/(m2*s) to W/m2
      units='W/m2'
    if len(data_v1.shape)<4:
      data_v0=f_v0[var_dict[var]]

    if var=="precip_liq_surf":
        data_v1=data_v1+f_v1['precip_ice_surf'][:]
        data_v1=data_v1*3600.*24.*1000.
        data_v0=data_v0*3600.*24.*1000.
        units='mm/d'

    if var in ['qv','qc','qr','qi']:
        data_flattened=water_path(data_v1,p_int)
        data_v1=data_flattened
        print(var)
        print(data_v1.shape)
        data_flattened=water_path(data_v0,Pi_v0)
        data_v0=data_flattened
        print(data_v0.shape)
        units='g/m2'

    if var in ['nc','nr','ni']:
        data_flattened=vertintegrated(data_v1,p_int)
        data_v1=data_flattened
        print(var)
        print(data_v1.shape)
        data_flattened=vertintegrated(data_v0,Pi_v0)
        data_v0=data_flattened
        print(data_v0.shape)
        units='1/m2'

    if var in ['LW_flux_up','LW_flux_dn','SW_flux_up','SW_flux_dn','T_mid','omega']:
        data_flattened=toplevel(data_v1)
        data_v1_top=data_flattened
        print(var)
        print(data_v1.shape)
        print(data_v1_top.shape)
        data_flattened=bottomlevel(data_v1)
        data_v1_bottom=data_flattened
        #process v0 data
        data_flattened=toplevel(data_v0)
        data_v0_top=data_flattened
        print(var)
        data_flattened=bottomlevel(data_v0)
        data_v0_bottom=data_flattened
        print(data_v0_bottom.shape)

    if var in ['horiz_winds']:
        u_v0=f_v0['U']
        v_v0=f_v0['V']
        u_v1=data_v1[:,:,0,:]
        v_v1=data_v1[:,:,1,:]
        # for U
        u_flattened=toplevel(u_v1)
        u_v1_top=u_flattened
        print(var)
        print(u_v1.shape)
        print(u_v1_top.shape)
        u_flattened=bottomlevel(u_v1)
        u_v1_bottom=u_flattened
        #process v0 u
        u_flattened=toplevel(u_v0)
        u_v0_top=u_flattened
        print(var)
        u_flattened=bottomlevel(u_v0)
        u_v0_bottom=u_flattened
        print(u_v0_bottom.shape)
        # for V
        v_flattened=toplevel(v_v1)
        v_v1_top=v_flattened
        print(var)
        print(v_v1.shape)
        print(v_v1_top.shape)
        v_flattened=bottomlevel(v_v1)
        v_v1_bottom=v_flattened
        #process v0 v
        v_flattened=toplevel(v_v0)
        v_v0_top=v_flattened
        print(var)
        v_flattened=bottomlevel(v_v0)
        v_v0_bottom=v_flattened
        print(v_v0_bottom.shape)
    
    #=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #FOR 2D VARIABLES:
    #=-=-=-=-=-=-=-=-=-=-=-=-=-=
    if data_v1.shape==(len(ti),len(lats_v1)):
      plot_2d(lats_v1,lons_v1,np.array(data_v1)[0,:],lats_v0,lons_v0,np.array(data_v0)[0,:],units)
    else:
      if var in ['horiz_winds']:
        var='U'
        plot_2d(lats_v1,lons_v1,np.array(u_v1_top)[0,:],lats_v0,lons_v0,np.array(u_v0_top)[0,:],units+' top')
        plot_2d(lats_v1,lons_v1,np.array(u_v1_bottom)[0,:],lats_v0,lons_v0,np.array(u_v0_bottom)[0,:],units+' bot')
        var='V'
        plot_2d(lats_v1,lons_v1,np.array(v_v1_top)[0,:],lats_v0,lons_v0,np.array(v_v0_top)[0,:],units+' top')
        plot_2d(lats_v1,lons_v1,np.array(v_v1_bottom)[0,:],lats_v0,lons_v0,np.array(v_v0_bottom)[0,:],units+' bot')
      else:    
        print('data shape is 3d')
        plot_2d(lats_v1,lons_v1,np.array(data_v1_top)[0,:],lats_v0,lons_v0,np.array(data_v0_top)[0,:],units+' top')
        plot_2d(lats_v1,lons_v1,np.array(data_v1_bottom)[0,:],lats_v0,lons_v0,np.array(data_v0_bottom)[0,:],units+' bot')

    #=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #FOR 3D VARIABLES:
    #=-=-=-=-=-=-=-=-=-=-=-=-=-=
    #handles both lev and ilev cases unlike: (len(ti),len(lats),len(levs)):
    #if len(data_v1.shape)==3: 
    #  plot_3d(data_v1,units)

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
