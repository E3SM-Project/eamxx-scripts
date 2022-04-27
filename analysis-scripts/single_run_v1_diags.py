#!/usr/bin/env python
"""
Explore first v1 output
"""

from netCDF4 import MFDataset, Dataset
import numpy as np
import pylab as pl
import diags_utils

#USER CHOICES:
#=========================
#File to dump plots into:
output_pdf_name = 'output/ne30_test.pdf'

#Directory to search for model outputs:
#fi_dir='/global/cfs/cdirs/e3sm/terai/SCREAM/v1_analysis/run4/'
#fi_dir='/p/lustre2/caldwep/e3sm_scratch/quartz/SMS_D_Ld1_P180x1.ne30_ne30.F2000SCREAMv1.quartz_intel.20220425_165547_zaftpl/run/'
fi_dir='/p/lustre2/caldwep/e3sm_scratch/syrah/SMS_Ld1_P180x1.ne30_ne30.F2000SCREAMv1.syrah_intel.20220426_174217_xzm7bp/run/'

fi_name='*INSTANT*0.nc'

#Mapping file:
#map_dir='/global/homes/z/zender/data/maps/' #on NERSC
map_dir='/usr/gdata/climdat/maps/' #on LC

map_fi='map_ne30np4_to_fv128x256_aave.20160301.nc'
#map_fi='map_ne4np4_to_fv25x48_aave.20170401.nc'

#Create a list of variables to ignore
#((makes calcs faster and allows skipping big vars that would blow memory)
ignore=['time','T_prev_micro_step','qv_prev_micro_step','micro_liq_ice_exchange','micro_vap_liq_exchange','micro_vap_ice_exchange','eff_radius_qc','eff_radius_qi','inv_qc_relvar','SW_flux_dn_dir','p_int','nc_activated','ni_activated','nc_nuceat_tend','sfc_alb_dir_nir','sfc_alb_dif_nir','aero_g_sw','aero_ssa_sw','aero_tau_lw']

#INIT PDF TO SAVE ALL FIGURES TO:
#=========================
import matplotlib.backends.backend_pdf
outfi = matplotlib.backends.backend_pdf.PdfPages(output_pdf_name)

#GET MAPPING TOOLS:
#=========================
map_f=Dataset(map_dir+'/'+map_fi)

#area-weighted ave = np.sum(var[0:ncol]*areas)/np.sum(areas)
areas=map_f['area_a'][:]

#lat and lon at cell centers
lats=map_f['yc_a'][:]
lons=map_f['xc_a'][:]

#lats_edges=map_f['yv_a'][:] #dims=[len(data),5] since corners can have 5 vertices?
#lons_edges=map_f['xv_a'][:]

map_f.close()

#GET ACTUAL DATA:
#==========================
#note trickery below: restart files are *INSTANT*.r.nc and screw up MFDataset...
f=MFDataset(fi_dir+'/'+fi_name)

#time is fractional days since start of simulation
ti=f['time'][0:-1]

#edge pressures needed for vert integrals
p_int=f.variables['p_int'][0:-1]

vars=f.variables.keys()

for var in vars:

  #This chunk of code just skips to next var if var is to be ignored
  if var in ignore:
      print('Skipping '+var)
      continue
  
  print('Handling '+var)  
  data=f[var][0:-1]
  units=f[var].units

  if var=="surf_latent_flux":
    data=data*2.5e6; #convert kg w/(m2*s) to W/m2
    units='W/m2'

  if var=="surf_mom_flux":
    var="surf_mom_flux_u"
    diags_utils.plot_2d(ti,data[:,:,0],var,units,lats,lons,areas,outfi)
    var="surf_mom_flux_v"
    diags_utils.plot_2d(ti,data[:,:,1],var,units,lats,lons,areas,outfi)

  elif var=="horiz_winds":
    var="U"
    diags_utils.plot_3d(ti,data[:,:,0,:],var,units,lats,lons,areas,outfi)
    var="V"
    diags_utils.plot_3d(ti,data[:,:,1,:],var,units,lats,lons,areas,outfi)

  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  #FOR 2D VARIABLES:
  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  elif data.shape==(len(ti),len(lats)):
    diags_utils.plot_2d(ti,data,var,units,lats,lons,areas,outfi)

  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  #FOR 3D VARIABLES:
  #=-=-=-=-=-=-=-=-=-=-=-=-=-=
  #handles both lev and ilev cases unlike: (len(ti),len(lats),len(levs)):
  elif len(data.shape)==3: 
    diags_utils.plot_3d(ti,data,var,units,lats,lons,areas,outfi)

  else:
    print(var+' has weird size. skipping')


#LOOP OVER DERIVED VARIABLES:
#=================
var_dict={'LWP':'qc','IWP':'qi','WVP':'qv'}
for var in ['LWP','IWP','WVP']:
  print('Handling '+var)  
  data=f.variables[var_dict[var]][0:-1]
  data_vint=diags_utils.water_path(data,p_int)
  diags_utils.plot_2d(ti,data_vint,var,'g/m2',lats,lons,areas,outfi)
  
var='TBOT'
print('Handling '+var)  
data=f.variables['T_mid'][0:-1,:,-1]
units=f.variables['T_mid'].units
diags_utils.plot_2d(ti,data,var,units,lats,lons,areas,outfi)

var='QBOT'
print('Handling '+var)  
data=f.variables['qv'][0:-1,:,-1]
units=f.variables['qv'].units
diags_utils.plot_2d(ti,data,var,units,lats,lons,areas,outfi)


#CLEAN UP/FINALIZE
#=================
outfi.close()
#pl.show()
print('*** DONE ***')

