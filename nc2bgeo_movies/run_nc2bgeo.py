#!/usr/bin/env python
"""
Load E3SM netcdf output and use nc2bgeo module
to convert it to .bgeo format for loading into
Houdini

Notes:
1. 256x512x128 regridded 3d data takes ~3 min to write to file
2. Subsetting the regridded data to 0.2% of its original volume 
   takes 4 min or so because checking whether each point is in range
   is slow (while writing less data is faster).
3. Processing native grid data is hard because you can't store many
   copies of the global 3d array in memory at once but looping over the 
   25M columns takes forever. Chopping off all levels above 100 mb (~30%
   of the levels) and not running interactively to aid trash collection(?)
   seems to work on cori-haswell.
4. If latslice is non-empty than lonslice needs to be as well.

"""

#IMPORT STUFF:
#==============
from netCDF4 import Dataset
import numpy as np
#import pylab as pl
import time
import nc2bgeo
"""
#INPUTS FOR A GLOBAL 3D 256X521 REGRIDDED VARIABLE
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=--=-=-==-=-=-

varname='CLDLIQ'
regridded=False
vertdim=True #if 3d data rather than just lat/lon
####outfile=varname+'_14-25Sx98-110W.bgeo'
#outfile=varname+'_native_20200216_timesnap1.bgeo'
latlonht=True #if true, dims=lat,lon,ht. Otherwise meters from center of earth

#var_fi='/global/cscratch1/sd/terai/SCREAM/ForBradCarvey/'\
#        'rgr_1800x3600_SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.2020-02-16-00000.nc'

#Only need the next file if regridded=False.
#native_grid_fi='/global/homes/z/zender/data/grids/ne1024pg2.nc'

#Only need the next 3 files if vertdim=True
#- - - - - - - - - - - - - - - - - - - - - - - -
#topo_fi='/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/'\
#           'USGS-gtopo30_ne1024np4pg2_16xconsistentSGH_20190528_converted.nc'
#topo_fi='/global/cscratch1/sd/terai/SCREAM/ForBradCarvey/'\
#            'rgr_1800x3600_USGS-gtopo30_ne1024np4pg2_16xconsistentSGH_20190528_converted.nc' #for regridded
#note: if regridded=True, no topo file exists yet. Create it via
#      ncremap -v PHIS -m map.nc ${topo_fi}.nc ${topo_fi}_regridded256x512.nc
#      where map.nc is extracted from ncdump -h ${var_fi}
#ps_fi='/global/cscratch1/sd/terai/SCREAM/ForBradCarvey/'\
#       'rgr_1800x3600_PS_SCREAMv0.eam.h4.2020-02-16-00000.nc'
#T_fi='/global/cscratch1/sd/terai/SCREAM/ForBradCarvey/'\
#      'rgr_1800x3600_T_SCREAMv0.eam.h6.2020-02-16-00000.nc'

timesnap=1 #index of time entry to write to .bgeo. Must be a scalar integer.

###latslice=[-25,-14]
###lonslice=[98,110]
#latslice=0
#lonslice=0
"""
toplev=35  #skip levels closer to model top than this. 35=100mb


#INPUTS FOR GLOBAL 3D NATIVE-GRID VARIABLE
#=============================================
varname='CLDLIQ'   #'CLDLIQ' on h7
regridded=False
vertdim=True #if 3d data rather than just lat/lon
outfile=varname+'_native_14-25Sx98-110W.bgeo'

timesnap=0 #index of time entry to write to .bgeo. Must be a scalar integer.
latslice=[-25,-14]
lonslice=[98,110]
#latslice=0
#lonslice=0
latlonht=True #if true, dims=lat,lon,ht. Otherwise meters from center of earth

time_on_file='2020-02-16-00000' 

var_fi='/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127/'\
        +'SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h7.'+time_on_file+'.nc'

#Next 3 files aren't used if vertdim=False; You shouldn't need to edit these file names.
#- - - - - - - - - - - - - - - - - - - - - - - -
topo_fi='/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/'\
           'USGS-gtopo30_ne1024np4pg2_16xconsistentSGH_20190528_converted.nc'

ps_fi='/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127/'\
        +'SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h4.'+time_on_file+'.nc'
T_fi='/global/cfs/cdirs/e3sm/terai/SCREAM/DYAMOND2/Output/20201127/'\
        +'SCREAMv0.SCREAM-DY2.ne1024pg2.20201127.eam.h6.'+time_on_file+'.nc'

#Only need the next file if regridded=False.
#native_grid_fi='/global/homes/z/zender/data/grids/ne1024pg2.nc' #is netCDF-4, doesn't work on compute nodes
native_grid_fi='/global/homes/p/petercal/py/e3sm/movies/ne1024pg2-cdf5.nc' #converted to cdf5 using nccopy -k; seems to work

#DEFINE CONSTANTS:
#============================================
Rd=287.15 #gas constant for dry air
grav=9.8  #gravity

#LOAD OUTPUT
#============================================
start=time.time()

f=Dataset(var_fi)

#Note data needs to be in numpy format and is expected 32 rather than 64 bit.
#That's what [] and .astype commands ensure.
if vertdim:
    var=f.variables[varname][timesnap,toplev:].astype(np.float32)
else:
    var=f.variables[varname][timesnap].astype(np.float32)

#Eventually will want to get time info automatically from file...
#ti=f.variables['time'][:].astype(np.float32)

if regridded:
    lat=f.variables['lat'][:].astype(np.float32)
    lon=f.variables['lon'][:].astype(np.float32)

f.close()

if not regridded:
    f=Dataset(native_grid_fi)
    lat=f.variables['grid_center_lat'][:].astype(np.float32)
    lon=f.variables['grid_center_lon'][:].astype(np.float32)
    #also has grid_corner_lat and grid_corner_lon if we
    #eventually want cell connectivity!
    f.close()

end=time.time()
print('Getting var took %f secs'%(end-start))

#COMPUTE HEIGHT DIMENSION FOR 3D OUTPUT:
#=======================================
if not vertdim:
    hts=0 #put 2d vars at sea level by default. May want to change this to top of atm for some vars.

else:
    
    #GET SURFACE HEIGHT
    #------------------
    start=time.time()

    f=Dataset(topo_fi)
    surf_z=f.variables['PHIS'][:].astype(np.float32)/grav #PHIS had units m2/s2
    f.close()
    end=time.time()
    print('Getting surf_z took %f sec'%(end-start))
    
    #GET SPATIALLY-VARYING PRESSURE FIELD:
    #------------------
    start=time.time()
        
    f=Dataset(ps_fi)
    ps=f.variables['PS'][timesnap].astype(np.float32)
    f.close()


    #COMBINE SURFACE HEIGHT AND PRESSURE W/ HYPSOMETRIC EQ TO GET HEIGHT
    #------------------   
    f=Dataset(T_fi)
    #note hack: should be *virtual* temperature here!
    T=f.variables['T'][timesnap,toplev:].astype(np.float32)
    hyam=f.variables['hyam'][toplev:].astype(np.float32)
    hybm=f.variables['hybm'][toplev:].astype(np.float32)
    f.close()

    P=nc2bgeo.make_3d_pres(ps,hyam,hybm)
    
    hts=np.zeros(var.shape,np.float32)
    hts[-1]=surf_z+10. #this is bottom edge of lowest cell. Adding 10m as a guess of center.
    #Can't avoid looping here since lev depends on the one above it.
    for i in range(len(hyam)-2,-1,-1):
        hts[i]= hts[i+1] + Rd*T[i]/grav*np.log(P[i+1]/P[i])
        
    end=time.time()
    print('Getting hts took %f sec'%(end-start))

#MAKE LAT AND LON THE SAME SIZE AS VAR FOR VECTOR OPERATIONS
#==============================
start=time.time()
if vertdim:
    if regridded:
        LAT,LON=nc2bgeo.resize_latlonlev(lat,lon,hts)
    else:
        LAT,LON=nc2bgeo.resize_ncollev(lat,lon,hts)
else: #if 2d
    if regridded:
        LAT,LON=nc2bgeo.resize_latlon(lat,lon)
    else:
        LAT=lat
        LON=lon

end=time.time()
print('Vectorizing lat/lon took %f sec'%(end-start))

#GET POINTS TO WRITE
#=======================
start=time.time()
if latlonht==False: #if output is x,y,z in meters with origin center of earth

    if latslice!=0 or lonslice!=0:
        #if need to subset array, need to iterate over indices
        X,Y,Z,V = nc2bgeo.subset_latlonht2xyz(LAT,LON,hts,var,latslice,lonslice)
    else:
        #if using all elements, can do a simple/fast array operation
        X,Y,Z = nc2bgeo.latlonht2xyz(LAT,LON,hts)
        V=var

else:
    if latslice==0 and lonslice==0:
        #recall that this function insists that lat and lon already have size of var
        X=LAT
        Y=LON
        V=var
        if np.isscalar(hts):
            Z=lev*np.ones(X.shape)
        else:
            Z=hts
    else:
        X,Y,Z,V=nc2bgeo.subset_latlonht(LAT,LON,hts,var,latslice,lonslice)
        
end=time.time()
print('Finding points to write took %f sec'%(end-start))
        
#DEBUGGING
#==============================
"""
#does hts look right?
pl.figure()
if regridded: #regridded peak ~5km=16k ft
    pl.pcolor(LON[0],LAT[0],surf_z);
    pl.colorbar();
    pl.title('Surf Ht (m)')
else:
    pl.plot(surf_z);
    pl.title('Surf Ht (m)')
pl.figure()
if regridded: #1st lev at 40km, last at surf.
    pl.plot(hts[:,0,0]);
    pl.title('heights')
else:
    pl.plot(hts[:,0]);
    pl.title('heights')
pl.show()
"""

#WRITE THE FILE
#==============================
start=time.time()
nc2bgeo.writeBgeoFile('/global/cscratch1/sd/terai/SCREAM/ForBradCarvey/'+outfile,varname,V,X,Y,Z)
end=time.time()
print('Writing file took %f sec'%(end-start))
