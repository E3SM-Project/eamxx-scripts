#!/bin/csh 
## ====================================================================
# Purpose is to do post-processing on E3SM/SCREAM output
#  to generate prescribed aerosol climatology.  

##### USER INPUT BEGIN

# Path to where output is located (most likely *.cam.h0* files) 
setenv datapath /global/cscratch1/sd/terai/E3SM_simulation/maint-1.0.F2010C5-CMIP6-LR.E3SMv1_PresAeroFinal_201102.ne30_oEC/run
# The case name
setenv casename maint-1.0.F2010C5-CMIP6-LR.E3SMv1_PresAeroFinal_201102.ne30_oEC

# Desired prefix name of aerosol generated output file
setenv outname mam4_0.9x1.2_L72_2000FSCREAMLR_c201110

# Location of mapfile to convert from ne30 to 1-degree FV
setenv mapfile /global/cfs/cdirs/e3sm/mapping/maps/map_ne30np4_to_0.9x1.25_aave_110121.nc

# Year start
setenv years 2011 
# Year end
setenv yeare 2015

##### USER INPUT END

set the_year = $years

# First pull out only relevant information pertaining to the prescribed aerosol files
echo "Making subset of data"
while ($the_year <= $yeare)
  set the_month = 1
  
  while ($the_month <= 12)
    set the_month_s = `printf "%02d" $the_month`
    set the_year_s = `printf "%04d" $the_year`
      
    # Generate subset files
    setenv thefile $datapath/{$casename}.cam.h0.{$the_year_s}-{$the_month_s}.nc
    setenv outfile $datapath/subset_{$the_year_s}_{$the_month_s}.nc

    ncks -O -v P0,PS,bc_a1_logm,bc_a1_logv,dst_a1_logm,dst_a1_logv,dst_a3_logm,dst_a3_logv,ncl_a1_logm,ncl_a1_logv,ncl_a2_logm,ncl_a2_logv,ncl_a3_logm,ncl_a3_logv,num_a1_logm,num_a1_logv,num_a2_logm,num_a2_logv,num_a3_logm,num_a3_logv,pom_a1_logm,pom_a1_logv,so4_a1_logm,so4_a1_logv,so4_a2_logm,so4_a2_logv,so4_a3_logm,so4_a3_logv,soa_a1_logm,soa_a1_logv,soa_a2_logm,soa_a2_logv,bc_a1,bc_a1DDF,bc_a1SFWET,bc_c1,bc_c1DDF,bc_c1SFWET,dst_a1,dst_a1DDF,dst_a1SFWET,dst_a3,dst_a3DDF,dst_a3SFWET,dst_c1,dst_c1DDF,dst_c1SFWET,dst_c3,dst_c3DDF,dst_c3SFWET,date,hyai,hyam,hybi,hybm,ilev,lev,lon,ncl_a1,ncl_a1DDF,ncl_a1SFWET,ncl_a2,ncl_a2DDF,ncl_a2SFWET,ncl_a3,ncl_a3DDF,ncl_a3SFWET,ncl_c1,ncl_c1DDF,ncl_c1SFWET,ncl_c2,ncl_c2DDF,ncl_c2SFWET,ncl_c3,ncl_c3DDF,ncl_c3SFWET,num_a1,num_a1DDF,num_a1SFWET,num_a2,num_a2DDF,num_a2SFWET,num_a3,num_a3DDF,num_a3SFWET,num_c1,num_c1DDF,num_c1SFWET,num_c2,num_c2DDF,num_c2SFWET,num_c3,num_c3DDF,num_c3SFWET,pom_a1,pom_a1DDF,pom_a1SFWET,pom_c1,pom_c1DDF,pom_c1SFWET,so4_a1,so4_a1DDF,so4_a1SFWET,so4_a2,so4_a2DDF,so4_a2SFWET,so4_a3,so4_a3DDF,so4_a3SFWET,so4_c1,so4_c1DDF,so4_c1SFWET,so4_c2,so4_c2DDF,so4_c2SFWET,so4_c3,so4_c3DDF,so4_c3SFWET,soa_a1,soa_a1DDF,soa_a1SFWET,soa_a2,soa_a2DDF,soa_a2SFWET,soa_c1,soa_c1DDF,soa_c1SFWET,soa_c2,soa_c2DDF,soa_c2SFWET,time,time_bnds,wat_a1,wat_a2,wat_a3 $thefile $outfile

    set the_month = `expr $the_month + 1`
  end
  set the_year = `expr $the_year + 1`
end

# Create climatological monthly averages
echo "Generating climatological monthly averages"
set the_month = 1
while ($the_month <= 12)

    set the_month_s = `printf "%02d" $the_month`
    
    ncra -O {$datapath}/subset*{$the_month_s}.nc {$datapath}/climo_${the_month_s}.nc 
    
    set the_month = `expr $the_month + 1`
    
end

rm subset*.nc

cd {$datapath}

echo "Changing time and time bounds in files"
# time needs to be changed for each file, to appease the prescribe aerosol code
ncap2 -O -s 'time(0)=31' -s 'date(0)=10201' -s 'time_bnds(0,0) = 0' -s 'time_bnds(0,1) = 31' climo_01.nc climo_01_out.nc
ncap2 -O -s 'time(0)=59' -s 'date(0)=10301' -s 'time_bnds(0,0) = 31' -s 'time_bnds(0,1) = 59' climo_02.nc climo_02_out.nc
ncap2 -O -s 'time(0)=90' -s 'date(0)=10401' -s 'time_bnds(0,0) = 59' -s 'time_bnds(0,1) = 90' climo_03.nc climo_03_out.nc
ncap2 -O -s 'time(0)=120' -s 'date(0)=10501' -s 'time_bnds(0,0) = 90' -s 'time_bnds(0,1) = 120' climo_04.nc climo_04_out.nc
ncap2 -O -s 'time(0)=151' -s 'date(0)=10601' -s 'time_bnds(0,0) = 120' -s 'time_bnds(0,1) = 151' climo_05.nc climo_05_out.nc
ncap2 -O -s 'time(0)=181' -s 'date(0)=10701' -s 'time_bnds(0,0) = 151' -s 'time_bnds(0,1) = 181' climo_06.nc climo_06_out.nc
ncap2 -O -s 'time(0)=212' -s 'date(0)=10801' -s 'time_bnds(0,0) = 181' -s 'time_bnds(0,1) = 212' climo_07.nc climo_07_out.nc
ncap2 -O -s 'time(0)=243' -s 'date(0)=10901' -s 'time_bnds(0,0) = 212' -s 'time_bnds(0,1) = 243' climo_08.nc climo_08_out.nc
ncap2 -O -s 'time(0)=273' -s 'date(0)=11001' -s 'time_bnds(0,0) = 243' -s 'time_bnds(0,1) = 273' climo_09.nc climo_09_out.nc
ncap2 -O -s 'time(0)=304' -s 'date(0)=11101' -s 'time_bnds(0,0) = 273' -s 'time_bnds(0,1) = 304' climo_10.nc climo_10_out.nc
ncap2 -O -s 'time(0)=334' -s 'date(0)=11201' -s 'time_bnds(0,0) = 304' -s 'time_bnds(0,1) = 334' climo_11.nc climo_11_out.nc
ncap2 -O -s 'time(0)=365' -s 'date(0)=20101' -s 'time_bnds(0,0) = 334' -s 'time_bnds(0,1) = 365' climo_12.nc climo_12_out.nc

# Now stitch the monthly files together into one file
echo "Stitching files together"
ncrcat -O climo_01_out.nc climo_02_out.nc climo_03_out.nc climo_04_out.nc climo_05_out.nc climo_06_out.nc climo_07_out.nc climo_08_out.nc climo_09_out.nc climo_10_out.nc climo_11_out.nc climo_12_out.nc {$outname}se.nc

# Clean up some stuff
echo "Cleaning up"
rm *_out.nc
rm climo_*.nc
rm subset*.nc

# convert to FV grid, which is needed for inputs to E3SM
echo "Converting from SE to FV"
ncks --map=$mapfile {$outname}se.nc {$outname}.nc

rm {$outname}se.nc

echo "Program done"
