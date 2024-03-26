import glob
import numpy as np
import os
import os.path
import re
from netCDF4 import Dataset
from cftime import num2date, date2num
from datetime import datetime, timedelta

# Input files
input_path = '/pscratch/sd/n/ndk/e3sm_scratch/pm-gpu/se66-may15-11145e5338/t.se66-may15-11145e5338.F2010-SCREAMv1.ne120pg2_ne120pg2.pm-gpu.n085t4xX.vth200/run'
input_file = 'output.scream.monthly.AVERAGE.nmonths_x1.????-??-??-00000.nc'

# Output files
output_path = '/pscratch/sd/t/terai/EAMxx/CessExp/ne120/pmcpu_cntl/run'
output_file = 't.se66-may15-11145e5338.F2010-SCREAMv1.ne120pg2_ne120pg2.pm-gpu.n085t4xX.vth200.eam.h0.YYYY-MM.nc'

# Create output directory
#os.makedirs(output_path)

# Lits of input files
input_files = sorted(glob.glob(os.path.join(input_path, input_file)))

# Loop over input files
for ifile in input_files:

    # Extract time stamp string
    fname = os.path.basename(ifile)
    tmp = re.findall(r"[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}[.]nc$", fname)
    if len(tmp) != 1:
       print("ERROR: Could not file time stamp in %s" % (fname))
       sys.exit("Input file name error")
    tstamp = tmp[0][:-3]
    YYYY = tstamp[0:4]
    MM = tstamp[5:7]
    DD = tstamp[8:10]
    SSSSS = tstamp[11:16]

    # Output file name
    ofile = output_file.replace("YYYY", YYYY).replace("MM", MM).replace("DD", DD).replace("SSSSS", SSSSS)
    print(ofile)

    # Open model input file
    model = Dataset(os.path.join(input_path, ifile), "r")

    # Open output file
    output = Dataset(os.path.join(output_path, ofile), "w")

    # Global attribute
    output.Conventions = "CF-1.8"

    # Copy dimensions
    dimensions = [model.dimensions[x] for x in ['time', 'ncol', 'lev', 'ilev', 'nbnd', 'dim2']]
    for d in dimensions:
        output.createDimension(d.name, (len(d) if not d.isunlimited() else None))

    # Copy time,tine_bnds,lat,lon,hyam,hybm,area
    variables = [model.variables[x] for x in ['time', 'time_bnds', 'lat', 'lon', 'hyam', 'hybm', 'area']]
    for v in variables:
        x = output.createVariable(v.name, v.datatype, v.dimensions)
        # copy variable attributes all at once via dictionary
        output[v.name].setncatts(model[v.name].__dict__)
        # copy data
        output[v.name][:] = model[v.name][:]

    # Add lev
    lev = output.createVariable('lev', np.float64, ('lev',))
    lev.long_name = "hybrid level at midpoints (1000*(A+B))"
    lev.units = "hPa"
    lev.positive = "down"
    lev.standard_name = "atmosphere_hybrid_sigma_pressure_coordinate"
    lev.formula_terms = "a: hyam b: hybm p0: P0 ps: PS"
    lev[:] = 1000.0*(output['hyam'][:] + output['hybm'][:])

    # Create output variables

    # PS
    PS = output.createVariable('PS', np.float32, ('time', 'ncol'))
    PS.units = 'Pa'
    PS.long_name = 'Surface pressure'
    PS.standard_name = 'surface_air_pressure'

    # PRECT
    PRECT = output.createVariable('PRECT', np.float32, ('time', 'ncol'))
    PRECT.units = 'm/s'
    PRECT.standard_name = 'precipitation_flux'
    PRECT.long_name = 'Total precipitation rate (liq + ice)'

    # PRECC
    PRECC = output.createVariable('PRECC', np.float32, ('time', 'ncol'))
    PRECC.units = 'm/s'
    PRECC.long_name = 'Convective precipitation rate (liq + ice)'

    # PRECL
    PRECL = output.createVariable('PRECL', np.float32, ('time', 'ncol'))
    PRECL.units = 'm/s'
    PRECL.long_name = 'Large-scale precipitation rate (liq + ice)'

    # PRECSC
    PRECSC = output.createVariable('PRECSC', np.float32, ('time', 'ncol'))
    PRECSC.units = 'm/s'
    PRECSC.long_name = 'Convective snow rate (water equivalent)'

    # PRECSL
    PRECSL = output.createVariable('PRECSL', np.float32, ('time', 'ncol'))
    PRECSL.units = 'm/s'
    PRECSL.long_name = 'Large-scale snow rate (water equivalent)'

    # QFLX
    QFLX = output.createVariable('QFLX', np.float32, ('time', 'ncol'))
    QFLX.units = 'kg/m2/s'
    QFLX.standard_name = 'water_evapotranspiration_flux'
    QFLX.long_name = 'Surface water flux'

    # SHFLX
    SHFLX = output.createVariable('SHFLX', np.float32, ('time', 'ncol'))
    SHFLX.units = 'W/m2'
    SHFLX.standard_name = 'surface_upward_sensible_heat_flux'
    SHFLX.long_name = 'Surface sensible heat flux'

    # LHFLX
    LHFLX = output.createVariable('LHFLX', np.float32, ('time', 'ncol'))
    LHFLX.units = 'W/m2'
    LHFLX.standard_name = 'surface_upward_latent_heat_flux'
    LHFLX.long_name = 'Surface latent heat flux'

    # PSL
    PSL = output.createVariable('PSL', np.float32, ('time', 'ncol'))
    PSL.units = 'Pa'
    PSL.long_name = 'Sea level pressure'
    PSL.standard_name = 'air_pressure_at_mean_sea_level'

    # TREFHT
    TREFHT = output.createVariable('TREFHT', np.float32, ('time', 'ncol'))
    TREFHT.units = 'K'
    TREFHT.long_name = 'Reference height temperature'
    TREFHT.standard_name = 'air_temperature'

    # TS
    TS = output.createVariable('TS', np.float32, ('time', 'ncol'))
    TS.units = 'K'
    TS.long_name = 'Surface temperature (radiative)'
    TS.standard_name = 'surface_temperature'

    # TMQ
    TMQ = output.createVariable('TMQ', np.float32, ('time', 'ncol'))
    TMQ.units = 'kg/m2'
    TMQ.long_name = 'Total (vertically integrated) precipitable water'
    TMQ.standard_name = 'atmosphere_mass_content_of_water_vapor'

    # TGCLDLWP
    TGCLDLWP = output.createVariable('TGCLDLWP', np.float32, ('time', 'ncol'))
    TGCLDLWP.units = 'kg/m2'
    TGCLDLWP.long_name = 'Total grid-box cloud liquid water path'
    TGCLDLWP.standard_name = 'atmosphere_mass_content_of_cloud_liquid_water'

    # SOLIN
    SOLIN = output.createVariable('SOLIN', np.float32, ('time', 'ncol'))
    SOLIN.units = 'W/m2'
    SOLIN.long_name = 'Solar insolation'

    # FSNTOA
    FSNTOA = output.createVariable('FSNTOA', np.float32, ('time', 'ncol'))
    FSNTOA.units = 'W/m2'
    FSNTOA.long_name = 'Net solar flux at top of atmosphere'

    # FSNTOAC
    FSNTOAC = output.createVariable('FSNTOAC', np.float32, ('time', 'ncol'))
    FSNTOAC.units = 'W/m2'
    FSNTOAC.long_name = 'Clearsky net solar flux at top of atmosphere'

    # FLUT
    FLUT = output.createVariable('FLUT', np.float32, ('time', 'ncol'))
    FLUT.units = 'W/m2'
    FLUT.long_name = 'Upwelling longwave flux at top of model'

    # FLUTC
    FLUTC = output.createVariable('FLUTC', np.float32, ('time', 'ncol'))
    FLUTC.units = 'W/m2'
    FLUTC.long_name = 'Clearsky upwelling longwave flux at top of model'

    # FSNT
    FSNT = output.createVariable('FSNT', np.float32, ('time', 'ncol'))
    FSNT.units = 'W/m2'
    FSNT.long_name = 'Net solar flux at top of model'

    # FLNT
    FLNT = output.createVariable('FLNT', np.float32, ('time', 'ncol'))
    FLNT.units = 'W/m2'
    FLNT.long_name = 'Net longwave flux at top of model'

    # SWCF
    SWCF = output.createVariable('SWCF', np.float32, ('time', 'ncol'))
    SWCF.units = 'W/m2'
    SWCF.long_name = 'Shortwave cloud forcing'

    # LWCF
    LWCF = output.createVariable('LWCF', np.float32, ('time', 'ncol'))
    LWCF.units = 'W/m2'
    LWCF.long_name = 'Longwave cloud forcing'

    # FSDS
    FSDS = output.createVariable('FSDS', np.float32, ('time', 'ncol'))
    FSDS.units = 'W/m2'
    FSDS.long_name = 'Downwelling solar flux at surface'
    FSDS.standard_name = 'surface_downwelling_shortwave_flux_in_air'

    # FSDSC
    FSDSC = output.createVariable('FSDSC', np.float32, ('time', 'ncol'))
    FSDSC.units = 'W/m2'
    FSDSC.long_name = 'Clearsky downwelling solar flux at surface'
    FSDSC.standard_name = 'surface_downwelling_shortwave_flux_in_air'

    # FSNS
    FSNS = output.createVariable('FSNS', np.float32, ('time', 'ncol'))
    FSNS.units = 'W/m2'
    FSNS.long_name = 'Net solar flux at surface'

    # FSNSC
    FSNSC = output.createVariable('FSNSC', np.float32, ('time', 'ncol'))
    FSNSC.units = 'W/m2'
    FSNSC.long_name = 'Clearsky net solar flux at surface'

    # FLDS
    FLDS = output.createVariable('FLDS', np.float32, ('time', 'ncol'))
    FLDS.units = 'W/m2'
    FLDS.long_name = 'Downwelling longwave flux at surface'
    FLDS.standard_name = 'surface_downwelling_longwave_flux_in_air'

    # FLNS
    FLNS = output.createVariable('FLNS', np.float32, ('time', 'ncol'))
    FLNS.units = 'W/m2'
    FLNS.long_name = 'Net longwave flux at surface'

    # FLNSC
    FLNSC = output.createVariable('FLNSC', np.float32, ('time', 'ncol'))
    FLNSC.units = 'W/m2'
    FLNSC.long_name = 'Clearsky net longwave flux at surface'

    # T
    T = output.createVariable('T', np.float32, ('time', 'lev', 'ncol'))
    T.units = 'K'
    T.long_name = 'Temperature'
    T.standard_name = 'air_temperature'

    # Q
    Q = output.createVariable('Q', np.float32, ('time', 'lev', 'ncol'))
    Q.units = 'kg/kg'
    Q.long_name = 'Specific humidity'

    # CLDLIQ
    CLDLIQ = output.createVariable('CLDLIQ', np.float32, ('time', 'lev', 'ncol'))
    CLDLIQ.units = 'kg/kg'
    CLDLIQ.long_name = 'Cloud liquid mixing ratio'

    # CLDICE
    CLDICE = output.createVariable('CLDICE', np.float32, ('time', 'lev', 'ncol'))
    CLDICE.units = 'kg/kg'
    CLDICE.long_name = 'Cloud ice mixing ratio'

    # CLOUD
    CLOUD = output.createVariable('CLOUD', np.float32, ('time', 'lev', 'ncol'))
    CLOUD.units = 'fraction'
    CLOUD.long_name = 'Cloud fraction'

    # RELHUM
    RELHUM = output.createVariable('RELHUM', np.float32, ('time', 'lev', 'ncol'))
    RELHUM.units = 'percent'
    RELHUM.long_name = 'Relative humidity'
  
    # U
    U = output.createVariable('U', np.float32, ('time', 'lev', 'ncol'))
    U.units = 'm/s'
    U.long_name = 'Zonal wind'
    U.standard_name = 'eastward_wind'

    # V
    V = output.createVariable('V', np.float32, ('time', 'lev', 'ncol'))
    V.units = 'm/s'
    V.long_name = 'Meridional wind'
    V.standard_name = 'northward_wind'

    # OMEGA
    OMEGA = output.createVariable('OMEGA', np.float32, ('time', 'lev', 'ncol'))
    OMEGA.units = 'Pa/s'
    OMEGA.long_name = 'Vertical velocity (pressure)'
    OMEGA.standard_name = 'lagrangian_tendency_of_air_pressure'

    # Z3
    Z3 = output.createVariable('Z3', np.float32, ('time', 'lev', 'ncol'))
    Z3.units = 'm'
    Z3.long_name = 'Geopotential Height (above sea level)'
    Z3.standard_name = 'geopotential_height'

    # Loop over time levels and copy data
    for i in range(len(model.variables['time'][:])):

        # Surface pressure
        PS[i] = model['ps'][i]

        # Precipitation variables
        # Already in m/2
        PRECT[i] = model['precip_liq_surf_mass_flux'][i] + model['precip_ice_surf_mass_flux'][i]
        PRECC[i] = 0.0
        PRECL[i] = model['precip_liq_surf_mass_flux'][i] + model['precip_ice_surf_mass_flux'][i]
        PRECSC[i] = 0.0
        PRECSL[i] = model['precip_ice_surf_mass_flux'][i]

        # Surface variables
        QFLX[i] = model['surf_evap'][i]
        SHFLX[i] = model['surf_sens_flux'][i]
        PSL[i] = model['SeaLevelPressure'][i]
        TREFHT[i] = model['T_2m'][i]
        TS[i] = model['surf_radiative_T'][i]

        # Other 2d variables
        TMQ[i] = model['VapWaterPath'][i]
        TGCLDLWP[i] = model['LiqWaterPath'][i] + model['RainWaterPath'][i]
        
        # TOA radiation
        SOLIN[i] = model['SW_flux_dn_at_top_of_model'][i]
        FSNTOA[i] = model['SW_flux_dn_at_top_of_model'][i] - model['SW_flux_up_at_top_of_model'][i]
        FSNTOAC[i] = model['SW_flux_dn_at_top_of_model'][i] - model['SW_clrsky_flux_up_at_top_of_model'][i]
        FLUT[i] = model['LW_flux_up_at_top_of_model'][i]
        FLUTC[i] = model['LW_clrsky_flux_up_at_top_of_model'][i]
        FSNT[i] = model['SW_flux_dn_at_top_of_model'][i] - model['SW_flux_up_at_top_of_model'][i]
        FLNT[i] = model['LW_flux_up_at_top_of_model'][i]
        SWCF[i] = - model['SW_flux_up_at_top_of_model'][i] + model['SW_clrsky_flux_up_at_top_of_model'][i]
        LWCF[i] = model['LW_clrsky_flux_up_at_top_of_model'][i] - model['LW_flux_up_at_top_of_model'][i]

        # Surface radiation
        FSDS[i] = model['SW_flux_dn_at_bot_of_model'][i]
        FSDSC[i] = model['SW_clrsky_flux_dn_at_bot_of_model'][i]
        FSNS[i] = model['SW_flux_dn_at_bot_of_model'][i] - model['SW_flux_up_at_bot_of_model'][i]
        FSNSC[i] = model['SW_clrsky_flux_dn_at_bot_of_model'][i] - model['SW_clrsky_flux_up_at_bot_of_model'][i]
        FLDS[i] = model['LW_flux_dn_at_bot_of_model'][i]
        FLNS[i] = - model['LW_flux_dn_at_bot_of_model'][i] + model['LW_flux_up_at_bot_of_model'][i]
        FLNSC[i] = - model['LW_clrsky_flux_dn_at_bot_of_model'][i] + model['LW_clrsky_flux_up_at_bot_of_model'][i]

        # 3d variables
        T[i,:,:] = np.swapaxes(model['T_mid'][i,:,:],0,1)
        Q[i,:,:] = np.swapaxes(model['qv'][i,:,:],0,1)
        CLDLIQ[i,:,:] = np.swapaxes(model['qc'][i,:,:],0,1)
        CLDICE[i,:,:] = np.swapaxes(model['qi'][i,:,:],0,1)
        CLOUD[i,:,:] = np.swapaxes(model['cldfrac_tot_for_analysis'][i,:,:],0,1)
        RELHUM[i,:,:] = np.swapaxes(model['RelativeHumidity'][i,:,:],0,1) * 100.0
        U[i,:,:] = np.swapaxes(np.squeeze(model['U'][i,:,:]),0,1)
        V[i,:,:] = np.swapaxes(np.squeeze(model['V'][i,:,:]),0,1)
        OMEGA[i,:,:] = np.swapaxes(model['omega'][i,:,:],0,1)
        Z3[i,:,:] =model['z_mid'][i,:,:]

    # Close current input, output files
    model.close()
    output.close()
