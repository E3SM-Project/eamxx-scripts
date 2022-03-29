#!/usr/bin/env python3
from subprocess import call
import os
from datetime import date

# Define input and output files
datestring = date.today().strftime('%Y%m%d')
res = "ne120np4L72"
#input_file = "/global/cscratch1/sd/bhillma/scream/cases/aarondonahue/scream_p3_sa_input_from_f90.ne30_ne30.F2010-SCREAM-LR.cori-knl_intel.32x32x2.20220311-1614/run/scream_p3_sa_input_from_f90.ne30_ne30.F2010-SCREAM-LR.cori-knl_intel.32x32x2.20220311-1614.eam.h1.0001-01-01-00000.nc"
input_file = "/global/cscratch1/sd/bhillma/scream/cases/aarondonahue/scream_p3_sa_input_from_f90.ne120_r0125_oRRS18to6v3.F2010-SCREAM-LR.cori-knl_intel.384x32x2.20220328/run/scream_p3_sa_input_from_f90.ne120_r0125_oRRS18to6v3.F2010-SCREAM-LR.cori-knl_intel.384x32x2.20220328.eam.h1.0001-01-01-00000.nc"
output_file = f'/global/cscratch1/sd/bhillma/scream/data/init/{res}/screami_{res}_{datestring}.nc'
#input_file = '/global/cscratch1/sd/bhillma/scream/cases/aarondonahue/scream_p3_sa_input_from_f90.ne4_ne4.F2010-SCREAM-LR.cori-knl_intel.1x32x2.20220315-1829/run/scream_p3_sa_input_from_f90.ne4_ne4.F2010-SCREAM-LR.cori-knl_intel.1x32x2.20220315-1829.eam.h1.0001-01-01-00000.nc'
#output_file = f'/global/cscratch1/sd/bhillma/scream/data/init/ne4np4L72/{datestring}/screami_ne4np4L72_{datestring}.nc'

# Make sure directory for output file exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Define required vars
spa_vars=(
    "CCN3",
    "AER_G_SW",
    "AER_SSA_SW",
    "AER_TAU_SW",
    "AER_TAU_LW",
)
homme_vars=(
    "hyam",
    "hybm",
    "hyai",
    "hybi",
    "P0",
    "dp_inHOMME" ,
    "phi_int_inHOMME",
    "qv_inHOMME",
    "v_inHOMME",
    "vtheta_dp_inHOMME",
    "w_i_inHOMME",
    "ps_inHOMME",
)
p3_vars=(
    "T_mid_inP3",
    "p_mid_inP3",
    "pseudo_density_inP3",
    "qc_inP3",
    "z_int_inP3",
    "inv_qc_relvar_inP3",
    "bm_inP3",
    "nc_inP3",
    "ni_inP3",
    "nr_inP3",
    "qi_inP3",
    "qm_inP3",
    "qr_inP3",
    "cldfrac_tot_inP3",
    "T_prev_micro_step_inP3",
    "nc_activated_inP3",
    "nc_nuceat_tend_inP3",
    "ni_activated_inP3",
    "qv_prev_micro_step_inP3",
)
shoc_vars=(
    #"T_mid_inSHOC",  # Already added in P3
    "cldfrac_liq_inSHOC",
    "eddy_diff_mom_inSHOC",
    "horiz_winds_inSHOC",
    "host_dx_inSHOC",
    "host_dy_inSHOC",
    "omega_inSHOC",
    "p_int_inSHOC",
    #"p_mid_inSHOC",
    "phis_inSHOC",
    #"pseudo_density_inSHOC",
    #"qc_inSHOC",
    #"qv_inSHOC",
    "sgs_buoy_flux_inSHOC",
    "surf_latent_flux_inSHOC",
    "surf_sens_flux_inSHOC",
    "surf_mom_flux_inSHOC",
    "tke_inSHOC",
    #"z_int_inSHOC",
    "z_mid_inSHOC",
)
rad_vars=(
    #"T_mid_inRAD",  # Already added in P3
    "eff_radius_qc_inRAD",
    "eff_radius_qi_inRAD",
    #"cldfrac_tot_inRAD",
    "sfc_alb_dir_vis_inRAD",
    "sfc_alb_dif_vis_inRAD",
    "sfc_alb_dir_nir_inRAD",
    "sfc_alb_dif_nir_inRAD",
    "surf_lw_flux_up_inRAD",
    "H2O_inRAD",
    "CO2_inRAD",
    "O3_inRAD",
    "N2O_inRAD",
    "CO_inRAD",
    "CH4_inRAD",
    "O2_inRAD",
    "N2_inRAD",
)
rename_vars = {
    "H2O_inRAD": 'h2o',
    "CO2_inRAD": 'co2',
    "O3_inRAD" : 'o3' ,
    "N2O_inRAD": 'n2o',
    "CO_inRAD" : 'co' ,
    "CH4_inRAD": 'ch4',
    "O2_inRAD" : 'o2' ,
    "N2_inRAD" : 'n2' ,
}

# Full list of required variables from v0 output file. Note that we are leaving
# SPA out, since the SPA data is in a separate input file.
required_vars = homme_vars + p3_vars + shoc_vars + rad_vars

# Create a mapping from SCREAMv0 output names to SCREAMv1 input names; using the
# branch aarondonahue/scream_p3_sa_input_from_f90, these all have names
# <var>_in<PROCESS>; i.e., qv_inHOMME. In the future, we might give these
# standard names. This naming convention was adopted in anticipation that we
# might use this to create inputs and outputs from each individual process.
screamv0_to_screamv1_names = {v: v.split('_in')[0] for v in required_vars}

# And some of these need to be renamed still
for v0, v1 in rename_vars.items():
    screamv0_to_screamv1_names[v0] = v1

# Copy all the required fields to a new file
required_vars_string = ','.join(required_vars)
print('Adding variables from input file to output file...')
call(f'ncks -O -v {required_vars_string} {input_file} {output_file}'.split(' '))

# Rename variables
print('Rename variables from SCREAMv0 output to SCREAMv1 input names...')
rename_string = ' '.join([f'-v {v0},{v1}' for v0, v1 in screamv0_to_screamv1_names.items() if v0 != v1])
call(f'ncrename -O {rename_string} {output_file}'.split(' '))

# Fix dimension ordering
print('Fix dimension ordering...')
call(f'ncpdq -O -a time,ncol,dim2,lev {output_file} {output_file}'.split(' '))

# Append another var identical to hybm but called pref_mid for SHOC
# TODO: why does SHOC not just use hybm?
print('Add a hybm->pref_mid var for SHOC...')
call(f'ncks -O -v hybm {input_file} pref_mid.nc'.split(' '))
call(f'ncrename -O -v hybm,pref_mid pref_mid.nc'.split(' '))
call(f'ncks -A -v pref_mid pref_mid.nc {output_file}'.split(' '))
os.remove('pref_mid.nc')
