#!/usr/bin/env python3
from subprocess import call
import os
datestring = '20220315'
#input_file = "/global/cscratch1/sd/bhillma/scream/cases/aarondonahue/scream_p3_sa_input_from_f90.ne30_ne30.F2010-SCREAM-LR.cori-knl_intel.32x32x2.20220311-1614/run/scream_p3_sa_input_from_f90.ne30_ne30.F2010-SCREAM-LR.cori-knl_intel.32x32x2.20220311-1614.eam.h1.0001-01-01-00000.nc"
#output_file = f'/global/cscratch1/sd/bhillma/scream/data/init/ne30np4L72/{datestring}/screami_ne30np4L72_{datestring}.nc'
input_file = '/global/cscratch1/sd/bhillma/scream/cases/aarondonahue/scream_p3_sa_input_from_f90.ne4_ne4.F2010-SCREAM-LR.cori-knl_intel.1x32x2.20220315-1829/run/scream_p3_sa_input_from_f90.ne4_ne4.F2010-SCREAM-LR.cori-knl_intel.1x32x2.20220315-1829.eam.h1.0001-01-01-00000.nc'
output_file = f'/global/cscratch1/sd/bhillma/scream/data/init/ne4np4L72/{datestring}/screami_ne4np4L72_{datestring}.nc'
os.makedirs(os.path.dirname(output_file), exist_ok=True)
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
rename_vars_in_orig=(
    "H2O_inRAD",
    "CO2_inRAD",
    "O3_inRAD",
    "N2O_inRAD",
    "CO_inRAD",
    "CH4_inRAD",
    "O2_inRAD",
    "N2_inRAD",
)
rename_vars_in_new=(
    "h2o",
    "co2",
    "o3",
    "n2o",
    "co",
    "ch4",
    "o2",
    "n2",
)

# Leave SPA out of it; separate input file
required_vars = homme_vars + p3_vars + shoc_vars + rad_vars

# Create a mapping from SCREAMv0 output names to SCREAMv1 input names
screamv0_to_screamv1_names = {v: v.split('_in')[0] for v in required_vars}

# And some of these need to be renamed still
for v0, v1 in zip(rename_vars_in_orig, rename_vars_in_new): 
    screamv0_to_screamv1_names[v0] = v1

# Only keep one version of each variable added?
#unique_names = set(v1 for v0, v1 in screamv0_to_screamv1_names.items())
#print(unique_names)
#exit()

# Make file with all these variables
required_vars_string = ','.join(required_vars)
print('Adding variables from input file to output file...')
call(f'ncks -O -v {required_vars_string} {input_file} {output_file}'.split(' '))

# Rename
print('Rename variables from SCREAMv0 output to SCREAMv1 input names...')
rename_string = ' '.join([f'-v {v0},{v1}' for v0, v1 in screamv0_to_screamv1_names.items() if v0 != v1])
call(f'ncrename -O {rename_string} {output_file}'.split(' '))

# Fix dimension ording
print('Fix dimension ordering...')
call(f'ncpdq -O -a time,ncol,dim2,lev {output_file} {output_file}'.split(' '))
