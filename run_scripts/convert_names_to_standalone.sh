#!/bin/sh

# Purpose: translate fields from SCREAMv0 history files to SCREAMv1 initial
# conditions
#
# Authors: Aaron Donahue (original version)
#          Benjamin Hillman (minor edits)

# Exit script if any command returns non-zero exit code
set -e

module load nco

res="ne30np4"  # ne2np4, ne4np4

dirout="${CSCRATCH}/scream/init_files/$res"
mkdir -p ${dirout}

# Path to SCREAMv0 history file from which to pull values
#filein="cori-knl.scream.standalone_init_data.ne2_ne2.eam.h1.0001-01.nc"
#filein="cori-knl.scream.standalone_init_data.ne4_ne4.eam.h1.0001-01.nc"
#filein="cori-knl.scream.standalone_init_data.ne30_ne30.eam.h1.0001-01.nc"
filein="/global/cscratch1/sd/bhillma/scream/cases/add-aquaplanet.ne30_ne30.F-SCREAM-LR-AQP1.cori-knl_intel.32x32x2.20211111-1350/run/add-aquaplanet.ne30_ne30.F-SCREAM-LR-AQP1.cori-knl_intel.32x32x2.20211111-1350.eam.h1.0001-01-01-00000.nc"

# SPA needs additional input file we do not have access to right now, so disable
# for the time being
do_spa=0
spain="spa_init_$res.nc"

# Standalone initial condition files we will generate
#file_forout="quartz.scream.PR_setup_init_vals.ne4_ne4.eam.h2.0001-01.nc"
homme_file="$dirout/homme_init_$res.nc"
p3_file="$dirout/p3_init_$res.nc"
shoc_file="$dirout/shoc_init_$res.nc"
rad_file="$dirout/rrtmgp_init_$res.nc"
shoc_cld_p3_file="$dirout/shoc_cld_p3_init_$res.nc"
shoc_cld_p3_rad_file="$dirout/shoc_cld_p3_rrtmgp_init_$res.nc"
homme_shoc_cld_p3_rad_file="$dirout/homme_shoc_cld_p3_rrtmgp_init_$res.nc"
homme_shoc_cld_spa_p3_rad_file="$dirout/homme_shoc_cld_spa_p3_rrtmgp_init_$res.nc"
shoc_cld_spa_p3_rad_file="$dirout/shoc_cld_spa_p3_rrtmgp_init_$res.nc"
spa_rad_file="$dirout/spa_rrtmgp_init_$res.nc"

# Baseline output files
p3_file_out="$dirout/p3_$res_baseline.nc"
shoc_file_out="$dirout/shoc_$res_baseline.nc"

spa_vars=(
"PS"
"CCN3"
"AER_G_SW"
"AER_SSA_SW"
"AER_TAU_SW"
"AER_TAU_LW"
)

homme_vars=(
"dp_inHOMME" \
"phi_int_inHOMME" \
"qv_inHOMME" \
"v_inHOMME" \
"vtheta_dp_inHOMME" \
"w_i_inHOMME" \
"ps_inHOMME" \
)

p3_vars=(
"T_mid_inP3" \
"p_mid_inP3" \
"pseudo_density_inP3" \
"qc_inP3" \
"qv_inP3" \
"z_int_inP3" \
"inv_qc_relvar_inP3" \
"bm_inP3" \
"nc_inP3" \
"ni_inP3" \
"nr_inP3" \
"qi_inP3" \
"qm_inP3" \
"qr_inP3" \
"cldfrac_tot_inP3" \
"T_prev_micro_step_inP3" \
"nc_activated_inP3" \
"nc_nuceat_tend_inP3" \
"ni_activated_inP3" \
"qv_prev_micro_step_inP3" \
)

shoc_vars=(
"T_mid_inSHOC" \
"cldfrac_liq_inSHOC" \
"eddy_diff_mom_inSHOC" \
"horiz_winds_inSHOC" \
"host_dx_inSHOC" \
"host_dy_inSHOC" \
"omega_inSHOC" \
"p_int_inSHOC" \
"p_mid_inSHOC" \
"phis_inSHOC" \
"pseudo_density_inSHOC" \
"qc_inSHOC" \
"qv_inSHOC" \
"sgs_buoy_flux_inSHOC" \
"surf_latent_flux_inSHOC" \
"surf_sens_flux_inSHOC" \
"surf_mom_flux_inSHOC" \
"tke_inSHOC" \
"z_int_inSHOC" \
"z_mid_inSHOC" \
)


shoc_2_rad_vars=(
"p_int_inSHOC" \
)
p3_2_rad_vars=(
"p_mid_inP3" \
"pseudo_density_inP3" \
"qc_inP3" \
"qv_inP3" \
"qi_inP3" \
)
rad_vars=(
"T_mid_inRAD" \
"eff_radius_qc_inRAD" \
"eff_radius_qi_inRAD" \
"cldfrac_tot_inRAD" \
"sfc_alb_dir_vis_inRAD" \
"sfc_alb_dif_vis_inRAD" \
"sfc_alb_dir_nir_inRAD" \
"sfc_alb_dif_nir_inRAD" \
"surf_lw_flux_up_inRAD" \
"H2O_inRAD" \
"CO2_inRAD" \
"O3_inRAD" \
"N2O_inRAD" \
"CO_inRAD" \
"CH4_inRAD" \
"O2_inRAD" \
"N2_inRAD" \
)
#"p_mid_inRAD" \
#"p_int_inRAD" \
#"qc_inRAD" \
#"qi_inRAD" \
#"pseudo_density_inRAD" \ This one is suspicious, and wrong - all zero in src file.

p3_vars_out=(
    "T_mid_outP3" \
    "qc_outP3" \
    "qv_outP3" \
    "bm_outP3" \
    "nc_outP3" \
    "ni_outP3" \
    "nr_outP3" \
    "qi_outP3" \
    "qm_outP3" \
    "qr_outP3" \
    "T_prev_micro_step_outP3" \
    "qv_prev_micro_step_outP3" \
    "eff_radius_qc_outP3" \
    "eff_radius_qi_outP3" \
    "liq_ice_exchange_outP3" \
    "vap_liq_exchange_outP3" \
    "vap_ice_exchange_outP3"
)
shoc_vars_out=(
    "T_mid_outSHOC" \
    "cldfrac_liq_outSHOC" \
    "eddy_diff_mom_outSHOC" \
    "horiz_winds_outSHOC" \
    "qc_outSHOC" \
    "qv_outSHOC" \
    "sgs_buoy_flux_outSHOC" \
    "tke_outSHOC" \
    "inv_qc_relvar_outSHOC" \
    "pbl_height_outSHOC"
)

rename_vars_in_orig=(
"H2O" \
"CO2" \
"O3" \
"N2O" \
"CO" \
"CH4" \
"O2" \
"N2" \
)
rename_vars_in_new=(
"h2o" \
"co2" \
"o3" \
"n2o" \
"co" \
"ch4" \
"o2" \
"n2" \
)
num_rename_vars_in=${#rename_vars_in_orig[@]}
#---------------------------------------------------------------------------------------#
#HOMME - only
echo "HOMME-only..."
file_out=$homme_file
if [ ! -e ${file_out} ]; then
    for v in "${homme_vars[@]}"
    do
      ncks -h -A -v $v $filein $file_out
      v_strip=${v%"_inHOMME"}
      ncrename -h -v $v,$v_strip $file_out
    done
    ncks -h -A -v lat $filein $file_out
    ncks -h -A -v lon $filein $file_out
    ncpdq -h -a ncol,lev $file_out tmp.nc
    mv tmp.nc $file_out
    ncpdq -h -a ncol,ilev $file_out tmp.nc
    mv tmp.nc $file_out
    ncpdq -h -a ncol,dim2 $file_out tmp.nc
    mv tmp.nc $file_out
else
    echo "${file_out} exists, skipping."
fi
#---------------------------------------------------------------------------------------#
#P3 - only
echo "P3-only..."
file_out=$p3_file
if [ ! -e ${file_out} ]; then
    for v in "${p3_vars[@]}"
    do
      ncks -h -A -v $v $filein $file_out
      v_strip=${v%"_inP3"}
      ncrename -h -v $v,$v_strip $file_out
    done
    ncks -h -A -v lat $filein $file_out
    ncks -h -A -v lon $filein $file_out
    ncpdq -h -a ncol,lev $file_out tmp.nc
    mv tmp.nc $file_out
    ncpdq -h -a ncol,ilev $file_out tmp.nc
    mv tmp.nc $file_out
    ncpdq -h -a ncol,dim2 $file_out tmp.nc
    mv tmp.nc $file_out
else
    echo "${file_out} exists, skipping."
fi
#--------------------------------------#
#ASD file_out=$p3_file_out
#ASD for v in "${p3_vars_out[@]}"
#ASD do
#ASD   ncks -A -v $v $file_forout $file_out
#ASD   v_strip=${v%"_outP3"}
#ASD   if [ $v_strip == "liq_ice_exchange" ]
#ASD   then
#ASD     v_strip="micro_liq_ice_exchange"
#ASD   elif [ $v_strip == "vap_liq_exchange" ]
#ASD   then
#ASD     v_strip="micro_vap_liq_exchange"
#ASD   elif [ $v_strip == "vap_ice_exchange" ]
#ASD   then
#ASD     v_strip="micro_vap_ice_exchange"
#ASD   fi
#ASD   ncrename -h -v $v,$v_strip $file_out
#ASD done
#ASD ncks -A -v lat $filein $file_out
#ASD ncks -A -v lon $filein $file_out
#ASD ncpdq -h -a ncol,lev $file_out tmp.nc
#ASD mv tmp.nc $file_out
#ASD ncpdq -h -a ncol,ilev $file_out tmp.nc
#ASD mv tmp.nc $file_out
#ASD ncpdq -h -a ncol,dim2 $file_out tmp.nc
#ASD mv tmp.nc $file_out

#---------------------------------------------------------------------------------------#
## SHOC - only
echo "SHOC-only..."
file_out=$shoc_file
if [ ! -e ${file_out} ]; then
    for v in "${shoc_vars[@]}"
    do
      ncks -h -A -v $v $filein $file_out
      v_strip=${v%"_inSHOC"}
      ncrename -h -v $v,$v_strip $file_out
    done
    ncks -h -A -v lat $filein $file_out
    ncks -h -A -v lon $filein $file_out
    ncks -h -A -v hybm $filein $file_out
    ncrename -h -v hybm,pref_mid $file_out
    ncpdq -h -a ncol,lev $file_out tmp.nc
    mv tmp.nc $file_out
    ncpdq -h -a ncol,ilev $file_out tmp.nc
    mv tmp.nc $file_out
    ncpdq -h -a ncol,dim2 $file_out tmp.nc
    mv tmp.nc $file_out
else
    echo "${file_out} exists, skipping."
fi
#--------------------------------------#
#ASD file_out=$shoc_file_out
#ASD for v in "${shoc_vars_out[@]}"
#ASD do
#ASD   ncks -A -v $v $file_forout $file_out
#ASD   v_strip=${v%"_outSHOC"}
#ASD   ncrename -h -v $v,$v_strip $file_out
#ASD done
#ASD ncks -A -v lat $filein $file_out
#ASD ncks -A -v lon $filein $file_out
#ASD ncks -A -v hybm $filein $file_out
#ASD ncrename -h -v hybm,pref_mid $file_out
#ASD ncpdq -h -a ncol,lev $file_out tmp.nc
#ASD mv tmp.nc $file_out
#ASD ncpdq -h -a ncol,ilev $file_out tmp.nc
#ASD mv tmp.nc $file_out
#ASD ncpdq -h -a ncol,dim2 $file_out tmp.nc
#ASD mv tmp.nc $file_out

#---------------------------------------------------------------------------------------#
## RAD - only
echo "RAD-only..."
file_out=$rad_file
for v in "${rad_vars[@]}"
do
  ncks -h -A -v $v $filein $file_out
  v_strip=${v%"_inRAD"}
  ncrename -h -v $v,$v_strip $file_out
done
# We need qv as well, TODO have a qv_inRAD variable.  For now, take the one from p3.
# pseudo_density from rad seems wrong too, fixing that by using p3 value.
for v in "${p3_2_rad_vars[@]}"
do
  ncks -h -A -v $v $filein $file_out
  v_strip=${v%"_inP3"}
  ncrename -h -v $v,$v_strip $file_out
done
for v in "${shoc_2_rad_vars[@]}"
do
  ncks -h -A -v $v $filein $file_out
  v_strip=${v%"_inSHOC"}
  ncrename -h -v $v,$v_strip $file_out
done
# Recall, TODO is to replace the above with an inRAD version.
ncks -h -A -v lat $filein $file_out
ncks -h -A -v lon $filein $file_out
ncpdq -h -a ncol,lev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,ilev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,ngas $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a lev,ngas $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,swband $file_out tmp.nc
mv tmp.nc $file_out
## rename gases to be lowercase
echo "Rename gases"
for ii in $(seq 0 $(expr $num_rename_vars_in - 1))
do
  v=${rename_vars_in_orig[$ii]}
  v_new=${rename_vars_in_new[$ii]}
  ncrename -h -O -v $v,$v_new $file_out
done
#---------------------------------------------------------------------------------------#
## SPA+RAD
echo "SPA+RAD..."
file_out=$spa_rad_file
cp $rad_file $file_out
for v in "${spa_vars[@]}"
do
  v_tmp="${v}_beg"
  if [ $do_spa -eq 1 ]; then ncks -h -A -v $v_tmp $spain $file_out; fi
  v_tmp="${v}_end"
  if [ $do_spa -eq 1 ]; then ncks -h -A -v $v_tmp $spain $file_out; fi
done
ncpdq -h -a ncol,lev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,ilev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,dim2 $file_out tmp.nc
mv tmp.nc $file_out
#---------------------------------------------------------------------------------------#
## SHOC+CLDFRAC+P3+RAD
echo "SHOC+CLDFRAC+P3+RAD..."
vars_set=()
file_out=$shoc_cld_p3_rad_file
# SHOC
for v in "${shoc_vars[@]}"
do
  ncks -h -A -v $v $filein $file_out
  v_strip=${v%"_inSHOC"}
  ncrename -h -O -v $v,$v_strip $file_out
  vars_set+=("$v_strip")
done
# P3 - note cloud frac doesn't have individual outputs
for v in "${p3_vars[@]}"
do
  v_strip=${v%"_inP3"}
  ix=$( printf "%s\n" "${vars_set[@]}" | grep -n -m 1 "^${v_strip}$" | cut -d ":" -f1 )
  if [[ -z $ix ]]
    then
      ncks -h -A -v $v $filein $file_out
      ncrename -h -v $v,$v_strip $file_out
      vars_set+=("$v_strip")
  fi
done
# RAD
echo "Rename rad vars..."
for v in "${rad_vars[@]}"
do
  v_strip=${v%"_inRAD"}
  ix=$( printf "%s\n" "${vars_set[@]}" | grep -n -m 1 "^${v_strip}$" | cut -d ":" -f1 )
  if [[ -z $ix ]]
    then
      ncks -h -A -v $v $filein $file_out
      ncrename -h -v $v,$v_strip $file_out
  fi
done
ncks -h -A -v lat $filein $file_out
ncks -h -A -v lon $filein $file_out
ncks -h -A -v hybm $filein $file_out
ncrename -h -v hybm,pref_mid $file_out
ncpdq -h -a ncol,lev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,ilev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,dim2 $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,ngas $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a lev,ngas $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,swband $file_out tmp.nc
mv tmp.nc $file_out
## rename gases to be lowercase
for ii in $(seq 0 $(expr $num_rename_vars_in - 1)) 
do
  v=${rename_vars_in_orig[$ii]}
  v_new=${rename_vars_in_new[$ii]}
  ncrename -h -O -v $v,$v_new $file_out
done
#---------------------------------------------------------------------------------------#
## HOMME+SHOC+CLDFRAC+P3+RAD
echo "HOMME+SHOC+CLDFRAC+P3+RAD..."
file_out=$homme_shoc_cld_p3_rad_file
cp $shoc_cld_p3_rad_file $file_out
# Add HOMME
for v in "${homme_vars[@]}"
do
  v_strip=${v%"_inHOMME"}
  ix=$( printf "%s\n" "${vars_set[@]}" | grep -n -m 1 "^${v_strip}$" | cut -d ":" -f1 )
  if [[ -z $ix ]]
    then
      ncks -h -A -v $v $filein $file_out
      ncrename -h -v $v,$v_strip $file_out
  fi
done
ncpdq -h -a ncol,lev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,ilev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,dim2 $file_out tmp.nc
mv tmp.nc $file_out
#---------------------------------------------------------------------------------------#
## HOMME+SHOC+CLDFRAC+SPA+P3+RAD
echo "HOMME+SHOC+CLDFRAC+SPA+P3+RAD..."
file_out=$homme_shoc_cld_spa_p3_rad_file
cp $homme_shoc_cld_p3_rad_file $file_out
for v in "${spa_vars[@]}"
do
  v_tmp="${v}_beg"
  if [ $do_spa -eq 1 ]; then ncks -h -A -v $v_tmp $spain $file_out; fi
  v_tmp="${v}_end"
  if [ $do_spa -eq 1 ]; then ncks -h -A -v $v_tmp $spain $file_out; fi
done
ncpdq -h -a ncol,lev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,ilev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,dim2 $file_out tmp.nc
mv tmp.nc $file_out
#---------------------------------------------------------------------------------------#
## SHOC+CLDFRAC+SPA+P3+RAD
echo "SHOC+CLDFRAC+SPA+P3+RAD..."
file_out=$shoc_cld_spa_p3_rad_file
cp $shoc_cld_p3_rad_file $file_out
for v in "${spa_vars[@]}"
do
  v_tmp="${v}_beg"
  if [ $do_spa -eq 1 ]; then ncks -h -A -v $v_tmp $spain $file_out; fi
  v_tmp="${v}_end"
  if [ $do_spa -eq 1 ]; then ncks -h -A -v $v_tmp $spain $file_out; fi
done
ncpdq -h -a ncol,lev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,ilev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,dim2 $file_out tmp.nc
mv tmp.nc $file_out

#---------------------------------------------------------------------------------------#
## SHOC+CLDFRAC+P3
echo "SHOC+CLDFRAC+P3..."
vars_set=()
file_out=$shoc_cld_p3_file
for v in "${shoc_vars[@]}"
do
  ncks -h -A -v $v $filein $file_out
  v_strip=${v%"_inSHOC"}
  ncrename -h -O -v $v,$v_strip $file_out
  vars_set+=("$v_strip")
done
for v in "${p3_vars[@]}"
do
  v_strip=${v%"_inP3"}
  ix=$( printf "%s\n" "${vars_set[@]}" | grep -n -m 1 "^${v_strip}$" | cut -d ":" -f1 )
  if [[ -z $ix ]]
    then
      ncks -h -A -v $v $filein $file_out
      ncrename -h -v $v,$v_strip $file_out
  fi
done
ncks -h -A -v lat $filein $file_out
ncks -h -A -v lon $filein $file_out
ncks -h -A -v hybm $filein $file_out
ncrename -h -v hybm,pref_mid $file_out
ncpdq -h -a ncol,lev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,ilev $file_out tmp.nc
mv tmp.nc $file_out
ncpdq -h -a ncol,dim2 $file_out tmp.nc
mv tmp.nc $file_out
