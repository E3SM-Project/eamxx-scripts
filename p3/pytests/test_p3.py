#!/usr/bin/env python
"""
This script contains first steps at exercising P3 from 
python. It assumes you've already compiled 
p3.so using f2py.
"""

#IMPORT STUFF:
#===============
import p3
import numpy as np
import pylab as pl


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=
def plot_p3_output(case,pres,qc,nc,qr,nr,th_old,th,qv_old,qv,\
                   qitot,qirim,nitot,birim,ssat,qc_init,nc_init,\
                   qr_init,nr_init,th_old_init,th_init,qv_old_init,\
                   qv_init,qitot_init,qirim_init,nitot_init,\
                   birim_init,ssat_init,prt_liq,prt_sol,\
                   diag_ze,diag_effc,diag_effi,diag_vmi,diag_di,\
                   diag_rhoi,diag_3d,prt_drzl,prt_rain,prt_crys,prt_snow,\
                   prt_grpl,prt_pell,prt_hail,prt_sndp):
    """
    This function just takes in all the variables that p3 changes and plots them.
    The input variable "case" is a string to specify in each figure title which run 
    is being tested. In order to see the plots from this call, 'pl.show()' must be issued.
    """

    # --- LOOP OVER INOUT VARS ---
    inout_vars=['qc','nc','qr','nr','th_old','th','qv_old','qv','qitot','qirim','nitot','birim','ssat']

    for k in range(len(inout_vars)):
        #MAKE A NEW FIGURE IF NEEDED
        if k%4==0:
            pl.figure()

        pl.subplot(2,2,k%4+1)
        pl.plot(eval(inout_vars[k]+'_init').squeeze(),pres.squeeze()/100.,'b-')
        pl.plot(eval(inout_vars[k]).squeeze(),pres.squeeze()/100.,'b--')
        pl.gca().invert_yaxis()
        pl.title(case+': '+inout_vars[k])

    # --- LOOP OVER OUTPUT VARS ---
    # removed 'diag_2d' because it has shape [1,1] which isn't a scalar but can't be plotted.
    out_vars=['prt_liq','prt_sol','diag_ze','diag_effc','diag_effi','diag_vmi','diag_di',\
              'diag_rhoi','diag_3d','prt_drzl','prt_rain','prt_crys','prt_snow',\
              'prt_grpl','prt_pell','prt_hail','prt_sndp']

    cnt=-1
    for k in range(len(inout_vars)):
        if eval(out_vars[k]+'.shape')==(1,): #if scalar
            print case+': '+out_vars[k]+' = %0.3e'%eval(out_vars[k])

        else:
            cnt+=1

            #MAKE A NEW FIGURE IF NEEDED
            if cnt%4==0:
                pl.figure()

            pl.subplot(2,2,cnt%4+1)
            pl.plot(eval(out_vars[k]).squeeze(),pres.squeeze()/100.,'b-')
            pl.gca().invert_yaxis()
            pl.title(case+': '+out_vars[k])

    return

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#===========================
#RUN INITIALIZATION ROUTINE:
#===========================
#data_dir='/project/projectdirs/acme/inputdata/atm/cam/physprops'
data_dir='/g/g11/caldwep/scream/p3'
ncat=1 #number of ice categories. We will always use 1.

#This just primes the p3 module so p3_main works right
p3.micro_p3.p3_init(data_dir,ncat)

#===========================
# LOOP OVER TEST CASES:
#===========================

cases = ['mixed'] #inactive,ice,liq,mixed
for case in cases:

    #LOAD ALL INPUT DATA NEEDED BY P3
    #---------------------------
    exec('from '+case+'_case_data import *')

    #MAKE COPIES OF INOUT VARIABLES BEFORE THEY GET UPDATED:
    #---------------------------
    qc_init=qc.copy()
    nc_init=nc.copy()
    qr_init=qr.copy()
    nr_init=nr.copy()
    th_old_init=th_old.copy()
    th_init=th.copy()
    qv_init=qv.copy()
    qv_old_init=qv_old.copy()
    qitot_init=qitot.copy()
    qirim_init=qirim.copy()
    nitot_init=nitot.copy()
    birim_init=birim.copy()
    ssat_init=ssat.copy()

    #CONFIRM ALL INOUT VARS HAVE FORTRAN BYTE ORDER:
    #---------------------------
    vars=['qc','nc','qr','nr','th_old','th','qv_old','qv','qitot','qirim',\
          'nitot','birim','ssat','uzpl','pres','dzq']
    for var in vars:
        exec('isF=np.isfortran('+var+')')
        if not isF:
            raise Exception(var+' is not Fortran contiguous')
    
    
    #ACTUALLY CALL P3_MAIN:
    #---------------------------
    print 'Running '+case

    prt_liq,prt_sol,diag_ze,diag_effc,diag_effi,diag_vmi,diag_di,diag_rhoi,diag_2d,diag_3d,prt_drzl,prt_rain,prt_crys,prt_snow,prt_grpl,prt_pell,prt_hail,prt_sndp = p3.micro_p3.p3_main(qc,nc,qr,nr,th_old,th,qv_old,qv,dt,qitot,qirim,nitot,birim,ssat,uzpl,pres,dzq,it,its,ite,kts,kte,n_diag_2d,n_diag_3d,log_predictnc,typediags_on,model,[ncat])

    print 'Completed '+case

    #PLOT THE RESULTS:
    #---------------------------
    plot_p3_output(case,pres,qc,nc,qr,nr,th_old,th,qv_old,qv,\
                   qitot,qirim,nitot,birim,ssat,qc_init,nc_init,\
                   qr_init,nr_init,th_old_init,th_init,qv_old_init,\
                   qv_init,qitot_init,qirim_init,nitot_init,\
                   birim_init,ssat_init,prt_liq,prt_sol,\
                   diag_ze,diag_effc,diag_effi,diag_vmi,diag_di,\
                   diag_rhoi,diag_3d,prt_drzl,prt_rain,prt_crys,prt_snow,\
                   prt_grpl,prt_pell,prt_hail,prt_sndp)


#=============================
#NOW DISPLAY ALL THE FIGURES:
#=============================
pl.show()

