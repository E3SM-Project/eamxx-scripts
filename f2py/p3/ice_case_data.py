"""
This module just contains all the variables needed to run p3_main in
a configuration where only ice microphysical processes should be active. 
This is done by setting ice mass and number to positive values and 
making the temperature 250 at the surface and decreasing to 150 K at TOA.
"""

#IMPORT NEEDED STUFF:
#======================
import numpy as np

# --- INPUT: DIMENSIONS ---
it= 1  #how many steps the model has taken. If first step, behaves differently to avoid uninitialized vals.
its=1  #lower index for column chunk. always 1
ite=1  #upper index for column chunk. use single column for tests.
kts=1  #1st vertical layer (model top)
kte=72 #top vertical layer (near surface)
ncat=1 #number of ice categories (will always be 1)
n_diag_2d = 1 #number of user-specified 2d diagnostic output fields (we aren't using any, just returns zeros)
n_diag_3d = 1 #number of user-specified 3d diagnostic output fields (we aren't using any, just returns zeros)
log_predictnc = True      # true = 2 moment, 1 = prescribed number conc
typediags_on = True       # compute in-depth precip types for output
model = 'GEM'             # model is used to decide if index=1 is surf. GEM has index=1 top like E3SM.

# --- INOUT VARIABLES ---
qc = np.zeros([ite,kte],np.float32)       #cloud liquid water mixing ratio    kg/kg
nc = np.ones([ite,kte],np.float32)*1e6    #cloud liquid drop number           #/kg
qr = np.zeros([ite,kte],np.float32)       #rain water mixing ratio            kg/kg
nr = np.ones([ite,kte],np.float32)*1e6    #rain drop number                   #/kg

#note: indexing for making some levels nonzero assumes there's only 1 column and 1 category.
qitot = np.zeros([ite,kte,ncat],np.float32)    #total ice mass mixing ratio   kg/kg
qitot[0,-20:-5,0]=1e-4*(1.-np.arange(15)/14.)  #max cld at ~700mb, decreasing to 0 at 900mb.
nitot = np.ones([ite,kte,ncat],np.float32)*1e6 #total ice number              #/kg
qirim = np.zeros([ite,kte,ncat],np.float32)    #rime ice mass mixing ratio    kg/kg
qirim[0,-20:-5,0]=1e-6*(1.-np.arange(15)/14.)  #max cld at ~700mb, decreasing to 0 at 900mb.
birim = np.zeros([ite,kte,ncat],np.float32)    #rime ice volume mixing ratio  m3/kg
birim[0,-20:-5,0]=1e-2                         #guess based on: m3/kg is 1/density and liquid water has a density of 1000 kg/m3
ssat = np.zeros([ite,kte],np.float32)          #supersaturation (qv - qs)     kg/kg

#qv goes to zero halfway through profile (to avoid condensate near model top)
tmp=-5e-4+np.reshape(1e-3/float(kte)\
                     *np.arange(kte),[1,kte]).astype(np.float32)
qv = np.where(tmp>0.,tmp,0.,)                                 # water vapor mixing ratio  kg kg-1
#qv[:,-20:-5]=5e-4  #make layer with qi saturated

#pres is actually an input variable, but needed here to compute theta.
pres = 100.+np.reshape(1e5/float(kte)\
                  *np.arange(kte),[1,kte]).astype(np.float32)  # pressure (note: must be >0 to avoid nans!)  Pa

#To get potential temperature, start by making absolute temperature
#vary between 150K at top of atmos and 273k at surface, then convert to potential temp.
T = (150.+np.reshape(110./float(kte)\
                *np.arange(kte),[1,kte])).astype(np.float32) 
Rd=287.    #Gas constant for dry air                  J/kg/K
cp=1004.   #heat constant of air at constant pressure J/kg
p0=100000. #reference pressure                        Pa
th=T*(p0/pres)**(Rd/cp)                                        # potential temperature           K


th_old = th.copy()                                           # theta @ beginning of timestep   K
qv_old = qv.copy()                                           # qv at beginning of timestep    kg kg-1

# --- INPUT VARIABLES ---
dt =np.float32(1800.)                                          # model time step                s
uzpl = np.zeros([ite,kte],np.float32)                          # vertical air velocity          m s-1

#compute vertical grid spacing dzq (in m) from pres and theta:
g=9.8      #gravity                                   m/s2
edge_pres=np.zeros([ite,kte+1],np.float32)
edge_pres[:,1:-1]=(pres[:,0:-1]+pres[:,1:])/2.
edge_pres[:,0]=max(0.,pres[:,0] - 0.5*(pres[:,1]-pres[:,0])/(1.-0.))
edge_pres[:,-1]=pres[:,-1] + 0.5*(pres[:,-1]-pres[:,-2])/1.0
dpres=edge_pres[:,1:]-edge_pres[:,0:-1]

#note dpres is >0 since index=1 is top of atmos. dz should also be pos, 
#so dropped minus sign in dp/dz=-rho*g
dzq  = Rd*T/(g*pres)*dpres  


# --- OUTPUT VARIABLES (included for reference) ---
#prt_liq   #liquid surface precip          m/s
#prt_sol     #frozen surface precip          m/s
#diag_ze     #equiv reflectivity             dBZ
#diag_effc   #liq effective radius           m
#diag_effi   #ice effective radius           m
#diag_vmi    #mass-weighted ice fall speed   m/s
#diag_di     #mean ice diameter              m
#diag_rhoi   #bulk density of ice            kg/m            
#diag_2d     #user-specified 2d outputs      N/A
#diag_3d     #user-specified 3d outputs      N/A
